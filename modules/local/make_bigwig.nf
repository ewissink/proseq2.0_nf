process MAKE_BIGWIG {
    tag   "$meta.id"
    label 'process_medium'
    conda "${projectDir}/env/environment.yml"

    publishDir "${params.outdir}/bam",    mode: params.publish_mode, pattern: "*.sort.bam"
    publishDir "${params.outdir}/bigwig", mode: params.publish_mode, pattern: "*.bw"
    publishDir "${params.outdir}/qc",     mode: params.publish_mode, pattern: "*.align.log"

    input:
    tuple val(meta), path(bam)
    path  chrom_info
    val   norm            // 'standard' (raw + rpm) | 'subsample' (.subsample track only)

    output:
    path "*.bw",                   emit: bw
    path "${meta.id}.align.log",   emit: log,     optional: true
    path "${meta.id}.map.metrics", emit: metrics, optional: true
    path "${meta.id}.sort.bam", optional: true

    script:
    def prefix = meta.id
    def single = meta.single_end ? 'true' : 'false'
    def se_out = meta.se_output ?: 'G'
    def rna5   = meta.rna5 ?: 'R1_5prime'
    def opp    = meta.opp  ? 'TRUE' : 'FALSE'
    def report = meta.report ?: 'rna_3prime'
    """
    set -eu

    # gen_bed <map5 TRUE|FALSE> <out.bed.gz> : single-base BED of one reported RNA end
    gen_bed() {
        local map5="\$1" out="\$2"
        if [ "${single}" = "true" ]; then
            if [ "${se_out}" = "G" ]; then
                bedtools bamtobed -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$5>0){print \$0}' | awk 'BEGIN{OFS="\\t"} (\$6=="+"){print \$1,\$2,\$2+1,\$4,\$5,\$6}; (\$6=="-"){print \$1,\$3-1,\$3,\$4,\$5,\$6}' | gzip > "\$out"
            else
                bedtools bamtobed -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$5>0){print \$0}' | awk 'BEGIN{OFS="\\t"} (\$6=="+"){print \$1,\$2,\$2+1,\$4,\$5,"-"}; (\$6=="-"){print \$1,\$3-1,\$3,\$4,\$5,"+"}' | gzip > "\$out"
            fi
        elif [ "${rna5}" = "R1_5prime" ]; then
            if [ "${opp}" = "FALSE" ]; then
                if [ "\$map5" = "TRUE" ]; then
                    bedtools bamtobed -bedpe -mate1 -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$9=="+"){print \$1,\$2,\$2+1,\$7,\$8,\$9}; (\$9=="-"){print \$1,\$3-1,\$3,\$7,\$8,\$9}' | gzip > "\$out"
                else
                    bedtools bamtobed -bedpe -mate1 -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$10=="-"){print \$1,\$6-1,\$6,\$7,\$8,\$9}; (\$10=="+"){print \$1,\$5,\$5+1,\$7,\$8,\$9}' | gzip > "\$out"
                fi
            else
                if [ "\$map5" = "TRUE" ]; then
                    bedtools bamtobed -bedpe -mate1 -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$9=="+"){print \$1,\$2,\$2+1,\$7,\$8,\$10}; (\$9=="-"){print \$1,\$3-1,\$3,\$7,\$8,\$10}' | gzip > "\$out"
                else
                    bedtools bamtobed -bedpe -mate1 -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$10=="-"){print \$1,\$6-1,\$6,\$7,\$8,\$10}; (\$10=="+"){print \$1,\$5,\$5+1,\$7,\$8,\$10}' | gzip > "\$out"
                fi
            fi
        else   # rna5 == R2_5prime
            if [ "${opp}" = "FALSE" ]; then
                if [ "\$map5" = "TRUE" ]; then
                    bedtools bamtobed -bedpe -mate1 -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$10=="+"){print \$1,\$5,\$5+1,\$7,\$8,\$10}; (\$10=="-"){print \$1,\$6-1,\$6,\$7,\$8,\$10}' | gzip > "\$out"
                else
                    bedtools bamtobed -bedpe -mate1 -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$9=="+"){print \$1,\$2,\$2+1,\$7,\$8,\$10}; (\$9=="-"){print \$1,\$3-1,\$3,\$7,\$8,\$10}' | gzip > "\$out"
                fi
            else
                if [ "\$map5" = "TRUE" ]; then
                    bedtools bamtobed -bedpe -mate1 -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$10=="+"){print \$1,\$5,\$5+1,\$7,\$8,\$9}; (\$10=="-"){print \$1,\$6-1,\$6,\$7,\$8,\$9}' | gzip > "\$out"
                else
                    bedtools bamtobed -bedpe -mate1 -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$9=="+"){print \$1,\$2,\$2+1,\$7,\$8,\$9}; (\$9=="-"){print \$1,\$3-1,\$3,\$7,\$8,\$9}' | gzip > "\$out"
                fi
            fi
        fi
    }

    # metrics_from <raw.bed.gz> : align.log counts + map.metrics (standard runs only)
    metrics_from() {
        local bed="\$1" MAP RC
        MAP=\$(zcat "\$bed" | wc -l)
        echo "Number of mappable reads:"                  >> ${prefix}.align.log
        echo \$MAP                                          >> ${prefix}.align.log
        zcat "\$bed" | grep "rRNA\\|chrM" -v | grep "_" -v | sort-bed - | gzip > _metrics.nr.rs.bed.gz
        RC=\$(zcat _metrics.nr.rs.bed.gz | wc -l)
        echo "Number of mappable reads (excluding rRNA):" >> ${prefix}.align.log
        echo \$RC                                           >> ${prefix}.align.log
        printf 'sample\\t%s\\nmappable\\t%s\\nmappable_no_rrna\\t%s\\n' "${prefix}" "\$MAP" "\$RC" > ${prefix}.map.metrics
        rm -f _metrics.nr.rs.bed.gz
    }

    # emit_track <raw.bed.gz> <endlabel> : fwd/rev bigWigs (raw + rpm, or .subsample)
    emit_track() {
        local bed="\$1" end="\$2" RC
        zcat "\$bed" | grep "rRNA\\|chrM" -v | grep "_" -v | sort-bed - | gzip > _track.nr.rs.bed.gz
        RC=\$(zcat _track.nr.rs.bed.gz | wc -l)
        bedtools genomecov -bg -i _track.nr.rs.bed.gz -g ${chrom_info} -strand + > _fwd.bedGraph
        bedtools genomecov -bg -i _track.nr.rs.bed.gz -g ${chrom_info} -strand - > _rev.noinv.bedGraph
        awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,-1*\$4}' _rev.noinv.bedGraph > _rev.bedGraph
        if [ "${norm}" = "subsample" ]; then
            bedGraphToBigWig _fwd.bedGraph ${chrom_info} ${prefix}\${end}_fwd.subsample.bw
            bedGraphToBigWig _rev.bedGraph ${chrom_info} ${prefix}\${end}_rev.subsample.bw
        else
            bedGraphToBigWig _fwd.bedGraph ${chrom_info} ${prefix}\${end}_fwd.bw
            bedGraphToBigWig _rev.bedGraph ${chrom_info} ${prefix}\${end}_rev.bw
            awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,\$4*1000*1000/'\$RC'/1}' _fwd.bedGraph > _fwd.rpm.bedGraph
            awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,\$4*1000*1000/'\$RC'/1}' _rev.bedGraph > _rev.rpm.bedGraph
            bedGraphToBigWig _fwd.rpm.bedGraph ${chrom_info} ${prefix}\${end}_fwd.rpm.bw
            bedGraphToBigWig _rev.rpm.bedGraph ${chrom_info} ${prefix}\${end}_rev.rpm.bw
        fi
        rm -f _track.nr.rs.bed.gz _fwd.bedGraph _rev.noinv.bedGraph _rev.bedGraph _fwd.rpm.bedGraph _rev.rpm.bedGraph
    }

    # standard runs write the align.log + a copy of the sorted BAM; subsample runs don't
    if [ "${norm}" != "subsample" ]; then
        echo "${prefix}" > ${prefix}.align.log
        if [ "${bam}" != "${prefix}.sort.bam" ]; then cp "${bam}" "${prefix}.sort.bam"; fi
    fi

    if [ "${single}" = "true" ]; then
        gen_bed TRUE ${prefix}.bed.gz
        [ "${norm}" != "subsample" ] && metrics_from ${prefix}.bed.gz || true
        emit_track ${prefix}.bed.gz ""
    elif [ "${report}" = "both" ]; then
        gen_bed TRUE  ${prefix}.5prime.bed.gz
        gen_bed FALSE ${prefix}.3prime.bed.gz
        [ "${norm}" != "subsample" ] && metrics_from ${prefix}.5prime.bed.gz || true
        emit_track ${prefix}.5prime.bed.gz "_5prime"
        emit_track ${prefix}.3prime.bed.gz "_3prime"
    else
        if [ "${report}" = "rna_5prime" ]; then M=TRUE; else M=FALSE; fi
        gen_bed \$M ${prefix}.bed.gz
        [ "${norm}" != "subsample" ] && metrics_from ${prefix}.bed.gz || true
        emit_track ${prefix}.bed.gz ""
    fi
    """

    stub:
    def prefix = meta.id
    def both   = (!meta.single_end && meta.report == 'both')
    def ends   = both ? ['_5prime', '_3prime'] : ['']
    """
    for end in ${ends.collect{ "\"${it}\"" }.join(' ')}; do
        if [ "${norm}" = "subsample" ]; then
            touch ${prefix}\${end}_fwd.subsample.bw ${prefix}\${end}_rev.subsample.bw
        else
            touch ${prefix}\${end}_fwd.bw ${prefix}\${end}_rev.bw ${prefix}\${end}_fwd.rpm.bw ${prefix}\${end}_rev.rpm.bw
        fi
    done
    if [ "${norm}" != "subsample" ]; then
        touch ${prefix}.align.log ${prefix}.sort.bam
        printf 'sample\\t%s\\nmappable\\t0\\nmappable_no_rrna\\t0\\n' "${prefix}" > ${prefix}.map.metrics
    fi
    """
}
