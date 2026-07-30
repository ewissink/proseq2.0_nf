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

    output:
    path "*.bw",                   emit: bw
    path "${meta.id}.align.log",   emit: log
    path "${meta.id}.map.metrics", emit: metrics
    path "${meta.id}.sort.bam"

    script:
    def prefix = meta.id
    def single = meta.single_end ? 'true' : 'false'
    def se_out = meta.se_output ?: 'G'
    def rna5   = meta.rna5 ?: 'R1_5prime'
    def opp    = meta.opp  ? 'TRUE' : 'FALSE'
    def report = meta.report ?: 'rna_3prime'
    """
    set -eu
    echo "${prefix}" > ${prefix}.align.log

    # keep a copy of the sorted BAM (same name) for publishing
    if [ "${bam}" != "${prefix}.sort.bam" ]; then cp "${bam}" "${prefix}.sort.bam"; fi

    # gen_bed <map5 TRUE|FALSE> <out.bed.gz> : single-base BED of one reported RNA end
    gen_bed() {
        local map5="\$1" out="\$2"
        if [ "${single}" = "true" ]; then
            if [ "${se_out}" = "G" ]; then   # report read 5' end, same strand
                bedtools bamtobed -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$5>0){print \$0}' | awk 'BEGIN{OFS="\\t"} (\$6=="+"){print \$1,\$2,\$2+1,\$4,\$5,\$6}; (\$6=="-"){print \$1,\$3-1,\$3,\$4,\$5,\$6}' | gzip > "\$out"
            else                             # report read 5' end, strand flipped
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

    # metrics_from <raw.bed.gz> : write align.log counts + map.metrics from one track
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

    # emit_track <raw.bed.gz> <label> : rRNA/chrM filter -> raw + RPM plus/minus bigWigs
    emit_track() {
        local bed="\$1" lab="\$2" RC
        zcat "\$bed" | grep "rRNA\\|chrM" -v | grep "_" -v | sort-bed - | gzip > _track.nr.rs.bed.gz
        RC=\$(zcat _track.nr.rs.bed.gz | wc -l)
        bedtools genomecov -bg -i _track.nr.rs.bed.gz -g ${chrom_info} -strand + > _p.bedGraph
        bedtools genomecov -bg -i _track.nr.rs.bed.gz -g ${chrom_info} -strand - > _m.noinv.bedGraph
        awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,-1*\$4}' _m.noinv.bedGraph > _m.bedGraph
        awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,\$4*1000*1000/'\$RC'/1}' _p.bedGraph > _p.rpm.bedGraph
        awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,\$4*1000*1000/'\$RC'/1}' _m.bedGraph > _m.rpm.bedGraph
        bedGraphToBigWig _p.rpm.bedGraph ${chrom_info} ${prefix}\${lab}_plus.rpm.bw
        bedGraphToBigWig _m.rpm.bedGraph ${chrom_info} ${prefix}\${lab}_minus.rpm.bw
        bedGraphToBigWig _p.bedGraph     ${chrom_info} ${prefix}\${lab}_plus.bw
        bedGraphToBigWig _m.bedGraph     ${chrom_info} ${prefix}\${lab}_minus.bw
        rm -f _track.nr.rs.bed.gz _p.bedGraph _m.noinv.bedGraph _m.bedGraph _p.rpm.bedGraph _m.rpm.bedGraph
    }

    # --- decide which track(s) to emit ---
    if [ "${single}" = "true" ]; then
        gen_bed TRUE ${prefix}.bed.gz          # SE always reports the read 5' end
        metrics_from ${prefix}.bed.gz
        emit_track   ${prefix}.bed.gz ""
    elif [ "${report}" = "both" ]; then
        gen_bed TRUE  ${prefix}.5prime.bed.gz
        gen_bed FALSE ${prefix}.3prime.bed.gz
        metrics_from ${prefix}.5prime.bed.gz   # read/pair count is the same for either end
        emit_track   ${prefix}.5prime.bed.gz "_5prime"
        emit_track   ${prefix}.3prime.bed.gz "_3prime"
    else
        if [ "${report}" = "rna_5prime" ]; then M=TRUE; else M=FALSE; fi
        gen_bed \$M ${prefix}.bed.gz
        metrics_from ${prefix}.bed.gz
        emit_track   ${prefix}.bed.gz ""
    fi
    """

    stub:
    def prefix = meta.id
    def both   = (!meta.single_end && meta.report == 'both')
    def labels = both ? ['_5prime', '_3prime'] : ['']
    """
    for lab in ${labels.collect{ "\"${it}\"" }.join(' ')}; do
        touch ${prefix}\${lab}_plus.bw ${prefix}\${lab}_minus.bw ${prefix}\${lab}_plus.rpm.bw ${prefix}\${lab}_minus.rpm.bw
    done
    touch ${prefix}.align.log ${prefix}.sort.bam
    printf 'sample\\t%s\\nmappable\\t0\\nmappable_no_rrna\\t0\\n' "${prefix}" > ${prefix}.map.metrics
    """
}
