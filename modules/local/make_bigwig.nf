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
    tuple val(meta), path("${meta.id}_plus.bw"), path("${meta.id}_minus.bw"), emit: bw
    path "${meta.id}_plus.rpm.bw"
    path "${meta.id}_minus.rpm.bw"
    path "${meta.id}.align.log",   emit: log
    path "${meta.id}.map.metrics", emit: metrics
    path "${meta.id}.sort.bam"

    script:
    def prefix = meta.id
    def single = meta.single_end ? 'true' : 'false'
    def se_out = meta.se_output ?: 'G'
    def rna5   = meta.rna5 ?: 'R1_5prime'
    def opp    = meta.opp  ? 'TRUE' : 'FALSE'
    def map5   = meta.map5 ? 'TRUE' : 'FALSE'
    """
    set -eu
    echo "${prefix}" > ${prefix}.align.log

    # keep a copy of the sorted BAM (same name) for publishing
    if [ "${bam}" != "${prefix}.sort.bam" ]; then cp "${bam}" "${prefix}.sort.bam"; fi

    # ------------------------------------------------------------------
    # 1) BAM -> single-base BED of the reported RNA end
    # ------------------------------------------------------------------
    if [ "${single}" = "true" ]; then
        if [ "${se_out}" = "G" ]; then   # GRO-seq: 5' end of RNA, same strand
            bedtools bamtobed -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$5 > 0){print \$0}' | awk 'BEGIN{OFS="\\t"} (\$6=="+"){print \$1,\$2,\$2+1,\$4,\$5,\$6}; (\$6=="-"){print \$1,\$3-1,\$3,\$4,\$5,\$6}' | gzip > ${prefix}.bed.gz
        else                             # PRO-seq: 3' end, strand flipped
            bedtools bamtobed -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$5 > 0){print \$0}' | awk 'BEGIN{OFS="\\t"} (\$6=="+"){print \$1,\$2,\$2+1,\$4,\$5,"-"}; (\$6=="-"){print \$1,\$3-1,\$3,\$4,\$5,"+"}' | gzip > ${prefix}.bed.gz
        fi
    else
        if [ "${rna5}" = "R1_5prime" ]; then
            if [ "${opp}" = "FALSE" ]; then
                if [ "${map5}" = "TRUE" ]; then
                    bedtools bamtobed -bedpe -mate1 -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$9=="+"){print \$1,\$2,\$2+1,\$7,\$8,\$9}; (\$9=="-"){print \$1,\$3-1,\$3,\$7,\$8,\$9}' | gzip > ${prefix}.bed.gz
                else
                    bedtools bamtobed -bedpe -mate1 -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$10=="-"){print \$1,\$6-1,\$6,\$7,\$8,\$9}; (\$10=="+"){print \$1,\$5,\$5+1,\$7,\$8,\$9}' | gzip > ${prefix}.bed.gz
                fi
            else
                if [ "${map5}" = "TRUE" ]; then
                    bedtools bamtobed -bedpe -mate1 -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$9=="+"){print \$1,\$2,\$2+1,\$7,\$8,\$10}; (\$9=="-"){print \$1,\$3-1,\$3,\$7,\$8,\$10}' | gzip > ${prefix}.bed.gz
                else
                    bedtools bamtobed -bedpe -mate1 -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$10=="-"){print \$1,\$6-1,\$6,\$7,\$8,\$10}; (\$10=="+"){print \$1,\$5,\$5+1,\$7,\$8,\$10}' | gzip > ${prefix}.bed.gz
                fi
            fi
        else   # rna5 == R2_5prime
            if [ "${opp}" = "FALSE" ]; then
                if [ "${map5}" = "TRUE" ]; then
                    bedtools bamtobed -bedpe -mate1 -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$10=="+"){print \$1,\$5,\$5+1,\$7,\$8,\$10}; (\$10=="-"){print \$1,\$6-1,\$6,\$7,\$8,\$10}' | gzip > ${prefix}.bed.gz
                else
                    bedtools bamtobed -bedpe -mate1 -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$9=="+"){print \$1,\$2,\$2+1,\$7,\$8,\$10}; (\$9=="-"){print \$1,\$3-1,\$3,\$7,\$8,\$10}' | gzip > ${prefix}.bed.gz
                fi
            else
                if [ "${map5}" = "TRUE" ]; then
                    bedtools bamtobed -bedpe -mate1 -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$10=="+"){print \$1,\$5,\$5+1,\$7,\$8,\$9}; (\$10=="-"){print \$1,\$6-1,\$6,\$7,\$8,\$9}' | gzip > ${prefix}.bed.gz
                else
                    bedtools bamtobed -bedpe -mate1 -i ${bam} 2> kill.warnings | awk 'BEGIN{OFS="\\t"} (\$9=="+"){print \$1,\$2,\$2+1,\$7,\$8,\$9}; (\$9=="-"){print \$1,\$3-1,\$3,\$7,\$8,\$9}' | gzip > ${prefix}.bed.gz
                fi
            fi
        fi
    fi

    # ------------------------------------------------------------------
    # 2) Drop rRNA / chrM / non-primary contigs; sort
    # ------------------------------------------------------------------
    MAP=\$(zcat ${prefix}.bed.gz | wc -l)
    echo "Number of mappable reads:"     >> ${prefix}.align.log
    echo \$MAP                            >> ${prefix}.align.log

    zcat ${prefix}.bed.gz | grep "rRNA\\|chrM" -v | grep "_" -v | sort-bed - | gzip > ${prefix}.nr.rs.bed.gz
    readCount=\$(zcat ${prefix}.nr.rs.bed.gz | wc -l)
    echo "Number of mappable reads (excluding rRNA):" >> ${prefix}.align.log
    echo \$readCount                                   >> ${prefix}.align.log

    printf 'sample\\t%s\\nmappable\\t%s\\nmappable_no_rrna\\t%s\\n' "${prefix}" "\$MAP" "\$readCount" > ${prefix}.map.metrics

    # ------------------------------------------------------------------
    # 3) bedGraph -> bigWig  (raw and RPM-normalized, plus/minus)
    # ------------------------------------------------------------------
    bedtools genomecov -bg -i ${prefix}.nr.rs.bed.gz -g ${chrom_info} -strand + > ${prefix}_plus.bedGraph
    bedtools genomecov -bg -i ${prefix}.nr.rs.bed.gz -g ${chrom_info} -strand - > ${prefix}_minus.noinv.bedGraph
    awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,-1*\$4}' ${prefix}_minus.noinv.bedGraph > ${prefix}_minus.bedGraph

    awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,\$4*1000*1000/'\$readCount'/1}' ${prefix}_plus.bedGraph  > ${prefix}_plus.rpm.bedGraph
    awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,\$4*1000*1000/'\$readCount'/1}' ${prefix}_minus.bedGraph > ${prefix}_minus.rpm.bedGraph

    bedGraphToBigWig ${prefix}_plus.rpm.bedGraph  ${chrom_info} ${prefix}_plus.rpm.bw
    bedGraphToBigWig ${prefix}_minus.rpm.bedGraph ${chrom_info} ${prefix}_minus.rpm.bw
    bedGraphToBigWig ${prefix}_plus.bedGraph      ${chrom_info} ${prefix}_plus.bw
    bedGraphToBigWig ${prefix}_minus.bedGraph     ${chrom_info} ${prefix}_minus.bw
    """

    stub:
    def prefix = meta.id
    """
    touch ${prefix}_plus.bw ${prefix}_minus.bw ${prefix}_plus.rpm.bw ${prefix}_minus.rpm.bw
    touch ${prefix}.align.log ${prefix}.sort.bam
    printf 'sample\\t%s\\nmappable\\t0\\nmappable_no_rrna\\t0\\n' "${prefix}" > ${prefix}.map.metrics
    """
}
