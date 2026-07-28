process SORTMERNA_PE {
    tag   "$meta.id"
    label 'process_medium'
    conda "${projectDir}/env/environment.yml"

    publishDir "${params.outdir}/qc/sortmerna", mode: params.publish_mode, pattern: "*.sortmerna.log"

    input:
    tuple val(meta), path(read1), path(read2)
    path  rrna_refs

    output:
    tuple val(meta), path("${meta.id}_nonrRNA_1.fastq.gz"), path("${meta.id}_nonrRNA_2.fastq.gz"), emit: reads
    path "${meta.id}.sortmerna.log",     emit: log
    path "${meta.id}.sortmerna.metrics", emit: metrics

    script:
    def prefix = meta.id
    """
    set -eu
    REF_ARGS=""
    for r in ${rrna_refs}; do REF_ARGS="\$REF_ARGS --ref \$r"; done

    # --paired_in: drop the whole pair if EITHER mate is rRNA; --out2: split fwd/rev
    sortmerna \$REF_ARGS --reads ${read1} --reads ${read2} --paired_in --out2 --threads ${task.cpus} \
        --workdir "\$PWD/smr" --fastx --aligned rRNA_reads --other non_rRNA_reads

    if ls non_rRNA_reads_fwd.f*q.gz >/dev/null 2>&1; then
        mv non_rRNA_reads_fwd.f*q.gz ${prefix}_nonrRNA_1.fastq.gz
        mv non_rRNA_reads_rev.f*q.gz ${prefix}_nonrRNA_2.fastq.gz
    else
        mv non_rRNA_reads_fwd.f*q ${prefix}_nonrRNA_1.fastq && gzip ${prefix}_nonrRNA_1.fastq
        mv non_rRNA_reads_rev.f*q ${prefix}_nonrRNA_2.fastq && gzip ${prefix}_nonrRNA_2.fastq
    fi
    mv rRNA_reads.log ${prefix}.sortmerna.log

    KEPT=\$(( \$(zcat ${prefix}_nonrRNA_1.fastq.gz | wc -l) / 4 ))
    printf 'sample\\t%s\\npost_rrna\\t%s\\n' "${prefix}" "\$KEPT" > ${prefix}.sortmerna.metrics
    """

    stub:
    def prefix = meta.id
    """
    echo '' | gzip > ${prefix}_nonrRNA_1.fastq.gz
    echo '' | gzip > ${prefix}_nonrRNA_2.fastq.gz
    touch ${prefix}.sortmerna.log
    printf 'sample\\t%s\\npost_rrna\\t0\\n' "${prefix}" > ${prefix}.sortmerna.metrics
    """
}
