process SORTMERNA_SE {
    tag   "$meta.id"
    label 'process_medium'
    conda "${projectDir}/env/environment.yml"

    publishDir "${params.outdir}/qc/sortmerna", mode: params.publish_mode, pattern: "*.sortmerna.log"

    input:
    tuple val(meta), path(reads)
    path  rrna_refs

    output:
    tuple val(meta), path("${meta.id}_nonrRNA.fastq.gz"), emit: reads
    path "${meta.id}.sortmerna.log",                       emit: log
    path "${meta.id}.sortmerna.metrics",                   emit: metrics

    script:
    def prefix = meta.id
    """
    set -eu
    REF_ARGS=""
    for r in ${rrna_refs}; do REF_ARGS="\$REF_ARGS --ref \$r"; done

    sortmerna \$REF_ARGS --reads ${reads} --threads ${task.cpus} \
        --workdir "\$PWD/smr" --fastx --aligned rRNA_reads --other non_rRNA_reads

    # normalize the kept (non-rRNA) reads to a gzipped file
    if ls non_rRNA_reads.f*q.gz >/dev/null 2>&1; then
        mv non_rRNA_reads.f*q.gz ${prefix}_nonrRNA.fastq.gz
    else
        mv non_rRNA_reads.f*q ${prefix}_nonrRNA.fastq && gzip ${prefix}_nonrRNA.fastq
    fi
    mv rRNA_reads.log ${prefix}.sortmerna.log

    KEPT=\$(( \$(zcat ${prefix}_nonrRNA.fastq.gz | wc -l) / 4 ))
    printf 'sample\\t%s\\npost_rrna\\t%s\\n' "${prefix}" "\$KEPT" > ${prefix}.sortmerna.metrics
    """

    stub:
    def prefix = meta.id
    """
    echo '' | gzip > ${prefix}_nonrRNA.fastq.gz
    touch ${prefix}.sortmerna.log
    printf 'sample\\t%s\\npost_rrna\\t0\\n' "${prefix}" > ${prefix}.sortmerna.metrics
    """
}
