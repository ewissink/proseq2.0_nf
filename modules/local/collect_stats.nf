process COLLECT_STATS {
    label 'process_low'
    conda "${projectDir}/env/environment.yml"

    publishDir "${params.outdir}/qc", mode: params.publish_mode

    input:
    path metrics

    output:
    path "proseq_read_stats_mqc.tsv", emit: mqc

    script:
    """
    collect_stats.py . > proseq_read_stats_mqc.tsv
    """

    stub:
    """
    collect_stats.py . > proseq_read_stats_mqc.tsv
    """
}
