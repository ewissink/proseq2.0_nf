process CAT_FASTQ_SE {
    tag   "$meta.id"
    label 'process_low'
    conda "${projectDir}/env/environment.yml"

    input:
    tuple val(meta), path(reads, stageAs: "input*/*")   // list of R1 files

    output:
    tuple val(meta), path("${meta.id}.merged.fastq.gz"), emit: reads

    script:
    // gzip streams concatenate directly; read order doesn't affect the output
    """
    cat ${reads} > ${meta.id}.merged.fastq.gz
    """

    stub:
    """
    cat ${reads} > ${meta.id}.merged.fastq.gz
    """
}
