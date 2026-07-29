process CAT_FASTQ_PE {
    tag   "$meta.id"
    label 'process_low'
    conda "${projectDir}/env/environment.yml"

    input:
    tuple val(meta), path(reads1, stageAs: "input1*/*"), path(reads2, stageAs: "input2*/*")

    output:
    tuple val(meta), path("${meta.id}.merged_R1.fastq.gz"), path("${meta.id}.merged_R2.fastq.gz"), emit: reads

    script:
    // R1s and R2s must be concatenated in the SAME order to keep mates paired
    """
    cat ${reads1} > ${meta.id}.merged_R1.fastq.gz
    cat ${reads2} > ${meta.id}.merged_R2.fastq.gz
    """

    stub:
    """
    cat ${reads1} > ${meta.id}.merged_R1.fastq.gz
    cat ${reads2} > ${meta.id}.merged_R2.fastq.gz
    """
}
