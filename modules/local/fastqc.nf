process FASTQC {
    tag   "${meta.id}:${stage}"
    label 'process_low'
    conda "${projectDir}/env/environment.yml"

    publishDir { "${params.outdir}/qc/fastqc/${stage}" }, mode: params.publish_mode

    input:
    tuple val(meta), path(reads)
    val   stage              // 'raw' or 'trim' — for tag / publish subdir only

    output:
    path "*_fastqc.zip",  emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc --threads ${task.cpus} --quiet ${reads}
    """

    stub:
    """
    touch ${meta.id}_${stage}_fastqc.zip ${meta.id}_${stage}_fastqc.html
    """
}
