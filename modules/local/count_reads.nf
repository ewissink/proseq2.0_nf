process COUNT_READS {
    tag   "$meta.id"
    label 'process_low'
    conda "${projectDir}/env/environment.yml"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), stdout, emit: counts   // stdout = mapped-read count

    script:
    """
    samtools view -c ${bam}
    """

    stub:
    """
    echo 1000
    """
}
