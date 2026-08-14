process MULTIQC {
    label 'process_low'
    conda "${projectDir}/env/environment.yml"

    publishDir "${params.outdir}/multiqc", mode: params.publish_mode

    input:
    path multiqc_files, stageAs: "?/*"

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data

    script:
    // default output is already multiqc_report.html + multiqc_data/; passing -n
    // renames the data dir to <name>_data, which broke the declared output.
    """
    multiqc -f .
    """

    stub:
    """
    mkdir -p multiqc_data
    touch multiqc_report.html
    """
}
