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
    """
    multiqc -f -n multiqc_report.html .
    """

    stub:
    """
    mkdir -p multiqc_data
    touch multiqc_report.html
    """
}
