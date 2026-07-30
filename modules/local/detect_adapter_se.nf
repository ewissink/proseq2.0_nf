process DETECT_ADAPTER_SE {
    tag   "$meta.id"
    label 'process_low'
    // own env: fastp's libdeflate conflicts with the pinned bedtools=2.28.0
    conda "conda-forge::python=3.11 bioconda::fastp=0.23.4"

    publishDir "${params.outdir}/qc/adapter", mode: params.publish_mode,
        pattern: "*.{fastp.json,adapter_check.txt}"

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta.id}.fastp.json",         emit: json
    path "${meta.id}.adapter_check.txt",  emit: check

    script:
    def prefix = meta.id
    """
    set -eu
    # fastp detects the adapter de-novo (report only; trimmed output is discarded)
    fastp -i ${reads} -o discard.fastq.gz --thread ${task.cpus} \
        --disable_quality_filtering --disable_length_filtering \
        --json ${prefix}.fastp.json --html ${prefix}.fastp.html 2> ${prefix}.fastp.log
    rm -f discard.fastq.gz

    adapter_check.py ${prefix}.fastp.json ${prefix} SE "${params.adapter_se}"
    """

    stub:
    def prefix = meta.id
    """
    printf '{"adapter_cutting":{"read1_adapter_sequence":"${params.adapter_se}"}}' > ${prefix}.fastp.json
    adapter_check.py ${prefix}.fastp.json ${prefix} SE "${params.adapter_se}"
    """
}
