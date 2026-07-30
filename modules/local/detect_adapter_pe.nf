process DETECT_ADAPTER_PE {
    tag   "$meta.id"
    label 'process_low'
    // own env: fastp's libdeflate conflicts with the pinned bedtools=2.28.0
    conda "conda-forge::python=3.11 bioconda::fastp=0.23.4"

    publishDir "${params.outdir}/qc/adapter", mode: params.publish_mode,
        pattern: "*.{fastp.json,adapter_check.txt}"

    input:
    tuple val(meta), path(read1), path(read2)

    output:
    path "${meta.id}.fastp.json",         emit: json
    path "${meta.id}.adapter_check.txt",  emit: check

    script:
    def prefix = meta.id
    // pipeline convention: R1's 3' adapter = --adapter2, R2's 3' adapter = --adapter1
    """
    set -eu
    fastp -i ${read1} -I ${read2} -o d1.fastq.gz -O d2.fastq.gz --detect_adapter_for_pe --thread ${task.cpus} \
        --disable_quality_filtering --disable_length_filtering \
        --json ${prefix}.fastp.json --html ${prefix}.fastp.html 2> ${prefix}.fastp.log
    rm -f d1.fastq.gz d2.fastq.gz

    adapter_check.py ${prefix}.fastp.json ${prefix} PE "${params.adapter2}" "${params.adapter1}"
    """

    stub:
    def prefix = meta.id
    """
    printf '{"adapter_cutting":{"read1_adapter_sequence":"${params.adapter2}","read2_adapter_sequence":"${params.adapter1}"}}' > ${prefix}.fastp.json
    adapter_check.py ${prefix}.fastp.json ${prefix} PE "${params.adapter2}" "${params.adapter1}"
    """
}
