process PREPROCESS_SE {
    tag   "$meta.id"
    label 'process_medium'
    conda "${projectDir}/env/environment.yml"

    publishDir "${params.outdir}/qc", mode: params.publish_mode,
        pattern: "*.{QC.log,prinseq-pcrDups.gd}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_QC_end.fastq.gz"), emit: reads
    path "${meta.id}.QC.log",                            emit: log
    path "${meta.id}.cutadapt.log",                      emit: cutadapt
    path "${meta.id}.reads.metrics",                     emit: metrics
    path "${meta.id}.prinseq-pcrDups.gd", optional: true

    script:
    def prefix = meta.id
    def n1     = (params.umi1 as int) + (params.add_b1 as int)
    def n2     = (params.umi2 as int) + (params.add_b2 as int)
    def dedup  = (((params.umi1 as int) > 0) || ((params.umi2 as int) > 0)) ? 'true' : 'false'
    def force  = params.force_deduplicate ? 'true' : 'false'
    """
    set -eu
    GD=${prefix}.prinseq-pcrDups.gd
    : > \$GD

    echo "${prefix}"                       >  ${prefix}.QC.log
    echo "Number of original input reads:" >> ${prefix}.QC.log
    echo \$(( \$(zcat ${reads} | wc -l) / 4 )) >> ${prefix}.QC.log

    # --- Remove 3' adapter; quality trim; strip 3' UMI/barcode ---
    cutadapt -a ${params.adapter_se} -e 0.10 --overlap 2 --output=${prefix}_trim.fastq --untrimmed-output=${prefix}_untrim.fastq ${reads} -j ${task.cpus} > ${prefix}.cutadapt.log
    cutadapt --cut -${n2} --minimum-length=10 ${prefix}_trim.fastq --output=${prefix}_trimN.fastq -q 20 -j ${task.cpus}
    cutadapt --minimum-length=10 ${prefix}_untrim.fastq --output=${prefix}_q20trim.fastq -q 20 -j ${task.cpus}

    cat ${prefix}_q20trim.fastq ${prefix}_trimN.fastq | paste - - - - | LC_ALL=C sort -k1,1 -S 2G | tr '\\t' '\\n' > ${prefix}_noadapt.fastq

    echo "Number of reads after adapter removal and QC:" >> ${prefix}.QC.log
    echo \$(( \$(cat ${prefix}_noadapt.fastq | wc -l) / 4 ))    >> ${prefix}.QC.log

    dedup_L=30
    if [ "${dedup}" = "true" ]; then
        # --- Remove PCR duplicates using the first \${dedup_L} nt ---
        seqtk seq -L \$dedup_L ${prefix}_noadapt.fastq > ${prefix}_l30.fastq
        prinseq-lite.pl -fastq ${prefix}_l30.fastq -derep 1 -out_format 3 -out_bad null -out_good ${prefix}_dedup -min_len 15 2>> \$GD
        awk '(NR%4==1){print substr(\$1,2)}' ${prefix}_dedup.fastq > ${prefix}_dedup.names
        seqtk subseq ${prefix}_noadapt.fastq ${prefix}_dedup.names > ${prefix}_dedup_withBarcode.fastq
        prinseq-lite.pl -trim_left ${n1} -fastq ${prefix}_dedup_withBarcode.fastq -out_format 3 -out_bad null -out_good ${prefix}_dedup_BarcodeRemoved 2>> \$GD
        prinseq-lite.pl -min_len 15 -fastq ${prefix}_dedup_BarcodeRemoved.fastq -out_format 3 -out_bad null -out_good ${prefix}_QC_end 2>> \$GD
        echo "Number of reads after PCR duplicate removal and QC:" >> ${prefix}.QC.log
        echo \$(( \$(cat ${prefix}_QC_end.fastq | wc -l) / 4 ))    >> ${prefix}.QC.log
    elif [ "${force}" = "true" ]; then
        # --- Force dedup on full-length reads (no UMI) ---
        prinseq-lite.pl -derep 1 -fastq ${prefix}_noadapt.fastq -out_format 3 -out_bad null -out_good ${prefix}_dedup_withBarcode 2>> \$GD
        prinseq-lite.pl -trim_left ${n1} -fastq ${prefix}_dedup_withBarcode.fastq -out_format 3 -out_bad null -out_good ${prefix}_dedup_BarcodeRemoved 2>> \$GD
        prinseq-lite.pl -min_len 15 -fastq ${prefix}_dedup_BarcodeRemoved.fastq -out_format 3 -out_bad null -out_good ${prefix}_QC_end 2>> \$GD
        echo "Number of reads after PCR duplicate removal and QC:" >> ${prefix}.QC.log
        echo \$(( \$(cat ${prefix}_QC_end.fastq | wc -l) / 4 ))    >> ${prefix}.QC.log
    else
        # --- No dedup: strip additional barcode only ---
        prinseq-lite.pl -trim_left ${params.add_b1} -fastq ${prefix}_noadapt.fastq -out_format 3 -out_bad null -out_good ${prefix}_BarcodeRemoved 2>> \$GD
        prinseq-lite.pl -min_len 15 -fastq ${prefix}_BarcodeRemoved.fastq -out_format 3 -out_bad null -out_good ${prefix}_QC_end 2>> \$GD
        echo "Number of reads after final QC:" >> ${prefix}.QC.log
        echo \$(( \$(cat ${prefix}_QC_end.fastq | wc -l) / 4 )) >> ${prefix}.QC.log
    fi

    # --- Optional dataset-wide length cutoff for mapping ---
    if [ "${params.map_length}" -ne 0 ]; then
        seqtk seq -L ${params.map_length} ${prefix}_QC_end.fastq > ${prefix}_QC_end.trimmed.fastq
        mv ${prefix}_QC_end.trimmed.fastq ${prefix}_QC_end.fastq
    fi

    gzip -nf ${prefix}_QC_end.fastq

    # --- machine-readable metrics for MultiQC ---
    IN=\$(( \$(zcat ${reads} | wc -l) / 4 ))
    PQ=\$(( \$(zcat ${prefix}_QC_end.fastq.gz | wc -l) / 4 ))
    printf 'sample\\t%s\\nseq_type\\tSE\\ninput\\t%s\\npass_qc\\t%s\\n' "${prefix}" "\$IN" "\$PQ" > ${prefix}.reads.metrics
    """

    stub:
    def prefix = meta.id
    """
    echo '' | gzip > ${prefix}_QC_end.fastq.gz
    touch ${prefix}.QC.log ${prefix}.cutadapt.log ${prefix}.prinseq-pcrDups.gd
    printf 'sample\\t%s\\nseq_type\\tSE\\ninput\\t0\\npass_qc\\t0\\n' "${prefix}" > ${prefix}.reads.metrics
    """
}
