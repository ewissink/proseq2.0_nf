process PREPROCESS_PE {
    tag   "$meta.id"
    label 'process_medium'
    conda "${projectDir}/env/environment.yml"

    publishDir "${params.outdir}/qc", mode: params.publish_mode,
        pattern: "*.{QC.log,prinseq-pcrDups.gd}"

    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta), path("${meta.id}_QC_end_1.fastq.gz"), path("${meta.id}_QC_end_2.fastq.gz"), emit: reads
    path "${meta.id}.QC.log",                              emit: log
    path "${meta.id}_R*.cutadapt.log",                     emit: cutadapt
    path "${meta.id}.reads.metrics",                       emit: metrics
    path "${meta.id}.prinseq-pcrDups.gd", optional: true

    script:
    def prefix = meta.id
    def n1     = (params.umi1 as int) + (params.add_b1 as int)   // removed from 3' of R2
    def n2     = (params.umi2 as int) + (params.add_b2 as int)   // removed from 3' of R1
    def dedup  = (((params.umi1 as int) > 0) || ((params.umi2 as int) > 0)) ? 'true' : 'false'
    def force  = params.force_deduplicate ? 'true' : 'false'
    """
    set -eu
    GD=${prefix}.prinseq-pcrDups.gd
    : > \$GD

    echo "Number of original input reads:"                     >  ${prefix}.QC.log
    echo "R1: \$(( \$(zcat ${read1} | wc -l) / 4 ))"           >> ${prefix}.QC.log
    echo "R2: \$(( \$(zcat ${read2} | wc -l) / 4 ))"           >> ${prefix}.QC.log

    # --- Remove 3' adapters (ADAPT2 from R1, ADAPT1 from R2) ---
    cutadapt -a ${params.adapter2} -e 0.10 --overlap 2 --output=${prefix}_trim_R1.fastq --untrimmed-output=${prefix}_untrim_R1.fastq ${read1} -j ${task.cpus} > ${prefix}_R1.cutadapt.log
    cutadapt -a ${params.adapter1} -e 0.10 --overlap 2 --output=${prefix}_trim_R2.fastq --untrimmed-output=${prefix}_untrim_R2.fastq ${read2} -j ${task.cpus} > ${prefix}_R2.cutadapt.log

    # R1: remove UMI2+ADD_B2 from 3' end, quality trim
    cutadapt --cut -${n2} --minimum-length=10 ${prefix}_trim_R1.fastq --output=${prefix}_trimN_R1.fastq -q 20 -j ${task.cpus}
    cutadapt --minimum-length=10 ${prefix}_untrim_R1.fastq --output=${prefix}_q20_R1.fastq -q 20 -j ${task.cpus}
    cat ${prefix}_q20_R1.fastq ${prefix}_trimN_R1.fastq | paste - - - - | LC_ALL=C sort -k1,1 -S 2G | tr '\\t' '\\n' > ${prefix}_noadapt_R1.fastq

    # R2: remove UMI1+ADD_B1 from 3' end, quality trim
    cutadapt --cut -${n1} --minimum-length=10 ${prefix}_trim_R2.fastq --output=${prefix}_trimN_R2.fastq -q 20 -j ${task.cpus}
    cutadapt --minimum-length=10 ${prefix}_untrim_R2.fastq --output=${prefix}_q20_R2.fastq -q 20 -j ${task.cpus}
    cat ${prefix}_q20_R2.fastq ${prefix}_trimN_R2.fastq | paste - - - - | LC_ALL=C sort -k1,1 -S 2G | tr '\\t' '\\n' > ${prefix}_noadapt_R2.fastq

    echo "Number of reads after adapter removal and QC (R1):" >> ${prefix}.QC.log
    echo \$(( \$(cat ${prefix}_noadapt_R1.fastq | wc -l) / 4 )) >> ${prefix}.QC.log
    echo "Number of reads after adapter removal and QC (R2):" >> ${prefix}.QC.log
    echo \$(( \$(cat ${prefix}_noadapt_R2.fastq | wc -l) / 4 )) >> ${prefix}.QC.log

    dedup_L=30
    if [ "${dedup}" = "true" ]; then
        # --- Dedup using the first \${dedup_L} nt of the pair ---
        seqtk seq -L \$dedup_L ${prefix}_noadapt_R1.fastq > ${prefix}_l30_R1.fastq
        seqtk seq -L \$dedup_L ${prefix}_noadapt_R2.fastq > ${prefix}_l30_R2.fastq
        prinseq-lite.pl -fastq ${prefix}_l30_R1.fastq -fastq2 ${prefix}_l30_R2.fastq -derep 1 -out_format 3 -out_bad null -out_good ${prefix}_dedup -min_len 15 2>> \$GD
        awk '(NR%4==1){print substr(\$1,2)}' ${prefix}_dedup_1.fastq > ${prefix}_dedup.names
        seqtk subseq ${prefix}_noadapt_R1.fastq ${prefix}_dedup.names > ${prefix}_dwb_1.fastq
        seqtk subseq ${prefix}_noadapt_R2.fastq ${prefix}_dedup.names > ${prefix}_dwb_2.fastq
        prinseq-lite.pl -trim_left ${n1} -fastq ${prefix}_dwb_1.fastq -out_format 3 -out_bad null -out_good ${prefix}_dbr_1 2>> \$GD
        prinseq-lite.pl -trim_left ${n2} -fastq ${prefix}_dwb_2.fastq -out_format 3 -out_bad null -out_good ${prefix}_dbr_2 2>> \$GD
        prinseq-lite.pl -min_len 15 -fastq ${prefix}_dbr_1.fastq -fastq2 ${prefix}_dbr_2.fastq -out_format 3 -out_bad null -out_good ${prefix}_QC_end 2>> \$GD
        echo "Number of paired reads after PCR duplicate removal and QC:" >> ${prefix}.QC.log
        echo \$(( \$(cat ${prefix}_QC_end_1.fastq | wc -l) / 4 ))         >> ${prefix}.QC.log
    elif [ "${force}" = "true" ]; then
        # --- Force dedup on full-length reads ---
        prinseq-lite.pl -derep 1 -fastq ${prefix}_noadapt_R1.fastq -fastq2 ${prefix}_noadapt_R2.fastq -out_format 3 -out_bad null -out_good ${prefix}_dwb 2>> \$GD
        prinseq-lite.pl -trim_left ${n1} -fastq ${prefix}_dwb_1.fastq -out_format 3 -out_bad null -out_good ${prefix}_dbr_1 2>> \$GD
        prinseq-lite.pl -trim_left ${n2} -fastq ${prefix}_dwb_2.fastq -out_format 3 -out_bad null -out_good ${prefix}_dbr_2 2>> \$GD
        prinseq-lite.pl -min_len 15 -fastq ${prefix}_dbr_1.fastq -fastq2 ${prefix}_dbr_2.fastq -out_format 3 -out_bad null -out_good ${prefix}_QC_end 2>> \$GD
        echo "Number of paired reads after PCR duplicate removal and QC:" >> ${prefix}.QC.log
        echo \$(( \$(cat ${prefix}_QC_end_1.fastq | wc -l) / 4 ))         >> ${prefix}.QC.log
    else
        # --- No dedup: strip additional barcode only ---
        prinseq-lite.pl -trim_left ${params.add_b1} -fastq ${prefix}_noadapt_R1.fastq -out_format 3 -out_bad null -out_good ${prefix}_br_1 2>> \$GD
        prinseq-lite.pl -trim_left ${params.add_b2} -fastq ${prefix}_noadapt_R2.fastq -out_format 3 -out_bad null -out_good ${prefix}_br_2 2>> \$GD
        prinseq-lite.pl -min_len 15 -fastq ${prefix}_br_1.fastq -fastq2 ${prefix}_br_2.fastq -out_format 3 -out_bad null -out_good ${prefix}_QC_end 2>> \$GD
        echo "Number of paired reads after final QC:"      >> ${prefix}.QC.log
        echo \$(( \$(cat ${prefix}_QC_end_1.fastq | wc -l) / 4 )) >> ${prefix}.QC.log
    fi

    # --- Optional dataset-wide length cutoff for mapping ---
    if [ "${params.map_length}" -ne 0 ]; then
        seqtk seq -L ${params.map_length} ${prefix}_QC_end_1.fastq > ${prefix}_QC_end_1.trimmed.fastq
        seqtk seq -L ${params.map_length} ${prefix}_QC_end_2.fastq > ${prefix}_QC_end_2.trimmed.fastq
        mv ${prefix}_QC_end_1.trimmed.fastq ${prefix}_QC_end_1.fastq
        mv ${prefix}_QC_end_2.trimmed.fastq ${prefix}_QC_end_2.fastq
    fi

    gzip -nf ${prefix}_QC_end_1.fastq ${prefix}_QC_end_2.fastq

    # --- machine-readable metrics for MultiQC ---
    IN=\$(( \$(zcat ${read1} | wc -l) / 4 ))
    PQ=\$(( \$(zcat ${prefix}_QC_end_1.fastq.gz | wc -l) / 4 ))
    printf 'sample\\t%s\\nseq_type\\tPE\\ninput\\t%s\\npass_qc\\t%s\\n' "${prefix}" "\$IN" "\$PQ" > ${prefix}.reads.metrics
    """

    stub:
    def prefix = meta.id
    """
    echo '' | gzip > ${prefix}_QC_end_1.fastq.gz
    echo '' | gzip > ${prefix}_QC_end_2.fastq.gz
    touch ${prefix}.QC.log ${prefix}_R1.cutadapt.log ${prefix}_R2.cutadapt.log ${prefix}.prinseq-pcrDups.gd
    printf 'sample\\t%s\\nseq_type\\tPE\\ninput\\t0\\npass_qc\\t0\\n' "${prefix}" > ${prefix}.reads.metrics
    """
}
