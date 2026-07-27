process INFER_STRAND {
    tag   "$meta.id"
    label 'process_low'
    conda "${projectDir}/env/environment.yml"

    publishDir "${params.outdir}/qc/strand", mode: params.publish_mode

    input:
    tuple val(meta), path(bam)
    path  gene_bed

    output:
    path "${meta.id}.infer_experiment.txt", emit: report
    path "${meta.id}.strand_check.txt",     emit: check

    script:
    def prefix   = meta.id
    def expected = meta.expected_sense ? 'sense' : 'antisense'
    def seqtype  = meta.single_end ? 'SE' : 'PE'
    """
    set -eu
    # infer_experiment needs a coordinate-sorted, indexed BAM (ours is name-sorted)
    samtools sort -@ ${task.cpus} -o ${prefix}.cs.bam ${bam}
    samtools index ${prefix}.cs.bam

    infer_experiment.py -i ${prefix}.cs.bam -r ${gene_bed} -s 400000 -q 20 \
        > ${prefix}.infer_experiment.txt 2> ${prefix}.rseqc.err || true

    strand_check.py ${prefix}.infer_experiment.txt ${expected} ${prefix} ${seqtype}

    rm -f ${prefix}.cs.bam ${prefix}.cs.bam.bai
    """

    stub:
    def prefix = meta.id
    """
    printf 'This is stub Data\\nFraction of reads explained by "++,--": 1.0\\n' > ${prefix}.infer_experiment.txt
    echo 'VERDICT: PASS - stub' > ${prefix}.strand_check.txt
    """
}
