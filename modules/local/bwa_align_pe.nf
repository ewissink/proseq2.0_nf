process BWA_ALIGN_PE {
    tag   "$meta.id"
    label 'process_high'
    conda "${projectDir}/env/environment.yml"

    input:
    tuple val(meta), path(read1), path(read2)
    path  index          // all bwa index files (share a common prefix)
    val   bwa_prefix     // basename of the index prefix

    output:
    tuple val(meta), path("${meta.id}.sort.bam"), emit: bam

    script:
    def prefix  = meta.id
    def use_aln = ((params.aligner ?: 'mem') == 'aln') ? 'true' : 'false'
    """
    set -eu

    if [ "${use_aln}" = "true" ]; then
        # BWA-backtrack
        bwa aln -t ${task.cpus} ${bwa_prefix} ${read1} > ${prefix}_sa1.sai
        bwa aln -t ${task.cpus} ${bwa_prefix} ${read2} > ${prefix}_sa2.sai
        bwa sampe -n 1 -f ${prefix}_end.sam ${bwa_prefix} ${prefix}_sa1.sai ${prefix}_sa2.sai ${read1} ${read2}
        samtools view -bf 0x2 -q 20 ${prefix}_end.sam | samtools sort -n -@ ${task.cpus} - > ${prefix}.sort.bam
        rm -f ${prefix}_sa1.sai ${prefix}_sa2.sai ${prefix}_end.sam
    else
        # default: BWA-MEM (properly-paired, q20)
        bwa mem -k 19 -t ${task.cpus} ${bwa_prefix} ${read1} ${read2} | samtools view -bf 0x2 -q 20 - | samtools sort -n -@ ${task.cpus} - > ${prefix}.sort.bam
    fi
    """

    stub:
    """
    touch ${meta.id}.sort.bam
    """
}
