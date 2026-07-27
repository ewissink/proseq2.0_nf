process BWA_ALIGN_SE {
    tag   "$meta.id"
    label 'process_high'
    conda "${projectDir}/env/environment.yml"

    input:
    tuple val(meta), path(reads)
    path  index          // all bwa index files (share a common prefix)
    val   bwa_prefix     // basename of the index prefix

    output:
    tuple val(meta), path("${meta.id}.sort.bam"), emit: bam

    script:
    def prefix  = meta.id
    def use_mem = (params.aligner ?: 'aln') == 'mem'
    def dreg    = params.dreg ? 'true' : 'false'
    def mem     = use_mem ? 'true' : 'false'
    """
    set -eu

    if [ "${dreg}" = "true" ]; then
        # dREG-compatible: bwa aln, keep all reads (q0)
        bwa aln -t ${task.cpus} ${bwa_prefix} ${reads} | bwa samse -n 1 -f ${prefix}_end.sam ${bwa_prefix} - ${reads}
        samtools view -b -q 0 ${prefix}_end.sam | samtools sort -n -@ ${task.cpus} - > ${prefix}.sort.bam
        rm -f ${prefix}_end.sam
    elif [ "${mem}" = "true" ]; then
        bwa mem -k 19 -t ${task.cpus} ${bwa_prefix} ${reads} | samtools view -b -q 20 - | samtools sort -n -@ ${task.cpus} - > ${prefix}.sort.bam
    else
        # default: bwa aln (backtrack), q20
        bwa aln -t ${task.cpus} ${bwa_prefix} ${reads} | bwa samse -n 1 -f ${prefix}_end.sam ${bwa_prefix} - ${reads}
        samtools view -b -q 20 ${prefix}_end.sam | samtools sort -n -@ ${task.cpus} - > ${prefix}.sort.bam
        rm -f ${prefix}_end.sam
    fi
    """

    stub:
    """
    touch ${meta.id}.sort.bam
    """
}
