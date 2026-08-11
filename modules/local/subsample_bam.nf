process SUBSAMPLE_BAM {
    tag   "$meta.id"
    label 'process_low'
    conda "${projectDir}/env/environment.yml"

    input:
    tuple val(meta), path(bam), val(count), val(min_count)

    output:
    tuple val(meta), path("${meta.id}.sub.bam"), emit: bam

    script:
    def prefix = meta.id
    def seed   = params.subsample_seed
    """
    set -eu
    # fraction to keep; the shallowest sample (count <= min) is kept whole.
    frac=\$(awk -v a=${min_count} -v b=${count} 'BEGIN{ printf "%.6f", (b>0 ? a/b : 1) }')
    if awk -v f=\$frac 'BEGIN{ exit !(f>=1) }'; then
        cp ${bam} ${prefix}.sub.bam
    else
        # samtools -s takes SEED.FRACTION as one float (integer part = seed)
        S=\$(awk -v s=${seed} -v f=\$frac 'BEGIN{ printf "%.6f", s+f }')
        samtools view -b -s \$S ${bam} > ${prefix}.sub.bam
    fi
    """

    stub:
    """
    cp ${bam} ${meta.id}.sub.bam
    """
}
