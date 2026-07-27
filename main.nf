#!/usr/bin/env nextflow
/*
 * proseq2.0-nf : Nextflow DSL2 port of the Danko-Lab proseq2.0 pipeline.
 * Preprocess & align Run-On sequencing (PRO/GRO/ChRO-seq) data and emit bigWigs.
 */
nextflow.enable.dsl = 2

include { PREPROCESS_SE } from './modules/local/preprocess_se'
include { PREPROCESS_PE } from './modules/local/preprocess_pe'
include { BWA_ALIGN_SE  } from './modules/local/bwa_align_se'
include { BWA_ALIGN_PE  } from './modules/local/bwa_align_pe'
include { MAKE_BIGWIG   } from './modules/local/make_bigwig'
include { FASTQC as FASTQC_RAW  } from './modules/local/fastqc'
include { FASTQC as FASTQC_TRIM } from './modules/local/fastqc'
include { COLLECT_STATS } from './modules/local/collect_stats'
include { MULTIQC       } from './modules/local/multiqc'

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------

// Reproduce the SE -G/-P and PE --RNA5/--RNA3 derivation logic from proseq2.0.bsh
def deriveStrand() {
    def rna5 = params.rna5
    def rna3 = params.rna3
    if (!rna5 && !rna3) rna5 = 'R1_5prime'
    if (rna3 == 'R1_5prime')      rna5 = 'R2_5prime'
    else if (rna3 == 'R2_5prime') rna5 = 'R1_5prime'
    if (rna5 != 'R1_5prime' && rna5 != 'R2_5prime')
        error "--rna5/--rna3 must resolve to R1_5prime or R2_5prime (got rna5=${rna5}, rna3=${rna3})"
    return [ rna5: rna5, rna3: rna3, map5: (params.map5 as boolean), opp: (params.opposite_strand as boolean) ]
}

def helpMessage() {
    log.info """
    proseq2.0-nf  —  Nextflow port of Danko-Lab proseq2.0

    Usage:
      nextflow run main.nf --input samplesheet.csv \\
          --bwa_index /path/to/index_prefix --chrom_info /path/to/chromInfo \\
          --outdir results -profile conda [options]

    Samplesheet (CSV, header required):
      sample,fastq_1,fastq_2
      sampleA,/data/A_R1.fastq.gz,/data/A_R2.fastq.gz   # paired-end
      sampleB,/data/B.fastq.gz,                          # single-end (empty fastq_2)

    Required:
      --input        Samplesheet CSV
      --bwa_index    Prefix of the BWA index (no .bwt suffix)
      --chrom_info   chromInfo table (chrom <TAB> size)

    Key options (see nextflow.config for all + defaults):
      SE:  --se_read {RNA_5prime|RNA_3prime}        (GRO-seq | PRO-seq)
      PE:  --rna5/--rna3 {R1_5prime|R2_5prime}  --map5 {true|false}  --opposite_strand {true|false}
      UMI: --umi1 N --umi2 N --add_b1 N --add_b2 N --force_deduplicate {true|false}
      Map: --aligner {aln|mem}  --dreg  --map_length N
    """.stripIndent()
}

// -----------------------------------------------------------------------------
// Workflow
// -----------------------------------------------------------------------------
workflow {
    if (params.help) { helpMessage(); return }

    // --- validate ---
    if (!params.input)      error "Missing --input (samplesheet CSV)."
    if (!params.bwa_index)  error "Missing --bwa_index (BWA index prefix)."
    if (!params.chrom_info) error "Missing --chrom_info (chromInfo table)."
    if (params.dreg && (params.aligner == 'mem'))
        error "--dreg is only compatible with bwa aln (do not combine with --aligner mem)."

    // --- reference channels ---
    ch_index      = Channel.fromPath("${params.bwa_index}*", checkIfExists: true).collect()
    def bwa_prefix = file(params.bwa_index).name
    ch_chrom_info = file(params.chrom_info, checkIfExists: true)

    def strand = deriveStrand()

    // --- read samplesheet -> [meta, reads] ---
    ch_input = Channel.fromPath(params.input, checkIfExists: true)
        | splitCsv(header: true)
        | map { row ->
            if (!row.sample)   error "Samplesheet row missing 'sample' column: ${row}"
            if (!row.fastq_1)  error "Samplesheet row '${row.sample}' missing fastq_1."
            def meta = [ id: row.sample.trim() ]
            def r1 = file(row.fastq_1.trim(), checkIfExists: true)
            if (row.fastq_2?.trim()) {
                meta.single_end = false
                meta.rna5 = strand.rna5; meta.rna3 = strand.rna3
                meta.map5 = strand.map5; meta.opp  = strand.opp
                return [ meta, [ r1, file(row.fastq_2.trim(), checkIfExists: true) ] ]
            } else {
                meta.single_end = true
                meta.se_output  = (params.se_read == 'RNA_5prime') ? 'G' : 'P'
                return [ meta, [ r1 ] ]
            }
        }

    ch_input
        .branch { meta, reads ->
            se: meta.single_end
            pe: !meta.single_end
        }
        .set { ch_branched }

    // --- single-end path ---
    ch_se_reads = ch_branched.se.map { meta, reads -> [ meta, reads[0] ] }
    PREPROCESS_SE( ch_se_reads )
    BWA_ALIGN_SE( PREPROCESS_SE.out.reads, ch_index, bwa_prefix )

    // --- paired-end path ---
    ch_pe_reads = ch_branched.pe.map { meta, reads -> [ meta, reads[0], reads[1] ] }
    PREPROCESS_PE( ch_pe_reads )
    BWA_ALIGN_PE(
        PREPROCESS_PE.out.reads.map { meta, r1, r2 -> [ meta, r1, r2 ] },
        ch_index, bwa_prefix
    )

    // --- bigWigs (shared) ---
    ch_bams = BWA_ALIGN_SE.out.bam.mix( BWA_ALIGN_PE.out.bam )
    MAKE_BIGWIG( ch_bams, ch_chrom_info )

    // --- QC: FastQC (raw + trimmed) and MultiQC ---
    ch_multiqc = Channel.empty()

    if (!params.skip_fastqc) {
        // raw reads: ch_input carries [meta, [reads]] for both SE and PE
        FASTQC_RAW( ch_input, 'raw' )

        // trimmed reads, normalized to [meta, [reads]]
        ch_trim = PREPROCESS_SE.out.reads.map { meta, r -> [ meta, [ r ] ] }
            .mix( PREPROCESS_PE.out.reads.map { meta, r1, r2 -> [ meta, [ r1, r2 ] ] } )
        FASTQC_TRIM( ch_trim, 'trim' )

        ch_multiqc = ch_multiqc.mix( FASTQC_RAW.out.zip, FASTQC_TRIM.out.zip )
    }

    ch_multiqc = ch_multiqc.mix(
        PREPROCESS_SE.out.cutadapt,
        PREPROCESS_PE.out.cutadapt
    )

    // read-count / dedup / mapping table as MultiQC custom content
    ch_metrics = PREPROCESS_SE.out.metrics
        .mix( PREPROCESS_PE.out.metrics, MAKE_BIGWIG.out.metrics )
    COLLECT_STATS( ch_metrics.collect() )
    ch_multiqc = ch_multiqc.mix( COLLECT_STATS.out.mqc )

    if (!params.skip_multiqc) {
        MULTIQC( ch_multiqc.collect() )
    }
}
