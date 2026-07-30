#!/usr/bin/env nextflow
/*
 * proseq2.0-nf : Nextflow DSL2 port of the Danko-Lab proseq2.0 pipeline.
 * Preprocess & align Run-On sequencing (PRO/GRO/ChRO-seq) data and emit bigWigs.
 */
nextflow.enable.dsl = 2

include { PREPROCESS_SE } from './modules/local/preprocess_se'
include { PREPROCESS_PE } from './modules/local/preprocess_pe'
include { CAT_FASTQ_SE  } from './modules/local/cat_fastq_se'
include { CAT_FASTQ_PE  } from './modules/local/cat_fastq_pe'
include { SORTMERNA_SE  } from './modules/local/sortmerna_se'
include { SORTMERNA_PE  } from './modules/local/sortmerna_pe'
include { BWA_ALIGN_SE  } from './modules/local/bwa_align_se'
include { BWA_ALIGN_PE  } from './modules/local/bwa_align_pe'
include { MAKE_BIGWIG   } from './modules/local/make_bigwig'
include { DETECT_ADAPTER_SE } from './modules/local/detect_adapter_se'
include { DETECT_ADAPTER_PE } from './modules/local/detect_adapter_pe'
include { FASTQC as FASTQC_RAW  } from './modules/local/fastqc'
include { FASTQC as FASTQC_TRIM } from './modules/local/fastqc'
include { COLLECT_STATS } from './modules/local/collect_stats'
include { INFER_STRAND  } from './modules/local/infer_strand'
include { MULTIQC       } from './modules/local/multiqc'

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------

// Resolve the library geometry (RNA-end centric) from --assay + low-level flags.
//   read_start : which RNA end sits at the read's (SE) / R1's (PE) 5' end.
//     rna_5prime = 5' methods (GRO-seq, PRO-cap, GRO-cap) -> read is SENSE to the gene
//     rna_3prime = 3' methods (PRO-seq, ChRO-seq)         -> read is ANTISENSE
//   report (PE): which RNA end to record (rna_5prime | rna_3prime | both); default = read_start.
//   antisense  : report on the opposite strand.
def resolveGeometry() {
    def assay = params.assay ? params.assay.toString().toUpperCase() : null
    def ASSAY_READ_START = [ GRO:'rna_5prime', PROCAP:'rna_5prime', GROCAP:'rna_5prime',
                             PRO:'rna_3prime', CHRO:'rna_3prime' ]
    if (assay && !ASSAY_READ_START.containsKey(assay))
        error "--assay must be one of GRO, PRO, ChRO, PROcap, GROcap (got '${params.assay}')."

    def read_start = params.read_start ?: (assay ? ASSAY_READ_START[assay] : 'rna_3prime')
    if (!(read_start in ['rna_5prime', 'rna_3prime']))
        error "--read_start must be rna_5prime or rna_3prime (got '${read_start}')."

    def report = params.report ?: read_start   // default: report the captured end
    if (!(report in ['rna_5prime', 'rna_3prime', 'both']))
        error "--report must be rna_5prime, rna_3prime or both (got '${report}')."

    return [ read_start: read_start, report: report, antisense: (params.antisense as boolean) ]
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
      (rows sharing a sample name are concatenated as technical replicates)

    Required:
      --input        Samplesheet CSV
      --bwa_index    Prefix of the BWA index (no .bwt suffix)
      --chrom_info   chromInfo table (chrom <TAB> size)

    Key options (see nextflow.config for all + defaults):
      Assay: --assay {GRO|PRO|ChRO|PROcap|GROcap}   preset for the geometry below
      Geom:  --read_start {rna_3prime|rna_5prime}   which RNA end is at the read/R1 5' end
                 rna_3prime = 3' methods (PRO/ChRO) | rna_5prime = 5' (GRO/*cap)  [default rna_3prime]
             --report {rna_5prime|rna_3prime|both}  PE only; which RNA end to record  [default = read_start]
             --antisense                            report on the opposite strand
             (--read_start/--report/--antisense override --assay)
      UMI: --umi1 N --umi2 N --add_b1 N --add_b2 N --force_deduplicate {true|false}
      Map: --aligner {aln|mem}  --dreg  --map_length N
      rRNA:--remove_rrna --rrna_refs FILE[,FILE...]   SortMeRNA pre-alignment rRNA depletion
      QC:  --gene_bed FILE   BED12 gene model -> RSeQC strand inference (validate & warn)
           --skip_adapter_detect   (fastp adapter detection, validate-and-warn; on by default)
           --skip_fastqc  --skip_multiqc
      Res: --max_cpus N   --max_memory '16.GB'   --max_time '24.h'
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
    if (params.remove_rrna && !params.rrna_refs)
        error "--remove_rrna requires --rrna_refs (rRNA FASTA file(s): comma-separated list or a glob)."

    // --- reference channels ---
    ch_index      = Channel.fromPath("${params.bwa_index}*", checkIfExists: true).collect()
    def bwa_prefix = file(params.bwa_index).name
    ch_chrom_info = file(params.chrom_info, checkIfExists: true)

    def geom = resolveGeometry()
    // Derived internals:
    //   primary read (SE read / R1) is SENSE to genes iff its 5' end is the RNA 5' end
    def expected_sense = (geom.read_start == 'rna_5prime')
    //   PE: which mate's 5' end is the RNA 5' end
    def rna5_mate = expected_sense ? 'R1_5prime' : 'R2_5prime'
    //   SE: report read 5' end same strand (G) if read sense, flipped (P) if antisense;
    //       --antisense toggles it
    def se_output = ((expected_sense) != geom.antisense) ? 'G' : 'P'
    log.info "Geometry: read_start=${geom.read_start} report=${geom.report} antisense=${geom.antisense}" +
             (params.assay ? " (from --assay ${params.assay})" : "")

    // --- read samplesheet -> one [meta, reads] per row ---
    ch_rows = Channel.fromPath(params.input, checkIfExists: true)
        | splitCsv(header: true)
        | map { row ->
            if (!row.sample)   error "Samplesheet row missing 'sample' column: ${row}"
            if (!row.fastq_1)  error "Samplesheet row '${row.sample}' missing fastq_1."
            def meta = [ id: row.sample.trim(), expected_sense: expected_sense ]
            def r1 = file(row.fastq_1.trim(), checkIfExists: true)
            if (row.fastq_2?.trim()) {
                meta.single_end = false
                meta.rna5   = rna5_mate
                meta.opp    = geom.antisense
                meta.report = geom.report
                return [ meta, [ r1, file(row.fastq_2.trim(), checkIfExists: true) ] ]
            } else {
                meta.single_end = true
                meta.se_output  = se_output
                return [ meta, [ r1 ] ]
            }
        }

    // --- group technical replicates: rows sharing `sample` are concatenated ---
    ch_grouped = ch_rows
        .map { meta, reads -> [ meta.id, meta, reads ] }
        .groupTuple()
        .map { id, metas, rlist ->
            if (metas.collect { it.single_end }.unique().size() > 1)
                error "Sample '${id}' mixes single-end and paired-end rows in the samplesheet."
            [ metas[0], rlist ]   // rlist = list of per-row read-lists
        }

    ch_grouped
        .branch { meta, rlist -> multi: rlist.size() > 1
                                 single: true }
        .set { ch_rep }

    // single-row samples pass straight through (no concat)
    ch_single = ch_rep.single.map { meta, rlist -> [ meta, rlist[0] ] }

    // multi-row samples: concatenate R1s (and R2s) in order
    ch_rep.multi
        .branch { meta, rlist -> se: meta.single_end
                                 pe: !meta.single_end }
        .set { ch_multi }
    CAT_FASTQ_SE( ch_multi.se.map { meta, rlist -> [ meta, rlist.collect { it[0] }.sort() ] } )
    CAT_FASTQ_PE( ch_multi.pe.map { meta, rlist -> [ meta, rlist.collect { it[0] }.sort(), rlist.collect { it[1] }.sort() ] } )

    ch_reads = ch_single
        .mix( CAT_FASTQ_SE.out.reads.map { meta, r      -> [ meta, [ r ] ] } )
        .mix( CAT_FASTQ_PE.out.reads.map { meta, r1, r2 -> [ meta, [ r1, r2 ] ] } )

    ch_reads
        .branch { meta, reads ->
            se: meta.single_end
            pe: !meta.single_end
        }
        .set { ch_branched }

    // --- preprocess ---
    ch_se_reads = ch_branched.se.map { meta, reads -> [ meta, reads[0] ] }
    PREPROCESS_SE( ch_se_reads )
    ch_pe_reads = ch_branched.pe.map { meta, reads -> [ meta, reads[0], reads[1] ] }
    PREPROCESS_PE( ch_pe_reads )

    // --- adapter detection on raw reads (validate-and-warn; never changes trimming) ---
    if (!params.skip_adapter_detect) {
        DETECT_ADAPTER_SE( ch_se_reads )
        DETECT_ADAPTER_PE( ch_pe_reads )
    }

    // --- optional rRNA depletion (SortMeRNA), else pass QC'd reads straight through ---
    def ch_se_to_align = PREPROCESS_SE.out.reads
    def ch_pe_to_align = PREPROCESS_PE.out.reads
    if (params.remove_rrna) {
        ch_rrna = Channel.fromPath(params.rrna_refs.tokenize(','), checkIfExists: true).collect()
        SORTMERNA_SE( PREPROCESS_SE.out.reads, ch_rrna )
        SORTMERNA_PE( PREPROCESS_PE.out.reads, ch_rrna )
        ch_se_to_align = SORTMERNA_SE.out.reads
        ch_pe_to_align = SORTMERNA_PE.out.reads
    }

    // --- align ---
    BWA_ALIGN_SE( ch_se_to_align, ch_index, bwa_prefix )
    BWA_ALIGN_PE( ch_pe_to_align, ch_index, bwa_prefix )

    // --- bigWigs (shared) ---
    ch_bams = BWA_ALIGN_SE.out.bam.mix( BWA_ALIGN_PE.out.bam )
    MAKE_BIGWIG( ch_bams, ch_chrom_info )

    // --- QC: FastQC (raw + trimmed) and MultiQC ---
    ch_multiqc = Channel.empty()

    if (!params.skip_fastqc) {
        // raw reads (post-merge): ch_reads carries [meta, [reads]] for SE and PE
        FASTQC_RAW( ch_reads, 'raw' )

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
    if (!params.skip_adapter_detect) {
        ch_multiqc = ch_multiqc.mix( DETECT_ADAPTER_SE.out.json, DETECT_ADAPTER_PE.out.json )
    }
    if (params.remove_rrna) {
        ch_multiqc = ch_multiqc.mix( SORTMERNA_SE.out.log, SORTMERNA_PE.out.log )
    }

    // --- strand inference / validation (opt-in via --gene_bed) ---
    if (params.gene_bed) {
        ch_gene_bed = file(params.gene_bed, checkIfExists: true)
        ch_align    = BWA_ALIGN_SE.out.bam.mix( BWA_ALIGN_PE.out.bam )
        INFER_STRAND( ch_align, ch_gene_bed )
        ch_multiqc = ch_multiqc.mix( INFER_STRAND.out.report )
    }

    // read-count / dedup / mapping table as MultiQC custom content
    ch_metrics = PREPROCESS_SE.out.metrics
        .mix( PREPROCESS_PE.out.metrics, MAKE_BIGWIG.out.metrics )
    if (params.remove_rrna) {
        ch_metrics = ch_metrics.mix( SORTMERNA_SE.out.metrics, SORTMERNA_PE.out.metrics )
    }
    COLLECT_STATS( ch_metrics.collect() )
    ch_multiqc = ch_multiqc.mix( COLLECT_STATS.out.mqc )

    if (!params.skip_multiqc) {
        MULTIQC( ch_multiqc.collect() )
    }
}
