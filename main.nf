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
include { FASTQC as FASTQC_RAW  } from './modules/local/fastqc'
include { FASTQC as FASTQC_TRIM } from './modules/local/fastqc'
include { COLLECT_STATS } from './modules/local/collect_stats'
include { INFER_STRAND  } from './modules/local/infer_strand'
include { MULTIQC       } from './modules/local/multiqc'

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------

// Expand the --assay preset into the SE/PE geometry. Individual flags override it.
//   GRO  : capture RNA 5' end  -> SE -G ; PE RNA 5' at R1 (rna3=R2_5prime)
//   PRO  : capture RNA 3' end  -> SE -P ; PE RNA 3' at R1 (rna3=R1_5prime)
//   ChRO : same geometry as PRO-seq (PE run-on)
def effectiveGeometry() {
    def assay = params.assay ? params.assay.toString().toUpperCase() : null
    if (assay && !(assay in ['GRO', 'PRO', 'CHRO']))
        error "--assay must be one of GRO, PRO, ChRO (got '${params.assay}')."
    def assay_se_read = (assay == 'PRO' || assay == 'CHRO') ? 'RNA_3prime' : 'RNA_5prime'
    def assay_rna3    = (assay == 'PRO' || assay == 'CHRO') ? 'R1_5prime'  : 'R2_5prime'
    return [
        se_read: params.se_read ?: (assay ? assay_se_read : 'RNA_3prime'),   // SE default: PRO (-P)
        rna3   : params.rna3    ?: (assay ? assay_rna3    : 'R1_5prime'),   // PE default: PRO/ChRO
    ]
}

// Reproduce the SE -G/-P and PE --RNA5/--RNA3 derivation logic from proseq2.0.bsh
def deriveStrand(String rna3_in) {
    def rna5 = params.rna5
    def rna3 = rna3_in
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
      (rows sharing a sample name are concatenated as technical replicates)

    Required:
      --input        Samplesheet CSV
      --bwa_index    Prefix of the BWA index (no .bwt suffix)
      --chrom_info   chromInfo table (chrom <TAB> size)

    Key options (see nextflow.config for all + defaults):
      Assay: --assay {GRO|PRO|ChRO}   preset for the SE/PE geometry below
                                      (default, no assay: PRO/ChRO geometry)
      SE:  --se_read {RNA_3prime|RNA_5prime}   PRO-seq (-P) | GRO-seq (-G)   [default RNA_3prime]
      PE:  --rna5/--rna3 {R1_5prime|R2_5prime}   [default --rna3 R1_5prime]
           --map5 {true|false}   --opposite_strand {true|false}
           (individual flags override --assay)
      UMI: --umi1 N --umi2 N --add_b1 N --add_b2 N --force_deduplicate {true|false}
      Map: --aligner {aln|mem}  --dreg  --map_length N
      rRNA:--remove_rrna --rrna_refs FILE[,FILE...]   SortMeRNA pre-alignment rRNA depletion
      QC:  --gene_bed FILE   BED12 gene model -> RSeQC strand inference (validate & warn)
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

    def geom   = effectiveGeometry()
    def strand = deriveStrand(geom.rna3)
    log.info "Geometry: SE se_read=${geom.se_read} | PE rna5=${strand.rna5} rna3=${strand.rna3} " +
             "map5=${strand.map5} opposite_strand=${strand.opp}" +
             (params.assay ? " (from --assay ${params.assay})" : "")

    // --- read samplesheet -> one [meta, reads] per row ---
    ch_rows = Channel.fromPath(params.input, checkIfExists: true)
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
                // primary read (R1) is sense to genes iff RNA 5' end is at R1's 5' end
                meta.expected_sense = (strand.rna5 == 'R1_5prime')
                return [ meta, [ r1, file(row.fastq_2.trim(), checkIfExists: true) ] ]
            } else {
                meta.single_end = true
                meta.se_output  = (geom.se_read == 'RNA_5prime') ? 'G' : 'P'
                // SE read is sense to genes for GRO (-G), antisense for PRO (-P)
                meta.expected_sense = (meta.se_output == 'G')
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
