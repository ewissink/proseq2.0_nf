# proseq2.0-nf — Nextflow port

A [Nextflow](https://www.nextflow.io/) DSL2 port of the Danko-Lab **proseq2.0**
pipeline for preprocessing and aligning Run-On sequencing (PRO-seq / GRO-seq /
ChRO-seq) data. It reproduces the logic of the original `proseq2.0.bsh` script,
decomposed into per-sample, resumable, parallelizable processes. Comparison 
showing the two pipelines have identical outputs are available [here](https://github.com/ewissink/proseq2.0_nf/tree/main/test).

The original bash script (`proseq2.0.bsh`) is untouched and still works; its
documentation is preserved in [`README_proseq2.0.md`](README_proseq2.0.md). This
Nextflow port lives alongside it.

## What it does

Per sample, the pipeline runs three stages (SE and PE variants):

1. **Preprocess** — remove 3′ adapters and quality-trim (`cutadapt`), strip
   UMI / additional barcodes, and (when UMIs are present or `--force_deduplicate`)
   remove PCR duplicates (`prinseq-lite.pl`, `seqtk`).
2. **Align** — map to the reference with `bwa` (`aln`/backtrack or `mem`),
   filter and name-sort with `samtools`.
3. **BigWig** — take the reported RNA end (5′ or 3′), drop rRNA/chrM/non-primary
   contigs, and write raw + RPM-normalized plus/minus `bigWig` files
   (`bedtools`, `bedops`, `bedGraphToBigWig`).

Plus optional QC (on by default):

- **FastQC** on raw input reads and again on the trimmed/QC'd reads.
- **MultiQC** aggregating FastQC + cutadapt trimming stats, plus a **proseq2.0
  read-summary table** (input → pass-QC → mappable → mappable-excl-rRNA per
  sample), into one `multiqc_report.html`.
- **Strand inference** (opt-in, `--gene_bed`) — RSeQC `infer_experiment.py` on a
  BAM subsample checks whether the data's orientation matches the configured
  strand flags and **warns** on a mismatch (see *Strand inference* below).
- **rRNA depletion** (opt-in, `--remove_rrna`) — SortMeRNA removes rRNA reads
  *before* alignment against reference rRNA FASTA(s), giving a clean `% rRNA`
  metric without needing rDNA in the genome index (see *rRNA removal* below).

## Requirements

- Nextflow ≥ 22.10 (tested on 26.04)
- One of: **conda/mamba**, **docker**, or **singularity** — or the tools on
  `$PATH` (see `env/environment.yml`): `cutadapt`, `seqtk`, `prinseq-lite.pl`,
  `bwa`, `samtools`, `bedtools`, `sort-bed` (bedops), `bedGraphToBigWig`.

## Inputs

**Samplesheet** (`--input`, CSV with a header). Leave `fastq_2` empty for
single-end samples. For technical replicates, include two rows with the
same sample name, and the fastq files will be concatenated in the pipeline.

```csv
sample,fastq_1,fastq_2
sampleA_PE,/data/sampleA_R1.fastq.gz,/data/sampleA_R2.fastq.gz
sampleB_SE,/data/sampleB.fastq.gz,
```

**Reference**:
- `--bwa_index` — the **prefix** of a `bwa index` (no `.bwt` suffix). All
  `PREFIX.*` files are staged automatically.
- `--chrom_info` — a 2-column `chrom<TAB>size` table.

The reference genome should not include any alternative haplotypes. To create your own BWA index and 
 chrom_info file for `mm10`, `hg38`, or `dm6`,use the the following
with the correct genome:
```bash
bin/setup_genome.sh <genome>
```

The original pipeline aligned to a reference genome that included the rDNA transcript, then
removed rRNA-matching reads after-the-fact. That is still possible here; however, rRNA can
instead be removed using SortMeRNA prior to alignment.

## Minimal usage

```bash
nextflow run main.nf -profile <conda/docker/singularity> \
    --input <samplesheet-path> \
    --bwa_index <bwa-path> \
    --chrom_info <chrom_info-path> \
    --outdir results \
    --assay <PRO/GRO/ChRO|PROcap|GROcap> 

```

`nextflow run main.nf --help` prints the option summary. Add `-resume` to reuse
cached results after an interruption.

## Assay presets

`--assay {GRO|PRO|ChRO|PROcap|GROcap}` sets the library geometry, which then sets presets for 
other pipeline flags. These are based on published methods for creating these types of libraries:
- GRO-seq: captures the active site of Pol II at the 5' end of the RNA (see [Core et al, 2008](https://www.science.org/doi/10.1126/science.1162228))
- GRO-cap: captures the TSS of the RNA at the 5' end (see [Core et al, 2014](https://www.nature.com/articles/ng.3142))
- PRO-/ and ChRO-seq: captures the active site of Pol II at the 3' end of the RNA (see [Kwak et al, 2013](https://www.science.org/doi/10.1126/science.1229386) and [Chu et al, 2018](https://www.nature.com/articles/s41588-018-0244-3))
- PRO-cap: captures the TSS of the RNA at the 5' end (see [Kwak et al, 2013](https://www.science.org/doi/10.1126/science.1229386))

| `--assay` | `--read_start` | `--report` | `--antisense` |
|-----------|------------------|-----------|--------- |
| `GRO`     | `rna_5prime` | `rna_5prime` | false |
| `GROcap`     | `rna_5prime` | `rna_5prime` |  false |
| `PRO`     | `rna_3prime` | `rna_3prime` | false |
| `PROcap`     | `rna_5prime` |`rna_5prime` |  false |
| `ChRO`    | `rna_3prime`  | `rna_3prime` | false |

Any explicit flag (`--re_start`, `--report`, `--antisense`) 
**overrides** the preset. 


## Options

| Parameter   |  Purpose         | Default            | Options |
|------------------|------------------|--------------------|---------|
| `--assay` |  set geometry for read reporting (see above for details) | none | `{GRO\|PRO\|ChRO\|GROcap\|PROcap`} | 
| `--read_start` | RNA end sits at the read's (SE) / R1's (PE) 5' end, based on assay. Overrides `--assay` preset | `{rna_3prime\|rna_5prime}` |
| `--report` | paired-end only; sets the RNA end to record in bigWigs; both is useful for PRO-cap where TSS and active site can be captured | `--read_start` parameter | `{rna_3prime\|rna_5prime\|both}`} |
| `--antisense` | report on the opposite (antisense) strand | `false` | `{false\|true}` |
| `--adapter_se` | Adapter sequence trimmed in single-end sequencing | `'TGGAATTCTCGGGTGCCAAGG'` | string of adapter sequence |
| `--adapter1` | Adapter sequence trimmed from R2 | `'GATCGTCGGACTGTAGAACTCTGAACG'` | `string` of adapter sequence |
| `--adapter2` | Adapter sequence trimmed from R1 | `'AGATCGGAAGAGCACACGTCTGAACTC'` | `string` of adapter sequence |
| `--umi1` | length of UMI on R1 | 0 | `integer` of UMI length |
| `--umi2` | length of UMI on R2 | 0 | `integer` of UMI length |
| `--add_b1` | extra barcode trimmed from R1 | 0 | `integer` of length to trim |
| `--add_b2` | extra barcode trimmed from R2 | 0 | `integer` of length to trim |
| `--force_deduplicate` | deduplicate in absence of UMIs | `false` | `{false\|true}` |
| `--aligner` |choice between BWA-backtrack (`aln`) or BWA-mem (`mem`) | SE->`aln`, PE->`mem`| `{aln\|mem}` |
| `--dreg` | use dREG-compatible parameters, only for SE | `false` | `{false\|true}` |
| `--map_length` | align whole read (0) or set a data-set wide length cutoff for mapping (other integer) | 0 | `integer` of read-length to align |
| `--skip_fastqc` | Don't run FastQC on raw and trimmed reads | `false` | `{false\|true}` |
| `--skip_multiqc` | Don't run MultiQC at end of run | `false` | `{false\|true}` |
| `--gene_bed` | when set, run RSeQC strand inference and warn if it doesn't match configured strand | `null` | path to BED12 gene model for your species |
| `--remove_rrna` | Run pre-alignment rRNA depletion | `false` | `{false\|true}` |
| `--rrna_refs` | Required when running `--remove_rna` | `null` | path to rRNA reference FASTA(s) |


Resource ceilings: `--max_cpus`, `--max_memory`, `--max_time`. On a laptop,
lower them (e.g. `--max_memory '8.GB' --max_cpus 4`).

## Outputs (`--outdir`)

```
results/
├── bam/       <sample>.sort.bam            (name-sorted)
├── bigwig/    <sample>_{plus,minus}.bw     (raw counts; minus strand negated)
│              <sample>_{plus,minus}.rpm.bw (RPM-normalized)
├── qc/        <sample>.QC.log, <sample>.align.log, <sample>.prinseq-pcrDups.gd
│              <sample>.cutadapt.log, proseq_read_stats_mqc.tsv
│   ├── fastqc/{raw,trim}/  per-sample FastQC reports
│   ├── strand/            <sample>.infer_experiment.txt, <sample>.strand_check.txt
│   └── sortmerna/         <sample>.sortmerna.log   (when --remove_rrna)
├── multiqc/   multiqc_report.html + multiqc_data/
└── pipeline_info/  execution report, timeline, trace, DAG
```



## Strand inference

The read orientation is fixed by the geometry above: the primary read (SE read,
or R1 in PE) is *sense* to genes iff the RNA 5′ end sits at that read's 5′ end.
Which end to *report* (`--map5`) and `--opposite_strand` are experiment choices
and can't be inferred from data.

Pass `--gene_bed <BED12>` to have RSeQC `infer_experiment.py` (run on a BAM
subsample) check the data against your configured geometry. Per sample you get
`qc/strand/<sample>.strand_check.txt` with a `VERDICT: PASS|WARN|UNDETERMINED`,
and the raw RSeQC output flows into the MultiQC report. A `WARN` means the data's
orientation disagrees with your flags (e.g. swapped R1/R2, or wrong assay) — the
pipeline still completes; it does not auto-change your settings.

## rRNA removal

By default the pipeline removes rRNA the way proseq2.0 does: *after* alignment, by
dropping reads on contigs named `rRNA`/`chrM`. Many assemblies (e.g. UCSC dm6/mm10/hg38)
do **not** include the nuclear rDNA arrays, so rRNA reads simply fail to map and
vanish as "unmapped".

For a portable, standard alternative, `--remove_rrna` runs **SortMeRNA** *before*
alignment to filter reads by sequence against reference rRNA FASTA(s):

```bash
nextflow run main.nf ... --remove_rrna --rrna_refs <rrna-db-path>
```

- `--rrna_refs` takes one or more FASTA files (comma-separated, or a glob). Use an
  rRNA database (e.g. the SortMeRNA/SILVA sets) or your organism's rRNA sequences. You
  can download a ready-made database that has rRNA sequences for eukaryotes and prokaryotes
  using:
  ```bash
  wget https://github.com/sortmerna/sortmerna/releases/download/v7.0.0/smr_v4.3_default_db.fasta.gz
  gunzip smr_v4.3_default_db.fasta.gz
  ```
- Paired-end uses `--paired_in` (a pair is dropped if *either* mate is rRNA).
- Adds a **Post rRNA-removal** column to the read-summary table and feeds the
  SortMeRNA log into MultiQC (`% rRNA`).

This runs between preprocessing and alignment; the post-alignment `chrM`/`rRNA`
coordinate filter still applies as a backstop.

## Profiles

`-profile conda | mamba | docker | singularity | slurm` (combine with commas,
e.g. `-profile singularity,slurm`). `conda`/`mamba` build the env from
`env/environment.yml`. For docker/singularity, set a container in a custom
config or `process.container`.

## Testing without data

Every process has a `stub:` block, so the channel topology can be verified with
no tools or reference data:

```bash
nextflow run main.nf -stub-run --input samplesheet.csv \
    --bwa_index dummy --chrom_info dummy.txt
```

## Not yet ported / differences from `proseq2.0.bsh`

- **Pre-mapping** (`--PREMAP_BWAIDX`, an undocumented flag that maps to a decoy
  index and keeps unmapped reads) is not implemented. Open an issue if you need it.
- Paired-end **singleton** reads from dedup are not carried forward (the original
  wrote them to disk but never used them downstream — behavior is unchanged).
- Read counts in `*.QC.log` are computed as `lines/4` rather than the original's
  `grep @ -c` (more accurate; the numbers may differ slightly).
- The original wrote all samples into one shared temp dir and one combined log;
  here each sample is isolated in its own Nextflow work dir with per-sample logs.

Alignment parameters, adapter/UMI handling, strand logic (all 8 PE cases + the
2 SE cases), rRNA/chrM filtering, and RPM normalization mirror the original.

## Cite

Chu T, Wang Z, Chou SP, Danko CG. *Discovering Transcriptional Regulatory 
Elements From Run-On and Sequencing Data Using the Web-Based dREG Gateway.* 
Curr Protoc Bioinformatics. 2019 Jun;66(1):e70. doi: 10.1002/cpbi.70. 

Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. 
*Nextflow enables reproducible computational workflows.* 
Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. 
