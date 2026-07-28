# Concordance tests

Validate the Nextflow port by running the **same subsampled reads** through both
`proseq2.0.bsh` (original) and the Nextflow port with matched parameters, then
diffing the bigWigs. Identical bigWigs ⇒ the port faithfully reproduces the
original.

This exercises the strand/end matrix (GRO 5′, PRO, PE report-5′, PE report-3′)
because those flags only change how reads are *reported* — the same reads test
every code path.

## 1. Environment

```bash
conda env create -f test/environment.yml
conda activate proseq2.0-nf-test
# nextflow must also be on PATH
```

The test env pins the **same tool versions** used by the pipeline, so both the
bash script and the NF port trim/align identically (a fair comparison). The runs
are forced single-threaded to remove bwa's multi-thread nondeterminism.

## 2. Genome (mouse or human)

If you already have a BWA index + chromInfo, use it. Otherwise build one (heavy,
one-time — mammalian `bwa index` ≈ 1 h, several GB):

```bash
test/setup_genome.sh mm10        # or: hg38
```

It prints the `--bwa-index` / `--chrom-info` paths to pass along. For *real*
analysis (not just concordance) add the rDNA sequence (GenBank U13369.1) to the
FASTA before indexing — see the note in `setup_genome.sh`.

## 3. Datasets

Edit `datasets.tsv` and fill in real SRA run accessions (the shipped ones are
placeholders). Suggested minimal matrix — pick mouse **or** human, stay consistent
with your index:

| name           | assay | layout | report | why |
|----------------|-------|--------|--------|-----|
| `gro_se`       | GRO   | SE     | 5prime | GRO-seq 5′ path (`-G`) |
| `pro_se`       | PRO   | SE     | 5prime | PRO-seq path (`-P`) |
| `prochro_pe_5p`| ChRO  | PE     | 5prime | PE report-5′ (`--map5 true`) |
| `prochro_pe_3p`| ChRO  | PE     | 3prime | PE report-3′ (`--map5 false`) — reuse the same PE accession |

Find accessions via GEO → "SRA Run Selector" (search the Danko lab or the assay +
organism). For PE, use libraries with **distinct i7 barcodes** (proseq2.0 requires
that).

## 4. Run

```bash
test/run_concordance.sh \
  --bwa-index  test/ref/mm10.fa.gz \
  --chrom-info test/ref/mm10.chromInfo \
  --reads 1000000 --seed 100
```

Per dataset it fetches + subsamples (fixed seed; PE mates subsampled with the same
seed so pairs stay matched), runs both pipelines into `test/out/<name>/{bash,nf}`,
and prints `PASS (identical)` or a `DIFF` summary per strand. Exit code is non-zero
if anything differs.

## Reading the results

- **PASS (identical)** on all strands → the port matches the original. 
- **DIFF** → it reports interval counts and how many are identical, and leaves the
  two `.bg` files so you can `diff` them. A handful of differing intervals usually
  means tool-version drift (e.g. cutadapt); large or systematic diffs mean a real
  logic difference worth chasing.
- The raw `_plus.bw` / `_minus.bw` are compared (RPM tracks scale from the same
  counts, so raw is the cleanest signal).

## Biological spot-check (optional)

Add `--gene_bed <BED12>` to a normal pipeline run to have the strand-inference step
confirm each library's orientation matches its assay flags (a refGene BED12 for
your genome works). See the main README's *Strand inference* section.
