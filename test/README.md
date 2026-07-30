# Concordance tests

Validation the Nextflow port by running the **same subsampled reads** through both
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
```

On bioHPC servers at Cornell, load the Nextflow module with
```bash
module load nextflow/25.4.3
```

The test env pins the **same tool versions** used by the pipeline, so both the
bash script and the NF port trim/align identically (a fair comparison). The runs
are forced single-threaded to remove bwa's multi-thread nondeterminism.

## 2. Genome (mouse or human)

If you already have a BWA index + chromInfo, use it. Otherwise build one:

```bash
bin/setup_genome.sh dm6         # Drosophila: ~140 Mb, indexes in ~1-2 min (recommended)
# test/setup_genome.sh mm10      # or mm10 / hg38 (mammalian bwa index ~1 h, several GB)
```

**dm6 (Drosophila) is the easiest target**  beacuase it is tiny/fast.

It prints the `--bwa-index` / `--chrom-info` paths to pass along. For *real*
analysis (not just concordance) add the rDNA sequence (GenBank U13369.1) to the
FASTA before indexing — see the note in `setup_genome.sh`.

## 3. Datasets

Edit `datasets.tsv` and fill in real SRA run accessions (the shipped ones are
placeholders). The ones chosen are from the following GEO accessions:
- Single-end GRO-seq: [Fuda et al, Mol Cell Bio 2012](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77607)
- Single-end PRO-seq: [Duarte et al, Genes Dev 2016](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38748)
- Paired-end PRO-seq: [Judd et al, Genes Dev 20206](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149332)


| name           | assay | layout | report | strand   | why |
|----------------|-------|--------|--------|----------|-----|
| `gro_se`       | GRO   | SE     | 5prime | same     | GRO-seq 5′ path (`-G`) |
| `pro_se`       | PRO   | SE     | 5prime | same     | PRO-seq path (`-P`) |
| `prochro_pe_5p`| ChRO  | PE     | 5prime | same     | PE report-5′ (`--map5 true`) |
| `prochro_pe_3p`| ChRO  | PE     | 3prime | same     | PE report-3′ (`--map5 false`) |
| `prochro_pe_opp`| ChRO | PE     | 5prime | opposite | opposite-strand path (`-s`) |

The last three rows can share **one** PE accession — reads are fetched/subsampled
once and cached, then reused. `3prime` and `opposite` are PE-only.


## 4. Run

```bash
test/run_concordance.sh \
  --bwa-index  test/ref/dm6.fa.gz \
  --chrom-info test/ref/dm6.chromInfo \
  --reads 1000000 --seed 100
```

Per dataset it fetches + subsamples (fixed seed; PE mates subsampled with the same
seed so pairs stay matched), runs both pipelines into `test/out/<name>/{bash,nf}`,
and prints `PASS (identical)` or a `DIFF` summary per strand. Exit code is non-zero
if anything differs.

Results are also written to **`test/out/concordance_mqc.tsv`** — a per-`sample_strand`
table (`IDENTICAL` / `DIFF` / `MISSING` + interval counts) printed at the end and
saved to disk. It's in MultiQC custom-content format, so pointing MultiQC at
`test/out/` renders it as a "proseq2.0-nf concordance" table.

## Results from testing above samples
| Assay     | Single/Paired-end | Strand   | Bigwig |Result     | 
|-----------|-------------------|----------|--------|-----------|
| GRO-seq   | SE                | same     | plus   | IDENTICAL |
| GRO-seq   | SE                | same     | minus  | IDENTICAL |
| PRO-seq   | SE                | same     | plus   | IDENTICAL |
| PRO-seq   | SE                | same     | minus  | IDENTICAL |
| PRO-seq   | PE                | same     | plus   | IDENTICAL |
| PRO-seq   | PE                | same     | minus  | IDENTICAL |
| PRO-seq   | PE                | same     | plus   | IDENTICAL |
| PRO-seq   | PE                | same     | minus  | IDENTICAL |
| PRO-seq   | PE                | opposite | plus   | IDENTICAL |
| PRO-seq   | PE                | opposite | minus  | IDENTICAL |

