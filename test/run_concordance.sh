#!/usr/bin/env bash
#
# run_concordance.sh — validate the Nextflow port against the original proseq2.0.bsh
#
# For each dataset in a TSV, this will:
#   1. fetch the SRA run and subsample to N reads (fixed seed; PE pairs stay in sync)
#   2. run the ORIGINAL proseq2.0.bsh  and  the Nextflow port, with matched
#      parameters, tool versions (same conda env), and a single thread
#   3. diff the resulting bigWigs and print PASS / DIFF per strand
#
# Prereqs (activate the test env first):
#   conda env create -f test/environment.yml && conda activate proseq2.0-nf-test
#   plus: nextflow on PATH, and a BWA index (+ chromInfo) — see test/README.md.
#
# Usage:
#   test/run_concordance.sh \
#       --bwa-index /ref/mm10/mm10.rRNA \
#       --chrom-info /ref/mm10/mm10.chromInfo \
#       [--datasets test/datasets.tsv] [--reads 1000000] [--seed 100] [--outdir test/out]
#
set -euo pipefail

# ---- locate repo root (this script lives in <repo>/test) ----
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/.." && pwd)"

# ---- defaults ----
DATASETS="$SCRIPT_DIR/datasets.tsv"
READS=1000000
SEED=100
OUTDIR="$SCRIPT_DIR/out"
BWA_INDEX=""
CHROM_INFO=""

# ---- args ----
while [ $# -gt 0 ]; do
  case "$1" in
    --bwa-index)  BWA_INDEX="$2"; shift 2 ;;
    --chrom-info) CHROM_INFO="$2"; shift 2 ;;
    --datasets)   DATASETS="$2"; shift 2 ;;
    --reads)      READS="$2"; shift 2 ;;
    --seed)       SEED="$2"; shift 2 ;;
    --outdir)     OUTDIR="$2"; shift 2 ;;
    -h|--help)    grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

[ -n "$BWA_INDEX" ]  || { echo "ERROR: --bwa-index is required" >&2; exit 1; }
[ -n "$CHROM_INFO" ] || { echo "ERROR: --chrom-info is required" >&2; exit 1; }
[ -f "$CHROM_INFO" ] || { echo "ERROR: chrom-info not found: $CHROM_INFO" >&2; exit 1; }

for tool in nextflow seqtk fasterq-dump bwa samtools bigWigToBedGraph; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: '$tool' not on PATH (activate test env)" >&2; exit 1; }
done

mkdir -p "$OUTDIR"
READSDIR="$OUTDIR/reads"; mkdir -p "$READSDIR"

# ---- map assay+layout(+report) to matched bash / NF flags ----
map_flags() {  # $1=assay $2=layout $3=report ; sets BASH_FLAGS, NF_FLAGS
  case "$1:$2" in
    GRO:SE)          BASH_FLAGS="-SE -G";               NF_FLAGS="--assay GRO"  ;;
    PRO:SE)          BASH_FLAGS="-SE -P";               NF_FLAGS="--assay PRO"  ;;
    GRO:PE)          BASH_FLAGS="-PE";                  NF_FLAGS="--assay GRO"  ;;
    PRO:PE|ChRO:PE)  BASH_FLAGS="-PE --RNA3=R1_5prime"; NF_FLAGS="--assay ChRO" ;;
    *) echo "ERROR: unsupported assay:layout '$1:$2'" >&2; return 1 ;;
  esac
  if [ "$3" = "3prime" ]; then
    [ "$2" = "PE" ] || { echo "ERROR: report=3prime is PE-only ($1:$2)" >&2; return 1; }
    BASH_FLAGS="$BASH_FLAGS -3"; NF_FLAGS="$NF_FLAGS --map5 false"
  fi
}

# ---- fetch + subsample (idempotent: skips if outputs exist) ----
subsample() {  # $1=name $2=layout $3=accession
  local name="$1" layout="$2" acc="$3"
  if [ "$layout" = "SE" ]; then
    [ -s "$READSDIR/${name}.fastq.gz" ] && { echo "  reads cached"; return; }
    fasterq-dump --split-files -O "$READSDIR" "$acc"
    seqtk sample -s"$SEED" "$READSDIR/${acc}.fastq" "$READS" | gzip > "$READSDIR/${name}.fastq.gz"
    rm -f "$READSDIR/${acc}.fastq"
  else
    [ -s "$READSDIR/${name}_R1.fastq.gz" ] && { echo "  reads cached"; return; }
    fasterq-dump --split-files -O "$READSDIR" "$acc"
    # SAME seed on both mates keeps pairs matched
    seqtk sample -s"$SEED" "$READSDIR/${acc}_1.fastq" "$READS" | gzip > "$READSDIR/${name}_R1.fastq.gz"
    seqtk sample -s"$SEED" "$READSDIR/${acc}_2.fastq" "$READS" | gzip > "$READSDIR/${name}_R2.fastq.gz"
    rm -f "$READSDIR/${acc}_1.fastq" "$READSDIR/${acc}_2.fastq"
  fi
}

# ---- run the original bash pipeline ----
run_bash() {  # $1=name $2=layout ; uses BASH_FLAGS
  local name="$1" layout="$2" o="$OUTDIR/$name/bash"
  rm -rf "$o"; mkdir -p "$o"
  ( cd "$READSDIR" && bash "$REPO/proseq2.0.bsh" $BASH_FLAGS \
      -i "$BWA_INDEX" -c "$CHROM_INFO" -I "$name" -T "$o/tmp" -O "$o" --thread 1 )
}

# ---- run the Nextflow port ----
run_nf() {  # $1=name $2=layout ; uses NF_FLAGS
  local name="$1" layout="$2" o="$OUTDIR/$name/nf"
  rm -rf "$o"; mkdir -p "$o"
  local sheet="$o/samplesheet.csv"
  echo "sample,fastq_1,fastq_2" > "$sheet"
  if [ "$layout" = "SE" ]; then
    echo "${name},${READSDIR}/${name}.fastq.gz," >> "$sheet"
  else
    echo "${name},${READSDIR}/${name}_R1.fastq.gz,${READSDIR}/${name}_R2.fastq.gz" >> "$sheet"
  fi
  ( cd "$o" && nextflow run "$REPO/main.nf" -profile conda $NF_FLAGS \
      --input "$sheet" --bwa_index "$BWA_INDEX" --chrom_info "$CHROM_INFO" \
      --outdir "$o/results" --max_cpus 1 --skip_fastqc --skip_multiqc )
}

# ---- compare one strand's bigWig; returns 0 on exact match ----
compare_bw() {  # $1=name $2=strand(plus|minus)
  local name="$1" s="$2"
  local a="$OUTDIR/$name/bash/${name}_${s}.bw"
  local b="$OUTDIR/$name/nf/results/bigwig/${name}_${s}.bw"
  if [ ! -f "$a" ] || [ ! -f "$b" ]; then
    printf "    %-6s MISSING (bash:%s nf:%s)\n" "$s" "$( [ -f "$a" ] && echo ok || echo -- )" "$( [ -f "$b" ] && echo ok || echo -- )"
    return 1
  fi
  local ta="$OUTDIR/$name/_${s}_bash.bg" tb="$OUTDIR/$name/_${s}_nf.bg"
  bigWigToBedGraph "$a" "$ta"; bigWigToBedGraph "$b" "$tb"
  if cmp -s "$ta" "$tb"; then
    printf "    %-6s PASS (identical)\n" "$s"; rm -f "$ta" "$tb"; return 0
  fi
  local la lb both
  la=$(wc -l < "$ta"); lb=$(wc -l < "$tb")
  both=$(comm -12 <(sort "$ta") <(sort "$tb") | wc -l)
  printf "    %-6s DIFF  bash=%s nf=%s intervals; %s identical -> inspect %s vs %s\n" \
    "$s" "$la" "$lb" "$both" "$ta" "$tb"
  return 1
}

# ---- main loop ----
echo "Repo:       $REPO"
echo "BWA index:  $BWA_INDEX"
echo "Reads/seed: $READS / $SEED"
echo "Datasets:   $DATASETS"
echo

fails=0
while read -r name assay layout accession report; do
  [ -z "${name:-}" ] && continue
  case "$name" in \#*) continue ;; esac
  report="${report:-5prime}"
  echo "=== $name  ($assay $layout report=$report, $accession) ==="
  map_flags "$assay" "$layout" "$report"
  echo "  bash: $BASH_FLAGS"
  echo "  nf:   $NF_FLAGS"
  subsample "$name" "$layout" "$accession"
  run_bash "$name" "$layout"
  run_nf   "$name" "$layout"
  echo "  concordance:"
  compare_bw "$name" plus  || fails=$((fails+1))
  compare_bw "$name" minus || fails=$((fails+1))
  echo
done < "$DATASETS"

if [ "$fails" -eq 0 ]; then
  echo "ALL bigWigs identical — port matches proseq2.0.bsh."
else
  echo "$fails strand comparison(s) differed — see the .bg files noted above."
  echo "(Small diffs can come from tool-version drift; large diffs mean a logic bug.)"
  exit 1
fi
