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

echo "== proseq2.0-nf concordance test =="

[ -n "$BWA_INDEX" ]  || { echo "ERROR: --bwa-index is required" >&2; exit 1; }
[ -n "$CHROM_INFO" ] || { echo "ERROR: --chrom-info is required" >&2; exit 1; }
[ -f "$CHROM_INFO" ] || { echo "ERROR: chrom-info not found: $CHROM_INFO" >&2; exit 1; }
[ -e "${BWA_INDEX}.bwt" ] || { echo "ERROR: BWA index not found at prefix '${BWA_INDEX}' (expected ${BWA_INDEX}.bwt). Build it with test/setup_genome.sh, or fix --bwa-index." >&2; exit 1; }

# The pipelines run from other working dirs, so make these absolute now.
CHROM_INFO="$(cd "$(dirname "$CHROM_INFO")" && pwd)/$(basename "$CHROM_INFO")"
BWA_INDEX="$(cd "$(dirname "$BWA_INDEX")" && pwd)/$(basename "$BWA_INDEX")"

for tool in nextflow seqtk fasterq-dump bwa samtools bigWigToBedGraph; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: '$tool' not on PATH — did you 'conda activate proseq2.0-nf-test'?" >&2; exit 1; }
done

mkdir -p "$OUTDIR"; OUTDIR="$(cd "$OUTDIR" && pwd)"
READSDIR="$OUTDIR/reads"; mkdir -p "$READSDIR"

# ---- map assay+layout(+report) to matched bash / NF flags ----
map_flags() {  # $1=assay $2=layout $3=report $4=strand ; sets BASH_FLAGS, NF_FLAGS
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
  if [ "$4" = "opposite" ]; then
    [ "$2" = "PE" ] || { echo "ERROR: strand=opposite is PE-only ($1:$2)" >&2; return 1; }
    BASH_FLAGS="$BASH_FLAGS -s"; NF_FLAGS="$NF_FLAGS --opposite_strand true"
  fi
}

# ---- fetch + subsample ----
# Subsampled reads are cached per accession+seed; each dataset row links its
# <name> files to that cache, so an accession reused across rows (e.g. PE 5'/3'/
# opposite) is fetched and subsampled only once. `accession` may be an SRR id
# (fetched) or a local prefix: <prefix>.fastq.gz (SE) / <prefix>_R{1,2}.fastq.gz (PE).
subsample() {  # $1=name $2=layout $3=accession
  local name="$1" layout="$2" acc="$3" key="${3##*/}"
  if [ "$layout" = "SE" ]; then
    local sub="$READSDIR/${key}.sub.fastq.gz"
    if [ ! -s "$sub" ]; then
      local raw
      if   [ -s "${acc}.fastq.gz" ]; then raw="${acc}.fastq.gz"
      elif [ -s "${acc}.fastq" ];    then raw="${acc}.fastq"
      else fasterq-dump --split-files -O "$READSDIR" "$acc"; raw="$READSDIR/${acc}.fastq"; fi
      seqtk sample -s"$SEED" "$raw" "$READS" | gzip > "$sub"
      case "$raw" in "$READSDIR/"*) rm -f "$raw" ;; esac   # only remove what we fetched
    else echo "  reads cached (${key})"; fi
    ln -sf "$(basename "$sub")" "$READSDIR/${name}.fastq.gz"
  else
    local s1="$READSDIR/${key}.sub_R1.fastq.gz" s2="$READSDIR/${key}.sub_R2.fastq.gz"
    if [ ! -s "$s1" ]; then
      local r1 r2
      if [ -s "${acc}_R1.fastq.gz" ]; then r1="${acc}_R1.fastq.gz"; r2="${acc}_R2.fastq.gz"
      else fasterq-dump --split-files -O "$READSDIR" "$acc"; r1="$READSDIR/${acc}_1.fastq"; r2="$READSDIR/${acc}_2.fastq"; fi
      # SAME seed on both mates keeps pairs matched
      seqtk sample -s"$SEED" "$r1" "$READS" | gzip > "$s1"
      seqtk sample -s"$SEED" "$r2" "$READS" | gzip > "$s2"
      case "$r1" in "$READSDIR/"*) rm -f "$r1" "$r2" ;; esac
    else echo "  reads cached (${key})"; fi
    ln -sf "$(basename "$s1")" "$READSDIR/${name}_R1.fastq.gz"
    ln -sf "$(basename "$s2")" "$READSDIR/${name}_R2.fastq.gz"
  fi
}

# ---- run the original bash pipeline ----
run_bash() {  # $1=name $2=layout ; uses BASH_FLAGS
  local name="$1" layout="$2" o="$OUTDIR/$name/bash"
  rm -rf "$o"; mkdir -p "$o"
  ( cd "$READSDIR" && bash "$REPO/proseq2.0.bsh" $BASH_FLAGS \
      -i "$BWA_INDEX" -c "$CHROM_INFO" -I "$name" -T "$o/tmp" -O "$o" --thread=1 )
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
  # No -profile conda: use the tools already on PATH (the activated test env),
  # so BOTH pipelines run identical tool versions and we skip a slow per-run
  # conda-env build. (For real analysis runs, use -profile conda/mamba instead.)
  ( cd "$o" && nextflow run "$REPO/main.nf" $NF_FLAGS \
      --input "$sheet" --bwa_index "$BWA_INDEX" --chrom_info "$CHROM_INFO" \
      --outdir "$o/results" --max_cpus 1 --skip_fastqc --skip_multiqc )
}

# ---- compare one strand's bigWig; returns 0 on exact match ----
# Match by CONTENT, not filename: the original leaves a spurious _QC / _QC_end in
# its output names (fragile name-mangling in proseq2.0.bsh), while the port names
# them <sample>_<strand>.bw. So glob each side's raw (non-.rpm) bigWig for the strand.
compare_bw() {  # $1=name $2=strand(plus|minus) ; appends a row to $SUMMARY
  local name="$1" s="$2" a b
  a=$(ls "$OUTDIR/$name/bash/"*_"${s}.bw" 2>/dev/null | grep -v '\.rpm\.bw$' | head -1)
  b=$(ls "$OUTDIR/$name/nf/results/bigwig/"*_"${s}.bw" 2>/dev/null | grep -v '\.rpm\.bw$' | head -1)
  if [ -z "$a" ] || [ -z "$b" ]; then
    printf "    %-6s MISSING (bash:%s nf:%s)\n" "$s" "$( [ -n "$a" ] && echo ok || echo -- )" "$( [ -n "$b" ] && echo ok || echo -- )"
    printf '%s\tMISSING\t0\t0\t0\n' "${name}_${s}" >> "$SUMMARY"
    return 1
  fi
  local ta="$OUTDIR/$name/_${s}_bash.bg" tb="$OUTDIR/$name/_${s}_nf.bg"
  bigWigToBedGraph "$a" "$ta"; bigWigToBedGraph "$b" "$tb"
  local la lb both
  la=$(wc -l < "$ta"); lb=$(wc -l < "$tb")
  if cmp -s "$ta" "$tb"; then
    printf "    %-6s PASS (identical)\n" "$s"
    printf '%s\tIDENTICAL\t%s\t%s\t%s\n' "${name}_${s}" "$la" "$lb" "$la" >> "$SUMMARY"
    rm -f "$ta" "$tb"; return 0
  fi
  both=$(comm -12 <(sort "$ta") <(sort "$tb") | wc -l)
  printf "    %-6s DIFF  bash=%s nf=%s intervals; %s identical -> inspect %s vs %s\n" \
    "$s" "$la" "$lb" "$both" "$ta" "$tb"
  printf '%s\tDIFF\t%s\t%s\t%s\n' "${name}_${s}" "$la" "$lb" "$both" >> "$SUMMARY"
  return 1
}

# ---- main loop ----
echo "Repo:       $REPO"
echo "BWA index:  $BWA_INDEX"
echo "Reads/seed: $READS / $SEED"
echo "Datasets:   $DATASETS"
echo

# Concordance report as MultiQC custom content (also a plain TSV you can open).
SUMMARY="$OUTDIR/concordance_mqc.tsv"
{
  echo "# id: 'proseq_concordance'"
  echo "# section_name: 'proseq2.0-nf concordance (port vs original)'"
  echo "# description: 'bigWig agreement between the Nextflow port and proseq2.0.bsh on identical subsampled reads.'"
  echo "# plot_type: 'table'"
  echo "# pconfig:"
  echo "#     id: 'proseq_concordance_table'"
  printf 'Sample\tResult\tbash intervals\tnf intervals\tshared\n'
} > "$SUMMARY"

fails=0
while read -r name assay layout accession report strand || [ -n "${name:-}" ]; do
  [ -z "${name:-}" ] && continue
  case "$name" in \#*) continue ;; esac
  report="${report:-5prime}"
  strand="${strand:-same}"
  echo "=== $name  ($assay $layout report=$report strand=$strand, $accession) ==="
  map_flags "$assay" "$layout" "$report" "$strand"
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

echo "Concordance summary  (report: $SUMMARY):"
grep -v '^#' "$SUMMARY" | column -t -s $'\t' | sed 's/^/  /'
echo

if [ "$fails" -eq 0 ]; then
  echo "ALL bigWigs identical — port matches proseq2.0.bsh."
else
  echo "$fails strand comparison(s) not identical — see the report / .bg files above."
  echo "(Small diffs can come from tool-version drift; large diffs mean a logic bug.)"
  exit 1
fi
