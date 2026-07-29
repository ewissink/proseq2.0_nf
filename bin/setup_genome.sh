#!/usr/bin/env bash
#
# setup_genome.sh — download a genome and build a BWA index + chromInfo.
#
# Usage: bin/setup_genome.sh {dm6|mm10|hg38} [OUTDIR]
#
# Alt haplotypes: for hg38 this pulls UCSC's "analysisSet" FASTA, which has the
# alternate haplotypes and patch scaffolds already removed (the assembly UCSC
# recommends for read alignment). Unplaced/unlocalized scaffolds are kept on
# purpose, so reads that belong there don't mis-map onto the main chromosomes.
# mm10 and dm6 have no alt-haplotype contigs, so their standard assemblies are
# used as-is.
#
# chromInfo is generated FROM the downloaded FASTA, so it always matches the
# exact set of sequences in the index (important for bedGraphToBigWig).
#
# dm6 (Drosophila) is tiny (~140 Mb) -> indexes in ~1-2 min. Mammalian
# `bwa index` takes ~1 h and several GB.
#
# rDNA caveat: none of these UCSC downloads include the nuclear rDNA array, so
# rRNA reads won't map/filter cleanly. For rRNA handling either append an rDNA
# contig (named with 'rRNA') to the FASTA before indexing, or use the pipeline's
# --remove_rrna (SortMeRNA).
#
set -euo pipefail

GENOME="${1:-}"
OUTDIR="${2:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/test/ref}"
[ -n "$GENOME" ] || { echo "Usage: $0 {dm6|mm10|hg38} [OUTDIR]" >&2; exit 1; }

UCSC="https://hgdownload.soe.ucsc.edu/goldenPath"
case "$GENOME" in
  dm6)  URL="$UCSC/dm6/bigZips/dm6.fa.gz";                        FA="dm6.fa.gz"                ;;
  mm10) URL="$UCSC/mm10/bigZips/mm10.fa.gz";                      FA="mm10.fa.gz"               ;;
  hg38) URL="$UCSC/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz"; FA="hg38.analysisSet.fa.gz" ;;
  *) echo "ERROR: genome must be dm6, mm10 or hg38" >&2; exit 1 ;;
esac

for tool in wget bwa; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: '$tool' not on PATH (activate the env)" >&2; exit 1; }
done

mkdir -p "$OUTDIR"; cd "$OUTDIR"

echo ">> Downloading $GENOME FASTA ($FA) ..."
[ -s "$FA" ] || wget -q --show-progress "$URL"

CHROMINFO="${FA%.fa.gz}.chromInfo"
echo ">> Generating chromInfo from the FASTA (matches the index exactly) ..."
gzip -dc "$FA" \
  | awk '/^>/{if(name!="")print name"\t"n; name=substr($1,2); n=0; next}{n+=length($0)} END{if(name!="")print name"\t"n}' \
  > "$CHROMINFO"
echo "   $(wc -l < "$CHROMINFO") sequences"

echo ">> Building BWA index (this is the slow part) ..."
[ -s "${FA}.bwt" ] || bwa index "$FA"

echo
echo "Done. Pass these to the pipeline (or the test harness):"
echo "  --bwa_index  $OUTDIR/$FA          (--bwa-index for run_concordance.sh)"
echo "  --chrom_info $OUTDIR/$CHROMINFO   (--chrom-info for run_concordance.sh)"
