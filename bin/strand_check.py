#!/usr/bin/env python3
"""Compare RSeQC infer_experiment.py output against the configured strand geometry.

The primary read (SE read, or R1 in PE) is *sense* to genes iff the RNA 5' end is
at that read's 5' end (proseq2.0 `-G` / `--rna5 R1_5prime`). `-5/-3` (--map5) and
`-s` (--opposite_strand) are reporting choices and are NOT validated here.

Usage: strand_check.py INFER_EXPERIMENT_TXT EXPECTED SAMPLE SEQTYPE
  EXPECTED : 'sense' | 'antisense'   (what the configured flags imply)
  SEQTYPE  : 'SE' | 'PE'
Writes SAMPLE.strand_check.txt and echoes the verdict to stderr.
"""
import sys
import re

infile, expected, sample, seqtype = sys.argv[1:5]

sense = antisense = failed = 0.0
text = open(infile).read()
m = re.search(r"failed to determine:\s*([0-9.]+)", text)
if m:
    failed = float(m.group(1))
for pattern, frac in re.findall(r'explained by "([^"]+)":\s*([0-9.]+)', text):
    f = float(frac)
    # '++,--' (SE) or '1++,1--,2+-,2-+' (PE) == primary read sense to genes
    if pattern.startswith("++") or "1++" in pattern:
        sense += f
    else:
        antisense += f

determined = sense + antisense
exp_label = "SENSE" if expected == "sense" else "ANTISENSE"

lines = [
    f"Sample: {sample} ({seqtype})",
    f"Configured expectation: primary read {exp_label} to genes",
    f"Detected: sense={sense:.3f} antisense={antisense:.3f} undetermined={failed:.3f}",
]

if determined < 0.10:
    verdict = "UNDETERMINED"
    note = ("Too few reads assigned to gene models — check that --gene_bed matches "
            "this genome/annotation.")
else:
    rel = sense / determined
    if rel >= 0.6:
        detected = "sense"
    elif rel <= 0.4:
        detected = "antisense"
    else:
        detected = "ambiguous"

    if detected == "ambiguous":
        verdict = "WARN"
        note = (f"Orientation ambiguous (sense = {rel*100:.0f}% of determined reads); "
                "library may be unstranded or the annotation is noisy.")
    elif detected == expected:
        verdict = "PASS"
        note = "Data orientation matches the configured strand settings."
    else:
        verdict = "WARN"
        fix = ("flip -G/-P (--se_read / --assay)" if seqtype == "SE"
               else "swap --rna5/--rna3 (or the --assay preset), then re-check -s")
        note = (f"Data looks {detected.upper()} but the flags imply {exp_label}. "
                f"Likely wrong strand config: {fix}.")

lines.append(f"VERDICT: {verdict} - {note}")
report = "\n".join(lines) + "\n"

with open(f"{sample}.strand_check.txt", "w") as fh:
    fh.write(report)
sys.stderr.write(report)
if verdict == "WARN":
    sys.stderr.write(f"\n*** STRAND WARNING [{sample}]: {note}\n")
