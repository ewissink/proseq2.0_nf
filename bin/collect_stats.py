#!/usr/bin/env python3
"""Merge per-sample *.metrics files into a single MultiQC custom-content table.

Each *.metrics file holds tab-separated key<TAB>value lines and must include a
`sample` key. Files for the same sample (e.g. one from preprocessing, one from
bigWig generation) are merged by sample id.

Usage: collect_stats.py [DIR]   (DIR defaults to '.')
"""
import sys
import os
import glob

metrics_dir = sys.argv[1] if len(sys.argv) > 1 else "."

# column key -> display name, in the order they should appear
COLUMNS = [
    ("seq_type",         "Seq type"),
    ("input",            "Input (reads/pairs)"),
    ("pass_qc",          "Pass QC"),
    ("mappable",         "Mappable"),
    ("mappable_no_rrna", "Mappable (excl. rRNA)"),
]

rows = {}
for path in glob.glob(os.path.join(metrics_dir, "*.metrics")):
    kv = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                kv[parts[0]] = parts[1]
    sample = kv.get("sample")
    if sample:
        rows.setdefault(sample, {}).update(kv)

header = [
    "# id: 'proseq_read_stats'",
    "# section_name: 'proseq2.0 read summary'",
    "# description: 'Read counts through preprocessing, QC and alignment.'",
    "# plot_type: 'table'",
    "# pconfig:",
    "#     id: 'proseq_read_stats_table'",
    "#     namespace: 'proseq2.0'",
]
print("\n".join(header))
print("\t".join(["Sample"] + [disp for _, disp in COLUMNS]))
for sample in sorted(rows):
    kv = rows[sample]
    print("\t".join([sample] + [str(kv.get(key, "")) for key, _ in COLUMNS]))
