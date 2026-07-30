#!/usr/bin/env python3
"""Compare fastp's de-novo-detected adapter against the configured --adapter_* .

Validate-and-warn only: writes SAMPLE.adapter_check.txt with a
PASS / WARN / UNDETECTED verdict and echoes it to stderr. Never changes trimming.

Usage: adapter_check.py FASTP_JSON SAMPLE SE|PE CONFIGURED_R1 [CONFIGURED_R2]
  CONFIGURED_R1 : adapter expected on read1's 3' end  (SE: --adapter_se; PE: --adapter2)
  CONFIGURED_R2 : adapter expected on read2's 3' end  (PE: --adapter1)
"""
import json
import sys

fastp_json, sample, layout, conf_r1 = sys.argv[1:5]
conf_r2 = sys.argv[5] if len(sys.argv) > 5 else None

with open(fastp_json) as fh:
    data = json.load(fh)
ac = data.get("adapter_cutting", {})


def detected(mate):
    # fastp key names have varied across versions; try the common ones
    for k in (f"{mate}_adapter_sequence", "adapter_sequence"):
        v = ac.get(k)
        if v and v not in ("", "no adapter"):
            return v
    return None


def verdict(det, conf):
    if not det:
        return "UNDETECTED", f"fastp found no adapter (low read-through?); configured={conf}"
    n = min(len(det), len(conf), 12)
    if n >= 6 and det[:n].upper() == conf[:n].upper():
        return "PASS", f"detected {det} matches configured {conf}"
    return "WARN", f"detected {det} does NOT match configured {conf}"


lines = [f"Sample: {sample} ({layout})"]
results = [("read1", detected("read1"), conf_r1)]
if layout == "PE":
    results.append(("read2", detected("read2"), conf_r2))

worst = "PASS"
for mate, det, conf in results:
    v, note = verdict(det, conf)
    lines.append(f"{mate}: VERDICT: {v} - {note}")
    if v == "WARN":
        worst = "WARN"
    elif v == "UNDETECTED" and worst != "WARN":
        worst = "UNDETECTED"

report = "\n".join(lines) + "\n"
with open(f"{sample}.adapter_check.txt", "w") as fh:
    fh.write(report)
sys.stderr.write(report)
if worst == "WARN":
    sys.stderr.write(f"\n*** ADAPTER WARNING [{sample}]: detected adapter disagrees with the "
                     f"configured --adapter_* — check your adapter settings.\n")
