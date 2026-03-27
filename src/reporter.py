"""
 Module 5: Report Generator

Takes a list of DiagnosticResults and writes two outputs:
  1. bmi_summary.tsv  — machine-readable TSV, one row per BUSCO
  2. bmi_report.txt   — human-readable plain-text report with recommendations

The full BMI pipeline will add HTML/Markdown output here.
This prototype demonstrates the core reporting structure.
"""

from __future__ import annotations
import csv
from datetime import datetime
from typing import Optional
from src.aggregator import DiagnosticResult


# ── TSV report ─────────

TSV_COLUMNS = [
    "BUSCO_ID",
    "Status",
    "Top_Reason_Code",
    "Confidence",
    "All_Reason_Codes",
    "Recommendation",
]


def write_tsv(results: list[DiagnosticResult], output_path: str) -> None:
    """
    Write a TSV summary: one row per BUSCO with top reason code and confidence.

    Args:
        results     : list of DiagnosticResult from aggregator.py
        output_path : path to write the .tsv file
    """
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TSV_COLUMNS, delimiter="\t")
        writer.writeheader()

        for result in results:
            all_codes = "|".join(rc.code for rc in result.reason_codes)
            writer.writerow({
                "BUSCO_ID":        result.busco_id,
                "Status":          result.status,
                "Top_Reason_Code": result.top_code,
                "Confidence":      f"{result.top_confidence:.2f}",
                "All_Reason_Codes": all_codes,
                "Recommendation":  result.recommendation,
            })


# ── Plain-text report ─────────

def write_text_report(
    results: list[DiagnosticResult],
    output_path: str,
    busco_version: int = 5,
    lineage: Optional[str] = None,
) -> None:
    """
    Write a human-readable report: overview statistics followed by
    per-BUSCO detail sections with ranked reason codes and recommendations.

    Args:
        results      : list of DiagnosticResult from aggregator.py
        output_path  : path to write the .txt report
        busco_version: BUSCO version detected by parser
        lineage      : lineage dataset name if available
    """
    # ── Summary statistics ──────
    total        = len(results)
    missing      = sum(1 for r in results if r.status == "Missing")
    fragmented   = sum(1 for r in results if r.status == "Fragmented")

    code_counts: dict[str, int] = {}
    for r in results:
        code_counts[r.top_code] = code_counts.get(r.top_code, 0) + 1

    with open(output_path, "w") as fh:
        # Header
        fh.write("=" * 72 + "\n")
        fh.write("  BUSCO-Missing Investigator (BMI) — Diagnostic Report\n")
        fh.write(f"  Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        if lineage:
            fh.write(f"  Lineage dataset: {lineage}\n")
        fh.write(f"  BUSCO version: {busco_version}\n")
        fh.write("=" * 72 + "\n\n")

        # Overview
        fh.write("OVERVIEW\n")
        fh.write("-" * 40 + "\n")
        fh.write(f"  Total investigated:  {total}\n")
        fh.write(f"  Missing:             {missing}\n")
        fh.write(f"  Fragmented:          {fragmented}\n\n")

        fh.write("  Reason code distribution:\n")
        for code, count in sorted(code_counts.items(), key=lambda x: -x[1]):
            bar = "█" * count
            fh.write(f"    {code:<25} {bar} ({count})\n")
        fh.write("\n")

        # Per-BUSCO details
        fh.write("PER-BUSCO DIAGNOSTICS\n")
        fh.write("-" * 40 + "\n\n")

        for result in results:
            fh.write(f"  BUSCO: {result.busco_id}  [{result.status}]\n")
            fh.write(f"  Top diagnosis: {result.top_code}  "
                     f"(confidence: {result.top_confidence:.2f})\n")

            if len(result.reason_codes) > 1:
                other = ", ".join(
                    f"{rc.code} ({rc.confidence:.2f})"
                    for rc in result.reason_codes[1:]
                )
                fh.write(f"  Other signals:  {other}\n")

            fh.write(f"\n  Recommendation:\n")
            # Word-wrap the recommendation at 65 chars
            words = result.recommendation.split()
            line  = "    "
            for word in words:
                if len(line) + len(word) + 1 > 68:
                    fh.write(line + "\n")
                    line = "    " + word
                else:
                    line = line + (" " if line != "    " else "") + word
            fh.write(line + "\n")
            fh.write("\n" + "  " + "-" * 60 + "\n\n")

        fh.write("=" * 72 + "\n")
        fh.write("  End of BMI report\n")
        fh.write("=" * 72 + "\n")
