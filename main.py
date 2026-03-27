"""
BMI Prototype Entry Point

Runs the full prototype pipeline:
  1. Parse BUSCO full_table.tsv        (Module 1: parser.py)
  2. Run genome-level diagnostics       (Module 2: diagnostics.py)
  3. Aggregate signals → reason codes   (Module 4: aggregator.py)
  4. Write TSV + text report            (Module 5: reporter.py)

Usage:
    python main.py --busco data/full_table.tsv
    python main.py --busco data/full_table.tsv --assembly data/assembly.fasta
    python main.py --busco data/full_table.tsv --assembly data/assembly.fasta --out results/
"""

from __future__ import annotations
import argparse
import os
import sys

from src.parser      import parse_busco_table
from src.diagnostics import load_assembly, check_contig_edge, check_assembly_gap
from src.aggregator  import aggregate
from src.reporter    import write_tsv, write_text_report


def run(busco_path: str, assembly_path: str | None, out_dir: str) -> None:
    """Run the full BMI prototype pipeline."""

    os.makedirs(out_dir, exist_ok=True)

    # ── Step 1: Parse BUSCO output ────────────────────────────────────────────
    print(f"[BMI] Parsing BUSCO table: {busco_path}")
    entries, version = parse_busco_table(busco_path)
    print(f"[BMI] Found {len(entries)} Missing/Fragmented BUSCOs (BUSCO v{version})")

    # ── Step 2: Load assembly (optional) ─────────────────────────────────────
    assembly: dict[str, str] = {}
    if assembly_path:
        print(f"[BMI] Loading assembly: {assembly_path}")
        assembly = load_assembly(assembly_path)
        print(f"[BMI] Loaded {len(assembly)} contigs")
    else:
        print("[BMI] No assembly provided — skipping genome-level diagnostics")

    # ── Step 3: Run diagnostics + aggregate ──────────────────────────────────
    print("[BMI] Running diagnostics...")
    results = []
    for entry in entries:
        evidence = []

        if assembly:
            evidence.append(check_contig_edge(entry, assembly))
            evidence.append(check_assembly_gap(entry, assembly))
        # In full BMI: BLAST search and coverage analysis would be added here

        result = aggregate(entry, evidence)
        results.append(result)

    # ── Step 4: Write reports ─────────────────────────────────────────────────
    tsv_path    = os.path.join(out_dir, "bmi_summary.tsv")
    report_path = os.path.join(out_dir, "bmi_report.txt")

    write_tsv(results, tsv_path)
    write_text_report(results, report_path, busco_version=version)

    print(f"[BMI] Done.")
    print(f"[BMI] Summary TSV → {tsv_path}")
    print(f"[BMI] Text report → {report_path}")

    # ── Print quick summary to terminal ──────────────────────────────────────
    print("\n" + "=" * 60)
    print("  BMI Quick Summary")
    print("=" * 60)
    print(f"  {'BUSCO ID':<18} {'Status':<12} {'Top Reason':<22} {'Conf'}")
    print("  " + "-" * 56)
    for r in results:
        print(f"  {r.busco_id:<18} {r.status:<12} {r.top_code:<22} {r.top_confidence:.2f}")
    print("=" * 60)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="BMI Prototype — BUSCO Missing Investigator"
    )
    parser.add_argument(
        "--busco", required=True,
        help="Path to BUSCO full_table.tsv"
    )
    parser.add_argument(
        "--assembly", default=None,
        help="Path to assembly FASTA (optional — enables genome-level diagnostics)"
    )
    parser.add_argument(
        "--out", default="results",
        help="Output directory for reports (default: results/)"
    )
    args = parser.parse_args()

    if not os.path.isfile(args.busco):
        print(f"[BMI] Error: BUSCO file not found: {args.busco}", file=sys.stderr)
        sys.exit(1)

    if args.assembly and not os.path.isfile(args.assembly):
        print(f"[BMI] Error: Assembly file not found: {args.assembly}", file=sys.stderr)
        sys.exit(1)

    run(args.busco, args.assembly, args.out)


if __name__ == "__main__":
    main()
