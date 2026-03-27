"""
Module 2: Genome-Level Diagnostics

Two evidence signals:
  1. Contig-edge proximity,,, is the BUSCO hit near a scaffold boundary?
  2. Assembly gap scanner,,, are there N-stretches at the BUSCO locus?

Both functions accept a BuscoEntry and the parsed assembly (a dict of
contig_name -> sequence string) and return an EvidenceResult.

These are the two fastest, assembly-only checks wiht no external tools needed.
BLAST and coverage analysis would be added here in the full BMI pipeline
as additional functions following the same pattern.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional
from src.parser import BuscoEntry


# ── Evidence result model 

@dataclass
class EvidenceResult:
    """
    Holds the raw evidence gathered for one diagnostic check on one BUSCO.

    signal_name   : short identifier for this evidence type
    fired         : True if the signal was triggered 
    confidence    : 0.0–1.0 confidence that this explains the missing BUSCO
    detail        : human-readable description of what was found
    """
    signal_name: str
    fired: bool
    confidence: float
    detail: str


# ── Assembly loader 

def load_assembly(fasta_path: str) -> dict[str, str]:
    """
    Parse a FASTA file into {contig_name: sequence} dict.
    Sequence is stored as uppercase string for consistent N-counting.

    Args:
        fasta_path: path to assembly FASTA file

    Returns:
        dict mapping contig name (first word after >) to full sequence
    """
    assembly: dict[str, str] = {}
    current_name: Optional[str] = None
    current_seq: list[str] = []

    with open(fasta_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                # Save previous contig
                if current_name is not None:
                    assembly[current_name] = "".join(current_seq).upper()
                # Start new contig — use first word of header as name
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    # Save last contig
    if current_name is not None:
        assembly[current_name] = "".join(current_seq).upper()

    return assembly


# ── Signal 1: Contig-edge proximity 

def check_contig_edge(
    entry: BuscoEntry,
    assembly: dict[str, str],
    edge_threshold_bp: int = 500,
) -> EvidenceResult:
    """
    Check whether the BUSCO hit is within `edge_threshold_bp` bases of
    either end of its contig. A hit near the edge suggests the gene may
    be split across a contig boundary — a classic assembly fragmentation signal.

    Args:
        entry             : BuscoEntry with coordinates
        assembly          : {contig_name: sequence} from load_assembly()
        edge_threshold_bp : distance from contig end to flag as "near edge"

    Returns:
        EvidenceResult with fired=True if near edge
    """
    signal = "CONTIG_EDGE"

    # Can only check if we have coordinates and the contig is in the assembly
    if not entry.has_coordinates:
        return EvidenceResult(
            signal_name=signal,
            fired=False,
            confidence=0.0,
            detail="No coordinates available (Missing entry — cannot check edge proximity).",
        )

    contig_seq = assembly.get(entry.sequence)
    if contig_seq is None:
        return EvidenceResult(
            signal_name=signal,
            fired=False,
            confidence=0.0,
            detail=f"Contig '{entry.sequence}' not found in assembly FASTA.",
        )

    contig_length = len(contig_seq)
    dist_from_start = entry.start
    dist_from_end   = contig_length - entry.end

    near_start = dist_from_start <= edge_threshold_bp
    near_end   = dist_from_end   <= edge_threshold_bp
    fired      = near_start or near_end

    if fired:
        which = "start" if near_start else "end"
        dist  = dist_from_start if near_start else dist_from_end
        # Confidence scales with how close to the edge:
        # 0 bp away → 1.0, edge_threshold_bp away → 0.5 — clamped to [0.5, 1.0]
        ratio = min(dist / edge_threshold_bp, 1.0)
        confidence = round(1.0 - ratio * 0.5, 2)
        detail = (
            f"BUSCO hit is {dist} bp from the contig {which} "
            f"(contig length: {contig_length} bp, threshold: {edge_threshold_bp} bp). "
            f"Gene may be split across a contig boundary."
        )
    else:
        confidence = 0.0
        detail = (
            f"BUSCO hit is {dist_from_start} bp from start, "
            f"{dist_from_end} bp from end of contig. "
            f"Not near any contig edge (threshold: {edge_threshold_bp} bp)."
        )

    return EvidenceResult(
        signal_name=signal,
        fired=fired,
        confidence=confidence,
        detail=detail,
    )


# ── Signal 2: Assembly gap scanner 

def check_assembly_gap(
    entry: BuscoEntry,
    assembly: dict[str, str],
    min_n_run: int = 10,
    flank_bp: int = 200,
) -> EvidenceResult:
    """
    Scan the BUSCO locus (plus flanking regions) for runs of Ns in the
    assembly sequence. Long N-stretches represent assembly gaps and can
    explain why a gene model was not recovered.

    Args:
        entry      : BuscoEntry with coordinates
        assembly   : {contig_name: sequence} from load_assembly()
        min_n_run  : minimum consecutive Ns to count as a gap
        flank_bp   : bases to scan on either side of the BUSCO locus

    Returns:
        EvidenceResult with fired=True if gaps are found in or near locus
    """
    signal = "ASSEMBLY_GAP"

    if not entry.has_coordinates:
        return EvidenceResult(
            signal_name=signal,
            fired=False,
            confidence=0.0,
            detail="No coordinates available — cannot scan for gaps.",
        )

    contig_seq = assembly.get(entry.sequence)
    if contig_seq is None:
        return EvidenceResult(
            signal_name=signal,
            fired=False,
            confidence=0.0,
            detail=f"Contig '{entry.sequence}' not found in assembly FASTA.",
        )

    # Define the window to scan: locus ± flank_bp, clamped to contig bounds
    scan_start = max(0, entry.start - flank_bp)
    scan_end   = min(len(contig_seq), entry.end + flank_bp)
    window     = contig_seq[scan_start:scan_end]

    # Find all N-runs of at least min_n_run length
    n_runs: list[tuple[int, int]] = []  # (start, length) in window coords
    i = 0
    while i < len(window):
        if window[i] == "N":
            run_start = i
            while i < len(window) and window[i] == "N":
                i += 1
            run_length = i - run_start
            if run_length >= min_n_run:
                n_runs.append((scan_start + run_start, run_length))
        else:
            i += 1

    fired = len(n_runs) > 0

    if fired:
        total_n = sum(length for _, length in n_runs)
        locus_length = entry.end - entry.start
        # Confidence based on how much of the locus is covered by Ns
        gap_fraction = min(total_n / max(locus_length, 1), 1.0)
        confidence   = round(0.4 + gap_fraction * 0.6, 2)  # range: 0.4–1.0
        runs_desc = ", ".join(
            f"{length} Ns at position {pos}" for pos, length in n_runs
        )
        detail = (
            f"Found {len(n_runs)} N-run(s) in/near BUSCO locus "
            f"(scanned {scan_start}–{scan_end}): {runs_desc}. "
            f"Total gap bases: {total_n}."
        )
    else:
        confidence = 0.0
        detail = (
            f"No N-runs of >= {min_n_run} bp found in scan window "
            f"({scan_start}–{scan_end}). No assembly gap signal."
        )

    return EvidenceResult(
        signal_name=signal,
        fired=fired,
        confidence=confidence,
        detail=detail,
    )
