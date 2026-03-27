"""
Module 4: Evidence Aggregator & Reason Code Scorer

Takes a BuscoEntry and a list of EvidenceResults from the diagnostic
modules and produces a DiagnosticResult: a ranked list of reason codes
with confidence scores and a plain-English recommendation.

This is the decision layer of BMI — it combines all signals into a
single, interpretable output per BUSCO gene.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional
from src.parser import BuscoEntry
from src.diagnostics import EvidenceResult


# ── Reason code definitions ──────

REASON_CODES: dict[str, str] = {
    "CONTIG_EDGE": (
        "The BUSCO hit is near a contig boundary — the gene may be split "
        "across contigs. Consider scaffolding with long reads or Hi-C data, "
        "or rerunning assembly with a long-read assembler."
    ),
    "ASSEMBLY_GAP": (
        "An N-gap is present at or near the BUSCO locus. This suggests an "
        "assembly gap disrupted the gene. Consider gap-filling tools such as "
        "TGS-GapCloser or LR_Gapcloser."
    ),
    "LOW_SCORE": (
        "The BUSCO alignment score is low. This may indicate a lineage dataset "
        "mismatch — try rerunning BUSCO with a closer lineage dataset for your "
        "species, or check for high sequence divergence at this locus."
    ),
    "NO_HOMOLOGY": (
        "No sequence coordinates were found — the gene may be genuinely absent "
        "from the assembly. Verify with BLAST against a closely related species, "
        "or check for true biological gene loss."
    ),
    "FRAGMENTED_SHORT": (
        "The BUSCO hit is very short relative to the expected gene length. "
        "The gene model may be truncated due to a premature stop codon, "
        "exon skipping, or annotation failure. Consider rerunning annotation."
    ),
    "INVESTIGATION_NEEDED": (
        "No strong diagnostic signal was detected with the available inputs. "
        "Providing an assembly FASTA, annotation GFF, and/or a coverage BAM "
        "would allow more specific diagnostics."
    ),
}

# Score thresholds — configurable in the full pipeline via nextflow.config
LOW_SCORE_THRESHOLD      = 50.0   
SHORT_HIT_FRACTION       = 0.3    
EXPECTED_BUSCO_LENGTH_BP = 500     


# ── Diagnostic result model ───────

@dataclass
class ReasonCode:
    """One reason code with its confidence score."""
    code: str
    confidence: float
    recommendation: str


@dataclass
class DiagnosticResult:
    """
    Full diagnostic output for one BUSCO gene.
    Contains ranked reason codes and the top recommendation.
    """
    busco_id: str
    status: str
    reason_codes: list[ReasonCode] = field(default_factory=list)
    evidence_summary: list[str]    = field(default_factory=list)

    @property
    def top_reason(self) -> Optional[ReasonCode]:
        """The highest-confidence reason code, or None if no signals fired."""
        return self.reason_codes[0] if self.reason_codes else None

    @property
    def top_code(self) -> str:
        return self.top_reason.code if self.top_reason else "INVESTIGATION_NEEDED"

    @property
    def top_confidence(self) -> float:
        return self.top_reason.confidence if self.top_reason else 0.0

    @property
    def recommendation(self) -> str:
        return (
            self.top_reason.recommendation
            if self.top_reason
            else REASON_CODES["INVESTIGATION_NEEDED"]
        )


# ── Aggregator ─────────

def aggregate(
    entry: BuscoEntry,
    evidence: list[EvidenceResult],
) -> DiagnosticResult:
    """
    Combine a BuscoEntry and its evidence signals into a DiagnosticResult.

    Logic:
      1. Check evidence signals from diagnostics.py (contig edge, gap, etc.)
      2. Check score-based signals directly on the BuscoEntry
      3. Check coordinate-based signals (no homology, very short hit)
      4. Rank all fired signals by confidence (descending)
      5. Fall back to INVESTIGATION_NEEDED if nothing fired

    Args:
        entry    : BuscoEntry from parser.py
        evidence : list of EvidenceResult from diagnostics.py functions

    Returns:
        DiagnosticResult with ranked reason codes and recommendations
    """
    fired_codes: list[ReasonCode] = []
    summary: list[str] = []

    # ── Evidence-signal reason codes ───────────
    for ev in evidence:
        summary.append(f"[{ev.signal_name}] fired={ev.fired} conf={ev.confidence:.2f} — {ev.detail}")
        if ev.fired:
            fired_codes.append(ReasonCode(
                code=ev.signal_name,
                confidence=ev.confidence,
                recommendation=REASON_CODES.get(ev.signal_name, REASON_CODES["INVESTIGATION_NEEDED"]),
            ))

    # ── Score-based reason codes ───────────
    if entry.score is not None and entry.score < LOW_SCORE_THRESHOLD:
        # Scale confidence: score=0 → 0.9, score=threshold → 0.5
        conf = round(0.9 - (entry.score / LOW_SCORE_THRESHOLD) * 0.4, 2)
        fired_codes.append(ReasonCode(
            code="LOW_SCORE",
            confidence=conf,
            recommendation=REASON_CODES["LOW_SCORE"],
        ))
        summary.append(
            f"[LOW_SCORE] fired=True conf={conf:.2f} — "
            f"BUSCO score {entry.score:.1f} is below threshold {LOW_SCORE_THRESHOLD}."
        )

    # ── Coordinate-based reason codes ─────

    # No coordinates at all → gene not found anywhere in assembly
    if entry.status == "Missing" and not entry.has_coordinates:
        fired_codes.append(ReasonCode(
            code="NO_HOMOLOGY",
            confidence=0.75,
            recommendation=REASON_CODES["NO_HOMOLOGY"],
        ))
        summary.append(
            "[NO_HOMOLOGY] fired=True conf=0.75 — "
            "Missing entry with no sequence coordinates."
        )

    # Very short hit relative to expected length to truncated gene model
    if entry.has_coordinates and entry.gene_length is not None:
        expected = entry.length if entry.length else EXPECTED_BUSCO_LENGTH_BP
        if expected > 0 and (entry.gene_length / expected) < SHORT_HIT_FRACTION:
            conf = round(
                0.6 + (1.0 - entry.gene_length / (expected * SHORT_HIT_FRACTION)) * 0.3,
                2,
            )
            conf = min(conf, 0.9)
            fired_codes.append(ReasonCode(
                code="FRAGMENTED_SHORT",
                confidence=conf,
                recommendation=REASON_CODES["FRAGMENTED_SHORT"],
            ))
            summary.append(
                f"[FRAGMENTED_SHORT] fired=True conf={conf:.2f} — "
                f"Gene span {entry.gene_length} bp is < {SHORT_HIT_FRACTION*100:.0f}% "
                f"of expected {expected} bp."
            )

    # ── Rank by confidence descending ─────────
    fired_codes.sort(key=lambda r: r.confidence, reverse=True)

    # ── Fallback if nothing fired ────────
    if not fired_codes:
        fired_codes.append(ReasonCode(
            code="INVESTIGATION_NEEDED",
            confidence=0.0,
            recommendation=REASON_CODES["INVESTIGATION_NEEDED"],
        ))
        summary.append(
            "[INVESTIGATION_NEEDED] No diagnostic signals fired with available inputs."
        )

    return DiagnosticResult(
        busco_id=entry.busco_id,
        status=entry.status,
        reason_codes=fired_codes,
        evidence_summary=summary,
    )
