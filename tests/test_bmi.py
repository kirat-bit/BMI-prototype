"""
Test suite for the BMI prototype

Covers:
  - parser.py         : BuscoEntry creation, v4/v5 detection, edge cases
  - diagnostics.py    : contig-edge detection, gap detection
  - aggregator.py     : reason code assignment, ranking, fallback
  - reporter.py       : TSV and text report writing

Run with:  pytest tests/ -v
"""

from __future__ import annotations
import os
import csv
import tempfile
import pytest

from src.parser      import parse_busco_table, parse_row, detect_busco_version, BuscoEntry
from src.diagnostics import load_assembly, check_contig_edge, check_assembly_gap
from src.aggregator  import aggregate, DiagnosticResult
from src.reporter    import write_tsv, write_text_report


# ── Fixtures ──────────────────────────────────────────────────────────────────

SAMPLE_TSV_CONTENT = """\
# BUSCO version is: 5.4.3
# The lineage dataset is: insecta_odb10
# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength
EOG090W00A9\tComplete\tcontig_001\t145230\t147890\t+\t987.3\t892
EOG090W00D5\tMissing\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A
EOG090W00C3\tFragmented\tcontig_004\t3820\t4210\t+\t43.2\t187
EOG090W00G8\tFragmented\tcontig_007\t99821\t99995\t+\t38.7\t94
"""

SAMPLE_FASTA_CONTENT = """\
>contig_001 length=300000
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>contig_004 length=5000
TTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAAGGGG
>contig_007 length=100500
AAAACCCCTTTTGGGGAAAACCCCTTTTGGGGAAAANNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNAAAAGGGGAAAACCCCTTTTGGGG
"""


@pytest.fixture
def sample_tsv(tmp_path):
    """Write sample TSV to a temp file and return its path."""
    p = tmp_path / "full_table.tsv"
    p.write_text(SAMPLE_TSV_CONTENT)
    return str(p)


@pytest.fixture
def sample_fasta(tmp_path):
    """Write sample FASTA to a temp file and return its path."""
    p = tmp_path / "assembly.fasta"
    p.write_text(SAMPLE_FASTA_CONTENT)
    return str(p)


@pytest.fixture
def sample_assembly(sample_fasta):
    """Return a loaded assembly dict from the sample FASTA."""
    return load_assembly(sample_fasta)


@pytest.fixture
def missing_entry():
    """A Missing BuscoEntry with no coordinates."""
    return BuscoEntry(
        busco_id="EOG090W00D5",
        status="Missing",
        sequence=None,
        start=None,
        end=None,
        score=None,
        length=None,
        busco_version=5,
    )


@pytest.fixture
def fragmented_near_edge():
    """A Fragmented BuscoEntry near the end of contig_007 (length ~100500)."""
    return BuscoEntry(
        busco_id="EOG090W00G8",
        status="Fragmented",
        sequence="contig_007",
        start=99821,
        end=99995,
        score=38.7,
        length=94,
        busco_version=5,
    )


@pytest.fixture
def fragmented_low_score():
    """A Fragmented BuscoEntry on contig_004 with a low BUSCO score."""
    return BuscoEntry(
        busco_id="EOG090W00C3",
        status="Fragmented",
        sequence="contig_004",
        start=3820,
        end=4210,
        score=43.2,
        length=187,
        busco_version=5,
    )


# ── parser.py tests ───────────────────────────────────────────────────────────

class TestParser:

    def test_parse_returns_only_missing_and_fragmented(self, sample_tsv):
        entries, _ = parse_busco_table(sample_tsv)
        assert all(e.status in ("Missing", "Fragmented") for e in entries)

    def test_complete_entries_are_excluded(self, sample_tsv):
        entries, _ = parse_busco_table(sample_tsv)
        ids = [e.busco_id for e in entries]
        assert "EOG090W00A9" not in ids  # Complete — should be excluded

    def test_missing_entry_has_no_coordinates(self, sample_tsv):
        entries, _ = parse_busco_table(sample_tsv)
        missing = next(e for e in entries if e.status == "Missing")
        assert missing.sequence is None
        assert missing.start is None
        assert missing.end is None
        assert not missing.has_coordinates

    def test_fragmented_entry_has_coordinates(self, sample_tsv):
        entries, _ = parse_busco_table(sample_tsv)
        fragmented = [e for e in entries if e.status == "Fragmented"]
        assert len(fragmented) == 2
        assert all(e.has_coordinates for e in fragmented)

    def test_version_detection_v5(self, sample_tsv):
        _, version = parse_busco_table(sample_tsv)
        assert version == 5

    def test_version_detection_v4(self, tmp_path):
        v4_content = "# BUSCO version is: 4.1.0\n# Busco id\tStatus\tSequence\tStart\tEnd\tScore\tLength\nEOG090W00X1\tMissing\tN/A\tN/A\tN/A\tN/A\tN/A\n"
        p = tmp_path / "v4_table.tsv"
        p.write_text(v4_content)
        entries, version = parse_busco_table(str(p))
        assert version == 4
        assert len(entries) == 1

    def test_gene_length_property(self, sample_tsv):
        entries, _ = parse_busco_table(sample_tsv)
        fragmented = next(e for e in entries if e.sequence == "contig_004")
        assert fragmented.gene_length == 4210 - 3820

    def test_empty_file_raises_value_error(self, tmp_path):
        p = tmp_path / "empty.tsv"
        p.write_text("# just a header\n")
        with pytest.raises(ValueError):
            parse_busco_table(str(p))

    def test_score_parsed_as_float(self, sample_tsv):
        entries, _ = parse_busco_table(sample_tsv)
        scored = [e for e in entries if e.score is not None]
        assert all(isinstance(e.score, float) for e in scored)


# ── diagnostics.py tests ──────────────────────────────────────────────────────

class TestDiagnostics:

    def test_load_assembly_keys(self, sample_assembly):
        assert "contig_001" in sample_assembly
        assert "contig_004" in sample_assembly
        assert "contig_007" in sample_assembly

    def test_load_assembly_sequence_uppercase(self, sample_assembly):
        seq = sample_assembly["contig_001"]
        assert seq == seq.upper()

    def test_contig_edge_fires_near_end(self, fragmented_near_edge, sample_assembly):
        result = check_contig_edge(fragmented_near_edge, sample_assembly, edge_threshold_bp=500)
        assert result.fired is True
        assert result.confidence > 0.5

    def test_contig_edge_does_not_fire_far_from_edge(self, sample_assembly):
        # contig_001 in sample FASTA is 60 bp — use coords well inside it
        # with a very small threshold so it doesn't trigger
        entry = BuscoEntry(
            busco_id="TEST",
            status="Fragmented",
            sequence="contig_001",
            start=20,
            end=40,
            score=100.0,
            length=20,
            busco_version=5,
        )
        result = check_contig_edge(entry, sample_assembly, edge_threshold_bp=5)
        assert result.fired is False
        assert result.confidence == 0.0

    def test_contig_edge_missing_entry_does_not_fire(self, missing_entry, sample_assembly):
        result = check_contig_edge(missing_entry, sample_assembly)
        assert result.fired is False

    def test_contig_edge_unknown_contig(self, sample_assembly):
        entry = BuscoEntry(
            busco_id="TEST",
            status="Fragmented",
            sequence="contig_999",
            start=100,
            end=200,
            score=50.0,
            length=100,
            busco_version=5,
        )
        result = check_contig_edge(entry, sample_assembly)
        assert result.fired is False

    def test_assembly_gap_detects_n_stretch(self, sample_assembly):
        """contig_007 has a long N-run embedded in its sequence."""
        entry = BuscoEntry(
            busco_id="TEST_GAP",
            status="Fragmented",
            sequence="contig_007",
            start=60,
            end=120,
            score=40.0,
            length=60,
            busco_version=5,
        )
        result = check_assembly_gap(entry, sample_assembly, min_n_run=5, flank_bp=50)
        assert result.fired is True

    def test_assembly_gap_no_n_in_clean_contig(self, sample_assembly):
        entry = BuscoEntry(
            busco_id="TEST_CLEAN",
            status="Fragmented",
            sequence="contig_001",
            start=100,
            end=200,
            score=100.0,
            length=100,
            busco_version=5,
        )
        result = check_assembly_gap(entry, sample_assembly, min_n_run=10, flank_bp=50)
        assert result.fired is False


# ── aggregator.py tests ───────────────────────────────────────────────────────

class TestAggregator:

    def test_missing_no_evidence_gives_no_homology(self, missing_entry):
        result = aggregate(missing_entry, evidence=[])
        assert result.top_code == "NO_HOMOLOGY"
        assert result.top_confidence == 0.75

    def test_low_score_fires(self, fragmented_low_score):
        result = aggregate(fragmented_low_score, evidence=[])
        codes = [rc.code for rc in result.reason_codes]
        assert "LOW_SCORE" in codes

    def test_reason_codes_ranked_by_confidence(self, missing_entry, sample_assembly):
        from src.diagnostics import check_contig_edge, check_assembly_gap, EvidenceResult
        ev = EvidenceResult(
            signal_name="CONTIG_EDGE",
            fired=True,
            confidence=0.9,
            detail="Near edge."
        )
        result = aggregate(missing_entry, evidence=[ev])
        confidences = [rc.confidence for rc in result.reason_codes]
        assert confidences == sorted(confidences, reverse=True)

    def test_fallback_investigation_needed(self):
        entry = BuscoEntry(
            busco_id="TEST",
            status="Fragmented",
            sequence="contig_001",
            start=50000,
            end=51000,
            score=200.0,  # high score → no LOW_SCORE
            length=1000,
            busco_version=5,
        )
        result = aggregate(entry, evidence=[])
        assert result.top_code == "INVESTIGATION_NEEDED"

    def test_result_has_recommendation_text(self, missing_entry):
        result = aggregate(missing_entry, evidence=[])
        assert isinstance(result.recommendation, str)
        assert len(result.recommendation) > 10


# ── reporter.py tests ─────────────────────────────────────────────────────────

class TestReporter:

    def _make_results(self) -> list[DiagnosticResult]:
        from src.aggregator import ReasonCode
        r1 = DiagnosticResult(
            busco_id="EOG001",
            status="Missing",
            reason_codes=[ReasonCode("NO_HOMOLOGY", 0.75, "Check for gene loss.")],
            evidence_summary=["[NO_HOMOLOGY] fired=True"],
        )
        r2 = DiagnosticResult(
            busco_id="EOG002",
            status="Fragmented",
            reason_codes=[ReasonCode("CONTIG_EDGE", 0.85, "Consider scaffolding.")],
            evidence_summary=["[CONTIG_EDGE] fired=True"],
        )
        return [r1, r2]

    def test_tsv_written_with_correct_columns(self, tmp_path):
        results = self._make_results()
        tsv_path = str(tmp_path / "summary.tsv")
        write_tsv(results, tsv_path)

        with open(tsv_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)

        assert len(rows) == 2
        assert rows[0]["BUSCO_ID"] == "EOG001"
        assert rows[0]["Top_Reason_Code"] == "NO_HOMOLOGY"
        assert rows[1]["Top_Reason_Code"] == "CONTIG_EDGE"

    def test_tsv_confidence_is_numeric_string(self, tmp_path):
        results = self._make_results()
        tsv_path = str(tmp_path / "summary.tsv")
        write_tsv(results, tsv_path)
        with open(tsv_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                float(row["Confidence"])  # should not raise

    def test_text_report_written(self, tmp_path):
        results = self._make_results()
        report_path = str(tmp_path / "report.txt")
        write_text_report(results, report_path)
        assert os.path.isfile(report_path)
        content = open(report_path).read()
        assert "EOG001" in content
        assert "EOG002" in content
        assert "OVERVIEW" in content

    def test_text_report_contains_reason_codes(self, tmp_path):
        results = self._make_results()
        report_path = str(tmp_path / "report.txt")
        write_text_report(results, report_path)
        content = open(report_path).read()
        assert "NO_HOMOLOGY" in content
        assert "CONTIG_EDGE" in content
