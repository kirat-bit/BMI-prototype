"""
parser.py — Module 1: BUSCO Output Parser

Reads a BUSCO full_table.tsv file and returns a list of BuscoEntry
objects for all Missing and Fragmented genes.

Supports BUSCO v4 and v5 output formats.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional
import re


# ── Data model ──────

@dataclass
class BuscoEntry:
    """One row from full_table.tsv representing a single BUSCO gene."""
    busco_id: str
    status: str                   
    sequence: Optional[str]       
    start: Optional[int]          
    end: Optional[int]            
    score: Optional[float]        
    length: Optional[int]        
    busco_version: int            

    @property
    def has_coordinates(self) -> bool:
        """True if this entry has sequence coordinates (Fragmented entries do)."""
        return self.sequence is not None

    @property
    def gene_length(self) -> Optional[int]:
        """Genomic span of the hit in base pairs."""
        if self.start is not None and self.end is not None:
            return self.end - self.start
        return None


# ── Parser ─────────

def detect_busco_version(header_lines: list[str]) -> int:
    """
    Detect BUSCO version from comment lines at the top of full_table.tsv.
    Returns 5 or 4. Defaults to 5 if version cannot be determined.
    """
    for line in header_lines:
        match = re.search(r"BUSCO version is:\s*(\d+)", line, re.IGNORECASE)
        if match:
            major = int(match.group(1))
            return major if major in (4, 5) else 5
    return 5  


def parse_row(row: list[str], version: int) -> Optional[BuscoEntry]:
    """
    Parse one data row from full_table.tsv into a BuscoEntry.

    BUSCO v5 columns: busco_id, status, sequence, start, end, strand, score, length
    BUSCO v4 columns: busco_id, status, sequence, start, end, score, length
    Returns None if the row is not Missing or Fragmented.
    """
    if len(row) < 7:
        return None

    status = row[1].strip()
    if status not in ("Missing", "Fragmented"):
        return None

    busco_id = row[0].strip()

    # Coordinates are N/A for Missing entries
    def parse_optional_int(value: str) -> Optional[int]:
        return int(value) if value not in ("N/A", "", ".") else None

    def parse_optional_float(value: str) -> Optional[float]:
        return float(value) if value not in ("N/A", "", ".") else None

    sequence = row[2].strip() if row[2].strip() not in ("N/A", "") else None

    if version == 5:
        # v5: busco_id, status, sequence, start, end, strand, score, length
        start  = parse_optional_int(row[3])
        end    = parse_optional_int(row[4])
        # row[5] is strand — skip
        score  = parse_optional_float(row[6]) if len(row) > 6 else None
        length = parse_optional_int(row[7])   if len(row) > 7 else None
    else:
        # v4: busco_id, status, sequence, start, end, score, length
        start  = parse_optional_int(row[3])
        end    = parse_optional_int(row[4])
        score  = parse_optional_float(row[5]) if len(row) > 5 else None
        length = parse_optional_int(row[6])   if len(row) > 6 else None

    return BuscoEntry(
        busco_id=busco_id,
        status=status,
        sequence=sequence,
        start=start,
        end=end,
        score=score,
        length=length,
        busco_version=version,
    )


def parse_busco_table(filepath: str) -> tuple[list[BuscoEntry], int]:
    """
    Read a BUSCO full_table.tsv and return:
      - list of BuscoEntry for all Missing and Fragmented genes
      - detected BUSCO version (4 or 5)

    Args:
        filepath: path to full_table.tsv

    Returns:
        (entries, version)

    Raises:
        FileNotFoundError: if the file does not exist
        ValueError: if the file appears to be malformed
    """
    header_lines: list[str] = []
    entries: list[BuscoEntry] = []
    row_count = 0

    with open(filepath, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")

            # Collect comment/header lines for version detection
            if line.startswith("#"):
                header_lines.append(line)
                continue

            # Skip blank lines
            if not line.strip():
                continue

            row = line.split("\t")
            row_count += 1

            # Detect version once we have headers, before processing rows
            if row_count == 1:
                version = detect_busco_version(header_lines)

            entry = parse_row(row, version)
            if entry is not None:
                entries.append(entry)

    if row_count == 0:
        raise ValueError(f"No data rows found in {filepath}. Is this a valid full_table.tsv?")

    return entries, version
