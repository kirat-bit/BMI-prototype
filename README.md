# bmi-prototype

A working prototype of the **BUSCO-Missing Investigator (BMI)** pipeline —
built as part of a Google Summer of Code 2026 proposal for EMBL-EBI (Ensembl).

BMI explains *why* genes are reported as Missing or Fragmented by BUSCO,
and provides ranked, confidence-scored reason codes with actionable recommendations.

This prototype implements:
**Module 1 — BUSCO Parser** (`src/parser.py`): reads `full_table.tsv`, supports BUSCO v4 and v5
**Module 2 — Genome Diagnostics** (`src/diagnostics.py`): contig-edge proximity and N-gap detection
**Module 4 — Evidence Aggregator** (`src/aggregator.py`): reason code assignment and confidence scoring
**Module 5 — Report Generator** (`src/reporter.py`): TSV summary and plain-text report

Built directly in response to mentor feedback: "sketch a simple prototype that parses
BUSCO results and attaches preliminary diagnostic labels."*

---

## Project structure

```
bmi-prototype/
├── src/
│   ├── parser.py       # Module 1: BUSCO output parser
│   ├── diagnostics.py  # Module 2: genome-level evidence signals
│   ├── aggregator.py   # Module 4: reason code scorer
│   └── reporter.py     # Module 5: TSV + text report generator
├── tests/
│   └── test_bmi.py     # 22 pytest tests across all modules
├── data/
│   ├── full_table.tsv  # Sample BUSCO output (synthetic)
│   └── assembly.fasta  # Sample assembly FASTA (synthetic)
├── main.py             # Pipeline entry point
├── requirements.txt
└── README.md
```

---

## Setup

```bash
git clone https://github.com/kirat-bit/bmi-prototype.git
cd bmi-prototype

python -m venv venv
source venv/bin/activate       # Mac/Linux
# or: .\venv\Scripts\activate  # Windows

pip install -r requirements.txt
```

---

## Usage

```bash
# BUSCO output only (assigns reason codes from scores + coordinates)
python main.py --busco data/full_table.tsv

# With assembly FASTA (enables contig-edge and gap detection)
python main.py --busco data/full_table.tsv --assembly data/assembly.fasta

# Custom output directory
python main.py --busco data/full_table.tsv --assembly data/assembly.fasta --out results/
```

### Example terminal output

```
[BMI] Parsing BUSCO table: data/full_table.tsv
[BMI] Found 6 Missing/Fragmented BUSCOs (BUSCO v5)
[BMI] Loading assembly: data/assembly.fasta
[BMI] Loaded 5 contigs
[BMI] Running diagnostics...
[BMI] Done.
[BMI] Summary TSV → results/bmi_summary.tsv
[BMI] Text report → results/bmi_report.txt

============================================================
  BMI Quick Summary
============================================================
  BUSCO ID           Status       Top Reason             Conf
  --------------------------------------------------------
  EOG090W00D5        Missing      NO_HOMOLOGY            0.75
  EOG090W00C3        Fragmented   LOW_SCORE              0.72
  EOG090W00G8        Fragmented   CONTIG_EDGE            0.90
  ...
============================================================
```

---

## Running tests

```bash
pytest tests/ -v
```

Expected output: **22 passed**

---

## Reason code vocabulary

| Code | Level | Meaning |
|------|-------|---------|
| `CONTIG_EDGE` | Genome | BUSCO hit near scaffold boundary |
| `ASSEMBLY_GAP` | Genome | N-stretch present at locus |
| `NO_HOMOLOGY` | Genome | Missing entry with no coordinates |
| `LOW_SCORE` | Genome | BUSCO score below threshold (lineage mismatch) |
| `FRAGMENTED_SHORT` | Annotation | Hit span very short — truncated gene model |
| `INVESTIGATION_NEEDED` | — | No signal detected with available inputs |

---

## Relevance to BMI

This prototype demonstrates the core diagnostic loop of the full BMI pipeline:

```
full_table.tsv → Parser → Diagnostics → Aggregator → Report
```

The full pipeline (proposed for GSoC 2026) will extend this with:
- BLAST-based homology search (Module 2b)
- BAM coverage analysis (Module 2c)
- GFF annotation overlap checking (Module 3)
- Nextflow workflow + Docker containerisation
- HTML/Markdown report with plots

---

## Author

**Jaskirat Singh**
[GitHub](https://github.com/kirat-bit) 