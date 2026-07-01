#!/usr/bin/env python3
"""
Prepare GEO submission metadata for Parse Biosciences scRNA-seq experiment.

Reads project metadata, FASTQ inventory, and processed output paths, then
writes a filled GEO metadata Excel file based on the GEO template structure.

Usage:
    python src/prep_geo_submission.py [--output OUTPUT]
"""
import argparse
import csv
import os
import re
import sys
from pathlib import Path

import openpyxl
from openpyxl.styles import Font, PatternFill

PROJECT_DIR = Path(__file__).resolve().parent.parent
METADATA_CSV = PROJECT_DIR / "metadata" / "metadata.csv"
FASTQ_DIR = PROJECT_DIR / "data" / "FASTQ"
CHECKSUMS_FILE = PROJECT_DIR / "data" / "checksums.md5"
PROCESSED_DIR = PROJECT_DIR / "output" / "parse_comb"
TEMPLATE_PATH = PROJECT_DIR / "docs" / "parse_template.xlsx"
DEFAULT_OUTPUT = PROJECT_DIR / "output" / "geo_submission_metadata.xlsx"

# ---- Experiment-level constants ----
ORGANISM = "Homo sapiens"
LIBRARY_STRATEGY = "scRNA-seq"
MOLECULE = "polyA RNA"
SINGLE_OR_PAIRED = "paired-end"
INSTRUMENT_MODEL = "Illumina NovaSeq X Plus"
GENOME_BUILD = "GRCh38 (hg38)"

EXTRACT_PROTOCOL = (
    "Human embryonic stem cells (hESCs) were differentiated toward an osteogenic "
    "lineage for 7 or 20 days. At the designated time point, cells were exposed to "
    "the indicated toxicant (bisphenol analogue or Zyn nicotine pouch extract) at the "
    "indicated concentration (IC10, IC25, or IC50) or left untreated. Cells were then "
    "fixed and processed following the Parse Biosciences Whole Transcriptome Kit (WTK) "
    "fixation protocol."
)

LIBRARY_CONSTRUCTION_PROTOCOL = (
    "Single-cell RNA-seq libraries were prepared using the Parse Biosciences Whole "
    "Transcriptome Kit (WT Mega, v3 chemistry) following the manufacturer's protocol. "
    "Fixed cells from all 30 samples were combinatorially barcoded and pooled across "
    "8 sublibraries for sequencing."
)

DATA_PROCESSING_STEPS = [
    (
        "Raw FASTQ files from each sublibrary were processed independently using "
        "split-pipe v1.6.0 (Parse Biosciences) in '--mode all' with the GRCh38 (hg38) "
        "reference transcriptome."
    ),
    (
        "Per-sublibrary outputs were combined using split-pipe v1.6.0 in '--mode comb' "
        "to generate per-sample demultiplexed DGE matrices."
    ),
    (
        "Processed data files per sample: count_matrix.mtx (sparse UMI count matrix in "
        "Matrix Market format), all_genes.csv (gene ID and gene name annotations), "
        "cell_metadata.csv (per-cell barcode and QC metadata)."
    ),
]

PROCESSED_FILES_FORMAT = (
    "count_matrix.mtx: sparse UMI count matrix (Matrix Market format, cells × genes); "
    "all_genes.csv: gene annotations (Ensembl ID, gene name); "
    "cell_metadata.csv: per-cell metadata including barcode, sample assignment, and QC metrics."
)

STUDY_TITLE_PLACEHOLDER = (
    "Single-cell transcriptomic profiling of bisphenol analogue and nicotine "
    "pouch extract effects on human osteogenic differentiation"
)

STUDY_SUMMARY_PLACEHOLDER = (
    "[TODO: Insert abstract here. Describe the study goal, experimental system "
    "(hESC osteogenic differentiation), toxicants tested (BPA, BPF, BPS; Zyn pouches), "
    "and key findings.]"
)

EXPERIMENTAL_DESIGN_PLACEHOLDER = (
    "Human embryonic stem cells were differentiated toward osteoblasts over 20 days. "
    "Cells were treated with bisphenol analogues (BPA, BPF, BPS) at IC10/IC25/IC50 "
    "concentrations at day 7 or day 20, or with Zyn nicotine pouch extracts "
    "(Chill, Wintergreen, Smooth, Original; 3 mg or 6 mg) at day 7. "
    "Untreated controls were collected at day 0, 7, and 20. Single-cell RNA-seq "
    "was performed using Parse Biosciences combinatorial barcoding across 8 "
    "sublibraries (30 samples, ~100,000 cells total)."
)


def load_metadata(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def load_checksums(path: Path) -> dict[str, str]:
    """Return {filename_basename: md5} from a standard md5sum file."""
    checksums: dict[str, str] = {}
    if not path.exists():
        return checksums
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            md5, relpath = line.split(None, 1)
            checksums[Path(relpath).name] = md5
    return checksums


def collect_fastq_runs(fastq_dir: Path) -> list[dict]:
    """
    Return list of dicts with keys: sublibrary, lane, R1, R2, I1, I2.
    Derived from the FASTQ directory structure.
    """
    runs: list[dict] = []
    for sub_dir in sorted(fastq_dir.iterdir()):
        if not sub_dir.is_dir():
            continue
        sublibrary = sub_dir.name
        lane_files: dict[str, dict[str, str]] = {}
        for fq in sorted(sub_dir.glob("*.fastq.gz")):
            m = re.search(r"_(L\d+)_(I1|I2|R1|R2)_\d+\.fastq\.gz$", fq.name)
            if not m:
                continue
            lane, read = m.group(1), m.group(2)
            lane_files.setdefault(lane, {})[read] = fq.name
        for lane, files in sorted(lane_files.items()):
            # Only emit rows that have at least R1 and R2
            if "R1" not in files or "R2" not in files:
                continue
            runs.append(
                {
                    "sublibrary": sublibrary,
                    "lane": lane,
                    "R1": files.get("R1", ""),
                    "R2": files.get("R2", ""),
                    "I1": files.get("I1", ""),
                    "I2": files.get("I2", ""),
                }
            )
    return runs


def collect_processed_files(processed_dir: Path, sample_title: str) -> list[str]:
    """Return basenames of processed DGE files for a sample."""
    dge_dir = processed_dir / sample_title / "DGE_filtered"
    if not dge_dir.exists():
        return []
    files = ["count_matrix.mtx", "all_genes.csv", "cell_metadata.csv"]
    return [f for f in files if (dge_dir / f).exists()]


def cell_type_to_tissue(cell_type: str) -> str:
    mapping = {
        "hESC": "embryoid body",
        "Osteoprogenitor": "bone marrow",
        "Osteoblast": "bone",
    }
    return mapping.get(cell_type, "hESC-derived osteogenic lineage")


def build_sample_title_display(row: dict) -> str:
    """Human-readable title for a GEO sample record."""
    parts = [row["sample_title"]]
    if row["treatment"] != "UNT":
        parts.append(f"{row['treatment']} {row['treatment_level']}")
    parts.append(f"day {row['day']}")
    return ", ".join(parts)


def build_geo_workbook(metadata: list[dict], fastq_runs: list[dict]) -> openpyxl.Workbook:
    wb = openpyxl.Workbook()
    wb.remove(wb.active)  # remove default sheet

    _write_study_sheet(wb)
    _write_samples_sheet(wb, metadata)
    _write_protocols_sheet(wb)
    _write_paired_end_sheet(wb, fastq_runs)
    _write_md5_sheet(wb, fastq_runs)

    return wb


def _section_header(ws, row: int, label: str):
    cell = ws.cell(row=row, column=1, value=label)
    cell.font = Font(bold=True, size=12)
    cell.fill = PatternFill("solid", fgColor="DDEBF7")


def _comment_row(ws, row: int, text: str):
    cell = ws.cell(row=row, column=1, value=f"# {text}")
    cell.font = Font(italic=True, color="595959")


def _write_study_sheet(wb: openpyxl.Workbook):
    ws = wb.create_sheet("STUDY")
    ws.column_dimensions["A"].width = 30
    ws.column_dimensions["B"].width = 100

    rows = [
        ("*title", STUDY_TITLE_PLACEHOLDER),
        ("*summary (abstract)", STUDY_SUMMARY_PLACEHOLDER),
        ("*experimental design", EXPERIMENTAL_DESIGN_PLACEHOLDER),
        ("contributor", "Sparks, Nicole"),
        ("contributor", "[TODO: Add additional contributors as Lastname, Firstname]"),
        ("supplementary file", ""),
    ]
    for i, (label, value) in enumerate(rows, start=1):
        ws.cell(row=i, column=1, value=label)
        ws.cell(row=i, column=2, value=value)


def _write_samples_sheet(wb: openpyxl.Workbook, metadata: list[dict]):
    ws = wb.create_sheet("SAMPLES")

    header = [
        "*library name",
        "*title",
        "*library strategy",
        "*organism",
        "**tissue",
        "**cell line",
        "**cell type",
        "treatment",
        "treatment level",
        "day",
        "experiment",
        "*molecule",
        "*single or paired-end",
        "*instrument model",
        "description",
        "raw file 1",
        "processed data file 1",
        "processed data file 2",
        "processed data file 3",
    ]
    ws.append(header)
    for cell in ws[1]:
        cell.font = Font(bold=True)
        cell.fill = PatternFill("solid", fgColor="BDD7EE")

    # All samples share the same sublibraries as raw files. Use the first
    # sublibrary (L001 of Sublibrary1) as the representative raw file name.
    # GEO reviewers understand that for combinatorial barcoding, all sublibraries
    # contribute to each biological sample; full file listings are in PAIRED-END.
    representative_raw = "See PAIRED-END EXPERIMENTS sheet; all 8 sublibraries shared"

    for row in metadata:
        title = row["sample_title"]
        processed = collect_processed_files(PROCESSED_DIR, title)
        # Pad to 3 entries (count_matrix.mtx, all_genes.csv, cell_metadata.csv)
        while len(processed) < 3:
            processed.append("")

        ws.append([
            title,                                      # *library name
            build_sample_title_display(row),            # *title
            LIBRARY_STRATEGY,                           # *library strategy
            ORGANISM,                                   # *organism
            cell_type_to_tissue(row["cell_type"]),      # **tissue
            "hESC-derived",                             # **cell line
            row["cell_type"],                           # **cell type
            row["treatment"],                           # treatment
            row["treatment_level"],                     # treatment level
            row["day"],                                 # day
            row["experiment"],                          # experiment
            MOLECULE,                                   # *molecule
            SINGLE_OR_PAIRED,                           # *single or paired-end
            INSTRUMENT_MODEL,                           # *instrument model
            row["experiment_explanation"],              # description
            representative_raw,                         # raw file 1
            processed[0],                               # processed data file 1
            processed[1],                               # processed data file 2
            processed[2],                               # processed data file 3
        ])

    # Auto-width columns
    for col in ws.columns:
        max_len = max((len(str(c.value or "")) for c in col), default=0)
        ws.column_dimensions[col[0].column_letter].width = min(max_len + 2, 60)


def _write_protocols_sheet(wb: openpyxl.Workbook):
    ws = wb.create_sheet("PROTOCOLS")
    ws.column_dimensions["A"].width = 35
    ws.column_dimensions["B"].width = 120

    rows = [
        ("growth protocol", ""),
        ("treatment protocol", (
            "Cells were treated with bisphenol analogues (BPA, BPF, BPS) at IC10, "
            "IC25, or IC50 concentrations, or with Zyn nicotine pouch extracts "
            "(Chill, Wintergreen, Smooth, Original) at 3 mg or 6 mg equivalents, "
            "for 24 hours prior to harvest. Untreated controls received vehicle only."
        )),
        ("*extract protocol", EXTRACT_PROTOCOL),
        ("*library construction protocol", LIBRARY_CONSTRUCTION_PROTOCOL),
    ]
    for label, value in rows:
        ws.append([label, value])

    ws.append(["*data processing step", DATA_PROCESSING_STEPS[0]])
    for step in DATA_PROCESSING_STEPS[1:]:
        ws.append(["data processing step", step])

    ws.append(["*genome build/assembly", GENOME_BUILD])
    ws.append(["*processed data files format and content", PROCESSED_FILES_FORMAT])

    for cell in ws["A"]:
        cell.font = Font(bold=True)


def _write_paired_end_sheet(wb: openpyxl.Workbook, fastq_runs: list[dict]):
    ws = wb.create_sheet("PAIRED-END EXPERIMENTS")

    header = ["file name 1 (R1)", "file name 2 (R2)", "file name 3 (I1)", "file name 4 (I2)", "sublibrary", "lane"]
    ws.append(header)
    for cell in ws[1]:
        cell.font = Font(bold=True)
        cell.fill = PatternFill("solid", fgColor="BDD7EE")

    for run in fastq_runs:
        ws.append([run["R1"], run["R2"], run["I1"], run["I2"], run["sublibrary"], run["lane"]])

    for col in ws.columns:
        max_len = max((len(str(c.value or "")) for c in col), default=0)
        ws.column_dimensions[col[0].column_letter].width = min(max_len + 2, 80)


def _write_md5_sheet(wb: openpyxl.Workbook, fastq_runs: list[dict]):
    ws = wb.create_sheet("MD5 Checksums")
    checksums = load_checksums(CHECKSUMS_FILE)

    header = ["file name", "MD5 checksum"]
    ws.append(header)
    for cell in ws[1]:
        cell.font = Font(bold=True)
        cell.fill = PatternFill("solid", fgColor="BDD7EE")

    for run in fastq_runs:
        for fname in [run["R1"], run["R2"], run["I1"], run["I2"]]:
            if fname:
                ws.append([fname, checksums.get(fname, "[TODO: compute]")])

    ws.column_dimensions["A"].width = 70
    ws.column_dimensions["B"].width = 40


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=DEFAULT_OUTPUT,
        help="Output Excel file path (default: output/geo_submission_metadata.xlsx)"
    )
    args = parser.parse_args()

    if not METADATA_CSV.exists():
        sys.exit(f"ERROR: metadata not found at {METADATA_CSV}")
    if not FASTQ_DIR.exists():
        sys.exit(f"ERROR: FASTQ directory not found at {FASTQ_DIR}")

    print(f"Loading metadata from {METADATA_CSV}")
    metadata = load_metadata(METADATA_CSV)
    print(f"  {len(metadata)} samples")

    print(f"Scanning FASTQ directory {FASTQ_DIR}")
    fastq_runs = collect_fastq_runs(FASTQ_DIR)
    print(f"  {len(fastq_runs)} FASTQ runs (sublibrary × lane)")

    print("Building GEO workbook...")
    wb = build_geo_workbook(metadata, fastq_runs)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    wb.save(args.output)
    print(f"Saved to {args.output}")

    # Summary of what needs manual review
    print("\n--- Fields requiring manual review ---")
    print("  STUDY sheet:")
    print("    *summary (abstract)  — insert final abstract text")
    print("    contributor          — add all co-authors")
    print("    supplementary file   — add path to combined AnnData/Seurat object if any")
    print("  SAMPLES sheet:")
    print("    raw file 1           — currently placeholder; update once FASTQ files are")
    print("                           registered in GEO FTP or note multiplexed structure")
    print("  PROTOCOLS sheet:")
    print("    growth protocol      — add hESC culture and differentiation protocol details")


if __name__ == "__main__":
    main()
