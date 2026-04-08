#!/usr/bin/env python3
"""
Generate a project report for Parse Biosciences scRNA-seq analysis driven by
split-pipe, with summaries from produced outputs.

USAGE:
    python generate_report.py [OPTIONS]

OPTIONS:
    --fastq-dir DIR           Path to FASTQ directory (sublibrary subfolders)
                              Default: data/FASTQ

    --parse-dir DIR           Path to parse_comb output directory
                              Default: output/parse_comb

    --metadata FILE           Path to sample metadata CSV
                              Default: metadata/metadata.csv

    --output FILE             Output PDF path
                              Default: SparN_ParseBio_Report.pdf

    --author NAME             Report author name
                              Default: Kevin Stachelek

OUTPUTS:
    - PDF report with project information, split-pipe pipeline details,
      per-sample QC statistics, study design summary, and references.

REQUIREMENTS:
    - reportlab: for PDF generation
    - Python 3.8+
"""

import os
import glob
import csv
import json
import re
from datetime import datetime
from pathlib import Path
from reportlab.lib.pagesizes import letter, landscape
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (
    BaseDocTemplate,
    PageTemplate,
    Frame,
    NextPageTemplate,
    Table,
    TableStyle,
    Paragraph,
    Spacer,
    PageBreak,
)
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT


class ParseSublibraryExtractor:
    """Scan FASTQ directory for Parse sublibrary subfolders and file sizes."""

    def __init__(self, fastq_dir: str):
        self.fastq_dir = fastq_dir
        self.sublibraries = {}
        self._scan()

    def _scan(self):
        if not os.path.isdir(self.fastq_dir):
            return
        for subdir in sorted(os.listdir(self.fastq_dir)):
            full = os.path.join(self.fastq_dir, subdir)
            if not os.path.isdir(full):
                continue
            r1_files = glob.glob(os.path.join(full, "*_R1_*.fastq.gz"))
            r2_files = glob.glob(os.path.join(full, "*_R2_*.fastq.gz"))
            total_gb = sum(self._size_gb(f) for f in r1_files + r2_files)
            self.sublibraries[subdir] = {
                'r1_count': len(r1_files),
                'r2_count': len(r2_files),
                'total_size_gb': total_gb,
            }

    @staticmethod
    def _size_gb(path: str) -> float:
        try:
            return os.path.getsize(path) / (1024 ** 3)
        except Exception:
            return 0.0

    def get_summary(self):
        n = len(self.sublibraries)
        total_gb = sum(s['total_size_gb'] for s in self.sublibraries.values())
        return {
            'n_sublibraries': n,
            'total_size_gb': total_gb,
            'generation_date': datetime.now().strftime("%B %d, %Y"),
        }


class ParseSampleSummary:
    """Load aggregated per-sample summary from output/parse_comb/agg_sample_summary.csv."""

    # Columns we display in the report table
    DISPLAY_COLS = [
        ('number_of_cells',          'Cells'),
        ('hg38_median_tscp_per_cell','Median Transcripts/Cell'),
        ('hg38_median_genes_per_cell','Median Genes/Cell'),
        ('mean_reads_per_cell',      'Mean Reads/Cell'),
        ('sequencing_saturation',    'Seq. Saturation'),
        ('hg38_fraction_reads_in_cells', 'Frac. Reads in Cells'),
    ]

    def __init__(self, parse_dir: str):
        self.parse_dir = parse_dir
        self.rows = {}      # sample -> dict of stats
        self.all_row = {}   # the 'all-sample' combined row
        self.splitpipe_version = self._detect_version()
        self._load()

    def _detect_version(self) -> str:
        log = os.path.join(self.parse_dir, 'split-pipe_v1_6_0.log')
        if os.path.isfile(log):
            return 'v1.6.0'
        # Try to find version from any analysis_process.json
        for f in glob.glob(os.path.join(self.parse_dir, '*', 'report', 'analysis_process.json')):
            try:
                with open(f) as fh:
                    d = json.load(fh)
                ver = d.get('header', {}).get('ver_number', '')
                if ver:
                    return f'v{ver}'
            except Exception:
                pass
        return 'unknown'

    def _load(self):
        path = os.path.join(self.parse_dir, 'agg_sample_summary.csv')
        if not os.path.isfile(path):
            return
        try:
            with open(path, newline='') as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    sample = row.get('sample', '').strip()
                    if sample == 'all-sample':
                        self.all_row = row
                    else:
                        self.rows[sample] = row
        except Exception as e:
            print(f"Warning: could not read {path}: {e}")

    def get_overall_stats(self):
        row = self.all_row
        if not row:
            return {}
        return {
            'total_cells': self._fmt_int(row.get('number_of_cells', '')),
            'total_reads': self._fmt_large(row.get('number_of_reads', '')),
            'n_samples': len(self.rows),
            'seq_saturation': self._fmt_pct(row.get('sequencing_saturation', '')),
            'valid_barcode_fraction': self._fmt_pct(row.get('valid_barcode_fraction', '')),
            'median_tscp': self._fmt_num(row.get('hg38_median_tscp_per_cell', '')),
            'median_genes': self._fmt_num(row.get('hg38_median_genes_per_cell', '')),
        }

    @staticmethod
    def _fmt_int(v):
        try:
            return f"{int(float(v)):,}"
        except Exception:
            return v or 'N/A'

    @staticmethod
    def _fmt_large(v):
        try:
            n = float(v)
            if n >= 1e9:
                return f"{n/1e9:.2f}B"
            if n >= 1e6:
                return f"{n/1e6:.1f}M"
            return f"{int(n):,}"
        except Exception:
            return v or 'N/A'

    @staticmethod
    def _fmt_pct(v):
        try:
            return f"{float(v)*100:.1f}%"
        except Exception:
            return v or 'N/A'

    @staticmethod
    def _fmt_num(v):
        try:
            f = float(v)
            return f"{f:,.1f}" if f != int(f) else f"{int(f):,}"
        except Exception:
            return v or 'N/A'

    def table_rows(self):
        """Return (header_row, data_rows) for the per-sample stats table."""
        header = ['Sample'] + [label for _, label in self.DISPLAY_COLS]
        data = []
        for sample in sorted(self.rows.keys()):
            row = self.rows[sample]
            cells = [sample]
            for col, _ in self.DISPLAY_COLS:
                val = row.get(col, '')
                if col == 'sequencing_saturation' or col == 'hg38_fraction_reads_in_cells':
                    cells.append(self._fmt_pct(val))
                elif col in ('number_of_cells', 'hg38_median_tscp_per_cell',
                             'hg38_median_genes_per_cell'):
                    cells.append(self._fmt_int(val))
                elif col == 'mean_reads_per_cell':
                    cells.append(self._fmt_large(val))
                else:
                    cells.append(val or 'N/A')
            data.append(cells)
        return header, data


class ParseMetadata:
    """Load project metadata CSV with experimental design columns."""

    def __init__(self, path: str):
        self.path = path
        self.rows = []
        self._load()

    def _load(self):
        if not self.path or not os.path.isfile(self.path):
            return
        try:
            with open(self.path, newline='') as fh:
                reader = csv.DictReader(fh)
                self.rows = [row for row in reader]
        except Exception as e:
            print(f"Warning: could not read metadata {self.path}: {e}")

    def get_experiment_groups(self):
        """Return dict of experiment_id -> list of sample rows."""
        groups = {}
        for row in self.rows:
            exp = row.get('experiment', '').strip()
            groups.setdefault(exp, []).append(row)
        return groups


class ReportGenerator:
    """Generate PDF report summarizing Parse scRNA-seq pipeline inputs and outputs."""

    def __init__(self, output_path, author, fastq_dir, parse_dir, metadata_path, workdir='.'):
        self.output_path = output_path
        self.author = author
        self.fastq_dir = fastq_dir
        self.parse_dir = parse_dir
        self.workdir = workdir
        self.sublibs = ParseSublibraryExtractor(fastq_dir)
        self.summary = self.sublibs.get_summary()
        self.parse_summary = ParseSampleSummary(parse_dir)
        self.metadata = ParseMetadata(metadata_path)

    def generate(self):
        left_margin = 0.75 * inch
        right_margin = 0.75 * inch
        top_margin = 0.75 * inch
        bottom_margin = 0.75 * inch

        pw, ph = letter
        lw, lh = landscape(letter)

        portrait_frame = Frame(
            left_margin, bottom_margin,
            pw - left_margin - right_margin,
            ph - top_margin - bottom_margin,
            id='portrait_frame'
        )
        landscape_frame = Frame(
            left_margin, bottom_margin,
            lw - left_margin - right_margin,
            lh - top_margin - bottom_margin,
            id='landscape_frame'
        )

        def on_portrait(canvas, doc):
            canvas.setPageSize(letter)

        def on_landscape(canvas, doc):
            canvas.setPageSize(landscape(letter))

        portrait_template = PageTemplate(id='Portrait', frames=[portrait_frame], onPage=on_portrait)
        landscape_template = PageTemplate(id='Landscape', frames=[landscape_frame], onPage=on_landscape)

        doc = BaseDocTemplate(
            self.output_path,
            pagesize=letter,
            rightMargin=right_margin,
            leftMargin=left_margin,
            topMargin=top_margin,
            bottomMargin=bottom_margin,
            pageTemplates=[portrait_template, landscape_template],
        )

        elements = []
        styles = getSampleStyleSheet()

        title_style = ParagraphStyle(
            'CustomTitle',
            parent=styles['Heading1'],
            fontSize=24,
            textColor=colors.HexColor('#1f4788'),
            spaceAfter=12,
            alignment=TA_CENTER,
            fontName='Helvetica-Bold'
        )
        heading_style = ParagraphStyle(
            'CustomHeading',
            parent=styles['Heading2'],
            fontSize=14,
            textColor=colors.HexColor('#1f4788'),
            spaceAfter=12,
            spaceBefore=12,
            fontName='Helvetica-Bold'
        )
        body_style = ParagraphStyle(
            'CustomBody',
            parent=styles['BodyText'],
            fontSize=10,
            alignment=TA_LEFT,
            spaceAfter=10,
            leading=14
        )
        cell_style_small = ParagraphStyle(
            'CellSmall',
            parent=styles['BodyText'],
            fontSize=8,
            leading=10,
            wordWrap='CJK'
        )

        # ── Title ───────────────────────────────────────────────────────────
        elements.append(Paragraph("Parse Biosciences scRNA-seq Project Report", title_style))
        elements.append(Paragraph("SparN Toxicant Exposure Study", ParagraphStyle(
            'Subtitle', parent=styles['BodyText'], fontSize=14,
            alignment=TA_CENTER, textColor=colors.HexColor('#555555'), spaceAfter=16
        )))
        elements.append(Spacer(1, 0.2 * inch))

        # ── Project Information ──────────────────────────────────────────────
        elements.append(Paragraph("Project Information", heading_style))
        overall = self.parse_summary.get_overall_stats()
        project_info = [
            ['Generation Date:', self.summary['generation_date']],
            ['Author:', self.author],
            ['Pipeline:', f"Parse Biosciences split-pipe {self.parse_summary.splitpipe_version}"],
            ['Reference Genome:', 'hg38 (GRCh38)'],
            ['Sublibraries:', str(self.summary['n_sublibraries'])],
            ['Total Input Data:', f"{self.summary['total_size_gb']:.1f} GB"],
        ]

        info_table = Table(project_info, colWidths=[2.7 * inch, 3.3 * inch])
        info_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (0, -1), colors.HexColor('#E8EEF7')),
            ('TEXTCOLOR', (0, 0), (-1, -1), colors.black),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
            ('TOPPADDING', (0, 0), (-1, -1), 8),
            ('GRID', (0, 0), (-1, -1), 1, colors.grey),
        ]))
        elements.append(info_table)
        elements.append(Spacer(1, 0.3 * inch))

        # ── Overall QC Summary ───────────────────────────────────────────────
        elements.append(Paragraph("Overall QC Summary", heading_style))
        summary_data = [
            ['Metric', 'Value'],
            ['Total Samples Processed', str(overall.get('n_samples', 'N/A'))],
            ['Total Cells Recovered', overall.get('total_cells', 'N/A')],
            ['Total Reads', overall.get('total_reads', 'N/A')],
            ['Median Transcripts per Cell', overall.get('median_tscp', 'N/A')],
            ['Median Genes per Cell', overall.get('median_genes', 'N/A')],
            ['Sequencing Saturation', overall.get('seq_saturation', 'N/A')],
            ['Valid Barcode Fraction', overall.get('valid_barcode_fraction', 'N/A')],
        ]

        summary_table = Table(summary_data, colWidths=[3.0 * inch, 2.5 * inch])
        summary_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1f4788')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
            ('TOPPADDING', (0, 0), (-1, -1), 8),
            ('GRID', (0, 0), (-1, -1), 1, colors.grey),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F0F0F0')]),
        ]))
        elements.append(summary_table)
        elements.append(Spacer(1, 0.3 * inch))

        # ── Pipeline Overview ────────────────────────────────────────────────
        elements.append(Paragraph("Pipeline Overview", heading_style))

        elements.append(Paragraph(
            f"<b>Parse Biosciences split-pipe {self.parse_summary.splitpipe_version}</b>",
            styles['Heading3']
        ))
        elements.append(Paragraph(
            "Raw FASTQ files were processed with the Parse Biosciences split-pipe pipeline. "
            "Split-pool ligation-based transcriptomics (SPLiT-seq) uses three rounds of combinatorial "
            "barcoding to label individual cells without requiring microfluidic capture. "
            "The pipeline demultiplexes reads by sublibrary, corrects barcodes using a whitelist "
            "(allowing up to one edit distance), maps reads to the reference genome using STAR, "
            "assigns reads to cells based on the UMI knee-point cutoff, and produces a digital "
            "gene expression (DGE) matrix per sample.",
            body_style
        ))

        elements.append(Paragraph("<b>Reference Genome and Annotation</b>", styles['Heading3']))
        elements.append(Paragraph(
            "Reads were aligned to the human genome (hg38/GRCh38) with the corresponding "
            "GENCODE annotation. Transcript mapping fraction and exonic fraction are reported "
            "per sublibrary.",
            body_style
        ))

        elements.append(Paragraph("<b>Cell Calling</b>", styles['Heading3']))
        elements.append(Paragraph(
            "Cells were called using a UMI knee-point algorithm applied to the transcript count "
            "distribution. The cell_tscp_cutoff represents the minimum transcript count threshold "
            "used to distinguish cells from empty droplets/beads.",
            body_style
        ))

        elements.append(Paragraph("<b>Downstream Analysis</b>", styles['Heading3']))
        elements.append(Paragraph(
            "Filtered DGE matrices were analyzed with Scanpy and scVI-tools. Integration across "
            "sublibraries and samples was performed with scVI/scANVI. Clustering used the Leiden "
            "algorithm. Additional integration was performed with Harmony for comparison.",
            body_style
        ))

        # ── Study Design ─────────────────────────────────────────────────────
        elements.append(Paragraph("Study Design", heading_style))
        elements.append(Paragraph(
            "Human embryonic stem cells (hESCs) were differentiated toward an osteoblast lineage "
            "and exposed to toxicant compounds at multiple concentrations and timepoints. "
            "Two independent experiments were conducted:",
            body_style
        ))

        if self.metadata.rows:
            exp_groups = self.metadata.get_experiment_groups()
            for exp_id, rows in sorted(exp_groups.items()):
                treatments = sorted({r.get('treatment', '') for r in rows if r.get('treatment')})
                days = sorted({r.get('day', '') for r in rows if r.get('day')}, key=lambda x: int(x) if x.isdigit() else 0)
                cell_types = sorted({r.get('cell_type', '') for r in rows if r.get('cell_type')})
                exp_label = "Experiment 1 (Bisphenols)" if exp_id == '1' else f"Experiment {exp_id} (Zyn Nicotine Pouches)" if exp_id == '2' else f"Experiment {exp_id}"
                elements.append(Paragraph(f"<b>{exp_label}</b>", styles['Heading3']))
                elements.append(Paragraph(
                    f"Treatments: {', '.join(treatments)} | "
                    f"Days: {', '.join(days)} | "
                    f"Cell Types: {', '.join(cell_types)}",
                    body_style
                ))

        # ── Sample Details (landscape) ───────────────────────────────────────
        elements.append(NextPageTemplate('Landscape'))
        elements.append(PageBreak())
        elements.append(Paragraph("Per-Sample QC Statistics", heading_style))
        elements.append(Paragraph(
            "Key QC metrics from the split-pipe combined analysis for each experimental sample. "
            "Values are from the combined (all-sublibrary) run per sample.",
            body_style
        ))
        elements.append(Spacer(1, 0.1 * inch))

        header, data_rows = self.parse_summary.table_rows()
        table_data = [header] + data_rows

        # Column widths: sample name wider, rest equal
        n_data_cols = len(header) - 1
        name_w = 1.6 * inch
        data_w = (9.5 * inch - name_w) / n_data_cols
        col_widths = [name_w] + [data_w] * n_data_cols

        stats_table = Table(table_data, colWidths=col_widths, repeatRows=1)
        stats_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1f4788')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('ALIGN', (0, 0), (0, -1), 'LEFT'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 8),
            ('FONTSIZE', (0, 1), (-1, -1), 7),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 5),
            ('TOPPADDING', (0, 0), (-1, -1), 5),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F5F5F5')]),
        ]))
        elements.append(stats_table)
        elements.append(Spacer(1, 0.2 * inch))

        # ── Sample Metadata Table ────────────────────────────────────────────
        elements.append(Paragraph("Sample Descriptions", heading_style))

        if self.metadata.rows:
            meta_header = ['#', 'Sample', 'Experiment', 'Day', 'Treatment', 'Concentration', 'Cell Type']
            meta_data = [meta_header]
            for row in self.metadata.rows:
                meta_data.append([
                    row.get('sample_number', ''),
                    row.get('sample_title', ''),
                    row.get('experiment', ''),
                    row.get('day', ''),
                    row.get('treatment', ''),
                    row.get('treatment_level', ''),
                    row.get('cell_type', ''),
                ])
            meta_col_widths = [0.4*inch, 1.4*inch, 1.0*inch, 0.5*inch, 1.0*inch, 1.7*inch, 1.5*inch]
            meta_table = Table(meta_data, colWidths=meta_col_widths, repeatRows=1)
            meta_table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1f4788')),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('ALIGN', (1, 1), (1, -1), 'LEFT'),
                ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, 0), 8),
                ('FONTSIZE', (0, 1), (-1, -1), 7),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
                ('TOPPADDING', (0, 0), (-1, -1), 4),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
                ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F5F5F5')]),
            ]))
            elements.append(meta_table)

        # ── References (portrait) ────────────────────────────────────────────
        elements.append(NextPageTemplate('Portrait'))
        elements.append(PageBreak())
        elements.append(Paragraph("References", heading_style))

        reference_style = ParagraphStyle(
            'Reference',
            parent=styles['BodyText'],
            fontSize=9,
            leftIndent=0.2 * inch,
            spaceAfter=8,
            leading=11,
            textColor=colors.black
        )

        references = [
            "Rosenberg et al. (2018). Single-cell profiling of the developing mouse brain and "
            "spinal cord with split-pool barcoding. Science 360(6385):176-182.",
            "Parse Biosciences. (2023). split-pipe: Single-cell RNA sequencing analysis pipeline. "
            "https://www.parsebiosciences.com",
            "Wolf et al. (2018). SCANPY: large-scale single-cell gene expression data analysis. "
            "Genome Biology 19:15.",
            "Lopez et al. (2018). Deep generative modeling for single-cell transcriptomics. "
            "Nature Methods 15:1053-1058.",
            "Traag et al. (2019). From Louvain to Leiden: guaranteeing well-connected communities. "
            "Scientific Reports 9:5233.",
            "Korsunsky et al. (2019). Fast, sensitive and accurate integration of single-cell data "
            "with Harmony. Nature Methods 16:1289-1296.",
            "Dobin et al. (2013). STAR: ultrafast universal RNA-seq aligner. "
            "Bioinformatics 29(1):15-21.",
        ]

        for i, ref in enumerate(references, 1):
            elements.append(Paragraph(f"<b>{i}.</b> {ref}", reference_style))

        doc.build(elements)
        print(f"Report generated successfully: {self.output_path}")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Generate a Parse Biosciences scRNA-seq project report'
    )
    parser.add_argument(
        '--fastq-dir',
        default='data/FASTQ',
        help='Path to FASTQ directory with sublibrary subfolders (default: data/FASTQ)'
    )
    parser.add_argument(
        '--parse-dir',
        default='output/parse_comb',
        help='Path to parse_comb output directory (default: output/parse_comb)'
    )
    parser.add_argument(
        '--metadata',
        default='metadata/metadata.csv',
        help='Path to sample metadata CSV (default: metadata/metadata.csv)'
    )
    parser.add_argument(
        '--output',
        default='SparN_ParseBio_Report.pdf',
        help='Output PDF path (default: SparN_ParseBio_Report.pdf)'
    )
    parser.add_argument(
        '--author',
        default='Kevin Stachelek',
        help='Report author name'
    )

    args = parser.parse_args()

    generator = ReportGenerator(
        output_path=args.output,
        author=args.author,
        fastq_dir=args.fastq_dir,
        parse_dir=args.parse_dir,
        metadata_path=args.metadata,
        workdir='.',
    )
    generator.generate()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
