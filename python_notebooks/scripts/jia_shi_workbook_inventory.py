#!/usr/bin/env python3
"""Inventory accessible Jia/Shi Excel workbooks and their encoded labels.

This is an accounting helper, not an overlap analysis. It answers:

- Which workbooks are physically accessible?
- Which sheets exist in each workbook?
- Which sheets have recognizable `gene`/`cluster` marker-table columns?
- Which cluster labels are encoded in those tables?
- Do curated crosswalk terms such as M2/M3/M4/M5/M6/M7 appear as cluster
  labels, gene symbols, sheet names, or anywhere in workbook text?
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Iterable

import openpyxl
import pandas as pd


PROJECT_ROOT_DEFAULT = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
DEFAULT_JIA_XLSX = "reference/science.adw1803_data_s9.xlsx"
DEFAULT_SHI_DIR = "reference/shi_2021_tables_s2_to_s9"
DEFAULT_RUN_LABEL = "jia_shi_workbook_inventory_v1"

CURATED_TERMS = [
    "M2",
    "M3",
    "M4",
    "M5",
    "M6",
    "M7",
    "pM1",
    "pM2",
    "pM3",
    "pM4",
    "pL1",
    "pL2",
    "pL3",
    "pC1",
    "pC2",
    "pC3",
    "LHX8",
    "ISL1",
    "NR2F1",
    "NR2F2",
    "EPHA5",
    "MEF2C",
    "LHX6",
    "NFIA",
    "CRABP1",
    "ANGPT2",
    "ETV1",
    "GBX2",
    "ZIC1",
    "ZFHX3",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=os.environ.get("PROJECT_ROOT", PROJECT_ROOT_DEFAULT))
    parser.add_argument("--jia-xlsx", default=DEFAULT_JIA_XLSX)
    parser.add_argument("--shi-dir", default=DEFAULT_SHI_DIR)
    parser.add_argument("--run-label", default=DEFAULT_RUN_LABEL)
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def resolve_path(path: str | Path, project_root: Path) -> Path:
    candidate = Path(path).expanduser()
    if not candidate.is_absolute():
        candidate = project_root / candidate
    return candidate.resolve()


def prepare_outdir(args: argparse.Namespace, project_root: Path) -> tuple[Path, Path, Path]:
    outdir = (
        Path(args.outdir).expanduser().resolve()
        if args.outdir
        else project_root / "results" / "jia_s9_shi_lineage_overlap" / args.run_label
    )
    if outdir.exists() and any(outdir.iterdir()):
        if not args.overwrite:
            raise FileExistsError(f"Output directory is not empty; use --overwrite: {outdir}")
        for path in outdir.iterdir():
            if path.is_file() or path.is_symlink():
                path.unlink()
            else:
                shutil.rmtree(path)
    table_dir = outdir / "tables"
    report_dir = outdir / "reports"
    table_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    return outdir, table_dir, report_dir


def normalize_token(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def detect_header_row(raw: pd.DataFrame) -> int | None:
    for idx, row in raw.iterrows():
        values = {str(value).strip().lower() for value in row.tolist() if not pd.isna(value)}
        if {"gene", "cluster"}.issubset(values):
            return int(idx)
    return None


def detect_metadata_header_row(raw: pd.DataFrame) -> int | None:
    for idx, row in raw.iterrows():
        values = {str(value).strip().lower() for value in row.tolist() if not pd.isna(value)}
        if {"cells", "major types"}.issubset(values):
            return int(idx)
    return None


def safe_read_excel(path: Path, sheet: str, header: int | None = None) -> pd.DataFrame | None:
    try:
        return pd.read_excel(path, sheet_name=sheet, header=header)
    except Exception:
        return None


def term_pattern(term: str) -> re.Pattern[str]:
    return re.compile(rf"(?<![A-Za-z0-9_]){re.escape(term)}(?![A-Za-z0-9_])", re.IGNORECASE)


def scan_workbook_text(path: Path, terms: Iterable[str]) -> dict[str, list[str]]:
    patterns = {term: term_pattern(term) for term in terms}
    hits = {term: [] for term in terms}
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    for sheet in workbook.worksheets:
        for row in sheet.iter_rows():
            for cell in row:
                value = cell.value
                if value is None:
                    continue
                text = str(value)
                for term, pattern in patterns.items():
                    if pattern.search(text):
                        locator = f"{sheet.title}!{cell.coordinate}"
                        if locator not in hits[term]:
                            hits[term].append(locator)
                        break
    workbook.close()
    return hits


def extract_pdf_text(path: Path) -> str:
    if not path.exists() or shutil.which("pdftotext") is None:
        return ""
    result = subprocess.run(
        ["pdftotext", str(path), "-"],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
    )
    return result.stdout if result.returncode == 0 else ""


def workbook_inventory(path: Path, workbook_role: str) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    workbook_rows: list[dict[str, object]] = []
    sheet_rows: list[dict[str, object]] = []
    term_rows: list[dict[str, object]] = []

    xl = pd.ExcelFile(path)
    text_hits = scan_workbook_text(path, CURATED_TERMS)
    workbook_rows.append(
        {
            "workbook_role": workbook_role,
            "path": str(path),
            "file_name": path.name,
            "file_size_bytes": path.stat().st_size,
            "n_sheets": len(xl.sheet_names),
            "sheet_names": ";".join(xl.sheet_names),
            "accessible": path.exists(),
        }
    )

    gene_values_by_term: dict[str, set[str]] = {term: set() for term in CURATED_TERMS}
    cluster_values_by_term: dict[str, set[str]] = {term: set() for term in CURATED_TERMS}
    sheets_by_term: dict[str, set[str]] = {term: set() for term in CURATED_TERMS}

    for sheet in xl.sheet_names:
        raw = safe_read_excel(path, sheet, header=None)
        if raw is None:
            continue
        n_rows, n_cols = raw.shape
        header_row = detect_header_row(raw)
        parsed_columns: list[str] = []
        n_data_rows = 0
        n_gene_rows = 0
        n_unique_genes = 0
        cluster_labels: list[str] = []
        n_cluster_labels = 0
        marker_like = False
        metadata_header_row = detect_metadata_header_row(raw)
        metadata_like = False
        metadata_columns: list[str] = []
        metadata_label_column = ""
        metadata_label_values: list[str] = []
        n_metadata_rows = 0
        if header_row is not None:
            parsed = safe_read_excel(path, sheet, header=header_row)
            if parsed is not None:
                parsed.columns = [str(col).strip() for col in parsed.columns]
                parsed_columns = list(parsed.columns)
                marker_like = "gene" in parsed_columns and "cluster" in parsed_columns
                if marker_like:
                    work = parsed.copy()
                    work["gene"] = work["gene"].map(normalize_token)
                    work["cluster"] = work["cluster"].map(normalize_token)
                    work = work.loc[work["gene"].ne("") & work["cluster"].ne("")]
                    n_data_rows = len(work)
                    genes = sorted(set(work["gene"].str.upper()))
                    clusters = sorted(set(work["cluster"]))
                    n_gene_rows = len(work["gene"])
                    n_unique_genes = len(genes)
                    cluster_labels = clusters
                    n_cluster_labels = len(clusters)
                    upper_genes = set(genes)
                    cluster_set = set(clusters)
                    for term in CURATED_TERMS:
                        if term.upper() in upper_genes:
                            gene_values_by_term[term].add(sheet)
                    if term in cluster_set:
                        cluster_values_by_term[term].add(sheet)
        if metadata_header_row is not None:
            metadata = safe_read_excel(path, sheet, header=metadata_header_row)
            if metadata is not None:
                metadata.columns = [str(col).strip() for col in metadata.columns]
                metadata_columns = list(metadata.columns)
                metadata_like = "Cells" in metadata_columns and "Major types" in metadata_columns
                if metadata_like:
                    n_metadata_rows = len(metadata.loc[metadata["Cells"].notna()])
                    metadata_label_column = "Major types"
                    metadata_label_values = sorted(
                        set(metadata["Major types"].dropna().astype(str).str.strip()) - {"Major types"}
                    )
        for term in CURATED_TERMS:
            if term_pattern(term).search(sheet):
                sheets_by_term[term].add(sheet)

        sheet_rows.append(
            {
                "workbook_role": workbook_role,
                "file_name": path.name,
                "sheet": sheet,
                "n_raw_rows": n_rows,
                "n_raw_cols": n_cols,
                "detected_header_row_1based": "" if header_row is None else header_row + 1,
                "parsed_columns": ";".join(parsed_columns),
                "marker_like_gene_cluster_table": marker_like,
                "metadata_like_cell_table": metadata_like,
                "metadata_columns": ";".join(metadata_columns),
                "metadata_label_column": metadata_label_column,
                "metadata_label_values": ";".join(metadata_label_values),
                "n_metadata_rows": n_metadata_rows,
                "n_data_rows": n_data_rows,
                "n_gene_rows": n_gene_rows,
                "n_unique_genes": n_unique_genes,
                "n_cluster_labels": n_cluster_labels,
                "cluster_labels": ";".join(cluster_labels),
            }
        )

    for term in CURATED_TERMS:
        term_rows.append(
            {
                "workbook_role": workbook_role,
                "file_name": path.name,
                "term": term,
                "appears_in_sheet_name": bool(sheets_by_term[term]),
                "sheet_name_hits": ";".join(sorted(sheets_by_term[term])),
                "appears_as_gene_symbol": bool(gene_values_by_term[term]),
                "gene_symbol_sheet_hits": ";".join(sorted(gene_values_by_term[term])),
                "appears_as_cluster_label": bool(cluster_values_by_term[term]),
                "cluster_label_sheet_hits": ";".join(sorted(cluster_values_by_term[term])),
                "appears_anywhere_in_workbook_cells": bool(text_hits[term]),
                "n_cell_hits": len(text_hits[term]),
                "example_cell_hits": ";".join(text_hits[term][:12]),
            }
        )

    return workbook_rows, sheet_rows, term_rows


def simple_markdown_table(df: pd.DataFrame, max_rows: int = 20) -> str:
    if df.empty:
        return "_No rows._"
    display = df.head(max_rows).copy()
    for col in display.columns:
        display[col] = display[col].map(lambda value: "" if pd.isna(value) else str(value).replace("|", "\\|"))
    lines = [
        "| " + " | ".join(display.columns) + " |",
        "| " + " | ".join(["---"] * len(display.columns)) + " |",
    ]
    for row in display.to_numpy().tolist():
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def write_report(
    report_path: Path,
    workbook_df: pd.DataFrame,
    sheet_df: pd.DataFrame,
    term_df: pd.DataFrame,
    pdf_term_df: pd.DataFrame,
) -> None:
    shi_sheet_df = sheet_df.loc[sheet_df["workbook_role"].eq("shi")]
    shi_term_df = term_df.loc[term_df["workbook_role"].eq("shi")]
    curated_m = ["M2", "M3", "M4", "M5", "M6", "M7"]
    curated_m_df = shi_term_df.loc[shi_term_df["term"].isin(curated_m)]
    parsed_tables = shi_sheet_df.loc[shi_sheet_df["marker_like_gene_cluster_table"].astype(bool)]
    metadata_tables = shi_sheet_df.loc[shi_sheet_df["metadata_like_cell_table"].astype(bool)]
    lines = [
        "# Jia/Shi Workbook Accessibility And Label Inventory",
        "",
        "## Bottom Line",
        "",
        "- The Jia S9 workbook and Shi S2-S9 Excel workbooks are accessible on NFS.",
        "- The previous overlap workflow was not saying the Shi files were absent.",
        "- It was saying that the curated lineage terms `M2/M3/M4/M5/M6/M7` are not encoded as literal `cluster` values in the parsed Shi marker tables.",
        "- This inventory separates physical file access, sheet access, parsed marker-table columns, gene symbols, cluster labels, and raw workbook cell text.",
        "",
        "## Accessible Workbooks",
        "",
        simple_markdown_table(workbook_df[["workbook_role", "file_name", "file_size_bytes", "n_sheets", "accessible"]], max_rows=20),
        "",
        "## Parsed Shi Marker Tables",
        "",
        simple_markdown_table(
            parsed_tables[
                [
                    "file_name",
                    "sheet",
                    "n_data_rows",
                    "n_unique_genes",
                    "n_cluster_labels",
                    "cluster_labels",
                ]
            ],
            max_rows=30,
        ),
        "",
        "## Parsed Shi Cell-Metadata Tables",
        "",
        simple_markdown_table(
            metadata_tables[
                [
                    "file_name",
                    "sheet",
                    "n_metadata_rows",
                    "metadata_columns",
                    "metadata_label_column",
                    "metadata_label_values",
                ]
            ],
            max_rows=20,
        ),
        "",
        "## Curated M-Term Accounting In Shi Excel Files",
        "",
        simple_markdown_table(
            curated_m_df[
                [
                    "file_name",
                    "term",
                    "appears_in_sheet_name",
                    "appears_as_gene_symbol",
                    "appears_as_cluster_label",
                    "appears_anywhere_in_workbook_cells",
                    "example_cell_hits",
                ]
            ],
            max_rows=80,
        ),
        "",
        "## Caption PDF Term Hits",
        "",
        simple_markdown_table(pdf_term_df, max_rows=80),
        "",
        "## Output Tables",
        "",
        "- `tables/workbook_inventory.tsv`",
        "- `tables/sheet_inventory.tsv`",
        "- `tables/term_inventory.tsv`",
        "- `tables/pdf_caption_term_inventory.tsv`",
    ]
    report_path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).expanduser().resolve()
    jia_xlsx = resolve_path(args.jia_xlsx, project_root)
    shi_dir = resolve_path(args.shi_dir, project_root)
    outdir, table_dir, report_dir = prepare_outdir(args, project_root)

    shi_workbooks = sorted(shi_dir.glob("science.abj6641_table_s*.xlsx"))
    workbooks = [(jia_xlsx, "jia")] + [(path, "shi") for path in shi_workbooks]
    missing = [path for path, _ in workbooks if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing workbook(s): " + ", ".join(str(path) for path in missing))

    workbook_rows: list[dict[str, object]] = []
    sheet_rows: list[dict[str, object]] = []
    term_rows: list[dict[str, object]] = []
    for path, role in workbooks:
        print(f"Inventorying {role}: {path}", flush=True)
        w_rows, s_rows, t_rows = workbook_inventory(path, role)
        workbook_rows.extend(w_rows)
        sheet_rows.extend(s_rows)
        term_rows.extend(t_rows)

    pdf_path = shi_dir / "science.abj6641_table_captions.pdf"
    pdf_text = extract_pdf_text(pdf_path)
    pdf_rows = []
    for term in CURATED_TERMS:
        hits = [match.start() for match in term_pattern(term).finditer(pdf_text)]
        snippets = []
        for start in hits[:5]:
            snippets.append(re.sub(r"\s+", " ", pdf_text[max(0, start - 50) : start + 80]).strip())
        pdf_rows.append(
            {
                "file_name": pdf_path.name,
                "term": term,
                "pdf_accessible": pdf_path.exists(),
                "pdftotext_available": shutil.which("pdftotext") is not None,
                "appears_in_pdf_text": bool(hits),
                "n_pdf_hits": len(hits),
                "example_pdf_snippets": " || ".join(snippets),
            }
        )

    workbook_df = pd.DataFrame(workbook_rows)
    sheet_df = pd.DataFrame(sheet_rows)
    term_df = pd.DataFrame(term_rows)
    pdf_term_df = pd.DataFrame(pdf_rows)

    workbook_df.to_csv(table_dir / "workbook_inventory.tsv", sep="\t", index=False)
    sheet_df.to_csv(table_dir / "sheet_inventory.tsv", sep="\t", index=False)
    term_df.to_csv(table_dir / "term_inventory.tsv", sep="\t", index=False)
    pdf_term_df.to_csv(table_dir / "pdf_caption_term_inventory.tsv", sep="\t", index=False)
    (table_dir / "inventory_manifest.json").write_text(
        json.dumps(
            {
                "jia_xlsx": str(jia_xlsx),
                "shi_dir": str(shi_dir),
                "outdir": str(outdir),
                "n_workbooks": len(workbook_df),
                "n_sheets": len(sheet_df),
                "n_term_rows": len(term_df),
            },
            indent=2,
        )
        + "\n"
    )
    write_report(
        report_dir / "jia_shi_workbook_inventory_report.md",
        workbook_df,
        sheet_df,
        term_df,
        pdf_term_df,
    )
    print(f"Wrote workbook inventory to {outdir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
