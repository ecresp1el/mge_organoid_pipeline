#!/usr/bin/env python3
"""
Convert Shi et al Table S3 (.xlsx) to a TSV without external Python dependencies.

This parser reads OOXML internals directly (zip + XML) and supports:
- shared strings
- inline strings
- numeric cells

Example:
  python scripts/05h_shi_table_s3_xlsx_to_tsv.py \
    --xlsx /path/to/science.abj6641_table_s3.xlsx \
    --sheet 1 \
    --out /path/to/science.abj6641_table_s3.tsv
"""

import argparse
import csv
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple


MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
DOC_REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
NS = {"m": MAIN_NS, "pr": PKG_REL_NS}


def col_ref_to_zero_based(cell_ref: str) -> int:
    m = re.match(r"^([A-Za-z]+)", cell_ref)
    if not m:
        raise ValueError(f"Invalid cell reference: {cell_ref}")
    letters = m.group(1).upper()
    idx = 0
    for ch in letters:
        idx = idx * 26 + (ord(ch) - ord("A") + 1)
    return idx - 1


def parse_shared_strings(zf: zipfile.ZipFile) -> List[str]:
    shared_path = "xl/sharedStrings.xml"
    if shared_path not in zf.namelist():
        return []
    root = ET.fromstring(zf.read(shared_path))
    out = []
    for si in root.findall(".//m:si", NS):
        parts = []
        for t in si.findall(".//m:t", NS):
            parts.append(t.text or "")
        out.append("".join(parts))
    return out


def discover_sheets(zf: zipfile.ZipFile) -> List[Tuple[str, str]]:
    workbook_root = ET.fromstring(zf.read("xl/workbook.xml"))
    rels_root = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))

    rid_to_target = {}
    for rel in rels_root.findall(".//pr:Relationship", NS):
        rid = rel.get("Id")
        target = rel.get("Target")
        if rid and target:
            if target.startswith("/"):
                norm = target.lstrip("/")
            else:
                norm = f"xl/{target}" if not target.startswith("xl/") else target
            rid_to_target[rid] = norm

    sheets = []
    for sheet in workbook_root.findall(".//m:sheets/m:sheet", NS):
        name = sheet.get("name", "")
        rid = sheet.get(f"{{{DOC_REL_NS}}}id")
        if not rid:
            continue
        target = rid_to_target.get(rid)
        if target:
            sheets.append((name, target))
    return sheets


def cell_value(cell: ET.Element, shared: List[str]) -> str:
    cell_type = cell.get("t", "")
    if cell_type == "inlineStr":
        # inline string may be in is/t (possibly multiple rich text runs)
        texts = [t.text or "" for t in cell.findall(".//m:is//m:t", NS)]
        return "".join(texts)

    v = cell.find("m:v", NS)
    if v is None or v.text is None:
        return ""
    raw = v.text

    if cell_type == "s":
        try:
            idx = int(raw)
            if 0 <= idx < len(shared):
                return shared[idx]
            return raw
        except ValueError:
            return raw
    if cell_type == "b":
        return "TRUE" if raw == "1" else "FALSE"
    return raw


def parse_sheet_rows(zf: zipfile.ZipFile, sheet_xml_path: str, shared: List[str]) -> List[List[str]]:
    root = ET.fromstring(zf.read(sheet_xml_path))
    rows_out = []

    for row in root.findall(".//m:sheetData/m:row", NS):
        row_map = {}
        max_col = -1

        fallback_col = 0
        for cell in row.findall("m:c", NS):
            ref = cell.get("r")
            if ref:
                col_idx = col_ref_to_zero_based(ref)
            else:
                col_idx = fallback_col
            fallback_col = col_idx + 1

            val = cell_value(cell, shared)
            row_map[col_idx] = val
            if col_idx > max_col:
                max_col = col_idx

        if max_col < 0:
            rows_out.append([])
            continue

        dense = [""] * (max_col + 1)
        for idx, val in row_map.items():
            dense[idx] = val
        rows_out.append(dense)

    return rows_out


def main() -> int:
    p = argparse.ArgumentParser(description="Convert an XLSX worksheet to TSV (no external deps).")
    p.add_argument("--xlsx", required=True, help="Path to .xlsx file")
    p.add_argument("--sheet", default="1", help="Sheet index (1-based) or sheet name (default: 1)")
    p.add_argument("--out", required=True, help="Output TSV path")
    p.add_argument("--list-sheets", action="store_true", help="Print sheet list and exit")
    args = p.parse_args()

    xlsx_path = Path(args.xlsx)
    out_path = Path(args.out)

    if not xlsx_path.exists():
        print(f"ERROR: XLSX not found: {xlsx_path}", file=sys.stderr)
        return 2

    with zipfile.ZipFile(xlsx_path, "r") as zf:
        sheets = discover_sheets(zf)
        if not sheets:
            print("ERROR: No worksheets found in workbook", file=sys.stderr)
            return 3

        if args.list_sheets:
            for i, (name, target) in enumerate(sheets, start=1):
                print(f"{i}\t{name}\t{target}")
            return 0

        chosen_name = None
        chosen_target = None

        if re.fullmatch(r"[0-9]+", args.sheet):
            idx = int(args.sheet)
            if idx < 1 or idx > len(sheets):
                print(f"ERROR: Sheet index out of range: {idx} (1..{len(sheets)})", file=sys.stderr)
                return 4
            chosen_name, chosen_target = sheets[idx - 1]
        else:
            for name, target in sheets:
                if name == args.sheet:
                    chosen_name, chosen_target = name, target
                    break
            if chosen_target is None:
                print(f"ERROR: Sheet name not found: {args.sheet}", file=sys.stderr)
                print("Available sheets:", file=sys.stderr)
                for i, (name, target) in enumerate(sheets, start=1):
                    print(f"  {i}: {name} ({target})", file=sys.stderr)
                return 5

        shared = parse_shared_strings(zf)
        rows = parse_sheet_rows(zf, chosen_target, shared)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for r in rows:
            writer.writerow(r)

    nrows = len(rows)
    ncols = max((len(r) for r in rows), default=0)
    print(f"WROTE_TSV={out_path}")
    print(f"SHEET={chosen_name}")
    print(f"NROWS={nrows}")
    print(f"NCOLS={ncols}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
