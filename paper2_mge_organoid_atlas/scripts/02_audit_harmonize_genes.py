#!/usr/bin/env python3
"""Audit and map canonical feature IDs without reading expression matrices."""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re

import h5py
import pandas as pd


STUDIES = [
    "varela_div30",
    "varela_div90",
    "walsh",
    "bershteyn_2025",
    "bershteyn_2023",
    "siebert_2026",
]
ENSEMBL_RE = re.compile(r"^ENSG[0-9]+(?:\.[0-9]+)?$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    return parser.parse_args()


def read_env(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line and not line.startswith("#") and "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def decode_array(values) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def read_h5ad_features(path: Path) -> pd.DataFrame:
    # Read only /var. X and layers are never accessed.
    with h5py.File(path, "r") as handle:
        var = handle["var"]
        index_key = var.attrs.get("_index", "_index")
        if isinstance(index_key, bytes):
            index_key = index_key.decode()
        feature_ids = decode_array(var[index_key][:])
        columns = {"canonical_feature_id": feature_ids}
        for name in ["source_feature_id", "feature_symbol"]:
            columns[name] = decode_array(var[name][:]) if name in var else feature_ids
    result = pd.DataFrame(columns)
    if not result["canonical_feature_id"].is_unique:
        raise ValueError(f"Canonical feature IDs are not unique in {path}")
    return result


def parse_gtf_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for match in re.finditer(r'(\S+)\s+"([^"]*)";', text):
        attributes[match.group(1)] = match.group(2)
    return attributes


def read_gencode_genes(path: Path) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = parse_gtf_attributes(fields[8])
            gene_id_versioned = attrs["gene_id"]
            rows.append(
                {
                    "common_gene_id": gene_id_versioned.split(".", 1)[0],
                    "gencode_gene_id_versioned": gene_id_versioned,
                    "common_gene_symbol": attrs.get("gene_name", ""),
                    "gencode_gene_type": attrs.get("gene_type", attrs.get("gene_biotype", "")),
                    "chromosome": fields[0],
                    "start": fields[3],
                    "end": fields[4],
                    "strand": fields[6],
                }
            )
    genes = pd.DataFrame(rows).drop_duplicates()
    duplicated = genes[genes.duplicated("common_gene_id", keep=False)]
    if not duplicated.empty:
        inconsistent = duplicated.groupby("common_gene_id").agg(
            symbols=("common_gene_symbol", "nunique"),
            types=("gencode_gene_type", "nunique"),
        )
        inconsistent = inconsistent[(inconsistent["symbols"] > 1) | (inconsistent["types"] > 1)]
        if not inconsistent.empty:
            raise ValueError("GENCODE has inconsistent duplicate versionless gene IDs")
    return genes.drop_duplicates("common_gene_id").reset_index(drop=True)


def split_multi(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []
    return [item.strip() for item in str(value).split("|") if item.strip()]


def build_hgnc_lookups(
    path: Path, gencode_ids: set[str]
) -> tuple[dict[str, set[str]], dict[str, set[str]], dict[str, set[str]], pd.DataFrame]:
    hgnc = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    approved = hgnc.loc[hgnc["status"].eq("Approved")].copy()
    approved_lookup: dict[str, set[str]] = defaultdict(set)
    previous_lookup: dict[str, set[str]] = defaultdict(set)
    alias_lookup: dict[str, set[str]] = defaultdict(set)
    for record in approved.to_dict(orient="records"):
        ids = {
            gene_id.split(".", 1)[0]
            for gene_id in split_multi(record.get("ensembl_gene_id", ""))
            if gene_id.split(".", 1)[0] in gencode_ids
        }
        if not ids:
            continue
        approved_lookup[record["symbol"]].update(ids)
        for symbol in split_multi(record.get("prev_symbol", "")):
            previous_lookup[symbol].update(ids)
        for symbol in split_multi(record.get("alias_symbol", "")):
            alias_lookup[symbol].update(ids)
    return approved_lookup, previous_lookup, alias_lookup, approved


def classify_namespace(feature_id: str) -> str:
    if ENSEMBL_RE.fullmatch(feature_id):
        return "ensembl_gene_id"
    if re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]*", feature_id):
        return "gene_symbol_like"
    return "other"


def mapping_result(
    feature_id: str,
    gencode_by_id: dict[str, dict[str, str]],
    gencode_by_symbol: dict[str, set[str]],
    hgnc_approved: dict[str, set[str]],
    hgnc_previous: dict[str, set[str]],
    hgnc_alias: dict[str, set[str]],
) -> dict[str, object]:
    namespace = classify_namespace(feature_id)
    candidates: set[str] = set()
    method = ""
    if namespace == "ensembl_gene_id":
        gene_id = feature_id.split(".", 1)[0]
        if gene_id in gencode_by_id:
            candidates = {gene_id}
            method = "ensembl_direct_gencode50"
        else:
            method = "ensembl_absent_from_gencode50"
    else:
        exact_gencode = gencode_by_symbol.get(feature_id, set())
        if exact_gencode:
            candidates = set(exact_gencode)
            method = "symbol_exact_gencode50"
        elif hgnc_approved.get(feature_id):
            candidates = set(hgnc_approved[feature_id])
            method = "symbol_exact_hgnc_approved"
        elif hgnc_previous.get(feature_id):
            candidates = set(hgnc_previous[feature_id])
            method = "symbol_previous_hgnc"
        elif hgnc_alias.get(feature_id):
            candidates = set(hgnc_alias[feature_id])
            method = "symbol_alias_hgnc"
        else:
            method = "no_reference_match"

    ordered = sorted(candidates)
    if len(ordered) == 1:
        common_id = ordered[0]
        reference = gencode_by_id[common_id]
        status = "mapped"
        common_symbol = reference["common_gene_symbol"]
        gene_type = reference["gencode_gene_type"]
    elif len(ordered) > 1:
        common_id = ""
        common_symbol = ""
        gene_type = ""
        status = "ambiguous"
    else:
        common_id = ""
        common_symbol = ""
        gene_type = ""
        status = "unresolved"
    return {
        "input_namespace": namespace,
        "mapping_status": status,
        "mapping_method": method,
        "candidate_count": len(ordered),
        "candidate_common_gene_ids": "|".join(ordered),
        "common_gene_id": common_id,
        "common_gene_symbol": common_symbol,
        "gencode_gene_type": gene_type,
    }


def overlap_row(label: str, left: set[str], right: set[str]) -> dict[str, object]:
    intersection = len(left & right)
    union = len(left | right)
    return {
        "representation": label,
        "n_left": len(left),
        "n_right": len(right),
        "n_intersection": intersection,
        "n_union": union,
        "jaccard": intersection / union if union else 1.0,
    }


def main() -> None:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    config = read_env(run_dir / "config/resolved.env")
    canonical = Path(config["CANONICAL_DIR"])
    tables = run_dir / "tables"
    provenance = run_dir / "provenance"
    gencode_path = run_dir / "reference" / f"gencode.v{config['GENCODE_RELEASE']}.annotation.gtf.gz"
    hgnc_path = run_dir / "reference/hgnc_complete_set.txt"

    print("[gene-audit] validating frozen reference files", flush=True)
    reference_manifest = pd.read_csv(provenance / "reference_manifest.tsv", sep="\t", dtype=str)
    for record in reference_manifest.to_dict(orient="records"):
        path = run_dir / "reference" / record["reference"]
        if int(record["size_bytes"]) != path.stat().st_size or record["sha256"] != sha256_file(path):
            raise ValueError(f"Reference checksum mismatch: {path}")

    canonical_manifest = pd.read_csv(
        run_dir / "config/canonical_dataset_manifest.tsv", sep="\t", dtype=str
    )
    if canonical_manifest["study_id"].tolist() != STUDIES:
        raise ValueError("Canonical dataset order/content is not the expected six-study manifest")
    checksum_manifest = pd.read_csv(
        run_dir / "config/canonical_file_checksums.tsv", sep="\t", dtype=str
    )

    print("[gene-audit] parsing GENCODE gene annotation", flush=True)
    genes = read_gencode_genes(gencode_path)
    genes.to_csv(tables / "gencode50_gene_identity_reference.tsv.gz", sep="\t", index=False)
    gencode_by_id = genes.set_index("common_gene_id").to_dict(orient="index")
    gencode_ids = set(gencode_by_id)
    gencode_by_symbol: dict[str, set[str]] = defaultdict(set)
    for record in genes.to_dict(orient="records"):
        if record["common_gene_symbol"]:
            gencode_by_symbol[record["common_gene_symbol"]].add(record["common_gene_id"])

    print("[gene-audit] parsing HGNC approved/previous/alias symbols", flush=True)
    hgnc_approved, hgnc_previous, hgnc_alias, hgnc = build_hgnc_lookups(hgnc_path, gencode_ids)
    reference_summary = pd.DataFrame(
        [
            {"metric": "gencode_release", "value": config["GENCODE_RELEASE"]},
            {"metric": "gencode_assembly", "value": config["GENCODE_ASSEMBLY"]},
            {"metric": "gencode_unique_versionless_gene_ids", "value": len(gencode_ids)},
            {"metric": "gencode_unique_gene_symbols", "value": len(gencode_by_symbol)},
            {"metric": "gencode_ambiguous_gene_symbols", "value": sum(len(v) > 1 for v in gencode_by_symbol.values())},
            {"metric": "hgnc_approved_records", "value": hgnc.shape[0]},
            {"metric": "hgnc_approved_symbols_with_gencode_id", "value": len(hgnc_approved)},
            {"metric": "hgnc_previous_symbols_with_gencode_id", "value": len(hgnc_previous)},
            {"metric": "hgnc_alias_symbols_with_gencode_id", "value": len(hgnc_alias)},
            {"metric": "common_gene_namespace", "value": config["COMMON_GENE_NAMESPACE"]},
        ]
    )
    reference_summary.to_csv(tables / "reference_summary.tsv", sep="\t", index=False)

    all_mappings: list[pd.DataFrame] = []
    namespace_rows: list[dict[str, object]] = []
    raw_sets: dict[str, set[str]] = {}
    print("[gene-audit] reading only H5AD /var feature metadata", flush=True)
    for record in canonical_manifest.to_dict(orient="records"):
        study_id = record["study_id"]
        checksum_row = checksum_manifest.loc[
            checksum_manifest["study_id"].eq(study_id) & checksum_manifest["format"].eq("h5ad")
        ]
        if checksum_row.shape[0] != 1:
            raise ValueError(f"Expected one H5AD checksum for {study_id}")
        checksum_record = checksum_row.iloc[0]
        h5ad = canonical / checksum_record["relative_path"]
        if h5ad.stat().st_size != int(checksum_record["size_bytes"]):
            raise ValueError(f"Canonical H5AD size mismatch for {study_id}")
        if sha256_file(h5ad) != checksum_record["sha256"]:
            raise ValueError(f"Canonical H5AD checksum mismatch for {study_id}")
        features = read_h5ad_features(h5ad)
        if features.shape[0] != int(record["n_features"]):
            raise ValueError(f"Feature count mismatch for {study_id}")
        if not features["canonical_feature_id"].equals(features["source_feature_id"]):
            raise ValueError(f"Canonical/source feature IDs unexpectedly differ in {study_id}")
        raw_sets[study_id] = set(features["canonical_feature_id"])
        mapped_rows = [
            mapping_result(
                feature_id,
                gencode_by_id,
                gencode_by_symbol,
                hgnc_approved,
                hgnc_previous,
                hgnc_alias,
            )
            for feature_id in features["canonical_feature_id"]
        ]
        mapped = pd.concat([features, pd.DataFrame(mapped_rows)], axis=1)
        mapped.insert(0, "study_id", study_id)
        all_mappings.append(mapped)
        namespace_counts = mapped["input_namespace"].value_counts()
        namespace_rows.append(
            {
                "study_id": study_id,
                "n_features": mapped.shape[0],
                "n_unique_feature_ids": mapped["canonical_feature_id"].nunique(),
                "n_ensembl_gene_ids": int(namespace_counts.get("ensembl_gene_id", 0)),
                "n_gene_symbol_like": int(namespace_counts.get("gene_symbol_like", 0)),
                "n_other": int(namespace_counts.get("other", 0)),
                "feature_id_namespace_conclusion": (
                    "mixed_gene_symbols_and_ensembl_gene_ids"
                    if namespace_counts.get("ensembl_gene_id", 0) > 0
                    and namespace_counts.get("gene_symbol_like", 0) > 0
                    else "gene_symbols"
                ),
            }
        )

    mapping = pd.concat(all_mappings, ignore_index=True)
    mapped_mask = mapping["mapping_status"].eq("mapped")
    duplicate_sizes = (
        mapping.loc[mapped_mask]
        .groupby(["study_id", "common_gene_id"])["canonical_feature_id"]
        .transform("size")
    )
    mapping["duplicate_common_identity_group_size"] = 0
    mapping.loc[mapped_mask, "duplicate_common_identity_group_size"] = duplicate_sizes.astype(int)
    mapping["duplicate_common_identity"] = mapping["duplicate_common_identity_group_size"].gt(1)
    mapping["strict_one_to_one_eligible"] = mapped_mask & ~mapping["duplicate_common_identity"]

    mapping.to_csv(tables / "feature_mapping_long.tsv.gz", sep="\t", index=False)
    pd.DataFrame(namespace_rows).to_csv(tables / "feature_namespace_by_dataset.tsv", sep="\t", index=False)
    (
        mapping.groupby(["study_id", "mapping_status", "mapping_method"], dropna=False)
        .size()
        .rename("n_features")
        .reset_index()
        .to_csv(tables / "mapping_summary_by_dataset.tsv", sep="\t", index=False)
    )
    mapping.loc[mapping["mapping_status"].eq("ambiguous")].to_csv(
        tables / "ambiguous_feature_mappings.tsv", sep="\t", index=False
    )
    mapping.loc[mapping["mapping_status"].eq("unresolved")].to_csv(
        tables / "unresolved_features.tsv", sep="\t", index=False
    )
    mapping.loc[mapping["duplicate_common_identity"]].sort_values(
        ["study_id", "common_gene_id", "canonical_feature_id"]
    ).to_csv(tables / "duplicate_common_gene_id_mappings.tsv", sep="\t", index=False)

    mapped_sets: dict[str, set[str]] = {}
    strict_sets: dict[str, set[str]] = {}
    per_study_by_gene: dict[str, dict[str, list[str]]] = {}
    for study_id in STUDIES:
        subset = mapping.loc[mapping["study_id"].eq(study_id)]
        mapped_subset = subset.loc[subset["mapping_status"].eq("mapped")]
        mapped_sets[study_id] = set(mapped_subset["common_gene_id"])
        strict_sets[study_id] = set(
            subset.loc[subset["strict_one_to_one_eligible"], "common_gene_id"]
        )
        per_study_by_gene[study_id] = (
            mapped_subset.groupby("common_gene_id")["canonical_feature_id"].apply(list).to_dict()
        )

    union_mapped = set().union(*mapped_sets.values())
    presence_rows = []
    for common_id in sorted(union_mapped):
        ref = gencode_by_id[common_id]
        row: dict[str, object] = {
            "common_gene_id": common_id,
            "common_gene_symbol": ref["common_gene_symbol"],
            "gencode_gene_type": ref["gencode_gene_type"],
        }
        for study_id in STUDIES:
            source_features = per_study_by_gene[study_id].get(common_id, [])
            row[f"{study_id}_present"] = bool(source_features)
            row[f"{study_id}_strict_one_to_one"] = common_id in strict_sets[study_id]
            row[f"{study_id}_source_features"] = "|".join(source_features)
        row["n_studies_present"] = sum(common_id in mapped_sets[s] for s in STUDIES)
        row["n_studies_strict_one_to_one"] = sum(common_id in strict_sets[s] for s in STUDIES)
        presence_rows.append(row)
    presence = pd.DataFrame(presence_rows)
    presence.to_csv(tables / "mapped_gene_presence_matrix.tsv.gz", sep="\t", index=False)

    mapped_intersection = set.intersection(*(mapped_sets[s] for s in STUDIES))
    strict_intersection = set.intersection(*(strict_sets[s] for s in STUDIES))
    raw_intersection = set.intersection(*(raw_sets[s] for s in STUDIES))
    presence.loc[presence["common_gene_id"].isin(mapped_intersection)].to_csv(
        tables / "six_way_gene_intersection_identity_level.tsv", sep="\t", index=False
    )
    presence.loc[presence["common_gene_id"].isin(strict_intersection)].to_csv(
        tables / "six_way_gene_intersection_strict.tsv", sep="\t", index=False
    )
    pd.DataFrame({"source_feature_id": sorted(raw_intersection)}).to_csv(
        tables / "six_way_raw_feature_id_intersection.tsv", sep="\t", index=False
    )

    overlap_rows = []
    for left_index, left in enumerate(STUDIES):
        for right in STUDIES[left_index + 1 :]:
            for representation, sets in [
                ("raw_feature_id", raw_sets),
                ("mapped_common_gene_id", mapped_sets),
                ("strict_one_to_one_common_gene_id", strict_sets),
            ]:
                row = overlap_row(representation, sets[left], sets[right])
                row.update({"left_study": left, "right_study": right})
                overlap_rows.append(row)
    pd.DataFrame(overlap_rows).to_csv(tables / "pairwise_gene_overlap.tsv", sep="\t", index=False)

    summary_rows = []
    for study_id in STUDIES:
        subset = mapping.loc[mapping["study_id"].eq(study_id)]
        summary_rows.append(
            {
                "study_id": study_id,
                "input_features": subset.shape[0],
                "mapped_features": int(subset["mapping_status"].eq("mapped").sum()),
                "mapped_unique_common_gene_ids": len(mapped_sets[study_id]),
                "ambiguous_features": int(subset["mapping_status"].eq("ambiguous").sum()),
                "unresolved_features": int(subset["mapping_status"].eq("unresolved").sum()),
                "features_in_duplicate_common_id_groups": int(subset["duplicate_common_identity"].sum()),
                "duplicate_common_gene_ids": int(
                    subset.loc[subset["duplicate_common_identity"], "common_gene_id"].nunique()
                ),
                "strict_one_to_one_common_gene_ids": len(strict_sets[study_id]),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(tables / "gene_mapping_outcome_by_dataset.tsv", sep="\t", index=False)
    intersection_summary = pd.DataFrame(
        [
            {"intersection_definition": "raw_exact_source_feature_id", "n_genes": len(raw_intersection)},
            {"intersection_definition": "mapped_common_gene_identity_in_all_six", "n_genes": len(mapped_intersection)},
            {
                "intersection_definition": "strict_one_to_one_common_gene_identity_in_all_six",
                "n_genes": len(strict_intersection),
            },
        ]
    )
    intersection_summary.to_csv(tables / "six_way_intersection_summary.tsv", sep="\t", index=False)

    validation = pd.DataFrame(
        [
            {"check": "six_expected_canonical_inputs", "status": "PASS"},
            {"check": "canonical_h5ad_checksums", "status": "PASS"},
            {"check": "reference_checksums", "status": "PASS"},
            {"check": "only_h5ad_var_read", "status": "PASS"},
            {"check": "canonical_inputs_modified", "status": "NO"},
            {"check": "expression_matrices_loaded", "status": "NO"},
            {"check": "expression_matrices_modified", "status": "NO"},
            {"check": "datasets_concatenated", "status": "NO"},
            {"check": "normalization_performed", "status": "NO"},
            {"check": "hvg_selection_performed", "status": "NO"},
            {"check": "integration_performed", "status": "NO"},
            {"check": "review_stop_enforced", "status": "YES"},
        ]
    )
    validation.to_csv(tables / "scope_and_validation_checks.tsv", sep="\t", index=False)

    completed = datetime.now(timezone.utc).astimezone().isoformat()
    readme = f"""# Step 02: feature/gene identifier audit and harmonization reports

Status: COMPLETE — STOP FOR REVIEW

This report-only step read feature metadata from the six frozen canonical H5AD
files. It did not read expression matrices and did not modify canonical inputs.
It did not concatenate datasets, normalize data, select HVGs, create combined
objects, or run integration.

## Mapping authority and common identity

The common identity is a versionless Ensembl human gene ID from GENCODE release
{config['GENCODE_RELEASE']} ({config['GENCODE_ASSEMBLY']}). Exact Ensembl IDs and
exact GENCODE gene symbols were mapped first. When no exact GENCODE symbol was
present, a frozen HGNC complete-set snapshot was used for approved, previous,
then alias symbols. A mapping was accepted only when the active tier produced
one GENCODE gene ID. Multi-candidate mappings are `ambiguous`; absent mappings
are `unresolved`. No case-insensitive or heuristic guessing was used.

## Six-way intersections

- Raw exact source feature IDs: {len(raw_intersection):,}
- Mapped common gene identities present in all six: {len(mapped_intersection):,}
- Strict one-to-one common identities in all six: {len(strict_intersection):,}

The identity-level intersection includes a common gene if it is mapped in all
six datasets, even if multiple source features in one dataset map to it. The
strict intersection excludes any common identity with duplicate source-feature
mappings in any dataset. Both are reports only; no matrix has been collapsed or
subset. Review the duplicate and ambiguity tables before choosing a downstream
matrix policy.

## Key tables

- `feature_namespace_by_dataset.tsv`
- `gene_mapping_outcome_by_dataset.tsv`
- `feature_mapping_long.tsv.gz`
- `ambiguous_feature_mappings.tsv`
- `unresolved_features.tsv`
- `duplicate_common_gene_id_mappings.tsv`
- `mapped_gene_presence_matrix.tsv.gz`
- `six_way_gene_intersection_identity_level.tsv`
- `six_way_gene_intersection_strict.tsv`
- `six_way_raw_feature_id_intersection.tsv`
- `pairwise_gene_overlap.tsv`
- `scope_and_validation_checks.tsv`

Reference files, source/reference checksums, exact code/configuration, SLURM
logs, Git state, and package checksums are included in this run directory.

Completed: {completed}

## Review stop

Do not create a common matrix or proceed to metadata harmonization until the
mapping, ambiguous aliases, unresolved features, duplicate common identities,
and strict-vs-identity-level intersection policy have been reviewed.
"""
    (run_dir / "README.md").write_text(readme, encoding="utf-8")

    package_rows = []
    for path in sorted(run_dir.rglob("*")):
        relative = path.relative_to(run_dir)
        if (
            path.is_file()
            and relative.parts[0] != "logs"
            and relative != Path("provenance/package_checksums.tsv")
            and relative != Path("SUCCESS.txt")
        ):
            package_rows.append(
                {
                    "relative_path": str(relative),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
    pd.DataFrame(package_rows).to_csv(
        provenance / "package_checksums.tsv", sep="\t", index=False
    )
    (run_dir / "SUCCESS.txt").write_text(
        "PASS\n"
        f"completed={completed}\n"
        "scope=feature_identifier_audit_and_mapping_reports_only\n"
        f"six_way_mapped_identity_intersection={len(mapped_intersection)}\n"
        f"six_way_strict_one_to_one_intersection={len(strict_intersection)}\n"
        "review_stop=YES\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": "PASS",
        "mapped_intersection": len(mapped_intersection),
        "strict_intersection": len(strict_intersection),
        "raw_intersection": len(raw_intersection),
    }), flush=True)


if __name__ == "__main__":
    main()
