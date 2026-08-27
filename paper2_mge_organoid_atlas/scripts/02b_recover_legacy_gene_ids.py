#!/usr/bin/env python3
"""Recover candidate legacy gene IDs without changing inputs or Step 02 outputs."""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
from pathlib import Path
import re

import pandas as pd


STUDIES = [
    "varela_div30",
    "varela_div90",
    "walsh",
    "bershteyn_2025",
    "bershteyn_2023",
    "siebert_2026",
]
RELEASES = ["27", "32", "35", "44", "50"]
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


def parse_gtf_attributes(text: str) -> dict[str, str]:
    return {
        match.group(1): match.group(2)
        for match in re.finditer(r'(\S+)\s+"([^"]*)";', text)
    }


def read_gencode_gene_rows(path: Path, release: str) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id", "").split(".", 1)[0]
            symbol = attrs.get("gene_name", "")
            if gene_id:
                rows.append(
                    {
                        "reference_release": release,
                        "common_gene_id": gene_id,
                        "gene_symbol": symbol,
                        "gene_type": attrs.get(
                            "gene_type", attrs.get("gene_biotype", "")
                        ),
                    }
                )
    return pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)


def read_cellranger_features(path: Path) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Malformed Cell Ranger feature row {line_number}: {path}")
            feature_type = fields[2] if len(fields) >= 3 else "Gene Expression"
            if feature_type != "Gene Expression":
                continue
            rows.append(
                {
                    "source_gene_id": fields[0].split(".", 1)[0],
                    "source_gene_symbol": fields[1],
                    "feature_type": feature_type,
                }
            )
    return pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)


def join_ids(values: set[str] | list[str]) -> str:
    return "|".join(sorted(set(values)))


def build_reference_stability(
    release_tables: dict[str, pd.DataFrame],
) -> tuple[pd.DataFrame, dict[str, dict[str, set[str]]], dict[str, set[str]], dict[str, str]]:
    by_release_symbol: dict[str, dict[str, set[str]]] = {}
    ids_by_release: dict[str, set[str]] = {}
    newest_symbol_by_id: dict[str, str] = {}
    all_symbols: set[str] = set()
    for release in RELEASES:
        table = release_tables[release]
        symbol_map: dict[str, set[str]] = defaultdict(set)
        for row in table.itertuples(index=False):
            if row.gene_symbol:
                symbol_map[row.gene_symbol].add(row.common_gene_id)
                all_symbols.add(row.gene_symbol)
            newest_symbol_by_id[row.common_gene_id] = row.gene_symbol
        by_release_symbol[release] = symbol_map
        ids_by_release[release] = set(table["common_gene_id"])

    rows: list[dict[str, object]] = []
    for symbol in sorted(all_symbols):
        observed_releases = [r for r in RELEASES if symbol in by_release_symbol[r]]
        union_ids = set().union(
            *(by_release_symbol[r][symbol] for r in observed_releases)
        )
        if len(union_ids) == 1 and len(observed_releases) >= 2:
            status = "historical_consensus_unique"
        elif len(union_ids) == 1:
            status = "historical_single_release_unique"
        else:
            status = "historical_ambiguous"
        row: dict[str, object] = {
            "gene_symbol": symbol,
            "n_releases_observed": len(observed_releases),
            "releases_observed": "|".join(observed_releases),
            "union_candidate_count": len(union_ids),
            "union_candidate_common_gene_ids": join_ids(union_ids),
            "historical_symbol_status": status,
        }
        for release in RELEASES:
            row[f"gencode_{release}_candidate_ids"] = join_ids(
                by_release_symbol[release].get(symbol, set())
            )
        rows.append(row)
    return pd.DataFrame(rows), by_release_symbol, ids_by_release, newest_symbol_by_id


def historical_candidate(
    feature_id: str,
    stability_by_symbol: dict[str, dict[str, object]],
    ids_by_release: dict[str, set[str]],
) -> dict[str, object]:
    if ENSEMBL_RE.fullmatch(feature_id):
        gene_id = feature_id.split(".", 1)[0]
        releases = [release for release in RELEASES if gene_id in ids_by_release[release]]
        if releases:
            return {
                "status": "mapped",
                "ids": {gene_id},
                "evidence": "historical_ensembl_id_present",
                "releases": releases,
            }
        return {"status": "unresolved", "ids": set(), "evidence": "absent_all_references", "releases": []}

    record = stability_by_symbol.get(feature_id)
    if record is None:
        return {"status": "unresolved", "ids": set(), "evidence": "absent_all_references", "releases": []}
    ids = set(str(record["union_candidate_common_gene_ids"]).split("|"))
    releases = str(record["releases_observed"]).split("|")
    status = str(record["historical_symbol_status"])
    if status in {"historical_consensus_unique", "historical_single_release_unique"}:
        return {"status": "mapped", "ids": ids, "evidence": status, "releases": releases}
    return {"status": "ambiguous", "ids": ids, "evidence": status, "releases": releases}


def package_checksums(run_dir: Path) -> pd.DataFrame:
    rows = []
    for path in sorted(run_dir.rglob("*")):
        relative = path.relative_to(run_dir)
        if (
            path.is_file()
            and relative.parts[0] != "logs"
            and relative != Path("provenance/package_checksums.tsv")
            and relative != Path("SUCCESS.txt")
        ):
            rows.append(
                {
                    "relative_path": str(relative),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    tables_dir = run_dir / "tables"
    reference_dir = run_dir / "reference"
    provenance_dir = run_dir / "provenance"
    config = read_env(run_dir / "config/resolved.env")
    evidence = pd.read_csv(
        provenance_dir / "source_evidence_registry.tsv", sep="\t", dtype=str
    ).fillna("")

    print("[02b] validating every registered evidence file", flush=True)
    evidence_paths: dict[str, Path] = {}
    for row in evidence.itertuples(index=False):
        path = Path(row.path)
        if not path.is_file():
            raise ValueError(f"Registered evidence file is missing: {path}")
        if path.stat().st_size != int(row.size_bytes) or sha256_file(path) != row.sha256:
            raise ValueError(f"Registered evidence checksum mismatch: {path}")
        evidence_paths[row.evidence_id] = path

    print("[02b] parsing five GENCODE releases", flush=True)
    release_tables: dict[str, pd.DataFrame] = {}
    for release in RELEASES:
        gene_table = read_gencode_gene_rows(evidence_paths[f"gencode{release}"], release)
        if gene_table.empty:
            raise ValueError(f"No genes parsed from GENCODE {release}")
        gene_table.to_csv(
            reference_dir / f"gencode_{release}_gene_identity.tsv.gz",
            sep="\t",
            index=False,
        )
        release_tables[release] = gene_table

    stability, _, ids_by_release, newest_symbol_by_id = build_reference_stability(
        release_tables
    )
    stability.to_csv(
        tables_dir / "historical_symbol_reference_stability.tsv.gz",
        sep="\t",
        index=False,
    )
    stability_by_symbol = stability.set_index("gene_symbol").to_dict(orient="index")

    reference_summary = []
    for release in RELEASES:
        table = release_tables[release]
        symbol_sizes = table.loc[table["gene_symbol"].ne("")].groupby("gene_symbol")[
            "common_gene_id"
        ].nunique()
        reference_summary.append(
            {
                "reference_release": release,
                "n_gene_rows": table.shape[0],
                "n_unique_versionless_ensembl_ids": table["common_gene_id"].nunique(),
                "n_unique_symbols": symbol_sizes.shape[0],
                "n_symbols_mapping_to_multiple_ids": int(symbol_sizes.gt(1).sum()),
            }
        )
    pd.DataFrame(reference_summary).to_csv(
        tables_dir / "historical_reference_summary.tsv", sep="\t", index=False
    )

    print("[02b] validating and comparing the Varela source feature tables", flush=True)
    varela_features: dict[str, pd.DataFrame] = {}
    varela_lookups: dict[str, dict[str, set[str]]] = {}
    for study in ["varela_div30", "varela_div90"]:
        frame = read_cellranger_features(evidence_paths[f"{study}_features"])
        lookup: dict[str, set[str]] = defaultdict(set)
        for row in frame.itertuples(index=False):
            lookup[row.source_gene_symbol].add(row.source_gene_id)
        varela_features[study] = frame
        varela_lookups[study] = lookup

    same_varela_table = varela_features["varela_div30"].equals(
        varela_features["varela_div90"]
    )
    source_table_summary = []
    for study, frame in varela_features.items():
        symbol_sizes = frame.groupby("source_gene_symbol")["source_gene_id"].nunique()
        source_table_summary.append(
            {
                "study_id": study,
                "n_gene_expression_rows": frame.shape[0],
                "n_unique_gene_ids": frame["source_gene_id"].nunique(),
                "n_unique_gene_symbols": frame["source_gene_symbol"].nunique(),
                "n_symbols_with_multiple_gene_ids": int(symbol_sizes.gt(1).sum()),
                "div30_div90_tables_exactly_equal": same_varela_table,
            }
        )
    pd.DataFrame(source_table_summary).to_csv(
        tables_dir / "varela_source_feature_table_summary.tsv", sep="\t", index=False
    )

    study_evidence = pd.DataFrame(
        [
            {
                "study_id": "varela_div30",
                "original_feature_table_found": "YES",
                "exact_reference_bundle_identified": "source feature table found; bundle name not required for row identity",
                "recovery_policy": "confirmed source feature-table symbol-to-Ensembl mapping first",
            },
            {
                "study_id": "varela_div90",
                "original_feature_table_found": "YES",
                "exact_reference_bundle_identified": "source feature table found; bundle name not required for row identity",
                "recovery_policy": "confirmed source feature-table symbol-to-Ensembl mapping first",
            },
            {
                "study_id": "walsh",
                "original_feature_table_found": "NO",
                "exact_reference_bundle_identified": "NO; GEO records Cell Ranger 6.1.2 and GRCh38 but not the exact reference bundle",
                "recovery_policy": "historical-reference consensus candidate only",
            },
            {
                "study_id": "bershteyn_2025",
                "original_feature_table_found": "NO",
                "exact_reference_bundle_identified": "NO; GEO records Cell Ranger 7/GRCh38 but supplies the Seurat object only",
                "recovery_policy": "historical-reference consensus candidate only",
            },
            {
                "study_id": "bershteyn_2023",
                "original_feature_table_found": "NO",
                "exact_reference_bundle_identified": "NO; GEO records Cell Ranger 3.1.0/GRCh38 but not the exact reference bundle",
                "recovery_policy": "historical-reference consensus candidate only",
            },
            {
                "study_id": "siebert_2026",
                "original_feature_table_found": "NO",
                "exact_reference_bundle_identified": "NO; local NeMO metadata has no original feature-table URL",
                "recovery_policy": "historical-reference consensus candidate only",
            },
        ]
    )
    study_evidence.to_csv(
        tables_dir / "study_source_evidence_and_recovery_policy.tsv", sep="\t", index=False
    )

    print("[02b] producing report-only candidate mappings", flush=True)
    mapping = pd.read_csv(
        evidence_paths["step02_mapping"], sep="\t", compression="gzip", dtype=str
    ).fillna("")
    if mapping["study_id"].drop_duplicates().tolist() != STUDIES:
        raise ValueError("Parent Step 02 mapping is not the expected six-study mapping")
    if mapping.duplicated(["study_id", "canonical_feature_id"]).any():
        raise ValueError("Parent Step 02 mapping contains duplicate study/feature rows")

    proposed_rows: list[dict[str, object]] = []
    for row in mapping.to_dict(orient="records"):
        study = row["study_id"]
        feature_id = row["canonical_feature_id"]
        step02_status = row["mapping_status"]
        step02_id = row["common_gene_id"]
        exact_source_ids: set[str] = set()
        exact_source_state = "not_available"
        if study in varela_lookups:
            exact_source_ids = varela_lookups[study].get(feature_id, set())
            if len(exact_source_ids) == 1:
                exact_source_state = "unique"
            elif len(exact_source_ids) > 1:
                exact_source_state = "ambiguous"
            else:
                exact_source_state = "absent"

        historical = historical_candidate(feature_id, stability_by_symbol, ids_by_release)
        proposed_status = "unresolved"
        proposed_ids: set[str] = set()
        evidence_label = "unresolved_after_all_references"
        evidence_tier = "5_unresolved"

        if len(exact_source_ids) == 1:
            proposed_status = "mapped"
            proposed_ids = exact_source_ids
            evidence_label = "confirmed_source_feature_table"
            evidence_tier = "1_confirmed_source_feature_table"
        elif step02_status == "mapped" and step02_id:
            proposed_status = "mapped"
            proposed_ids = {step02_id}
            evidence_label = "step02_current_reference"
            if exact_source_state == "ambiguous":
                evidence_label += "_fallback_source_symbol_ambiguous"
            elif exact_source_state == "absent":
                evidence_label += "_fallback_source_symbol_absent"
            evidence_tier = "2_step02_current_reference"
        elif historical["status"] == "mapped":
            proposed_status = "mapped"
            proposed_ids = set(historical["ids"])
            evidence_label = str(historical["evidence"])
            if evidence_label in {
                "historical_consensus_unique",
                "historical_ensembl_id_present",
            }:
                evidence_tier = "3_historical_reference_consensus"
            else:
                evidence_tier = "4_historical_single_release"
        elif historical["status"] == "ambiguous" or exact_source_state == "ambiguous":
            proposed_status = "ambiguous"
            proposed_ids = set(historical["ids"]) | exact_source_ids
            evidence_label = (
                "confirmed_source_symbol_ambiguous"
                if exact_source_state == "ambiguous"
                else "historical_reference_ambiguous"
            )
            evidence_tier = "5_ambiguous"

        proposed_id = next(iter(proposed_ids)) if proposed_status == "mapped" else ""
        if step02_status != "mapped" and proposed_status == "mapped":
            mapping_action = "recovered_candidate"
        elif (
            step02_status == "mapped"
            and proposed_status == "mapped"
            and step02_id != proposed_id
        ):
            mapping_action = "reassigned_by_confirmed_source_table"
        elif step02_status == proposed_status and step02_id == proposed_id:
            mapping_action = "unchanged"
        else:
            mapping_action = "remains_without_unique_mapping"

        proposed_rows.append(
            {
                **row,
                "source_feature_table_match_state": exact_source_state,
                "source_feature_table_candidate_ids": join_ids(exact_source_ids),
                "historical_reference_status": historical["status"],
                "historical_reference_evidence": historical["evidence"],
                "historical_reference_releases": "|".join(historical["releases"]),
                "historical_reference_candidate_ids": join_ids(historical["ids"]),
                "proposed_mapping_status": proposed_status,
                "proposed_common_gene_id": proposed_id,
                "proposed_common_gene_symbol": newest_symbol_by_id.get(proposed_id, ""),
                "proposed_candidate_common_gene_ids": join_ids(proposed_ids),
                "proposed_mapping_evidence": evidence_label,
                "proposed_evidence_tier": evidence_tier,
                "mapping_action_vs_step02": mapping_action,
                "application_status": "REPORT_ONLY_NOT_APPLIED",
            }
        )

    proposed = pd.DataFrame(proposed_rows)
    mapped_mask = proposed["proposed_mapping_status"].eq("mapped")
    group_sizes = (
        proposed.loc[mapped_mask]
        .groupby(["study_id", "proposed_common_gene_id"])["canonical_feature_id"]
        .transform("size")
    )
    proposed["proposed_duplicate_group_size"] = 0
    proposed.loc[mapped_mask, "proposed_duplicate_group_size"] = group_sizes.astype(int)
    proposed["proposed_duplicate_common_identity"] = proposed[
        "proposed_duplicate_group_size"
    ].gt(1)
    proposed["proposed_strict_one_to_one_eligible"] = mapped_mask & ~proposed[
        "proposed_duplicate_common_identity"
    ]
    proposed.to_csv(
        tables_dir / "recovery_candidate_mapping_long.tsv.gz", sep="\t", index=False
    )
    proposed.loc[proposed["mapping_action_vs_step02"].ne("unchanged")].to_csv(
        tables_dir / "recovery_changes_and_unresolved.tsv.gz", sep="\t", index=False
    )
    proposed.loc[proposed["proposed_mapping_status"].eq("unresolved")].to_csv(
        tables_dir / "unresolved_after_recovery.tsv", sep="\t", index=False
    )
    proposed.loc[proposed["proposed_mapping_status"].eq("ambiguous")].to_csv(
        tables_dir / "ambiguous_after_recovery.tsv", sep="\t", index=False
    )
    proposed.loc[proposed["proposed_duplicate_common_identity"]].sort_values(
        ["study_id", "proposed_common_gene_id", "canonical_feature_id"]
    ).to_csv(
        tables_dir / "proposed_duplicate_common_gene_id_mappings.tsv",
        sep="\t",
        index=False,
    )

    outcome_rows = []
    proposed_sets: dict[str, set[str]] = {}
    strict_sets: dict[str, set[str]] = {}
    for study in STUDIES:
        subset = proposed.loc[proposed["study_id"].eq(study)]
        mapped = subset.loc[subset["proposed_mapping_status"].eq("mapped")]
        proposed_sets[study] = set(mapped["proposed_common_gene_id"])
        strict_sets[study] = set(
            subset.loc[
                subset["proposed_strict_one_to_one_eligible"],
                "proposed_common_gene_id",
            ]
        )
        outcome_rows.append(
            {
                "study_id": study,
                "input_features": subset.shape[0],
                "step02_mapped_features": int(subset["mapping_status"].eq("mapped").sum()),
                "newly_recovered_candidate_features": int(
                    subset["mapping_action_vs_step02"].eq("recovered_candidate").sum()
                ),
                "confirmed_source_table_mapped_features": int(
                    subset["proposed_mapping_evidence"].eq(
                        "confirmed_source_feature_table"
                    ).sum()
                ),
                "historical_consensus_mapped_features": int(
                    subset["proposed_evidence_tier"].eq(
                        "3_historical_reference_consensus"
                    ).sum()
                ),
                "historical_single_release_mapped_features": int(
                    subset["proposed_evidence_tier"].eq(
                        "4_historical_single_release"
                    ).sum()
                ),
                "reassigned_by_confirmed_source_table": int(
                    subset["mapping_action_vs_step02"].eq(
                        "reassigned_by_confirmed_source_table"
                    ).sum()
                ),
                "proposed_mapped_features": mapped.shape[0],
                "proposed_unique_common_gene_ids": len(proposed_sets[study]),
                "remaining_ambiguous_features": int(
                    subset["proposed_mapping_status"].eq("ambiguous").sum()
                ),
                "remaining_unresolved_features": int(
                    subset["proposed_mapping_status"].eq("unresolved").sum()
                ),
                "proposed_duplicate_common_gene_ids": int(
                    subset.loc[
                        subset["proposed_duplicate_common_identity"],
                        "proposed_common_gene_id",
                    ].nunique()
                ),
                "proposed_strict_one_to_one_ids": len(strict_sets[study]),
                "application_status": "REPORT_ONLY_NOT_APPLIED",
            }
        )
    outcome = pd.DataFrame(outcome_rows)
    outcome.to_csv(
        tables_dir / "recovery_outcome_by_dataset.tsv", sep="\t", index=False
    )

    identity_intersection = set.intersection(*(proposed_sets[s] for s in STUDIES))
    strict_intersection = set.intersection(*(strict_sets[s] for s in STUDIES))
    original_identity = set(
        mapping.loc[
            mapping.groupby("common_gene_id")["study_id"].transform("nunique").eq(6)
            & mapping["common_gene_id"].ne(""),
            "common_gene_id",
        ]
    )
    original_strict = set.intersection(
        *(
            set(
                mapping.loc[
                    mapping["study_id"].eq(study)
                    & mapping["strict_one_to_one_eligible"].eq("True"),
                    "common_gene_id",
                ]
            )
            for study in STUDIES
        )
    )
    intersection_summary = pd.DataFrame(
        [
            {
                "intersection_definition": "step02_mapped_identity_level",
                "n_genes": len(original_identity),
                "application_status": "PARENT_REPORT_ONLY",
            },
            {
                "intersection_definition": "step02_strict_one_to_one",
                "n_genes": len(original_strict),
                "application_status": "PARENT_REPORT_ONLY",
            },
            {
                "intersection_definition": "02b_proposed_identity_level",
                "n_genes": len(identity_intersection),
                "application_status": "REPORT_ONLY_NOT_APPLIED",
            },
            {
                "intersection_definition": "02b_proposed_strict_one_to_one",
                "n_genes": len(strict_intersection),
                "application_status": "REPORT_ONLY_NOT_APPLIED",
            },
        ]
    )
    intersection_summary.to_csv(
        tables_dir / "proposed_six_way_intersection_summary.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        {"proposed_common_gene_id": sorted(identity_intersection)}
    ).to_csv(
        tables_dir / "proposed_six_way_intersection_identity_level.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(
        {"proposed_common_gene_id": sorted(strict_intersection)}
    ).to_csv(
        tables_dir / "proposed_six_way_intersection_strict.tsv", sep="\t", index=False
    )

    validation = pd.DataFrame(
        [
            {"check": "registered_evidence_checksums", "status": "PASS"},
            {"check": "five_historical_gencode_releases_parsed", "status": "PASS"},
            {"check": "varela_source_feature_tables_exactly_equal", "status": "PASS" if same_varela_table else "FAIL"},
            {"check": "canonical_inputs_opened", "status": "NO"},
            {"check": "canonical_inputs_modified", "status": "NO"},
            {"check": "parent_step02_outputs_modified", "status": "NO"},
            {"check": "expression_matrices_loaded", "status": "NO"},
            {"check": "datasets_concatenated", "status": "NO"},
            {"check": "normalization_or_hvg_selection", "status": "NO"},
            {"check": "integration_performed", "status": "NO"},
            {"check": "candidate_mappings_applied", "status": "NO"},
            {"check": "review_stop_enforced", "status": "YES"},
        ]
    )
    if validation["status"].eq("FAIL").any():
        raise ValueError("One or more validation checks failed")
    validation.to_csv(
        tables_dir / "scope_and_validation_checks.tsv", sep="\t", index=False
    )

    completed = datetime.now(timezone.utc).astimezone().isoformat()
    recovered_total = int(
        proposed["mapping_action_vs_step02"].eq("recovered_candidate").sum()
    )
    readme = f"""# Step 02b: legacy gene-ID recovery audit

Status: COMPLETE — REPORT ONLY — STOP FOR REVIEW

This extension asks whether unresolved or ambiguous source feature names can be
connected to stable versionless Ensembl gene IDs using stronger historical
evidence. It does not alter the six canonical inputs or the completed Step 02
mapping package. No expression matrix was opened, concatenated, normalized,
filtered, or integrated.

## Evidence policy

Varela DIV30 and DIV90 have original project Cell Ranger `features.tsv.gz`
files. A unique exact symbol match in those tables is labeled
`confirmed_source_feature_table`. The two tables are checksum-identical and
their parsed Gene Expression rows are exactly equal.

The exact original feature tables/reference bundles were not found for Walsh,
Bershteyn 2025, Bershteyn 2023, or Siebert 2026. For those datasets, GENCODE
27, 32, 35, 44, and 50 are used only to generate candidate identities:

- `historical_consensus_unique`: the exact symbol maps to one identical
  versionless Ensembl ID across at least two releases where it occurs;
- `historical_single_release_unique`: found uniquely in only one tested
  release, and therefore lower confidence;
- `historical_reference_ambiguous`: multiple Ensembl IDs remain possible;
- `absent_all_references`: no exact recovery was found.

Historical consensus is not represented as proof that a particular study used
that reference bundle. Every candidate row is marked `REPORT_ONLY_NOT_APPLIED`.

## Reported effect

- Candidate features newly assigned beyond Step 02: {recovered_total:,}
- Parent Step 02 identity-level six-way intersection: {len(original_identity):,}
- 02b proposed identity-level six-way intersection: {len(identity_intersection):,}
- Parent Step 02 strict six-way intersection: {len(original_strict):,}
- 02b proposed strict six-way intersection: {len(strict_intersection):,}

These proposed intersections are audit outputs, not approved feature sets. The
per-study evidence tiers, changed rows, remaining ambiguous/unresolved names,
and duplicate mappings must be reviewed before any mapping is adopted.

## Key tables

- `study_source_evidence_and_recovery_policy.tsv`
- `historical_reference_summary.tsv`
- `historical_symbol_reference_stability.tsv.gz`
- `varela_source_feature_table_summary.tsv`
- `recovery_candidate_mapping_long.tsv.gz`
- `recovery_changes_and_unresolved.tsv.gz`
- `recovery_outcome_by_dataset.tsv`
- `ambiguous_after_recovery.tsv`
- `unresolved_after_recovery.tsv`
- `proposed_duplicate_common_gene_id_mappings.tsv`
- `proposed_six_way_intersection_summary.tsv`
- `scope_and_validation_checks.tsv`

The source-evidence registry records exact paths, sizes, SHA-256 checksums, and
the claim each source is allowed to support. Parsed identity tables for every
GENCODE release are frozen under `reference/` to keep the report reproducible
without copying five large GTF files into the package.

Completed: {completed}

## Review stop

Do not rewrite canonical objects, create a common matrix, or proceed to Step 03
until the evidence tiers and proposed duplicate/intersection policy are
reviewed explicitly.
"""
    (run_dir / "README.md").write_text(readme, encoding="utf-8")
    package_checksums(run_dir).to_csv(
        provenance_dir / "package_checksums.tsv", sep="\t", index=False
    )
    (run_dir / "SUCCESS.txt").write_text(
        "PASS\n"
        f"completed={completed}\n"
        "scope=legacy_gene_id_recovery_report_only\n"
        f"newly_recovered_candidate_features={recovered_total}\n"
        f"proposed_six_way_identity_intersection={len(identity_intersection)}\n"
        f"proposed_six_way_strict_intersection={len(strict_intersection)}\n"
        "canonical_inputs_modified=NO\n"
        "candidate_mappings_applied=NO\n"
        "review_stop=YES\n",
        encoding="utf-8",
    )
    print(
        {
            "status": "PASS",
            "newly_recovered_candidate_features": recovered_total,
            "proposed_identity_intersection": len(identity_intersection),
            "proposed_strict_intersection": len(strict_intersection),
        },
        flush=True,
    )


if __name__ == "__main__":
    main()
