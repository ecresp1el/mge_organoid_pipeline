#!/usr/bin/env python3
"""Modular GSE94641 reference-PCA and kNN label-transfer implementation.

The full four-age published reference defines the feature set and PCA. E15.5
cells are the primary neighbor pool; all ages provide a separately labeled
context result. Every published label field is transferred independently and
as an exact tuple. PCDH19 classifications are never read or used as features.
"""

from __future__ import print_function

import csv
import gzip
import hashlib
import json
import os
import platform
import re
import shutil
import sys
import tempfile
import urllib.request
from collections import Counter, OrderedDict

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.sparse as sparse


class MappingError(RuntimeError):
    pass


def sha256_file(path, chunk_size=8 * 1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=header,
            lineterminator="\n", extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in header})


def read_tsv(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def decode(values):
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


class LockedGeneAnnotation(object):
    """Materialize and parse the exact NCBI mouse Entrez/Ensembl crosswalk."""

    ENSEMBL_PATTERN = re.compile(r"(?:^|\|)Ensembl:(ENSMUSG\d+)(?:\||$)")

    def __init__(self, reference_root, definition):
        self.path = os.path.join(reference_root, "source_files", definition["name"])
        self.definition = definition

    def _verify(self, path):
        return (
            os.path.isfile(path)
            and os.path.getsize(path) == int(self.definition["bytes"])
            and sha256_file(path) == self.definition["sha256"]
        )

    def materialize(self):
        if self._verify(self.path):
            return self.path
        if os.path.exists(self.path):
            raise MappingError("Existing NCBI gene annotation fails its checksum lock")
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        descriptor, temporary = tempfile.mkstemp(prefix=".gene-info.", dir=os.path.dirname(self.path))
        os.close(descriptor)
        try:
            request = urllib.request.Request(
                self.definition["url"],
                headers={"User-Agent": "paper3-pcdh19-gse94641-mapping/1.0"},
            )
            with urllib.request.urlopen(request) as response, open(temporary, "wb") as output:
                shutil.copyfileobj(response, output, 1024 * 1024)
                output.flush()
                os.fsync(output.fileno())
            if not self._verify(temporary):
                raise MappingError("Downloaded NCBI gene annotation fails its checksum lock")
            os.replace(temporary, self.path)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
        return self.path

    def parse(self):
        mapping = {}
        with gzip.open(self.materialize(), "rt", encoding="utf-8", errors="strict", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                matches = self.ENSEMBL_PATTERN.findall(row["dbXrefs"])
                if len(set(matches)) != 1:
                    continue
                mapping[row["GeneID"]] = {
                    "ensembl_gene_id": matches[0],
                    "gene_symbol": row["Symbol"],
                    "gene_description": row["description"],
                }
        return mapping


class ReferencePreprocessor(object):
    """Build a log-library-normalized full-age reference and reference PCA."""

    def __init__(self, reference_root, validation_root, annotation, configuration):
        self.reference_root = reference_root
        self.validation_root = validation_root
        self.annotation = annotation
        self.configuration = configuration

    def metadata(self):
        rows = read_tsv(os.path.join(self.validation_root, "GSE94641_cell_metadata.tsv"))
        if len(rows) != int(self.configuration["expected_reference_cells"]):
            raise MappingError("Unexpected validated reference cell count")
        return rows

    def shared_gene_table(self, query_gene_ids):
        identifiers = read_tsv(os.path.join(self.validation_root, "GSE94641_gene_identifiers.tsv"))
        query_index = {gene_id: index for index, gene_id in enumerate(query_gene_ids)}
        seen_ensembl = Counter()
        candidates = []
        for row in identifiers:
            record = self.annotation.get(row["ID_REF"])
            if record is None or record["ensembl_gene_id"] not in query_index:
                continue
            seen_ensembl[record["ensembl_gene_id"]] += 1
            candidates.append({
                "entrez_gene_id": row["ID_REF"],
                "ensembl_gene_id": record["ensembl_gene_id"],
                "gene_symbol": record["gene_symbol"],
                "gene_description": record["gene_description"],
                "query_feature_index": query_index[record["ensembl_gene_id"]],
            })
        shared = [row for row in candidates if seen_ensembl[row["ensembl_gene_id"]] == 1]
        if len(shared) < int(self.configuration["minimum_shared_genes"]):
            raise MappingError("Insufficient one-to-one reference/query gene overlap: {}".format(len(shared)))
        return shared

    def expression(self, metadata, shared):
        manifest = read_tsv(os.path.join(self.validation_root, "GSE94641_expression_file_manifest.tsv"))
        by_gsm = {row["geo_accession"]: row for row in manifest}
        shared_index = {row["entrez_gene_id"]: index for index, row in enumerate(shared)}
        matrix = np.zeros((len(metadata), len(shared)), dtype=np.float32)
        library_sizes = np.zeros(len(metadata), dtype=np.float64)
        for cell_index, cell in enumerate(metadata):
            record = by_gsm.get(cell["geo_accession"])
            if record is None:
                raise MappingError("Reference metadata/expression mismatch")
            path = os.path.join(self.reference_root, record["relative_path"])
            with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
                for _ in range(3):
                    next(handle)
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    count = int(row["count"])
                    library_sizes[cell_index] += count
                    position = shared_index.get(row["ID_REF"])
                    if position is not None:
                        matrix[cell_index, position] = count
            if library_sizes[cell_index] <= 0:
                raise MappingError("Zero reference library size")
        target = float(self.configuration["normalization_target"])
        matrix = np.log1p(matrix * (target / library_sizes[:, None])).astype(np.float32)
        return matrix, library_sizes

    def fit_pca(self, matrix, shared):
        variances = matrix.var(axis=0, dtype=np.float64)
        count = min(int(self.configuration["selected_features"]), matrix.shape[1])
        selected = np.argsort(variances, kind="mergesort")[-count:][::-1]
        selected_matrix = matrix[:, selected].astype(np.float64)
        center = selected_matrix.mean(axis=0)
        centered = selected_matrix - center
        _, singular_values, vt = np.linalg.svd(centered, full_matrices=False)
        dimensions = min(int(self.configuration["pca_dimensions"]), vt.shape[0])
        loadings = vt[:dimensions].T
        coordinates = centered.dot(loadings)
        eigenvalues = (singular_values[:dimensions] ** 2) / float(centered.shape[0] - 1)
        explained = eigenvalues / np.sum((singular_values ** 2) / float(centered.shape[0] - 1))
        selected_rows = []
        for rank, index in enumerate(selected, 1):
            row = dict(shared[index])
            row.update({
                "selected_feature_rank": rank,
                "reference_log_variance": "{:.12g}".format(variances[index]),
                "query_selected_row_index": rank - 1,
            })
            selected_rows.append(row)
        return {
            "selected_shared_indices": selected,
            "query_feature_indices": np.asarray([shared[index]["query_feature_index"] for index in selected], dtype=np.int64),
            "center": center,
            "loadings": loadings,
            "coordinates": coordinates,
            "explained_variance_ratio": explained,
            "selected_rows": selected_rows,
        }


class QuerySampleLoader(object):
    """Load one Cell Ranger sample, normalize it, project it, and reuse its UMAP."""

    def __init__(self, sample_root, expected_features, normalization_target):
        self.sample_root = sample_root
        self.expected_features = expected_features
        self.normalization_target = normalization_target

    def paths(self, sample_id):
        root = os.path.join(self.sample_root, "per_sample_outs", sample_id)
        return {
            "matrix": os.path.join(root, "sample_filtered_feature_bc_matrix.h5"),
            "umap": os.path.join(root, "analysis", "umap", "gene_expression_2_components", "projection.csv"),
        }

    def feature_contract(self, sample_id):
        path = self.paths(sample_id)["matrix"]
        with h5py.File(path, "r") as handle:
            matrix = handle["matrix"]
            ids = decode(matrix["features"]["id"][:])
            names = decode(matrix["features"]["name"][:])
            types = decode(matrix["features"]["feature_type"][:])
        if len(ids) != self.expected_features or len(set(ids)) != len(ids):
            raise MappingError("Unexpected or duplicate query features in {}".format(sample_id))
        if set(types) != {"Gene Expression"}:
            raise MappingError("Unexpected feature types in {}".format(sample_id))
        return ids, names

    def project(self, sample_id, expected_gene_ids, pca):
        paths = self.paths(sample_id)
        with h5py.File(paths["matrix"], "r") as handle:
            group = handle["matrix"]
            ids = decode(group["features"]["id"][:])
            if ids != expected_gene_ids:
                raise MappingError("Query feature identity/order mismatch in {}".format(sample_id))
            barcodes = decode(group["barcodes"][:])
            shape = tuple(int(value) for value in group["shape"][:])
            matrix = sparse.csc_matrix(
                (group["data"][:], group["indices"][:], group["indptr"][:]),
                shape=shape,
                dtype=np.float32,
            )
        library_sizes = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
        if np.any(library_sizes <= 0):
            raise MappingError("Zero query library size in {}".format(sample_id))
        repeats = np.diff(matrix.indptr)
        matrix.data = np.log1p(
            matrix.data.astype(np.float64)
            * self.normalization_target
            / np.repeat(library_sizes, repeats)
        ).astype(np.float32)
        selected = matrix[pca["query_feature_indices"], :].transpose().tocsr()
        coordinates = selected.dot(pca["loadings"])
        coordinates = np.asarray(coordinates) - pca["center"].dot(pca["loadings"])
        del matrix, selected
        umap = self._umap(paths["umap"])
        if set(umap) != set(barcodes):
            raise MappingError("Vendor UMAP/barcode identity mismatch in {}".format(sample_id))
        umap_coordinates = np.asarray([umap[barcode] for barcode in barcodes], dtype=np.float64)
        return barcodes, coordinates, umap_coordinates, library_sizes, paths

    @staticmethod
    def _umap(path):
        result = {}
        with open(path, "r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                result[row["Barcode"]] = (float(row["UMAP-1"]), float(row["UMAP-2"]))
        return result


class KNNLabelTransfer(object):
    """Find fixed-k neighbors and transfer published labels by majority vote."""

    def __init__(self, reference_coordinates, metadata, indices, k, fields):
        self.coordinates = np.asarray(reference_coordinates[indices], dtype=np.float64)
        self.metadata = [metadata[index] for index in indices]
        self.k = min(int(k), len(indices))
        self.fields = fields
        self.reference_norm = np.sum(self.coordinates ** 2, axis=1)

    @staticmethod
    def composite(row):
        return "extendedphenotype={}|subtype={}|nonneuralcelltype={}".format(
            row["extendedphenotype"], row["subtype"], row["nonneuralcelltype"]
        )

    @staticmethod
    def _vote(values, distances):
        counts = Counter(values)
        maximum = max(counts.values())
        tied = set(value for value, count in counts.items() if count == maximum)
        winner = next(value for value, _ in zip(values, distances) if value in tied)
        return winner, maximum / float(len(values))

    def predict(self, query_coordinates, chunk_size=5000):
        outputs = []
        for start in range(0, query_coordinates.shape[0], chunk_size):
            query = np.asarray(query_coordinates[start:start + chunk_size], dtype=np.float64)
            distances = (
                np.sum(query ** 2, axis=1)[:, None]
                + self.reference_norm[None, :]
                - 2.0 * query.dot(self.coordinates.T)
            )
            np.maximum(distances, 0.0, out=distances)
            neighbors = np.argpartition(distances, self.k - 1, axis=1)[:, :self.k]
            for local_index in range(query.shape[0]):
                order = neighbors[local_index]
                order = order[np.argsort(distances[local_index, order], kind="mergesort")]
                neighbor_rows = [self.metadata[index] for index in order]
                neighbor_distances = np.sqrt(distances[local_index, order])
                result = {
                    "nearest_distance": float(neighbor_distances[0]),
                    "mean_k_distance": float(np.mean(neighbor_distances)),
                }
                for field in self.fields:
                    winner, confidence = self._vote(
                        [row[field] for row in neighbor_rows], neighbor_distances
                    )
                    result[field] = winner
                    result[field + "_confidence"] = confidence
                winner, confidence = self._vote(
                    [self.composite(row) for row in neighbor_rows], neighbor_distances
                )
                result["composite_label"] = winner
                result["composite_label_confidence"] = confidence
                if "age" not in self.fields:
                    age, age_confidence = self._vote(
                        [row["age"] for row in neighbor_rows], neighbor_distances
                    )
                    result["age"] = age
                    result["age_confidence"] = age_confidence
                outputs.append(result)
        return outputs


class MappingPlotter(object):
    """Create only the requested reference/query sanity figures."""

    BROAD_COLORS = {
        "postmitotic_immature_neuron": "#4C78A8",
        "proliferating_neural_progenitor": "#F58518",
        "not_assigned_neural_state": "#9E9E9E",
    }

    @staticmethod
    def broad_state(value):
        return {
            "maturing": "postmitotic_immature_neuron",
            "proliferating": "proliferating_neural_progenitor",
            "NA": "not_assigned_neural_state",
        }.get(value, "published_state_{}".format(value))

    def reference_pca(self, path, coordinates, metadata):
        fig, axes = plt.subplots(1, 3, figsize=(16, 5))
        panels = [("age", "Published age"), ("extendedphenotype", "Published broad state"), ("subtype", "Published subtype")]
        for axis, (field, title) in zip(axes, panels):
            values = [row[field] for row in metadata]
            labels = sorted(set(values))
            cmap = plt.get_cmap("tab20")
            for index, label in enumerate(labels):
                mask = np.asarray([value == label for value in values])
                axis.scatter(coordinates[mask, 0], coordinates[mask, 1], s=16, alpha=0.8, color=cmap(index % 20), label=label)
            axis.set_title(title)
            axis.set_xlabel("Reference PC1")
            axis.set_ylabel("Reference PC2")
            axis.legend(fontsize=7, frameon=False, markerscale=1.4)
        fig.suptitle("GSE94641 full E11.5-E17.5 reference PCA (published labels)")
        fig.tight_layout()
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)

    def query_umap_labels(self, path, plot_records):
        fig, axes = plt.subplots(3, 4, figsize=(18, 13))
        for axis, record in zip(axes.ravel(), plot_records):
            coordinates = record["umap"]
            labels = record["broad_state"]
            for label in sorted(set(labels)):
                mask = np.asarray(labels == label)
                axis.scatter(
                    coordinates[mask, 0], coordinates[mask, 1], s=0.45,
                    alpha=0.6, rasterized=True,
                    color=self.BROAD_COLORS.get(label, "#54A24B"), label=label,
                )
            axis.set_title(record["sample_id"])
            axis.set_xticks([]); axis.set_yticks([])
        handles, labels = axes.ravel()[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False)
        fig.suptitle("Existing per-sample Cell Ranger UMAP: E15.5-focused transferred broad state")
        fig.tight_layout(rect=(0, 0.04, 1, 0.97))
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)

    def query_umap_confidence(self, path, plot_records):
        fig, axes = plt.subplots(3, 4, figsize=(18, 13))
        scatter = None
        for axis, record in zip(axes.ravel(), plot_records):
            scatter = axis.scatter(
                record["umap"][:, 0], record["umap"][:, 1],
                c=record["confidence"], vmin=0, vmax=1, cmap="viridis",
                s=0.45, alpha=0.7, rasterized=True,
            )
            axis.set_title(record["sample_id"])
            axis.set_xticks([]); axis.set_yticks([])
        fig.colorbar(scatter, ax=axes.ravel().tolist(), fraction=0.02, pad=0.02, label="Composite-label neighbor-vote fraction")
        fig.suptitle("Existing per-sample Cell Ranger UMAP: E15.5 transfer confidence")
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)

    def proportions(self, path, sample_counts, sample_order):
        labels = sorted(set(label for counter in sample_counts.values() for label in counter))
        bottoms = np.zeros(len(sample_order))
        fig, axis = plt.subplots(figsize=(13, 6))
        for label in labels:
            values = np.asarray([
                sample_counts[sample][label] / float(sum(sample_counts[sample].values()))
                for sample in sample_order
            ])
            axis.bar(sample_order, values, bottom=bottoms, color=self.BROAD_COLORS.get(label, None), label=label)
            bottoms += values
        axis.set_ylabel("Fraction of cells")
        axis.set_ylim(0, 1)
        axis.tick_params(axis="x", rotation=45)
        axis.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
        axis.set_title("E15.5-focused GSE94641 transferred broad-state proportions")
        fig.tight_layout()
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)


class OutputPublisher(object):
    """Verify existing output or atomically publish a newly staged package."""

    @staticmethod
    def verify_existing(root):
        if not os.path.isdir(root):
            return False
        manifest = os.path.join(root, "output_manifest.tsv")
        if not os.path.isfile(manifest):
            raise MappingError("Existing mapping output has no manifest")
        for row in read_tsv(manifest):
            path = os.path.join(root, row["relative_path"])
            if (not os.path.isfile(path) or os.path.getsize(path) != int(row["bytes"])
                    or sha256_file(path) != row["sha256"]):
                raise MappingError("Existing mapping output fails manifest verification")
        return True

    @staticmethod
    def write_manifest(root):
        rows = []
        for directory, _, filenames in os.walk(root):
            for filename in sorted(filenames):
                path = os.path.join(directory, filename)
                relative = os.path.relpath(path, root)
                if relative == "output_manifest.tsv":
                    continue
                rows.append({
                    "relative_path": relative,
                    "bytes": os.path.getsize(path),
                    "sha256": sha256_file(path),
                })
        rows.sort(key=lambda row: row["relative_path"])
        write_tsv(os.path.join(root, "output_manifest.tsv"), ["relative_path", "bytes", "sha256"], rows)


class GSE94641LabelTransferWorkflow(object):
    """Coordinate full-age PCA, E15.5-primary transfer, plotting, and publication."""

    OUTPUT_HEADER = [
        "sample_id", "submitted_sample_name", "genotype", "sex", "design_group", "cell_barcode",
        "vendor_umap_1", "vendor_umap_2", "query_library_umi",
        "GSE94641_label", "GSE94641_label_confidence", "GSE94641_broad_state",
        "GSE94641_extendedphenotype", "GSE94641_extendedphenotype_confidence",
        "GSE94641_subtype", "GSE94641_subtype_confidence",
        "GSE94641_nonneuralcelltype", "GSE94641_nonneuralcelltype_confidence",
        "GSE94641_nearest_distance", "GSE94641_mean_k_distance", "GSE94641_reference",
        "GSE94641_all_age_label", "GSE94641_all_age_label_confidence",
        "GSE94641_all_age_extendedphenotype", "GSE94641_all_age_extendedphenotype_confidence",
        "GSE94641_all_age_subtype", "GSE94641_all_age_subtype_confidence",
        "GSE94641_all_age_nonneuralcelltype", "GSE94641_all_age_nonneuralcelltype_confidence",
        "GSE94641_all_age_neighbor_age", "GSE94641_all_age_neighbor_age_confidence",
        "GSE94641_all_age_nearest_distance", "GSE94641_all_age_mean_k_distance",
    ]

    def __init__(self, configuration, paths):
        self.configuration = configuration
        self.paths = paths
        self.plotter = MappingPlotter()

    def _validate_locks(self):
        manifest = os.path.join(self.paths["validation_root"], "output_manifest.tsv")
        if sha256_file(manifest) != self.configuration["reference_validation_manifest_sha256"]:
            raise MappingError("Reference validation manifest checksum mismatch")
        if sha256_file(self.paths["sample_key"]) != self.configuration["sample_key_sha256"]:
            raise MappingError("Sample key checksum mismatch")

    def _samples(self):
        with open(self.paths["sample_key"], "r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        rows.sort(key=lambda row: int(row["technical_sample_id"].rsplit("-", 1)[1]))
        if len(rows) != int(self.configuration["expected_query_samples"]):
            raise MappingError("Unexpected query sample count")
        return rows

    @staticmethod
    def _reference_rows(metadata, coordinates, library_sizes, primary_age):
        rows = []
        for index, cell in enumerate(metadata):
            row = OrderedDict(cell)
            row["reference_library_read_count"] = int(library_sizes[index])
            row["primary_e15_5_neighbor_pool"] = str(cell["age"] == primary_age).lower()
            for component in range(coordinates.shape[1]):
                row["reference_PC{}".format(component + 1)] = "{:.10g}".format(coordinates[index, component])
            rows.append(row)
        return rows

    def run(self):
        output_root = self.paths["output_root"]
        if OutputPublisher.verify_existing(output_root):
            print("verified existing mapping package: {}".format(output_root))
            return
        self._validate_locks()
        samples = self._samples()
        loader = QuerySampleLoader(
            self.paths["query_root"], int(self.configuration["expected_query_features"]),
            float(self.configuration["normalization_target"]),
        )
        query_gene_ids, query_gene_names = loader.feature_contract(samples[0]["technical_sample_id"])
        annotation_loader = LockedGeneAnnotation(self.paths["reference_root"], self.configuration["gene_annotation"])
        annotation = annotation_loader.parse()
        reference = ReferencePreprocessor(
            self.paths["reference_root"], self.paths["validation_root"], annotation, self.configuration
        )
        metadata = reference.metadata()
        shared = reference.shared_gene_table(query_gene_ids)
        matrix, reference_library_sizes = reference.expression(metadata, shared)
        pca = reference.fit_pca(matrix, shared)
        del matrix

        primary_indices = np.asarray([
            index for index, row in enumerate(metadata)
            if row["age"] == self.configuration["primary_reference_age"]
        ], dtype=np.int64)
        if len(primary_indices) != int(self.configuration["expected_primary_reference_cells"]):
            raise MappingError("Unexpected E15.5 primary reference count")
        all_indices = np.arange(len(metadata), dtype=np.int64)
        fields = list(self.configuration["transferred_fields"])
        primary_mapper = KNNLabelTransfer(
            pca["coordinates"], metadata, primary_indices,
            self.configuration["primary_k"], fields,
        )
        all_age_mapper = KNNLabelTransfer(
            pca["coordinates"], metadata, all_indices,
            self.configuration["all_age_context_k"], fields + ["age"],
        )

        parent = os.path.dirname(output_root)
        os.makedirs(parent, exist_ok=True)
        staging = tempfile.mkdtemp(prefix=".gse94641-transfer.", dir=parent)
        try:
            figures = os.path.join(staging, "figures")
            os.makedirs(figures)
            self.plotter.reference_pca(
                os.path.join(figures, "gse94641_reference_pca_published_labels.png"),
                pca["coordinates"], metadata,
            )
            reference_rows = self._reference_rows(
                metadata, pca["coordinates"], reference_library_sizes,
                self.configuration["primary_reference_age"],
            )
            reference_header = list(reference_rows[0].keys())
            write_tsv(os.path.join(staging, "gse94641_reference_cells_with_pca.tsv"), reference_header, reference_rows)
            write_tsv(
                os.path.join(staging, "gse94641_shared_selected_genes.tsv"),
                ["selected_feature_rank", "entrez_gene_id", "ensembl_gene_id", "gene_symbol", "gene_description", "query_feature_index", "query_selected_row_index", "reference_log_variance"],
                pca["selected_rows"],
            )
            write_tsv(
                os.path.join(staging, "gse94641_pca_variance.tsv"),
                ["component", "explained_variance_ratio"],
                ({"component": index + 1, "explained_variance_ratio": "{:.12g}".format(value)} for index, value in enumerate(pca["explained_variance_ratio"])),
            )

            combined_path = os.path.join(staging, "gse94641_query_cell_label_transfer.tsv")
            plot_records = []
            broad_counts = OrderedDict()
            summary_counts = Counter()
            input_manifest = []
            with open(combined_path, "w", encoding="utf-8", newline="") as combined_handle:
                writer = csv.DictWriter(combined_handle, delimiter="\t", fieldnames=self.OUTPUT_HEADER, lineterminator="\n")
                writer.writeheader()
                total_cells = 0
                for sample in samples:
                    sample_id = sample["technical_sample_id"]
                    print("mapping {}".format(sample_id))
                    current_ids, _ = loader.feature_contract(sample_id)
                    if current_ids != query_gene_ids:
                        raise MappingError("Query feature contract differs in {}".format(sample_id))
                    barcodes, query_pca, umap, library_sizes, input_paths = loader.project(sample_id, query_gene_ids, pca)
                    primary = primary_mapper.predict(query_pca)
                    all_age = all_age_mapper.predict(query_pca)
                    broad_states = np.asarray([
                        self.plotter.broad_state(prediction["extendedphenotype"])
                        for prediction in primary
                    ], dtype=object)
                    confidences = np.asarray([
                        prediction["composite_label_confidence"] for prediction in primary
                    ], dtype=np.float64)
                    plot_records.append({
                        "sample_id": sample_id, "umap": umap,
                        "broad_state": broad_states, "confidence": confidences,
                    })
                    broad_counts[sample_id] = Counter(broad_states)
                    for index, barcode in enumerate(barcodes):
                        first = primary[index]
                        context = all_age[index]
                        row = {
                            "sample_id": sample_id,
                            "submitted_sample_name": sample["submitted_sample_name"],
                            "genotype": sample["genotype"],
                            "sex": sample["sex"],
                            "design_group": sample["design_group"],
                            "cell_barcode": barcode,
                            "vendor_umap_1": "{:.10g}".format(umap[index, 0]),
                            "vendor_umap_2": "{:.10g}".format(umap[index, 1]),
                            "query_library_umi": int(library_sizes[index]),
                            "GSE94641_label": first["composite_label"],
                            "GSE94641_label_confidence": "{:.6f}".format(first["composite_label_confidence"]),
                            "GSE94641_broad_state": broad_states[index],
                            "GSE94641_extendedphenotype": first["extendedphenotype"],
                            "GSE94641_extendedphenotype_confidence": "{:.6f}".format(first["extendedphenotype_confidence"]),
                            "GSE94641_subtype": first["subtype"],
                            "GSE94641_subtype_confidence": "{:.6f}".format(first["subtype_confidence"]),
                            "GSE94641_nonneuralcelltype": first["nonneuralcelltype"],
                            "GSE94641_nonneuralcelltype_confidence": "{:.6f}".format(first["nonneuralcelltype_confidence"]),
                            "GSE94641_nearest_distance": "{:.10g}".format(first["nearest_distance"]),
                            "GSE94641_mean_k_distance": "{:.10g}".format(first["mean_k_distance"]),
                            "GSE94641_reference": "GSE94641_E15.5_k{}".format(self.configuration["primary_k"]),
                            "GSE94641_all_age_label": context["composite_label"],
                            "GSE94641_all_age_label_confidence": "{:.6f}".format(context["composite_label_confidence"]),
                            "GSE94641_all_age_extendedphenotype": context["extendedphenotype"],
                            "GSE94641_all_age_extendedphenotype_confidence": "{:.6f}".format(context["extendedphenotype_confidence"]),
                            "GSE94641_all_age_subtype": context["subtype"],
                            "GSE94641_all_age_subtype_confidence": "{:.6f}".format(context["subtype_confidence"]),
                            "GSE94641_all_age_nonneuralcelltype": context["nonneuralcelltype"],
                            "GSE94641_all_age_nonneuralcelltype_confidence": "{:.6f}".format(context["nonneuralcelltype_confidence"]),
                            "GSE94641_all_age_neighbor_age": context["age"],
                            "GSE94641_all_age_neighbor_age_confidence": "{:.6f}".format(context["age_confidence"]),
                            "GSE94641_all_age_nearest_distance": "{:.10g}".format(context["nearest_distance"]),
                            "GSE94641_all_age_mean_k_distance": "{:.10g}".format(context["mean_k_distance"]),
                        }
                        writer.writerow(row)
                        summary_counts[(sample_id, "extendedphenotype", first["extendedphenotype"])] += 1
                        summary_counts[(sample_id, "subtype", first["subtype"])] += 1
                        summary_counts[(sample_id, "nonneuralcelltype", first["nonneuralcelltype"])] += 1
                        summary_counts[(sample_id, "broad_state", broad_states[index])] += 1
                    total_cells += len(barcodes)
                    for role, path in input_paths.items():
                        input_manifest.append({
                            "sample_id": sample_id, "input_role": role,
                            "path": path, "bytes": os.path.getsize(path),
                            "sha256": sha256_file(path),
                        })
                    del query_pca, primary, all_age
            if total_cells != int(self.configuration["expected_query_cells"]):
                raise MappingError("Unexpected total mapped query cells: {}".format(total_cells))

            summary_rows = []
            for (sample_id, field, label), count in sorted(summary_counts.items()):
                denominator = sum(value for (sid, fld, _), value in summary_counts.items() if sid == sample_id and fld == field)
                summary_rows.append({
                    "sample_id": sample_id, "transferred_field": field,
                    "transferred_label": label, "cells": count,
                    "fraction_of_sample": "{:.8f}".format(count / float(denominator)),
                })
            write_tsv(
                os.path.join(staging, "gse94641_transferred_label_summary_by_sample.tsv"),
                ["sample_id", "transferred_field", "transferred_label", "cells", "fraction_of_sample"],
                summary_rows,
            )
            write_tsv(
                os.path.join(staging, "query_input_manifest.tsv"),
                ["sample_id", "input_role", "path", "bytes", "sha256"], input_manifest,
            )
            configuration_rows = [
                {"parameter": key, "value": json.dumps(value, sort_keys=True) if isinstance(value, (dict, list)) else value}
                for key, value in self.configuration.items()
            ]
            configuration_rows.extend([
                {"parameter": "actual_shared_one_to_one_genes", "value": len(shared)},
                {"parameter": "actual_selected_features", "value": len(pca["selected_rows"])},
                {"parameter": "actual_query_cells", "value": total_cells},
                {"parameter": "pcdh19_classification_features_used", "value": "false"},
            ])
            write_tsv(os.path.join(staging, "mapping_configuration.tsv"), ["parameter", "value"], configuration_rows)
            write_tsv(
                os.path.join(staging, "software_environment.tsv"),
                ["component", "version"],
                [
                    {"component": "python", "version": platform.python_version()},
                    {"component": "numpy", "version": np.__version__},
                    {"component": "scipy", "version": scipy.__version__},
                    {"component": "h5py", "version": h5py.__version__},
                    {"component": "matplotlib", "version": matplotlib.__version__},
                ],
            )
            self.plotter.query_umap_labels(
                os.path.join(figures, "query_vendor_umap_gse94641_e15_5_broad_state.png"), plot_records
            )
            self.plotter.query_umap_confidence(
                os.path.join(figures, "query_vendor_umap_gse94641_e15_5_confidence.png"), plot_records
            )
            self.plotter.proportions(
                os.path.join(figures, "gse94641_e15_5_transferred_broad_state_proportions_by_sample.png"),
                broad_counts, [sample["technical_sample_id"] for sample in samples],
            )
            OutputPublisher.write_manifest(staging)
            os.replace(staging, output_root)
            print("published mapping package: {}".format(output_root))
        finally:
            if os.path.isdir(staging):
                shutil.rmtree(staging)
