"""Ingestion, the immutable preprocessed baseline, and safe export."""

from batchdx.cyto.io.dataset import (
    Dataset,
    DatasetValidationError,
    load_dataset,
    write_candidate_output,
)
from batchdx.cyto.io.hashing import file_sha256, frame_hash, matrix_hash
from batchdx.cyto.io.manifests import (
    ManifestError,
    read_marker_table,
    read_processing_manifest,
    read_sample_table,
)

__all__ = [
    "Dataset",
    "DatasetValidationError",
    "ManifestError",
    "file_sha256",
    "frame_hash",
    "load_dataset",
    "matrix_hash",
    "read_marker_table",
    "read_processing_manifest",
    "read_sample_table",
    "write_candidate_output",
]
