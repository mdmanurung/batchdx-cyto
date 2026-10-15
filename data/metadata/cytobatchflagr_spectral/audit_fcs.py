#!/usr/bin/env python3
"""Audit every downloaded spectral-control FCS file against the frozen record."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import flowio

MANIFEST = Path(__file__).with_name("source_files.tsv")


def hashes(path: Path) -> tuple[str, str, int]:
    md5 = hashlib.md5()  # noqa: S324 - source metadata supplies MD5
    sha256 = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            md5.update(block)
            sha256.update(block)
            size += len(block)
    return md5.hexdigest(), sha256.hexdigest(), size


def write_tsv(path: Path, header: tuple[str, ...], rows: list[tuple]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    source = args.data_dir.resolve()
    out = args.out_dir.resolve()
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if row["kind"] == "fcs"]
    if len(rows) != 72:
        raise ValueError(f"expected all 72 FCS source rows, found {len(rows)}")

    inventory: list[tuple] = []
    channel_layouts: dict[tuple[tuple[str, str], ...], int] = Counter()
    file_layouts: list[tuple[str, tuple[tuple[str, str], ...]]] = []
    for row in rows:
        path = source / row["file"]
        if not path.is_file():
            raise FileNotFoundError(f"missing source FCS: {path}")
        md5, sha256, size = hashes(path)
        if (md5, size) != (row["md5"], int(row["bytes"])):
            raise ValueError(f"source FCS differs from frozen manifest: {path}")
        data = flowio.FlowData(str(path))
        layout = tuple(
            (
                str((data.channels.get(i) or data.channels.get(str(i)) or {}).get("pnn", "")),
                str((data.channels.get(i) or data.channels.get(str(i)) or {}).get("pns", "") or ""),
            )
            for i in range(1, int(data.channel_count) + 1)
        )
        if not all(name for name, _ in layout) or len({name for name, _ in layout}) != len(layout):
            raise ValueError(f"missing or duplicate channel names: {path}")
        if int(data.event_count) <= 0:
            raise ValueError(f"empty FCS: {path}")
        channel_layouts[layout] += 1
        file_layouts.append((row["file"], layout))
        inventory.append((row["file"], row["batch"], row["control"], row["sample_id"],
                          int(data.event_count), int(data.channel_count), sha256, md5))
        print(f"audited {row['file']}", flush=True)

    counts = Counter((batch, control) for _, batch, control, *_ in inventory)
    batches = sorted({batch for batch, _ in counts})
    controls = sorted({control for _, control in counts})
    complete_design = (
        len(batches) == 18 and controls == ["tcs3", "tcs4", "tcu1", "tcu2"]
        and all(counts[(batch, control)] == 1 for batch in batches for control in controls)
    )
    one_layout = len(channel_layouts) == 1

    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "fcs_inventory.tsv",
              ("file", "batch", "control", "sample_id", "events", "channels", "sha256", "md5"),
              inventory)
    ordered_layouts = sorted(channel_layouts.items(), key=lambda item: (-item[1], item[0]))
    layout_ids = {layout: i for i, (layout, _) in enumerate(ordered_layouts, 1)}
    write_tsv(out / "file_layouts.tsv", ("file", "layout"),
              [(file, layout_ids[layout]) for file, layout in file_layouts])
    largest_layout = ordered_layouts[0][0]
    eligible_files = {file for file, layout in file_layouts if layout == largest_layout}
    eligible_rows = [row for row in inventory if row[0] in eligible_files]
    eligible_batches = sorted({row[1] for row in eligible_rows})
    eligible_counts = Counter((row[1], row[2]) for row in eligible_rows)
    subset_complete = (
        len(eligible_batches) >= 2
        and all(eligible_counts[(batch, control)] == 1
                for batch in eligible_batches for control in controls)
    )
    write_tsv(out / "technical_control_subset.tsv",
              ("file", "batch", "control", "sample_id", "events", "channels", "sha256", "md5"),
              eligible_rows if subset_complete else [])
    channel_rows = [
        (layout_number, i, name, label, n_files)
        for layout_number, (layout, n_files) in enumerate(ordered_layouts, 1)
        for i, (name, label) in enumerate(layout, 1)
    ]
    write_tsv(out / "channel_layouts.tsv",
              ("layout", "channel_index", "pnn", "pns", "n_files"), channel_rows)
    report = {
        "source_record": "10.5281/zenodo.15388817",
        "files_audited": len(inventory),
        "events_total": sum(row[4] for row in inventory),
        "batches": batches,
        "controls": controls,
        "files_per_batch_control": {f"{b}:{c}": counts[(b, c)] for b in batches for c in controls},
        "channel_layout_count": len(channel_layouts),
        "technical_control_design": "eligible_for_review" if complete_design and one_layout else "ineligible_as_is",
        "technical_control_reason": (
            "one measured file for each of four controls in each of 18 batches; one channel layout"
            if complete_design and one_layout else
            "the complete 18-batch four-control grid has multiple channel layouts"
            if complete_design else "incomplete batch/control grid"
        ),
        "homogeneous_subset_design": "eligible_for_review" if subset_complete else "ineligible",
        "homogeneous_subset_layout": layout_ids[largest_layout],
        "homogeneous_subset_batches": eligible_batches,
        "homogeneous_subset_files": len(eligible_rows) if subset_complete else 0,
        "excluded_from_subset_files": len(inventory) - len(eligible_rows),
        "donor_disjoint_design": "unestablished",
        "donor_disjoint_reason": "published Table_2.csv has no donor identity",
        "biological_preservation_design": "unestablished",
        "biological_preservation_reason": "no donor or outcome mapping established from the source metadata",
    }
    (out / "design_audit.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in ("files_audited", "events_total",
                                               "channel_layout_count", "technical_control_design",
                                               "homogeneous_subset_files")},
                     indent=2))
    return 0 if (complete_design and one_layout) or subset_complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
