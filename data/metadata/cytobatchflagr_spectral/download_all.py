#!/usr/bin/env python3
"""Download every file in the frozen Zenodo manifest to a directory outside Git."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import tempfile
import urllib.request
from pathlib import Path

MANIFEST = Path(__file__).with_name("source_files.tsv")
REPO_ROOT = Path(__file__).resolve().parents[3]


def md5_and_size(path: Path) -> tuple[str, int]:
    digest = hashlib.md5()  # noqa: S324 - comparison with the source's MD5, not security
    size = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
            size += len(chunk)
    return digest.hexdigest(), size


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dest", type=Path, required=True, help="directory outside Git")
    args = parser.parse_args()
    dest = args.dest.resolve()
    if dest == REPO_ROOT or REPO_ROOT in dest.parents:
        parser.error("destination must be outside this Git repository")
    dest.mkdir(parents=True, exist_ok=True)

    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 74 or sum(row["kind"] == "fcs" for row in rows) != 72:
        raise ValueError("frozen manifest must contain all 72 FCS and 2 attachments")

    for row in rows:
        name = row["file"]
        if Path(name).name != name:
            raise ValueError(f"unsafe manifest filename: {name!r}")
        target = dest / name
        expected = row["md5"], int(row["bytes"])
        if target.exists():
            if md5_and_size(target) != expected:
                raise ValueError(f"existing file does not match frozen manifest: {target}")
            print(f"verified {name}", flush=True)
            continue
        fd, temporary = tempfile.mkstemp(prefix=f".{name}.", suffix=".part", dir=dest)
        try:
            with os.fdopen(fd, "wb") as output, urllib.request.urlopen(row["url"], timeout=120) as response:
                for chunk in iter(lambda: response.read(1024 * 1024), b""):
                    output.write(chunk)
            temp_path = Path(temporary)
            if md5_and_size(temp_path) != expected:
                raise ValueError(f"download does not match frozen manifest: {name}")
            os.replace(temp_path, target)
            print(f"downloaded and verified {name}", flush=True)
        finally:
            Path(temporary).unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
