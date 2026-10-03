#!/usr/bin/env bash
# Clone every upstream resource at the commit pinned in config/upstream_lock.yaml.
#
# The checkouts are what the Stage 0 audit reads and what the CytoBatchFlagR adapter
# source()s.  Nothing here installs an R package: cyCombine and CytoNorm pull
# dependencies from CRAN and Bioconductor and need the container image in containers/.

set -euo pipefail

TARGET="${1:-.upstream}"
LOCK="${2:-config/upstream_lock.yaml}"

if [[ ! -f "$LOCK" ]]; then
  echo "lock manifest not found: $LOCK" >&2
  exit 2
fi

mkdir -p "$TARGET"

fetch () {
  local name="$1" url="$2" commit="$3"
  local dir="$TARGET/$name"
  if [[ -d "$dir/.git" ]]; then
    echo "== $name: already present, fetching $commit"
    git -C "$dir" fetch --quiet origin "$commit" 2>/dev/null || git -C "$dir" fetch --quiet origin
  else
    echo "== $name: cloning $url"
    git clone --quiet "$url" "$dir"
  fi
  git -C "$dir" checkout --quiet "$commit"
  echo "   HEAD = $(git -C "$dir" rev-parse HEAD)"
}

# Kept in sync with config/upstream_lock.yaml by
# tests/integration/test_stage0_audit.py::test_fetch_script_pins_the_same_commits_as_the_lock.
fetch cyCombine        https://github.com/biosurf/cyCombine.git        fcd42dbbe92bfa6bffa6094bf9b4de536e71c249
fetch CytoNorm         https://github.com/saeyslab/CytoNorm.git        278f352d8edffd7a5d3e6b0b76c402d991a1e7a9
fetch CytoBatchFlagR   https://github.com/AndorfLab/CytoBatchFlagR.git ad21e3b104126783044351eb11fa24622a0da2f7
fetch task_cyto_batch_integration \
      https://github.com/openproblems-bio/task_cyto_batch_integration.git \
      37151479aaa89edf37ae23068de8fbb41a1b7bb4
fetch spiscy           https://github.com/emilier17/spiscy.git         a6e55c8f722112bf01d5364996b1867844efbd02

echo
echo "checkouts are in $TARGET"
echo "next: python scripts/audit_upstream.py --upstream-dir $TARGET --out docs/audit/audit_python.json"
echo "      Rscript scripts/audit_upstream.R docs/audit/audit_r.json $TARGET"
