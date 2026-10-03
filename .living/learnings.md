# Learnings

Append-only log of gotchas, surprises, and insights.

**Entry template:** copy from `skills/core/templates/learning-entry.md` (includes Category, What happened, Why it matters, Resolution, Tags fields). The `**Tags**:` line is consumed by `generate_index.py --summary-heuristic` to build the cluster summary in INDEX.md — use them.

### 2026-10-03 — Full package tests need the isolated R runtime and writable caches

- **Category:** test environment
- **What happened:** A sandboxed full pytest run reported backend failures when `Rscript` was absent from `PATH`; Snakemake also tried to write below the read-only home cache. With `renv/bin` on `PATH` and caches in `/tmp`, cyCombine still needed a BiocParallel localhost port that the sandbox denied. The isolated cyCombine test passed with local socket access.
- **Why it matters:** These failures can be mistaken for changes in correction algorithms or workflow behavior.
- **Resolution:** Run integration checks with the repository's `renv/bin`, writable `XDG_CACHE_HOME`, `MPLCONFIGDIR` and `NUMBA_CACHE_DIR`, and local socket access for BiocParallel.
- **Tags:** testing, R, Snakemake, sandbox
- **mitigation_type:** convention
- **structural_mitigation_candidate:** A test launcher could set the local runtime and caches when this environment recurs; it cannot grant socket access inside a restricted sandbox.
