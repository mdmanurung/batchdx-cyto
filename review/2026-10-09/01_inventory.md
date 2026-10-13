# 01 Inventory (generated 2026-10-09)

## Source modules (src/batchdx/cyto)

| module | files | LOC |
|---|---|---|
| (top-level) | 8 | 1960 |
| adapters | 7 | 1480 |
| decisions | 2 | 411 |
| design | 2 | 503 |
| diagnostics | 6 | 933 |
| io | 7 | 1227 |
| metrics | 8 | 3487 |
| preprocess | 5 | 1669 |
| reporting | 4 | 749 |
| sampling | 3 | 1201 |
| schemas | 3 | 636 |

## R / workflow scripts

| file | LOC |
|---|---|
| workflow/Snakefile | 271 |
| workflow/scripts/cycombine_correct.R | 111 |
| workflow/scripts/cycombine_detect.R | 148 |
| workflow/scripts/cycombine_evaluate.R | 127 |
| workflow/scripts/cytobatchflagr_detect.R | 210 |
| workflow/scripts/cytonorm_correct.R | 228 |
| workflow/scripts/lib/binmat.R | 76 |
| workflow/scripts/peacoqc_qc.R | 235 |

## Tests per directory

| dir | test files | test functions |
|---|---|---|
| tests/unit | 7 | 125 |
| tests/integration | 11 | 134 |
| tests/scientific | 6 | 105 |
| tests/testthat | 1 | 13 |

## Gate ↔ test map (scripts/gates.yaml)

| stage | tier | title | #tests | blocked_reason |
|---|---|---|---|---|
| X | A | Cross-stage contracts (typed, versioned, round-tripped across the R boundary) | 10 |  |
| 0 | A | Audit upstream resources and establish contracts | 8 |  |
| 1 | A | Runnable package and minimal Snakemake workflow | 32 |  |
| 2 | A | Ingestion, immutable data, preprocessing machinery (not yet in a run) | 64 |  |
| 3 | A | Design audit and correction eligibility | 12 |  |
| 4 | B | Freeze evaluation populations, sampling and splits | 26 |  |
| 5 | A | Detection adapters and diagnosis report | 19 |  |
| 6 | A | Correction backends and output validation | 33 |  |
| 7 | A | Native and common comparative metrics | 52 |  |
| 8 | B | Biological preservation and uncertainty | 18 |  |
| 9 | B | Candidate assessment and recommendation rules | 13 |  |
| 10 | A | Reporting and safe exports | 34 |  |
| 11 | C | Snakemake execution, performance, deployment | 11 |  |
| 12 | C | Scientific validation, documentation, release readiness | 4 |  |
