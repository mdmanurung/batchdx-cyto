# batchdx-cyto

Batch-effect **diagnosis**, **correction** and **evaluation** for mass cytometry,
conventional flow cytometry and spectral flow cytometry, orchestrated by Snakemake.

Built for cross-cohort immune studies, where biological differences and technical
batches are often correlated and where the honest answer is frequently *"this contrast
cannot be corrected"*.

> **Status: alpha, partially blocked.** The Python core, the design audit, the native
> diagnostics, the common-partition metrics, the decision layer, the report and the
> Snakemake workflow run end to end and are tested. The **R-backed adapters (cyCombine,
> CytoNorm, CytoBatchFlagR) cannot be executed** in the environment this was developed
> in, because CRAN and Bioconductor were unreachable. Those gates are recorded as
> **BLOCKED**, not passed. See [`docs/implementation_status.md`](docs/implementation_status.md)
> (generated from an executed test report) and [`docs/blocked.md`](docs/blocked.md).

---

## What it does

1. **Diagnoses** batch-associated differences and asks whether a *technical*
   interpretation is actually supported by the design.
2. **Runs eligible correction methods independently**, always including an identity
   candidate, because "no correction" is a valid result.
3. **Evaluates** technical improvement and, separately, what the correction did to the
   data's structure.
4. **Reports**: one offline HTML file, structured JSON/CSV records, and corrected
   datasets with an explicit scale.

Three things it deliberately refuses to do:

* **It will not turn association into attribution.** A batch contrast is classified as
  `reference_supported`, `biologically_stratified_conditional`, `descriptive_only` or
  `non_identifiable`, and perfect confounding blocks automatic recommendations for the
  affected contrast. Adding the confounded factor as a covariate does not fix it, and
  the audit says so.
* **It will not invent evidence.** An unavailable dependency, an invalid comparison, an
  absent control or a failed method gets an explicit status and a reason. A missing
  metric is never a zero, and a near-zero denominator is `undefined`, not an epsilon.
* **It will not force a winner.** Below a declared design floor it abstains and produces
  descriptive output only. With uncalibrated thresholds it produces *review*, never
  automatic approval.

## Install

```bash
pip install -e '.[dev]'          # from a checkout
batchdx-cyto info                # what this build actually ships
```

Python ≥ 3.10. The R backends additionally need an R environment; see
[`containers/Dockerfile.r-backends`](containers/Dockerfile.r-backends) (**not built or
tested in this build**) and `docs/decisions/0004-dependency-isolation.md`.

## Five-minute tour

```bash
batchdx-cyto init my_project --example paired_reference
cd my_project

batchdx-cyto validate  --config config/config.yaml
batchdx-cyto diagnose  --config config/config.yaml --outdir results/run1
batchdx-cyto benchmark --config config/config.yaml --outdir results/run1
open results/run1/report/report.html
```

`init` writes everything the configuration references — synthetic data, the three
manifests and the gating threshold table — so the example runs with no further setup
and no developer-specific paths. Three examples ship: `paired_reference`,
`no_reference` and `confounded` (see [`config/README.md`](config/README.md)).

### Running Snakemake directly

The CLI delegates to Snakemake; it is not a second workflow engine. Anything the CLI
does you can do yourself:

```bash
snakemake --snakefile "$(batchdx-cyto info --json | python -c 'import json,sys; print(json.load(sys.stdin)["workflow"])')" \
    --cores 4 \
    --config config=config/config.yaml outdir=results/run1 base_dir=. mode=benchmark
```

`mode` is one of `diagnose`, `benchmark`, `correct` (which also needs
`candidate=<id>`) or `report`.

## The three commands

| Command | What it does |
|---|---|
| `diagnose` | design audit, detectors, diagnosis-only report. No correction runs. |
| `benchmark` | every configured candidate, scored on one frozen partition, then assessed under a published policy. |
| `correct` | applies **one explicitly named** candidate. A correction is never chosen for you. |

Plus `init`, `validate`, `report` (regenerates the HTML without rerunning anything) and
`audit` (runs the Stage 0 upstream audit).

## How the scientific guarantees are enforced

| Guarantee | Where it lives |
|---|---|
| The identity candidate is always benchmarked | `config.py` rejects a configuration without a `backend: none` candidate |
| The baseline is immutable | `Dataset.expression()` returns copies; `assert_baseline_unchanged()` hashes it |
| Evaluation labels are frozen before any candidate runs | `sampling/partitions.py`; the tier and confidence ride on every metric record |
| Clustering labels can never alone support a recommendation | `decisions/policy.py` |
| Correction fitting is separated from held-out evaluation | not yet — see Stage 4 in the status file |
| Cells are not independent replicates | the design audit and the composition test both work at sample level |
| Perfect confounding is not correctable | `design/audit.py::_perfectly_confounded`, and the eligibility entry is marked non-overridable |
| A failure never returns unchanged data as success | `adapters/base.py::failure_result`, tested in Stage 6 |
| A locally written metric never claims to be an upstream component | `metrics/registry.py` namespaces, asserted by `scripts/audit_upstream.py` |

## Upstream tools

Pinned by commit in [`config/upstream_lock.yaml`](config/upstream_lock.yaml), audited by
two executable scripts, and interpreted in [`docs/upstream_audit.md`](docs/upstream_audit.md).

| Tool | Commit | Licence | Integration |
|---|---|---|---|
| [cyCombine](https://github.com/biosurf/cyCombine) 0.3.0 | `fcd42dbb` | MIT | subprocess; ComBat route pinned, the count model refused |
| [CytoNorm](https://github.com/saeyslab/CytoNorm) 2.0.12 | `278f352d` | GPL (≥2) | declared, **deferred** |
| [CytoBatchFlagR](https://github.com/AndorfLab/CytoBatchFlagR) (was `cytoFlagR`) | `ad21e3b1` | **GPL-3.0** | scripts `source()`d unmodified at arm's length; nothing vendored |
| [OpenProblems cyto batch integration](https://github.com/openproblems-bio/task_cyto_batch_integration) | `37151479` | MIT | reference only; metrics deferred |
| [SPiSCy](https://github.com/emilier17/spiscy) | `a6e55c8f` | MIT | reference reading only |

```bash
scripts/fetch_upstream.sh .upstream
batchdx-cyto audit --upstream-dir .upstream
```

## Tests and status

```bash
pytest -q                                              # everything
pytest -q -m "not slow"                                # skip the Snakemake runs
pytest --json-report --json-report-file=build/report.json
python scripts/generate_status.py --report build/report.json
```

`docs/implementation_status.md` is **generated** from that report against the gate
manifest in `scripts/gates.yaml`. A status line without an executed test cannot exist,
and a gate with a recorded external blocker never reads as passed.

## Documentation

* [`docs/implementation_status.md`](docs/implementation_status.md) — generated status per stage
* [`docs/usage.md`](docs/usage.md) — inputs, reference requirements, confounding, interpretation, limits
* [`docs/upstream_audit.md`](docs/upstream_audit.md) — what the pinned sources actually say
* [`docs/validation_results.md`](docs/validation_results.md) — what the scientific suite establishes, and what it does not
* [`docs/blocked.md`](docs/blocked.md) — every blocked item with its command and error
* [`docs/decisions/`](docs/decisions/) — architecture and licence decisions

## Repository history

This repository previously held a cookiecutter Python project template. Those files are
retained untouched (`cookiecutter.json`, `hooks/`, `pixi.toml`,
`{{ cookiecutter.__repo_name }}/`, `tests/test_template_bake.py`,
`tests/test_post_gen_hook.py`) and are unrelated to the package documented here.
`pyproject.toml` scopes pytest to `tests/unit`, `tests/integration` and
`tests/scientific`.

## Licence

MIT. Third-party components keep their own licences — see
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md). Note that CytoBatchFlagR is GPL-3.0
and CytoNorm is GPL (≥2); neither is redistributed here, and neither is copied into this
package.
