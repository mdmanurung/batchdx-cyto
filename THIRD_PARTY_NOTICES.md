# Third-party notices

batchdx-cyto is MIT licensed. It **does not redistribute** any of the components below;
each is fetched or installed by the user and keeps its own licence. No source from any
of them is copied into this package.

## R components invoked as separate programs

| Component | Licence | How it is used | Redistributed here? |
|---|---|---|---|
| [cyCombine](https://github.com/biosurf/cyCombine) @ `fcd42dbb` | MIT | called from a separate `Rscript` process | no |
| [CytoNorm](https://github.com/saeyslab/CytoNorm) @ `278f352d` | GPL (≥ 2) | called from a separate `Rscript` process; fitted on the training split's controls, applied to a held-out split | no |
| [CytoBatchFlagR](https://github.com/AndorfLab/CytoBatchFlagR) @ `ad21e3b1` | **GPL-3.0** | its scripts are `source()`d **unmodified** from a checkout the user supplies, in a separate `Rscript` process | no |

The CytoBatchFlagR case is the one that constrains design rather than only paperwork.
GPL-3 is copyleft: copying or refactoring those scripts into this MIT package would make
the result a derivative work. The decision recorded in
`docs/decisions/0003-cytobatchflagr-integration.md` is therefore to integrate at arm's
length and to leave the obvious engineering improvement — turning the scripts into a
proper R package — undone here. If that becomes necessary, the correct route is a
separate GPL-3 package.

`workflow/scripts/cytobatchflagr_detect.R` is original work by this project: it loads
the upstream file and calls its functions. It contains no upstream code.

## Referenced, and partly executed

| Component | Licence | Relationship |
|---|---|---|
| [OpenProblems cytometry batch integration](https://github.com/openproblems-bio/task_cyto_batch_integration) @ `37151479` | MIT | seven of its nine metric components are re-executed from the pinned Viash `helper.py` / `helper.R` and registered in this package's `compatibility` namespace with `upstream_verified=True`: `emd_mean_ct_horiz`, `emd_mean_ct_vert`, `ratio_consistent_peaks`, `average_batch_r2_ct`, `iLisi`, `cLisi` and `flowsom_mean_mapping_similarity`. The other two remain deferred with a stated obstacle. See `docs/blocked.md` B-007 (RESOLVED 2026-09-29) for what that resolution did and did not settle, and `config/upstream_lock.yaml` for the pinned commit. **No upstream source is copied into this package** — each helper is imported from the user's checkout by path and its output compared against ours under a recorded tolerance (`docs/compatibility_matrix.md` lists the tolerance per component and the test that enforces it) |
| [SPiSCy](https://github.com/emilier17/spiscy) @ `a6e55c8f` | MIT | reference reading for workflow layout; no code or contract used |

## Python dependencies

Declared in `pyproject.toml` and installed by the user's package manager; each retains its
own licence. Core: numpy, pandas, scipy, anndata, h5py, pydantic, typer, jinja2, PyYAML,
flowio, pyarrow, matplotlib.

Optional extras:

| Extra | Packages | Why it is optional |
|---|---|---|
| `workflow` | snakemake | only needed to execute the workflow; the package's stages import none of it |
| `metrics` | scikit-learn, scib-metrics | three of the six verified `compatibility` components import these from the *pinned upstream helper*, so without them those rows cannot run their own parity comparison |
| `tests`, `dev` | pytest, pytest-json-report, ruff (plus everything above) | the test and lint toolchain; `dev` is the union |

Run `pip-licenses` or equivalent against your resolved environment for the exact set and
versions you deployed.
