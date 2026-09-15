# Third-party notices

batchdx-cyto is MIT licensed. It **does not redistribute** any of the components below;
each is fetched or installed by the user and keeps its own licence. No source from any
of them is copied into this package.

## R components invoked as separate programs

| Component | Licence | How it is used | Redistributed here? |
|---|---|---|---|
| [cyCombine](https://github.com/biosurf/cyCombine) @ `fcd42dbb` | MIT | called from a separate `Rscript` process | no |
| [CytoNorm](https://github.com/saeyslab/CytoNorm) @ `278f352d` | GPL (≥ 2) | would be called from a separate process; the backend is deferred | no |
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

## Referenced but not invoked

| Component | Licence | Relationship |
|---|---|---|
| [OpenProblems cytometry batch integration](https://github.com/openproblems-bio/task_cyto_batch_integration) @ `37151479` | MIT | its metric inventory was read to reconcile names; no implementation was extracted (see `docs/blocked.md` B-007) |
| [SPiSCy](https://github.com/emilier17/spiscy) @ `a6e55c8f` | MIT | reference reading for workflow layout; no code or contract used |

## Python dependencies

Declared in `pyproject.toml` and installed by the user's package manager: numpy, pandas,
scipy, anndata, h5py, pydantic, typer, jinja2, PyYAML, flowio, pyarrow, matplotlib and
(optionally) snakemake. Each retains its own licence; run `pip-licenses` or equivalent
against your resolved environment for the exact set and versions you deployed.
