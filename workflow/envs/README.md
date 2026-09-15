# Environments

Fine-grained environment files live here so a rule pulls only what it needs.

**Status in this build:** the Python side is pinned through `pyproject.toml` and the
lock file of whatever installer you use. The R side is *not* reproducibly pinnable from
a conda environment file alone, because cyCombine and CytoNorm are distributed through
GitHub rather than CRAN or Bioconductor — see `docs/decisions/0004-dependency-isolation.md`
and `containers/`.

`scripts/fetch_upstream.sh` materialises the pinned commits; the container image is what
makes the transitive CRAN and Bioconductor dependencies reproducible.
