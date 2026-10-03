# Environments & Installations

## Primary Environment

- **Manager**: Python virtual environment (`.venv/` in this checkout)
- **Python version**: 3.12.13 in this checkout; package requires Python >=3.10
- **Workflow**: Snakemake, invoked by the `batchdx-cyto` CLI

### Setup from scratch

```bash
python -m venv .venv
.venv/bin/python -m pip install -e '.[dev]'
```

## Dependencies

Python dependencies are declared in `pyproject.toml`. The optional metric dependencies are in its `metrics` extra. Pinned upstream repositories are recorded in `config/upstream_lock.yaml`.

The R backends use the isolated environment built by `scripts/setup_r_env.sh`; that script pins R 4.4 and installs its required packages. See `docs/decisions/0004-dependency-isolation.md` before rebuilding it.

## System Dependencies

`scripts/setup_r_env.sh` uses micromamba for the R environment. The Snakemake workflow also has environment definitions under `workflow/envs/`.
