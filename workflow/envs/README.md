# Environments

**Two environment files live here, and neither is consumed by the workflow.** They are
human-run recipes: something reads them, but it is a person, not Snakemake.

| File | What it creates | Consumed by |
|---|---|---|
| `python.yaml` | every Python-side dependency the rules in `workflow/Snakefile` need | nothing; you run `conda env create` yourself |
| `openproblems.yaml` | the same, plus `scikit-learn` / `scib-metrics`, which three of the six verified `compatibility` components need to run their own parity comparison | nothing; same |

That "consumed by: nothing" is a measured statement, not an oversight. The Snakefile
contains no `conda:` directive on any of its nine rules, no `--use-conda` appears in
`cli.py` or in the test harness, and no CI step reads either file. Wiring them into the
rules is a central change to `workflow/Snakefile`, tracked as an open item rather than
claimed here.

```bash
conda env create -f workflow/envs/python.yaml
conda activate batchdx-cyto
pip install -e '.[workflow]'      # the package itself is not on PyPI

# only needed to re-run the upstream-parity comparisons:
conda env create -f workflow/envs/openproblems.yaml
pip install -e '.[metrics,dev]'
```

**There is no per-rule environment file in this build, and there is not one coming from
this directory.** An earlier version of this README claimed that "fine-grained
environment files live here so a rule pulls only what it needs"; nothing did. The claim
was wrong, and it is stated here in the negative on purpose so that it is not repeated.

Making it per-rule would mean adding `conda:` directives to the rules, which is a
central change to `workflow/Snakefile` rather than a file added beside it. It would also
buy less than it sounds: the rules share one process per candidate, the shared modules
are the analytical core, and per-rule environments would have to be solved and rebuilt
for every distinct (rule, platform) pair.

**The R side is not pinnable from a conda environment file at all.** cyCombine and
CytoNorm are distributed through GitHub rather than CRAN or Bioconductor, so their
transitive dependencies are not resolvable by a solver. That is what
`scripts/fetch_upstream.sh` and the container image in `containers/` are for — see
`docs/decisions/0004-dependency-isolation.md`. Rules that shell out to R ask for the
interpreter by name (`upstream.rscript`, `upstream.checkout_dir`) and record
`not_available` with a reason when it is absent, rather than failing the run.
