# Examples

Every example is materialised by the CLI, data included, so nothing here depends on a
path that exists only on a developer's machine:

```bash
batchdx-cyto init my_project --example paired_reference
batchdx-cyto init my_project --example no_reference
batchdx-cyto init my_project --example confounded
```

Each writes `data/preprocessed.h5ad`, `config/samples.tsv`, `config/markers.tsv`,
`config/processing.json`, `config/thresholds.tsv`, `config/config.yaml` and a `README.md`
with the exact commands to run. The configuration written is identical to the
corresponding file in `config/` (see `config/README.md`).

| Example | Design | The point of it |
|---|---|---|
| `paired_reference` | 3 batches, a measured reference aliquot in each, balanced case/control, technical offset on CD4 and HLA-DR | the reference-supported path end to end |
| `no_reference` | 3 batches, overlapping biology, no measured technical replicate | evidence is conditional, not reference-supported |
| `confounded` | 2 batches, perfectly confounded with group, no reference | the system must abstain |

The data are synthetic fixtures from `batchdx.cyto.simulate`. They exercise the
workflow; they are not a scientific claim about real cytometry, and no real-data example
is included (see `docs/blocked.md` B-008).
