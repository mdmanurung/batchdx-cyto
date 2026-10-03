# Shipped configurations

Three complete, runnable configurations, one per design situation the package has to
handle honestly:

| File | Design | What it exercises |
|---|---|---|
| `example_paired_reference.yaml` | three batches, a measured reference aliquot in each, balanced case/control biology | the reference-supported path end to end |
| `example_no_reference.yaml` | three batches, biologically overlapping, no measured technical replicate | conditional evidence; the design floor with `require_shared_reference: false` |
| `example_confounded_diagnostic_only.yaml` | two batches, perfectly confounded with biological group | abstention: non-identifiable contrasts and a below-floor design |

Each one is materialised together with its data and manifests by:

```
batchdx-cyto init --example paired_reference my_project
batchdx-cyto init --example no_reference     my_project
batchdx-cyto init --example confounded       my_project
```

`init` writes a copy of the corresponding file above into `my_project/config/config.yaml`
— the same configuration key for key, with the same values, but a one-line comment header
instead of the multi-line one in `config/` (everything below the header is byte for byte
the same). It also writes `data/preprocessed.h5ad`, `config/samples.tsv`,
`config/markers.tsv`, `config/processing.json` and `config/thresholds.tsv`. Every path
in the configuration is relative to the project directory, so nothing depends on a
developer's machine.

`upstream_lock.yaml` is not a run configuration: it is the Stage 0 pinned manifest of
upstream repositories, commits, licences and the exact functions the adapters call.
