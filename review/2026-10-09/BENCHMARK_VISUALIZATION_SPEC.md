# Benchmark ranking visualization for the first prototype

Date: 2026-10-10. Status: required by user; implementation pending.
This specifies presentation and descriptive ranking of the nine required
OpenProblems cytometry scores. It does not establish calibrated recommendations.

## Reference and intended result

The user requested the same benchmark ranking visualization style as
OpenProblems/scIB. The primary references checked on 2026-10-10 are:

- [OpenProblems cytometry leaderboard](https://www.openproblems.bio/benchmarks/cyto_batch_integration/):
  control-anchored scores, ordering by scaled overall mean, and per-metric
  normalization plots with raw-score axes and control anchors.
- [OpenProblems overview](https://www.openproblems.bio/): method-by-metric/dataset
  funkyheatmap and scaling relative to control methods.
- [scIB plotting implementation](https://github.com/yoseflab/scib-metrics/blob/main/src/scib_metrics/benchmark/_core.py):
  grouped metric columns, coloured numeric circles, summary bars, and rows sorted
  by total. The installed implementation was also inspected. Its default total
  weights batch correction by 0.4 and biological conservation by 0.6; those are
  scIB defaults, not the cytometry leaderboard's aggregation rule.

Keep metric calculations at the existing pinned OpenProblems commit
`37151479aaa89edf37ae23068de8fbb41a1b7bb4`. The current website is a visual
reference, not authorization to replace that metric pin. Before implementation,
pin and inspect the upstream reporting/scaling implementation separately: the
local task delegates rendering to an unpopulated `common` submodule. Website
documentation does not prove exact aggregation or missing-value parity.

## Required outputs

1. **Ranked benchmark table.** Rows are method/configuration candidates, including
   identity; columns are all nine named scores with grouped headers. Show numeric
   scores inside coloured circles, batch and biology summary bars, overall bar,
   numeric rank, and coverage/status. Keep failed and unsupported candidates in
   the table below the completely evaluated rows. The static layout follows scIB;
   the report permits sorting and inspecting raw values without external assets.
2. **Normalization QC.** One compact panel per metric: candidate and benchmark
   control points on a common higher-is-better scaled axis, labelled 0/1 anchors,
   raw values and original direction, and explicit out-of-range indicators.
   Controls used for scaling are artificial benchmark control methods, distinct
   from physical control samples that upstream evaluation may exclude.
3. **Inspectable exports.** The offline HTML report, standalone SVG/PDF/PNG,
   and CSV/JSON come from the same canonical results. Raw values, normalized
   values, ranks, coverage, anchors, weights, comparison scope, and statuses are
   accessible without relying on colour or hover. Rerendering runs no corrections.

The optional multi-dataset view uses the same method ordering plus dataset
summary columns once multiple full-suite datasets exist. A single real dataset
is sufficient for the scoped prototype; no additional correction backend is
required solely to populate the figure.

## Column mapping

Use these presentation groups. This grouping is a package design choice grounded
in metric descriptions, not an assertion that upstream declares scIB categories.
FlowSOM assesses consistency of cell-type mapping across technical replicates.

| Group | Score IDs | Original direction |
|---|---|---|
| Batch / replicate consistency | `emd_mean_ct_horiz`, `emd_mean_ct_vert`, `average_batch_r2_ct` | lower |
| Batch / replicate consistency | `iLisi`, `flowsom_mean_mapping_similarity` | higher |
| Biological preservation | `cLisi`, `ratio_consistent_peaks`, `functional_marker_preservation_wilcoxon` | higher |
| Biological preservation | `functional_marker_preservation_cohens_d` | lower |

Category bars summarize their declared group; the overall bar averages the nine
metric scores with equal metric weight. Display that formula in the caption.
Category bars do not introduce 40/60 weighting. Normalized upstream cLISI already
has higher-is-better direction; do not invert it as though it were raw cLISI.

## Scaling and ranking contract

- Compute benchmark controls on the same eligible input and two-split contract
  as the candidates. Preserve their method IDs, seeds, source pin, raw results,
  and metric-specific anchor selection. Identity alone cannot establish two
  distinct best/worst anchors. Benchmark controls are diagnostic references,
  not selectable corrections; never export a synthetic positive-control output
  as the user's corrected dataset.
- For a verified linear control-anchor mapping, scaled score is
  `(raw - worst_anchor) / (best_anchor - worst_anchor)`. The same formula handles
  lower-is-better metrics because the denominator is negative. Verify the actual
  reference rule, including anchor selection, against the pinned renderer before
  claiming OpenProblems parity; do not substitute candidate-dependent min/max
  scaling. Preserve raw and unclipped scaled values; clamp only drawing with a
  visible overflow indicator. Equal/non-finite/missing anchors give an explicit
  scaling-unavailable result. Do not create an overall rank from unavailable
  scaling.
- Within one dataset, the descriptive overall score is the arithmetic mean of
  all nine valid, control-scaled metric scores. For multiple datasets, use equal
  dataset weight after each dataset's complete nine-score mean. Record the fixed
  comparison dataset set and formula before inspecting method results. This is
  the package's explicit aggregation contract; verify reference multi-dataset
  parity before claiming equivalent upstream overall scores.
- Only candidates complete on that fixed scope receive an overall rank. Do not
  use a missing-value-skipping mean or silently rank methods on different metrics
  or datasets. Partial category summaries disclose coverage; incomplete categories
  have no summary bar. All raw valid cells remain inspectable. A failed score is
  never zero and never removed to improve the candidate's rank.
- Rank descending using full precision and competition ties (1, 1, 3). Display
  ties explicitly; stable method/configuration IDs order tied rows without
  inventing a preference. Per-metric ranks follow declared score directions and
  disclose the evaluable comparison set. Method filtering changes only the
  displayed rows, not anchors or original benchmark ranks.
- Identify each benchmark by inputs, donor/sample coverage, correction masks,
  splits, seeds, metrics pin, implementations, controls, and aggregation version.
  Group independently executed repetitions before ranking by a declared rule;
  do not treat the two required metric splits as independent benchmark repeats.
  No fabricated uncertainty interval is shown for a single execution.
- A leading descriptive rank does not pass protection, applicability, failure,
  or calibrated decision gates automatically. Keep recommendation status visible
  separately from benchmark rank.

## Layout and accessibility

Show best complete candidates first with identity clearly labelled. Use clean
short metric headings, full IDs in the accessible table/details, original
direction arrows, and a legend saying larger scaled values are better. Numerical
labels and bar lengths carry the values independently of colour. Distinguish
failed, unsupported, undefined, unavailable, and missing states with text badges
and neutral hatching; reserve coloured score glyphs for valid numbers. Use a
colour-vision-accessible palette and readable text contrast. Show coverage as
`9/9`, the evaluated dataset scope, and the scaling/aggregation caption beside
the figure. SVG IDs and metadata must follow the existing stable rendering helper.

The attached [layout preview](prototype-evidence/benchmark-layout-preview.png)
uses invented candidates and synthetic already-scaled values to demonstrate the
table. Its values and ranks are not measurements or evidence about any method.
It is reproduced by [the review script](prototype-evidence/render_benchmark_preview.py).

## Work placement and acceptance checks

This extends P08 and P20 in [PROTOTYPE_TASKS.md](PROTOTYPE_TASKS.md). P06–P09
provide scores; P13 provides split-correct corrections; P16–P17 provide the
eligible measured demonstration. P08 also needs the benchmark controls and
scaling provenance. Estimate 1–2 engineer-days for presentation/export work
once canonical scores exist, plus verified control/scaling integration. This
effort is provisional and not an executed timing measurement.

Reuse `src/batchdx/cyto/reporting/plots.py`, `html.py`, and the existing template
and canonical metric contracts. Matplotlib is already a core dependency;
`plottable` is installed through the metrics environment and can be reused if
declared for the supported install. Do not add a plotting framework or make
computing scIB RNA metrics a dependency of the cytometry figure.

Acceptance must cover:

1. Every configured candidate and all nine metric IDs appear in figure and data
   exports, including identity, failures, and a candidate with only eight scores.
2. Lower/higher directions and already-normalized cLISI give correct order;
   fixed anchors stay fixed when an unrelated candidate is added or filtered.
3. Zero-width/missing anchors, non-finite results, out-of-range valid scores,
   ties, and missing required results are represented truthfully; missing evidence
   cannot promote a candidate into the complete overall ranking.
4. Plotted totals/ranks reproduce CSV/JSON from full-precision values and a
   documented fixed comparison scope. Normalization panels expose raw anchors.
5. An actual eligible measured-data run supplies all nine scores and the rank
   visualization. Synthetic layout snapshots alone do not close this gate.
6. Offline sorting/details and accessible value table work; labels are escaped
   and redacted; identical results rerender to stable SVG/HTML without correction
   execution. Exported vectors remain readable at intended figure size.

No package plotting code or scientific benchmarking was implemented by this
specification update. The preview is a review artifact, not a product renderer.
The [visualization review receipt](prototype-evidence/2026-10-10-benchmark-visualization.json)
records successful preview assertions, repeat-byte-stable SVG, scoped lint, local
link checks, and unchanged reviewed source hashes. Mycelium strict structure
validation also passed. These checks establish artifact consistency only.
