# Documentation index

| File | What it is |
|---|---|
| [`implementation_status.md`](implementation_status.md) | **generated** from an executed test report; the authority on what passed |
| [`usage.md`](usage.md) | installation, input preparation, all three modes, reference requirements, confounding, metric interpretation, calibration, failure recovery, memory, privacy, extension interfaces |
| [`upstream_audit.md`](upstream_audit.md) | what the pinned upstream sources actually say, including corrections to the plan's premises |
| [`validation_results.md`](validation_results.md) | what the scientific suite establishes and what it does **not**, including the synthetic-only scope limit and the exact scope of the multi-seed sweep |
| [`blocked.md`](blocked.md) | every blocked item with the exact command, the error, and the affected stage |
| [`compatibility_matrix.md`](compatibility_matrix.md) | what is actually tested, declared only, or untested, with the command that verifies each row — Python versions, dependency versions and bounds, R and the R backends, and the statement that no real cytometry data ships |
| [`decisions/`](decisions/) | architecture and licence decision records |
| [`plans/`](plans/) | implementation plans and progress records, tracked deliberately: the remaining work, the tranches executing it, and what each tranche verified |
| [`releases/0.1.0-alpha.md`](releases/0.1.0-alpha.md) | draft Alpha release text, held until the checklist and final approval |
| `audit/` | captured Stage 0 audit artifacts (produced by `batchdx-cyto audit`) |

Two files at the repository root belong in this index's reading order even though they
are not under `docs/`:

| File | What it is |
|---|---|
| [`../CHANGELOG.md`](../CHANGELOG.md) | every notable change, with known issues stated rather than omitted |
| [`../RELEASE.md`](../RELEASE.md) | the release checklist: what is verified, what is declared-only, and what must change before the next release is honest |

Read `implementation_status.md` and `blocked.md` before relying on any capability, and
`compatibility_matrix.md` before relying on any version claim.

**On `plans/`.** The roadmap and the tranche records are tracked rather than kept local,
for one reason: several of their entries record *corrections* to earlier claims that
turned out to be false. A plan whose premise was wrong is more useful tracked than
deleted, because the next reader needs to know that the evidence was re-checked rather
than assumed. Their contents are plans and findings, not authority — the authority is
`implementation_status.md` (generated), `compatibility_matrix.md` (tested vs declared)
and `blocked.md` (what still fails).
