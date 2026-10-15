# Handoff — whole-package review → Biomni review-of-review → completion roadmap

This is the historical 2026-10-09 review snapshot. For the later authorized
commit/push work and current validation scope, see the
[2026-10-14 journal](.living/log/2026-10-14-001-batchdx-cyto.md).

**Date:** 2026-10-09 · **Branch:** alpha-ci-20261003 @ b313006 (local, unpushed) · **Status:** Phase 1 done; Phase 2 (Biomni) 2/6 done, 3 running, 1 queued; Phase 3 not started

## Goal
Review all of batchdx-cyto, have Biomni check the review, then write `review/2026-10-09/ROADMAP.md` to a COMPLETE
cytometry batch-correction package (mass/flow/spectral; full diagnostics + correction). Analysis only: no source edits,
no commits, no pushes. Write only inside `review/2026-10-09/`.

## Next action
Poll Biomni project `prj_014ByJfTRi8HUzQHkF0FqKRK` (`list_tasks`, then `wait_for_next_update` on each running task) and
start B6 (message is in the original prompt: adversarial review of the review; files below) as soon as a slot frees.

## State
- Uncommitted: all of `review/` (new) plus pre-existing doc/data/gate edits (not mine). Nothing committed.
- Phase 1 files: `review/2026-10-09/{00_README,01_inventory,02_capability_matrix}.md`, `lenses/L1..L10_*.md` (159 findings: 2 blocker, 55 major, 79 minor, 23 nit).
- Biomni (files uploaded: `findings_short.md` = file_0131NoX3oRnRf3b0GIvcsGnn, `capability_matrix_condensed.md` = file_011cEbScD4PQFHH41x2Cxen6):
  - B1 literature `tsk_011yPcBj12PsEbaBcq8UZwsY` DONE; B2 missing methods `tsk_014hmFy5M7qWwzX66pNllWHH` DONE. Summaries in `review/2026-10-09/biomni/B1_*.md`, `B2_*.md`.
  - B3 metrics `tsk_012ZxPdWEFfpmHNFRpMOdTVM`, B4 public data `tsk_014M1q1iS6pl0jTjis5Wj5yX`, B5 stats `tsk_011DcDdn8jhKq3ZqEyZamv9Z`: RUNNING (~35-55 min each, started ~14:50).
  - B6 adversarial: NOT STARTED (429 while 3 tasks run).
- Claude lens agents are all finished.

## Locked decisions
- Biomni pass approved: run all 6 tasks, docs-only upload, paths/user/email stripped (user, 2026-10-09).
- I drive Biomni tasks myself, not via subagents (cheaper on the Claude limit); uploads are condensed rewrites, not the full files (upload tool is inline-text only).
- Biomni output is untrusted; check each claim against repo or primary source before it enters `03_adjudication.md` (ACCEPTED / REJECTED / REVISED).

## Dead ends — do not redo
- Launching 11 Claude subagents at once hit the usage limit twice and lost work. Launch ≤3-5; resume with SendMessage (agents keep context).
- Biomni `start_new_task` returns 429 when 3 tasks are running; wait, don't retry in a loop.
- Biomni result files cannot be downloaded via MCP; only the final text comes back. Full reports are on the web UI.

## Read first
- `review/2026-10-09/02_capability_matrix.md` — gaps ranked; top-10 list at the end.
- `review/2026-10-09/biomni/B2_missing_methods.md` — conflicts with B1 to resolve (fdaNorm, OpenProblems counts).
- `lenses/L4_evaluation_metrics.md` L4-01/02, `L1` L1-01/02, `L7` L7-01 — metrics not wired into the pipeline (3 reviewers agree).
- `src/batchdx/cyto/metrics/openproblems.py:363-364` — verified: EMD rows get baseline `panel`, not `candidate`.
- `.mycelium/last-session.md`, `.living/decisions.md` D-5 — earlier session context (SVG fix b313006).

## Open items
- Adjudication conflicts: fdaNorm (main session: 0 mentions in current flowStats manual, so B1 right, B2 wrong); OpenProblems 20 methods/6 metrics (B1) vs 21/8 (B2, reviewers) — re-check via GitHub API; CytoNorm 2.0 is a separate control-free paper (B1) — matters for L2-01; B2's CytoVI/CytoMDS/BPS/HDCytoData claims unverified.
- Refresh the capability matrix against L7-L10 (it was built from L1-L6 only).
- Unverified Claude-lens claims worth a quick check before the roadmap leans on them: L10-01 (Jinja autoescape off; `html.py:533` consistent), L3-02 (Kruskal 2v2 can't reach p<0.05; B3 will test), L6-05, L7-11.
- Phase 3: write `03_adjudication.md`, then `ROADMAP.md` (definition of complete, ranked gaps, milestones M0..Mn with days/gates, critical path, ≤5 decisions for the user, first 3 actions <2h each); trace each item to an accepted finding.
- Mycelium post-action protocol (learnings/decisions) after the roadmap; push of b313006 still needs the user's explicit OK.
