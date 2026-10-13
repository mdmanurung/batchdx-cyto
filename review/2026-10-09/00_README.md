# Review of batchdx-cyto — 2026-10-09

**Goal:** assess the package against a *complete* cytometry batch-correction package
(mass, flow, spectral; full diagnostics; full correction capacity) and produce a roadmap.

**Repo state reviewed:** branch `alpha-ci-20261003`, HEAD `b313006` (local, unpushed),
with uncommitted doc/data/gate edits present in the worktree. Status record =
`RELEASE.md`, `scripts/gates.yaml`, executed reports. Older roadmaps are historical.

**Constraints:** analysis only; no source edits/commits/pushes; every finding cites
`path:line` or an executed artifact.

## Method
1. Phase 0 orient → `01_inventory.md`.
2. Phase 1: one reviewer per lens (L1–L10) → `lenses/`, plus `02_capability_matrix.md`.
3. Phase 2: Biomni review-of-review (B1–B6) → `biomni/`, `03_adjudication.md`. **Gated on user OK.**
4. Phase 3: `ROADMAP.md`.

## Finding ID scheme
`L<lens>-NN`, schema: `ID | lens | severity | claim vs reality | evidence | fix | effort (days) | confidence`.
Severity: blocker / major / minor / nit.

## File index
- `01_inventory.md` — modules, LOC, tests, gate map
- `lenses/L1..L10_*.md` — lens reviews
- `02_capability_matrix.md` — target vs implemented vs real-data vs gap
- `biomni/`, `03_adjudication.md`, `ROADMAP.md` — later phases
