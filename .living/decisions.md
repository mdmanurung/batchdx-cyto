# Decision Log

Append-only log of non-obvious decisions and their rationale.

**Entry template:** copy from `skills/core/templates/decision-log-entry.md` (includes Context, Decision, Alternatives considered, Rationale, Consequences, Tags fields).

### 2026-10-03 — Add Mycelium without moving application files

- **Context:** This is an existing Python/R package with imports, packaged workflow assets, tests, and substantial uncommitted work.
- **Decision:** Initialize the Mycelium living layer and manifests in place. Keep `src/`, `workflow/`, `docs/`, `config/`, and existing project files at their current paths.
- **Rationale:** Moving those files would break package paths and workflow references without improving the requested initialization.
- **Consequences:** The Mycelium manifests begin as registries for future analyses and data; the existing package layout remains authoritative.
- **Tags:** mycelium, repository-layout

### 2026-10-03 — Use Python 3.12 as the Alpha CI status environment

- **Context:** The workflow environment pins Python 3.12, while the prior CI job selected 3.11 and the generated status came from a local 3.12 run.
- **Decision:** Select Python 3.12 in the primary CI job and run the R exchange suite there with the job's Python interpreter.
- **Rationale:** The status comparison needs an executed environment matching the local report, and cross-language tests must not silently skip.
- **Consequences:** Python 3.10, 3.11 and 3.13 remain declared-only until an executed matrix establishes support; the CI job remains unverified until it runs.
- **Tags:** ci, release, python
