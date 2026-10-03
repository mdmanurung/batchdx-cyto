# Decision Log

Append-only log of non-obvious decisions and their rationale.

**Entry template:** copy from `skills/core/templates/decision-log-entry.md` (includes Context, Decision, Alternatives considered, Rationale, Consequences, Tags fields).

### 2026-10-03 — Add Mycelium without moving application files

- **Context:** This is an existing Python/R package with imports, packaged workflow assets, tests, and substantial uncommitted work.
- **Decision:** Initialize the Mycelium living layer and manifests in place. Keep `src/`, `workflow/`, `docs/`, `config/`, and existing project files at their current paths.
- **Rationale:** Moving those files would break package paths and workflow references without improving the requested initialization.
- **Consequences:** The Mycelium manifests begin as registries for future analyses and data; the existing package layout remains authoritative.
- **Tags:** mycelium, repository-layout
