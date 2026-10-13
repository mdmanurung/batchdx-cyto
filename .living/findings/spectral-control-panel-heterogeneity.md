## F-001: The public spectral controls use six distinct FCS channel layouts
**Status:** preliminary
**Claim:** The 72 FCS files in Zenodo record 15388817 contain 27,022,765 events and form a complete four-control grid across 18 batches, but their headers define six distinct channel layouts. Batches 13–18 form a 24-file subset with one exact 33-channel layout.
**Implications:** A common-channel workflow cannot pool all 18 batches without an explicit marker mapping. The six-batch subset can be evaluated for technical-control endpoints; donor and biological outcomes remain unidentified.
**Tags:** cytometry, spectral-flow, panels, validation

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-10-03 | `data/metadata/cytobatchflagr_spectral/audit_fcs.py` | Zenodo 15388817, 72 checksum-verified FCS | batchdx-cyto | Six layouts; 24 files in batches 13–18 have one exact layout and all four controls | supports |

### Open Questions
- Can the other layouts be harmonized by antigen identity while preserving comparable measurement scales?
- Which technical endpoints are identifiable from the six-batch subset without donor IDs?
