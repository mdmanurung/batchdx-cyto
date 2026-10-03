# CytoBatchFlagR spectral controls: access and design audit

**State (2026-10-03): all source files downloaded and checksum verified; a 24-file homogeneous technical-control subset is eligible for workflow review. No correction result or biological claim has been validated.**

## Source and rights

- Public record: [Zenodo 15388817](https://zenodo.org/records/15388817), DOI `10.5281/zenodo.15388817`, titled *CytoBatchFlagR: Spectral flow cytometry dataset*.
- Creators listed by the record: Zachary Koenig, Sandra Andorf, Tamara Tilburgs, and Shruti Eswar. Attribute the record and DOI in derived work.
- The record API declares `cc-by-4.0`. Confirm any third-party rights and the terms attached to the downloaded files before redistribution. This project does not commit raw FCS.
- `source_files.tsv` freezes all 74 record files with sizes, MD5 checksums, and direct URLs. Its 72 FCS filenames match all 72 rows in the record's `Table_2.csv`; `source_Table_2.csv` is an unchanged copy with MD5 `734391c8d456039bf9d4dabefec898bc`.
- `source_Table_1.pdf` is the unchanged antibody table, MD5 `f0ba3fbe9eaaed5b8118102f14f5b02b`. It lists 25 antigens and fluorophores. `channel_layouts.tsv` records measured FCS channel names and labels; a final assay-marker manifest still needs review.

## Design known from source metadata

The source description says the FCS files contain FlowJo-gated live `CD45+CD14-` single cells. `Table_2.csv` lists four controls (`tcu1`, `tcu2`, `tcs3`, `tcs4`) in each of 18 batches: 72 files total. The description calls `tcu1` and `tcu2` unstimulated controls and `tcs3` and `tcs4` stimulated controls. The table's `sample_id` is a batch/control label; it does not establish a biological donor ID.

The local FCS audit read and checksum verified all 72 FCS files, totaling **27,022,765 events**. There are **six distinct channel layouts**, and each batch uses one layout. The largest exact-layout group is batches **13–18**: 24 files, one of each control per batch, all with the same 33 channels. `technical_control_subset.tsv` freezes that subset. The other **48 files are excluded from this first technical-control analysis solely because their channel layouts differ**; they remain in the full source inventory. The full 18-batch set is ineligible as-is for a common-marker pipeline. Cross-panel harmonization would need an explicit marker mapping and separate validation.

**Role decision:** the 24-file, six-batch subset is eligible for a technical-control workflow review, subject to a marker manifest and the raw-FCS preprocessing implementation. It cannot support donor-disjoint or biological-preservation claims because donor and outcome identities are not established. No marker exclusions have been made.

## Download and validation recipe

Run from the repository root, supplying a directory outside Git with at least 4 GB free:

```bash
python data/metadata/cytobatchflagr_spectral/download_all.py \
  --dest /path/outside/git/cytobatchflagr_spectral
```

The script downloads every manifest row, verifies bytes and MD5, and atomically moves each verified file into place. Re-running skips files with valid hashes and refuses an existing file with the wrong hash. Audit the downloaded FCS files with:

```bash
python data/metadata/cytobatchflagr_spectral/audit_fcs.py \
  --data-dir /path/outside/git/cytobatchflagr_spectral \
  --out-dir data/metadata/cytobatchflagr_spectral
```

`fcs_inventory.tsv` freezes each FCS SHA-256, event count and channel count. `file_layouts.tsv`, `channel_layouts.tsv`, `technical_control_subset.tsv`, and `design_audit.json` record the panel split and eligibility decision. The audit fails if any downloaded FCS differs from the source manifest or if no complete homogeneous multi-batch control subset exists. Raw FCS remain outside Git.

## Access check

On 2026-10-03, the Zenodo record API, `Table_1.pdf`, and `Table_2.csv` were reachable. FlowRepository's HTTPS certificate failed validation from this host, so FR-FCM-Z247 could not yet be assessed here. ImmPort's SDY2015 study page returned HTTP 200, but file access and terms were not yet audited.
