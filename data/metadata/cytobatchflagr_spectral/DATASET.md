# CytoBatchFlagR spectral controls: access and design audit

**State (2026-10-03): metadata frozen; raw FCS not downloaded; eligibility not decided.**

## Source and rights

- Public record: [Zenodo 15388817](https://zenodo.org/records/15388817), DOI `10.5281/zenodo.15388817`, titled *CytoBatchFlagR: Spectral flow cytometry dataset*.
- Creators listed by the record: Zachary Koenig, Sandra Andorf, Tamara Tilburgs, and Shruti Eswar. Attribute the record and DOI in derived work.
- The record API declares `cc-by-4.0`. Confirm any third-party rights and the terms attached to the downloaded files before redistribution. This project does not commit raw FCS.
- `source_files.tsv` freezes all 74 record files with sizes, MD5 checksums, and direct URLs. Its 72 FCS filenames match all 72 rows in the record's `Table_2.csv`; `source_Table_2.csv` is an unchanged copy with MD5 `734391c8d456039bf9d4dabefec898bc`.
- `source_Table_1.pdf` is the unchanged antibody table, MD5 `f0ba3fbe9eaaed5b8118102f14f5b02b`. It lists 25 antigens and fluorophores. The correspondence to FCS channel names remains unverified.

## Design known from source metadata

The source description says the FCS files contain FlowJo-gated live `CD45+CD14-` single cells. `Table_2.csv` lists four controls (`tcu1`, `tcu2`, `tcs3`, `tcs4`) in each of 18 batches: 72 files total. The description calls `tcu1` and `tcu2` unstimulated controls and `tcs3` and `tcs4` stimulated controls. The table's `sample_id` is a batch/control label; it does not establish a biological donor ID. Marker names, channel aliases, and event counts remain unverified until the FCS headers and files are read.

**Provisional role:** candidate multi-batch technical-control dataset if FR-FCM-Z247 cannot be obtained. Its repeated control labels support a batch-control audit. It cannot yet support donor-disjoint or biological-preservation claims. No sample or marker exclusions have been made.

## Download and validation recipe

Run from the repository root, supplying a directory outside Git with at least 4 GB free:

```bash
python data/metadata/cytobatchflagr_spectral/download_all.py \
  --dest /path/outside/git/cytobatchflagr_spectral
```

The script downloads every manifest row, verifies bytes and MD5, and atomically moves each verified file into place. Re-running skips files with valid hashes and refuses an existing file with the wrong hash. After download, record the actual local checksums, FCS header marker mapping, event counts, and design-audit output here before declaring eligibility.

## Access check

On 2026-10-03, the Zenodo record API, `Table_1.pdf`, and `Table_2.csv` were reachable. FlowRepository's HTTPS certificate failed validation from this host, so FR-FCM-Z247 could not yet be assessed here. ImmPort's SDY2015 study page returned HTTP 200, but file access and terms were not yet audited.
