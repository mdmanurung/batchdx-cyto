# B2 missing-methods scan (Biomni task tsk_014hmFy5M7qWwzX66pNllWHH)

Status: completed 2026-10-09. Full report `B2_missing_methods.md` (19 KB) lives in the Biomni project
(https://biomni.phylo.bio/projects/prj_014ByJfTRi8HUzQHkF0FqKRK/tasks/tsk_014hmFy5M7qWwzX66pNllWHH). This file keeps Biomni's
final summary only. UNTRUSTED until adjudicated.

## Biomni's summary
Top must-haves "missing from the matrix":
- CytoVI: purpose-built cytometry VAE, Nat Methods 2026, fit-once/apply-to-new-batch (claimed as the real answer to the matrix's "scVI" target).
- diffcyt-style DA/DS spike-in evaluation on HDCytoData: ground-truth over-correction test.
- PVCA + variancePartition: graded confounding, replacing the binary flag.
- PCR comparison.
- Public benchmark-dataset bundle: HDCytoData, cyCombine DFCI, FR-FCM-Z4KT / Z3WR.
- Spillover spreading matrix (SSM) for flow/spectral instrument tracking.
Should-haves: C2ST / Batch Probing Score, CytoMDS, NMI/ARI + isolated-label F1/ASW, replicate CV + pseudobulk concordance,
flowClean, cytofQC, CytOpT (optimal transport), MMD-ResNet.
Skip: SAUCIE, DeepCyTOF, optimalFlow, Scanorama/BBKNN/Conos/scMerge, Symphony/scArches/scPoli, scIB cell-cycle/HVG,
flowQ/QUALIFIER (removed from Bioconductor), premessa, CytoPipeline, energy distance; "watch" tier of 2025-26 preprints
(CytoOne, dioscRi, UVAE, RBET).
Target-column flags: limma removeBatchEffect is a mis-target (no cytometry validation); OpenProblems cyto task real but
in development; embedding integration conflicts with the same-shape-matrix contract. Also claimed: LISI attribution
confirmed; spillR (Bioinformatics 2024, LGPL-3) is distinct from Chevrier 2018; Front Immunol 2021 spectral workflow
DOI 10.3389/fimmu.2021.768113.
Ecosystem gaps (no tool exists): CyTOF carryover metric, open CS&T/PMT-voltage tracking, tandem-dye lot QC, and no
spectral-specific batch-correction algorithm or benchmark in 2018-2026.

## Conflicts with B1 / main-session checks (to resolve in adjudication)
- fdaNorm: B2 says "in flowStats, maintained"; B1 says removed; main session found 0 mentions in the current flowStats
  manual -> B2 WRONG on this point.
- OpenProblems counts: B2 repeats the reviewers' "21 methods / 8 metrics"; B1 measured 20 (+3 controls) / 6 via GitHub API.
  Not yet re-checked by the main session.
- Several B2 items (CytoVI Nat Methods 2026, CytoMDS, Batch Probing Score, HDCytoData accessions) are not yet verified
  by the main session; mark roadmap uses of them as "to verify before committing".
