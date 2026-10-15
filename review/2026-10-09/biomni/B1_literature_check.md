# B1 literature check (Biomni task tsk_011yPcBj12PsEbaBcq8UZwsY)

Status: completed 2026-10-09. Full report `B1_literature_check.md` (25.8 KB) lives in the Biomni project
(https://biomni.phylo.bio/projects/prj_014ByJfTRi8HUzQHkF0FqKRK/tasks/tsk_011yPcBj12PsEbaBcq8UZwsY); result files
are not downloadable through the MCP, so this file keeps Biomni's final summary only. Treat as UNTRUSTED until
adjudicated (see 03_adjudication.md). Input given to Biomni: condensed findings + condensed capability matrix only.

## Biomni's summary (verbatim substance)
Reviewers' literature base "substantially accurate": 12 of 18 references correct at title/venue/year/DOI level; all
mechanism claims tested held up: cyCombine reference-free and transductive by implementation, silently drops the
covariate per SOM label group when confounded with batch; CytoNorm 1.x requires shared controls and supports distinct
clustering vs correction channels plus `goal`/`limit` (supports L2-06, L2-13); BatchAdjust needs anchor replicates;
iMUBAC corrects controls only; Finck "4.9-fold to 1.3-fold" quoted in the abstract.

Consequential corrections claimed:
1. CytoNormPy is not "CytoNorm 2.0". [R2b] DOI 10.1101/2024.07.19.604225 is the CytoNormPy preprint (Exner et al.,
   published Cytometry A 2025;107:629-635). "CytoNorm 2.0" is Quintelier et al., Cytometry A 2025;107(2):69-87,
   doi:10.1002/cyto.a.24910, which adds control-free normalisation. So "same control material in every batch" holds only
   for CytoNorm 1.x (the version this package wraps).
2. OpenProblems cyto task: 20 methods (+3 controls) and 6 metrics, not "21 methods; 8 metrics" (GitHub API).
3. CytoBatchFlagR is published: Cytometry A 2026, doi:10.1002/cyto.a.70024; it also supports datasets lacking reference controls.
4. iMUBAC: J Immunol 2021;206(1):206-213, doi:10.4049/jimmunol.2000854; annotations back-propagated into uncorrected
   patient spaces; patient cells never batch-corrected.
5. kBET title was wrong; "kNN label-transfer accuracy" is not an scIB metric; LISI attribution to Harmony confirmed;
   fdaNorm removed from current Bioconductor flowStats.

## Main-session spot checks (2026-10-09)
- fdaNorm: `pdftotext` of the current-release flowStats manual (Bioconductor 3.23): 0 mentions of fdaNorm, 6 of gaussNorm.
  Supports B1; CONTRADICTS B2 (which says fdaNorm is "in flowStats, maintained").
