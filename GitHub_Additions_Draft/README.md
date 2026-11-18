# TBI LCP1/CORO1A study materials — local additions draft

Repository: https://github.com/younng367/analysis-code

This directory adds study-specific files to the existing 17 analysis templates. It is staged locally and has not been pushed, released or assigned a DOI. The submitted minor-revision manuscript is the reference for scientific content. The current animal data include the latest author-supplied LCP1 IHC values; older animal-results archives are superseded.

## Contents

- `data/network`: all 16 candidates, original cytoHubba exports, full ranks, export order and saved confidence sensitivity results. Absent original scores are NA; recomputed isolated-node scores belong to the separate sensitivity analysis.
- `data/markers`: cell-type markers and the mechanics, immune-signature and ECM panels used in analysis. Check third-party source terms before redistribution.
- `data/animals`: current analysis values, IHC field records and WB ROI records. `source_row` is an assay-specific row index, not a newly inferred cross-assay animal identifier.
- `data/metadata` and `data/trajectory`: library metadata and saved coordinates/graph inputs.
- `results`: saved bulk, pseudobulk, immune-correlation, GSEA, ECM, trajectory and animal statistics. Directory names describe analyses; old Table numbers in retained filenames are provenance labels, not the current journal supplement numbering.
- `archived_analysis`: actual analysis scripts. These retain their original directory conventions and are not yet a portable end-to-end workflow. Earlier diagnostic calculations within these sources must not be substituted for the manuscript's declared primary selection rule.
- `environment`: recorded software versions; workstation username prefixes were redacted.
- `study_parameters.json`: parameters transcribed from the manuscript and scripts.

## Numerical checks that can be run from this directory

From a terminal, with Python 3 and R installed:

```sh
python3 scripts/check_network_top3.py .
Rscript scripts/check_animal_statistics.R . ./check_output
```

The first check verifies all 30 archived scores and the fixed-three intersection. The second reproduces the current six animal comparisons and checks every reported numeric column against the supplied summary at tolerance 1e-10. These checks do not reconstruct sequencing preprocessing or image acquisition.

## Public source data

- Bulk cortex: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79441
- Single-nucleus: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277487
- Spatial: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236171

Raw sequencing downloads need not be duplicated here. A complete public reproduction package still needs the exact accession-to-file/sample mapping, preprocessing/annotation instructions, and either the processed inputs below or runnable instructions to recreate them:

- Gene-by-sample integer counts used by the bulk script (`BULK/bulk/GSE79441_gene_exp_filtered.csv`).
- The annotated microglial object `Micro_ECM.h5ad` with raw counts, library IDs, states and score annotations; the marker object `Microglia_DEG.Rdata`.
- The archived Monocle object `monocle/micro_monocle.qs`, or a documented reconstruction from public counts and annotations.
- The selected spatial object `GSE236171/old/612_crop_st.h5ad`, including the original region definitions and coordinates, or a documented reconstruction.
- Original CellChat and Milo analysis inputs/scripts/settings needed for their reported panels.

Some derived matrices and graph inputs are already included for numerical inspection. They are not a substitute for documenting the preceding processing steps. WB/IHC source images and measurement records are supplied separately with the journal submission; no reviewer report, response letter or ethics document belongs in this repository.

## Release preparation

Finish the input mapping and portable run instructions; validate outputs; confirm redistribution permissions and choose suitable licenses. Then create a real tagged GitHub Release and a persistent versioned archive (for example a Zenodo release DOI). No release number, DOI or complete-reproduction claim has been invented in this draft.
