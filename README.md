# TTFields Sensitivity Study

Computational modelling of Tumour Treating Fields (Optune) for glioblastoma:
sensitivity to tissue conductivity and tumour location in an anatomically
realistic FEM head model.

**Author:** Abdel Rahman Mousa
**Supervisor:** Prof. William Holderbaum
**Module:** BI3RP3 Final-Year Research Project, BSc Biomedical Engineering, University of Reading
**Year:** 2025–2026

## Project summary

This repository contains the MATLAB scripts, run inputs, and exported analysis
tables produced for a final-year research project on TTFields modelling.

The project uses a hierarchical forward-modelling workflow:
1. Simplified 1D MATLAB models of layered conductors (skin → skull → CSF → GM → WM).
2. A conceptual 2D finite-element MATLAB model with a tumour inclusion.
3. Anatomically realistic 3D ROAST simulations on the MNI152 template
   (FT7 → T8 montage, ±47.12 mA total injection).
4. A paired comparison of two conductivity assignments: ROAST defaults
   versus a Fan-like literature-matched set (Fan et al., 2023, Table 2).
5. Whole-brain metric extraction (median, P95) for GM, WM, GM+WM and CSF.
6. A spherical tumour-region-of-interest sweep across five lateral positions
   and three axial depths (15 ROI conditions per conductivity run).

The contribution is a structured sensitivity study quantifying how a single
reported TTFields dosimetric number can shift under plausible modelling
choices.

## Repository structure

- `01_simplified_models/` — 1D and 2D simplified MATLAB models.
- `02_roast_analysis/` — analysis scripts applied to ROAST output volumes
  (tumour ROI sweep, run-bundling/figure pipeline).
- `03_run_inputs/` — saved conductivity input files for the Fan-like run.
- `04_tables/` — CSV exports underlying all reported numerical results.

The full development record (dated entries, troubleshooting, supervisor
meetings) is held separately in an electronic notebook linked from the
dissertation.

## Environment

- MATLAB R2022a (later versions showed compatibility issues with the
  ROAST workflow during testing).
- ROAST v3.0 (https://www.parralab.org/roast/)
- SPM12
- Gmsh
- GetDP

## Reproduction

The workflow assumes a working ROAST installation with the MNI152 template
(`MNI152_T1_1mm.nii`).

1. Run a baseline simulation with ROAST defaults using the FT7 → T8 montage
   at ±47.12 mA total injection.
2. Run a second simulation with the Fan-like conductivity values listed in
   `03_run_inputs/run_inputs_20260211T182002_summary.txt`, keeping all other
   parameters identical.
3. Apply the analysis scripts in `02_roast_analysis/` to the resulting
   `*_emag.nii` volumes and the SPM12 segmentation mask
   (`MNI152_T1_1mm_ras_T1orT2_SPM_masks.nii`) to reproduce the tables in
   `04_tables/`.

The tumour ROI is defined as the intersection of a sphere (radius 20 voxels
= 20 mm) with the GM/WM mask:

```matlab
[X,Y,Z]   = ndgrid(1:size(M,1), 1:size(M,2), 1:size(M,3));
sphere    = (X-cx).^2 + (Y-cy).^2 + (Z-cz).^2 <= r^2;
tumourROI = sphere & ((M==1) | (M==2));
```

ROI centres used in the extended sweep: Base (93,104,z), X+ (118,104,z),
X− (68,104,z), Y+ (93,129,z), Y− (93,79,z), at z = 68, 83, 98.

## Reference

Fan, Y. et al. (2023). Electric field simulation for Tumor Treating Fields
using FEM. *Proceedings of the 2023 8th International Conference on
Biomedical Signal and Image Processing (ICBIP '23).* doi:10.1145/3613307.3613328

## Licence

MIT — see `LICENSE`.
