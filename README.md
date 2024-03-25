# CellPro: Cell Proximity Analysis between Cell Types Using IMC Data of LUAD Patients

## Introduction
This repository contains the analysis of cell intensity, fraction, immune fraction, cell-cell proximity, and corresponding R codes computed using LUAD IMC data extracted from Mark Sorin et al. (PMID: 36725934).

### Contents:
1. **Main.RData**: Includes the following variables:
   - `IMC_clinical`: Clinical data and survival records of 416 LUAD patients.
   - `CellPositionList`: Coordinates of each cell in every IMC image.
   - `IMCsize`: Pixel dimensions on the X and Y axes of each IMC image.
   - `IMC_proximity`: Proximity scores of 169 cell pairs in each sample.
   - `IMC_intensity`: Intensity (or density) of 16 cell types in each sample.
   - `IMC_fraction`: Cell fractions of 16 cell types in each sample.
   - `IMC_fraction_immune`: Immune cell fractions of 14 immune cell types in each sample.
   - `IMC_proximityPermut100Z`: Normalized proximity scores of 169 cell pairs calculated using 100 permutations in each sample.
   - `M2ProximalCellCountList`: Cell counts of M2-proximal and M2-distal cells in each sample.
   - `CoxUniProximity169`: Results of Univariable Cox regression using proximity of 169 cell pairs from 416 patients, including P-values and Hazard ratios (HR).

2. **IMC_CellFraction_Intensity_calculation.R**: R code for calculating intensity and fraction of each cell type in every IMC sample, either in all cells or immune cells.

3. **IMC_CellCellProximity_calculation.R**: R code for calculating X⮕Y proximity between each cell type in every IMC sample and NULL proximities using 100 permutations. Results of each permutation may slightly vary.

4. **M2proximal_cellcount_calculation.R**: R code for categorizing cell types as M2-proximal or M2-distal based on their closest distances to M2 Macrophages using cutoffs from 10 pixels to 100 pixels, and for calculating their cell counts.

5. **M2proximal_distal_SurvivAlanalysis.R**: R code for calculating the cell fraction of M2-distal and M2-proximal cell types using previously calculated cell counts and performing survival analysis.

6. **Hierarchical_clustering_of_10CCI.R**: R code for hierarchical clustering of patients based on 10 X⮕Y proximities significantly associated with overall survival.
