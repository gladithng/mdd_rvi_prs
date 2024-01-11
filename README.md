# Brain-based and genetic risk scores for depression
<p>Scripts used for the following paper: "Comparing personalized brain-based and genetic risk scores for major depressive disorder in large population samples of adults and adolescents"</p>  
  
Published in [European Psychiatry](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9393914/)

## prep directory
* prep_abcd_disease_antidep.R: determine subjects' neurological, psychiatric and antidepressant-use status (in ABCD)
* prep_abcd_mdd_covariates.R: extract depression phenotype and compile necessary covariates (in ABCD)
* prep_abcd_unrel_whiteonly.R: extract IDs for neurologically healthy, unrelated (determined through family_ID) and white-only (determined through k-means clustering of genetic principal components) subjects
* prep_abcd_imaging_RVI.R: tidy up (QC) neuroimaging data for each neuroimaging measure and reformat to match ENIGMA tables (e.g., matching tracts because different atlases were used) (in ABCD)
* prep_abcd_enigma_forRVI.R: reformat ENIGMA tables to suit ABCD naming conventions (in ABCD)
* prep_GS_imaging_enigma_forRVI.R: tidy up (QC) GS imaging data and reformat ENIGMA tables (in GS)
* prep_GS_unrel_disease.R: extract IDs for unrelated (determined through family ID) subjects and to determine disease status (neurological and psychiatric) (in GS)
* prep_GS_phenos_covariates.R: compile MDD phenotype and necessary covariates (in GS)

## analysis directory 
* analy_rvi_calculation_abcd.R: calculate MDD-RVI for each neuroimaging measure using the RVIpkg R package
* analy_rviprs_regression.R: conduct regression analysis to determine association between MDD phenotype and MDD-RVI or MDD-PRS
* analy_r2_aic_calculation.R: calculate variance explained in MDD (R^2) and model fit (AIC) by MDD-RVI and MDD-PRS when considered individually or additively
* analy_r2_rviprs_change_cbcl: to determine if the change in severity of depressive symptoms over two years is associated with MDD-RVI or MDD-PRS (in ABCD)

## functions directory 
(for larger functions; smaller functions are within the respective scripts) 
* gs_rename_enigma_FUN.R: to rename/abbreviate ENIGMA imaging variables to suit GS naming convention
* abcd_rename_enigma_FUN.R: to rename/abbreviate ENIGMA imaging variables to suit ABCD naming convention
