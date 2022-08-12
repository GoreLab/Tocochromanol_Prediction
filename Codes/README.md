## Overview 
This folder has all R codes to perform the analysis in "paper_title". All data analyses were performed in R version 3.5.0 in linux OS (high-spec machine for computation), while a few post-hoc visualization (such as making plots based on the predicted values) were executed in R version 4.1.2 in windows OS (a laptop). All codes were developed by Ryokei Takana, expect for the "0.1-MyFun_BoxcoxFunction.R" that was originally developed by Christine Diepenbrock. 

## Codes to create data files from raw data files for the downstream statistical analyses
#### 0.1-MyFun_BoxcoxFunction.R
This R code has a function for the Box-Cox transformation, including an optimization of the convinient lambda value. This function was used in Wu et al. (in prep) for the tocochromanol phenotypes in Ames panel, and we applied the same method to the phenotypes in the Goodman panel.

#### 0.3-Functions_for_MTM.R
This R code has a function to run the multi-trait model, which is used in "4.4-Mtm_UseSomeGenes.R"

#### 1.1-MakeAmesPhenoData.R / 1.2-MakeGbPhenoData.R
These two R codes create files of tocochromanol phenotypes: transformed and untransformed BLUE (for Ames panel) or BLUP (for the Goodman panel) values, convinient lambda values, and the sample ID correspondence between phenotype and genotype data.

#### 1.3-Ibs.R
This code was developed fro the identification of overlapping lines between the Ames and Goodman panels. Lines having nearly identical line names between the two panels were paired, and their IBS was calculated. If IBS < 0.8, they were regarded as different lines. See Materials and Method for more details.

#### 1.4-MakeCvFold.R / 1.6-MakeCvFold_Exp.R
The former code "1.4-MakeCvFold.R" creates a fold (i.e., randomly created K groups, where K = 5 in this study) for the cross-validation with X replications (X = 10 in this study) for the Ames and the Goodman panels. Similarly, "1.6-MakeCvFold_Exp.R" creates the cross-validation fold for the 545 accessions for the transcriptome-based prediction analysis.

#### 1.5-MakeGenReMat.R
This code calculates genomic relationship matrix for (1) all accessions either in the Ames or Goodman panel, (2) accessions in the Ames panel, (3) accessions in the Goodman panel, and (4) the 545 Ames accessions for the transcriptome-based predcition analysis. The genomic relationship matrix was calculated in VanRaden's method 1 (VanRaden, 2008).


## Codes for the two whole-genome prediction models (GBLUP and BayesB)
#### 2.1-GenPreFit_SinglePop_trans.R / 2.2-PredictionAccuracy_trans.R
The former code "2.1-GenPreFit_SinglePop_trans.R fits the BayesB model on the training population (either Ames or Goodman panel) and estimates the marker effects using BGLR package. The latter code "2.2-PredictionAccuracy_trans.R" uses the estimated marker effects (and grand mean) to calculate the predicted values in the test population, and apply the back-transformation (i.e., inverse Box-Cox transfomration) with the convinient lambda values of the training population for each tocochromanol phenotype.

Note that the name of the latter code is misleading, as it does not calculate the prediction accuracy. Originally, the code was implemented to calculate the predictive ability, but we decided to evaluate the predictive ability using the non-overlapping lines.

#### 2.3-GenPreFit_SinglePop_trans_GBLUP.R
This code uses GBLUP to predict (1) the Goodman panel from the Ames panel and (2) the Ames panel from the Goodman panel. The predicted genotypic values were back-transformed as done in "2.2-PredictionAccuracy_trans.R" for the result from BayesB.

#### 3.1-GenPreCv_SinglePop_trans.R / 3.2-GenPreCv_SinglePop_trans_GBLUP.R
These codes perform the cross-validation: the former code "3.1-GenPreCv_SinglePop_trans.R" uses the BayesB model, while the latter uses the GBLUP model.


## Codes for the transcriptome-based prediction
#### 4.1-GenExpPred_trans.R
This code perform the cross-validation with five models: (1) GBLUP as baseline, (2) TBLUP with all 22137 genes, (3) TBLUP with 111 a priori pathway genes, (4) GBLUP + TBLUP with all 22137 genes, and (5) GBLUP + TBLUP with 111 a priori pathway genes. 

#### 4.2-GenExpFit_trans.R
This code fits the two GBLUP + TBLUP models to the entire 545 lines, and calculates the regression coefficients of the genes according to the formula shown in Zhang et al. (2021). The estimated regression coefficients were then used for the principle component analysis etc.

#### 4.3-MakeDataForMtm.R / 4.4-Mtm_UseSomeGenes.R
These two codes applies the multi-trait model as one of our transcriptome-based prediction models. The first script "4.3-MakeDataForMtm.R" creates an R object that includes minimal data for the multi-trait prediction (this step is not necessary to be separately done, but I made this code for simplisity). The latter code applies the multi-trait model using the MTM package.


## Codes for the multi-kernel GBLUP based on the NAM JL-QTL
#### 5.1-MultiKernel_Across_trans.R
This code applies the multi-kernel GBLUP model for the two predcition scenarios: (1) from the Ames to Goodman and (2) from the Goodman to Ames. 

#### 5.2-MultiKernel_Cv_trans
This code applies the multi-kernel GBLUP model for the cross-validation within Ames or Goodman


## Codes for summarizing the result: calculate predictive ability, make figures and tables
#### 99.01-MakeSummaryFiles.R
This code makes files of the predicted tocohromanol phenotypes (i.e., back-transformed values), using the result files from the previous codes. 

#### 99.02-MakeFigAndTables_UseNamQtl_NonOverlap.R
This code makes figures and tables of the multi-kernel prediction analysis using NAM QTL. [OLD CODE AND NOT USED FOR THE FINAL FIGURES]

#### 99.03-MakeFigsAndTables_ExprPred.R
This code makes figures and tables of the prediction analysis using transcripts [OLD CODE AND NOT USED FOR THE FINAL FIGURES]

#### 99.04-PCA_of_Importance_ExprPred.R
This code performs PCA for the regression coefficients of the transcripts and visualize the result in the PCA bi-plot [OLD CODE AND NOT USED FOR THE FINAL FIGURES]

#### 99.05-Panel_plot_Importance_ExprPred.R
This code visualizes the importance of NAM QTL with their PVE for each tocochromanol phenotype [This makes a prototype of a figure]

#### 999-MakeFigs.R
This code generates updated figures for publication (e.g., chage colors, sizes, legends)

## Codes for genotype processing
#### Ames_282_GP_imputation_filtering_pruning.sh
This code uses in the genotypes of both Ames and Goodman panel genotyped with GBS, and imputed onto Hapmap v3.2.1 independently for two panels. Afterwards, the two genotypes are merged and filtered based on imputation quality and MAF, and pruned based on LD of sliding windows.

