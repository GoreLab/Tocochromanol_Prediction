## Overview 
This folder has all R codes to perform the analysis in "paper_title". All data analyses were performed in R version 3.5.0 in linux OS (high-spec machine for computation), while a few post-hoc visualization (such as making plots based on the predicted values) were executed in R version 4.1.2 in windows OS. All codes were developed by Ryokei Takana, expect for the "0.1-MyFun_BoxcoxFunction.R" that was originally developed by Christine Diepenbrock.


#### 0.1-MyFun_BoxcoxFunction.R
This R code has a function for the Box-Cox transformation, including an optimization of the convinient lambda value. This function was used in Wu et al. (in prep) for the tocochromanol phenotypes in Ames panel, and therfore we applied the same method to the phenotypes in the 282 panel.

#### 1.1-MakeAmesPhenoData.R and 1.2-MakeGbPhenoData.R
These two R codes creates files of tocochromanol phenotypes: transformed and untransformed BLUE (for Ames panel) or BLUP (for the 282 panel) values, convinient lambda values, and the sample ID correspondence between phenotype and genotype data.

#### 1.4-MakeCvFold.R
This code makes a fold (i.e., randomly created K groups, where K = 5 in this study) for the cross-validation with X replication (X = 10 in this study). This cross-validation fold was created for each of the Ames and 282 panels, and used for all tocochromanol phenotypes and different models. Note that the transcriptome-based prediction was performed on a subset of the Ames panel and therefore we used another R code (1.6-MakeCvFold_Exp.R) for that analysis.

#### 1.5 MakeGenReMat.R


