## Overview 
This folder has all R codes to perform the analysis in "paper_title". All data analyses were performed in R version 3.5.0 in linux OS (high-spec machine for computation), while a few post-hoc visualization (such as making plots based on the predicted values) were executed in R version 4.1.2 in windows OS. All codes were developed by Ryokei Takana, expect for the "0.1-MyFun_BoxcoxFunction.R" that was originally developed by Christine Diepenbrock.


#### 0.1-MyFun_BoxcoxFunction.R
This R code has a function for the Box-Cox transformation, including an optimization of the convinient lambda value. This function was used in Wu et al. (in prep) for the tocochromanol phenotypes in Ames panel, and therfore we applied the same method to the phenotypes in the 282 panel.

#### 1.1-MakeAmesPhenoData.R / 1.2-MakeGbPhenoData.R
These two R codes create files of tocochromanol phenotypes: transformed and untransformed BLUE (for Ames panel) or BLUP (for the 282 panel) values, convinient lambda values, and the sample ID correspondence between phenotype and genotype data.

#### 1.4-MakeCvFold.R / 1.6-MakeCvFold_Exp.R
The former code "1.4-MakeCvFold.R" creates a fold (i.e., randomly created K groups, where K = 5 in this study) for the cross-validation with X replications (X = 10 in this study) for the Ames and the 282 panels. Similarly, "1.6-MakeCvFold_Exp.R" creates the cross-validation fold for the 545 accessions for the transcriptome-based prediction analysis.

#### 1.5 MakeGenReMat.R
This code calculates genomic relationship matrix for (1) all accessions either in the Ames or 282 panel, (2) accessions in the Ames panel, (3) accessions in the 282 panel, and (4) the 545 Ames accessions for the transcriptome-based predcition analysis. The genomic relationship matrix was calculated in VanRaden's method 1 (VanRaden, 2008).


