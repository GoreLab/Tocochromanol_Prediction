# Genomic and/or expression-based prediction

# set.seed
set.seed(123)

# source
library(data.table)
library(gdata)
library(BGLR)
library(ggplot2)
library(Hmisc)

# param
nIter <-  60000
burnIn <- 40000
thin <- 20
trait.all <- c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3",
               "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")

# mkdir
dir.save <- "RESULT/4.2-GenExpFit_trans"
dir.log <- "RESULT/4.2-GenExpFit_trans/Logfile"
dir.create(dir.log, recursive = T)



# ----------------------------------------------------------------------------------------------- #
# ------------------------------------------ load data ------------------------------------------ #
# ----------------------------------------------------------------------------------------------- #
# genomic relationship matrix
GenReMat <- as.matrix(read.csv("RESULT/1.5-MakeGenReMat/GenReMat_ForExpData.csv", row.names = 1))

# phenotype data (trans)
PhenoData <- read.csv("RESULT/1.1-MakeAmesPhenoData/AmesPheno_trans.csv", stringsAsFactors = F)
PhenoData <- PhenoData[match(rownames(GenReMat), PhenoData$ID), ]

# load box-cox parameter
BoxCoxParam.ames <- read.csv("RESULT/1.1-MakeAmesPhenoData/Ames_BoxCoxParam.csv")

# gene expression data (all + vte7 expression)
Vte7ExprData <- fread("RAWDATA/ExprData_Ames/BLUE_vte7.csv", data.table = F)
ExpRawData <- fread("RAWDATA/ExprData_Ames/expression_BLUE_final_v1.1_B73.csv", data.table = F)

# merge vte7 with others
m <- match(ExpRawData$Accession_ID, Vte7ExprData$Genotype)
ExpRawData[["Zm00001d006778"]] <- Vte7ExprData$BLUE[m]
ExpMat <- as.matrix(ExpRawData[, -1])
rownames(ExpMat) <- gsub("_", "", ExpRawData$Accession_ID)
ExpMat <- ExpMat[match(rownames(GenReMat), rownames(ExpMat)), ]
ExpMat.scaled <- scale(ExpMat, center = T, scale = T) # scale each column

# a priori candidate gene 
CandGeneList <- read.xls("RAWDATA/tocochromanol_all_candidate_genes_combined_DW_20190624_with_two_genes.xlsx")
cand.gene.id <- unique(as.character(CandGeneList$RefGen_v4.Gene.ID))

# check the order of accessions
all.equal(rownames(ExpMat.scaled), rownames(GenReMat)) # OK!
all.equal(rownames(ExpMat.scaled), PhenoData$ID) # OK!



# ---------------------------------------------------------------------------------------------- #
# ------------------------ relationship matrix based on expression data ------------------------ #
# ---------------------------------------------------------------------------------------------- #
# 1. Use all genes
ExpReMat.all <- tcrossprod(ExpMat.scaled) / ncol(ExpMat.scaled)

# 2. Use a priori candidate genes
cand.gene.id.obs <- intersect(cand.gene.id, colnames(ExpMat.scaled))
ExpMat.cand.gene <- ExpMat.scaled[, cand.gene.id.obs]
ExpReMat.cand <- tcrossprod(ExpMat.cand.gene) / ncol(ExpMat.cand.gene)

# 3. Calculte gen-inverse (use it later)
Ginv.ExpMat.scaled <- MASS::ginv(ExpMat.scaled)
Ginv.ExpMat.cand.scaled <- MASS::ginv(ExpMat.cand.gene)



# ---------------------------------------------------------------------------------------------- #
# -------------------------- Expression-based prediction coefficient  -------------------------- #
# ---------------------------------------------------------------------------------------------- #
# for each trait
for ( i in 1:length(trait.all) ) {
   # trait
   trait <- trait.all[i]
   pheno.vec <- PhenoData[[trait]]
   names(pheno.vec) <- PhenoData$ID
   
   # model setup etc
   ETA <- list("E" = list(K = ExpReMat.all, model = "RKHS"))
   fm <- BGLR(y = pheno.vec, ETA = ETA, nIter = nIter, burnIn = burnIn, thin = thin,
              saveAt = paste0(dir.log, "/fm_ExpAll_", trait, "_"), verbose = F)
   saveRDS(fm, file = paste0(dir.log, "/fm_ExpAll_", trait, ".Rdata"))
   
   # calculate coefficients
   coef <- as.numeric(Ginv.ExpMat.scaled %*% fm$ETA$E$u)
   names(coef) <- colnames(ExpMat.scaled)
   
   # save data: coef and gene ID for all genes
   df.save <- data.frame("GeneID" = names(coef), "coef" = coef)
   fwrite(df.save, file = paste0(dir.save, "/Coef_ExpAll_", trait, ".csv"))
}



# ---------------------------------------------------------------------------------------------- #
# ------------------------------- G + E prediction coefficient  -------------------------------- #
# ---------------------------------------------------------------------------------------------- #
# for each trait
for ( i in 1:length(trait.all) ) {
   # trait
   trait <- trait.all[i]
   pheno.vec <- PhenoData[[trait]]
   names(pheno.vec) <- PhenoData$ID
   
   # model setup etc
   ETA <- list("G" = list(K = GenReMat, model = "RKHS"),
               "E" = list(K = ExpReMat.all, model = "RKHS"))
   fm <- BGLR(y = pheno.vec, ETA = ETA, nIter = nIter, burnIn = burnIn, thin = thin,
              saveAt = paste0(dir.log, "/fm_GBLUP+ExpAll_", trait, "_"), verbose = F)
   saveRDS(fm, file = paste0(dir.log, "/fm_GBLUP+ExpAll_", trait, ".Rdata"))
   
   # calculate coefficients
   coef <- as.numeric(Ginv.ExpMat.scaled %*% fm$ETA$E$u)
   names(coef) <- colnames(ExpMat.scaled)
   
   # save data: coef and gene ID for all genes
   df.save <- data.frame("GeneID" = names(coef), "coef" = coef)
   fwrite(df.save, file = paste0(dir.save, "/Coef_GBLUP+ExpAll_", trait, ".csv"))
}



# ---------------------------------------------------------------------------------------------- #
# ------------------- Expression-based prediction coefficient (candidates)  -------------------- #
# ---------------------------------------------------------------------------------------------- #
# for each trait
for ( i in 1:length(trait.all) ) {
   # trait
   trait <- trait.all[i]
   pheno.vec <- PhenoData[[trait]]
   names(pheno.vec) <- PhenoData$ID
   
   # model setup etc
   ETA <- list("E" = list(K = ExpReMat.cand, model = "RKHS"))
   fm <- BGLR(y = pheno.vec, ETA = ETA, nIter = nIter, burnIn = burnIn, thin = thin,
              saveAt = paste0(dir.log, "/fm_ExpCand_", trait, "_"), verbose = F)
   saveRDS(fm, file = paste0(dir.log, "/fm_ExpCand_", trait, ".Rdata"))
   
   # calculate coefficients
   coef <- as.numeric(Ginv.ExpMat.cand.scaled %*% fm$ETA$E$u)
   names(coef) <- colnames(ExpMat.cand.gene)
   
   # save data: coef and gene ID for all genes
   df.save <- data.frame("GeneID" = names(coef), "coef" = coef)
   fwrite(df.save, file = paste0(dir.save, "/Coef_ExpCand_", trait, ".csv"))
}



# ---------------------------------------------------------------------------------------------- #
# ------------------------ G + E prediction coefficient (candidates)  -------------------------- #
# ---------------------------------------------------------------------------------------------- #
# for each trait
for ( i in 1:length(trait.all) ) {
   # trait
   trait <- trait.all[i]
   pheno.vec <- PhenoData[[trait]]
   names(pheno.vec) <- PhenoData$ID
   
   # model setup etc
   ETA <- list("G" = list(K = GenReMat, model = "RKHS"),
               "E" = list(K = ExpReMat.cand, model = "RKHS"))
   fm <- BGLR(y = pheno.vec, ETA = ETA, nIter = nIter, burnIn = burnIn, thin = thin,
              saveAt = paste0(dir.log, "/fm_GBLUP+ExpCand_", trait, "_"), verbose = F)
   saveRDS(fm, file = paste0(dir.log, "/fm_GBLUP+ExpCand_", trait, ".Rdata"))
   
   # calculate coefficients
   coef <- as.numeric(Ginv.ExpMat.cand.scaled %*% fm$ETA$E$u)
   names(coef) <- colnames(ExpMat.cand.gene)
   
   # save data: coef and gene ID for all genes
   df.save <- data.frame("GeneID" = names(coef), "coef" = coef)
   fwrite(df.save, file = paste0(dir.save, "/Coef_GBLUP+ExpCand_", trait, ".csv"))
}






# # ---------------------------------------------------------------------------- #
# # library
# library(rrBLUP)
# 
# # trait
# trait <- trait.all[1]
# pheno.vec <- PhenoData[[trait]]
# names(pheno.vec) <- PhenoData$ID
# 
# # use GBLUP and then calculate coefficients
# ms1 <- mixed.solve(y = pheno.vec, K = ExpReMat.all)
# coef1 <- Ginv.ExpMat.scaled %*% ms1$u
# 
# # use rrBLUP
# ms2 <- mixed.solve(y = pheno.vec, Z = ExpMat.scaled, K = diag(ncol(ExpMat.scaled)))
# coef2 <- ms2$u
# 
# # identical
# plot(coef1, coef2, pch = 20)
# abline(0, 1)
# 
# # the other way (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8495923/)
# coef1.v2 <- (1 / ncol(ExpMat.scaled)) * t(ExpMat.scaled) %*% MASS::ginv(ExpReMat.all) %*% ms1$u
# 
# # identical
# plot(coef1, coef1.v2)
# abline(0, 1)
# # ---------------------------------------------------------------------------- #





