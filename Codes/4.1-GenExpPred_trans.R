# Genomic and/or expression-based prediction

# set seed
set.seed(123)

# source
library(data.table)
library(gdata)
library(BGLR)
library(ggplot2)

# param
args <- commandArgs(trailingOnly=T)
trait <- args[1]
n.rep <- 10
n.fold <- 5
nIter <- 12000
burnIn <- 8000

# mkdir
dir.save <- "RESULT/4.1-GenExpPred_trans"
dir.log <- "RESULT/4.1-GenExpPred_trans/Logfile"
dir.create(dir.log, recursive = T)

# myFun for CV
MyFun.CrossValid <- function(pheno.vec, ETA, CvFold.mat, n.rep, n.fold, ...) {
   # This is the object to save predicted values
   yPred.mat <- CvFold.mat
   for ( r in 1:n.rep ) {
      # cross validation fold
      CvNum <- CvFold.mat[, r]
      
      # vector to save predicted values
      yPred <- rep(NA, times = length(pheno.vec))
      for ( k in 1:n.fold ) {
         # show status
         print(paste0("Cross validation: ", k, "-th fold of ", n.fold, " folds in the ",  r, "-th rep in ", n.rep, " reps."))
         
         # mask test phenotype
         y.train <- pheno.vec
         tf.mask <- CvNum == k
         y.train[tf.mask] <- NA
         
         # regression
         fm <- BGLR(y = y.train, ETA = ETA, nIter = nIter, burnIn = burnIn, saveAt = paste0(dir.log, "/fm_"), verbose = F)
         yPred[tf.mask] <- fm$yHat[tf.mask]
      }
      
      # save result
      yPred.mat[, r] <- yPred
   }
   # return
   return(yPred.mat)
}

# myFun inverse box-cox
MyFun.InvBoxCox <- function(pred.vec, train, ...) {
   lambda <- BoxCoxParam.ames$lambda[BoxCoxParam.ames$trait == trait]
   const <- BoxCoxParam.ames$const[BoxCoxParam.ames$trait == trait]
   tmp.vec <- as.numeric(pred.vec)
   if (lambda == 0) { tmp.vec.original <- exp(tmp.vec) } else { tmp.vec.original <- tmp.vec ^ (1/lambda) }
   tmp.vec.original <- tmp.vec.original - const
   return(tmp.vec.original)   
}



# ----------------------------------------------------------------------------------------------- #
# ------------------------------------------ load data ------------------------------------------ #
# ----------------------------------------------------------------------------------------------- #
# genomic relationship matrix
GenReMat <- as.matrix(read.csv("RESULT/1.5-MakeGenReMat/GenReMat_ForExpData.csv", row.names = 1))

# phenotype data (trans)
PhenoData <- read.csv("RESULT/1.1-MakeAmesPhenoData/AmesPheno_trans.csv", stringsAsFactors = F)
PhenoData <- PhenoData[match(rownames(GenReMat), PhenoData$ID), ]
pheno.vec <- PhenoData[[trait]]
names(pheno.vec) <- PhenoData$ID

# phenotype data (un-trans)
PhenoData.ut <- read.csv("RESULT/1.1-MakeAmesPhenoData/AmesPheno.csv", stringsAsFactors = F)
PhenoData.ut <- PhenoData.ut[match(rownames(GenReMat), PhenoData.ut$ID), ]
pheno.vec.ut <- PhenoData.ut[[trait]]
names(pheno.vec.ut) <- PhenoData.ut$ID

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

# CvFold
CvFold <- read.csv("RESULT/1.6-MakeCvFold_Exp/CvFold_Exp.csv", row.names = 1, stringsAsFactors = F)
CvFold.mat <- as.matrix(CvFold[, -1])
rownames(CvFold.mat) <- CvFold$ID
rm(CvFold)

# a priori candidate gene 
CandGeneList <- read.xls("RAWDATA/tocochromanol_all_candidate_genes_combined_DW_20190624_with_two_genes.xlsx")
cand.gene.id <- unique(as.character(CandGeneList$RefGen_v4.Gene.ID))

# check the order of accessions
all.equal(rownames(ExpMat.scaled), rownames(GenReMat)) # OK!
all.equal(rownames(ExpMat.scaled), names(pheno.vec)) # OK!
all.equal(rownames(ExpMat.scaled), names(pheno.vec)) # OK!

# ---------------------------------------------------------------------------------------------- #
# ------------------------ relationship matrix based on expression data ------------------------ #
# ---------------------------------------------------------------------------------------------- #
# 1. Use all genes
ExpReMat.all <- tcrossprod(ExpMat.scaled) / ncol(ExpMat.scaled)

# 2. Use a priori candidate genes
cand.gene.id.obs <- intersect(cand.gene.id, colnames(ExpMat.scaled))
ExpMat.cand.gene <- ExpMat.scaled[, cand.gene.id.obs]
ExpReMat.cand <- tcrossprod(ExpMat.cand.gene) / ncol(ExpMat.cand.gene)


# -------------------------------------------------------------------------------------------------------------- #
# ------------------------------------------- 1. GBLUP as a baseline ------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------- #
# cross validation
ETA <- list("G" = list("K" = GenReMat, model = "RKHS"))
yPred.mat <- MyFun.CrossValid(pheno.vec, ETA = ETA, CvFold.mat, n.rep, n.fold)
yPred.mat.InvBoxCox <- apply(yPred.mat, 2, FUN = MyFun.InvBoxCox)

# save CV result!
df.yPred <- data.frame("ID" = rownames(yPred.mat),
                       "obs" = pheno.vec.ut,
                       yPred.mat.InvBoxCox)
file.out <- paste0(dir.save, "/pred_", trait, "_GBLUP.csv")
write.csv(df.yPred, file.out, row.names = F)

# correlation
cor.01 <- cor(df.yPred$obs, df.yPred[, 3:ncol(df.yPred)], use = "pair")



# ------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------- 2. Expression data: use all ---------------------------------------- #
# ------------------------------------------------------------------------------------------------------------- #
# model setup etc
ETA <- list("E" = list(K = ExpReMat.all, model = "RKHS"))
yPred.mat <- MyFun.CrossValid(pheno.vec, ETA = ETA, CvFold.mat, n.rep, n.fold)
yPred.mat.InvBoxCox <- apply(yPred.mat, 2, FUN = MyFun.InvBoxCox)

# save CV result!
df.yPred <- data.frame("ID" = rownames(yPred.mat),
                       "obs" = pheno.vec.ut,
                       yPred.mat.InvBoxCox)
file.out <- paste0(dir.save, "/pred_", trait, "_ExpBLUP.csv")
write.csv(df.yPred, file.out, row.names = F)

# correlation
cor.02 <- cor(df.yPred$obs, df.yPred[, 3:ncol(df.yPred)], use = "pair")



# ------------------------------------------------------------------------------------------------------------- #
# -------------------------------------- 3. Expression data: use subset --------------------------------------- #
# ------------------------------------------------------------------------------------------------------------- #
# model setup etc
ETA <- list("E" = list(K = ExpReMat.cand, model = "RKHS"))
yPred.mat <- MyFun.CrossValid(pheno.vec, ETA = ETA, CvFold.mat, n.rep, n.fold)
yPred.mat.InvBoxCox <- apply(yPred.mat, 2, FUN = MyFun.InvBoxCox)

# save CV result!
df.yPred <- data.frame("ID" = rownames(yPred.mat),
                       "obs" = pheno.vec.ut,
                       yPred.mat.InvBoxCox)
file.out <- paste0(dir.save, "/pred_", trait, "_ExpBLUP_UseCandGene.csv")
write.csv(df.yPred, file.out, row.names = F)

# correlation
cor.03 <- cor(df.yPred$obs, df.yPred[, 3:ncol(df.yPred)], use = "pair")



# -------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------- 4. GBLUP + EBLUP ---------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------- #
# cross validation
ETA <- list("G" = list(K = GenReMat, model = "RKHS"),
            "E" = list(K = ExpReMat.all, model = "RKHS"))
yPred.mat <- MyFun.CrossValid(pheno.vec, ETA = ETA, CvFold.mat, n.rep, n.fold)
yPred.mat.InvBoxCox <- apply(yPred.mat, 2, FUN = MyFun.InvBoxCox)

# save CV result!
df.yPred <- data.frame("ID" = rownames(yPred.mat),
                       "obs" = pheno.vec.ut,
                       yPred.mat.InvBoxCox)
file.out <- paste0(dir.save, "/pred_", trait, "_GBLUP+ExpBLUP.csv")
write.csv(df.yPred, file.out, row.names = F)

# correlation
cor.04 <- cor(df.yPred$obs, df.yPred[, 3:ncol(df.yPred)], use = "pair")




# -------------------------------------------------------------------------------------------------------------- #
# ------------------------------------------ 5. GBLUP + EBLUP_sub ---------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------- #
# cross validation
ETA <- list("G" = list(K = GenReMat, model = "RKHS"),
            "E" = list(K = ExpReMat.cand, model = "RKHS"))
yPred.mat <- MyFun.CrossValid(pheno.vec, ETA = ETA, CvFold.mat, n.rep, n.fold)
yPred.mat.InvBoxCox <- apply(yPred.mat, 2, FUN = MyFun.InvBoxCox)

# save CV result!
df.yPred <- data.frame("ID" = rownames(yPred.mat),
                       "obs" = pheno.vec.ut,
                       yPred.mat.InvBoxCox)
file.out <- paste0(dir.save, "/pred_", trait, "_GBLUP+ExpBLUP_UseCandGene.csv")
write.csv(df.yPred, file.out, row.names = F)

# correlation
cor.05 <- cor(df.yPred$obs, df.yPred[, 3:ncol(df.yPred)], use = "pair")
