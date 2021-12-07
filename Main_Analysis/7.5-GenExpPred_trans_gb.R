# Genomic and/or expression-based prediction

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

# ###########################################################
# trait <- "d.T"
# n.rep <- 10
# n.fold <- 5
# nIter <-  600
# burnIn <- 400
# ###########################################################

# mkdir
dir.save <- "RESULT/7.5-GenExpPred_trans_gb"
dir.log <- "RESULT/7.5-GenExpPred_trans_gb/Logfile"
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
GenReMat <- as.matrix(read.csv("RESULT/7.3-MakeGenReMat_ForGbExprData/GenReMat_gb_expr.csv", row.names = 1))

# phenotype data (trans)
PhenoData <- read.csv("RESULT/1.2-MakeGbPhenoData/GbPheno_trans.csv", stringsAsFactors = F)
PhenoData <- PhenoData[match(rownames(GenReMat), PhenoData$ID), ]
pheno.vec <- PhenoData[[trait]]
names(pheno.vec) <- PhenoData$ID

# phenotype data (un-trans)
PhenoData.ut <- read.csv("RESULT/1.2-MakeGbPhenoData/GbPheno.csv", stringsAsFactors = F)
PhenoData.ut <- PhenoData.ut[match(rownames(GenReMat), PhenoData.ut$ID), ]
pheno.vec.ut <- PhenoData.ut[[trait]]
names(pheno.vec.ut) <- PhenoData.ut$ID

# load box-cox parameter
BoxCoxParam.ames <- read.csv("RESULT/1.2-MakeGbPhenoData/Gb_BoxCoxParam.csv")

# load Gb info data
InfoData <- read.csv("RESULT/1.2-MakeGbPhenoData/GbInfo.csv")

# gene expression data
ExpRawData <- fread("RESULT/7.4-Diagnosis_ExprData/ExpressionData.Filtered.and.Imputed.csv")
ExpMat <- as.matrix(ExpRawData[, 2:ncol(ExpRawData)])
rownames(ExpMat) <- ExpRawData$Sample_ID
m <- match(rownames(ExpMat), InfoData$Name.Lipka)
rownames(ExpMat) <- InfoData$GBS.taxa[m]
ExpMat <- ExpMat[match(rownames(GenReMat), rownames(ExpMat)), ]
ExpMat.scaled <- scale(ExpMat, center = T, scale = T) # scale each column

# make CvFold
set.seed(637)
n.fold <- 5; n.rep <- 10
myfun.MakeCvmat <- function(N, n.fold, n.rep) {
	CvMat <- matrix(NA, nr = N, nc = n.rep)
	for (r in 1:n.rep) {
		num <- rep(1:n.fold, length.out = N)
		rand.num <- sample(num, replace = FALSE)
		CvMat[, r] <- rand.num
	}
	colnames(CvMat) <- paste0("rep", formatC(1:n.rep, width = 2, flag = "0"))
	return(CvMat)
}
CvMat <- myfun.MakeCvmat(N = nrow(ExpMat.scaled), n.fold = n.fold, n.rep = n.rep)
df.CvMat <- data.frame("ID" = rownames(ExpMat.scaled), CvMat)
write.csv(df.CvMat, file = paste0(dir.save, "/CvFold.csv")) # save

# rename the CV data (to reuse the previous code)
CvFold.mat <- as.matrix(df.CvMat[, -1])
rownames(CvFold.mat) <- df.CvMat$ID

# a priori candidate gene 
CandGeneList <- read.xls("RAWDATA/tocochromanol_all_candidate_genes_combined_DW_20190624.xlsx")
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

num <- which(apply(ExpMat, 2, sd) == 0)
apply(ExpMat[, num], 2, sum)



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
