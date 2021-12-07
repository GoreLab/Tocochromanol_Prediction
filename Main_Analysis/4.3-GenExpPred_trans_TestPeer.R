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


# mkdir
dir.save <- "RESULT/4.3-GenExpPred_trans_TestPeer"
dir.log <- "RESULT/4.3-GenExpPred_trans_TestPeer/Logfile"
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

# loop for all K
for ( K in 1:25 ) {
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
   
   # CvFold
   CvFold <- read.csv("RESULT/1.6-MakeCvFold_Exp/CvFold_Exp.csv", row.names = 1, stringsAsFactors = F)
   CvFold.mat <- as.matrix(CvFold[, -1])
   rownames(CvFold.mat) <- CvFold$ID
   rm(CvFold)
   
   # a priori candidate gene 
   CandGeneList <- read.xls("RAWDATA/tocochromanol_all_candidate_genes_combined_DW_20190624.xlsx")
   cand.gene.id <- unique(as.character(CandGeneList$RefGen_v4.Gene.ID))
   
   
   # --- Use PEER instead of BLUEs
   # load gene expression data
   ExpRawData <- fread("RAWDATA/ExprData_Ames/expression_BLUE_final_v1.1_B73.csv")
   ExpMat <- as.matrix(ExpRawData[, -1])
   rownames(ExpMat) <- gsub("_", "", ExpRawData$Accession_ID)
   ExpMat <- ExpMat[match(rownames(GenReMat), rownames(ExpMat)), ]
   
   # load PEER result K = 25
   w.df <- fread("RAWDATA/ExprData_Ames/PEER/PeerResult_Use25Fact_v1.1_B73_weights.txt", data.table = F)
   f.df <- fread("RAWDATA/ExprData_Ames/PEER/PeerResult_Use25Fact_v1.1_B73_factors.txt", data.table = F)
   w.mat <- as.matrix(w.df[, -1]); rownames(w.mat) <- w.df[, 1]
   f.mat <- as.matrix(f.df[, -1]); rownames(f.mat) <- f.df[, 1]
   
   
   # my function to calculate peer residuals by using 1~K factors
   MyFun.Calc.PEER.resid <- function(K, ExpMat, f.mat, w.mat,...) {
      fw.mat <- f.mat[, 1:K, drop = F] %*% t(w.mat[, 1:K, drop = F])
      rownames(fw.mat) <- gsub("_", "", rownames(fw.mat))
      m <- match(rownames(ExpMat), rownames(fw.mat))
      fw.mat.sorted <- fw.mat[m, ]
      PEER.resid.mat <- ExpMat - fw.mat.sorted
      return(PEER.resid.mat)
   }
   PEER.resid <- MyFun.Calc.PEER.resid(K, ExpMat, f.mat, w.mat)
   PEER.resid.scaled <- scale(PEER.resid, center = T, scale = T) # scale each column
   rm(list = c("w.df", "f.df", "w.mat", "f.mat")); gc(); gc()
   rm(list = c("ExpRawData", "ExpMat")); gc(); gc()
   
   # check the order of accessions
   all.equal(rownames(PEER.resid.scaled), rownames(GenReMat)) # OK!
   all.equal(rownames(PEER.resid.scaled), names(pheno.vec)) # OK!
   all.equal(rownames(PEER.resid.scaled), names(pheno.vec)) # OK!
   
   # ---------------------------------------------------------------------------------------------- #
   # ------------------------ relationship matrix based on expression data ------------------------ #
   # ---------------------------------------------------------------------------------------------- #
   # 1. Use all genes
   ExpReMat.all <- tcrossprod(PEER.resid.scaled) / ncol(PEER.resid.scaled)
   
   # 2. Use a priori candidate genes
   cand.gene.id.obs <- intersect(cand.gene.id, colnames(PEER.resid.scaled))
   ExpMat.cand.gene <- PEER.resid.scaled[, cand.gene.id.obs]
   ExpReMat.cand <- tcrossprod(ExpMat.cand.gene) / ncol(ExpMat.cand.gene)
   
   
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
   file.out <- paste0(dir.save, "/pred_", trait, "_ExpBLUP_UsePeer_K", K, ".csv")
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
   file.out <- paste0(dir.save, "/pred_", trait, "_ExpBLUP_UseCandGene_UsePeer_K", K, ".csv")
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
   file.out <- paste0(dir.save, "/pred_", trait, "_GBLUP+ExpBLUP_UsePeer_K", K, ".csv")
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
   file.out <- paste0(dir.save, "/pred_", trait, "_GBLUP+ExpBLUP_UseCandGene_UsePeer_K", K, ".csv")
   write.csv(df.yPred, file.out, row.names = F)
   
   # correlation
   cor.05 <- cor(df.yPred$obs, df.yPred[, 3:ncol(df.yPred)], use = "pair")
}






# # test the function
# r.df <- fread("RAWDATA/ExprData_Ames/PEER/PeerResult_Use25Fact_v1.1_B73_residuals.txt", data.table = F)
# r.mat <- as.matrix(r.df[, -1]); rownames(r.mat) <- r.df[, 1]
# PEER.resid.mat <- MyFun.Calc.PEER.resid(K = 25, ExpMat, f.mat, w.mat)
# rownames(r.mat) <- gsub("_", "", rownames(r.mat))
# m <- match(rownames(PEER.resid.mat), rownames(r.mat))
# diff.vec <- as.numeric(PEER.resid.mat - r.mat[m, ])
# hist(diff.vec, breaks = 100)
