# multi-kernel regression

#
source("0.2-MultiKernel_MyFunToSelectSnps.R")
library(data.table)
library(gdata)
library(rrBLUP)
library(BGLR)

# param
args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]
pop <- args[2]
bp.width.qtl <- as.integer(args[3]) # set 0 if you do not use this
bp.width.gene <- as.integer(args[4]) # set 0 if you do not use this
name.file <- args[5]
nIter <-  12000
burnIn <- 8000
thin <- 5
n.rep <- 10
n.fold <- 5

# ####
# trait <- "d.T3"
# pop <- "gb"
# bp.width.qtl <- 1000000
# bp.width.gene <- 0
# name.file <- "exampletest"
# nIter <-  100
# burnIn <- 10
# thin <- 5
# n.rep <- 2
# n.fold <- 2
# ###


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

# myFun for Inv-Box-COx
MyFun.InvBoxCox <- function(pred.vec, ...) {
   lambda <- BoxCoxParam$lambda[BoxCoxParam$trait == trait]
   const <- BoxCoxParam$const[BoxCoxParam$trait == trait]
   tmp.vec <- as.numeric(pred.vec)
   if (lambda == 0) { tmp.vec.original <- exp(tmp.vec) } else { tmp.vec.original <- tmp.vec ^ (1/lambda) }
   tmp.vec.original <- tmp.vec.original - const
   return(tmp.vec.original)   
}

# mkdir
dir.save <- "RESULT/5.2-MultiKernel_Cv_trans"
dir.log <- "RESULT/5.2-MultiKernel_Cv_trans/Logfile"
dir.create(dir.log, recursive = T)



# ----------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------- load genotype data  ---------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------------- #
# load geotype & map
indv.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.indv", header = FALSE)
geno.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012")
geno.mat <- as.matrix(geno.tmp[, -1]) # 1st col is ID.
rownames(geno.mat) <- indv.tmp$V1
map.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.pos")
map <- data.frame("chr" = map.tmp$V1, "pos" = map.tmp$V2)
rm(indv.tmp); rm(geno.tmp); rm(map.tmp); gc(); gc() # remove non-essential object

# knwon QTLs
GenData <- read.csv("RAWDATA/NAM_QTL/qtl.data.csv", stringsAsFactors = FALSE)
colnames(GenData)[1] <- "QTL.ID"

# known genes
TocoGene <- read.xls("RAWDATA/NAM_QTL/TocoGene.xlsx", stringsAsFactors = FALSE)
colnames(TocoGene)[9:ncol(TocoGene)] <- c("a.T", "d.T", "g.T", "Total.Tocopherols",
                                          "a.T3", "d.T3", "g.T3", "Total.Tocotrienols",
                                          "Total.Tocochromanols") # rename

# Get SNPs: choose appopriate method
if (bp.width.qtl != 0 & bp.width.gene == 0) {
   print("Use only QTLs as there is no specification of window size for a priori genes.")
   SnpNumList <- MyFun.GetSnp.QTL(trait = trait, bp.width.qtl = bp.width.qtl, GenData = GenData, map = map)
}
if (bp.width.qtl == 0 & bp.width.gene != 0) {
   print("Use only a priori genes as there is no specification of window size for QTLs.")
   SnpNumList <- MyFun.GetSnp.GENE(trait = trait, bp.width.gene = bp.width.gene, TocoGene = TocoGene, map = map)
}
if (bp.width.qtl == 0 & bp.width.gene == 0) {
   print("Use support intervals as there is no specification of window size.")
   SnpNumList <- MyFun.GetSnp.SI(trait = trait, GenData = GenData, map = map)
}
if (bp.width.qtl != 0 & bp.width.gene != 0) {
   print("Use both QTLs and a priori genes: hybrid method.")
   SnpNumList <- MyFun.GetSnp.HYBRID(trait = trait, bp.width.qtl = bp.width.qtl, bp.width.gene = bp.width.gene, 
                                     TocoGene = TocoGene, GenData = GenData, map = map)
}

# save the object of SNPs
saveRDS(SnpNumList, file = paste0(dir.save, "/SnpNumList_", name.file, "_", trait, ".Rdata"))



# ------------------------------------------------------------------------------------------------------------------ #
# ---------------------------------------------- load phenotype data  ---------------------------------------------- #
# ------------------------------------------------------------------------------------------------------------------ #
# filenames to load
if ( pop == "ames" ) {
   Input.CvFold <- "RESULT/1.4-MakeCvFold/CvFold_ames.csv"
   Input.BoxCox <- "RESULT/1.1-MakeAmesPhenoData/Ames_BoxCoxParam.csv"
   input.Pheno.tr <- "RESULT/1.1-MakeAmesPhenoData/AmesPheno_trans.csv"
   input.Pheno.ut <- "RESULT/1.1-MakeAmesPhenoData/AmesPheno.csv"
}
if ( pop == "gb" ) {
   Input.CvFold <- "RESULT/1.4-MakeCvFold/CvFold_gb.csv"
   Input.BoxCox <- "RESULT/1.2-MakeGbPhenoData/Gb_BoxCoxParam.csv"
   input.Pheno.tr <- "RESULT/1.2-MakeGbPhenoData/GbPheno_trans.csv"
   input.Pheno.ut <- "RESULT/1.2-MakeGbPhenoData/GbPheno.csv"
}

# CvFold
CvFold <- read.csv(Input.CvFold, row.names = 1, stringsAsFactors = F)
CvFold.mat <- as.matrix(CvFold[, -1])
rownames(CvFold.mat) <- CvFold$ID
rm(CvFold)

# load box-cox parameter
BoxCoxParam <- read.csv(Input.BoxCox)

# phenotype data (trans/un-trans)
PhenoData <- read.csv(input.Pheno.tr, stringsAsFactors = F)
PhenoData.ut <- read.csv(input.Pheno.ut, stringsAsFactors = F)



# --------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------- Cross Validation  ---------------------------------------------- #
# --------------------------------------------------------------------------------------------------------------- #
# take subset of genotype
pheno.vec <- PhenoData[[trait]]
pheno.vec.ut <- PhenoData.ut[[trait]]
m <- match(PhenoData$ID, rownames(geno.mat))

# make kernel -> Model object
ETA <- list()
for (i in 1:length(SnpNumList)) {
   GenReMat.i <- A.mat(X = geno.mat[m, SnpNumList[[i]]] - 1, min.MAF = NULL, shrink = FALSE)
   ETA[[i]] <- list("K" = GenReMat.i, "model" = "RKHS")
}
names(ETA) <- names(SnpNumList)

# cross validation
yPred.mat <- MyFun.CrossValid(pheno.vec, ETA = ETA, CvFold.mat, n.rep, n.fold)
yPred.mat.InvBoxCox <- apply(yPred.mat, 2, FUN = MyFun.InvBoxCox)

# save CV result!
df.yPred <- data.frame("ID" = rownames(yPred.mat),
                       "obs" = pheno.vec.ut,
                       yPred.mat.InvBoxCox)
file.out <- paste0(dir.save, "/pred_", pop, "_", trait, "_", name.file, ".csv")
write.csv(df.yPred, file.out, row.names = F)



