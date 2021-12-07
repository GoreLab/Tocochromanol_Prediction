# Genomic and/or expression-based prediction
# with functional features in the Hapmap population

# source
library(data.table)
library(gdata)
library(BGLR)
library(ggplot2)
library(stringr)

# param
args <- commandArgs(trailingOnly=T)
trait <- args[1]
pop <- args[2]
name.feature <- args[3]
n.rep <- 10
n.fold <- 5
nIter <- 6000
burnIn <- 4000

# mkdir
dir.save <- "RESULT/6.2-UseFeatures_cv"
dir.log <- "RESULT/6.2-UseFeatures_cv/Logfile"
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
         b <- Sys.time()
         
         # mask test phenotype
         y.train <- pheno.vec
         tf.mask <- CvNum == k
         y.train[tf.mask] <- NA
         
         # regression
         fm <- BGLR(y = y.train, ETA = ETA, nIter = nIter, burnIn = burnIn, saveAt = paste0(dir.log, "/fm_"), verbose = F)
         yPred[tf.mask] <- fm$yHat[tf.mask]
         
         # print time duration
         a <- Sys.time()
         print(a-b)
      }
      
      # save result
      yPred.mat[, r] <- yPred
   }
   # return
   return(yPred.mat)
}

# myFun inverse box-cox
MyFun.InvBoxCox <- function(pred.vec, train, ...) {
   lambda <- BoxCoxParam$lambda[BoxCoxParam$trait == trait]
   const <- BoxCoxParam$const[BoxCoxParam$trait == trait]
   tmp.vec <- as.numeric(pred.vec)
   if (lambda == 0) { tmp.vec.original <- exp(tmp.vec) } else { tmp.vec.original <- tmp.vec ^ (1/lambda) }
   tmp.vec.original <- tmp.vec.original - const
   return(tmp.vec.original)   
}


# ----------------------------------------------------------------------------------------------- #
# ------------------------------------------ load data ------------------------------------------ #
# ----------------------------------------------------------------------------------------------- #
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


# ------------------------------------------------------------------------------------------ #
# ------------------------ relationship matrix based on the feature ------------------------ #
# ------------------------------------------------------------------------------------------ #
pop2 <- str_to_title(pop)
# 0. model setup
if ( name.feature == "Prox_100kb" ) {
   GenReMat_prox <- as.matrix(read.csv(paste0("RESULT/1.7-UseFeatures/GenReMat_For", pop2, "Data_proximity_100kb.csv"), row.names = 1))
   GenReMat_non <- as.matrix(read.csv(paste0("RESULT/1.7-UseFeatures/GenReMat_For", pop2, "Data_non_proximity_100kb.csv"), row.names = 1))
   ETA <- list("G1" = list(K = GenReMat_prox, model = "RKHS"),
               "G2" = list(K = GenReMat_non, model = "RKHS"))
}
if ( name.feature == "Prox_10kb" ) {
   GenReMat_prox <- as.matrix(read.csv(paste0("RESULT/1.7-UseFeatures/GenReMat_For", pop2, "Data_proximity_10kb.csv"), row.names = 1))
   GenReMat_non <- as.matrix(read.csv(paste0("RESULT/1.7-UseFeatures/GenReMat_For", pop2, "Data_non_proximity_10kb.csv"), row.names = 1))
   ETA <- list("G1" = list(K = GenReMat_prox, model = "RKHS"),
               "G2" = list(K = GenReMat_non, model = "RKHS")) # model
}
if ( name.feature == "Prox_1kb" ) {
   GenReMat_prox <- as.matrix(read.csv(paste0("RESULT/1.7-UseFeatures/GenReMat_For", pop2, "Data_proximity_1kb.csv"), row.names = 1))
   GenReMat_non <- as.matrix(read.csv(paste0("RESULT/1.7-UseFeatures/GenReMat_For", pop2, "Data_non_proximity_1kb.csv"), row.names = 1))
   ETA <- list("G1" = list(K = GenReMat_prox, model = "RKHS"),
               "G2" = list(K = GenReMat_non, model = "RKHS")) # model
}
if ( name.feature == "MnaseShoot" ) {
   GenReMat_hot <- as.matrix(read.csv(paste0("RESULT/1.7-UseFeatures/GenReMat_For", pop2, "Data_MNASE_SHOOT_hotspot.csv"), row.names = 1))
   GenReMat_non <- as.matrix(read.csv(paste0("RESULT/1.7-UseFeatures/GenReMat_MNASE_For", pop2, "Data_SHOOT_non_hotspot.csv"), row.names = 1))
   ETA <- list("G1" = list(K = GenReMat_hot, model = "RKHS"),
               "G2" = list(K = GenReMat_non, model = "RKHS")) # model
}
if ( name.feature == "MnaseRoot" ) {
   GenReMat_hot <- as.matrix(read.csv(paste0("RESULT/1.7-UseFeatures/GenReMat_For", pop2, "Data_MNASE_ROOT_hotspot.csv"), row.names = 1))
   GenReMat_non <- as.matrix(read.csv(paste0("RESULT/1.7-UseFeatures/GenReMat_For", pop2, "Data_MNASE_ROOT_non_hotspot.csv"), row.names = 1))
   ETA <- list("G1" = list(K = GenReMat_hot, model = "RKHS"),
               "G2" = list(K = GenReMat_non, model = "RKHS")) # model
}
if ( name.feature == "GERP" ) {
   GenReMat_GERP_posi <- as.matrix(read.csv(paste0("RESULT/1.7-UseFeatures/GenReMat_For", pop2, "Data_positive_GERP.csv"), row.names = 1))
   GenReMat_GERP_nega <- as.matrix(read.csv(paste0("RESULT/1.7-UseFeatures/GenReMat_For", pop2, "Data_nagative_GERP.csv"), row.names = 1))
   ETA <- list("G1" = list(K = GenReMat_GERP_posi, model = "RKHS"),
               "G2" = list(K = GenReMat_GERP_nega, model = "RKHS")) # model
}
if ( name.feature == "RecRate" ) {
   GenReMat_RecRate_high <- as.matrix(read.csv(paste0("RESULT/1.7-UseFeatures/GenReMat_For", pop2, "Data_high_recombination.csv"), row.names = 1))
   GenReMat_RecRate_middle <- as.matrix(read.csv(paste0("RESULT/1.7-UseFeatures/GenReMat_For", pop2, "Data_middle_recombination.csv"), row.names = 1))
   GenReMat_RecRate_low <- as.matrix(read.csv(paste0("RESULT/1.7-UseFeatures/GenReMat_For", pop2, "Data_low_recombination.csv"), row.names = 1))
   ETA <- list("G1" = list(K = GenReMat_RecRate_high, model = "RKHS"),
               "G2" = list(K = GenReMat_RecRate_middle, model = "RKHS"),
               "G3" = list(K = GenReMat_RecRate_low, model = "RKHS")) # model
}
if ( name.feature == "MAF" ) {
   GenReMat_MAF_high <- as.matrix(read.csv(paste0("RESULT/1.7-UseFeatures/GenReMat_For", pop2, "Data_high_maf.csv"), row.names = 1))
   GenReMat_MAF_middle <- as.matrix(read.csv(paste0("RESULT/1.7-UseFeatures/GenReMat_For", pop2, "Data_middle_maf.csv"), row.names = 1))
   GenReMat_MAF_low <- as.matrix(read.csv(paste0("RESULT/1.7-UseFeatures/GenReMat_For", pop2, "Data_low_maf.csv"), row.names = 1))
   ETA <- list("G1" = list(K = GenReMat_MAF_high, model = "RKHS"),
               "G2" = list(K = GenReMat_MAF_middle, model = "RKHS"),
               "G3" = list(K = GenReMat_MAF_low, model = "RKHS")) # model
}


# ---------------------------------------------------------------------------------------------- #
# -------------------------------------- main calculation -------------------------------------- #
# ---------------------------------------------------------------------------------------------- #
# match IDs
name.sort <- rownames(ETA$G1$K)
for ( i in 1:length(ETA) ) {
   m <- match(name.sort, rownames(ETA[[i]]$K))
   K.sorted <- ETA[[i]]$K[m, m]
   ETA[[i]]$K <- K.sorted
}
m <- match(name.sort, PhenoData$ID)
PhenoData.sort <- PhenoData[m, ]
m <- match(name.sort, PhenoData.ut$ID)
PhenoData.ut.sort <- PhenoData.ut[m, ]
m <- match(name.sort, rownames(CvFold.mat))
CvFold.mat.sort <- CvFold.mat[m, ]
all(PhenoData.sort$ID == PhenoData.ut.sort$ID) # ok
all(PhenoData.sort$ID == rownames(CvFold.mat.sort)) #ok
all(PhenoData.sort$ID == rownames(ETA$G1$K)) # ok

# cross validation
file.out <- paste0(dir.save, "/pred_", trait, "_", pop, "_GBLUP_Use_", name.feature, "_cv.csv") # filename to save
pheno.vec.sort <- PhenoData.sort[[trait]]
pheno.vec.ut.sort <- PhenoData.ut.sort[[trait]]
yPred.mat <- MyFun.CrossValid(pheno.vec = pheno.vec.sort,
                              ETA = ETA, 
                              CvFold.mat = CvFold.mat.sort, 
                              n.rep = n.rep, 
                              n.fold = n.fold) # run CV
yPred.mat.InvBoxCox <- apply(yPred.mat, 2, FUN = MyFun.InvBoxCox) # back-transform
df.yPred <- data.frame("ID" = rownames(yPred.mat), "obs" = pheno.vec.ut.sort, yPred.mat.InvBoxCox) # data to save
m <- match(PhenoData$ID, df.yPred$ID)
df.yPred.ord <- df.yPred[m, ] # re-order
write.csv(df.yPred.ord, file.out, row.names = F) # save
