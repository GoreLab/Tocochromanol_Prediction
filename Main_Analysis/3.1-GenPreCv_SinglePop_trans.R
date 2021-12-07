# Genomic Prediction, Cross valiation

# Library & source
library(data.table)
library(BGLR)

# variables
args <- commandArgs(trailingOnly=T)
trait <- args[1]
pop <- args[2] # "ames" or "gb"
model <- args[3] # "BRR" or "BayesB" (no 'GBLUP' option!)
nIter <-  12000
burnIn <- 8000
thin <- 5
n.rep <- 10

# ###########################################################
# trait <- "a.T"
# pop <- "ames" # "ames" or "gb"
# model <- "BRR" # "BRR" or "BayesB" (no 'GBLUP' option!)
# nIter <-  60
# burnIn <- 40
# thin <- 5
# n.rep <- 10
# ###########################################################


# -------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------- load data  ---------------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #
# Make Folder to save result
name.dir <- paste0("RESULT/3.1-GenPreCv_SinglePop_trans")
dir.create(name.dir, recursive = TRUE)
name.dir.logfile <- paste0(name.dir, "/Logfile")
dir.create(name.dir.logfile)

# load geotype & map
indv.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.indv", header = FALSE)
geno.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012")
geno.mat <- as.matrix(geno.tmp[, -1]) # 1st col is ID.
rownames(geno.mat) <- indv.tmp$V1
map.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.pos")
map <- data.frame("chr" = map.tmp$V1, "pos" = map.tmp$V2)
rm(indv.tmp); rm(geno.tmp); rm(map.tmp); gc(); gc() # remove non-essential object

# load phenotye 
if ( pop == "ames" ) { PhenoData <- read.csv("RESULT/1.1-MakeAmesPhenoData/AmesPheno_trans.csv", stringsAsFactors = F) }
if ( pop == "gb" ) { PhenoData <- read.csv("RESULT/1.2-MakeGbPhenoData/GbPheno_trans.csv", stringsAsFactors = F) }
pheno <- PhenoData[[trait]]

# load CV set
CvFold <- read.csv(file = paste0("RESULT/1.4-MakeCvFold/CvFold_", pop, ".csv"), row.names = 1)
CvNum <- as.matrix(CvFold[, -1])
rownames(CvNum) <- CvFold$ID


# -------------------------------------------------------------------------------------------------------- #
# ------------------------------ data processing before genomic prediction  ------------------------------ #
# -------------------------------------------------------------------------------------------------------- #
# get genotype data for the training population
m <- match(PhenoData$ID, rownames(geno.mat))
geno <- geno.mat[m ,]

# remove genotype without any polymorphism
tmp <- apply(geno, 2, sum)
tf <- (tmp == 0) | (tmp == (nrow(geno) * 2))
score <- geno[, !tf]

# build model
if ( model == "BayesB" ) { ETA <- list("MRK" = list(X = score, model = "BayesB", saveEffects = FALSE)) }
if ( model == "BRR" ) { ETA <- list("MRK" = list(X = score, model = "BRR", saveEffects = FALSE)) }


# ------------------------------------------------------------------------------- #
# ------------------------------ cross validation  ------------------------------ #
# ------------------------------------------------------------------------------- #
# CV params
n.fold <- length(unique(CvNum[, 1]))

# object to save result
yPred <- matrix(NA, nr = nrow(CvNum), nc = ncol(CvNum), dimnames = dimnames(CvNum))
# Run Cv
for ( r in 1:n.rep ) {
   # Cv
   cv.num <- CvNum[, r]
   
   for ( f in 1: n.fold ) {
      # clock
      b <- Sys.time()
      
      # make input vector
      y <- pheno
      y[cv.num == f] <- NA # mask phenotype
      
      # fitting via BGLR
      fm <- BGLR(y = y, ETA = ETA,
                 nIter = nIter, burnIn = burnIn, thin = thin,
                 saveAt = paste0(name.dir.logfile, "/fm_"), verbose = FALSE)
      saveRDS(fm, file = paste0(name.dir.logfile, "/fm_rep=", r, "_fold=", f, ".Rdata"))

      # keep result
      yPred[cv.num == f, r] <- fm$yHat[cv.num == f]
      
      # show current status
      print(c(r, f))
      
      # clock
      a <- Sys.time()
      print(a-b)
   }
}

# save result
df.save <- data.frame("ID" = rownames(yPred), "obs" = pheno, yPred)
write.csv(df.save, file = paste0("RESULT/3.1-GenPreCv_SinglePop_trans/yPred_Trait=", trait, "_Train=", pop, "_Model=", model, ".csv"))




