# Evaluate Prediction Accuracy Across Populations

# source
library(data.table)

# params
trait.all <- c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3",
               "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")
train.all <- c("gb", "ames")
model.all <- c("BRR", "BayesB")

# file to save
dir.save <- "RESULT/2.2-PredictionAccuracy_trans"

# mkdir to save
dir.create(dir.save)

# make a data.frame to save result
df.save <- expand.grid(trait.all, train.all, model.all)
colnames(df.save) <- c("trait", "train", "model")
df.save$trait <- as.character(df.save$trait)
df.save$train <- as.character(df.save$train)
df.save$model <- as.character(df.save$model)

# load box-cox parameter
BoxCoxParam.ames <- read.csv("RESULT/1.1-MakeAmesPhenoData/Ames_BoxCoxParam.csv")
BoxCoxParam.gb <- read.csv("RESULT/1.2-MakeGbPhenoData/Gb_BoxCoxParam.csv")

# loop
for ( i in 1:nrow(df.save) ) {
   # params
   trait <- df.save$trait[i]
   train <- df.save$train[i]
   model <- df.save$model[i]
   
   # -------------------------------------------------------------------------------------------------------- #
   # ------------------------------------------ Load Genotype data ------------------------------------------ #
   # -------------------------------------------------------------------------------------------------------- #
   # load geotype & map
   indv.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.indv", header = FALSE)
   geno.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012")
   geno.mat <- as.matrix(geno.tmp[, -1]) # 1st col is ID.
   rownames(geno.mat) <- indv.tmp$V1
   map.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.pos")
   map <- data.frame("chr" = map.tmp$V1, "pos" = map.tmp$V2)
   rm(indv.tmp); rm(geno.tmp); rm(map.tmp); gc(); gc() # remove non-essential object
   
   # ------------------------------------------------------------------------------------------------------ #
   # ------------------------------------------ Load MCMC result ------------------------------------------ #
   # ------------------------------------------------------------------------------------------------------ #
   # dir to load MCMC result
   name.dir <- paste0("RESULT/2.1-GenPreFit_SinglePop_trans/Trait=", trait, "_Train=", train, "_Model=", model)
   
   # load marker effects
   name.file <- paste0(name.dir, "/Logfile/fm_ETA_MRK_b.bin")
   num <- readBin(con = name.file, what = "numeric", n = 2) # number of chains & number of variables
   mrk.mcmc <- readBin(con = name.file, what = "numeric", n = 2 + prod(num))
   mrk.mcmc <- mrk.mcmc[3:length(mrk.mcmc)] # remove first two
   mrk.mcmc.mat <- t(matrix(mrk.mcmc, nr = num[2], nc = num[1]))
   
   # load location param (mu)
   name.file <- paste0(name.dir, "/Logfile/fm_mu.dat")
   mu <- read.table(file = name.file)
   mu <- as.numeric(tail(unlist(mu), num[1])) # the last num[1] is the samples after burn-in
   
   
   # ------------------------------------------------------------------------------------------------ #
   # ------------------------------------------ Load polym ------------------------------------------ #
   # ------------------------------------------------------------------------------------------------ #
   # polym info
   name.file <- paste0(name.dir, "/PolymCheck.csv")
   polym.info <- read.csv(name.file)
    
   # ----------------------------------------------------------------------------------------- #
   # --------------------------------- calc predicted values --------------------------------- #
   # ----------------------------------------------------------------------------------------- #
   # predcition
   geno.mat.new <- geno.mat[, polym.info$Polym] # focus on the test set & polym mrk
   GenValue.mcmc <- tcrossprod(x = geno.mat.new, y = mrk.mcmc.mat)
   PredValue <- apply(GenValue.mcmc, 1, mean) + mean(mu)
   
   
   # ----------------------------------------------------------------------------------------------- #
   # --------------------------------------- inverse box-cox --------------------------------------- #
   # ----------------------------------------------------------------------------------------------- #
   if ( train == "ames" ) {
      lambda <- BoxCoxParam.ames$lambda[BoxCoxParam.ames$trait == trait]
      const <- BoxCoxParam.ames$const[BoxCoxParam.ames$trait == trait]
   } else {
      lambda <- BoxCoxParam.gb$lambda[BoxCoxParam.gb$trait == trait]
      const <- BoxCoxParam.gb$const[BoxCoxParam.gb$trait == trait]
   }
   tmp.vec <- as.numeric(PredValue)
   if (lambda == 0) { tmp.vec.original <- exp(tmp.vec) } else { tmp.vec.original <- tmp.vec ^ (1/lambda) }
   tmp.vec.original <- tmp.vec.original - const
   PredValue <- tmp.vec.original
   
   
   # save
   filename.save <- paste0("RESULT/2.2-PredictionAccuracy_trans/yPred_trait=", trait, "_train=", train, "_model=", model, ".csv")
   df <- data.frame("ID" = rownames(geno.mat.new), "pred" = PredValue)
   write.csv(df, filename.save, row.names = F)
}

