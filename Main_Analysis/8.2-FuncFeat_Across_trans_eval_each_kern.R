# multi-kernel regression

#
library(data.table)
library(gdata)
library(rrBLUP)
library(BGLR)

# param
args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]
nIter <- 60000
burnIn <- 40000
thin <- 20

#####################
trait <- "a.T"
nIter <- 600
burnIn <- 400
thin <- 5
#####################

# create dir to save result
folder.save <- paste0("RESULT/8.2-FuncFeat_Across_trans_eval_each_kern/")
dir.create(folder.save, recursive = T)
dir.create(paste0(folder.save, "/LOGFILE"), recursive = T)

# load geotype & map (not used, but its names are used in this code)
indv.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.indv", header = FALSE)
N <- nrow(indv.tmp)
AccNames <- indv.tmp$V1

# load phenotye 
AmesPheno <- read.csv("RESULT/1.1-MakeAmesPhenoData/AmesPheno_trans.csv", stringsAsFactors = F)
GbPheno <- read.csv("RESULT/1.2-MakeGbPhenoData/GbPheno_trans.csv", stringsAsFactors = F)

# loop for all functional features
name.feature.all <- c("Prox_100kb", "Prox_10kb", "Prox_1kb", "GERP", "RecRate", "MAF")
for ( k in 1:length(name.feature.all) ) {
   # load relationship matrix
   name.feature <- name.feature.all[k]
   
   # 0. model setup
   if ( name.feature == "Prox_100kb" ) {
      GenReMat_prox <- as.matrix(read.csv("RESULT/1.7-UseFeatures/GenReMat_proximity_100kb.csv", row.names = 1))
      GenReMat_non <- as.matrix(read.csv("RESULT/1.7-UseFeatures/GenReMat_non_proximity_100kb.csv", row.names = 1))
      ETA <- list("prox_100kb" = list(K = GenReMat_prox, model = "RKHS"),
                  "dist_100kb" = list(K = GenReMat_non, model = "RKHS"))
   }
   if ( name.feature == "Prox_10kb" ) {
      GenReMat_prox <- as.matrix(read.csv("RESULT/1.7-UseFeatures/GenReMat_proximity_10kb.csv", row.names = 1))
      GenReMat_non <- as.matrix(read.csv("RESULT/1.7-UseFeatures/GenReMat_non_proximity_10kb.csv", row.names = 1))
      ETA <- list("prox_10kb" = list(K = GenReMat_prox, model = "RKHS"),
                  "dist_10kb" = list(K = GenReMat_non, model = "RKHS")) # model
   }
   if ( name.feature == "Prox_1kb" ) {
      GenReMat_prox <- as.matrix(read.csv("RESULT/1.7-UseFeatures/GenReMat_proximity_1kb.csv", row.names = 1))
      GenReMat_non <- as.matrix(read.csv("RESULT/1.7-UseFeatures/GenReMat_non_proximity_1kb.csv", row.names = 1))
      ETA <- list("prox_1kb" = list(K = GenReMat_prox, model = "RKHS"),
                  "dist_1kb" = list(K = GenReMat_non, model = "RKHS")) # model
   }
   if ( name.feature == "GERP" ) {
      GenReMat_GERP_posi <- as.matrix(read.csv("RESULT/1.7-UseFeatures/GenReMat_positive_GERP.csv", row.names = 1))
      GenReMat_GERP_nega <- as.matrix(read.csv("RESULT/1.7-UseFeatures/GenReMat_nagative_GERP.csv", row.names = 1))
      ETA <- list("GERP_posi" = list(K = GenReMat_GERP_posi, model = "RKHS"),
                  "GERP_nega" = list(K = GenReMat_GERP_nega, model = "RKHS")) # model
   }
   if ( name.feature == "RecRate" ) {
      GenReMat_RecRate_high <- as.matrix(read.csv("RESULT/1.7-UseFeatures/GenReMat_high_recombination.csv", row.names = 1))
      GenReMat_RecRate_middle <- as.matrix(read.csv("RESULT/1.7-UseFeatures/GenReMat_middle_recombination.csv", row.names = 1))
      GenReMat_RecRate_low <- as.matrix(read.csv("RESULT/1.7-UseFeatures/GenReMat_low_recombination.csv", row.names = 1))
      ETA <- list("RecRate_high" = list(K = GenReMat_RecRate_high, model = "RKHS"),
                  "RecRate_mid" = list(K = GenReMat_RecRate_middle, model = "RKHS"),
                  "RecRate_low" = list(K = GenReMat_RecRate_low, model = "RKHS")) # model
   }
   if ( name.feature == "MAF" ) {
      GenReMat_MAF_high <- as.matrix(read.csv("RESULT/1.7-UseFeatures/GenReMat_high_maf.csv", row.names = 1))
      GenReMat_MAF_middle <- as.matrix(read.csv("RESULT/1.7-UseFeatures/GenReMat_middle_maf.csv", row.names = 1))
      GenReMat_MAF_low <- as.matrix(read.csv("RESULT/1.7-UseFeatures/GenReMat_low_maf.csv", row.names = 1))
      ETA <- list("MAF_high" = list(K = GenReMat_MAF_high, model = "RKHS"),
                  "MAF_mid" = list(K = GenReMat_MAF_middle, model = "RKHS"),
                  "MAF_low" = list(K = GenReMat_MAF_low, model = "RKHS")) # model
   }
   
   # 1. Ames -> Gb prediction
   df.each <- data.frame("GBS.ID" = AccNames)
   for ( i in 1:length(ETA) ) {
      # build i-th model
      ETA.i <- ETA[i]
      y <- setNames(object = rep(NA, times = N), nm = AccNames)
      y[AmesPheno$ID] <- AmesPheno[[trait]]
      fm <- BGLR(y = y, ETA = ETA.i, burnIn = burnIn, nIter = nIter, thin = thin, verbose = F,
                 saveAt = paste0(folder.save, "/LOGFILE/fm_", trait, "_UseAmes_", names(ETA.i), "_"))
      # model
      saveRDS(object = fm, file = paste0(folder.save, "/LOGFILE/fm_", trait, "_UseAmes_", names(ETA.i), ".Rdata"))
      
      # prediction
      df.each[[names(ETA.i)]] <- fm$yHat
      gc(); gc()
   }
   f <- paste0(folder.save, "/PredRes_AmesToGb_EachKernel_", name.feature, "_", trait, "_trans.csv")
   write.csv(df.each, file = f, row.names = F)
   
   # 2. Gb -> Ames prediction
   df.each <- data.frame("GBS.ID" = AccNames, "obs" = y)
   for ( i in 1:length(ETA) ) {
      # build i-th model
      ETA.i <- ETA[i]
      y <- setNames(object = rep(NA, times = N), nm = AccNames)
      y[GbPheno$ID] <- GbPheno[[trait]]
      fm <- BGLR(y = y, ETA = ETA.i, burnIn = burnIn, nIter = nIter, thin = thin, verbose = F,
                 saveAt = paste0(folder.save, "/LOGFILE/fm_", trait, "_UseGb_", names(ETA.i), "_"))
      # model
      saveRDS(object = fm, file = paste0(folder.save, "/LOGFILE/fm_", trait, "_UseGb_", names(ETA.i), ".Rdata"))
      
      # prediction
      df.each[[names(ETA.i)]] <- fm$yHat
      gc(); gc()
   }
   f <- paste0(folder.save, "/PredRes_GbToAmes_EachKernel_", name.feature, "_", trait, "_trans.csv")
   write.csv(df.each, file = f, row.names = F)
   
   # clean-up again
   gc(); gc()
   
   print(name.feature)
}


