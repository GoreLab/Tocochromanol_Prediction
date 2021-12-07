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

# #################
# trait = "d.T"
# nIter = 1200
# burnIn = 400
# thin = 5
# #################

# create dir to save result
folder.save <- paste0("RESULT/6.3-UseFeatures_AcrossPopPred_trans_12M/")
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
      GenReMat_prox <- as.matrix(read.csv("RESULT/1.7-UseFeatures/Use_12M/GenReMat_use_12M_SNP_PROX.100Kb_Prox.csv", row.names = 1))
      GenReMat_non <- as.matrix(read.csv("RESULT/1.7-UseFeatures/Use_12M/GenReMat_use_12M_SNP_PROX.100Kb_Dist.csv", row.names = 1))
      ETA <- list("G1" = list(K = GenReMat_prox, model = "RKHS"),
                  "G2" = list(K = GenReMat_non, model = "RKHS"))
   }
   if ( name.feature == "Prox_10kb" ) {
      GenReMat_prox <- as.matrix(read.csv("RESULT/1.7-UseFeatures/Use_12M/GenReMat_use_12M_SNP_PROX.10Kb_Prox.csv", row.names = 1))
      GenReMat_non <- as.matrix(read.csv("RESULT/1.7-UseFeatures/Use_12M/GenReMat_use_12M_SNP_PROX.10Kb_Dist.csv", row.names = 1))
      ETA <- list("G1" = list(K = GenReMat_prox, model = "RKHS"),
                  "G2" = list(K = GenReMat_non, model = "RKHS")) # model
   }
   if ( name.feature == "Prox_1kb" ) {
      GenReMat_prox <- as.matrix(read.csv("RESULT/1.7-UseFeatures/Use_12M/GenReMat_use_12M_SNP_PROX.1Kb_Prox.csv", row.names = 1))
      GenReMat_non <- as.matrix(read.csv("RESULT/1.7-UseFeatures/Use_12M/GenReMat_use_12M_SNP_PROX.1Kb_Dist.csv", row.names = 1))
      ETA <- list("G1" = list(K = GenReMat_prox, model = "RKHS"),
                  "G2" = list(K = GenReMat_non, model = "RKHS")) # model
   }
   if ( name.feature == "GERP" ) {
      GenReMat_GERP_posi <- as.matrix(read.csv("RESULT/1.7-UseFeatures/Use_12M/GenReMat_use_12M_SNP_GERP_Positive.csv", row.names = 1))
      GenReMat_GERP_nega <- as.matrix(read.csv("RESULT/1.7-UseFeatures/Use_12M/GenReMat_use_12M_SNP_GERP_Negative.csv", row.names = 1))
      ETA <- list("G1" = list(K = GenReMat_GERP_posi, model = "RKHS"),
                  "G2" = list(K = GenReMat_GERP_nega, model = "RKHS")) # model
   }
   if ( name.feature == "RecRate" ) {
      GenReMat_RecRate_high <- as.matrix(read.csv("RESULT/1.7-UseFeatures/Use_12M/GenReMat_use_12M_SNP_REC.RATE_High.csv", row.names = 1))
      GenReMat_RecRate_middle <- as.matrix(read.csv("RESULT/1.7-UseFeatures/Use_12M/GenReMat_use_12M_SNP_REC.RATE_Middle.csv", row.names = 1))
      GenReMat_RecRate_low <- as.matrix(read.csv("RESULT/1.7-UseFeatures/Use_12M/GenReMat_use_12M_SNP_REC.RATE_Low.csv", row.names = 1))
      ETA <- list("G1" = list(K = GenReMat_RecRate_high, model = "RKHS"),
                  "G2" = list(K = GenReMat_RecRate_middle, model = "RKHS"),
                  "G3" = list(K = GenReMat_RecRate_low, model = "RKHS")) # model
   }
   if ( name.feature == "MAF" ) {
      GenReMat_MAF_high <- as.matrix(read.csv("RESULT/1.7-UseFeatures/Use_12M/GenReMat_use_12M_SNP_MAF_High.csv", row.names = 1))
      GenReMat_MAF_middle <- as.matrix(read.csv("RESULT/1.7-UseFeatures/Use_12M/GenReMat_use_12M_SNP_MAF_Middle.csv", row.names = 1))
      ETA <- list("G1" = list(K = GenReMat_MAF_high, model = "RKHS"),
                  "G2" = list(K = GenReMat_MAF_middle, model = "RKHS")) # model
   }
   
   # 1. Ames -> Gb prediction
   # build model
   y <- setNames(object = rep(NA, times = N), nm = AccNames)
   y[AmesPheno$ID] <- AmesPheno[[trait]]
   fm <- BGLR(y = y, ETA = ETA, burnIn = burnIn, nIter = nIter, thin = thin, verbose = F,
              saveAt = paste0(folder.save, "/LOGFILE/fm_", trait, "_UseAmes_", name.feature, "_"))
   saveRDS(object = fm, file = paste0(folder.save, "/LOGFILE/fm_", trait, "_UseAmes_", name.feature, ".Rdata"))
   # prediction
   df.each <- data.frame("GBS.ID" = names(y), "obs" = y)
   for (i in 1:length(fm$ETA)){
      pred.i <- fm$mu + fm$ETA[[i]]$u
      df.each[[names(fm$ETA)[[i]]]] <- pred.i
   }
   df.each$UseAll <- fm$yHat
   write.csv(df.each, 
             file = paste0(folder.save, "/PredRes_AmesToGb_MultiKernel_", name.feature, "_", trait, "_trans.csv"),
             row.names = F)
   gc(); gc()
   
   # 2. Gb -> Ames prediction
   # build model
   y <- setNames(object = rep(NA, times = N), nm = AccNames)
   y[GbPheno$ID] <- GbPheno[[trait]]
   fm <- BGLR(y = y, ETA = ETA, burnIn = burnIn, nIter = nIter, thin = thin, verbose = F, 
              saveAt = paste0(folder.save, "/LOGFILE/fm_", trait, "_UseGb_"))
   saveRDS(object = fm, file = paste0(folder.save, "/LOGFILE/fm_", trait, "_UseGb_", name.feature, ".Rdata"))
   # prediction
   df.each <- data.frame("GBS.ID" = names(y), "obs" = y)
   for (i in 1:length(fm$ETA)){
      pred.i <- fm$mu + fm$ETA[[i]]$u
      df.each[[names(fm$ETA)[[i]]]] <- pred.i
   }
   df.each$UseAll <- fm$yHat
   write.csv(df.each, 
             file = paste0(folder.save, "/PredRes_GbToAmes_MultiKernel_", name.feature, "_", trait, "_trans.csv"),
             row.names = F)
   
   gc(); gc()
   print(name.feature)
}



