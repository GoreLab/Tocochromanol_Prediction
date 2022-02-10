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
bp.width.qtl <- as.integer(args[2]) # set 0 if you do not use this
bp.width.gene <- as.integer(args[3]) # set 0 if you do not use this
name.file <- args[4]
nIter <-  60000
burnIn <- 40000
thin <- 20

# create dir to save result
folder.save <- paste0("RESULT/5.1-MultiKernel_Across_trans/")
dir.create(folder.save, recursive = T)
dir.create(paste0(folder.save, "/LOGFILE"), recursive = T)

# load geotype & map
indv.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.indv", header = FALSE)
geno.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012")
geno.mat <- as.matrix(geno.tmp[, -1]) # 1st col is ID.
rownames(geno.mat) <- indv.tmp$V1
map.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.pos")
map <- data.frame("chr" = map.tmp$V1, "pos" = map.tmp$V2)
rm(indv.tmp); rm(geno.tmp); rm(map.tmp); gc(); gc() # remove non-essential object

# knwon QTLs
GenData <- read.csv("RAWDATA/NAM_QTL/qtl.data.from.di.csv", stringsAsFactors = FALSE)
colnames(GenData)[1] <- "QTL.ID"

# known genes
TocoGene <- read.xls("RAWDATA/NAM_QTL/TocoGene.xlsx", stringsAsFactors = FALSE)
colnames(TocoGene)[9:ncol(TocoGene)] <- c("a.T", "d.T", "g.T", "Total.Tocopherols",
                                          "a.T3", "d.T3", "g.T3", "Total.Tocotrienols",
                                          "Total.Tocochromanols") # rename

# load phenotye 
AmesPheno <- read.csv("RESULT/1.1-MakeAmesPhenoData/AmesPheno_trans.csv", stringsAsFactors = F)
GbPheno <- read.csv("RESULT/1.2-MakeGbPhenoData/GbPheno_trans.csv", stringsAsFactors = F)


# Get SNPs: choose appropriate method
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
saveRDS(SnpNumList, file = paste0(folder.save, "/SnpNumList_", name.file, "_", trait, ".Rdata"))

# make kernel -> Model object
ETA <- list()
for (i in 1:length(SnpNumList)) {
   GenReMat.i <- A.mat(X = geno.mat[, SnpNumList[[i]]] - 1, min.MAF = NULL, shrink = FALSE)
   ETA[[i]] <- list("K" = GenReMat.i, "model" = "RKHS")
}
names(ETA) <- names(SnpNumList)

# ---------------------------------------------------------------------------------------------------------- #
# ---------------------------------------- 1. Ames -> Gb prediction ---------------------------------------- #
# ---------------------------------------------------------------------------------------------------------- #
# build model
y <- setNames(object = rep(NA, times = nrow(geno.mat)), nm = rownames(geno.mat))
y[AmesPheno$ID] <- AmesPheno[[trait]]
fm <- BGLR(y = y, ETA = ETA, burnIn = burnIn, nIter = nIter, thin = thin, verbose = F,
           saveAt = paste0(folder.save, "/LOGFILE/fm_", trait, "_UseAmes_"))
saveRDS(object = fm, file = paste0(folder.save, "/LOGFILE/fm_", trait, "_UseAmes.Rdata"))

# prediction
df.each <- data.frame("GBS.ID" = names(y), "obs" = y)
for (i in 1:length(fm$ETA)){
   pred.i <- fm$mu + fm$ETA[[i]]$u
   df.each[[names(fm$ETA)[[i]]]] <- pred.i
}
df.each$UseAll <- fm$yHat
write.csv(df.each, 
          file = paste0(folder.save, "/PredRes_AmesToGb_MultiKernel_", name.file, "_", trait, "_trans.csv"),
          row.names = F)

# ---------------------------------------------------------------------------------------------------------- #
# ---------------------------------------- 2. Gb -> Ames prediction ---------------------------------------- #
# ---------------------------------------------------------------------------------------------------------- #
# build model
y <- setNames(object = rep(NA, times = nrow(geno.mat)), nm = rownames(geno.mat))
y[GbPheno$ID] <- GbPheno[[trait]]
fm <- BGLR(y = y, ETA = ETA, burnIn = burnIn, nIter = nIter, thin = thin, verbose = F, 
           saveAt = paste0(folder.save, "/LOGFILE/fm_", trait, "_UseGb_"))
saveRDS(object = fm, file = paste0(folder.save, "/LOGFILE/fm_", trait, "_UseGb.Rdata"))

# prediction
df.each <- data.frame("GBS.ID" = names(y), "obs" = y)
for (i in 1:length(fm$ETA)){
   pred.i <- fm$mu + fm$ETA[[i]]$u
   df.each[[names(fm$ETA)[[i]]]] <- pred.i
}
df.each$UseAll <- fm$yHat
write.csv(df.each, 
          file = paste0(folder.save, "/PredRes_GbToAmes_MultiKernel_", name.file, "_", trait, "_trans.csv"),
          row.names = F)


