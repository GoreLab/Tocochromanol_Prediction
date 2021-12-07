# multi-kernel regression

#
library(data.table)
library(gdata)
library(rrBLUP)
library(BGLR)

# param
nIter <-  60000
burnIn <- 40000
thin <- 20

# create dir to save result
folder.save <- paste0("RESULT/8.4-MultiKernel_Across_trans_eval_each_kern_WithDiffWinSize/")
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

# large effect genes (PVE > 5%)
NamLargeGenes <- read.xls("MultiTraitModel/NamLargeGenes.xlsx")

# Di's pathway gene list
df.pathway.Di <- read.xls("RAWDATA/tocochromanol_all_candidate_genes_combined_DW_20190624.xlsx")

# ---------------------------------------------------------------------------- #
# make a data frame
df.all <- NULL
for ( i in 1:nrow(NamLargeGenes) ) {
   tr <- NamLargeGenes$Trait[i]
   gene.name <- NamLargeGenes$GeneName[i]
   id <- TocoGene[TocoGene$Gene == gene.name, "qtl.id"]
   df.i <- GenData[(GenData$Trait == tr) & (GenData$QTL.ID == id), ]
   df.all <- rbind(df.all, df.i)
}
df.all$Gene <- NamLargeGenes$GeneName
df.all$Gene.ID <- NamLargeGenes$Gene
m <- match(df.all$Gene.ID, df.pathway.Di$RefGen_v4.Gene.ID)
df.all$Start.v4 <- as.numeric(gsub(",", "", df.pathway.Di[m, "start_v4"]))
df.all$End.v4 <- as.numeric(gsub(",", "", df.pathway.Di[m, "end_v4"]))
df.all <- df.all[, c("QTL.ID", "Trait", "Gene", "Gene.ID",
                     "Chr", "Start.v4", "End.v4", "PVE",
                     "Peak.Pos.v4", "SupInt.L.Pos.v4", "SupInt.R.Pos.v4",
                     "Common.SupInt.L.Pos.v4", "Common.SupInt.R.Pos.v4",
                     "Peak.Pos.v2", "SupInt.L.Pos.v2", "SupInt.R.Pos.v2",
                     "Common.SupInt.L.Pos.v2", "Common.SupInt.R.Pos.v2")]
for ( i in 1:nrow(df.all) ) {
   df.i <- df.all[i, ]
   png(filename = paste0(folder.save, "/FigA", i, "-WindowSize.png"), 
       width = 700, height = 500)
   GENE <- df.i$Gene
   TRAIT <- df.i$Trait
   CSI.L <- df.i$Common.SupInt.L.Pos.v4 / 1000000
   CSI.R <- df.i$Common.SupInt.R.Pos.v4 / 1000000
   SI.L <- df.i$SupInt.L.Pos.v4 / 1000000
   SI.R <- df.i$SupInt.R.Pos.v4 / 1000000
   GENE.L <- df.i$Start.v4 / 1000000
   GENE.R <- df.i$End.v4 / 1000000
   QTL.peak <- df.i$Peak.Pos.v4 / 1000000
   R <- range(c(CSI.L, CSI.R, SI.L, SI.R, 
                GENE.L, GENE.R, QTL.peak, 
                QTL.peak - 1, QTL.peak + 1))
   plot(x = R, y = c(-1, 1), type = "n", bty = "n", yaxt = "n",
        xlab = "Position (Mbp)", ylab = "", main = paste0(GENE, " for ", TRAIT))
   polygon(x = c(GENE.L, GENE.R, GENE.R, GENE.L),
           y = c(0.75, 0.75, 1.00, 1.00), 
           col = "orange", border = "red", lwd = 2)
   polygon(x = c(SI.L, SI.R, SI.R, SI.L),
           y = c(0.25, 0.25, 0.50, 0.50), 
           col = "gray60", lwd = 2)
   polygon(x = c(CSI.L, CSI.R, CSI.R, CSI.L),
           y = c(-0.25, -0.25, 0.00, 0.00), 
           col = "gray40", lwd = 2)
   abline(v = QTL.peak, lty = 2)
   points(x = QTL.peak, y = -0.50, cex = 1.5, pch = 24, bg = "green", lwd = 2)
   arrows(x0 = QTL.peak - 0.25, x1 = QTL.peak + 0.25,
          y0 = -0.65, y1 = -0.65, length = 0.05, angle = 90, code = 3)
   arrows(x0 = QTL.peak - 1, x1 = QTL.peak + 1,
          y0 = -0.80, y1 = -0.80, length = 0.05, angle = 90, code = 3)
   dev.off()
}






################################################################################

















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
myFun <- function(x){ strsplit(x, "\\(")[[1]][1] }
names(SnpNumList) <- sapply(X = names(SnpNumList), FUN = myFun, USE.NAMES = F)

# make kernel -> Model object
ETA <- list()
for (i in 1:length(SnpNumList)) {
   GenReMat.i <- A.mat(X = geno.mat[, SnpNumList[[i]]] - 1, min.MAF = NULL, shrink = FALSE)
   ETA[[i]] <- list("K" = GenReMat.i, "model" = "RKHS")
}
names(ETA) <- names(SnpNumList)
gc(); gc()

# ---------------------------------------------------------------------------- #
# ----- loop for all QTLs ---------------------------------------------------- #
# ---------------------------------------------------------------------------- #
df.each.use.Ames <- df.each.use.282 <- data.frame("GBS.ID" = rownames(geno.mat))
for ( i.qtl in 1:length(SnpNumList) ) {
   # numbers to indicate SNPs for the i-th QTL
   name.qtl.i <- names(SnpNumList)[i.qtl]
   
   # make kernel -> Model object
   ETA.i <- ETA[-i.qtl]
   
   # -----  1. Ames -> Gb prediction
   # build model
   y <- setNames(object = rep(NA, times = nrow(geno.mat)), nm = rownames(geno.mat))
   y[AmesPheno$ID] <- AmesPheno[[trait]]
   fm <- BGLR(y = y, ETA = ETA, burnIn = burnIn, nIter = nIter, thin = thin, verbose = F,
              saveAt = paste0(folder.save, "/LOGFILE/fm_", trait, "_", name.qtl.i, "_UseAmes_"))
   df.each.use.Ames[[name.qtl.i]] <- fm$yHat
   
   # -----  2. Gb -> Ames prediction
   y <- setNames(object = rep(NA, times = nrow(geno.mat)), nm = rownames(geno.mat))
   y[GbPheno$ID] <- GbPheno[[trait]]
   fm <- BGLR(y = y, ETA = ETA, burnIn = burnIn, nIter = nIter, thin = thin, verbose = F, 
              saveAt = paste0(folder.save, "/LOGFILE/fm_", trait, "_", name.qtl.i, "_UseGb_"))
   df.each.use.282[[name.qtl.i]] <- fm$yHat
   
   print(paste0(name.qtl.i, ": ", i, "/", length(SnpNumList)))
}
# 
# # write out the result
# write.csv(df.each.use.Ames, 
#           file = paste0(folder.save, "/PredRes_AmesToGb_RmEachKernel_", name.file, "_", trait, "_trans.csv"),
#           row.names = F)
# write.csv(df.each.use.282, 
#           file = paste0(folder.save, "/PredRes_GbToAmes_RmEachKernel_", name.file, "_", trait, "_trans.csv"),
#           row.names = F)
# 

