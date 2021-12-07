# make genomic relationship matrix

# library
library(data.table)
library(rrBLUP)

# load geotype & map
indv.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.indv", header = FALSE)
geno.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012")
geno.mat <- as.matrix(geno.tmp[, -1]) # 1st col is ID.
rownames(geno.mat) <- indv.tmp$V1
rm(indv.tmp); rm(geno.tmp); gc(); gc() # remove non-essential object

# load map with features
map <- fread("RESULT/1.7-UseFeatures/map_with_features.csv", data.table = F)

# myfunction -- batch calcualtion of A.mat
myFun.calc.G <- function(X.012, batch = 3000, ...) {
   N <- nrow(X.012)
   P <- ncol(X.012)
   n <- ceiling(P / batch)
   c.all <- rep(0, n)
   WWt.all <- matrix(0, nr = N, nc = N)
   for (i in 1:n) {
      b <- Sys.time()
      print(paste0("start calculation for mini-batch ", i , " of ", n, " batches."))
      start <- (i - 1) * batch + 1
      end <- min(i * batch,  P)
      X.mini <- X.012[, start:end]
      vec.tmp <- apply(X.mini, 2, sum)
      allele.freq.vec <- vec.tmp / (2 * N)
      allele.freq.mat <- t(allele.freq.vec) %x% rep(1, N)
      c.mini <- 2 * sum(allele.freq.vec * (1 - allele.freq.vec))
      W <- X.mini - 2 * allele.freq.mat
      WWt.mini <- tcrossprod(W)
      a <- Sys.time()
      c.all[i] <- c.mini
      WWt.all <- WWt.all + WWt.mini
      print(a-b)
   }
   # G matrix
   G <- WWt.all / sum(c.all)
   return(G)
}


#######################################################################################################
####### For the expression data ######################################################################
#######################################################################################################
# For the accessions with expression data
ExpData <- fread("RAWDATA/ExprData_Ames/expression_BLUE_final_v1.1_B73.csv")
ExpId <- gsub("_", "", ExpData$Accession_ID)
ExpId.common <- intersect(rownames(geno.mat), ExpId)

# (A) GERP score
tf <- map$GERP == 0; table(tf)
G.GERP.nega <- myFun.calc.G(X.012 = geno.mat[ExpId.common, tf], batch = 1000)
G.GERP.posi <- myFun.calc.G(X.012 = geno.mat[ExpId.common, !tf], batch = 1000)
write.csv(G.GERP.nega, "RESULT/1.7-UseFeatures/GenReMat_ForExpData_nagative_GERP.csv")
write.csv(G.GERP.posi, "RESULT/1.7-UseFeatures/GenReMat_ForExpData_positive_GERP.csv")
rm(G.GERP.nega); rm(G.GERP.posi)
gc(); gc()

# (B) MNASE SHOOT
tf <- map$MNASE_SHOOT; table(tf)
G.MNASE.SHOOT.hotspot <- myFun.calc.G(X.012 = geno.mat[ExpId.common, tf], batch = 1000)
G.MNASE.SHOOT.non <- myFun.calc.G(X.012 = geno.mat[ExpId.common, !tf], batch = 1000)
write.csv(G.MNASE.SHOOT.hotspot, "RESULT/1.7-UseFeatures/GenReMat_ForExpData_MNASE_SHOOT_hotspot.csv")
write.csv(G.MNASE.SHOOT.non, "RESULT/1.7-UseFeatures/GenReMat_MNASE_ForExpData_SHOOT_non_hotspot.csv")
rm(G.MNASE.SHOOT.hotspot); rm(G.MNASE.SHOOT.non)
gc(); gc()

# (C) MNASE ROOT
tf <- map$MNASE_ROOT; table(tf)
G.MNASE.ROOT.hotspot <- myFun.calc.G(X.012 = geno.mat[ExpId.common, tf], batch = 1000)
G.MNASE.ROOT.non <- myFun.calc.G(X.012 = geno.mat[ExpId.common, !tf], batch = 1000)
write.csv(G.MNASE.ROOT.hotspot, "RESULT/1.7-UseFeatures/GenReMat_ForExpData_MNASE_ROOT_hotspot.csv")
write.csv(G.MNASE.ROOT.non, "RESULT/1.7-UseFeatures/GenReMat_ForExpData_MNASE_ROOT_non_hotspot.csv")
rm(G.MNASE.ROOT.hotspot); rm(G.MNASE.ROOT.non)
gc(); gc()

# (D) recombination
tf.low <- map$REC.RATE * 1000000 * 100 <= 0.45; table(tf.low)
tf.high <- map$REC.RATE * 1000000 * 100 > 1.65; table(tf.high)
tf.middle <- !(tf.low | tf.high); table(tf.middle)
G.low.rec <- myFun.calc.G(X.012 = geno.mat[ExpId.common, tf.low], batch = 1000)
G.middle.rec <- myFun.calc.G(X.012 = geno.mat[ExpId.common, tf.middle], batch = 1000)
G.high.rec <- myFun.calc.G(X.012 = geno.mat[ExpId.common, tf.high], batch = 1000)
write.csv(G.low.rec, "RESULT/1.7-UseFeatures/GenReMat_ForExpData_low_recombination.csv")
write.csv(G.middle.rec, "RESULT/1.7-UseFeatures/GenReMat_ForExpData_middle_recombination.csv")
write.csv(G.high.rec, "RESULT/1.7-UseFeatures/GenReMat_ForExpData_high_recombination.csv")
rm(G.low.rec); rm(G.middle.rec); rm(G.high.rec)
gc(); gc()

# (E) proximity 1Kb
tf <- map$PROX.1Kb; table(tf)
G.prox.1kb <- myFun.calc.G(X.012 = geno.mat[ExpId.common, tf], batch = 1000)
G.non.prox.1kb <- myFun.calc.G(X.012 = geno.mat[ExpId.common, !tf], batch = 1000)
write.csv(G.prox.1kb, "RESULT/1.7-UseFeatures/GenReMat_ForExpData_proximity_1kb.csv")
write.csv(G.non.prox.1kb, "RESULT/1.7-UseFeatures/GenReMat_ForExpData_non_proximity_1kb.csv")
rm(G.prox.1kb); rm(G.non.prox.1kb)
gc(); gc()

# (F) proximity 10Kb
tf <- map$PROX.10Kb; table(tf)
G.prox.10kb <- myFun.calc.G(X.012 = geno.mat[ExpId.common, tf], batch = 1000)
G.non.prox.10kb <- myFun.calc.G(X.012 = geno.mat[ExpId.common, !tf], batch = 1000)
write.csv(G.prox.10kb, "RESULT/1.7-UseFeatures/GenReMat_ForExpData_proximity_10kb.csv")
write.csv(G.non.prox.10kb, "RESULT/1.7-UseFeatures/GenReMat_ForExpData_non_proximity_10kb.csv")
rm(G.prox.10kb); rm(G.non.prox.10kb)
gc(); gc()

# (G) proximity 100Kb
tf <- map$PROX.100Kb; table(tf)
G.prox.100kb <- myFun.calc.G(X.012 = geno.mat[ExpId.common, tf], batch = 1000)
G.non.prox.100kb <- myFun.calc.G(X.012 = geno.mat[ExpId.common, !tf], batch = 1000)
write.csv(G.prox.100kb, "RESULT/1.7-UseFeatures/GenReMat_ForExpData_proximity_100kb.csv")
write.csv(G.non.prox.100kb, "RESULT/1.7-UseFeatures/GenReMat_ForExpData_non_proximity_100kb.csv")
rm(G.prox.100kb); rm(G.non.prox.100kb)
gc(); gc()


# (H) MAF
tf.low <- map$MAF <= 0.01; table(tf.low)
tf.high <- map$MAF > 0.05; table(tf.high)
tf.middle <- !(tf.low | tf.high); table(tf.middle)
G.low.maf <- myFun.calc.G(X.012 = geno.mat[ExpId.common, tf.low], batch = 1000)
G.middle.maf <- myFun.calc.G(X.012 = geno.mat[ExpId.common, tf.middle], batch = 1000)
G.high.maf <- myFun.calc.G(X.012 = geno.mat[ExpId.common, tf.high], batch = 1000)
write.csv(G.low.maf, "RESULT/1.7-UseFeatures/GenReMat_ForExpData_low_maf.csv")
write.csv(G.middle.maf, "RESULT/1.7-UseFeatures/GenReMat_ForExpData_middle_maf.csv")
write.csv(G.high.maf, "RESULT/1.7-UseFeatures/GenReMat_ForExpData_high_maf.csv")
rm(G.low.maf); rm(G.middle.maf); rm(G.high.maf)
gc(); gc()



#######################################################################################################
####### For the Ames data #############################################################################
#######################################################################################################
# For the Ames accessions
AmesData <- fread("RESULT/1.1-MakeAmesPhenoData/AmesPheno_trans.csv")
AmesId <- AmesData$ID
AmesId.common <- intersect(rownames(geno.mat), AmesId)

# (A) GERP score
tf <- map$GERP == 0; table(tf)
G.GERP.nega <- myFun.calc.G(X.012 = geno.mat[AmesId.common, tf], batch = 1000)
G.GERP.posi <- myFun.calc.G(X.012 = geno.mat[AmesId.common, !tf], batch = 1000)
write.csv(G.GERP.nega, "RESULT/1.7-UseFeatures/GenReMat_ForAmesData_nagative_GERP.csv")
write.csv(G.GERP.posi, "RESULT/1.7-UseFeatures/GenReMat_ForAmesData_positive_GERP.csv")
rm(G.GERP.nega); rm(G.GERP.posi)
gc(); gc()

# (B) MNASE SHOOT
tf <- map$MNASE_SHOOT; table(tf)
G.MNASE.SHOOT.hotspot <- myFun.calc.G(X.012 = geno.mat[AmesId.common, tf], batch = 1000)
G.MNASE.SHOOT.non <- myFun.calc.G(X.012 = geno.mat[AmesId.common, !tf], batch = 1000)
write.csv(G.MNASE.SHOOT.hotspot, "RESULT/1.7-UseFeatures/GenReMat_ForAmesData_MNASE_SHOOT_hotspot.csv")
write.csv(G.MNASE.SHOOT.non, "RESULT/1.7-UseFeatures/GenReMat_MNASE_ForAmesData_SHOOT_non_hotspot.csv")
rm(G.MNASE.SHOOT.hotspot); rm(G.MNASE.SHOOT.non)
gc(); gc()

# (C) MNASE ROOT
tf <- map$MNASE_ROOT; table(tf)
G.MNASE.ROOT.hotspot <- myFun.calc.G(X.012 = geno.mat[AmesId.common, tf], batch = 1000)
G.MNASE.ROOT.non <- myFun.calc.G(X.012 = geno.mat[AmesId.common, !tf], batch = 1000)
write.csv(G.MNASE.ROOT.hotspot, "RESULT/1.7-UseFeatures/GenReMat_ForAmesData_MNASE_ROOT_hotspot.csv")
write.csv(G.MNASE.ROOT.non, "RESULT/1.7-UseFeatures/GenReMat_ForAmesData_MNASE_ROOT_non_hotspot.csv")
rm(G.MNASE.ROOT.hotspot); rm(G.MNASE.ROOT.non)
gc(); gc()

# (D) recombination
tf.low <- map$REC.RATE * 1000000 * 100 <= 0.45; table(tf.low)
tf.high <- map$REC.RATE * 1000000 * 100 > 1.65; table(tf.high)
tf.middle <- !(tf.low | tf.high); table(tf.middle)
G.low.rec <- myFun.calc.G(X.012 = geno.mat[AmesId.common, tf.low], batch = 1000)
G.middle.rec <- myFun.calc.G(X.012 = geno.mat[AmesId.common, tf.middle], batch = 1000)
G.high.rec <- myFun.calc.G(X.012 = geno.mat[AmesId.common, tf.high], batch = 1000)
write.csv(G.low.rec, "RESULT/1.7-UseFeatures/GenReMat_ForAmesData_low_recombination.csv")
write.csv(G.middle.rec, "RESULT/1.7-UseFeatures/GenReMat_ForAmesData_middle_recombination.csv")
write.csv(G.high.rec, "RESULT/1.7-UseFeatures/GenReMat_ForAmesData_high_recombination.csv")
rm(G.low.rec); rm(G.middle.rec); rm(G.high.rec)
gc(); gc()

# (E) proximity 1Kb
tf <- map$PROX.1Kb; table(tf)
G.prox.1kb <- myFun.calc.G(X.012 = geno.mat[AmesId.common, tf], batch = 1000)
G.non.prox.1kb <- myFun.calc.G(X.012 = geno.mat[AmesId.common, !tf], batch = 1000)
write.csv(G.prox.1kb, "RESULT/1.7-UseFeatures/GenReMat_ForAmesData_proximity_1kb.csv")
write.csv(G.non.prox.1kb, "RESULT/1.7-UseFeatures/GenReMat_ForAmesData_non_proximity_1kb.csv")
rm(G.prox.1kb); rm(G.non.prox.1kb)
gc(); gc()

# (F) proximity 10Kb
tf <- map$PROX.10Kb; table(tf)
G.prox.10kb <- myFun.calc.G(X.012 = geno.mat[AmesId.common, tf], batch = 1000)
G.non.prox.10kb <- myFun.calc.G(X.012 = geno.mat[AmesId.common, !tf], batch = 1000)
write.csv(G.prox.10kb, "RESULT/1.7-UseFeatures/GenReMat_ForAmesData_proximity_10kb.csv")
write.csv(G.non.prox.10kb, "RESULT/1.7-UseFeatures/GenReMat_ForAmesData_non_proximity_10kb.csv")
rm(G.prox.10kb); rm(G.non.prox.10kb)
gc(); gc()

# (G) proximity 100Kb
tf <- map$PROX.100Kb; table(tf)
G.prox.100kb <- myFun.calc.G(X.012 = geno.mat[AmesId.common, tf], batch = 1000)
G.non.prox.100kb <- myFun.calc.G(X.012 = geno.mat[AmesId.common, !tf], batch = 1000)
write.csv(G.prox.100kb, "RESULT/1.7-UseFeatures/GenReMat_ForAmesData_proximity_100kb.csv")
write.csv(G.non.prox.100kb, "RESULT/1.7-UseFeatures/GenReMat_ForAmesData_non_proximity_100kb.csv")
rm(G.prox.100kb); rm(G.non.prox.100kb)
gc(); gc()

# (H) MAF
tf.low <- map$MAF <= 0.01; table(tf.low)
tf.high <- map$MAF > 0.05; table(tf.high)
tf.middle <- !(tf.low | tf.high); table(tf.middle)
G.low.maf <- myFun.calc.G(X.012 = geno.mat[AmesId.common, tf.low], batch = 1000)
G.middle.maf <- myFun.calc.G(X.012 = geno.mat[AmesId.common, tf.middle], batch = 1000)
G.high.maf <- myFun.calc.G(X.012 = geno.mat[AmesId.common, tf.high], batch = 1000)
write.csv(G.low.maf, "RESULT/1.7-UseFeatures/GenReMat_ForAmesData_low_maf.csv")
write.csv(G.middle.maf, "RESULT/1.7-UseFeatures/GenReMat_ForAmesData_middle_maf.csv")
write.csv(G.high.maf, "RESULT/1.7-UseFeatures/GenReMat_ForAmesData_high_maf.csv")
rm(G.low.maf); rm(G.middle.maf); rm(G.high.maf)
gc(); gc()



#######################################################################################################
####### For the 282 data #############################################################################
#######################################################################################################
# For the Ames accessions
GbData <- fread("RESULT/1.2-MakeGbPhenoData/GbPheno_trans.csv")
GbId <- GbData$ID
GbId.common <- intersect(rownames(geno.mat), GbId)

# (A) GERP score
tf <- map$GERP == 0; table(tf)
G.GERP.nega <- myFun.calc.G(X.012 = geno.mat[GbId.common, tf], batch = 1000)
G.GERP.posi <- myFun.calc.G(X.012 = geno.mat[GbId.common, !tf], batch = 1000)
write.csv(G.GERP.nega, "RESULT/1.7-UseFeatures/GenReMat_ForGbData_nagative_GERP.csv")
write.csv(G.GERP.posi, "RESULT/1.7-UseFeatures/GenReMat_ForGbData_positive_GERP.csv")
rm(G.GERP.nega); rm(G.GERP.posi)
gc(); gc()

# (B) MNASE SHOOT
tf <- map$MNASE_SHOOT; table(tf)
G.MNASE.SHOOT.hotspot <- myFun.calc.G(X.012 = geno.mat[GbId.common, tf], batch = 1000)
G.MNASE.SHOOT.non <- myFun.calc.G(X.012 = geno.mat[GbId.common, !tf], batch = 1000)
write.csv(G.MNASE.SHOOT.hotspot, "RESULT/1.7-UseFeatures/GenReMat_ForGbData_MNASE_SHOOT_hotspot.csv")
write.csv(G.MNASE.SHOOT.non, "RESULT/1.7-UseFeatures/GenReMat_MNASE_ForGbData_SHOOT_non_hotspot.csv")
rm(G.MNASE.SHOOT.hotspot); rm(G.MNASE.SHOOT.non)
gc(); gc()

# (C) MNASE ROOT
tf <- map$MNASE_ROOT; table(tf)
G.MNASE.ROOT.hotspot <- myFun.calc.G(X.012 = geno.mat[GbId.common, tf], batch = 1000)
G.MNASE.ROOT.non <- myFun.calc.G(X.012 = geno.mat[GbId.common, !tf], batch = 1000)
write.csv(G.MNASE.ROOT.hotspot, "RESULT/1.7-UseFeatures/GenReMat_ForGbData_MNASE_ROOT_hotspot.csv")
write.csv(G.MNASE.ROOT.non, "RESULT/1.7-UseFeatures/GenReMat_ForGbData_MNASE_ROOT_non_hotspot.csv")
rm(G.MNASE.ROOT.hotspot); rm(G.MNASE.ROOT.non)
gc(); gc()

# (D) recombination
tf.low <- map$REC.RATE * 1000000 * 100 <= 0.45; table(tf.low)
tf.high <- map$REC.RATE * 1000000 * 100 > 1.65; table(tf.high)
tf.middle <- !(tf.low | tf.high); table(tf.middle)
G.low.rec <- myFun.calc.G(X.012 = geno.mat[GbId.common, tf.low], batch = 1000)
G.middle.rec <- myFun.calc.G(X.012 = geno.mat[GbId.common, tf.middle], batch = 1000)
G.high.rec <- myFun.calc.G(X.012 = geno.mat[GbId.common, tf.high], batch = 1000)
write.csv(G.low.rec, "RESULT/1.7-UseFeatures/GenReMat_ForGbData_low_recombination.csv")
write.csv(G.middle.rec, "RESULT/1.7-UseFeatures/GenReMat_ForGbData_middle_recombination.csv")
write.csv(G.high.rec, "RESULT/1.7-UseFeatures/GenReMat_ForGbData_high_recombination.csv")
rm(G.low.rec); rm(G.middle.rec); rm(G.high.rec)
gc(); gc()

# (E) proximity 1Kb
tf <- map$PROX.1Kb; table(tf)
G.prox.1kb <- myFun.calc.G(X.012 = geno.mat[GbId.common, tf], batch = 1000)
G.non.prox.1kb <- myFun.calc.G(X.012 = geno.mat[GbId.common, !tf], batch = 1000)
write.csv(G.prox.1kb, "RESULT/1.7-UseFeatures/GenReMat_ForGbData_proximity_1kb.csv")
write.csv(G.non.prox.1kb, "RESULT/1.7-UseFeatures/GenReMat_ForGbData_non_proximity_1kb.csv")
rm(G.prox.1kb); rm(G.non.prox.1kb)
gc(); gc()

# (F) proximity 10Kb
tf <- map$PROX.10Kb; table(tf)
G.prox.10kb <- myFun.calc.G(X.012 = geno.mat[GbId.common, tf], batch = 1000)
G.non.prox.10kb <- myFun.calc.G(X.012 = geno.mat[GbId.common, !tf], batch = 1000)
write.csv(G.prox.10kb, "RESULT/1.7-UseFeatures/GenReMat_ForGbData_proximity_10kb.csv")
write.csv(G.non.prox.10kb, "RESULT/1.7-UseFeatures/GenReMat_ForGbData_non_proximity_10kb.csv")
rm(G.prox.10kb); rm(G.non.prox.10kb)
gc(); gc()

# (G) proximity 100Kb
tf <- map$PROX.100Kb; table(tf)
G.prox.100kb <- myFun.calc.G(X.012 = geno.mat[GbId.common, tf], batch = 1000)
G.non.prox.100kb <- myFun.calc.G(X.012 = geno.mat[GbId.common, !tf], batch = 1000)
write.csv(G.prox.100kb, "RESULT/1.7-UseFeatures/GenReMat_ForGbData_proximity_100kb.csv")
write.csv(G.non.prox.100kb, "RESULT/1.7-UseFeatures/GenReMat_ForGbData_non_proximity_100kb.csv")
rm(G.prox.100kb); rm(G.non.prox.100kb)
gc(); gc()

# (H) MAF
tf.low <- map$MAF <= 0.01; table(tf.low)
tf.high <- map$MAF > 0.05; table(tf.high)
tf.middle <- !(tf.low | tf.high); table(tf.middle)
G.low.maf <- myFun.calc.G(X.012 = geno.mat[GbId.common, tf.low], batch = 1000)
G.middle.maf <- myFun.calc.G(X.012 = geno.mat[GbId.common, tf.middle], batch = 1000)
G.high.maf <- myFun.calc.G(X.012 = geno.mat[GbId.common, tf.high], batch = 1000)
write.csv(G.low.maf, "RESULT/1.7-UseFeatures/GenReMat_ForGbData_low_maf.csv")
write.csv(G.middle.maf, "RESULT/1.7-UseFeatures/GenReMat_ForGbData_middle_maf.csv")
write.csv(G.high.maf, "RESULT/1.7-UseFeatures/GenReMat_ForGbData_high_maf.csv")
rm(G.low.maf); rm(G.middle.maf); rm(G.high.maf)
gc(); gc()

