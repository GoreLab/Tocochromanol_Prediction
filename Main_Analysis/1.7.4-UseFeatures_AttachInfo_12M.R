# library
library(data.table)
library(ape)

# mkdir
dir.create("RESULT/1.7-UseFeatures/Use_12M", recursive = TRUE)

# ---------------------------------------------------------------------------- #
# ----- 1. feature data from Ramstein et al., 2020 & GFF file
# ---------------------------------------------------------------------------- #
# load features
dat.feature <- readRDS("RAWDATA/FuncFeat/functional_annotations.rds")
gc(); gc()

# retain the necessary ones
dat.feature <- dat.feature[, c("CHROM", "POS", "maf", "gerp_non_negative_score", "r")]
gc(); gc()

# load gene annotation (GFF file)
GFF.v4 <- fread("RAWDATA/FuncFeat/Zea_mays.B73_RefGen_v4.59_anno.csv", data.table = F)
GFF.v4.gene <- GFF.v4[GFF.v4$type == "gene", ]
rm(GFF.v4)
gc(); gc()



# ---------------------------------------------------------------------------- #
# ----- 2. attach feature & save map for each chromosome
# ---------------------------------------------------------------------------- #
for ( chr in 1:10 ) {
   # clean-up
   gc(); gc()
   
   # load map file
   f <- paste0("RAWDATA/GRM_12M_SNP/NumericalGenotypeData/merged_ames_282", chr, "_maf.012.pos")
   map <- fread(f, data.table = F)
   colnames(map) <- c("chr", "pos")
   
   # feature data of the i-th chromosome
   dat.feature.i <- dat.feature[dat.feature$CHROM == chr, ]
   GFF.v4.gene.i <- GFF.v4.gene[GFF.v4.gene$chr == chr, ]
   
   # match
   setdiff(map$pos, dat.feature.i$POS)
   m <- match(map$pos, dat.feature.i$POS)
   sum(is.na(m))
   
   # attach
   map$MAF <- dat.feature.i$maf[m]
   map$GERP <- dat.feature.i$gerp_non_negative_score[m]
   map$REC.RATE <- dat.feature.i$r[m]
   
   # add gene proximity
   win.all <- c(1000, 10000, 100000)
   for (w in 1:length(win.all)) {
      win <- win.all[w]
      tf.in.all <- rep(NA, nrow(map))
      for (i in 1:nrow(map)) {
         pos.i <- map$pos[i]
         tf01 <- (GFF.v4.gene.i$start - win) < pos.i
         tf02 <- pos.i < (GFF.v4.gene.i$end + win)
         tf.vec.in <- tf01 & tf02
         tf.in.all[i] <- sum(tf.vec.in) > 0
      }
      if (win == 1000) { map$PROX.1Kb <- tf.in.all }
      if (win == 10000) { map$PROX.10Kb <- tf.in.all }
      if (win == 100000) { map$PROX.100Kb <- tf.in.all }
      gc(); gc()
   }
   
   # save
   f <- paste0("RESULT/1.7-UseFeatures/Use_12M/map_with_feature_chr", chr, ".csv")
   fwrite(map, file = f)
   
   # print
   print(chr)
}
