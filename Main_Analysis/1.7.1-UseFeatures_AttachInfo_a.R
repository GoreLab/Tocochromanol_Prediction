# make genomic relationship matrix

# library
library(data.table)
library(ape)

# mkdir
dir.create("RESULT/1.7-UseFeatures")

# load geotype & map
indv.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.indv", header = FALSE)
geno.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012")
geno.mat <- as.matrix(geno.tmp[, -1]) # 1st col is ID.
rownames(geno.mat) <- indv.tmp$V1
map.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.pos")
map <- data.frame("chr" = map.tmp$V1, "pos" = map.tmp$V2)
rm(indv.tmp); rm(geno.tmp); rm(map.tmp); gc(); gc() # remove non-essential object

# load gene annotation (GFF file)
GFF.v4 <- read.gff("RAWDATA/FuncFeat/Zea_mays.AGPv4.40.gff3.gz")
GFF.v4.gene <- GFF.v4[GFF.v4$type == "gene", ]
rm(GFF.v4)

# load features
dat.feature <- readRDS("RAWDATA/FuncFeat/functional_annotations.rds")
head(dat.feature)
gc(); gc()

# ###################################################################
# # check unique matching
# P <- nrow(map); tmp.all <- rep(NA, P); count <- 0
# for ( i in 1:10 ) {
#    # i-th chromosome
#    map.i <- map[map$chr == i, ]
#    dat.feature.i <- dat.feature[dat.feature$CHROM == i, ]
#    dat.feature.i.pos <- dat.feature.i$POS
#    gc(); gc()
# 
#    for ( j in 1:nrow(map.i) ) {
#       a <- Sys.time()
#       # j-th position in i-th chromosome
#       count <- count + 1
#       pos.ij <- map.i$pos[j]
#       tf <- dat.feature.i.pos == pos.ij
#       tmp.all[count] <- sum(tf)
#       b <- Sys.time()
#       print(b-a)
#    }
#    print(i)
# }
# table(tmp.all)
# saveRDS(tmp.all, file = "tmp.rds")

###################################################################
fature.attach <- NULL
for ( i in c(10, 1:9) ) {
   # i-th chromosome
   map.i <- map[map$chr == i, ]
   pos.i <- map.i$pos
   dat.feature.i <- dat.feature[dat.feature$CHROM == i, ]
   dat.feature.i.pos <- dat.feature.i$POS
   m <- match(pos.i, dat.feature.i.pos)
   dat.feature.i.m <- dat.feature.i[m, ]
   if (i == 10) {
      fature.attach <- dat.feature.i.m
   } else {
      fature.attach <- rbind(fature.attach, dat.feature.i.m)
   }
   gc(); gc()
   print(i)
}
all(map$chr == fature.attach$CHROM)
all(map$pos == fature.attach$POS)
tmp <- cbind(map, fature.attach)
fwrite(tmp, "RESULT/1.7-UseFeatures/map_a.csv")
gc(); gc()
# 
# 
# ###################################################################
# # Use 1Kb region
# tf.in.all <- rep(NA, nrow(map))
# for (i in 1:nrow(map)) {
#    b <- Sys.time()
#    chr.i <- map$chr[i]
#    pos.i <- map$pos[i]
#    GFF.v4.gene.i <- GFF.v4.gene[GFF.v4.gene$seqid == chr.i, ]
#    tf01 <- (GFF.v4.gene.i$start - 1000) < pos.i
#    tf02 <- pos.i < (GFF.v4.gene.i$end + 1000)
#    tf.vec.in <- tf01 & tf02
#    tf.in.all[i] <- sum(tf.vec.in) > 0
#    a <- Sys.time()
#    print(i); print(a -b)
# }
# map$gene.proximity.1Kb <- tf.in.all
# gc(); gc()
# 
# 
# ###################################################################
# # Use 10Kb region
# tf.in.all <- rep(NA, nrow(map))
# for (i in 1:nrow(map)) {
#    chr.i <- map$chr[i]
#    pos.i <- map$pos[i]
#    GFF.v4.gene.i <- GFF.v4.gene[GFF.v4.gene$seqid == chr.i, ]
#    tf01 <- (GFF.v4.gene.i$start - 10000) < pos.i
#    tf02 <- pos.i < (GFF.v4.gene.i$end + 10000)
#    tf.vec.in <- tf01 & tf02
#    tf.in.all[i] <- sum(tf.vec.in) > 0
#    gc(); gc()
# }
# map$gene.proximity.10Kb <- tf.in.all
# 
# 
# 
# ###################################################################
# # Use 100Kb region
# tf.in.all <- rep(NA, nrow(map))
# for (i in 1:nrow(map)) {
#    chr.i <- map$chr[i]
#    pos.i <- map$pos[i]
#    GFF.v4.gene.i <- GFF.v4.gene[GFF.v4.gene$seqid == chr.i, ]
#    tf01 <- (GFF.v4.gene.i$start - 100000) < pos.i
#    tf02 <- pos.i < (GFF.v4.gene.i$end + 100000)
#    tf.vec.in <- tf01 & tf02
#    tf.in.all[i] <- sum(tf.vec.in) > 0
#    gc(); gc()
# }
# map$gene.proximity.100Kb <- tf.in.all
# gc(); gc()
# 
# 
# ###################################################################
# fwrite(map, "RESULT/1.7-UseFeatures/map_b.csv")
# 
# 
# ###################################################################
# map.a <- fread("RESULT/1.7-UseFeatures/map_a.csv")
# map.b <- fread("RESULT/1.7-UseFeatures/map_b.csv")
# tmp <- cbind(map.a, map.b)
# df.save <- data.frame("chr" = tmp$chr,
#                       "pos" = tmp$pos,
#                       "MAF" = tmp$maf,
#                       "GERP" = tmp$gerp_non_negative_score,
#                       "MNASE_SHOOT" = tmp$mnase_hotspot_shoots,
#                       "MNASE_ROOT" = tmp$mnase_hotspot_roots,
#                       "REC.RATE" = tmp$r,
#                       "PROX.1Kb" = tmp$gene.proximity.1Kb,
#                       "PROX.10Kb" = tmp$gene.proximity.10Kb,
#                       "PROX.100Kb" = tmp$gene.proximity.100Kb)
# fwrite(df.save, "RESULT/1.7-UseFeatures/map_with_features.csv")
# 
