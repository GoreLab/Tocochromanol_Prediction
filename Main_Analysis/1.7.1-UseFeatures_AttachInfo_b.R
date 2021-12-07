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

###################################################################
# Use 1Kb region
tf.in.all <- rep(NA, nrow(map))
for (i in 1:nrow(map)) {
   chr.i <- map$chr[i]
   pos.i <- map$pos[i]
   GFF.v4.gene.i <- GFF.v4.gene[GFF.v4.gene$seqid == chr.i, ]
   tf01 <- (GFF.v4.gene.i$start - 1000) < pos.i
   tf02 <- pos.i < (GFF.v4.gene.i$end + 1000)
   tf.vec.in <- tf01 & tf02
   tf.in.all[i] <- sum(tf.vec.in) > 0
}
map$gene.proximity.1Kb <- tf.in.all
gc(); gc()


###################################################################
# Use 10Kb region
tf.in.all <- rep(NA, nrow(map))
for (i in 1:nrow(map)) {
   chr.i <- map$chr[i]
   pos.i <- map$pos[i]
   GFF.v4.gene.i <- GFF.v4.gene[GFF.v4.gene$seqid == chr.i, ]
   tf01 <- (GFF.v4.gene.i$start - 10000) < pos.i
   tf02 <- pos.i < (GFF.v4.gene.i$end + 10000)
   tf.vec.in <- tf01 & tf02
   tf.in.all[i] <- sum(tf.vec.in) > 0
}
map$gene.proximity.10Kb <- tf.in.all
gc(); gc()



###################################################################
# Use 100Kb region
tf.in.all <- rep(NA, nrow(map))
for (i in 1:nrow(map)) {
   chr.i <- map$chr[i]
   pos.i <- map$pos[i]
   GFF.v4.gene.i <- GFF.v4.gene[GFF.v4.gene$seqid == chr.i, ]
   tf01 <- (GFF.v4.gene.i$start - 100000) < pos.i
   tf02 <- pos.i < (GFF.v4.gene.i$end + 100000)
   tf.vec.in <- tf01 & tf02
   tf.in.all[i] <- sum(tf.vec.in) > 0
}
map$gene.proximity.100Kb <- tf.in.all
gc(); gc()


###################################################################
fwrite(map, "RESULT/1.7-UseFeatures/map_b.csv")



