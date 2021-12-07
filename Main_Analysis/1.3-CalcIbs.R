# source
library(data.table)

# mkdir 
dir.create("RESULT/1.3-CalcIbs")

# load geotype & map
indv.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.indv", header = FALSE)
geno.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012")
geno.mat <- as.matrix(geno.tmp[, -1]) # 1st col is ID.
rownames(geno.mat) <- indv.tmp$V1
map.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.pos")
map <- data.frame("chr" = map.tmp$V1, "pos" = map.tmp$V2)
rm(indv.tmp); rm(geno.tmp); rm(map.tmp); gc(); gc() # remove non-essential object

# calculate pairwise IBS
P <- ncol(geno.mat)
score.mat <- geno.mat - 1
XtX <- tcrossprod(score.mat)
IbsMat <- (XtX + P) / (2 * P)

# save the matrix
write.csv(IbsMat, "RESULT/1.3-CalcIbs/IbsMat.csv")
