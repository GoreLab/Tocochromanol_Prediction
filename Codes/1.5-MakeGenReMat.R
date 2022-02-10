# make genomic relationship matrix

# library
library(data.table)
library(rrBLUP)

# mkdir
dir.create("RESULT/1.5-MakeGenReMat")

# load geotype & map
indv.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.indv", header = FALSE)
geno.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012")
geno.mat <- as.matrix(geno.tmp[, -1]) # 1st col is ID.
rownames(geno.mat) <- indv.tmp$V1
rm(indv.tmp); rm(geno.tmp); gc(); gc() # remove non-essential object

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



# # test
# G1 <- myFun.calc.G(X.012 = geno.mat[, 1:2000], batch = 1000)
# G2 <- A.mat(geno.mat[, 1:2000] - 1)
# identical(round(G1, 8), round(G2, 8))

# (1) For all accessions
G <- myFun.calc.G(X.012 = geno.mat, batch = 1000)
write.csv(G, "RESULT/1.5-MakeGenReMat/GenReMat_all.csv")

# (2) For Ames accessions
AmesPheno <- read.csv("RESULT/1.1-MakeAmesPhenoData/AmesPheno.csv", stringsAsFactors = F)
AmesId <- AmesPheno$ID
G <- myFun.calc.G(X.012 = geno.mat[AmesId, ], batch = 1000)
write.csv(G, "RESULT/1.5-MakeGenReMat/GenReMat_ames.csv")

# (3) For Gb accessions
GbPheno <- read.csv("RESULT/1.2-MakeGbPhenoData/GbPheno.csv", stringsAsFactors = F)
GbId <- GbPheno$ID
G <- myFun.calc.G(X.012 = geno.mat[GbId, ], batch = 1000)
write.csv(G, "RESULT/1.5-MakeGenReMat/GenReMat_gb.csv")

# # (4) For the accessions with expression data
# ExpData <- fread("RAWDATA/expression_BLUE_filtered_v1.csv")
# ExpId <- gsub("_", "", ExpData$Accession_ID)
# ExpId.common <- intersect(rownames(geno.mat), ExpId)
# G <- myFun.calc.G(X.012 = geno.mat[ExpId.common, ], batch = 1000)
# write.csv(G, "RESULT/1.5-MakeGenReMat/GenReMat_ForExpData.csv")

# (4) For the accessions with expression data
ExpData <- fread("RAWDATA/ExprData_Ames/expression_BLUE_final_v1.1_B73.csv")
ExpId <- gsub("_", "", ExpData$Accession_ID)
ExpId.common <- intersect(rownames(geno.mat), ExpId) # all 545 retained
G <- myFun.calc.G(X.012 = geno.mat[ExpId.common, ], batch = 1000)
write.csv(G, "RESULT/1.5-MakeGenReMat/GenReMat_ForExpData.csv")


# ---------------------------------------------------------------------------- #
ExpData <- fread("RAWDATA/ExprData_Ames/expression_BLUE_final_v1.1_B73.csv")
ExpId <- gsub("_", "", ExpData$Accession_ID)
ExpId.common <- intersect(rownames(geno.mat), ExpId) # all 545 retained
X <- geno.mat[ExpId.common, ]
x <- apply(X, 2, sum)
af <- x / (2 * nrow(X))
maf <- rep(NA, length(af))
for ( i in 1:length(af) ) {
   maf[i] <- min(af[i], 1 - af[i])
}
hist(maf)
sum(maf == 0)
sum(maf < 0.01)
sum((0 < maf) & (maf <= 0.01))
sum(maf > 0.01)


# ---------------------------------------------------------------------------- #

AmesPheno <- read.csv("RESULT/1.1-MakeAmesPhenoData/AmesPheno.csv", stringsAsFactors = F)
AmesId <- AmesPheno$ID

