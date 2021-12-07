# make genomic relationship matrix

# library
library(data.table)
library(rrBLUP)

# mkdir
dir.create("RESULT/7.3-MakeGenReMat_ForGbExprData")

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


# For Gb accessions with expression data
GbPhenoInfo <- read.csv("RESULT/1.2-MakeGbPhenoData/GbInfo.csv")
GbExprInfo <- read.csv("RESULT/7.1-MakeExprData/GbExprDataInfo_unique.csv", stringsAsFactors = F)
m <- match(GbExprInfo$Name.Lipka.Manual, GbPhenoInfo$Name.Lipka)
GbExprInfo$GBS.taxa <- GbPhenoInfo$GBS.taxa[m]
GbId <- GbExprInfo$GBS.taxa
G <- myFun.calc.G(X.012 = geno.mat[GbId, ], batch = 1000)
write.csv(G, "RESULT/7.3-MakeGenReMat_ForGbExprData/GenReMat_gb_expr.csv")


