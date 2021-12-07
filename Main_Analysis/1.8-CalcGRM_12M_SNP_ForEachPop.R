# packages
library(data.table)

# mkdir
dir.save <- "RESULT/1.8-GRM_12M_EachPop"
dir.create(dir.save, recursive = T)

# objects
n.chr <- 10

# function to calculate important info for GRM
myFun.calc.GRM.info <- function(X.012, batch = 3000, ...) {
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
   # return list
   ResList <- list("Numerator" = WWt.all, "Denominator" = sum(c.all))
   return(ResList)
}

# just in case -- count the number of SNPs in each chromosome
df.n.snp <- data.frame("chr" = 1:n.chr, "n.SNP" = NA)
for ( i in 1:n.chr ) {
   f <- paste0("RAWDATA/GRM_12M_SNP/NumericalGenotypeData/merged_ames_282", i, "_maf.012.pos")
   df.pos <- fread(f, data.table = F)
   n.snp <- nrow(df.pos)
   df.n.snp$n.SNP[i] <- n.snp
}

# get id
df.id <- read.table("RAWDATA/GRM_12M_SNP/NumericalGenotypeData/merged_ames_2821_maf.012.indv")
id.vec <- df.id$V1

# phenotyupe data to get ID
AmesPheno <- read.csv("RESULT/1.1-MakeAmesPhenoData/AmesPheno.csv", stringsAsFactors = F)
AmesId <- AmesPheno$ID
m.ames <- match(AmesId, id.vec)
GbPheno <- read.csv("RESULT/1.2-MakeGbPhenoData/GbPheno.csv", stringsAsFactors = F)
GbId <- GbPheno$ID
m.gb <- match(GbId, id.vec)

# ---------------------------------------------------------------------------- #
# ----- Ames 
# ---------------------------------------------------------------------------- #
# calculate numerator and denominator of GRM
ResList <- list()
for ( i in 1:n.chr ) {
   # load i-th chromosome
   f <- paste0("RAWDATA/GRM_12M_SNP/NumericalGenotypeData/merged_ames_282", i, "_maf.012")
   GenoData <- fread(f, data.table = F)
   GenoMat <- as.matrix(GenoData[, -1]) # first column is ID
   rm("GenoData"); gc(); gc() # garbage collection
   
   # check the number of SNPs
   checker <- df.n.snp$n.SNP[i] == ncol(GenoMat)
   if ( checker == FALSE ) { print("Warning: Genotype file has wrong # of SNPs"); break }
   
   # calculate numerator and denominator of GRM in this chromosome
   ResList[[i]] <- myFun.calc.GRM.info(X.012 = GenoMat[m.ames, ], batch = 5000)
   gc(); gc()
}

# calculate GRM
WWt.all <- matrix(0, nr = length(m.ames), nc = length(m.ames))
c.all <- 0
for ( i in 1:length(ResList) ) {
   WWt.mini <- ResList[[i]][[1]]
   c.mini <- ResList[[i]][[2]]
   WWt.all <- WWt.all + WWt.mini
   c.all <- c.all + c.mini
}
GRM <- WWt.all / c.all
rownames(GRM) <- colnames(GRM) <- id.vec[m.ames]

# write the GRM
write.csv(GRM, paste0(dir.save, "/GenReMat_use_12M_SNP_Ames.csv"))


# ---------------------------------------------------------------------------- #
# ----- 282 
# ---------------------------------------------------------------------------- #
# calculate numerator and denominator of GRM
ResList <- list()
for ( i in 1:n.chr ) {
   # load i-th chromosome
   f <- paste0("RAWDATA/GRM_12M_SNP/NumericalGenotypeData/merged_ames_282", i, "_maf.012")
   GenoData <- fread(f, data.table = F)
   GenoMat <- as.matrix(GenoData[, -1]) # first column is ID
   rm("GenoData"); gc(); gc() # garbage collection
   
   # check the number of SNPs
   checker <- df.n.snp$n.SNP[i] == ncol(GenoMat)
   if ( checker == FALSE ) { print("Warning: Genotype file has wrong # of SNPs"); break }
   
   # calculate numerator and denominator of GRM in this chromosome
   ResList[[i]] <- myFun.calc.GRM.info(X.012 = GenoMat[m.gb, ], batch = 5000)
   gc(); gc()
}

# calculate GRM
WWt.all <- matrix(0, nr = length(m.gb), nc = length(m.gb))
c.all <- 0
for ( i in 1:length(ResList) ) {
   WWt.mini <- ResList[[i]][[1]]
   c.mini <- ResList[[i]][[2]]
   WWt.all <- WWt.all + WWt.mini
   c.all <- c.all + c.mini
}
GRM <- WWt.all / c.all
rownames(GRM) <- colnames(GRM) <- id.vec[m.gb]

# write the GRM
write.csv(GRM, paste0(dir.save, "/GenReMat_use_12M_SNP_Gb.csv"))






