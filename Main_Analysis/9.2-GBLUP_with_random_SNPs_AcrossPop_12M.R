# packages
library(data.table)
library(rrBLUP)
library(BGLR)

# parameter
args <- commandArgs(trailingOnly = TRUE)
p <- as.numeric(args[1])
r <- as.numeric(args[2])
nIter <-  60000
burnIn <- 40000
thin <- 20
seed <- 20210629
n.chr.all <- 10
set.seed(seed)

# function to use
myFun.calc.GRM.info <- function(X.012, batch = 3000, ...) {
   N <- nrow(X.012)
   P <- ncol(X.012)
   n <- ceiling(P / batch)
   c.all <- rep(0, n)
   WWt.all <- matrix(0, nr = N, nc = N)
   for (i in 1:n) {
      start <- (i - 1) * batch + 1
      end <- min(i * batch,  P)
      X.mini <- X.012[, start:end]
      vec.tmp <- apply(X.mini, 2, sum)
      allele.freq.vec <- vec.tmp / (2 * N)
      allele.freq.mat <- t(allele.freq.vec) %x% rep(1, N)
      c.mini <- 2 * sum(allele.freq.vec * (1 - allele.freq.vec))
      W <- X.mini - 2 * allele.freq.mat
      WWt.mini <- tcrossprod(W)
      c.all[i] <- c.mini
      WWt.all <- WWt.all + WWt.mini
   }
   # return list
   ResList <- list("Numerator" = WWt.all, "Denominator" = sum(c.all))
   return(ResList)
}

# create dir to save result
folder.save <- paste0("RESULT/9.2-GBLUP_with_random_SNPs_AcrossPop_12M")
dir.create(folder.save, recursive = T)
dir.create(paste0(folder.save, "/LOGFILE"), recursive = T)

# print
m.out <- paste0("Start GBLUP with random SNPs: #SNP = ", as.integer(p))
print(m.out)

# load ID of accessions
df.id <- read.table("RAWDATA/GRM_12M_SNP/NumericalGenotypeData/merged_ames_2821_maf.012.indv") # get id
id.vec <- df.id$V1

# load phenotype 
AmesPheno <- read.csv("RESULT/1.1-MakeAmesPhenoData/AmesPheno_trans.csv", stringsAsFactors = F)
GbPheno <- read.csv("RESULT/1.2-MakeGbPhenoData/GbPheno_trans.csv", stringsAsFactors = F)

# ---------------------------------------------------------------------------- #
# ----- Calculate GRM with random X SNPs
# ---------------------------------------------------------------------------- #
# print
m.out <- paste0("Start the random sampling of ", as.integer(p), " SNPs"); print(m.out)
print(m.out)

# get the number of SNPs in each chromosome
df.n.snp <- data.frame("chr" = 1:n.chr.all, "n.SNP" = NA)
for ( chr in 1:n.chr.all ) {
   # load SNP info with features
   map.w.f.chr <- fread(paste0("RESULT/1.7-UseFeatures/Use_12M/map_with_feature_chr", chr, ".csv"), data.table = F)
   n.snp <- nrow(map.w.f.chr)
   df.n.snp$n.SNP[chr] <- n.snp
}
df.n.snp$start <- c(1, (cumsum(df.n.snp$n.SNP) + 1)[1:(n.chr.all - 1)])
df.n.snp$end <- cumsum(df.n.snp$n.SNP)

# Get random SNPs
P <- sum(df.n.snp$n.SNP)
SnpNumList <- list()
for ( i in 1:r ) {
   SnpNumList[[i]] <- list()
   rand.num <- sort(sample(x = 1:P, size = p, replace = F))
   for ( chr in 1:nrow(df.n.snp) ) {
      s <- df.n.snp$start[chr]
      e <- df.n.snp$end[chr]
      tf <- (s <= rand.num) & (rand.num <= e)
      rand.num.chr <- rand.num[tf]
      rand.num.chr.new <- rand.num.chr - s + 1
      SnpNumList[[i]][[paste0("chr", chr)]] <- rand.num.chr.new
   }
}
names(SnpNumList) <- paste0("Rep", formatC(1:r, width = 2, flag = 0))
saveRDS(SnpNumList, file = paste0(folder.save, "/SnpNumList_", as.integer(p), ".Rdata")) # save

# print
m.out <- paste0("End the random sampling"); print(m.out)

# print
m.out <- paste0("Start the calculation of GRM"); print(m.out)

# make kernel and write them
dir.log <- paste0(folder.save, "/LogFile_", as.integer(p))
dir.create(dir.log, recursive = T)
for ( i in 1:r ) {
   # print
   m.out <- paste0("Start Rep.", i, " out of ", r, " reps"); print(m.out)
   # object to save n/d for GRM
   ResList <- list()
   for ( chr in 1:n.chr.all ) {
      # random SNP
      num.snp.i.chr <- SnpNumList[[i]][[chr]]
      # print
      m.out <- paste0("Loading the genotype data for chr.", chr); print(m.out)
      # load i-th chromosome
      f <- paste0("RAWDATA/GRM_12M_SNP/NumericalGenotypeData/merged_ames_282", chr, "_maf.012")
      GenoData <- fread(f, data.table = F)
      GenoMat <- as.matrix(GenoData[, -1]) # first column is ID
      rm("GenoData"); gc(verbose = F); gc(verbose = F)
      # print
      m.out <- paste0("End loading"); print(m.out)
      m.out <- paste0("We use ", as.integer(length(num.snp.i.chr)), "SNPs for chr.", chr, " in rep.", i); print(m.out)
      m.out <- paste0("Calculate GRM by using the batch-mode function"); print(m.out)
      # calculate numerator and denominator of GRM in this chromosome
      b <- Sys.time()
      ResList[[chr]] <- myFun.calc.GRM.info(X.012 = GenoMat[, num.snp.i.chr], batch = 5000)
      gc(verbose = F); gc(verbose = F)
      a <- Sys.time()
      # print
      m.out <- paste0("GRM calculation for chr.", chr, " in rep.", i, " Done"); print(m.out)
      print(a - b)
   }
   # calculate GRM
   WWt.all <- matrix(0, nr = nrow(GenoMat), nc = nrow(GenoMat))
   c.all <- 0
   for ( j in 1:length(ResList) ) {
      WWt.mini <- ResList[[j]][[1]]
      c.mini <- ResList[[j]][[2]]
      WWt.all <- WWt.all + WWt.mini
      c.all <- c.all + c.mini
   }
   GRM <- WWt.all / c.all
   rownames(GRM) <- colnames(GRM) <- id.vec
   # save GRM
   write.csv(GRM, paste0(dir.log, "/GRM_with_", as.integer(p), "_SNPs_rep", i, ".csv"))
}
rm("GenoMat"); gc(verbose = F); gc(verbose = F) # clean-up

# print
m.out <- paste0("All GRM calculations are done. Start across-population GP"); print(m.out)

# ---------------------------------------------------------------------------- #
# ----- loop for all reps ---------------------------------------------------- #
# ---------------------------------------------------------------------------- #
trait.all <- c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3",
               "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")
for ( trait in trait.all ) {
   df.each.use.Ames <- df.each.use.282 <- data.frame("GBS.ID" = id.vec)
   for ( j in 1:length(SnpNumList) ) {
      # numbers to indicate SNPs for the i-th random set of SNPs
      name.j <- names(SnpNumList)[j]
      
      # relationship matrix
      f.in <- paste0(dir.log, "/GRM_with_", as.integer(p), "_SNPs_rep", j, ".csv")
      GenReMat <- as.matrix(read.csv(f.in, row.names = 1))
      
      # -----  1. Ames -> Gb prediction
      # build model
      y <- setNames(object = rep(NA, times = length(id.vec)), nm = id.vec)
      y[AmesPheno$ID] <- AmesPheno[[trait]]
      fm <- BGLR(y = y, ETA = list("K" = list("K" = GenReMat, model = "RKHS")), 
                 nIter = nIter, burnIn = burnIn, thin = thin, verbose = F,
                 saveAt = paste0(folder.save, "/LOGFILE/fm_"))
      df.each.use.Ames[[name.j]] <- fm$yHat
      
      # -----  2. Gb -> Ames prediction
      y <- setNames(object = rep(NA, times = length(id.vec)), nm = id.vec)
      y[GbPheno$ID] <- GbPheno[[trait]]
      fm <- BGLR(y = y, ETA = list("K" = list("K" = GenReMat, model = "RKHS")), 
                 nIter = nIter, burnIn = burnIn, thin = thin, verbose = F, 
                 saveAt = paste0(folder.save, "/LOGFILE/fm_"))
      df.each.use.282[[name.j]] <- fm$yHat
   }
   
   # write out the result
   write.csv(df.each.use.Ames, 
             file = paste0(folder.save, "/PredRes_AmesToGb_RandomSnps_", 
                           as.integer(p), "_", trait, "_trans.csv"),
             row.names = F)
   write.csv(df.each.use.282, 
             file = paste0(folder.save, "/PredRes_GbToAmes_RandomSnps_", 
                           as.integer(p), "_", trait, "_trans.csv"),
             row.names = F)
}






