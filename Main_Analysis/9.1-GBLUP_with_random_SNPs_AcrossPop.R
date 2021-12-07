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
set.seed(seed)

# create dir to save result
folder.save <- paste0("RESULT/9.1-GBLUP_with_random_SNPs_AcrossPop")
dir.create(folder.save, recursive = T)
dir.create(paste0(folder.save, "/LOGFILE"), recursive = T)

# load genotype & map
indv.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.indv", header = FALSE)
geno.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012")
geno.mat <- as.matrix(geno.tmp[, -1]) # 1st col is ID.
names.genotypes <- rownames(geno.mat) <- indv.tmp$V1
map.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.pos")
map <- data.frame("chr" = map.tmp$V1, "pos" = map.tmp$V2)
rm(indv.tmp); rm(geno.tmp); rm(map.tmp); gc(); gc() # remove non-essential object

# load phenotype 
AmesPheno <- read.csv("RESULT/1.1-MakeAmesPhenoData/AmesPheno_trans.csv", stringsAsFactors = F)
GbPheno <- read.csv("RESULT/1.2-MakeGbPhenoData/GbPheno_trans.csv", stringsAsFactors = F)

# Get random SNPs
P <- ncol(geno.mat)
SnpNumList <- list()
for ( i in 1:r ) {
   rand.num <- sort(sample(x = 1:P, size = p, replace = F))
   SnpNumList[[i]] <- rand.num
}
names(SnpNumList) <- paste0("Rep", formatC(1:r, width = 2, flag = 0))

# save the object of SNPs
saveRDS(SnpNumList, file = paste0(folder.save, "/SnpNumList_", as.integer(p), ".Rdata"))

# make kernel and write them
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
dir.log <- paste0(folder.save, "/LogFile_", as.integer(p))
dir.create(dir.log, recursive = T)
for ( i in 1:r ) {
   num.snp.i <- SnpNumList[[i]]
   GRM.i <- myFun.calc.G(X.012 = geno.mat[, num.snp.i], batch = 100)
   write.csv(GRM.i, paste0(dir.log, "/GRM_with_", as.integer(p), "_SNPs_rep", i, ".csv"))
}
rm(geno.mat); gc(); gc() # clean-up


# ---------------------------------------------------------------------------- #
# ----- loop for all reps ---------------------------------------------------- #
# ---------------------------------------------------------------------------- #
trait.all <- c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3",
               "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")
for ( trait in trait.all ) {
   df.each.use.Ames <- df.each.use.282 <- data.frame("GBS.ID" = names.genotypes)
   for ( i.qtl in 1:length(SnpNumList) ) {
      # numbers to indicate SNPs for the i-th random set of SNPs
      name.qtl.i <- names(SnpNumList)[i.qtl]
      
      # relationship matrix
      f.in <- paste0("RESULT/9.1-GBLUP_with_random_SNPs_AcrossPop/LogFile_", as.integer(p), 
                     "/GRM_with_", as.integer(p), "_SNPs_rep", i.qtl, ".csv")
      GenReMat <- as.matrix(read.csv(f.in, row.names = 1))
      
      # -----  1. Ames -> Gb prediction
      # build model
      y <- setNames(object = rep(NA, times = length(names.genotypes)), nm = names.genotypes)
      y[AmesPheno$ID] <- AmesPheno[[trait]]
      fm <- BGLR(y = y, ETA = list("K" = list("K" = GenReMat, model = "RKHS")), 
                 nIter = nIter, burnIn = burnIn, thin = thin, verbose = F,
                 saveAt = paste0(folder.save, "/LOGFILE/fm_"))
      df.each.use.Ames[[name.qtl.i]] <- fm$yHat
      
      # -----  2. Gb -> Ames prediction
      y <- setNames(object = rep(NA, times = length(names.genotypes)), nm = names.genotypes)
      y[GbPheno$ID] <- GbPheno[[trait]]
      fm <- BGLR(y = y, ETA = list("K" = list("K" = GenReMat, model = "RKHS")), 
                 nIter = nIter, burnIn = burnIn, thin = thin, verbose = F, 
                 saveAt = paste0(folder.save, "/LOGFILE/fm_"))
      df.each.use.282[[name.qtl.i]] <- fm$yHat
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
