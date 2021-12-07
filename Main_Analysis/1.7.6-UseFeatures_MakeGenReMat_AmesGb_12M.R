# library
library(data.table)

# mkdir
dir.save <- "RESULT/1.7-UseFeatures/Use_12M"
dir.create(dir.save, recursive = TRUE)

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

# setup
df.comb.raw <- data.frame("Feature" = c("MAF", "MAF",
                                        "REC.RATE", "REC.RATE", "REC.RATE",
                                        "GERP", "GERP",
                                        "PROX.1Kb", "PROX.1Kb", 
                                        "PROX.10Kb", "PROX.10Kb",
                                        "PROX.100Kb", "PROX.100Kb"),
                          "Group" = c("Middle", "High",
                                      "Low", "Middle", "High",
                                      "Positive", "Negative",
                                      "Prox", "Dist", 
                                      "Prox", "Dist", 
                                      "Prox", "Dist"),
                          stringsAsFactors = F)
df.comb.ames <- cbind(df.comb.raw, "Pop" = "Ames")
df.comb.gb <- cbind(df.comb.raw, "Pop" = "Gb")
df.comb <- rbind.data.frame(df.comb.ames, df.comb.gb)
chr.all <- 1:10

# get ID of all
df.id <- read.table("RAWDATA/GRM_12M_SNP/NumericalGenotypeData/merged_ames_2821_maf.012.indv") # get id
id.vec <- df.id$V1

# get ID of Ames
AmesData <- fread("RESULT/1.1-MakeAmesPhenoData/AmesPheno_trans.csv")
AmesId <- AmesData$ID
AmesId.common <- intersect(id.vec, AmesId)

# get ID of Gb
GbData <- fread("RESULT/1.2-MakeGbPhenoData/GbPheno_trans.csv")
GbId <- GbData$ID
GbId.common <- intersect(id.vec, GbId)

# loop for all
for ( i in 1:nrow(df.comb) ) {
   # parameters
   feat <- df.comb$Feature[i]
   grp <- df.comb$Group[i]
   pop <- df.comb$Pop[i]
   if ( pop == "Ames" ) { m <- match(AmesId.common, id.vec) }
   if ( pop == "Gb" ) { m <- match(GbId.common, id.vec) }
   
   # calculate numerator and denominator of GRM
   ResList <- list()
   for ( chr in chr.all ) {
      # map info
      f <- paste0("RESULT/1.7-UseFeatures/Use_12M/map_with_feature_chr", chr, ".csv")
      map <- fread(file = f, data.table = F)
      
      # load i-th chromosome
      f <- paste0("RAWDATA/GRM_12M_SNP/NumericalGenotypeData/merged_ames_282", chr, "_maf.012")
      GenoData <- fread(f, data.table = F)
      GenoMat <- as.matrix(GenoData[, -1]) # first column is ID
      rm("GenoData"); gc(); gc() # garbage collection
      
      # get indicator
      vec.feat <- map[[feat]]
      if ( feat == "MAF" & grp == "Middle") { tf <- (-0.01 < vec.feat) & vec.feat <= 0.05 }
      if ( feat == "MAF" & grp == "High") { tf <- (0.05 < vec.feat) & (vec.feat <= 0.5 ) }
      if ( feat == "REC.RATE" & grp == "Low") { tf <- vec.feat * 1000000 * 100 <= 0.45 }
      if ( feat == "REC.RATE" & grp == "Middle") { tf <- (0.45 < vec.feat * 1000000 * 100) & (vec.feat * 1000000 * 100 <= 1.65) }
      if ( feat == "REC.RATE" & grp == "High") { tf <- 1.65 < vec.feat * 1000000 * 100 }
      if ( feat == "GERP" & grp == "Positive") { tf <- 0 < vec.feat }
      if ( feat == "GERP" & grp == "Negative") { tf <- vec.feat <= 0 }
      if ( feat == "PROX.1Kb" & grp == "Prox") { tf <- vec.feat }
      if ( feat == "PROX.1Kb" & grp == "Dist") { tf <- !vec.feat }
      if ( feat == "PROX.10Kb" & grp == "Prox") { tf <- vec.feat }
      if ( feat == "PROX.10Kb" & grp == "Dist") { tf <- !vec.feat }
      if ( feat == "PROX.100Kb" & grp == "Prox") { tf <- vec.feat }
      if ( feat == "PROX.100Kb" & grp == "Dist") { tf <- !vec.feat }
      
      # print
      print(paste0("We use ", sum(tf), "SNPs for ", feat, "-", grp, " in chr", chr))
      
      # calculate numerator and denominator of GRM in this chromosome
      ResList[[chr]] <- myFun.calc.GRM.info(X.012 = GenoMat[m, tf], batch = 5000)
      gc(); gc()
   }
   
   # calculate GRM
   WWt.all <- matrix(0, nr = length(m), nc = length(m))
   c.all <- 0
   for ( j in 1:length(ResList) ) {
      WWt.mini <- ResList[[j]][[1]]
      c.mini <- ResList[[j]][[2]]
      WWt.all <- WWt.all + WWt.mini
      c.all <- c.all + c.mini
   }
   GRM <- WWt.all / c.all
   rownames(GRM) <- colnames(GRM) <- id.vec[m]
   
   # write the GRM
   f.save <- paste0(dir.save, "/GenReMat_use_12M_SNP_", feat, "_", grp, "_For", pop, ".csv")
   write.csv(GRM, f.save)
}


