# make fold of cross validation

# set seed
set.seed(637)

# mkdir
dir.create("RESULT/1.4-MakeCvFold")

# params
n.fold <- 5
n.rep <- 10

# fucntion 
myfun.MakeCvmat <- function(N, n.fold, n.rep) {
   CvMat <- matrix(NA, nr = N, nc = n.rep)
   for (r in 1:n.rep) {
      num <- rep(1:n.fold, length.out = N)
      rand.num <- sample(num, replace = FALSE)
      CvMat[, r] <- rand.num
   }
   colnames(CvMat) <- paste0("rep", formatC(1:n.rep, width = 2, flag = "0"))
   return(CvMat)
}

# 1. Ames
pheno <- read.csv("RESULT/1.1-MakeAmesPhenoData/AmesPheno.csv", stringsAsFactors = F)
CvMat <- myfun.MakeCvmat(N = nrow(pheno), n.fold = n.fold, n.rep = n.rep)
df.CvMat <- data.frame("ID" = pheno$ID, CvMat)
write.csv(df.CvMat, file = "RESULT/1.4-MakeCvFold/CvFold_ames.csv")


# 2. Gb
pheno <- read.csv("RESULT/1.2-MakeGbPhenoData/GbPheno.csv", stringsAsFactors = F)
CvMat <- myfun.MakeCvmat(N = nrow(pheno), n.fold = n.fold, n.rep = n.rep)
df.CvMat <- data.frame("ID" = pheno$ID, CvMat)
write.csv(df.CvMat, file = "RESULT/1.4-MakeCvFold/CvFold_gb.csv")



