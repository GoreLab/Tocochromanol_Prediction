# make fold of cross validation

# set seed
set.seed(143)

# mkdir
dir.create("RESULT/1.6-MakeCvFold_Exp")

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

# Make fold for Exp.data 
GenReMat <- as.matrix(read.csv("RESULT/1.5-MakeGenReMat/GenReMat_ForExpData.csv", row.names = 1))
CvMat <- myfun.MakeCvmat(N = nrow(GenReMat), n.fold = n.fold, n.rep = n.rep)
df.CvMat <- data.frame("ID" = rownames(GenReMat), CvMat)
write.csv(df.CvMat, file = "RESULT/1.6-MakeCvFold_Exp/CvFold_Exp.csv")
