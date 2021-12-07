# 
library(data.table)
library(ggplot2)

# 
dir.save <- "RESULT/7.4-Diagnosis_ExprData"
dir.create(dir.save)

# params
cutoff.nonzero <- 0.5 # at least 50% samples should be non-zero
cutoff.MAD <- 100

# load data
ExprDat <- fread("RESULT/7.1-MakeExprData/GbExprData.csv", data.table = F)
ExprMat <- as.matrix(ExprDat[, 3:ncol(ExprDat)])
rownames(ExprMat) <- ExprDat$Name.Lipka

# 1. gene filtering
pr.nonzero.rlog <- 1 - apply(ExprMat == 0, 2, mean)
tf <- cutoff.nonzero <= pr.nonzero.rlog
ExprMat <- ExprMat[, tf] # filter

# MAD calculation
MAD.matrix <- matrix(NA, nr = nrow(ExprMat), nc = ncol(ExprMat),
										 dimnames = dimnames(ExprMat))
for ( i in 1:ncol(ExprMat) ) {
	x <- ExprMat[, i]
	thres.i <- mad(x, constant = 1)
	mad.i <- abs(x - median(x)) / thres.i
	MAD.matrix[, i] <- mad.i
}

# save MAD matrix
MAD.Dat <- data.frame("Sample_ID" = ExprDat$Name.Lipka, MAD.matrix)
fwrite(MAD.Dat, file = paste0(dir.save, "/MAD_matrix.csv"))

# remove observations when MAD > cutoff
bool.mat <- cutoff.MAD < MAD.matrix
ExprMat.RmOut <- ExprMat
ExprMat.RmOut[bool.mat] <- NA
print(sum(bool.mat))

# write the number of outliers
n.rm <- apply(bool.mat, 2, sum)
tf <- nrow(ExprMat.RmOut) * 0.1 < n.rm
tab.n.rm <- table(n.rm)
df.save.01 <- data.frame("GeneID" = names(n.rm),
												 "n.outlier" = n.rm,
												 "flag" = c("keep", "remove")[as.numeric(tf)+1],
												 row.names = NULL)
fwrite(df.save.01, file = paste0(dir.save, "/number_of_outliers.csv"))

# write summary table of the outlier removal
df.save.02 <- data.frame(tab.n.rm)
colnames(df.save.02) <- c("Number_of_removed_outliers", "Freq")
fwrite(df.save.02, file = paste0(dir.save, "/summary_table.csv"))

# when more than than 10% of the samples were identified as outliers, remove the gene
tf <- nrow(ExprMat.RmOut) * 0.1 < n.rm
genes.remove <- colnames(ExprMat.RmOut)[tf]
ExprMat.RmOut.RmGene <- ExprMat.RmOut[, !tf]

# imputation by using median
tf <- apply(is.na(ExprMat.RmOut.RmGene), 2, sum) != 0
genes.with.na <- colnames(ExprMat.RmOut.RmGene)[tf]
ExprMat.Imp <- ExprMat.RmOut.RmGene
for ( j in 1:length(genes.with.na) ) {
	gene.j <- genes.with.na[j]
	expr.vec <- ExprMat.RmOut.RmGene[, gene.j]
	med.j <- median(expr.vec, na.rm = T)
	ExprMat.Imp[is.na(expr.vec), gene.j] <- med.j
}

# write the outlier-removed & median-imputed matrix
ExprDat.Imp <- data.frame("Sample_ID" = ExprDat$Name.Lipka, ExprMat.Imp)
fwrite(ExprDat.Imp, file = paste0(dir.save, "/ExpressionData.Filtered.and.Imputed.csv"))
