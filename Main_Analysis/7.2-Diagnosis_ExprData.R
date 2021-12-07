# 
library(data.table)
library(ggplot2)
library(gdata)

# 
dir.save <- "RESULT/7.2-Diagnosis_ExprData"
dir.create(dir.save)

#
AmesExpr <- fread("RAWDATA/ExprData_Ames/expression_BLUE_final_v1.1_B73.csv", data.table = F)
AmesExprMat <- as.matrix(AmesExpr[, 2:ncol(AmesExpr)])
rownames(AmesExprMat) <- AmesExpr$Accession_ID

#
GbExpr <- fread("RESULT/7.1-MakeExprData/GbExprData.csv", data.table = F)
GbExprMat <- as.matrix(GbExpr[, 3:ncol(GbExpr)])
rownames(GbExprMat) <- GbExpr$Name.Lipka

#
AvgExpr.per.sample.ames <- apply(AmesExprMat, 1, mean)
AvgExpr.per.sample.gb <- apply(GbExprMat, 1, mean)
df.fig <- rbind.data.frame(data.frame("Data" = "Ames",
																			"Avg.Expr" = AvgExpr.per.sample.ames),
													 data.frame("Data" = "Gb",
													 					 "Avg.Expr" = AvgExpr.per.sample.gb))
p <- ggplot(df.fig, aes(x = Data, y = Avg.Expr))
p <- p + geom_boxplot()
p <- p + xlab("Dataset")
p <- p + ylab("Average expression per sample")
p <- p + ggtitle("Compare the expression data in Ames and the 282 panel")
ggsave(filename = paste0(dir.save, "/Fig01-BoxplotA.png"), p, width = 6, height = 6)

#
AvgExpr.per.gene.ames <- apply(AmesExprMat, 2, mean)
AvgExpr.per.gene.gb <- apply(GbExprMat, 2, mean)
df.fig <- rbind.data.frame(data.frame("Data" = "Ames",
																			"Avg.Expr" = AvgExpr.per.gene.ames),
													 data.frame("Data" = "Gb",
													 					 "Avg.Expr" = AvgExpr.per.gene.gb))
p <- ggplot(df.fig, aes(x = Data, y = Avg.Expr))
p <- p + geom_boxplot()
p <- p + xlab("Dataset")
p <- p + ylab("Average expression per gene")
p <- p + ggtitle("Compare the expression data in Ames and the 282 panel")
ggsave(filename = paste0(dir.save, "/Fig02-BoxplotB.png"), p, width = 6, height = 6)

# 
pr.zero.per.gene.gb <- apply(GbExprMat == 0, 2, mean)
df.n.rm <- data.frame("pr.zero" = seq(0, 0.9, by = 0.1),
											"n.gene.rm" = NA)
for ( i in 1:nrow(df.n.rm) ) {
	cutoff <- df.n.rm$pr.zero[i]
	df.n.rm$n.gene.rm[i] <- sum(cutoff < pr.zero.per.gene.gb)
}
df.n.rm$pr.zero.prcnt <- paste0(df.n.rm$pr.zero * 100, "%")
p <- ggplot(df.n.rm, aes(x = pr.zero.prcnt, y = n.gene.rm))
p <- p + geom_bar(stat = "identity")
p <- p + geom_text(aes(label = n.gene.rm), vjust = -1)
p <- p + xlab("threshold: % of zero samples")
p <- p + ylab("# of removed genes")
p <- p + ylim(c(0, 2500))
p <- p + ggtitle("Gene filtering on the all 202 samples (in the 282 panel)")
ggsave(filename = paste0(dir.save, "/Fig03-Barplot.png"), p, width = 6, height = 4)


