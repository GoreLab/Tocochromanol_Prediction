# package
library(MTM)
library(gdata)
source("0.3-Functions_for_MTM.R")

# mkdir
dir.save <- "RESULT/4.4-Mtm_UseSomeGenes"
dir.create(dir.save, recursive = TRUE)

# load data
load("RESULT/4.3-MakeDataForMtm/MultiTraitPredictionData.RData")

# arg
method <- "MultiNamLargeGenes"
GeneList <- read.xls("RAWDATA/NamLargeGenes.xlsx")
GeneList$Trait <- as.character(GeneList$Trait)
GeneList$Gene <- as.character(GeneList$Gene)
GeneList$GeneName <- as.character(GeneList$GeneName)

#
trait.all <- as.character(unique(GeneList$Trait))
for ( tr in 1:length(trait.all) ) {
	# target trait
	trait <- trait.all[tr]
	
	# phenotype data
	pheno.raw <- data.frame("ID" = PhenoData$ID,
													"Y" = PhenoData[[trait]])
	covar.mat <- ExpMat.scaled[, GeneList$Gene[GeneList$Trait == trait], drop = F]
	colnames(covar.mat) <- GeneList$GeneName[GeneList$Trait == trait]
	pheno <- cbind(pheno.raw, covar.mat)
	rownames(pheno) <- NULL
	
	# data frame to save the result
	yPred.cv <- matrix(NA, nrow = nrow(CvMat), nc = ncol(CvMat))
	rownames(yPred.cv) <- rownames(CvMat)
	colnames(yPred.cv) <- colnames(CvMat)
	
	# Cross validation
	n.rep <- ncol(CvMat); n.fold <- max(CvMat)
	for ( i in 1:n.rep ) {
		for ( j in 1:n.fold ) {
			# time print
			b <- Sys.time()
			
			# mask test data
			test.name <- rownames(CvMat)[CvMat[, i] == j]
			
			# prediction
			yPred <- myFun.MTM.v2(PhenoData = pheno, test = test.name, K = G)
			yPred.cv[test.name, i] <- yPred
			
			# time print
			a <- Sys.time()
			print(a-b)
		}
	}
	
	# write
	f <- paste0(dir.save, "/CvRes_", trait, "_", method, ".csv")
	write.csv(yPred.cv, file = f)
	
	print(trait)
}


