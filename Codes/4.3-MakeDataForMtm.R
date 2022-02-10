# packages
library(data.table)

# mkdir
dir.save <- "RESULT/4.3-MakeDataForMtm"
dir.create(dir.save, recursive = TRUE)

# ---------------------------------------------------------------------------- #
# genomic relationship matrix
f <- "RESULT/1.5-MakeGenReMat/GenReMat_ForExpData.csv"
G <- as.matrix(read.csv(f, row.names = 1))

# # gene expression data
# f <- "RAWDATA/ExprData_Ames/expression_BLUE_final_v1.1_B73.csv"
# ExpRawData <- fread(f, data.table = F)
# ExpMat <- as.matrix(ExpRawData[, -1])
# rownames(ExpMat) <- gsub("_", "", ExpRawData$Accession_ID)
# ExpMat <- ExpMat[match(rownames(G), rownames(ExpMat)), ]
# ExpMat.scaled <- scale(ExpMat, center = T, scale = T) # scale each column
# 
# 
# gene expression data (all + vte7 expression)
Vte7ExprData <- fread("RAWDATA/ExprData_Ames/BLUE_vte7.csv", data.table = F)
ExpRawData <- fread("RAWDATA/ExprData_Ames/expression_BLUE_final_v1.1_B73.csv", data.table = F)

# merge vte7 with others
m <- match(ExpRawData$Accession_ID, Vte7ExprData$Genotype)
ExpRawData[["Zm00001d006778"]] <- Vte7ExprData$BLUE[m]
ExpMat <- as.matrix(ExpRawData[, -1])
rownames(ExpMat) <- gsub("_", "", ExpRawData$Accession_ID)
ExpMat <- ExpMat[match(rownames(G), rownames(ExpMat)), ]
ExpMat.scaled <- scale(ExpMat, center = T, scale = T) # scale each column


# Ames phenotype data (trans)
f <- "RESULT/1.1-MakeAmesPhenoData/AmesPheno_trans.csv"
PhenoData <- read.csv(f, stringsAsFactors = F)
PhenoData <- PhenoData[match(rownames(G), PhenoData$ID), ]

# Ames Phenotype data (un-trans)
f <- "RESULT/1.1-MakeAmesPhenoData/AmesPheno.csv"
PhenoData.ut <- read.csv(f, stringsAsFactors = F)
PhenoData.ut <- PhenoData.ut[match(rownames(G), PhenoData.ut$ID), ]

# CV fold
f <- "RESULT/1.6-MakeCvFold_Exp/CvFold_Exp.csv"
CvFold <- read.csv(f, row.names = 1)
CvMat <- as.matrix(CvFold[, -1])
rownames(CvMat) <- CvFold[, 1]

# box-cox trans
f <- "RESULT/1.1-MakeAmesPhenoData/Ames_BoxCoxParam.csv"
BoxCoxParam.ames <- read.csv(f)

# make an R data object for convenience
save(G, ExpMat.scaled, PhenoData, PhenoData.ut, CvMat, BoxCoxParam.ames,
		 file = paste0(dir.save, "/MultiTraitPredictionData.RData"))


# # ---------------------------------------------------------------------------- #
# dim(ExpMat.scaled)
# dim(G)
# dim(PhenoData)
# dim(PhenoData.ut)
# dim(CvMat)
# dim(BoxCoxParam.ames)
# 
# lambda <- BoxCoxParam.ames$lambda[BoxCoxParam.ames$trait == trait]
# const <- BoxCoxParam.ames$const[BoxCoxParam.ames$trait == trait]
# 
# pheno.vec <- PhenoData[[trait]]
# names(pheno.vec) <- PhenoData$ID
# 
# pheno.vec.ut <- PhenoData.ut[[trait]]
# names(pheno.vec.ut) <- PhenoData.ut$ID









