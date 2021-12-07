# Genomic Prediction, Fit model & prediction

# Library & source
library(data.table)
library(BGLR)
library(ggplot2)

# variables
args <- commandArgs(trailingOnly=T)
trait <- args[1]
pop <- args[2] # "ames" or "gb"
model <- "GBLUP"
nIter <-  60000
burnIn <- 40000
thin <- 20

# ###########################################################
# trait <- "a.T3"
# pop <- "gb" # "ames" or "gb"
# model <- "GBLUP"
# nIter <-  600
# burnIn <- 400
# thin <- 4
# ###########################################################

# -------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------- load data  ---------------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #
# Make Folder to save result
name.dir <- paste0("RESULT/2.4-GenPreFit_SinglePop_trans_GBLUP_12M/Trait=", trait, "_Train=", pop, "_Model=", model)
dir.create(name.dir, recursive = TRUE)
name.dir.logfile <- paste0(name.dir, "/Logfile")
dir.create(name.dir.logfile)

# load GRM
GRM.df <- fread("RAWDATA/GRM_12M_SNP/GenReMat_use_12M_SNP.csv", data.table = F)
GRM <- as.matrix(GRM.df[, -1])
colnames(GRM) <- rownames(GRM) <- GRM.df[, 1]

# load phenotye 
if ( pop == "ames" ) { PhenoData <- read.csv("RESULT/1.1-MakeAmesPhenoData/AmesPheno_trans.csv", stringsAsFactors = F) }
if ( pop == "gb" ) { PhenoData <- read.csv("RESULT/1.2-MakeGbPhenoData/GbPheno_trans.csv", stringsAsFactors = F) }
pheno <- PhenoData[[trait]]

# load box-cox parameter
BoxCoxParam.ames <- read.csv("RESULT/1.1-MakeAmesPhenoData/Ames_BoxCoxParam.csv")
BoxCoxParam.gb <- read.csv("RESULT/1.2-MakeGbPhenoData/Gb_BoxCoxParam.csv")


# -------------------------------------------------------------------------------------------------------- #
# ---------------------------------------- run genomic prediction  --------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #
# make a phenotype vector
y <- rep(NA, nrow(GRM))
names(y) <- rownames(GRM)
m <- match(PhenoData$ID, names(y))
y[m] <- PhenoData[[trait]]

# fitting via BGLR
ETA <- list("G" = list(K = GRM, model = "RKHS"))
fm <- BGLR(y = y, ETA = ETA,
					 nIter = nIter, burnIn = burnIn, thin = thin,
					 saveAt = paste0(name.dir.logfile, "/fm_"), verbose = FALSE)
saveRDS(fm, file = paste0(name.dir, "/fm.Rdata")) # Save the result

# See result of fitting (1)
name.fig <- paste0(name.dir, "/Fig001-ObsFit.jpeg") # name of the figure
obs <- fm$y
fit <- fm$yHat
lim <- range(c(obs, fit))
COR <- round(cor(obs, fit), 3)
RMSE <- round(sqrt(mean((obs - fit) ^ 2)), 3)
HERIT <- round(var(fit) / var(obs), 3)
title.fig <- paste0("Fitting Accuracy; r = ", COR, ", RMSE = ", RMSE, ", G.HERIT = ", HERIT)
df.fig <- data.frame("obs" = obs, "fit" = fit)
p <- ggplot(data = df.fig, aes(x = obs, y = fit))
p <- p + geom_abline(slope = 1, intercept = 0, lty = 2, color = "gray")
p <- p + xlim(lim) + ylim(lim)
p <- p + ggtitle(label = title.fig)
p <- p + geom_point()
ggsave(p, filename = name.fig)

# See result of fitting (2)
name.fig <- paste0(name.dir, "/Fig002-FitSd.jpeg") # name of the figure
title.fig <- "Fitted value and its posterior StDev"
StDev <- fm$SD.yHat
df.fig <- data.frame("fit" = fit, "StDev" = StDev)
p <- ggplot(data = df.fig, aes(x = fit, y = StDev))
p <- p + ggtitle(label = title.fig)
p <- p + geom_point()
ggsave(p, filename = name.fig)

# Check convergence (3)
name.fig <- paste0(name.dir, "/Fig004-TracePlot_varE.jpeg") # varE
name.dat <- paste0(name.dir, "/Logfile/fm_varE.dat")
varE.mcmc <- scan(name.dat)
iter <- seq(from = thin, to = nIter, by = thin)
sample.burn <- burnIn
df.fig <- data.frame("iter" = iter, "varE" = varE.mcmc)
p <- ggplot(data = df.fig, aes(x = iter, y = varE))
p <- p + geom_point() + geom_line()
p <- p + geom_vline(xintercept = sample.burn, lty = 3, lwd = 1)
ggsave(p, file = name.fig)


# ----------------------------------------------------------------------------------------------- #
# --------------------------------------- inverse box-cox --------------------------------------- #
# ----------------------------------------------------------------------------------------------- #
PredValue <- fm$yHat
if ( pop == "ames" ) {
	lambda <- BoxCoxParam.ames$lambda[BoxCoxParam.ames$trait == trait]
	const <- BoxCoxParam.ames$const[BoxCoxParam.ames$trait == trait]
} else {
	lambda <- BoxCoxParam.gb$lambda[BoxCoxParam.gb$trait == trait]
	const <- BoxCoxParam.gb$const[BoxCoxParam.gb$trait == trait]
}
tmp.vec <- as.numeric(PredValue)
if (lambda == 0) { tmp.vec.original <- exp(tmp.vec) } else { tmp.vec.original <- tmp.vec ^ (1/lambda) }
tmp.vec.original <- tmp.vec.original - const
PredValue <- tmp.vec.original

# save
filename.save <- paste0("RESULT/2.4-GenPreFit_SinglePop_trans_GBLUP_12M/yPred_trait=", trait, "_train=", pop, "_model=", model, ".csv")
df <- data.frame("ID" = names(fm$yHat), "pred" = PredValue)
write.csv(df, filename.save, row.names = F)



