# Genomic Prediction, Fit model

# Library & source
library(data.table)
library(BGLR)
library(ggplot2)

# variables
args <- commandArgs(trailingOnly=T)
trait <- args[1]
pop <- args[2] # "ames" or "gb"
model <- args[3] # "BRR" or "BayesB" (no 'GBLUP' option!)
nIter <-  60000
burnIn <- 40000
thin <- 20

# ###########################################################
# trait <- "a.T"
# pop <- "gb" # "ames" or "gb"
# model <- "BayesB" # "BRR" or "BayesB" (no 'GBLUP' option!)
# nIter <-  600
# burnIn <- 400
# thin <- 4
# ###########################################################


# -------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------- load data  ---------------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #
# Make Folder to save result
name.dir <- paste0("RESULT/2.1-GenPreFit_SinglePop_trans/Trait=", trait, "_Train=", pop, "_Model=", model)
dir.create(name.dir, recursive = TRUE)
name.dir.logfile <- paste0(name.dir, "/Logfile")
dir.create(name.dir.logfile)

# load geotype & map
indv.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.indv", header = FALSE)
geno.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012")
geno.mat <- as.matrix(geno.tmp[, -1]) # 1st col is ID.
rownames(geno.mat) <- indv.tmp$V1
map.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.pos")
map <- data.frame("chr" = map.tmp$V1, "pos" = map.tmp$V2)
rm(indv.tmp); rm(geno.tmp); rm(map.tmp); gc(); gc() # remove non-essential object

# load phenotye 
if ( pop == "ames" ) { PhenoData <- read.csv("RESULT/1.1-MakeAmesPhenoData/AmesPheno_trans.csv", stringsAsFactors = F) }
if ( pop == "gb" ) { PhenoData <- read.csv("RESULT/1.2-MakeGbPhenoData/GbPheno_trans.csv", stringsAsFactors = F) }
pheno <- PhenoData[[trait]]



# -------------------------------------------------------------------------------------------------------- #
# ---------------------------------------- run genomic prediction  --------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #
# fitting via BGLR
ETA <- list("MRK" = list(X = score, model = model, saveEffects = TRUE))
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


# Check convergence (1)
if ( model == "BayesB" ) {
   iter <- seq(from = thin, to = nIter, by = thin)
   sample.burn <- burnIn
   name.fig <- paste0(name.dir, "/Fig003-TracePlot_probIn.jpeg") # probIn param in BayesB
   name.dat <- paste0(name.dir, "/Logfile/fm_ETA_MRK_parBayesB.dat")
   parBayesB.mcmc <- read.table(name.dat, header = T)
   df.fig <- data.frame("iter" = iter, "probIn" = parBayesB.mcmc$probIn)
   p <- ggplot(data = df.fig, aes(x = iter, y = probIn))
   p <- p + geom_point() + geom_line()
   p <- p + geom_vline(xintercept = sample.burn, lty = 3, lwd = 1)
   ggsave(p, file = name.fig)
}

# Check convergence (2)
if ( model == "BayesB" ) {
   iter <- seq(from = thin, to = nIter, by = thin)
   sample.burn <- burnIn
   name.fig <- paste0(name.dir, "/Fig004-TracePlot_scale.jpeg") # scale param in BayesB
   df.fig <- data.frame("iter" = iter, "scale" = parBayesB.mcmc$scale)
   p <- ggplot(data = df.fig, aes(x = iter, y = scale))
   p <- p + geom_point() + geom_line()
   p <- p + geom_vline(xintercept = sample.burn, lty = 3, lwd = 1)
   ggsave(p, file = name.fig)
}

# Check convergence (3)
if ( model == "BRR" ) {
   iter <- seq(from = thin, to = nIter, by = thin)
   sample.burn <- burnIn
   name.fig <- paste0(name.dir, "/Fig005-TracePlot_varB.jpeg")
   name.dat <- paste0(name.dir, "/Logfile/fm_ETA_MRK_varB.dat") # varB in BRR
   varB.mcmc <- scan(name.dat)
   df.fig <- data.frame("iter" = iter, "varB" = varB.mcmc)
   p <- ggplot(data = df.fig, aes(x = iter, y = varB))
   p <- p + geom_point() + geom_line()
   p <- p + geom_vline(xintercept = sample.burn, lty = 3, lwd = 1)
   ggsave(p, file = name.fig)
}


# Check convergence (4)
name.fig <- paste0(name.dir, "/Fig006-TracePlot_varE.jpeg") # varE
name.dat <- paste0(name.dir, "/Logfile/fm_varE.dat")
varE.mcmc <- scan(name.dat)
df.fig <- data.frame("iter" = iter, "varE" = varE.mcmc)
p <- ggplot(data = df.fig, aes(x = iter, y = varE))
p <- p + geom_point() + geom_line()
p <- p + geom_vline(xintercept = sample.burn, lty = 3, lwd = 1)
ggsave(p, file = name.fig)
