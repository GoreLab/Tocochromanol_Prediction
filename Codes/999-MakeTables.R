# Make tables for the manuscript
# Table 1, S1, S3, and S4 were made separately, as they are independent from the result of the prediction analysis 

# 
dir.save <- "RESULT/Tables_and_Figures"
dir.create(dir.save, recursive = T)

# ---------------------------------------------------------------------------- #
# ----- Table S2
# ---------------------------------------------------------------------------- #
# define objects
tr.all <- c("a.T", "d.T", "g.T",
            "a.T3", "d.T3", "g.T3",
            "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")
mod.all <- c("250K", "1M", "SI")

# Get map for the 341K SNP set
df.map <- read.table("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.pos")
colnames(df.map) <- c("chr", "pos")

# make table S2
df.S2 <- NULL
for ( i in 1:length(tr.all) ) {
  for (j in 1:length(mod.all) ) {
    tr <- tr.all[i]
    mod <- mod.all[j]
    dat <- readRDS(paste0("RESULT/5.1-MultiKernel_Across_trans/SnpNumList_", mod, "_", tr, ".Rdata"))
    dat.qtl.only <- dat[names(dat) != "NonQtl"]
    n.SNPs.vec <- sapply(dat.qtl.only, length, USE.NAMES = F)
    width.all <- rep(NA, length(dat.qtl.only))
    for ( k in 1:length(dat.qtl.only) ) {
      R <- range(df.map$pos[dat.qtl.only[[k]]])
      width.all[k] <- (R[2] - R[1]) / 1000000
    }
    df.ij <- data.frame("Phenotype" = tr, "Window" = mod, "n.QTL" = length(dat.qtl.only),
                        "Avg.n.SNP" = round(mean(n.SNPs.vec), 1), "Min.n.SNP" = min(n.SNPs.vec), "Max.n.SNP" = max(n.SNPs.vec),
                        "Avg.Mbp" = round(mean(width.all), 2), "Min.Mbp" = round(min(width.all), 2), "Max.Mbp" = round(max(width.all), 2))
    df.S2 <- rbind(df.S2, df.ij)
  }
}
write.csv(df.S2, paste0(dir.save, "/Table_S2.csv"), row.names = F)

# Clean-up
rm(list = ls()[ls() != "dir.save"])


# ---------------------------------------------------------------------------- #
# ----- Table S5
# ---------------------------------------------------------------------------- #
# package
library(data.table)
library(reshape2)
myFun.FormatImpRate <- function(x) {
  xf <- formatC(x, digits = 2, format = "f")
  if ( x == 0 ) { x.new <- paste0("±", xf) }
  if ( x > 0 ) { x.new <- paste0("+", xf) }
  if ( x < 0 ) { x.new <- xf }
  x.new <- paste0("(", x.new, ")")
  return(x.new)
}

# define object
tr.all <- c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3", "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")

# load phenotype data
df.ames.pheno <- fread("RESULT/1.1-MakeAmesPhenoData/AmesPheno.csv", data.table = F)
df.good.pheno <- fread("RESULT/1.2-MakeGbPhenoData/GbPheno.csv", data.table = F)

# Overlap
df.ames.ovlp <- read.csv("RESULT/AmesOverlapResult.csv")
df.good.ovlp <- read.csv("RESULT/GoodmanOverlapResult.csv")
id.ames.test <- df.ames.ovlp$Accession.Number[df.ames.ovlp$Overlap %in% c("No", "No (IBS < 0.8)")]
id.good.test <- df.good.ovlp$GBS.Sample[df.good.ovlp$Overlap %in% c("No", "No (IBS < 0.8)")]

# ##########################
# # replicate the old result
# tmp <- read.csv("/Users/ryokeitanaka/Desktop/TocoWork/Summary_of_1702_lines.csv")
# tmp.ames <- tmp[tmp$Panel == "Ames", ]
# id.ames.test <- tmp.ames$ID[tmp.ames$Overlap == "No"]
# tmp.good <- tmp[tmp$Panel == "282", ]
# id.good.test <- tmp.good$ID[tmp.good$Overlap == "No"]
# ##########################

# (1) From Ames to Goodman
df.pred.base <- read.csv("RESULT/99-SummaryFiles/BaselineModels_AmesTo282_BackTransPredVal.csv", check.names = FALSE)
df.pred.MK <- read.csv("RESULT/99-SummaryFiles/UseNamQtl_AmesTo282_BackTransPredVal.csv", check.names = FALSE)
checker <- all(colnames(df.pred.base) == colnames(df.pred.MK))
if ( checker ) { df.pred <- rbind(df.pred.base, df.pred.MK) }
df.PA <- data.frame("Trait" = df.pred$Trait, "Model" = df.pred$Model, "PA" = NA)
for ( i in 1:nrow(df.pred) ) {
  tr <- df.pred$Trait[i]
  mod <- df.pred$Model[i]
  pred <- unlist(df.pred[i, id.good.test])
  obs <- df.good.pheno[[tr]][match(id.good.test, df.good.pheno$ID)]
  df.PA$PA[i] <- cor(pred, obs, use = "p")
}
df.PA$IMP <- NA
for ( i in 1:nrow(df.PA) ) {
  tr <- df.PA$Trait[i]
  pa <- df.PA$PA[i]
  pa.base <- df.PA$PA[(df.PA$Trait == tr) & (df.PA$Model == "GBLUP")]
  df.PA$IMP[i] <- 100 * ((pa / pa.base) - 1)
}
df.PA$Trait <- factor(df.PA$Trait, levels = tr.all)
df.PA$Model <- factor(df.PA$Model, levels = c("GBLUP", "BayesB", "250K", "1M", "SI"))
df.PA.wide <- dcast(data = df.PA, formula = Model ~ Trait, value.var = "PA")
df.IMP.wide <- dcast(data = df.PA, formula = Model ~ Trait, value.var = "IMP")
mat.PA.wide <- df.PA.wide[, -1]; rownames(mat.PA.wide) <- df.PA.wide$Model
mat.IMP.wide <- df.IMP.wide[, -1]; rownames(mat.IMP.wide) <- df.IMP.wide$Model
mat.PA.wide <- cbind(mat.PA.wide, "Average" = apply(mat.PA.wide, 1, mean)) 
mat.IMP.wide <- cbind(mat.IMP.wide, "Average" = apply(mat.IMP.wide, 1, mean))
mat.PA.wide <- round(mat.PA.wide, 2)
mat.IMP.wide <- round(mat.IMP.wide, 2)
mat.S5.01 <- mat.PA.wide
for ( i in 1:nrow(mat.S5.01) ) {
  for ( j in 1:ncol(mat.S5.01) ) {
    mat.S5.01[i, j] <- paste(formatC(mat.PA.wide[i, j], digits = 2, format = "f"),
                             myFun.FormatImpRate(mat.IMP.wide[i, j]))    
  }
}

# (2) From Goodman to Ames
df.pred.base <- read.csv("RESULT/99-SummaryFiles/BaselineModels_282ToAmes_BackTransPredVal.csv", check.names = FALSE)
df.pred.MK <- read.csv("RESULT/99-SummaryFiles/UseNamQtl_282ToAmes_BackTransPredVal.csv", check.names = FALSE)
checker <- all(colnames(df.pred.base) == colnames(df.pred.MK))
if ( checker ) { df.pred <- rbind(df.pred.base, df.pred.MK) }
df.PA <- data.frame("Trait" = df.pred$Trait, "Model" = df.pred$Model, "PA" = NA)
for ( i in 1:nrow(df.pred) ) {
  tr <- df.pred$Trait[i]
  mod <- df.pred$Model[i]
  pred <- unlist(df.pred[i, id.ames.test])
  obs <- df.ames.pheno[[tr]][match(id.ames.test, df.ames.pheno$ID)]
  df.PA$PA[i] <- cor(pred, obs, use = "p")
}
df.PA$IMP <- NA
for ( i in 1:nrow(df.PA) ) {
  tr <- df.PA$Trait[i]
  pa <- df.PA$PA[i]
  pa.base <- df.PA$PA[(df.PA$Trait == tr) & (df.PA$Model == "GBLUP")]
  df.PA$IMP[i] <- 100 * ((pa / pa.base) - 1)
}
df.PA$Trait <- factor(df.PA$Trait, levels = tr.all)
df.PA$Model <- factor(df.PA$Model, levels = c("GBLUP", "BayesB", "250K", "1M", "SI"))
df.PA.wide <- dcast(data = df.PA, formula = Model ~ Trait, value.var = "PA")
df.IMP.wide <- dcast(data = df.PA, formula = Model ~ Trait, value.var = "IMP")
mat.PA.wide <- df.PA.wide[, -1]; rownames(mat.PA.wide) <- df.PA.wide$Model
mat.IMP.wide <- df.IMP.wide[, -1]; rownames(mat.IMP.wide) <- df.IMP.wide$Model
mat.PA.wide <- cbind(mat.PA.wide, "Average" = apply(mat.PA.wide, 1, mean)) 
mat.IMP.wide <- cbind(mat.IMP.wide, "Average" = apply(mat.IMP.wide, 1, mean))
mat.PA.wide <- round(mat.PA.wide, 2)
mat.IMP.wide <- round(mat.IMP.wide, 2)
mat.S5.02 <- mat.PA.wide
for ( i in 1:nrow(mat.S5.02) ) {
  for ( j in 1:ncol(mat.S5.02) ) {
    mat.S5.02[i, j] <- paste(formatC(mat.PA.wide[i, j], digits = 2, format = "f"), 
                             myFun.FormatImpRate(mat.IMP.wide[i, j]))    
  }
}

# (3) Within Ames
df.pred.base <- read.csv("RESULT/99-SummaryFiles/BaselineModels_WithinAmes_BackTransPredVal.csv", check.names = FALSE)
df.pred.MK <- read.csv("RESULT/99-SummaryFiles/UseNamQtl_WithinAmes_BackTransPredVal.csv", check.names = FALSE)
checker <- all(colnames(df.pred.base) == colnames(df.pred.MK))
if ( checker ) { df.pred <- rbind(df.pred.base, df.pred.MK) }
df.PA <- data.frame("Trait" = df.pred$Trait, "Model" = df.pred$Model, "Rep" = df.pred$Rep, "PA" = NA)
for ( i in 1:nrow(df.pred) ) {
  tr <- df.pred$Trait[i]
  mod <- df.pred$Model[i]
  rep <- df.pred$Rep[i]
  pred <- unlist(df.pred[i, 4:ncol(df.pred)])
  obs <- df.ames.pheno[[tr]][match(names(pred), df.ames.pheno$ID)]
  df.PA$PA[i] <- cor(pred, obs, use = "p")
}
df.PA$IMP <- NA
for ( i in 1:nrow(df.PA) ) {
  tr <- df.PA$Trait[i]
  rep <- df.PA$Rep[i]
  pa <- df.PA$PA[i]
  pa.base <- df.PA$PA[(df.PA$Trait == tr) & (df.PA$Rep == rep) & (df.PA$Model == "GBLUP")]
  df.PA$IMP[i] <- 100 * ((pa / pa.base) - 1)
}
df.PA$Trait <- factor(df.PA$Trait, levels = tr.all)
df.PA$Model <- factor(df.PA$Model, levels = c("GBLUP", "BayesB", "250K", "1M", "SI"))
df.PA.wide <- dcast(data = df.PA, formula = Model ~ Trait, value.var = "PA", fun.aggregate = mean)
df.IMP.wide <- dcast(data = df.PA, formula = Model ~ Trait, value.var = "IMP", fun.aggregate = mean)
mat.PA.wide <- df.PA.wide[, -1]; rownames(mat.PA.wide) <- df.PA.wide$Model
mat.IMP.wide <- df.IMP.wide[, -1]; rownames(mat.IMP.wide) <- df.IMP.wide$Model
mat.PA.wide <- cbind(mat.PA.wide, "Average" = apply(mat.PA.wide, 1, mean)) 
mat.IMP.wide <- cbind(mat.IMP.wide, "Average" = apply(mat.IMP.wide, 1, mean))
mat.PA.wide <- round(mat.PA.wide, 2)
mat.IMP.wide <- round(mat.IMP.wide, 2)
mat.S5.03 <- mat.PA.wide
for ( i in 1:nrow(mat.S5.03) ) {
  for ( j in 1:ncol(mat.S5.03) ) {
    mat.S5.03[i, j] <- paste(formatC(mat.PA.wide[i, j], digits = 2, format = "f"),
                             myFun.FormatImpRate(mat.IMP.wide[i, j]))    
  }
}

# (4) Within Goodman
df.pred.base <- read.csv("RESULT/99-SummaryFiles/BaselineModels_Within282_BackTransPredVal.csv", check.names = FALSE)
df.pred.MK <- read.csv("RESULT/99-SummaryFiles/UseNamQtl_Within282_BackTransPredVal.csv", check.names = FALSE)
checker <- all(colnames(df.pred.base) == colnames(df.pred.MK))
if ( checker ) { df.pred <- rbind(df.pred.base, df.pred.MK) }
df.PA <- data.frame("Trait" = df.pred$Trait, "Model" = df.pred$Model, "Rep" = df.pred$Rep, "PA" = NA)
for ( i in 1:nrow(df.pred) ) {
  tr <- df.pred$Trait[i]
  mod <- df.pred$Model[i]
  rep <- df.pred$Rep[i]
  pred <- unlist(df.pred[i, 4:ncol(df.pred)])
  obs <- df.good.pheno[[tr]][match(names(pred), df.good.pheno$ID)]
  df.PA$PA[i] <- cor(pred, obs, use = "p")
}
df.PA$IMP <- NA
for ( i in 1:nrow(df.PA) ) {
  tr <- df.PA$Trait[i]
  rep <- df.PA$Rep[i]
  pa <- df.PA$PA[i]
  pa.base <- df.PA$PA[(df.PA$Trait == tr) & (df.PA$Rep == rep) & (df.PA$Model == "GBLUP")]
  df.PA$IMP[i] <- 100 * ((pa / pa.base) - 1)
}
df.PA$Trait <- factor(df.PA$Trait, levels = tr.all)
df.PA$Model <- factor(df.PA$Model, levels = c("GBLUP", "BayesB", "250K", "1M", "SI"))
df.PA.wide <- dcast(data = df.PA, formula = Model ~ Trait, value.var = "PA", fun.aggregate = mean)
df.IMP.wide <- dcast(data = df.PA, formula = Model ~ Trait, value.var = "IMP", fun.aggregate = mean)
mat.PA.wide <- df.PA.wide[, -1]; rownames(mat.PA.wide) <- df.PA.wide$Model
mat.IMP.wide <- df.IMP.wide[, -1]; rownames(mat.IMP.wide) <- df.IMP.wide$Model
mat.PA.wide <- cbind(mat.PA.wide, "Average" = apply(mat.PA.wide, 1, mean)) 
mat.IMP.wide <- cbind(mat.IMP.wide, "Average" = apply(mat.IMP.wide, 1, mean))
mat.PA.wide <- round(mat.PA.wide, 2)
mat.IMP.wide <- round(mat.IMP.wide, 2)
mat.S5.04 <- mat.PA.wide
for ( i in 1:nrow(mat.S5.04) ) {
  for ( j in 1:ncol(mat.S5.04) ) {
    mat.S5.04[i, j] <- paste(formatC(mat.PA.wide[i, j], digits = 2, format = "f"),
                             myFun.FormatImpRate(mat.IMP.wide[i, j]))    
  }
}

# (5) merge the four result 
mat.S5 <- rbind(cbind("Scenario" = "From Ames to Goodman", "Model" = rownames(mat.S5.01), mat.S5.01),
                cbind("Scenario" = "From Goodman to Ames", "Model" = rownames(mat.S5.02), mat.S5.02),
                cbind("Scenario" = "Within Ames", "Model" = rownames(mat.S5.03), mat.S5.03),
                cbind("Scenario" = "Within Goodman", "Model" = rownames(mat.S5.04), mat.S5.04))
write.csv(mat.S5, file = paste0(dir.save, "/Table_S5_updated.csv"), row.names = F)

# Clean-up
rm(list = ls()[ls() != "dir.save"])

# ---------------------------------------------------------------------------- #
# ----- Table S6
# ---------------------------------------------------------------------------- #
# package
library(data.table)

# object
tr.all <- c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3", "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")
model.all <- c("250K", "1M", "SI")

# load NAM PVE 
NamData.Full <- read.csv("RAWDATA/NAM_QTL/qtl.data.from.di.csv")
colnames(NamData.Full)[1] <- "QTL.ID"
NamData <- NamData.Full[NamData.Full$Trait != "PC-8", ]
NamData$NAME3 <- paste0("QTL", NamData$QTL.ID)

# load phenotype data
df.ames.pheno <- fread("RESULT/1.1-MakeAmesPhenoData/AmesPheno.csv", data.table = F)
df.good.pheno <- fread("RESULT/1.2-MakeGbPhenoData/GbPheno.csv", data.table = F)

# Overlap
df.ames.ovlp <- read.csv("RESULT/AmesOverlapResult.csv")
df.good.ovlp <- read.csv("RESULT/GoodmanOverlapResult.csv")
id.ames.test <- df.ames.ovlp$Accession.Number[df.ames.ovlp$Overlap %in% c("No", "No (IBS < 0.8)")]
id.good.test <- df.good.ovlp$GBS.Sample[df.good.ovlp$Overlap %in% c("No", "No (IBS < 0.8)")]

# # replicate the old result
# tmp <- read.csv("/Users/ryokeitanaka/Desktop/TocoWork/Summary_of_1702_lines.csv")
# tmp.ames <- tmp[tmp$Panel == "Ames", ]
# id.ames.test <- tmp.ames$ID[tmp.ames$Overlap == "No"]
# tmp.good <- tmp[tmp$Panel == "282", ]
# id.good.test <- tmp.good$ID[tmp.good$Overlap == "No"]

# (1) From Ames to Goodman
df.pred.all <- read.csv("RESULT/99-SummaryFiles/PostHoc_UseNamQtl_AmesTo282_EvalEachKernel.csv", check.names = F)
df.pred <- df.pred.all[df.pred.all$QTL != "UseAll", ]
df.pcor.Ames.to.Good <- NULL
tr.all <- c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3", "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")
model.all <- c("250K", "1M", "SI")
for ( i in 1:length(tr.all) ) {
  for ( j in 1:length(model.all) ) {
    tr <- tr.all[i]
    win <- model.all[j]
    df.pred.sub <- df.pred[(df.pred$Trait == tr) & (df.pred$Window == win), ]
    mat.pred <- t(as.matrix(df.pred.sub[, id.good.test]))
    colnames(mat.pred) <- df.pred.sub$QTL
    obs.vec <- df.good.pheno[[tr]][match(id.good.test, df.good.pheno$ID)]
    df.fit.lm <- data.frame("y" = obs.vec, mat.pred)
    res.aov <- aov(y ~ ., df.fit.lm)
    res.type3 <- car::Anova(res.aov, type = "III")
    m <- match(colnames(df.fit.lm)[-1], rownames(res.type3))
    pval.vec <- res.type3$`Pr(>F)`[m]
    names(pval.vec) <- rownames(res.type3)[m]
    res.pcor <- ppcor::pcor(df.fit.lm[!is.na(df.fit.lm$y), ], method = "pearson")
    pcor.vec <- res.pcor$estimate[-1, 1]
    df.sub <- data.frame(df.pred.sub[, 1:3], "pcor" = pcor.vec)
    df.pcor.Ames.to.Good <- rbind(df.pcor.Ames.to.Good, df.sub)
  }
}
df.pcor.Ames.to.Good <- df.pcor.Ames.to.Good[df.pcor.Ames.to.Good$QTL != "NonQtl", ] # remove non-QTL
for ( i in 1:nrow(df.pcor.Ames.to.Good) ) {
  pve <- NamData[(NamData$NAME3 == df.pcor.Ames.to.Good$QTL[i]) & (NamData$Trait == df.pcor.Ames.to.Good$Trait[i]), "PVE"]
  df.pcor.Ames.to.Good$PVE[i] <- pve # get PVE
}

# (2) From Ames to Goodman
df.pred.all <- read.csv("RESULT/99-SummaryFiles/PostHoc_UseNamQtl_282ToAmes_EvalEachKernel.csv", check.names = F)
df.pred <- df.pred.all[df.pred.all$QTL != "UseAll", ]
df.pcor.Good.to.Ames <- NULL
tr.all <- c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3", "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")
model.all <- c("250K", "1M", "SI")
for ( i in 1:length(tr.all) ) {
  for ( j in 1:length(model.all) ) {
    tr <- tr.all[i]
    win <- model.all[j]
    df.pred.sub <- df.pred[(df.pred$Trait == tr) & (df.pred$Window == win), ]
    mat.pred <- t(as.matrix(df.pred.sub[, id.ames.test]))
    colnames(mat.pred) <- df.pred.sub$QTL
    obs.vec <- df.ames.pheno[[tr]][match(id.ames.test, df.ames.pheno$ID)]
    df.fit.lm <- data.frame("y" = obs.vec, mat.pred)
    res.aov <- aov(y ~ ., df.fit.lm)
    res.type3 <- car::Anova(res.aov, type = "III")
    m <- match(colnames(df.fit.lm)[-1], rownames(res.type3))
    pval.vec <- res.type3$`Pr(>F)`[m]
    names(pval.vec) <- rownames(res.type3)[m]
    res.pcor <- ppcor::pcor(df.fit.lm[!is.na(df.fit.lm$y), ], method = "pearson")
    pcor.vec <- res.pcor$estimate[-1, 1]
    df.sub <- data.frame(df.pred.sub[, 1:3], "pcor" = pcor.vec)
    df.pcor.Good.to.Ames <- rbind(df.pcor.Good.to.Ames, df.sub)
  }
}
df.pcor.Good.to.Ames <- df.pcor.Good.to.Ames[df.pcor.Good.to.Ames$QTL != "NonQtl", ] # remove non-QTL
for ( i in 1:nrow(df.pcor.Good.to.Ames) ) {
  pve <- NamData[(NamData$NAME3 == df.pcor.Good.to.Ames$QTL[i]) & (NamData$Trait == df.pcor.Good.to.Ames$Trait[i]), "PVE"]
  df.pcor.Good.to.Ames$PVE[i] <- pve # get PVE
}

# make table S6
df.S6 <- expand.grid("Window" = c("250K", "1M", "SI"), 
                     "Scenario" = c("From Ames to Goodman", "From Goodman to Ames"), 
                     "PVE" = c("Low", "High"),
                     stringsAsFactors = F)[, c(2, 1, 3)]
df.S6$P.value <- df.S6$Correlation <- NA
for ( i in 1:nrow(df.S6) ) {
  sc <- df.S6$Scenario[i]
  win <- df.S6$Window[i]
  pve <- df.S6$PVE[i]
  if ( sc == "From Ames to Goodman" ) {
    df.sub <- df.pcor.Ames.to.Good[df.pcor.Ames.to.Good$Window == win, ]  
  } else if ( sc == "From Goodman to Ames" ) {
    df.sub <- df.pcor.Good.to.Ames[df.pcor.Good.to.Ames$Window == win, ]
  }
  if ( pve == "Low" ) {
    df.sub <- df.sub[df.sub$PVE < 0.05, ]
  } else if ( pve == "High" ) {
    df.sub <- df.sub[df.sub$PVE >= 0.05, ]
  }
  df.S6$Correlation[i] <- cor(df.sub$pcor, df.sub$PVE)
  df.S6$P.value[i] <- cor.test(df.sub$pcor, df.sub$PVE)$p.value
}
df.S6$Correlation <- round(df.S6$Correlation, 2)

# write
write.csv(df.S6, file = paste0(dir.save, "/Table_S6_updated.csv"), row.names = F)

# Clean-up
rm(list = ls()[ls() != "dir.save"])


# ---------------------------------------------------------------------------- #
# ----- Table S7
# ---------------------------------------------------------------------------- #
# package
library(data.table)
library(reshape2)
myFun.FormatImpRate <- function(x) {
  xf <- formatC(x, digits = 2, format = "f")
  if ( x == 0 ) { x.new <- paste0("±", xf) }
  if ( x > 0 ) { x.new <- paste0("+", xf) }
  if ( x < 0 ) { x.new <- xf }
  x.new <- paste0("(", x.new, ")")
  return(x.new)
}

# define object
tr.all <- c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3", "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")

# load phenotype data
df.ames.pheno <- fread("RESULT/1.1-MakeAmesPhenoData/AmesPheno.csv", data.table = F)

#
df.pred <- read.csv("RESULT/99-SummaryFiles/TrsctPred_WithinAmes_BackTransPredVal.csv", check.names = FALSE)
df.PA <- data.frame("Trait" = df.pred$Trait, "Model" = df.pred$Model, "Rep" = df.pred$Rep, "PA" = NA)
for ( i in 1:nrow(df.pred) ) {
  tr <- df.pred$Trait[i]
  mod <- df.pred$Model[i]
  rep <- df.pred$Rep[i]
  pred <- unlist(df.pred[i, 4:ncol(df.pred)])
  obs <- df.ames.pheno[[tr]][match(names(pred), df.ames.pheno$ID)]
  df.PA$PA[i] <- cor(pred, obs, use = "p")
}
df.PA$IMP <- NA
for ( i in 1:nrow(df.PA) ) {
  tr <- df.PA$Trait[i]
  rep <- df.PA$Rep[i]
  pa <- df.PA$PA[i]
  pa.base <- df.PA$PA[(df.PA$Trait == tr) & (df.PA$Rep == rep) & (df.PA$Model == "GBLUP")]
  df.PA$IMP[i] <- 100 * ((pa / pa.base) - 1)
}
df.PA$Trait <- factor(df.PA$Trait, levels = tr.all)
df.PA$Model <- factor(df.PA$Model, levels = c("GBLUP", "TBLUP", "TBLUP_sub", "GBLUP+TBLUP", "GBLUP+TBLUP_sub", "MultiNamLargeGenes"))
df.PA.wide <- dcast(data = df.PA, formula = Model ~ Trait, value.var = "PA", fun.aggregate = mean)
df.IMP.wide <- dcast(data = df.PA, formula = Model ~ Trait, value.var = "IMP", fun.aggregate = mean)
mat.PA.wide <- df.PA.wide[, -1]; rownames(mat.PA.wide) <- df.PA.wide$Model
mat.IMP.wide <- df.IMP.wide[, -1]; rownames(mat.IMP.wide) <- df.IMP.wide$Model
mat.PA.wide <- cbind(mat.PA.wide, "Average" = apply(mat.PA.wide, 1, mean)) 
mat.IMP.wide <- cbind(mat.IMP.wide, "Average" = apply(mat.IMP.wide, 1, mean))
mat.PA.wide <- round(mat.PA.wide, 2)
mat.IMP.wide <- round(mat.IMP.wide, 2)
mat.S7 <- mat.PA.wide
for ( i in 1:nrow(mat.S7) ) {
  for ( j in 1:ncol(mat.S7) ) {
    mat.S7[i, j] <- paste(formatC(mat.PA.wide[i, j], digits = 2, format = "f"),
                          myFun.FormatImpRate(mat.IMP.wide[i, j]))    
  }
}
write.csv(mat.S7, file = paste0(dir.save, "/Table_S7.csv"), row.names = T)

