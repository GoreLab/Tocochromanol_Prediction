# package
library(data.table)
library(ggplot2)
library(extrafont)
library(inlmisc)
library(plyr)
library(reshape2)
library(car)
library(ggrepel)
library(scales)
library(ggpubr)
library(RColorBrewer)
library(ggtext)
library(readxl)

# function to be used
see <- function(x){x[1:5, 1:5]}
MyFun.load <- function(file, ...){
  dat <- fread(file, data.table = F)
  tf.ames <- all(AmesPheno$ID %in% colnames(dat))
  tf.across <- !("Rep" %in% colnames(dat))
  if ( tf.ames == T & tf.across == T ) { s <- "From 282 to Ames" }
  if ( tf.ames == F & tf.across == T ) { s <- "From Ames to 282" }
  if ( tf.ames == T & tf.across == F ) { s <- "Within Ames" }
  if ( tf.ames == F & tf.across == F ) { s <- "Within 282" }
  dat$Scenario <- s
  return(dat)   
}
myFun.Calc.PA <- function(df.pred, df.obs, pred.id.all) {
  # check the column names of the data
  checker <- !all(df.obs$ID %in% colnames(df.pred))
  if ( checker ) {
    out <- NULL
  } else {
    # split the data frame of the predicted values
    m <- match(df.obs$ID, colnames(df.pred))
    mat.pred.vals <- as.matrix(df.pred[, m])
    df.pred.info <- df.pred[, !(colnames(df.pred) %in% pred.id.all)]
    
    # calculate correlation
    cor.vec <- rep(NA, times = nrow(df.pred))
    for ( i in 1:nrow(df.pred) ) {
      tr <- df.pred.info$Trait[i]
      pred <- mat.pred.vals[i, ]
      obs <- df.obs[[tr]]
      cor.vec[i] <- cor(obs, pred, use = "p")
    }
    
    # add column of the correlation (predictive ability)
    df.pred.info$PA <- cor.vec
    out <- df.pred.info
  }
  
  # return
  return(out)
} # revised
myFun.Calc.PPA <- function(df.pred, df.obs) {
  # check the column names of the data
  checker <- !all(df.obs$ID %in% colnames(df.pred))
  if ( checker ) {
    out <- NULL
  } else {
    win.all <- unique(df.pred$Window)
    df.res.all <- NULL
    for (w in 1:length(win.all)) {
      win <- win.all[w]
      tf <- df.pred$Window == win
      df <- df.pred[tf, ]
      df <- df[df$QTL != "UseAll", ]
      df.res <- NULL
      for ( i in 1:length(unique(df$Trait)) ) {
        tr <- unique(df$Trait)[i]
        df.i <- df[df$Trait == tr, ]
        pred.mat <- t(df.i[, df.obs$ID])
        colnames(pred.mat) <- df.i$QTL
        obs.vec <- df.obs[, tr]
        df.fit.lm <- data.frame("y" = obs.vec, 
                                pred.mat)
        res.aov <- aov(y ~ ., df.fit.lm)
        res.type3 <- car::Anova(res.aov, type = "III")
        m <- match(colnames(df.fit.lm)[-1], rownames(res.type3))
        pval.vec <- res.type3$`Pr(>F)`[m]
        names(pval.vec) <- rownames(res.type3)[m]
        res.pcor <- ppcor::pcor(df.fit.lm[!is.na(df.fit.lm$y), ], method = "pearson")
        pcor.vec <- res.pcor$estimate[-1, 1]
        cheker.i <- all(df.i$QTL == names(pcor.vec)) & all(df.i$QTL == names(pval.vec))
        if ( cheker.i ) {
          df.res.i <- data.frame(df.i[, c("Scenario", "Trait", "Window", "QTL")],
                                 "Partial.Cor" = pcor.vec,
                                 "P.value" = pval.vec)
        } else {
          df.res.i <- NULL
        }
        df.res <- rbind.data.frame(df.res, df.res.i)
      }
      df.res.all <- rbind.data.frame(df.res.all, df.res)
    }
    out <- df.res.all
  }
  # return
  return(out)
}
myFun.Calc.SD <- function(df.PA) {
  if ( !("Rep" %in% colnames(df.PA)) ) {
    df.PA$SD <- 0
    out <- df.PA
  } else {
    factor.names <- colnames(df.PA)[!(colnames(df.PA) %in% c("PA", "Rep"))]
    formula <- paste0("PA ~ ", paste0(factor.names, collapse = " + "))
    df.avg <- aggregate(eval(parse(text=formula)), df.PA, mean)
    df.sd <- aggregate(eval(parse(text=formula)), df.PA, sd)
    colnames(df.sd)[colnames(df.sd) == "PA"] <- "SD"
    out <- merge.data.frame(df.avg, df.sd)
  }
  return(out)
}
myFun.Attach.SD.bar <- function(df.fig) {
  df.fig$SD.bar.min <- df.fig$PA - df.fig$SD
  df.fig$SD.bar.max <- df.fig$PA + df.fig$SD
  return(df.fig)
}
myFun.ChangeTraitNames <- function(x) {
  x <- as.character(x)
  x[x == "a.T"] <- "\u03B1T"
  x[x == "g.T"] <- "\u03B3T"
  x[x == "d.T"] <- "\u03B4T"
  x[x == "a.T3"] <- "\u03B1T3"
  x[x == "g.T3"] <- "\u03B3T3"
  x[x == "d.T3"] <- "\u03B4T3"
  x[x == "Total.Tocopherols"] <- "\u03A3T"
  x[x == "Total.Tocotrienols"] <- "\u03A3T3"
  x[x == "Total.Tocochromanols"] <- "\u03A3TT3"
  x <- factor(x, levels = c("\u03B1T", "\u03B4T", "\u03B3T",
                            "\u03B1T3", "\u03B4T3", "\u03B3T3",
                            "\u03A3T", "\u03A3T3", "\u03A3TT3"))
  return(x)
}
myFun.ChangeModeltNames.NamQtl <- function(x) {
  x[x == "250K"] <- "MK-GBLUP\n(± 250 kb)"
  x[x == "1M"] <- "MK-GBLUP\n(± 1 Mb)"
  x[x == "SI"] <- "MK-GBLUP\n(SI)"
  x <- factor(x, levels = c("GBLUP", 
                            "BayesB", 
                            "MK-GBLUP\n(± 250 kb)", 
                            "MK-GBLUP\n(± 1 Mb)", 
                            "MK-GBLUP\n(SI)"))
  return(x)
}
myFun.Calc.Imp <- function(df, base.model) {
  imp <- rep(NA, nrow(df))
  for ( i in 1:nrow(df) ) {
    tr.i <- df$Trait[i]
    sc.i <- df$Scenario[i]
    rep.i <- df$Rep[i]
    pa.i <- df$PA[i]
    if ( is.na(rep.i) ) {
      tf <- (df$Trait == tr.i) & (df$Scenario == sc.i) & (df$Model == base.model)
      base.r <- df$PA[tf]
    } else {
      tf <- (df$Trait == tr.i) & (df$Scenario == sc.i) & (df$Rep == rep.i) & (df$Model == base.model)
      base.r <- df$PA[tf]
    }
    imp[i] <- 100 * (pa.i / base.r - 1) 
  }
  return(imp)
}
myFun.ChangeModeltNames.Expr <- function(x) {
  x[x == "GBLUP"] <- "GBLUP\n(Baseline)"
  x[x == "TBLUP"] <- "TRM.all"
  x[x == "TBLUP_sub"] <- "TRM.cand"
  x[x == "GBLUP+TBLUP"] <- "GRM+TRM.all"
  x[x == "GBLUP+TBLUP_sub"] <- "GRM+TRM.cand"
  x[x == "MultiNamLargeGenes"] <- "Multi-trait GBLUP"
  x <- factor(x, levels = c("GBLUP\n(Baseline)", "TRM.all", 
                            "TRM.cand", "GRM+TRM.all",
                            "GRM+TRM.cand", "Multi-trait GBLUP"))
  return(x)
}

# mkdir
dir.save <- "RESULT/Tables_and_Figures"
dir.create(dir.save, recursive = T)

# load phenotype data
f.ames <- "RESULT/1.1-MakeAmesPhenoData/AmesPheno.csv" # Ames Phenotype data
f.gb <- "RESULT/1.2-MakeGbPhenoData/GbPheno.csv" # Gb Phenotype data
AmesPheno <- fread(f.ames, data.table = F)
GbPheno <- fread(f.gb, data.table = F)
trait.all <- c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3",
               "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")

# load overlap & create non-overlapped phenotype data
df.ames.ovlp <- read.csv("RESULT/AmesOverlapResult.csv")
df.good.ovlp <- read.csv("RESULT/GoodmanOverlapResult.csv")
id.ames.test <- df.ames.ovlp$Accession.Number[df.ames.ovlp$Overlap %in% c("No", "No (IBS < 0.8)")]
id.good.test <- df.good.ovlp$GBS.Sample[df.good.ovlp$Overlap %in% c("No", "No (IBS < 0.8)")]
GbPheno.unique <- GbPheno[GbPheno$ID %in% id.good.test, ]
AmesPheno.unique <- AmesPheno[AmesPheno$ID %in% id.ames.test, ]

# setup for figure
cols <- GetColors(n = 13)
cols.QTL <- c("gray", "black", cols[10:12])

# ---------------------------------------------------------------------------- #
# ----- Figure: PA on the MK-GBLUP using NAM-QTL 
# ---------------------------------------------------------------------------- #
list.all <- list(
  BaseL.Ato2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_AmesTo282_BackTransPredVal.csv"), 
                             df.obs = GbPheno.unique, pred.id.all = GbPheno$ID),
  BaseL.2toA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_282ToAmes_BackTransPredVal.csv"), 
                             df.obs = AmesPheno.unique, pred.id.all = AmesPheno$ID),
  BaseL.AtoA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_WithinAmes_BackTransPredVal.csv"), 
                             df.obs = AmesPheno, pred.id.all = AmesPheno$ID),
  BaseL.2to2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_Within282_BackTransPredVal.csv"), 
                             df.obs = GbPheno, pred.id.all = GbPheno$ID),
  UseNQ.Ato2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_AmesTo282_BackTransPredVal.csv"), 
                             df.obs = GbPheno.unique, pred.id.all = GbPheno$ID),
  UseNQ.2toA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_282ToAmes_BackTransPredVal.csv"), 
                             df.obs = AmesPheno.unique, pred.id.all = AmesPheno$ID),
  UseNQ.AtoA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_WithinAmes_BackTransPredVal.csv"), 
                             df.obs = AmesPheno, pred.id.all = AmesPheno$ID),
  UseNQ.2to2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_Within282_BackTransPredVal.csv"), 
                             df.obs = GbPheno, pred.id.all = GbPheno$ID)
)
list.all.SD <- lapply(list.all, myFun.Calc.SD)
df.fig <- ldply(list.all.SD, data.frame)
df.fig$Model <- myFun.ChangeModeltNames.NamQtl(df.fig$Model)
df.fig$Trait <- myFun.ChangeTraitNames(df.fig$Trait)
df.fig$Scenario[df.fig$Scenario == "Within Ames"] <- "'Within Ames ('*italic('n')*' = 1,462)'"
df.fig$Scenario[df.fig$Scenario == "From Ames to 282"] <- "'From Ames ('*italic('n')*' = 1,462) to Goodman ('*italic('n')*' = 73)'"
df.fig$Scenario[df.fig$Scenario == "From 282 to Ames"] <- "'From Goodman ('*italic('n')*' = 242) to Ames ('*italic('n')*' = 1,283)'"
df.fig$Scenario[df.fig$Scenario == "Within 282"] <- "'Within Goodman ('*italic('n')*' = 242)'"
df.fig$Scenario <- factor(df.fig$Scenario, 
                          levels = c("'Within Ames ('*italic('n')*' = 1,462)'", 
                                     "'From Ames ('*italic('n')*' = 1,462) to Goodman ('*italic('n')*' = 73)'", 
                                     "'From Goodman ('*italic('n')*' = 242) to Ames ('*italic('n')*' = 1,283)'", 
                                     "'Within Goodman ('*italic('n')*' = 242)'"))
df.fig <- myFun.Attach.SD.bar(df.fig)
p <- ggplot(df.fig, aes(x = Trait, y = PA, fill = Model))
p <- p + facet_wrap(~ Scenario, labeller = label_parsed)
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
                           ymax = SD.bar.max), 
                       position = position_dodge(.9), width=.5)
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times"))
p <- p + xlab("Phenotype")
p <- p + ylab("Predictive ability")
p <- p + scale_fill_manual(values = cols.QTL)
p <- p + theme(legend.key.height= unit(1, 'cm'))
ggplot2::ggsave(filename = paste0(dir.save, "/Fig1a_small_TwoSideBar_Updated.eps"), 
                plot = p, 
                device = cairo_ps, 
                dpi = 300, 
                width = 12 * 0.7,
                height = 6 * 0.6)

# ---------------------------------------------------------------------------- #
# ----- Figure: improvement of PA on the MK-GBLUP using NAM-QTL 
# ---------------------------------------------------------------------------- #
list.all <- list(
  BaseL.Ato2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_AmesTo282_BackTransPredVal.csv"), 
                             df.obs = GbPheno.unique, pred.id.all = GbPheno$ID),
  BaseL.2toA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_282ToAmes_BackTransPredVal.csv"), 
                             df.obs = AmesPheno.unique, pred.id.all = AmesPheno$ID),
  BaseL.AtoA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_WithinAmes_BackTransPredVal.csv"), 
                             df.obs = AmesPheno, pred.id.all = AmesPheno$ID),
  BaseL.2to2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_Within282_BackTransPredVal.csv"), 
                             df.obs = GbPheno, pred.id.all = GbPheno$ID),
  UseNQ.Ato2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_AmesTo282_BackTransPredVal.csv"), 
                             df.obs = GbPheno.unique, pred.id.all = GbPheno$ID),
  UseNQ.2toA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_282ToAmes_BackTransPredVal.csv"), 
                             df.obs = AmesPheno.unique, pred.id.all = AmesPheno$ID),
  UseNQ.AtoA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_WithinAmes_BackTransPredVal.csv"), 
                             df.obs = AmesPheno, pred.id.all = AmesPheno$ID),
  UseNQ.2to2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_Within282_BackTransPredVal.csv"), 
                             df.obs = GbPheno, pred.id.all = GbPheno$ID)
)
df.fig <- ldply(list.all, data.frame)
df.fig$PA <- myFun.Calc.Imp(df = df.fig, base.model = "GBLUP")
df.fig <- myFun.Calc.SD(df.fig)
df.fig <- df.fig[df.fig$Model != "GBLUP", ]
df.fig$Model <- myFun.ChangeModeltNames.NamQtl(df.fig$Model)
df.fig$Trait <- myFun.ChangeTraitNames(df.fig$Trait)
df.fig$Scenario[df.fig$Scenario == "Within Ames"] <- "'Within Ames ('*italic('n')*' = 1,462)'"
df.fig$Scenario[df.fig$Scenario == "From Ames to 282"] <- "'From Ames ('*italic('n')*' = 1,462) to Goodman ('*italic('n')*' = 73)'"
df.fig$Scenario[df.fig$Scenario == "From 282 to Ames"] <- "'From Goodman ('*italic('n')*' = 242) to Ames ('*italic('n')*' = 1,283)'"
df.fig$Scenario[df.fig$Scenario == "Within 282"] <- "'Within Goodman ('*italic('n')*' = 242)'"
df.fig$Scenario <- factor(df.fig$Scenario, 
                          levels = c("'Within Ames ('*italic('n')*' = 1,462)'", 
                                     "'From Ames ('*italic('n')*' = 1,462) to Goodman ('*italic('n')*' = 73)'", 
                                     "'From Goodman ('*italic('n')*' = 242) to Ames ('*italic('n')*' = 1,283)'", 
                                     "'Within Goodman ('*italic('n')*' = 242)'"))
df.fig <- myFun.Attach.SD.bar(df.fig)
p <- ggplot(df.fig, aes(x = Trait, y = PA, fill = Model))
p <- p + facet_wrap(~ Scenario, labeller = label_parsed)
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
                           ymax = SD.bar.max), 
                       position = position_dodge(.9), width=.5)
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times"))
p <- p + xlab("Phenotype")
p <- p + ylab("Improvement rate (%) relative to GBLUP")
p <- p + scale_fill_manual(values = cols.QTL[-1])
p <- p + theme(legend.key.height= unit(1, 'cm'))
ggplot2::ggsave(filename = paste0(dir.save, "/Fig1b_small_TwoSideBar_Updated.eps"), 
                plot = p, 
                device = cairo_ps, 
                dpi = 300, 
                width = 12 * 0.7,
                height = 6 * 0.6)

# ---------------------------------------------------------------------------- #
# ----- Figure: PA from the Transcriptome based prediction
# ---------------------------------------------------------------------------- #
pred <- read.csv("RESULT/99-SummaryFiles/TrsctPred_WithinAmes_BackTransPredVal.csv")
df.cor <- pred[, 1:3]
pred.mat <- as.matrix(pred[, 4:ncol(pred)])
m <- match(colnames(pred.mat), AmesPheno$ID)
AmesPheno.545 <- AmesPheno[m, ]
df.cor$PA <- NA
for ( i in 1:nrow(df.cor) ) {
  tr.i <- df.cor$Trait[i]
  obs <- AmesPheno.545[[tr.i]]
  pred <- pred.mat[i, ]
  df.cor$PA[i] <- cor(obs, pred, use = "p")
}
df.fig <- myFun.Calc.SD(df.cor)
df.fig <- myFun.Attach.SD.bar(df.fig)
df.fig$Trait <- myFun.ChangeTraitNames(df.fig$Trait)
df.fig$Model <- myFun.ChangeModeltNames.Expr(df.fig$Model)
df.fig <- df.fig[df.fig$Model != "GBLUP with covariates\n(Large-eff. genes)", ]

#
cols <- c("gray", brewer.pal(n = 5, name = "YlOrRd"))
p <- ggplot(df.fig, aes(x = Trait, y = PA, fill = Model))
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
                           ymax = SD.bar.max), 
                       position = position_dodge(.9), width=.5)
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times"))
p <- p + xlab("Phenotype")
p <- p + ylab("Predictive ability")
p <- p + theme(legend.key.height= unit(1, 'cm'))
p <- p + scale_fill_manual(values = cols)
ggplot2::ggsave(filename = paste0(dir.save, "/Fig2a_small_TwoSideBar.eps"), 
                plot = p, 
                device = cairo_ps, 
                dpi = 300, 
                width = 12 * 0.7,
                height = 5 * 0.7)

# ---------------------------------------------------------------------------- #
# ----- Figure: Imp.Rate from the Transcriptome based prediction
# ---------------------------------------------------------------------------- #
pred <- read.csv("RESULT/99-SummaryFiles/TrsctPred_WithinAmes_BackTransPredVal.csv")
df <- pred[, 1:3]
pred.mat <- as.matrix(pred[, 4:ncol(pred)])
m <- match(colnames(pred.mat), AmesPheno$ID)
AmesPheno.545 <- AmesPheno[m, ]
df$PA <- NA
for ( i in 1:nrow(df) ) {
  tr.i <- df$Trait[i]
  obs <- AmesPheno.545[[tr.i]]
  pred <- pred.mat[i, ]
  df$PA[i] <- cor(obs, pred, use = "p")
}
df$PA.IMP <- NA
for ( i in 1:nrow(df) ) {
  tr.i <- df$Trait[i]
  rep.i <- df$Rep[i]
  pa.i <- df$PA[i]
  tf <- (df$Trait == tr.i) & (df$Rep == rep.i) & (df$Model == "GBLUP")
  base.r <- df$PA[tf]
  df$PA.IMP[i] <- 100 * (pa.i / base.r - 1)
}
df.cor <- df[, c("Trait", "Model", "Rep", "PA.IMP")]
colnames(df.cor)[colnames(df.cor) == "PA.IMP"] <- "PA"
df.fig <- myFun.Calc.SD(df.cor)
df.fig <- myFun.Attach.SD.bar(df.fig)
df.fig$Trait <- myFun.ChangeTraitNames(df.fig$Trait)
df.fig$Model <- myFun.ChangeModeltNames.Expr(df.fig$Model)
df.fig <- df.fig[df.fig$Model != "GBLUP\n(Baseline)", ]
cols <- brewer.pal(n = 5, name = "YlOrRd")
p <- ggplot(df.fig, aes(x = Trait, y = PA, fill = Model))
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
                           ymax = SD.bar.max), 
                       position = position_dodge(.9), width=.5)
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times"))
p <- p + xlab("Phenotype")
p <- p + ylab("Improvement Rate (%) relative to GBLUP")
p <- p + scale_fill_manual(values = cols)
p <- p + theme(legend.key.height= unit(1, 'cm'))
ggplot2::ggsave(filename = paste0(dir.save, "/Fig2b_small_TwoSideBar_Updated.eps"), 
                plot = p, 
                device = cairo_ps, 
                dpi = 300, 
                width = 12 * 0.7,
                height = 5 * 0.7)

# ---------------------------------------------------------------------------- #
# ----- Figure: PCA biplot
# ---------------------------------------------------------------------------- #
# list of the 126 genes
list.cand <- as.data.frame(read_xlsx("RAWDATA/tocochromanol_all_candidate_genes_combined_DW_20190624_with_two_genes.xlsx"))
id.126 <- unique(list.cand$RefGen_v4.Gene.ID)
df.cand <- data.frame("id" = id.126)
df.cand$name <- list.cand$Gene.Name[match(df.cand$id, list.cand$RefGen_v4.Gene.ID)]

# show the two arodeH2
df.cand[df.cand$id %in% c("Zm00001d014734", "Zm00001d014737"), ]

# list of the 13 genes
df.gene.13 <- as.data.frame(read_xlsx("RAWDATA/CandidateGeneNames.xlsx"))
df.gene.13$Gene <- paste0("italic(", df.gene.13$Gene, ")")
df.gene.13$Gene[df.gene.13$`Gene ID` == "Zm00001d014734"] <- "italic(arodeH2)[Zm00001d014734]" # manually re-name
df.gene.13$Gene[df.gene.13$`Gene ID` == "Zm00001d014737"] <- "italic(arodeH2)[Zm00001d014737]" # manually re-name

# add sds (as this was in top10% for at least one trait in NAM)
df.gene.13 <- rbind.data.frame(df.gene.13, c("Zm00001d027694", "italic(sds)"))

# remove vte5
df.gene.13 <- df.gene.13[!(df.gene.13$`Gene ID` %in% "Zm00001d001896"), ]

# setup
trait.all <- c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3",
               "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")
model <- "GTBLUP.CandGenes" # G + T_cand model
path <- "RESULT/4.2-GenExpFit_trans/Coef_GBLUP+ExpCand_"

# set seed (PCA have a random choice of direction)
set.seed(139)

# Load reg. coef. data 
coef.mat <- NULL
for ( i in 1:length(trait.all) ) {
  trait <- trait.all[i]
  file.in <- paste0(path, trait, ".csv")
  df.coef <- read.csv(file.in)
  coef.mat <- cbind(coef.mat, df.coef$coef)
}
colnames(coef.mat) <- trait.all
rownames(coef.mat) <- df.coef$GeneID

# PCA
M <- coef.mat[, c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3")] # Use 6 traits
colnames(M) <- as.character(myFun.ChangeTraitNames(colnames(M)))
pr.res <- prcomp(M, center = TRUE, scale. = TRUE)

# make plot (figure 4)
pcobj = pr.res
choices = 1:2
scale = 1
obs.scale = 1 - scale
var.scale = scale
circle.prob = 0.57
varname.adjust = 1.5
labels.size = 3

nobs.factor <- sqrt(nrow(pcobj$x) - 1)
d <- pcobj$sdev
u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
v <- pcobj$rotation
df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, FUN = "*"))
v <- sweep(v, 2, d^var.scale, FUN = "*")
df.v <- as.data.frame(v[, choices])
names(df.u) <- c("xvar", "yvar")
names(df.v) <- names(df.u)
df.u <- df.u * nobs.factor
r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2)) ^ (1/4)
v.scale <- rowSums(v^2)
df.v <- r * df.v/sqrt(max(v.scale))
u.axis.labs <- paste("Standardized PC", choices, sep = "")
u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained variance)", 100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
df.v$varname <- rownames(v)
df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)

#
m <- match(df.gene.13$`Gene ID`, rownames(df.u))
df.fig <- cbind(df.gene.13, df.u[m, 1:2])
df.u$grp <- "A"
df.u$grp[rownames(df.u) %in% df.gene.13$`Gene ID`] <- "B"

#
g <- ggplot(data = df.u, aes(x = xvar, y = yvar, col = grp))
g <- g + geom_vline(xintercept = 0, lwd = 0.5, lty = 2) + geom_hline(yintercept = 0, lwd = 0.5, lty = 2)
g <- g + geom_point()
g <- g + scale_color_manual(values = c("black", "red"))
g <- g + geom_point(data = df.u[df.u$grp == "B", ], aes(x = xvar, y = yvar), color = "red")
g <- g + geom_segment(data = df.v, 
                      aes(x = 0, y = 0, xend = xvar, yend = yvar), 
                      arrow = arrow(length = unit(1/2, "picas")), 
                      color = "darkred")
g <- g + xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) 
g <- g + theme(panel.border = element_rect(fill = NA, size = 0.5))
g <- g + theme(panel.background = element_rect(fill = "white"))
g <- g + coord_equal()
g <- g + theme(text = element_text(family = "Times New Roman"))
g <- g + theme(legend.position = "none")
g <- g + theme(panel.grid.major = element_line(colour="gray92", size = 0.5),
               panel.grid.minor = element_line(colour="gray92", size = 0.5))
ggplot2::ggsave(filename = paste0(dir.save, "/Fig4.eps"), 
                plot = g, 
                device = cairo_ps, 
                dpi = 300, 
                width = 7 * 0.8,
                height = 6 * 0.8)
write.csv(df.v, file = paste0(dir.save, "/Fig4_arrows.csv"), row.names = F)
write.csv(df.u, file = paste0(dir.save, "/Fig4_points.csv"), row.names = T)


# ---------------------------------------------------------------------------- #
# ----- A Table to make figure3
# ---------------------------------------------------------------------------- #
trait.all <- c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3",
               "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")

# NAM data
df.NAM.genes <- as.data.frame(read_excel("RAWDATA/NAM_QTL/TocoGene_NAM.xlsx"))
df.NAM.genes <- df.NAM.genes[df.NAM.genes$Gene != "dxs3", ] # remove dxs3 as it only hits for PC-8
df.NAM.genes$RefGen_v4 <- gsub(" ", "", df.NAM.genes$RefGen_v4)
df.NAM.genes.simple <- data.frame("Gene" = df.NAM.genes$Gene,
                                  "ID" = df.NAM.genes$RefGen_v4,
                                  df.NAM.genes[, trait.all],
                                  "ceeQTL" = df.NAM.genes$eQTL)

# Load reg. coef. data 
coef.mat <- NULL
for ( i in 1:length(trait.all) ) {
  trait <- trait.all[i]
  file.in <- paste0("RESULT/4.2-GenExpFit_trans/Coef_GBLUP+ExpAll_", trait, ".csv")
  df.coef <- read.csv(file.in)
  coef.mat <- cbind(coef.mat, df.coef$coef)
}
colnames(coef.mat) <- trait.all
rownames(coef.mat) <- df.coef$GeneID

# calculate importance
abs.coef.mat <- abs(coef.mat)
imp.rank.mat <- apply((-1) * abs.coef.mat, 2, rank)
imp.rank.percent.mat <- imp.rank.mat / nrow(imp.rank.mat)

# for the NAM genes 
df.NAM.genes.simple.10 <- df.NAM.genes.simple[df.NAM.genes.simple$ID %in% rownames(imp.rank.percent.mat), ]
imp.rank.percent.mat.10 <- imp.rank.percent.mat[df.NAM.genes.simple.10$ID, ]
rownames(imp.rank.percent.mat.10) <- df.NAM.genes.simple.10$Gene

# write
write.csv(round(imp.rank.percent.mat.10, 3), file = paste0(dir.save, "/Fig3_left.csv"))
write.csv(df.NAM.genes.simple.10, file = paste0(dir.save, "/Fig3_right.csv"), row.names = F)
