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
myFun.Calc.PA <- function(df.pred, df.obs) {
   # check the column names of the data
   checker <- !all(df.obs$ID %in% colnames(df.pred))
   if ( checker ) {
      out <- NULL
   } else {
      # split the data frame of the predicted values
      m <- match(df.obs$ID, colnames(df.pred))
      mat.pred.vals <- as.matrix(df.pred[, m])
      df.pred.info <- df.pred[, !(colnames(df.pred) %in% df.obs$ID)]
      
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
}
myFun.Calc.PPA <- function(df.pred, df.obs, model = "MultiKernel") {
   # check the column names of the data
   checker <- !all(df.obs$ID %in% colnames(df.pred))
   if ( checker ) {
      out <- NULL
   } else {
      win.all <- unique(df.pred$Window)
      df.res.all <- NULL
      for (w in 1:length(win.all)) {
         win <- win.all[w]
         tf <- (df.pred$Model == model) & (df.pred$Window == win)
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
               df.res.i <- data.frame(df.i[, c("Scenario", "Trait", "Model", "Window", "QTL")],
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
   df.fig$SD.bar.max[df.fig$PA * df.fig$SD.bar.max < 0] <- 0
   df.fig$SD.bar.min[df.fig$PA * df.fig$SD.bar.min < 0] <- 0
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
   x[x == "BRR_341K"] <- "BRR"
   x[x == "BayesB_341K"] <- "BayesB"
   x[x == "GBLUP_341K"] <- "GBLUP"
   x[x == "250K"] <- "MK-GBLUP\n(250 Kb)"
   x[x == "1M"] <- "MK-GBLUP\n(1 Mb)"
   x[x == "SI"] <- "MK-GBLUP\n(SI)"
   x <- factor(x, levels = c("GBLUP", 
                             "BayesB", 
                             "MK-GBLUP\n(250 Kb)", 
                             "MK-GBLUP\n(1 Mb)", 
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
myFun.FormatImpRate <- function(X) {
   X.new <- matrix("", nr = nrow(X), nc = ncol(X))
   for (i in 1:nrow(X)) {
      for (j in 1:ncol(X)) {
         x <- X[i, j]
         xf <- formatC(x, digits = 2, format = "f")
         if ( x == 0 ) { x.new <- paste0("±", xf) }
         if ( x > 0 ) { x.new <- paste0("+", xf) }
         if ( x < 0 ) { x.new <- xf }
         X.new[i, j] <- x.new
      }
   }
   return(X.new)
}
MyFun.ConvertQtlNames <- function(x){
   flag <- grep("\\(", x)
   if  ( length(flag) == 0 ) {
      name.new <- x
   } else {
      name.new <- paste0("italic(", gsub(")", "", strsplit(x, "\\(")[[1]][2]), ")")
   }
   return(name.new)
}
myFun.ChangeModeltNames.Expr <- function(x) {
   x[x == "GBLUP"] <- "GBLUP\n(Baseline)"
   x[x == "TBLUP"] <- "TBLUP\n(All genes)"
   x[x == "TBLUP_sub"] <- "TBLUP\n(Cand. genes)"
   x[x == "GBLUP+TBLUP"] <- "GBLUP+TBLUP\n(All genes)"
   x[x == "GBLUP+TBLUP_sub"] <- "GBLUP+TBLUP\n(Cand. genes)"
   x[x == "MultiNamLargeGenes"] <- "Multi-trait model\n(Large-eff. genes)"
   x[x == "MultiNamLargeGenesCovar"] <- "GBLUP with covariates\n(Large-eff. genes)"
   x <- factor(x, levels = c("GBLUP\n(Baseline)", "TBLUP\n(All genes)", 
                             "TBLUP\n(Cand. genes)", "GBLUP+TBLUP\n(All genes)",
                             "GBLUP+TBLUP\n(Cand. genes)", "Multi-trait model\n(Large-eff. genes)",
                             "GBLUP with covariates\n(Large-eff. genes)"))
   return(x)
}
myFun.ChangeModelNames.Feat <- function(x) {
   x[x == "GBLUP_12M"] <- "GBLUP"
   x[x == "GERP"] <- "GERP"
   x[x == "RecRate"] <- "Rec. rate"
   x[x == "MAF"] <- "MAF"
   x[x == "Prox_1Kb"] <- "Prox.\n(1 Kb)"
   x[x == "Prox_10Kb"] <- "Prox.\n(10 Kb)"
   x[x == "Prox_100Kb"] <- "Prox.\n(100 Kb)"
   x <- factor(x, levels = c("GBLUP", "GERP", "Rec. rate", "MAF",
                             "Prox.\n(1 Kb)", "Prox.\n(10 Kb)", "Prox.\n(100 Kb)"))
   return(x)
}

# mkdir
dir.save <- "RESULT/99-Tables_and_Figures"
dir.create(dir.save, recursive = T)

# load phenotype data
f.ames <- "RESULT/1.1-MakeAmesPhenoData/AmesPheno.csv" # Ames Phenotype data
f.gb <- "RESULT/1.2-MakeGbPhenoData/GbPheno.csv" # Gb Phenotype data
AmesPheno <- fread(f.ames, data.table = F)
GbPheno <- fread(f.gb, data.table = F)
trait.all <- c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3",
               "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")
trait.six <- c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3")
trait.non.aT <- c("d.T", "g.T", "a.T3", "d.T3", "g.T3",
                  "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")

# setup for figure
loadfonts(device = "win")
cols <- GetColors(n = 13)
cols.9 <- cols[c(2:4, 6:8, 10:12)]
cols.QTL <- c("gray", "black", cols[10:12])
cols.Expr <- c("gray", cols[c(4, 3, 8, 6, 12, 11)])
cols.GF <- c("gray", cols[c(3, 5, 8, 10:12)])

# ---------------------------------------------------------------------------- #
# ----- Table: Summary of the NAM-QTL model
# ---------------------------------------------------------------------------- #
# map for the 341K SNP set
f <- "RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.pos"
df.map <- read.table(f)
colnames(df.map) <- c("chr", "pos")
# knwon QTLs (for validation)
GenData <- read.csv("RAWDATA/NAM_QTL/qtl.data.from.di.csv", stringsAsFactors = FALSE)
colnames(GenData)[1] <- "QTL.ID"
# loop for making the table
win.all <- c("250K", "1M", "SI")
df.save <- expand.grid("Trait" = trait.all, "Window" = win.all)
df.save$n.kernel <- df.save$n.QTL <- NA
df.save$Max.n.SNP <- df.save$Min.n.SNP <- df.save$Avg.n.SNP <- NA
df.save$Max.Width <- df.save$Min.Width <- df.save$Avg.Width <- NA
for ( i in 1:nrow(df.save) ) {
   tr <- df.save$Trait[i]
   win <- df.save$Window[i]
   SnpList <- readRDS(file = paste0("RESULT/5.1-MultiKernel_Across_trans/SnpNumList_", win, "_", tr, ".Rdata"))
   n.snp.vec <- width.vec <- c()
   for ( k in 1:length(SnpList) ) {
      qtl <- names(SnpList)[k]
      if (qtl != "NonQtl") {
         SnpVec <- SnpList[[k]]
         df.map.qtl <- df.map[SnpVec, ]
         n.snp <- nrow(df.map.qtl)
         width <- (max(df.map.qtl$pos) - min(df.map.qtl$pos)) / 1000000
         n.snp.vec <- c(n.snp.vec, n.snp)
         width.vec <- c(width.vec, width)
         # validation
         chr <- unique(df.map.qtl$chr)
         if ( length(chr) != 1 ) { print("WRONG.a"); break }
         tf <- (GenData$Trait == tr) & (GenData$NAME == qtl)
         chr.GenData <- as.integer(gsub("Chr", "", GenData$Chr[tf]))
         if ( chr != chr.GenData ) { print("WRONG.b"); break }
      }
   }
   # save
   df.save$n.QTL[i] <- length(SnpList) - 1
   df.save$n.kernel[i] <- length(SnpList)
   df.save$Max.n.SNP[i] <- max(n.snp.vec)
   df.save$Min.n.SNP[i] <- min(n.snp.vec)
   df.save$Avg.n.SNP[i] <- mean(n.snp.vec)
   df.save$Max.Width[i] <- max(width.vec)
   df.save$Min.Width[i] <- min(width.vec)
   df.save$Avg.Width[i] <- mean(width.vec)
}
o <- order(df.save$Trait)
df.save.ord <- df.save[o, ]
write.csv(df.save.ord, file = paste0(dir.save, "/Summary_MKGBLUP_using_NAM_QTL.csv"), row.names = F)


# ---------------------------------------------------------------------------- #
# ----- Figure: PA on the MK-GBLUP using NAM-QTL 
# ---------------------------------------------------------------------------- #
list.all <- list(
   BaseL.Ato2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_AmesTo282_BackTransPredVal.csv"), df.obs = GbPheno),
   BaseL.2toA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_282ToAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   BaseL.AtoA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_WithinAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   BaseL.2to2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_Within282_BackTransPredVal.csv"), df.obs = GbPheno),
   UseNQ.Ato2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_AmesTo282_BackTransPredVal.csv"), df.obs = GbPheno),
   UseNQ.2toA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_282ToAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   UseNQ.AtoA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_WithinAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   UseNQ.2to2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_Within282_BackTransPredVal.csv"), df.obs = GbPheno)
)
list.all.SD <- lapply(list.all, myFun.Calc.SD)
df.fig <- ldply(list.all.SD, data.frame)
df.fig <- df.fig[!(df.fig$Model %in% c("BRR_341K", "GBLUP_12M")), ]
df.fig$Model <- myFun.ChangeModeltNames.NamQtl(df.fig$Model)
df.fig$Trait <- myFun.ChangeTraitNames(df.fig$Trait)
df.fig$Scenario <- factor(df.fig$Scenario, 
                          levels = c("Within Ames", "From Ames to 282", 
                                     "From 282 to Ames", "Within 282"))
df.fig <- myFun.Attach.SD.bar(df.fig)
p <- ggplot(df.fig, aes(x = Model, y = PA, fill = Trait))
p <- p + facet_wrap(~ Scenario)
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
                           ymax = SD.bar.max), 
                       position = position_dodge(.9), width=.5)
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times New Roman"))
p <- p + ylab("Predictive ability")
p <- p + scale_fill_manual(values = cols.9)
ggsave(p, file = paste0(dir.save, "/Figure_PA_UseNamQtl.png"), width = 9, height = 6)

# Ver.2
p <- ggplot(df.fig, aes(x = Trait, y = PA, fill = Model))
p <- p + facet_wrap(~ Scenario)
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
                           ymax = SD.bar.max), 
                       position = position_dodge(.9), width=.5)
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times New Roman"))
p <- p + ylab("Predictive ability")
p <- p + scale_fill_manual(values = cols.QTL)
ggsave(p, file = paste0(dir.save, "/Figure_PA_UseNamQtl_Ver2.png"), width = 12, height = 6)


# ---------------------------------------------------------------------------- #
# ----- Figure: improvement of PA on the MK-GBLUP using NAM-QTL 
# ---------------------------------------------------------------------------- #
list.all <- list(
   BaseL.Ato2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_AmesTo282_BackTransPredVal.csv"), df.obs = GbPheno),
   BaseL.2toA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_282ToAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   BaseL.AtoA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_WithinAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   BaseL.2to2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_Within282_BackTransPredVal.csv"), df.obs = GbPheno),
   UseNQ.Ato2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_AmesTo282_BackTransPredVal.csv"), df.obs = GbPheno),
   UseNQ.2toA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_282ToAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   UseNQ.AtoA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_WithinAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   UseNQ.2to2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_Within282_BackTransPredVal.csv"), df.obs = GbPheno)
)
df.all <- ldply(list.all, data.frame)
df.all$PA <- myFun.Calc.Imp(df = df.all, base.model = "GBLUP_341K")
df.all <- myFun.Calc.SD(df.all)
df.fig <- df.all[!(df.all$Model %in% c("BRR_341K", "GBLUP_341K", "GBLUP_12M")), ]
df.fig$Model <- myFun.ChangeModeltNames.NamQtl(df.fig$Model)
df.fig$Trait <- myFun.ChangeTraitNames(df.fig$Trait)
df.fig$Scenario <- factor(df.fig$Scenario, 
                          levels = c("Within Ames", "From Ames to 282", 
                                     "From 282 to Ames", "Within 282"))
df.fig <- myFun.Attach.SD.bar(df.fig)
cols <- GetColors(n = 13); cols.9 <- cols[c(2:4, 6:8, 10:12)]
p <- ggplot(df.fig, aes(x = Model, y = PA, fill = Trait))
p <- p + facet_wrap(~ Scenario)
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
                           ymax = SD.bar.max), 
                       position = position_dodge(.9), width=.5)
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times New Roman"))
p <- p + ylab("Improvement rate (%) from GBLUP")
p <- p + scale_fill_manual(values = cols.9)
ggsave(p, file = paste0(dir.save, "/Figure_PA_Imp_UseNamQtl.png"), width = 9, height = 6)

# Ver.2
p <- ggplot(df.fig, aes(x = Trait, y = PA, fill = Model))
p <- p + facet_wrap(~ Scenario)
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
                           ymax = SD.bar.max), 
                       position = position_dodge(.9), width=.5)
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times New Roman"))
p <- p + ylab("Improvement rate (%) from GBLUP")
p <- p + scale_fill_manual(values = cols.QTL[-1])
ggsave(p, file = paste0(dir.save, "/Figure_PA_Imp_UseNamQtl_ver2.png"), width = 12, height = 6)


# ---------------------------------------------------------------------------- #
# ----- Table: Result of the MK-GBLUP using NAM-QTL
# ---------------------------------------------------------------------------- #
list.all <- list(
   BaseL.Ato2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_AmesTo282_BackTransPredVal.csv"), df.obs = GbPheno),
   BaseL.2toA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_282ToAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   BaseL.AtoA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_WithinAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   BaseL.2to2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_Within282_BackTransPredVal.csv"), df.obs = GbPheno),
   UseNQ.Ato2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_AmesTo282_BackTransPredVal.csv"), df.obs = GbPheno),
   UseNQ.2toA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_282ToAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   UseNQ.AtoA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_WithinAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   UseNQ.2to2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseNamQtl_Within282_BackTransPredVal.csv"), df.obs = GbPheno)
)
list.all.SD <- lapply(list.all, myFun.Calc.SD)
df.avg <- ldply(list.all.SD, data.frame)
df.imp <- ldply(list.all, data.frame)
df.imp$PA <- myFun.Calc.Imp(df = df.imp, base.model = "GBLUP_341K")
df.imp <- myFun.Calc.SD(df.imp)
colnames(df.avg)[colnames(df.avg) == "PA"] <- "AVG.PA"
colnames(df.imp)[colnames(df.imp) == "PA"] <- "AVG.PA.IMP"
df.long <- df.avg[, c("Trait", "Model", "Scenario", "AVG.PA")]
df.imp <- df.imp[, c("Trait", "Model", "Scenario", "AVG.PA.IMP")]
for ( i in 1:nrow(df.long) ) {
   tf <- (df.imp$Trait == df.long$Trait[i]) & 
      (df.imp$Model == df.long$Model[i]) & 
      (df.imp$Scenario == df.long$Scenario[i])
   imp.i <- df.imp$AVG.PA.IMP[tf]
   df.long$AVG.PA.IMP[i] <- imp.i
}
df.long <- df.long[!(df.long$Model %in% c("BRR_341K", "GBLUP_12M")), ]
df.long$Model <- factor(df.long$Model, levels = c("GBLUP_341K", "BayesB_341K", "250K", "1M", "SI"))
df.long$Scenario <- factor(df.long$Scenario, 
                           levels = c("From Ames to 282", "From 282 to Ames", "Within Ames", "Within 282"))
df.long$Trait <- factor(df.long$Trait, 
                        levels = c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3",
                                   "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols"))
df.PA <- dcast(df.long, Scenario + Model  ~ Trait, value.var = "AVG.PA")
df.PA$Avg.Six <- apply(df.PA[, trait.six], 1, mean)
df.PA$Avg.non.aT <- apply(df.PA[, trait.non.aT], 1, mean)
df.PA$Avg <- apply(df.PA[, trait.all], 1, mean)
df.PA.IMP <- dcast(df.long, Scenario + Model  ~ Trait, value.var = "AVG.PA.IMP")
df.PA.IMP$Avg.Six <- apply(df.PA.IMP[, trait.six], 1, mean)
df.PA.IMP$Avg.non.aT <- apply(df.PA.IMP[, trait.non.aT], 1, mean)
df.PA.IMP$Avg <- apply(df.PA.IMP[, trait.all], 1, mean)
M1 <- formatC(as.matrix(df.PA[, 3:ncol(df.PA)]), digits = 2, format = "f")
M2 <- myFun.FormatImpRate(as.matrix(df.PA.IMP[, 3:ncol(df.PA)]))
M <- matrix("", nr = nrow(M1), nc = ncol(M1))
for (j in 1:nrow(M1)) {
   for (k in 1:ncol(M1)){
      M[j, k] <- paste0(M1[j, k], " (", M2[j, k], ")")
   }
}
colnames(M) <- colnames(M1)
mat.save <- cbind(as.matrix(df.PA[, 1:2]), M)
write.csv(mat.save, file = paste0(dir.save, "/Table_PA_UseNamQtl.csv"), row.names = F)

# ---------------------------------------------------------------------------- #
# ----- Figure & Table: PA on the MK-GBLUP using NAM-QTL vs PVE
# ---------------------------------------------------------------------------- #
# load NAM PVE 
NamData.Full <- read.csv("RAWDATA/NAM_QTL/qtl.data.from.di.csv")
colnames(NamData.Full)[1] <- "QTL.ID"
NamData <- NamData.Full[NamData.Full$Trait != "PC-8", ]
NamData$NAME3 <- paste0("QTL", NamData$QTL.ID)
# map for the 341K SNP set
f <- "RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.pos"
df.map <- read.table(f)
colnames(df.map) <- c("chr", "pos")
# load prediction results
list.all <- list(
   UseNQ.Ato2 = myFun.Calc.PPA(df.pred = MyFun.load("RESULT/99-SummaryFiles/PostHoc_UseNamQtl_AmesTo282_EvalEachKernel.csv"), df.obs = GbPheno, model = "MultiKernel"),
   UseNQ.2toA = myFun.Calc.PPA(df.pred = MyFun.load("RESULT/99-SummaryFiles/PostHoc_UseNamQtl_282ToAmes_EvalEachKernel.csv"), df.obs = AmesPheno, model = "MultiKernel")
)
df.all <- ldply(list.all, data.frame)
# assign the number of SNPs
df.all$n.SNP <- NA
for ( i in 1:nrow(df.all) ) {
   tr.i <- df.all$Trait[i]
   win.i <- df.all$Window[i]
   qtl.i <- df.all$QTL[i]
   SnpList <- readRDS(file = paste0("RESULT/5.1-MultiKernel_Across_trans/SnpNumList_", win.i, "_", tr.i, ".Rdata"))
   names(SnpList) <- sapply(names(SnpList), function(x){strsplit(x, "\\(")[[1]][1]}, USE.NAMES = F)
   df.all$n.SNP[i] <- length(SnpList[[qtl.i]])
}
# get PVE & gene names
for ( i in 1:nrow(df.all) ) {
   qtl <- df.all$QTL[i]
   trait <- df.all$Trait[i]
   if (qtl != "NonQtl" ) {
      pve <- NamData[(NamData$NAME3 == qtl) & (NamData$Trait == trait), "PVE"]
      name <- NamData[(NamData$NAME3 == qtl) & (NamData$Trait == trait), "NAME"]
   } else {
      pve <- NA
      name <- qtl
   }
   df.all$PVE[i] <- pve
   df.all$NAME[i] <- name
}
# save the data frame
df.save <- df.all[, c("Scenario", "Trait", "Window", "NAME", "n.SNP", "Partial.Cor", "P.value", "PVE")]
write.csv(df.save, file = paste0(dir.save, "/Table_Partial_Correlation_and_PVE.csv"), row.names = F)
# make figures (1)
df.all$NAME <- sapply(df.all$NAME, MyFun.ConvertQtlNames, USE.NAMES = F)
df.fig <- df.all[, c("Scenario", "Trait", "Window", "NAME", "Partial.Cor", "P.value", "PVE")]
df.fig <- df.fig[df.fig$NAME != "NonQtl", ]
df.fig$Trait <- myFun.ChangeTraitNames(df.fig$Trait)
df.fig$Scenario <- factor(df.fig$Scenario, 
                          levels = c("From Ames to 282", "From 282 to Ames"))
cols <- GetColors(n = 13); cols.9 <- cols[c(2:4, 6:8, 10:12)]
df.fig$Size <- 1
df.fig$Size[df.fig$PVE > 0.1] <- 2
win.all <- as.character(unique(df.fig$Window))
for ( i in 1:length(win.all) ) {
   set.seed(100)
   win <- win.all[i]
   df.fig.sub <- df.fig[df.fig$Window == win, ]
   p <- ggplot(df.fig.sub, aes(x = PVE, y = Partial.Cor, color = Trait, size = Size))
   p <- p + facet_wrap(~ Scenario)
   p <- p + geom_point()
   p <- p + geom_text_repel(data = subset(df.fig.sub, PVE > 0.1),
                            aes(label = NAME),
                            box.padding = unit(0.15, "lines"),
                            point.padding = unit(0.15, "lines"),
                            parse = TRUE,
                            family = "Times New Roman")
   p <- p + theme_bw()
   p <- p + theme(text = element_text(family = "Times New Roman"))
   p <- p + scale_color_manual(values = cols.9)
   p <- p + xlab("PVE in NAM population")
   p <- p + ylab("Predictive ability (Partial correlation)")
   p <- p + scale_size(range = c(1, 2.5), guide = F)
   p <- p + guides(colour = guide_legend(override.aes = list(size = 2)))
   ggsave(filename = paste0(dir.save, "/Figure_PA_with_PVE_WindowSize_", win, ".png"), p, width = 9, height = 5)
}
# make figures (2)
df.fig <- df.all[, c("Scenario", "Trait", "Window", "NAME", "Partial.Cor", "n.SNP", "PVE")]
df.fig <- df.fig[df.fig$NAME != "NonQtl", ]
df.fig$Trait <- myFun.ChangeTraitNames(df.fig$Trait)
df.fig$Scenario <- factor(df.fig$Scenario, 
                          levels = c("From Ames to 282", "From 282 to Ames"))
df.fig$Window <- factor(df.fig$Window, levels = c("250K", "1M", "SI"))

cols <- GetColors(n = 13); cols.9 <- cols[c(2:4, 6:8, 10:12)]
p <- ggplot(df.fig, aes(x = n.SNP, y = Partial.Cor))
p <- p + facet_grid(Scenario ~ Window, scales = "free_x")
p <- p + geom_point(aes(col = Trait))
p <- p + geom_smooth(method = "lm", color = "black", formula = "y ~ x")
p <- p + stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                  size = 3)
p <- p + xlab("# of SNPs")
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times New Roman"))
p <- p + scale_color_manual(values = cols.9)
p <- p + ylab("Predictive ability (Partial correlation)")
ggsave(filename = paste0(dir.save, "/Figure_Partial_Correlation_and_the_number_of_SNPs.png"), p, 
       width = 9, height = 5)


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
cols <- GetColors(n = 13); cols.9 <- cols[c(2:4, 6:8, 10:12)]
p <- ggplot(df.fig, aes(x = Model, y = PA, fill = Trait))
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
                           ymax = SD.bar.max), 
                       position = position_dodge(.9), width=.5)
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times New Roman"))
p <- p + ylab("Predictive ability")
p <- p + scale_fill_manual(values = cols.9)
ggsave(p, file = paste0(dir.save, "/Figure_PA_TrsctPred.png"), width = 9, height = 5)

#
p <- ggplot(df.fig, aes(x = Trait, y = PA, fill = Model))
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
                           ymax = SD.bar.max), 
                       position = position_dodge(.9), width=.5)
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times New Roman"))
p <- p + ylab("Predictive ability")
p <- p + scale_fill_manual(values = cols.Expr)
ggsave(p, file = paste0(dir.save, "/Figure_PA_TrsctPred_Ver2.png"), width = 12, height = 5)


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

cols <- GetColors(n = 13); cols.9 <- cols[c(2:4, 6:8, 10:12)]
p <- ggplot(df.fig, aes(x = Model, y = PA, fill = Trait))
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
                           ymax = SD.bar.max), 
                       position = position_dodge(.9), width=.5)
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times New Roman"))
p <- p + ylab("Improvement Rate (%) from GBLUP")
p <- p + scale_fill_manual(values = cols.9)
ggsave(p, file = paste0(dir.save, "/Figure_PA_Imp_TrsctPred.png"), width = 9, height = 5)

#
p <- ggplot(df.fig, aes(x = Trait, y = PA, fill = Model))
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
                           ymax = SD.bar.max), 
                       position = position_dodge(.9), width=.5)
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times New Roman"))
p <- p + ylab("Improvement Rate (%) from GBLUP")
p <- p + scale_fill_manual(values = cols.Expr[-1])
ggsave(p, file = paste0(dir.save, "/Figure_PA_Imp_TrsctPred_Ver2.png"), width = 12, height = 5)


# ---------------------------------------------------------------------------- #
# ----- Table: Transcriptome based prediction
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
df$Model <- factor(df$Model, levels = c("GBLUP", "TBLUP", "TBLUP_sub",
                                        "GBLUP+TBLUP", "GBLUP+TBLUP_sub",
                                        "MultiNamLargeGenes",
                                        "MultiNamLargeGenesCovar"))
df.01 <- aggregate(PA ~ Trait + Model, df[, c("Model", "Trait", "Rep", "PA")], mean)
df.02 <- aggregate(PA.IMP ~ Trait + Model, df[, c("Model", "Trait", "Rep", "PA.IMP")], mean)
df.long <- merge.data.frame(df.01, df.02)
df.PA <- dcast(df.long, Model  ~ Trait, value.var = "PA")
df.PA$Avg.Six <- apply(df.PA[, trait.six], 1, mean)
df.PA$Avg.non.aT <- apply(df.PA[, trait.non.aT], 1, mean)
df.PA$Avg <- apply(df.PA[, trait.all], 1, mean)
df.PA.IMP <- dcast(df.long, Model  ~ Trait, value.var = "PA.IMP")
df.PA.IMP$Avg.Six <- apply(df.PA.IMP[, trait.six], 1, mean)
df.PA.IMP$Avg.non.aT <- apply(df.PA.IMP[, trait.non.aT], 1, mean)
df.PA.IMP$Avg <- apply(df.PA.IMP[, trait.all], 1, mean)
M1 <- formatC(as.matrix(df.PA[, 2:ncol(df.PA)]), digits = 2, format = "f")
M2 <- myFun.FormatImpRate(as.matrix(df.PA.IMP[, 2:ncol(df.PA)]))
M <- matrix("", nr = nrow(M1), nc = ncol(M1))
for (j in 1:nrow(M1)) {
   for (k in 1:ncol(M1)){
      M[j, k] <- paste0(M1[j, k], " (", M2[j, k], ")")
   }
}
colnames(M) <- colnames(M1)
mat.save <- cbind(as.matrix(df.PA[, 1]), M)
write.csv(mat.save, file = paste0(dir.save, "/Table_PA_TrsctPred.csv"), row.names = F)


# ---------------------------------------------------------------------------- #
# ----- Figure: PEER
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
df <- df[df$Model %in% c("TBLUP", "TBLUP_sub", "GBLUP+TBLUP", "GBLUP+TBLUP_sub"), ]
df.00 <- data.frame("Trait" = df$Trait,
                    "Model" = df$Model,
                    "n.PEER.Factor" = 0,
                    "Rep" = df$Rep,
                    "PA" = df$PA)
pred.peer <- read.csv("RESULT/99-SummaryFiles/PostHoc_TrsctPred_WithinAmes_EvalNumPeerFact_BackTransPredVal.csv")
df <- pred.peer[, 1:4]
pred.mat <- as.matrix(pred.peer[, 5:ncol(pred.peer)])
m <- match(colnames(pred.mat), AmesPheno$ID)
AmesPheno.545 <- AmesPheno[m, ]
df$PA <- NA
for ( i in 1:nrow(df) ) {
   tr.i <- df$Trait[i]
   obs <- AmesPheno.545[[tr.i]]
   pred <- pred.mat[i, ]
   df$PA[i] <- cor(obs, pred, use = "p")
}
df.long <- rbind.data.frame(df.00, df)
df.fig <- myFun.Calc.SD(df.long)
df.fig$Trait <- myFun.ChangeTraitNames(df.fig$Trait)
df.fig$Model[df.fig$Model == "TBLUP"] <- "TBLUP (All genes)"
df.fig$Model[df.fig$Model == "TBLUP_sub"] <- "TBLUP (Cand. genes)"
df.fig$Model[df.fig$Model == "GBLUP+TBLUP"] <- "GBLUP+TBLUP (All genes)"
df.fig$Model[df.fig$Model == "GBLUP+TBLUP_sub"] <- "GBLUP+TBLUP (Cand. genes)"
df.fig$Model <- factor(df.fig$Model, levels = c("TBLUP (All genes)", "TBLUP (Cand. genes)",
                                                "GBLUP+TBLUP (All genes)", "GBLUP+TBLUP (Cand. genes)"))
df.fig <- myFun.Attach.SD.bar(df.fig)
p <- ggplot(df.fig, aes(x = n.PEER.Factor, y = PA))
p <- p + facet_grid(Model ~ Trait, scales = "free_y") 
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
                           ymax = SD.bar.max), 
                       position = position_dodge(.9), width=.5)
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times New Roman"))
p <- p + ylab("Predictive ability")
p <- p + xlab("# of PEER factors")
ggsave(p, file = paste0(dir.save, "/Figure_PA_TrsctPred_num_PEER.png"), width = 15, height = 9)


# ---------------------------------------------------------------------------- #
# ----- Figure: PA on the prediction using genomic features 
# ---------------------------------------------------------------------------- #
list.all <- list(
   BaseL.Ato2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_AmesTo282_BackTransPredVal.csv"), df.obs = GbPheno),
   BaseL.2toA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_282ToAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   BaseL.AtoA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_WithinAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   BaseL.2to2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_Within282_BackTransPredVal.csv"), df.obs = GbPheno),
   UseGF.Ato2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseFeat_12M_AmesTo282_BackTransPredVal.csv"), df.obs = GbPheno),
   UseGF.2toA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseFeat_12M_282ToAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   UseGF.AtoA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseFeat_12M_WithinAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   UseGF.2to2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseFeat_12M_Within282_BackTransPredVal.csv"), df.obs = GbPheno)
)
list.all.SD <- lapply(list.all, myFun.Calc.SD)
df.fig <- ldply(list.all.SD, data.frame)
df.fig <- df.fig[!(df.fig$Model %in% c("BRR_341K", "BayesB_341K", "GBLUP_341K")), ]
df.fig <- myFun.Attach.SD.bar(df.fig)
df.fig$Model <- myFun.ChangeModelNames.Feat(df.fig$Model)
df.fig$Trait <- myFun.ChangeTraitNames(df.fig$Trait)
df.fig$Scenario <- factor(df.fig$Scenario, 
                          levels = c("Within Ames", "From Ames to 282", 
                                     "From 282 to Ames", "Within 282"))
cols <- GetColors(n = 13); cols.9 <- cols[c(2:4, 6:8, 10:12)]
p <- ggplot(df.fig, aes(x = Model, y = PA, fill = Trait))
p <- p + facet_wrap(~ Scenario)
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
                           ymax = SD.bar.max), 
                       position = position_dodge(.9), width=.5)
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times New Roman"))
p <- p + ylab("Predictive ability")
p <- p + scale_fill_manual(values = cols.9)
ggsave(p, file = paste0(dir.save, "/Figure_PA_UseFeat.png"), width = 9, height = 6)

#
p <- ggplot(df.fig, aes(x = Trait, y = PA, fill = Model))
p <- p + facet_wrap(~ Scenario)
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
                           ymax = SD.bar.max), 
                       position = position_dodge(.9), width=.5)
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times New Roman"))
p <- p + ylab("Predictive ability")
p <- p + scale_fill_manual(values = cols.GF)
ggsave(p, file = paste0(dir.save, "/Figure_PA_UseFeat_Ver2.png"), width = 12, height = 6)



# ---------------------------------------------------------------------------- #
# ----- Figure: improvement of PA on the prediction using genomic features 
# ---------------------------------------------------------------------------- #
list.all <- list(
   BaseL.Ato2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_AmesTo282_BackTransPredVal.csv"), df.obs = GbPheno),
   BaseL.2toA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_282ToAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   BaseL.AtoA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_WithinAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   BaseL.2to2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_Within282_BackTransPredVal.csv"), df.obs = GbPheno),
   UseGF.Ato2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseFeat_12M_AmesTo282_BackTransPredVal.csv"), df.obs = GbPheno),
   UseGF.2toA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseFeat_12M_282ToAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   UseGF.AtoA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseFeat_12M_WithinAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   UseGF.2to2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseFeat_12M_Within282_BackTransPredVal.csv"), df.obs = GbPheno)
)
df.all <- ldply(list.all, data.frame)
df.all$PA <- myFun.Calc.Imp(df = df.all, base.model = "GBLUP_12M")
df.all <- myFun.Calc.SD(df.all)
df.all <- df.all[!(df.all$Model %in% c("BayesB_341K", "BRR_341K", "GBLUP_341K")), ]
df.all$Model <- myFun.ChangeModelNames.Feat(df.all$Model)
df.all$Trait <- myFun.ChangeTraitNames(df.all$Trait)
df.all$Scenario <- factor(df.all$Scenario, 
                          levels = c("Within Ames", "From Ames to 282", 
                                     "From 282 to Ames", "Within 282"))
df.all <- myFun.Attach.SD.bar(df.all)
df.fig <- df.all[!(df.all$Model %in% c("GBLUP")), ]
cols <- GetColors(n = 13); cols.9 <- cols[c(2:4, 6:8, 10:12)]
p <- ggplot(df.fig, aes(x = Model, y = PA, fill = Trait))
p <- p + facet_wrap(~ Scenario)
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
                           ymax = SD.bar.max), 
                       position = position_dodge(.9), width=.5)
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times New Roman"))
p <- p + ylab("Improvement rate (%) from GBLUP")
p <- p + scale_fill_manual(values = cols.9)
ggsave(p, file = paste0(dir.save, "/Figure_PA_Imp_UseFeat.png"), width = 9, height = 6)

#
p <- ggplot(df.fig, aes(x = Trait, y = PA, fill = Model))
p <- p + facet_wrap(~ Scenario)
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
                           ymax = SD.bar.max), 
                       position = position_dodge(.9), width=.5)
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times New Roman"))
p <- p + ylab("Improvement rate (%) from GBLUP")
p <- p + scale_fill_manual(values = cols.GF[-1])
ggsave(p, file = paste0(dir.save, "/Figure_PA_Imp_UseFeat_Ver2.png"), width = 12, height = 6)



# ---------------------------------------------------------------------------- #
# ----- Figure for each genomic feature
# ---------------------------------------------------------------------------- #
list.all <- list(
   BaseL.Ato2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_AmesTo282_BackTransPredVal.csv"), df.obs = GbPheno),
   BaseL.2toA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_282ToAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   UseGF.Ato2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/PostHoc_UseFeat_AmesTo282_EvalEachKernel_BackTransPredVal.csv"), df.obs = GbPheno),
   UseGF.2toA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/PostHoc_UseFeat_282ToAmes_EvalEachKernel_BackTransPredVal.csv"), df.obs = AmesPheno)
)
df.all <- ldply(list.all, data.frame)
df.all <- df.all[df.all$Model %in% c("GBLUP_12M", "SingleKernel"), ]
df.all$Feature[is.na(df.all$Feature)] <- "GBLUP"
df.all <- df.all[, c("Scenario", "Feature", "Group", "Trait", 'PA')]
df.all$Group[is.na(df.all$Group)] <- "GBLUP"
df.all <- df.all[, c("Scenario", "Group", "Trait", "PA")]
df.all$Rep <- NA
colnames(df.all)[colnames(df.all) == "Group"] <- "Model"
df.all$Imp <- myFun.Calc.Imp(df = df.all, base.model = "GBLUP")
df.all$Scenario <- factor(df.all$Scenario, 
                          levels = c("From Ames to 282", "From 282 to Ames"))
df.all$Trait <- myFun.ChangeTraitNames(df.all$Trait)
df.all$Model[df.all$Model == "GBLUP"] <- "GBLUP (Baseline)"
df.all$Model[df.all$Model == "GERP_posi"] <- "GERP: positive"
df.all$Model[df.all$Model == "GERP_nega"] <- "GERP: negative"
df.all$Model[df.all$Model == "MAF_high"] <- "MAF: High"
df.all$Model[df.all$Model == "MAF_mid"] <- "MAF: Middle"
df.all$Model[df.all$Model == "RecRate_high"] <- "Rec. rate: High"
df.all$Model[df.all$Model == "RecRate_mid"] <- "Rec. rate: Middle"
df.all$Model[df.all$Model == "RecRate_low"] <- "Rec. rate: Low"
df.all$Model[df.all$Model == "prox_1kb"] <- "Proximal (1 Kb)"
df.all$Model[df.all$Model == "dist_1kb"] <- "Distal (1 Kb)"
df.all$Model[df.all$Model == "prox_10kb"] <- "Proximal (10 Kb)"
df.all$Model[df.all$Model == "dist_10kb"] <- "Distal (10 Kb)"
df.all$Model[df.all$Model == "prox_100kb"] <- "Proximal (100 Kb)"
df.all$Model[df.all$Model == "dist_100kb"] <- "Distal (100 Kb)"
df.all$Model <- factor(df.all$Model, levels = rev(c("GBLUP (Baseline)",
                                                    "GERP: positive", 
                                                    "GERP: negative",
                                                    "MAF: High", 
                                                    "MAF: Middle", 
                                                    "Rec. rate: High",
                                                    "Rec. rate: Middle",
                                                    "Rec. rate: Low",
                                                    "Proximal (1 Kb)",
                                                    "Proximal (10 Kb)",
                                                    "Proximal (100 Kb)",
                                                    "Distal (1 Kb)",
                                                    "Distal (10 Kb)",
                                                    "Distal (100 Kb)")))
p <- ggplot(df.all, aes(x = Trait, y = Model, fill=Imp))
p <- p + facet_wrap(~Scenario)
p <- p + geom_tile(color = "black", lwd = 0.5, linetype = 1)
p <- p + geom_text(aes(label = formatC(PA, digits = 2, format = "f")),
                   size = 2.5, family = "Times New Roman")
p <- p + scale_fill_gradient2(low = "#075AFF",
                              mid = "#FFFFFF",
                              high = "red",
                              midpoint = 0,
                              name = "Improvement rate\nfrom GBLUP (%)")
p <- p + theme_minimal()
p <- p + theme(panel.background = element_blank(), 
               axis.ticks.y = element_blank(),
               axis.ticks.x = element_blank(),
               panel.grid.major =  element_blank(),
               panel.grid.minor =  element_blank())
p <- p + theme(text = element_text(family = "Times New Roman"))
p <- p + coord_fixed()
ggsave(paste0(dir.save, "/Figrue_PA_SingleKernel_Feat.png"), p,
       width = 10, height = 7)


# ---------------------------------------------------------------------------- #
# ----- Table: Result of the MK-GBLUP using genomic features
# ---------------------------------------------------------------------------- #
list.all <- list(
   BaseL.Ato2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_AmesTo282_BackTransPredVal.csv"), df.obs = GbPheno),
   BaseL.2toA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_282ToAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   BaseL.AtoA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_WithinAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   BaseL.2to2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/BaselineModels_Within282_BackTransPredVal.csv"), df.obs = GbPheno),
   UseGF.Ato2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseFeat_12M_AmesTo282_BackTransPredVal.csv"), df.obs = GbPheno),
   UseGF.2toA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseFeat_12M_282ToAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   UseGF.AtoA = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseFeat_12M_WithinAmes_BackTransPredVal.csv"), df.obs = AmesPheno),
   UseGF.2to2 = myFun.Calc.PA(df.pred = MyFun.load("RESULT/99-SummaryFiles/UseFeat_12M_Within282_BackTransPredVal.csv"), df.obs = GbPheno)
)
list.all.SD <- lapply(list.all, myFun.Calc.SD)
df.avg <- ldply(list.all.SD, data.frame)
df.imp <- ldply(list.all, data.frame)
df.imp$PA <- myFun.Calc.Imp(df = df.imp, base.model = "GBLUP_12M")
df.imp <- myFun.Calc.SD(df.imp)
colnames(df.avg)[colnames(df.avg) == "PA"] <- "AVG.PA"
colnames(df.imp)[colnames(df.imp) == "PA"] <- "AVG.PA.IMP"
df.long <- df.avg[, c("Trait", "Model", "Scenario", "AVG.PA")]
df.imp <- df.imp[, c("Trait", "Model", "Scenario", "AVG.PA.IMP")]
for ( i in 1:nrow(df.long) ) {
   tf <- (df.imp$Trait == df.long$Trait[i]) & 
      (df.imp$Model == df.long$Model[i]) & 
      (df.imp$Scenario == df.long$Scenario[i])
   imp.i <- df.imp$AVG.PA.IMP[tf]
   df.long$AVG.PA.IMP[i] <- imp.i
}
df.long <- df.long[!(df.long$Model %in% c("BayesB_341K", "BRR_341K", "GBLUP_341K")), ]
df.long$Model <- factor(df.long$Model, levels = c("GBLUP_12M", "GERP", "RecRate", "MAF",
                                                  "Prox_1Kb", "Prox_10Kb", "Prox_100Kb"))
df.long$Scenario <- factor(df.long$Scenario, 
                           levels = c("From Ames to 282", "From 282 to Ames", "Within Ames", "Within 282"))
df.long$Trait <- factor(df.long$Trait, 
                        levels = c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3",
                                   "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols"))
df.PA <- dcast(df.long, Scenario + Model  ~ Trait, value.var = "AVG.PA")
df.PA$Avg.Six <- apply(df.PA[, trait.six], 1, mean)
df.PA$Avg.non.aT <- apply(df.PA[, trait.non.aT], 1, mean)
df.PA$Avg <- apply(df.PA[, trait.all], 1, mean)
df.PA.IMP <- dcast(df.long, Scenario + Model  ~ Trait, value.var = "AVG.PA.IMP")
df.PA.IMP$Avg.Six <- apply(df.PA.IMP[, trait.six], 1, mean)
df.PA.IMP$Avg.non.aT <- apply(df.PA.IMP[, trait.non.aT], 1, mean)
df.PA.IMP$Avg <- apply(df.PA.IMP[, trait.all], 1, mean)
M1 <- formatC(as.matrix(df.PA[, 3:ncol(df.PA)]), digits = 2, format = "f")
M2 <- myFun.FormatImpRate(as.matrix(df.PA.IMP[, 3:ncol(df.PA)]))
M <- matrix("", nr = nrow(M1), nc = ncol(M1))
for (j in 1:nrow(M1)) {
   for (k in 1:ncol(M1)){
      M[j, k] <- paste0(M1[j, k], " (", M2[j, k], ")")
   }
}
colnames(M) <- colnames(M1)
mat.save <- cbind(as.matrix(df.PA[, 1:2]), M)
write.csv(mat.save, file = paste0(dir.save, "/Table_PA_UseFeat.csv"), row.names = F)
