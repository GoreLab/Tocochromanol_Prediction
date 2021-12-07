# package
library(data.table)
library(ggplot2)
library(extrafont)
library(inlmisc)
library(plyr)
library(reshape2)
library(car)
library(ggrepel)
library(grid)

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
myFun.ChangeNames.Feat <- function(x) {
   x[x == "GBLUP_12M"] <- "GBLUP"
   x[x == "GERP"] <- "GERP"
   x[x == "RecRate"] <- "Recombination rate"
   x[x == "MAF"] <- "MAF"
   x[x == "Prox_1Kb"] <- "Gene Prox. (1 Kb)"
   x[x == "Prox_10Kb"] <- "Gene Prox. (10 Kb)"
   x[x == "Prox_100Kb"] <- "Gene Prox. (100 Kb)"
   x <- factor(x, levels = c("GBLUP", "GERP", "Recombination rate", "MAF",
                             "Gene Prox. (1 Kb)", "Gene Prox. (10 Kb)", "Gene Prox. (100 Kb)"))
   return(x)
}
myFun.ChangeNames.Feat.v2 <- function(x) {
   x[x == "GBLUP_12M"] <- "GBLUP"
   x[x == "GERP"] <- "GERP"
   x[x == "RecRate"] <- "Recombination rate"
   x[x == "MAF"] <- "MAF"
   x[x == "Prox_1Kb"] <- "Gene Proximity"
   x[x == "Prox_10Kb"] <- "Gene Proximity"
   x[x == "Prox_100Kb"] <- "Gene Proximity"
   x <- factor(x, levels = c("GBLUP", "GERP",
                             "Recombination rate", "MAF",
                             "Gene Proximity"))
   return(x)
}
myFun.ChangeNames.Grp <- function(x) {
   x[x == "GERP_posi"] <- "Positive"
   x[x == "GERP_nega"] <- "Negative"
   x[x == "MAF_high"] <- "High"
   x[x == "MAF_mid"] <- "Middle"
   x[x == "RecRate_high"] <- "High"
   x[x == "RecRate_mid"] <- "Middle"
   x[x == "RecRate_low"] <- "Low"
   x[x == "prox_1kb"] <- "Proximal"
   x[x == "dist_1kb"] <- "Distal"
   x[x == "prox_10kb"] <- "Proximal"
   x[x == "dist_10kb"] <- "Distal"
   x[x == "prox_100kb"] <- "Proximal"
   x[x == "dist_100kb"] <- "Distal"
   x <- factor(x, levels = c("Negative", "Positive",
                             "Low", "Middle", "High",
                             "Distal", "Proximal"))
   return(x)
}
myFun.ChangeNames.Grp.v2 <- function(x) {
   x[x == "GERP_posi"] <- "Positive"
   x[x == "GERP_nega"] <- "Negative"
   x[x == "MAF_high"] <- "High"
   x[x == "MAF_mid"] <- "Middle"
   x[x == "RecRate_high"] <- "High"
   x[x == "RecRate_mid"] <- "Middle"
   x[x == "RecRate_low"] <- "Low"
   x[x == "prox_1kb"] <- "Proximal\n(1 Kb)"
   x[x == "dist_1kb"] <- "Distal\n(1 Kb)"
   x[x == "prox_10kb"] <- "Proximal\n(10 Kb)"
   x[x == "dist_10kb"] <- "Distal\n(10 Kb)"
   x[x == "prox_100kb"] <- "Proximal\n(100 Kb)"
   x[x == "dist_100kb"] <- "Distal\n(100 Kb)"
   x <- factor(x, levels = c("Negative", "Positive",
                             "Low", "Middle", "High",
                             "Distal\n(100 Kb)", "Distal\n(10 Kb)", "Distal\n(1 Kb)",
                             "Proximal\n(100 Kb)", "Proximal\n(10 Kb)", "Proximal\n(1 Kb)"))
   return(x)
}
myFun.Calc.PPA.GF <- function(df.pred, df.obs, model = "MultiKernel") {
   
   # check the column names of the data
   checker <- !all(df.obs$ID %in% colnames(df.pred))
   if ( checker ) {
      out <- NULL
   } else {
      feat.all <- unique(df.pred$Feature)
      df.res.all <- NULL
      for (f in 1:length(feat.all)) {
         feat <- feat.all[f]
         tf <- (df.pred$Model == model) & (df.pred$Feature == feat)
         df <- df.pred[tf, ]
         df <- df[df$Group != "UseAll", ]
         df.res <- NULL
         for ( i in 1:length(unique(df$Trait)) ) {
            tr <- unique(df$Trait)[i]
            df.i <- df[df$Trait == tr, ]
            pred.mat <- t(df.i[, df.obs$ID])
            colnames(pred.mat) <- df.i$Group
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
            cheker.i <- all(df.i$Group == names(pcor.vec)) & all(df.i$Group == names(pval.vec))
            if ( cheker.i ) {
               df.res.i <- data.frame(df.i[, c("Scenario", "Trait", "Model", "Feature", "Group")],
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
   return(out)
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
cols <- GetColors(n = 13); cols.9 <- cols[c(2:4, 6:8, 10:12)]



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
df.all$Group[is.na(df.all$Group)] <- "GBLUP"
df.all <- df.all[, c("Scenario", "Feature", "Group", "Trait", "PA")]
colnames(df.all)[colnames(df.all) == "Group"] <- "Model"
df.all$Rep <- NA
df.all$Imp <- myFun.Calc.Imp(df = df.all, base.model = "GBLUP")
df.all$Trait <- myFun.ChangeTraitNames(df.all$Trait)
df.fig <- df.all[df.all$Feature != "GBLUP", ]
df.fig$Scenario <- factor(df.fig$Scenario, 
                          levels = c("From Ames to 282", "From 282 to Ames"))
#
df.fig.A <- df.fig
df.fig.A$Feature <- myFun.ChangeNames.Feat.v2(df.fig.A$Feature)
df.fig.A$Grp <- myFun.ChangeNames.Grp.v2(df.fig.A$Model)
p <- ggplot(df.fig.A, aes(x = Grp, y = Imp, colour = Trait, group = Trait))
p <- p + facet_grid(Scenario ~ Feature, scales = "free_x", space = "free_x")
p <- p + ylab("Improvement rate (%) from GBLUP")
p <- p + geom_point(size = 2) + geom_line(lwd = .5)
p <- p + geom_hline(yintercept = 0, lty = 2)
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times New Roman"))
p <- p + scale_color_manual(values = cols.9)
p <- p + xlab("Group")
ggsave(p, file = paste0(dir.save, "/Figrue_PaImp_SingleKernel_Feat.png"), width = 9, height = 5)

# ---------------------------------------------------------------------------- #
# ----- Partial correlation
# ---------------------------------------------------------------------------- #
# load prediction results
list.all <- list(
   UseGF.Ato2 = myFun.Calc.PPA.GF(df.pred = MyFun.load("RESULT/99-SummaryFiles/PostHoc_UseFeat_AmesTo282_EvalEachKernel_BackTransPredVal.csv"), df.obs = GbPheno),
   UseGF.2toA = myFun.Calc.PPA.GF(df.pred = MyFun.load("RESULT/99-SummaryFiles/PostHoc_UseFeat_282ToAmes_EvalEachKernel_BackTransPredVal.csv"), df.obs = AmesPheno)
)
df.all <- ldply(list.all, data.frame)
# save the data frame
df.save <- df.all[, c("Scenario", "Trait", "Feature", "Group", "Partial.Cor", "P.value")]
write.csv(df.save, file = paste0(dir.save, "/Table_Partial_Correlation_UseFeat.csv"), row.names = F)












# #
# df.fig.01 <- df.fig[df.fig$Scenario == "From Ames to 282", ]
# df.fig.01$Feature <- myFun.ChangeNames.Feat(df.fig.01$Feature)
# df.fig.01$Grp <- myFun.ChangeNames.Grp(df.fig.01$Model)
# p <- ggplot(df.fig.01, aes(x = Grp, y = Imp, colour = Trait, group = Trait))
# p <- p + facet_grid(~ Feature, scales = "free_x", space = "free_x")
# p <- p + ylab("Improvement rate (%) from GBLUP")
# p <- p + geom_point(size = 2) + geom_line(lwd = .5)
# p <- p + geom_hline(yintercept = 0, lty = 2)
# p <- p + theme_bw()
# p <- p + theme(text = element_text(family = "Times New Roman"))
# p <- p + scale_color_manual(values = cols.9)
# p <- p + xlab("Group")
# ggsave(p, file = paste0(dir.save, "/Figrue_PaImp_SingleKernel_Feat.png"), width = 9, height = 5)
# 
# 
# #
# df.fig.02 <- df.fig[df.fig$Scenario == "From Ames to 282", ]
# df.fig.02$Feature <- myFun.ChangeNames.Feat.v2(df.fig.02$Feature)
# df.fig.02$Grp <- myFun.ChangeNames.Grp.v2(df.fig.02$Model)
# p <- ggplot(df.fig.02, aes(x = Grp, y = Imp, colour = Trait, group = Trait))
# p <- p + facet_grid(~ Feature, scales = "free_x", space = "free_x")
# p <- p + ylab("Improvement rate (%) from GBLUP")
# p <- p + geom_point(size = 2) + geom_line(lwd = .5)
# p <- p + geom_hline(yintercept = 0, lty = 2)
# p <- p + theme_bw()
# p <- p + theme(text = element_text(family = "Times New Roman"))
# p <- p + scale_color_manual(values = cols.9)
# p <- p + xlab("Group")
# ggsave(p, file = paste0(dir.save, "/Figrue_PaImp_SingleKernel_Feat_v2.png"), width = 9, height = 5)
# 
