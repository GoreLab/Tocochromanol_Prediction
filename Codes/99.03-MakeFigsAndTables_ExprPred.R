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
myFun.ChangeModeltNames.Expr <- function(x) {
	x[x == "GBLUP"] <- "GBLUP\n(Baseline)"
	x[x == "TBLUP"] <- "TBLUP\n(All genes)"
	x[x == "TBLUP_sub"] <- "TBLUP\n(Cand. genes)"
	x[x == "GBLUP+TBLUP"] <- "GBLUP+TBLUP\n(All genes)"
	x[x == "GBLUP+TBLUP_sub"] <- "GBLUP+TBLUP\n(Cand. genes)"
	x[x == "MultiNamLargeGenes"] <- "Multi-trait model\n(Large-eff. genes)"
	x <- factor(x, levels = c("GBLUP\n(Baseline)", "TBLUP\n(All genes)", 
							  "TBLUP\n(Cand. genes)", "GBLUP+TBLUP\n(All genes)",
							  "GBLUP+TBLUP\n(Cand. genes)", "Multi-trait model\n(Large-eff. genes)"))
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
cols.Expr <- c("gray", cols[8:12])
windowsFonts(Times=windowsFont("Times New Roman"))

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
cols <- GetColors(n = 13); cols.9 <- cols[c(2:4, 6:8, 10:12)]
p <- ggplot(df.fig, aes(x = Model, y = PA, fill = Trait))
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
						   ymax = SD.bar.max), 
					   position = position_dodge(.9), width=.5)
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times"))
p <- p + ylab("Predictive ability")
p <- p + scale_fill_manual(values = cols.9)
ggsave(p, file = paste0(dir.save, "/Figure_PA_TrsctPred_ver1.png"), width = 9, height = 5)

#
p <- ggplot(df.fig, aes(x = Trait, y = PA, fill = Model))
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
						   ymax = SD.bar.max), 
					   position = position_dodge(.9), width=.5)
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times"))
p <- p + ylab("Predictive ability")
p <- p + scale_fill_manual(values = cols.Expr)
p <- p + theme(legend.key.height= unit(1, 'cm'))
ggsave(p, file = paste0(dir.save, "/Figure_PA_TrsctPred_ver2.png"), width = 12, height = 5)


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
p <- p + theme(text = element_text(family = "Times"))
p <- p + ylab("Improvement Rate (%) from GBLUP")
p <- p + scale_fill_manual(values = cols.9)
ggsave(p, file = paste0(dir.save, "/Figure_PA_Imp_TrsctPred_ver1.png"), width = 9, height = 5)

#
p <- ggplot(df.fig, aes(x = Trait, y = PA, fill = Model))
p <- p + geom_errorbar(aes(ymin = SD.bar.min, 
						   ymax = SD.bar.max), 
					   position = position_dodge(.9), width=.5)
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black")
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times"))
p <- p + ylab("Improvement Rate (%) from GBLUP")
p <- p + scale_fill_manual(values = cols.Expr[-1])
p <- p + theme(legend.key.height= unit(1, 'cm'))
ggsave(p, file = paste0(dir.save, "/Figure_PA_Imp_TrsctPred_ver2.png"), width = 12, height = 5)


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
										"MultiNamLargeGenes"))
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
mat.save <- mat.save[1:6, ] # remove one
write.csv(mat.save, file = paste0(dir.save, "/Table_PA_TrsctPred.csv"), row.names = F)

