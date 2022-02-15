# See the importance (reg. coef.) of each transcript

# Packages
# library(devtools)
# install_github("vqv/ggbiplot")
library(ggplot2)
library(ggbiplot)
library(readxl)
library(reshape2)
library(data.table)
library(extrafont)
library(inlmisc)
library(plyr)
library(car)
library(ggrepel)
library(scales)
library(ggpubr)
library(RColorBrewer)

# setup for figure
loadfonts(device = "win")
windowsFonts(Times=windowsFont("Times New Roman"))

# functions
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

# mkdir
dir.save <- "RESULT/99-Tables_and_Figures"
dir.create(dir.save, recursive = TRUE)

# setup for figure
cols <- GetColors(n = 13)
cols.9 <- cols[c(2:4, 6:8, 10:12)]
loadfonts(device = "win")
windowsFonts(Times = windowsFont("Times New Roman"))

# list of the 126 genes
list.cand <- as.data.frame(read_xlsx("RAWDATA/tocochromanol_all_candidate_genes_combined_DW_20190624_with_two_genes.xlsx"))
id.126 <- unique(list.cand$RefGen_v4.Gene.ID)
df.cand <- data.frame("id" = id.126)
df.cand$name <- list.cand$Gene.Name[match(df.cand$id, list.cand$RefGen_v4.Gene.ID)]

# show the two arodeH2
df.cand[df.cand$id %in% c("Zm00001d014734", "Zm00001d014737"), ]

# list of the 13 genes
df.gene.13 <- as.data.frame(read_xlsx("RAWDATA/CandidateGeneNames.xlsx"))
df.gene.13$Gene[df.gene.13$`Gene ID` == "Zm00001d014737"] <- "arodeH2-b" # manually curate

# get PVE etc
dat.MKGBLUP <- read.csv("RESULT/99-Tables_and_Figures/Table_Partial_Correlation_and_PVE.csv")
dat.NAM <- read.csv("RAWDATA/NAM_QTL/qtl.data.from.di.csv")
dat.NAM.2 <- as.data.frame(read_xlsx("RAWDATA/NAM_QTL/TocoGene.xlsx"))

# manually confirm that samt1 is in 29th NAM JL-QTL
df.cand[df.cand$id == "Zm00001d017937", ]
list.cand[list.cand$RefGen_v4.Gene.ID == "Zm00001d017937", ]
qtl.id.tmp <- unique(dat.NAM$QTL.ID[dat.NAM$Chr == "Chr5"])
dat.NAM.tmp <- dat.NAM[match(qtl.id.tmp, dat.NAM$QTL.ID), ]
dat.NAM.tmp[(dat.NAM.tmp$Common.SupInt.L.Pos.v4 < 210385310) & (210401948 < dat.NAM.tmp$Common.SupInt.R.Pos.v4), ]
dat.NAM[dat.NAM$QTL.ID == 29, ]

# update for samt1
dat.NAM$NAME[dat.NAM$NAME == "QTL29"] <- "QTL29(samt1)"
dat.MKGBLUP$NAME[dat.MKGBLUP$NAME == "QTL29"] <- "QTL29(samt1)"


# ---------------------------------------------------------------------------- #
# ----- Compare NAM and TBLUP 
# ---------------------------------------------------------------------------- #
# mkdir
dir.save.2 <- paste0(dir.save, "/Importance")
dir.create(dir.save.2)

# make a table
trait.all <- c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3",
							 "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")
qtl.all <- paste0("QTL", unique(dat.NAM$QTL.ID))
df.all <- expand.grid("Trait" = trait.all, 
											"QTL" = qtl.all,
											stringsAsFactors = F)
df.all$GENE <- df.all$PVE <- df.all$QTL.NAME <- NA
for ( i in 1:nrow(df.all) ) {
	trait <- df.all$Trait[i]
	qtl <- df.all$QTL[i]
	qtl.num <- as.numeric(gsub("QTL", "", qtl))
	tf <- (dat.NAM$Trait == trait) & (dat.NAM$QTL.ID == qtl.num)
	if ( sum(tf) == 1 ) {
		df.all$PVE[i] <- dat.NAM$PVE[tf]
		df.all$QTL.NAME[i] <- dat.NAM$NAME[tf]
		if ( length(grep("\\(", dat.NAM$NAME[tf])) == 1 ) {
			df.all$GENE[i] <-gsub(")", "", strsplit(dat.NAM$NAME[tf], "\\(")[[1]][2])
		}
	}
}
qtl.resolved <- unique(df.all$QTL[!is.na(df.all$GENE)])
df.NAM.resolved <- df.all[df.all$QTL %in% qtl.resolved, ]
for ( i in 1:nrow(df.NAM.resolved) ) {
	df.NAM.resolved$GENE[i] <- setdiff(unique(df.NAM.resolved$GENE[df.NAM.resolved$QTL %in% df.NAM.resolved$QTL[i]]), NA)
	df.NAM.resolved$QTL.NAME[i] <- setdiff(unique(df.NAM.resolved$QTL.NAME[df.NAM.resolved$QTL %in% df.NAM.resolved$QTL[i]]), NA)
}
df.NAM.resolved <- df.NAM.resolved[, c("GENE", "QTL", "QTL.NAME", "Trait", "PVE")]

# get gene id etc.
df.NAM.resolved$GENE.ID.v2 <- NA
df.NAM.resolved$GENE.ID.v4 <- NA
df.NAM.resolved$eQTL <- NA
for ( i in 1:nrow(df.NAM.resolved) ) {
	gene <- df.NAM.resolved$GENE[i]
	if ( gene != "samt1" ) {
		tf <- dat.NAM.2$Gene %in% gene
		df.NAM.resolved$GENE.ID.v2[i] <- dat.NAM.2$RefGen_v2[tf]
		df.NAM.resolved$GENE.ID.v4[i] <- dat.NAM.2$RefGen_v4[tf]
		df.NAM.resolved$eQTL[i] <- dat.NAM.2$eQTL[tf]
	} else {
		df.NAM.resolved$GENE.ID.v2[i] <- NA
		df.NAM.resolved$GENE.ID.v4[i] <- "Zm00001d017937"
		df.NAM.resolved$eQTL[i] <- "Yes" # check later if this is fine
	}
}

# get partial correlation in MK-GBLUP
df.NAM.resolved$p.cor.from.Ames.to.282 <- NA
df.NAM.resolved$p.cor.from.282.to.Ames <- NA
for ( i in 1:nrow(df.NAM.resolved) ) {
	qtl.name <- df.NAM.resolved$QTL.NAME[i]
	trait <- df.NAM.resolved$Trait[i]
	tf <- (dat.MKGBLUP$NAME %in% qtl.name) & (dat.MKGBLUP$Trait %in% trait)
	if ( sum(tf) != 0 ) {
		dat <- dat.MKGBLUP[tf, ]
		dat <- dat[dat$Window == "SI", ]
		df.NAM.resolved$p.cor.from.Ames.to.282[i] <- dat$Partial.Cor[dat$Scenario == "From Ames to 282"]
		df.NAM.resolved$p.cor.from.282.to.Ames[i] <- dat$Partial.Cor[dat$Scenario == "From 282 to Ames"]
	}
}

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

# convert the regression coefficients to importance %
coef.rank.mat <- coef.mat
for ( k in 1:ncol(coef.rank.mat) ) {
	coef.vec <- coef.rank.mat[, k]
	importance.vec <- abs(coef.vec)
	importance.rank.vec <- rank((-1) * importance.vec)
	coef.rank.mat[, k] <- importance.rank.vec
}

# add result to the data frame
df.NAM.resolved$Regression.Coef <- NA
df.NAM.resolved$Importance.rank <- NA
for ( i in 1:nrow(df.NAM.resolved) ) {
	gene.id <- df.NAM.resolved$GENE.ID.v4[i]
	trait <- df.NAM.resolved$Trait[i]
	if ( gene.id %in% rownames(coef.rank.mat) ) {
		df.NAM.resolved$Importance.rank[i] <- coef.rank.mat[gene.id, trait]
		df.NAM.resolved$Regression.Coef[i] <- coef.mat[gene.id, trait]
	}
}

# save the table
write.csv(df.NAM.resolved, 
					file = paste0(dir.save.2, "/SummaryTable.csv"), 
					row.names = F)

# ---------------------------------------------------------------------------- #
# ----- Make figure 2
# ---------------------------------------------------------------------------- #
# setup the data frame
df.fig <- df.NAM.resolved[, c("GENE", "Trait", "PVE", "eQTL",
							  "p.cor.from.Ames.to.282", "p.cor.from.282.to.Ames", "Importance.rank")]
colnames(df.fig)[colnames(df.fig) == "GENE"] <- "Gene"
colnames(df.fig)[colnames(df.fig) == "p.cor.from.Ames.to.282"] <- "MKGBLUP.Ames.to.282"
colnames(df.fig)[colnames(df.fig) == "p.cor.from.282.to.Ames"] <- "MKGBLUP.282.to.Ames"
colnames(df.fig)[colnames(df.fig) == "Importance.rank"] <- "GTBLUP.Importance.Rank"

# group by PVE
df.fig$PVE.group <- NA
df.fig$PVE.group[is.na(df.fig$PVE)] <- "n.s."
df.fig$PVE.group[!is.na(df.fig$PVE) & (df.fig$PVE < 0.05)] <- "PVE < 5%"
df.fig$PVE.group[!is.na(df.fig$PVE) & (df.fig$PVE > 0.05)] <- "PVE > 5%"
df.fig$PVE.group <- factor(df.fig$PVE.group,
						   levels = c("n.s.", "PVE < 5%", "PVE > 5%"))

# group by importance
df.fig$GTBLUP.Importance.Rank.group <- NA 
df.fig$GTBLUP.Importance.Rank.group[is.na(df.fig$GTBLUP.Importance.Rank)] <- "n.a."
df.fig$GTBLUP.Importance.Rank.group[!is.na(df.fig$GTBLUP.Importance.Rank) & (df.fig$GTBLUP.Importance.Rank <= floor(nrow(coef.rank.mat) * 0.10))] <- "top 10%"
df.fig$GTBLUP.Importance.Rank.group[!is.na(df.fig$GTBLUP.Importance.Rank) & (df.fig$GTBLUP.Importance.Rank <= floor(nrow(coef.rank.mat) * 0.01))] <- "top 1%"
df.fig$GTBLUP.Importance.Rank.group[!is.na(df.fig$GTBLUP.Importance.Rank) & (df.fig$GTBLUP.Importance.Rank > floor(nrow(coef.rank.mat) * 0.10))] <- "n.s."
df.fig$GTBLUP.Importance.Rank.group <- factor(df.fig$GTBLUP.Importance.Rank.group,
													 levels = c("n.a.", "n.s.", "top 10%", "top 1%"))

# change trait name
df.fig$Trait <- myFun.ChangeTraitNames(df.fig$Trait)

# change gene order
df.fig$Gene <- factor(df.fig$Gene,
					  levels = rev(c("arodeH2", "dxs2", "hggt1", "sds2", "vte3", "hppd1", "vte4",
					  		   "por1", "por2", "snare", "ltp", "phd", "fbn", "samt1", "vte7")))

# make a figure
df.fig.PVE <- df.fig
df.fig.PVE$pesudo.class <- "PVE evaluated in NAM"
p <- ggplot(df.fig.PVE, aes(x = Trait, y = Gene, fill = PVE.group))
p <- p + facet_wrap(~ pesudo.class)
p <- p + geom_tile(colour = "black", size = 1, width = 0.9, height = 0.9)
p <- p + scale_fill_manual(values = c("white", brewer.pal(5, "YlOrRd")[c(2, 5)]))
p <- p + theme_minimal()
p <- p + theme(panel.grid = element_line(color = "white"))
p <- p + theme(legend.position = "bottom")
p <- p + theme(axis.text.y = element_text(face="italic"))
p <- p + theme(text = element_text(family = "Times"))
p <- p + guides(fill=guide_legend(title="PVE group"))
p <- p + theme(strip.text.x = element_text(size = 15))
ggsave(filename = paste0(dir.save.2, "/Tile_PVE.png"), p, width = 5, height = 8)

# make a figure
df.fig.Imp <- df.fig
df.fig.Imp$pesudo.class <- "Importance in GBLUP+TBLUP"
p <- ggplot(df.fig.Imp, aes(x = Trait, y = Gene, fill = GTBLUP.Importance.Rank.group))
p <- p + facet_wrap(~ pesudo.class)
p <- p + geom_tile(colour = "black", size = 1, width = 0.9, height = 0.9)
p <- p + scale_fill_manual(values = c("white", brewer.pal(5, "YlOrRd")[c(1, 3, 5)]))
p <- p + theme_minimal()
p <- p + theme(panel.grid = element_line(color = "white"))
p <- p + theme(legend.position = "bottom")
p <- p + theme(axis.text.y = element_text(face="italic"))
p <- p + theme(text = element_text(family = "Times"))
p <- p + theme(strip.text.x = element_text(size = 15))
p <- p + guides(fill=guide_legend(title="Importance ranking"))
ggsave(filename = paste0(dir.save.2, "/Tile_Importance.png"), p, width = 5, height = 8)

# make a figure
df.fig.eQTL <- aggregate(eQTL ~ Gene, df.fig, unique)
df.fig.eQTL$x <- "eQTL"
df.fig.eQTL$pesudo.class <- "eQTL"
p <- ggplot(df.fig.eQTL, aes(x = x, y = Gene, fill = eQTL))
p <- p + facet_wrap(~ pesudo.class)
p <- p + geom_tile(colour = "black", size = 1, width = 0.9, height = 0.9)
p <- p + scale_fill_manual(values = c("white", brewer.pal(5, "YlOrRd")[5]))
p <- p + theme_minimal()
p <- p + xlab("")
p <- p + theme(panel.grid = element_line(color = "white"))
p <- p + theme(legend.position = "bottom")
p <- p + theme(axis.text.y = element_text(face="italic"))
p <- p + theme(text = element_text(family = "Times"))
p <- p + theme(strip.text.x = element_text(size = 15))
ggsave(filename = paste0(dir.save.2, "/Tile_ceeQTL.png"), p, width = 3, height = 8)


# # make a figure
# p <- ggplot(df.fig, aes(x = Trait, y = Gene, fill = MKGBLUP.Ames.to.282))
# p <- p + geom_tile(colour = "black", size = 1, width = 0.9, height = 0.9)
# p <- p + scale_fill_gradient(low = "yellow", high = "red", na.value = "white")
# p <- p + theme_minimal()
# p <- p + theme(panel.grid = element_line(color = "white"))
# p <- p + theme(legend.position = "bottom")
# ggsave(filename = paste0(dir.save.2, "/TmpFig2.png"), p, width = 5, height = 8)
# 
# # make a figure
# p <- ggplot(df.fig, aes(x = Trait, y = Gene, fill = MKGBLUP.282.to.Ames))
# p <- p + geom_tile(colour = "black", size = 1, width = 0.9, height = 0.9)
# p <- p + scale_fill_gradient(low = "yellow", high = "red", na.value = "white")
# p <- p + theme_minimal()
# p <- p + theme(panel.grid = element_line(color = "white"))
# p <- p + theme(legend.position = "bottom")
# ggsave(filename = paste0(dir.save.2, "/TmpFig3.png"), p, width = 5, height = 8)



# hist(df.fig$Importance.rank[is.na(df.fig$PVE)], breaks = 10)
# hist(df.fig$Importance.rank[!is.na(df.fig$PVE)], breaks = 10)
# hist(df.fig$Importance.rank[!is.na(df.fig$PVE) & (df.fig$PVE > 0.05)], breaks = 10)
# 
# # 
# observed.rank <- sort(df.fig$Importance.rank[is.na(df.fig$PVE)])
# expected.rank <- seq(from = 1, to = nrow(coef.mat), length.out = length(observed.rank))
# plot(x = observed.rank, y = expected.rank, 
# 		 xlim = c(0, nrow(coef.mat)), ylim = c(0, nrow(coef.mat)),
# 		 pch = 20)
# abline(0, 1, lty = 2)
# 
# 
# # 
# observed.rank <- sort(df.fig$Importance.rank[!is.na(df.fig$PVE)])
# expected.rank <- seq(from = 1, to = nrow(coef.mat), length.out = length(observed.rank))
# plot(x = observed.rank, y = expected.rank, 
# 		 xlim = c(0, nrow(coef.mat)), ylim = c(0, nrow(coef.mat)),
# 		 pch = 20)
# abline(0, 1, lty = 2)

#




# # ---------------------------------------------------------------------------- #
# # ----- Make figure 1
# # ---------------------------------------------------------------------------- #
# df.fig <- df.NAM.resolved[!is.na(df.NAM.resolved$Importance.rank), ]
# 
# # e.g.)
# df.fig$Group.based.on.NAM.PVE <- NA
# df.fig$Group.based.on.NAM.PVE[is.na(df.fig$PVE)] <- "not significant in NAM\n(n = 48 gene-trait pairs)"
# df.fig$Group.based.on.NAM.PVE[!is.na(df.fig$PVE) & (df.fig$PVE < 0.05)] <- "weakly significant in NAM (PVE < 5%)\n(n = 30 gene-trait pairs)"
# df.fig$Group.based.on.NAM.PVE[!is.na(df.fig$PVE) & (df.fig$PVE > 0.05)] <- "significant in NAM (PVE > 5%)\n(n = 21 gene-trait pairs)"
# df.fig$Group.based.on.NAM.PVE <- factor(df.fig$Group.based.on.NAM.PVE,
# 																				levels = c("not significant in NAM\n(n = 48 gene-trait pairs)",
# 																									 "weakly significant in NAM (PVE < 5%)\n(n = 30 gene-trait pairs)",
# 																									 "significant in NAM (PVE > 5%)\n(n = 21 gene-trait pairs)"))
# p <- ggplot(df.fig, aes(x = Group.based.on.NAM.PVE, y = Importance.rank))
# p <- p + geom_dotplot(binaxis = "y", stackdir = "center")
# p <- p + ggtitle("Relationship between the importance in GBLUP+TBLUP and NAM JL-QTL result")
# ggsave(p, filename = paste0(dir.save.2, "/FigX.png"), width = 8, height = 5)
