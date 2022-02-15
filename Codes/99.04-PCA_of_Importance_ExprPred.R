# See the importance (reg. coef.) of each transcript

# Packages
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
my.ggbiplot <- function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
												 obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
												 ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
												 alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
												 varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
												 ...) 
{
	library(ggplot2)
	library(plyr)
	library(scales)
	library(grid)
	stopifnot(length(choices) == 2)
	if (inherits(pcobj, "prcomp")) {
		nobs.factor <- sqrt(nrow(pcobj$x) - 1)
		d <- pcobj$sdev
		u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
		v <- pcobj$rotation
	}
	else if (inherits(pcobj, "princomp")) {
		nobs.factor <- sqrt(pcobj$n.obs)
		d <- pcobj$sdev
		u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
		v <- pcobj$loadings
	}
	else if (inherits(pcobj, "PCA")) {
		nobs.factor <- sqrt(nrow(pcobj$call$X))
		d <- unlist(sqrt(pcobj$eig)[1])
		u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
		v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
																									1]), FUN = "/")
	}
	else if (inherits(pcobj, "lda")) {
		nobs.factor <- sqrt(pcobj$N)
		d <- pcobj$svd
		u <- predict(pcobj)$x/nobs.factor
		v <- pcobj$scaling
		d.total <- sum(d^2)
	}
	else {
		stop("Expected a object of class prcomp, princomp, PCA, or lda")
	}
	choices <- pmin(choices, ncol(u))
	df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
															FUN = "*"))
	v <- sweep(v, 2, d^var.scale, FUN = "*")
	df.v <- as.data.frame(v[, choices])
	names(df.u) <- c("xvar", "yvar")
	names(df.v) <- names(df.u)
	if (pc.biplot) {
		df.u <- df.u * nobs.factor
	}
	r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
	v.scale <- rowSums(v^2)
	df.v <- r * df.v/sqrt(max(v.scale))
	if (obs.scale == 0) {
		u.axis.labs <- paste("standardized PC", choices, 
												 sep = "")
	}
	else {
		u.axis.labs <- paste("PC", choices, sep = "")
	}
	u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
																						100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
	if (!is.null(labels)) {
		df.u$labels <- labels
	}
	if (!is.null(groups)) {
		df.u$groups <- groups
	}
	if (varname.abbrev) {
		df.v$varname <- abbreviate(rownames(v))
	}
	else {
		df.v$varname <- rownames(v)
	}
	df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
	df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
	g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
		ylab(u.axis.labs[2]) + coord_equal()
	if (var.axes) {
		if (circle) {
			theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
																								length = 50))
			circle <- data.frame(xvar = r * cos(theta), yvar = r * 
													 	sin(theta))
			g <- g + geom_path(data = circle, color = muted("white"), 
												 size = 1/2, alpha = 1/3)
		}
		g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
																					 xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 
																					 																											 "picas")), color = muted("red"))
	}
	if (!is.null(df.u$labels)) {
		if (!is.null(df.u$groups)) {
			g <- g + geom_text(aes(label = labels, color = groups), 
												 size = labels.size)
		}
		else {
			g <- g + geom_text(aes(label = labels), size = labels.size)
		}
	}
	else {
		if (!is.null(df.u$groups)) {
			g <- g + geom_point(aes(color = groups), alpha = alpha)
		}
		else {
			g <- g + geom_point(alpha = alpha)
		}
	}
	if (!is.null(df.u$groups) && ellipse) {
		theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
		circle <- cbind(cos(theta), sin(theta))
		ell <- ddply(df.u, "groups", function(x) {
			if (nrow(x) <= 2) {
				return(NULL)
			}
			sigma <- var(cbind(x$xvar, x$yvar))
			mu <- c(mean(x$xvar), mean(x$yvar))
			ed <- sqrt(qchisq(ellipse.prob, df = 2))
			data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
											 mu, FUN = "+"), groups = x$groups[1])
		})
		names(ell)[1:2] <- c("xvar", "yvar")
		g <- g + geom_path(data = ell, aes(color = groups, group = groups))
	}
	if (var.axes) {
		g <- g + geom_text(data = df.v, aes(label = varname, 
											x = xvar, y = yvar, angle = angle, hjust = hjust), 
						   color = "darkred", size = varname.size, family = "Times")
	}
	return(g)
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
df.gene.13$Gene[df.gene.13$`Gene ID` == "Zm00001d014737"] <- "arodeH2-b" # manually re-name


# ---------------------------------------------------------------------------- #
# ----- PCA using regression coefficients
# ---------------------------------------------------------------------------- #
# mkdir
dir.save.2 <- paste0(dir.save, "/PCA")
dir.create(dir.save.2)

# setup
trait.all <- c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3",
							 "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")
TraitSets <- c("UseAll", "UseSix")
Models <- c("TBLUP.AllGenes", "GTBLUP.AllGenes", "TBLUP.CandGenes", "GTBLUP.CandGenes")
df.all.param <- expand.grid("TraitSet" = TraitSets, "Model" = Models, "Path"= NA, stringsAsFactors = F)
df.all.param$Path[df.all.param$Model == "TBLUP.AllGenes"] <-  "RESULT/4.2-GenExpFit_trans/Coef_ExpAll_"
df.all.param$Path[df.all.param$Model == "GTBLUP.AllGenes"] <-  "RESULT/4.2-GenExpFit_trans/Coef_GBLUP+ExpAll_"
df.all.param$Path[df.all.param$Model == "TBLUP.CandGenes"] <-  "RESULT/4.2-GenExpFit_trans/Coef_ExpCand_"
df.all.param$Path[df.all.param$Model == "GTBLUP.CandGenes"] <-  "RESULT/4.2-GenExpFit_trans/Coef_GBLUP+ExpCand_"

# loop to make figures
for ( k in 1:nrow(df.all.param) ) {
	# set seed (PCA have a ranom choice of direction)
	set.seed(123)
	
	# k-th setup
	traitset <- df.all.param$TraitSet[k]
	model <- df.all.param$Model[k]
	path <- df.all.param$Path[k]
	fig.out <- paste0(dir.save.2, "/PCA_biplot_", model, "_", traitset, ".png")
	
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
	
	# use 9 traits or 6 traits?
	if ( traitset == "UseAll" ) { M <- coef.mat }
	if ( traitset == "UseSix" ) { M <- coef.mat[, c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3")] }
	
	# PCA
	colnames(M) <- as.character(myFun.ChangeTraitNames(colnames(M)))
	pr.res <- prcomp(M, center = TRUE, scale. = TRUE)
	
	# PC1 and PC2 for the 13 genes
	m <- match(df.gene.13$`Gene ID`, rownames(pr.res$x))
	df.fig <- cbind(df.gene.13, 
									"PC1" = scale(pr.res$x[, 1])[m], 
									"PC2" = scale(pr.res$x[, 2])[m])
	p <- my.ggbiplot(pr.res)
	p <- p + geom_text_repel(df.fig, 
													 mapping = aes(x = PC1, y = PC2, label = Gene),
													 family = "Times", fontface = "italic", color = "red")
	p <- p + geom_point(df.fig, 
											mapping = aes(x = PC1, y = PC2), color = "red")
	p <- p + theme(legend.position = "none")
	p <- p + theme(text = element_text(family = "Times"))
	p$layers <- c(p$layers[[2]], p$layers[[1]], p$layers[[3]], p$layers[[4]], p$layers[[5]])
	ggsave(filename = fig.out, p, width = 10, height = 6)
}

