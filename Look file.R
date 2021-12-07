dat <- data.table::fread("Glb1_f300k.hmp_formatted.txt", data.table = F)
y <- rnorm(nrow(dat)) # this should be your own phenotypic values
pval.all <- c()
for ( i in 2:ncol(dat)) {
	x <- as.factor(dat[, i])
	y.wo.na <- y[!is.na(x)]
	x.wo.na <- x[!is.na(x)]
	m0 <- lm(y.wo.na ~ 1)
	m1 <- lm(y.wo.na ~ x.wo.na)
	pval <- anova(m0, m1)[2, "Pr(>F)"]
	pval.all <- c(pval.all, pval)
}


