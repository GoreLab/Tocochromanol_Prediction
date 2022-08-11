# use sommer to estimate covariance matrices & use MME
myFun.sommer.new <- function(PhenoData, test, K) {
	# make training data
	PhenoData.train <- PhenoData
	train <- setdiff(PhenoData.train$ID, test)
	PhenoData.train <- PhenoData.train[match(train, PhenoData.train$ID), ]
	K.train <- K[train, train]
	
	# fit the model
	PhenoData.train$ID <- factor(PhenoData.train$ID)
	mod <- mmer(cbind(Y1, Y2) ~ 1,
							random =~ vs(ID, Gu = K.train, Gtc = unsm(2)),
							rcov =~ vs(units, Gtc = unsm(2)),
							data = PhenoData.train, 
							method = "EMMA",
							verbose = TRUE, 
							date.warning = FALSE,
							tolparinv = 1)
	
	# predicted values
	y1 <- PhenoData$Y1
	if ( sum(is.na(PhenoData.train$Y2)) > 0 ) {
		num.na.Y2 <- which(is.na(PhenoData.train$Y2))
		y2 <- PhenoData.train$Y2[-num.na.Y2]
		m <- match(PhenoData.train$ID[-num.na.Y2], PhenoData$ID)
	} else {
		y2 <- PhenoData.train$Y2
		m <- match(PhenoData.train$ID, PhenoData$ID)
	}
	y <- c(y1, y2) # y vector
	Sigma.G <- mod$sigma$`u:ID`
	G <- Sigma.G %x% K # G matrix
	Z1 <- diag(nrow(PhenoData))
	Z2 <- Z1[m, ]
	Z3 <- rbind(cbind(Z1, matrix(0, nr = nrow(Z1), nc = ncol(Z1))),
							cbind(matrix(0, nr = nrow(Z2), nc = ncol(Z2)), Z2))
	Z <- Z3 # Z matrix
	Sigma.R <- mod$sigma$`u:units`
	sigma.e <- rep(diag(Sigma.R), c(length(y1), length(y2)))
	R <- diag(sigma.e) # R matrix
	X <- matrix(0, nr = length(y), nc = 2)
	X[1:length(y1), 1] <- 1
	X[(length(y1)+1):nrow(X), 2] <- 1  # X matrix
	R.inv <- tryCatch( R.inv <- MASS::ginv(R), error=function(e){} )
	G.inv <- tryCatch( G.inv <- MASS::ginv(G), error=function(e){} )
	tf <- !is.null(R.inv) & !is.null(G.inv)
	if ( tf ) {
		M11 <- t(X) %*% R.inv %*% X
		M12 <- t(X) %*% R.inv %*% Z
		M21 <- t(Z) %*% R.inv %*% X
		M22 <- t(Z) %*% R.inv %*% Z + G.inv
		M <- rbind(cbind(M11, M12), cbind(M21, M22)) # left-side matrix
		Y1 <- t(X) %*% R.inv %*% y
		Y2 <- t(Z) %*% R.inv %*% y
		Y <- c(Y1, Y2) # right-side vector
		MME.sol <- MASS::ginv(M) %*% Y
		est.beta <- MME.sol[1:2]
		est.u <- MME.sol[3:length(MME.sol)]
		u.mat <- matrix(est.u, nc = 2)
		beta.mat <- matrix(rep(est.beta, each = nrow(u.mat)), nc = 2)
		yPred.mat <- beta.mat + u.mat
		rownames(yPred.mat) <- PhenoData$ID
		y.pred <- yPred.mat[test, 2]
	} else {
		y.pred <- setNames(rep(NA, length(test)), test)
	}
	# return
	return(y.pred)
}

# use sommer: this is not an optimal way
myFun.sommer.old <- function(PhenoData, test, K) {
	# make training data
	PhenoData.train <- PhenoData
	train <- setdiff(PhenoData.train$ID, test)
	PhenoData.train$Y2[PhenoData.train$ID %in% test] <- NA

	# fit the model
	PhenoData.train$ID <- factor(PhenoData.train$ID)
	mod <- mmer(cbind(Y1, Y2) ~ 1,
							random =~ vs(ID, Gu = K, Gtc = unsm(2)),
							rcov =~ vs(units, Gtc = unsm(2)),
							data = PhenoData.train, 
							method = "EMMA",
							verbose = TRUE, 
							date.warning = FALSE,
							tolparinv = 1)
	y.pred <- mod$Beta$Estimate[2] + mod$U$`u:ID`$Y2[test]
	
	# return
	return(y.pred)
}

# use MTM
myFun.MTM <- function(PhenoData, test, K, dir.log = "LOGFILE") {
	# make folder to save log files
	dir.create(dir.log, recursive = TRUE)
	
	# make training data
	PhenoData.train <- PhenoData
	Y.train <- as.matrix(PhenoData.train[, -1])
	rownames(Y.train) <- PhenoData.train$ID
	Y.train[test, 2] <- NA # mask test data
	
	# fit the model
	fm <- MTM(Y = Y.train,
						K = list(list(K = K, COV = list(type = 'UN', df0 = 2, S0 = diag(2)))),
						resCov = list(type = 'DIAG', S0 = rep(1, 2), df0 = rep(1, 2)),
						nIter = 6000, burnIn = 4000, thin = 5, saveAt = paste0(dir.log, "/fm_"))
	y.pred <- fm$YHat[, 2][test]
	
	# return
	return(y.pred)
}

# use BME
myFun.BME <- function(PhenoData, test, K) {
	# re-order
	PhenoData$ID <- as.factor(PhenoData$ID)
	o <- order(PhenoData$ID)
	PhenoData <- PhenoData[o, ]
	K <- K[o, o]
	
	# mask test data
	PhenoData.train <- PhenoData
	PhenoData.train$Y2[PhenoData.train$ID %in% test] <- NA
	
	# run BME
	LG <- cholesky(K)
	ZG <- model.matrix(~ 0 + as.factor(PhenoData.train$ID))
	Z.G <- ZG %*% LG
	Y <- as.matrix(PhenoData.train[, -1])
	fm <- BME(Y = Y, Z1 = Z.G, nIter = 6000, burnIn = 4000, thin = 5, bs = 50)
	yPred <- fm$yHat[, 2]
	names(yPred) <- PhenoData.train$ID
	y.pred <- yPred[test]
	
	# return
	return(y.pred)
}

# GBLUP
myFun.GBLUP <- function (PhenoData, test, K, dir.log = "LOGFILE") {
	# make input vector
	y <- PhenoData$Y2
	names(y) <- PhenoData$ID
	y[test] <- NA # mask phenotype
	
	# model 
	covar.mat <- matrix(PhenoData$Y1, nc = 1)
	ETA <- list("G" = list(K = K, model = "RKHS"),
							"Covariate" = list(X = covar.mat, model = "FIXED"))
	
	# fitting via BGLR
	fm <- BGLR(y = y, ETA = ETA,
						 nIter = 12000, burnIn = 8000, thin = 5,
						 saveAt = paste0(dir.log, "/fm_"), verbose = FALSE)
	y.pred <- fm$yHat[test]
	
	# return
	return(y.pred)
}

# GBLUP with covariate
myFun.GBLUP.v2 <- function (PhenoVec, test, K, covar.mat, dir.log = "LOGFILE") {
	# make input vector
	PhenoVec[test] <- NA # mask phenotype
	
	# model 
	ETA <- list("G" = list(K = K, model = "RKHS"),
							"Covariate" = list(X = covar.mat, model = "FIXED"))
	
	# fitting via BGLR
	fm <- BGLR(y = PhenoVec, ETA = ETA,
						 nIter = 12000, burnIn = 8000, thin = 5,
						 saveAt = paste0(dir.log, "/fm_"), verbose = FALSE)
	y.pred <- fm$yHat[test]
	
	# return
	return(y.pred)
}

# use MTM with multiple genes
myFun.MTM.v2 <- function(PhenoData, test, K, dir.log = "LOGFILE") {
	# make folder to save log files
	dir.create(dir.log, recursive = TRUE)
	
	# make training data
	PhenoData.train <- PhenoData
	Y.train <- as.matrix(PhenoData.train[, -1])
	rownames(Y.train) <- PhenoData.train$ID
	Y.train[test, "Y"] <- NA # mask test data
	
	# fit the model
	p <- ncol(Y.train)
	fm <- MTM(Y = Y.train,
						K = list(list(K = K, COV = list(type = 'UN', df0 = p, S0 = diag(p)))),
						resCov = list(type = 'DIAG', S0 = rep(1, p), df0 = rep(1, p)),
						nIter = 12000, burnIn = 8000, thin = 5, saveAt = paste0(dir.log, "/fm_"))
	y.pred <- fm$YHat[, "Y"][test]
	
	# return
	return(y.pred)
}

