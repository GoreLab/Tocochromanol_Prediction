#### NOTE BY RYOEI TANAKA
#### This is a code From Di, originally written by the previous postdoc. Original name of the R code was 's1_boxcox_function.R'
#### I did not change the calculation, but I modified the code so that I can specify the folder/filename to save the plots.



# boxcox_function

#library(asreml)
library(MASS)
boxcox.fitting <- function(pheno, curr.trait, filename1 = NULL, filename2 = NULL) {
    
#########################################################################################
# Apply the Box-Cox transformation for each trait anc check if lambda is different from 1
# If lambda is not different from 1, then the trait is unchanged
# if it is, then both the untransformed and transformed values are output.

    #### REMOVE ZEROS and NEGATIVE VALUES !!! Zeros or negative values are not defined for some of the transformations
    # pheno <- mer_t
    # pheno[ which(pheno[curr.trait] <= 0) , curr.trait] <- NA # Test with logF

    #########################

    ## rmtrait <- which(is.na(eval(parse(text=paste('pheno$', curr.trait, sep='')))))
    ## if ( length(rmtrait > 0) ) {
    ##     tmppheno <- pheno[-rmtrait,]
    ## } else {
    ##     tmppheno <- pheno
    ## }
    ## rmgeno <- which(is.na(eval(parse(text=paste('tmppheno$', PREDICTION, sep='')))))
    ## if ( length(rmgeno > 0) ) {
    ##     tmppheno <- tmppheno[-rmgeno,]
    ## }
    
    # ---------------------- begin: written by RT ---------------------- #
    if ( is.null(filename1) ) { filename.pdf1 <- paste(curr.trait, "_lambda_values.jpeg", sep="") }
    if ( is.null(filename2) ) { filename.pdf1 <- paste(curr.trait,"_Distrib_After_Box_plot.jpeg",sep="") }
    # ----------------------   end: written by RT ---------------------- #

    # Run the Box-Cox procedure
    bxc <- paste( "boxcox(", curr.trait, "~", "1", ",data = pheno, lambda = seq(-2, 2, .05))") ### modified by Xiaowei Li
    jpeg(filename1) ### modified by RT
    the.boxcox <- eval(parse( text = bxc))
    # extract the value for lambda
    the.lambda <- with(the.boxcox, x[which.max(y)])

    # construct confidence 95% confidence interval
    the.interval <- range( the.boxcox$x[the.boxcox$y > max(the.boxcox$y) - qchisq(0.95,1)/2] )
    interval.test <- function(x, low, hi) (x >= low & x <= hi)
    # Check if 1 is in the interval
    is.one.in.interval <- interval.test(1, the.interval[1], the.interval[2])
    if (is.one.in.interval == TRUE) {
        #newvalue <- pheno[curr.trait]
        lambda.to.use <- 1
        abline(v=1, col='red')
    # If not, run the boxcox for the next two adjacent convenient values and take the one that maximizes the likelihood
    } else {
        convenient.lambdas <- c(-2, -1, -0.5, 0, 0.5, 1, 2)
        diff = abs(convenient.lambdas - the.lambda)
        # get the 2 smallest values
        aaa = sort(diff,partial=2)[1:2]
        idx.smallest.diff = which( diff %in% aaa)
        convenient.lambdas.to.test = convenient.lambdas[idx.smallest.diff]
        bxc1 <- paste( "boxcox(", curr.trait, "~", "1", ", plotit = FALSE, data = pheno, lambda = seq(", paste(convenient.lambdas.to.test[1],",",convenient.lambdas.to.test[2]), ", length=2))")
        the.boxcox1 <- eval(parse( text = bxc1))
        # get the closest.convenient.lambda
        lambda.to.use <- with(the.boxcox1, x[which.max(y)])
        lambda.not.use <- with(the.boxcox1, x[which.min(y)])
        abline(v=lambda.to.use, col='red')
        abline(h=the.boxcox1$y[which(the.boxcox1$x == lambda.to.use)], col='red')
        abline(v=lambda.not.use, col='blue', lty=2, lwd=0.5)
        abline(h=the.boxcox1$y[which(the.boxcox1$x == lambda.not.use)], col='blue', lty=2, lwd=0.5)
        title(main='Blue; convenient lambda not used -- Red: convenient lambda used', cex.main=1)
    }
    dev.off()

    if(lambda.to.use == -2) {
        newvalue <- 1/(pheno[curr.trait]^2)
    } else if (lambda.to.use == -1) {
        newvalue <- 1/(pheno[curr.trait])
    } else if( lambda.to.use == -0.5) {
        newvalue <- 1/(sqrt(pheno[curr.trait]))
    } else if (lambda.to.use == 0) {
        newvalue <- log(pheno[curr.trait])
    } else if (lambda.to.use == 0.5) {
        newvalue <- sqrt(pheno[curr.trait])
    # this is for the case where 1 is not in the CI but it is the closest convenient value
    } else if (lambda.to.use == 1) {
        newvalue <- pheno[curr.trait]
    } else { # (lambda.to.use == 2)
        newvalue <- pheno[curr.trait]^2
    }

    # plot the transformed (and non-transform) trait distribution
    jpeg(filename2 , width=720, height=560) ### modified by RT
    opar = par() 
    layout(matrix(c(1:4), nrow=2, ncol=2),1,1)
    newvalue <- as.matrix(newvalue)
    # boxplot
    boxplot(newvalue, main = paste('best lambda=', round(the.lambda,3), '; convenient lambda used=', lambda.to.use, sep="") )
    # qqplots for data
    qqnorm(newvalue)
    qqline(newvalue, col='red', lwd=2, lty=2)
    # histogram view
    hist(as.matrix(pheno[curr.trait]), breaks=20, freq=FALSE, main='Old distribution')
    curve(dnorm(x, mean=mean(as.matrix(pheno[curr.trait]), na.rm=TRUE ), sd=sd(as.matrix(pheno[curr.trait]), na.rm=TRUE)), add=TRUE, col="green", lwd=2)
    hist(newvalue, breaks=20, freq=FALSE, xlim= c(min(newvalue, na.rm=TRUE), max(newvalue, na.rm=TRUE)), main = 'New distribution')
    curve(dnorm(x, mean=mean(newvalue, na.rm=TRUE), sd=sd(newvalue, na.rm=TRUE)), add=TRUE, col="green", lwd=2)
    par(opar)
    dev.off()

    ### To write the values in file, remember that we have remove the NAs at the begining.

    # newvalue[which(is.na(pheno[curr.trait]))] <- NA
    return(list(TRT=newvalue, LAMBDA=lambda.to.use))
}
