#' rigorous Lasso for Linear Models: Inference
#'
#' Estimation and inference of selected (low-dimensional) coefficients in a high-dimensional linear model.
#'
#' The functions estimates selected (low-dimensional) coefficients in a high-dimensional linear model.
#' An application is e.g. estimation of a treatment effect \eqn{\alpha_0} in a
#' setting of high-dimensional controls. The so-called post-double-selection
#' estimator is implemented. The idea is to select variable by regression of
#' the outcome variable on the control variables and the treatment variable on
#' the control variables. The final estimation is done by a regression of the
#' outcome on the treatment effect and the union of the selected variables in
#' the first two steps. The resulting estimator for \eqn{\alpha_0} is normal
#' distributed which allows inference on the treatment effect.
#'
#' @param x matrix of regressor variables serving as controls and potential
#' treatments
#' @param y outcome variable (vector or matrix)
#' @param formula an element of class \code{formula}
#' @param data data.frame in connection with \code{formula}
#' @param index vector of integers, logicals or variables names indicating the position (column) of
#' variables (integer case), logical yes or no (TRUE or FALSE) or the variable names of x which should be used for inference / as treatment variables.
#' @param I3 logical vector with same length as the number of controls;
#' indicates if variables (TRUE) should be included in any case.
#' @param \dots parameters passed to the \code{\link{rlasso}} function.
#' @return The function returns an object of class \code{rlassoLM} with the following entries: \item{coefficients}{vector with estimated
#' values of the coefficients for each selected variable} \item{se}{standard error (vector)}
#' \item{t}{t-statistic} \item{pval}{p-value} \item{samplesize}{sample size of the data set} \item{I}{index of variables of the union of the lasso regressions}
#' @references A. Belloni, V. Chernozhukov, C. Hansen (2014). Inference on
#' treatment effects after selection among high-dimensional controls. The
#' Review of Economic Studies 81(2), 608--650.
#' @keywords Estimation Inference Treatment effect High-dimensional controls
#' @export
#' @rdname rlassoLM

rlassoLM <- function(x, ...)
  UseMethod("rlassoLM") # definition generic function

#' @rdname rlassoLM
#' @export

rlassoLM.default <- function(x, y, index=c(1:ncol(x)), I3=NULL, ...) {

  if(is.logical(index)){
    k <- p1 <- sum(index)
  } else {
    k <- p1 <- length(index)
  }
  n <- dim(x)[1]
  # preprocessing index
  # numerischer Vektor
  if (is.numeric(index)){
    index <- as.integer(index)
    stopifnot(all(index<=ncol(x)) && length(index)<=ncol(x))
  } else {
    # logical Vektor
    if (is.logical(index)){
      stopifnot(length(index)==ncol(x))
      index <- which(index==T)
    } else {
      # character Vektor
      if(is.character(index)){
        stopifnot(all(is.element(index,colnames(x))))
        index <- which(is.element(colnames(x),index))
      } else {
        stop("argument index has an invalid type")
      }
    }
  }

  # check validity of I3
  I3ind <- which(I3==T)
  if (length(intersect(index, I3ind)!=0)) stop("I3 and index must not overlap!")
  if (is.null(colnames(x))) colnames(x) <- paste("V", 1:dim(x)[2], sep="")
  coefficients <- as.vector(rep(NA_real_,k))
  se <-  rep(NA_real_,k)
  t <-  rep(NA_real_,k)
  pval <-  rep(NA_real_,k)
  lasso.regs <- vector("list", k)
  reside <- matrix(NA, nrow=n, ncol=p1)
  residv <- matrix(NA, nrow=n, ncol=p1)
  names(coefficients) <- names(se) <- names(t) <- names(pval) <- names(lasso.regs) <- colnames(reside) <- colnames(residv) <-colnames(x)[index]

  for (i in 1:k)
  {
    d <- x[,index[i], drop=FALSE]
    Xt <- x[,-index[i], drop=FALSE]

    lasso.regs[[i]] <- try(col <- rlassoLMone(Xt,y,d, I3=I3, ...))
    if(class(lasso.regs[[i]]) == "try-error") {
      next
    }
    else {
      coefficients[i] <- col$alpha
      se[i] <- col$se
      t[i] <- col$t
      pval[i] <- col$pval
      reside[,i] <- col$residuals$epsilon
      residv[,i] <- col$residuals$v
    }
  }
  residuals <- list(e=reside, v=residv)
  res <- list(coefficients=coefficients, se=se, t=t, pval=pval, lasso.regs=lasso.regs, index=I, call=match.call(), samplesize=n, residuals=residuals)
  class(res) <- "rlassoLM"
  return(res)
}

#' @rdname rlassoLM
#' @export

rlassoLM.formula <- function(formula, data, index, I3=NULL, ...) {
  # TBD
}


#' @rdname rlassoLM
#' @param d variable for which inference is conducted (treatment variable)
#' @export

rlassoLMone <- function(x, y, d, I3=NULL,  ...) {
  d <- as.matrix(d, ncol=1)
  y <- as.matrix(y, ncol=1)
  kx <- dim(x)[2]
  if (is.null(colnames(d))) colnames(d) <- "d1"
  if (is.null(colnames(x)) & !is.null(x)) colnames(x) <- paste("x", 1:kx, sep="")
  I1 <- rlasso(x, d, ...)$index
  I2 <- rlasso(x, y, ...)$index
  n <- dim(x)[1]

  if (is.logical(I3)) {
    I <- I1+I2+I3
    I <- as.logical(I)
  } else {
    I <- I1+I2
    I <- as.logical(I)
  }
  if (sum(I)==0) {I <- NULL}
  x <- cbind(d,x[,I, drop=FALSE])
  reg1 <- lm(y~x)
  alpha <- coef(reg1)[2]
  xi<-reg1$residuals*sqrt(n/(n-sum(I)-1))
  if(is.null(I)) {reg2 <- lm(d~1)}
  if(!is.null(I)) {reg2 <- lm(d~x[,-1, drop=FALSE])}
  v <- reg2$residuals
  var <- 1/n*1/mean(v^2)*mean(v^2*xi^2)*1/mean(v^2)
  se <- sqrt(var)
  tval <- alpha/sqrt(var)
  pval <- 2*pnorm(-abs(tval))
  if (is.null(I)) {
    no.selected <- 1
  } else {
    no.selected <- 0
  }
  res <- list(epsilon = xi, v=v)
  return(list(alpha=unname(alpha), se=drop(se), t=unname(tval), pval=unname(pval), no.selected=no.selected, coefficients=coef(reg1), residuals=res))
}


# HC.lasso.mult <- function(x, y, d, I3=NULL,  ...) {
#   y <- as.matrix(y, ncol=1)
#   d <- as.matrix(d)
#   if (is.null(colnames(d))) colnames(d) <- paste("d", 1:dim(d)[2], sep="")
#
#   kd <- dim(d)[2]
#   kx <- dim(x)[2]
#   n <- dim(x)[1]
#
#   I2 <- lasso(x, y, ...)$index
#   I1 <- vector("logical", length=kx)
#   for (i in 1:kd) {
#   Id <- lasso(x, d[,i], ...)$index
#   I1 <- I1 + Id
#   }
#
#   if (is.logical(I3)) {
#     I <- I1+I2+I3
#     I <- as.logical(I)
#   } else {
#     I <- I1+I2
#     I <- as.logical(I)
#   }
#
#   if (sum(I)==0) {I <- NULL}
#   x <- cbind(d,x[,I, drop=FALSE])
#   reg1 <- lm(y~x)
#   xi<-reg1$residuals*sqrt(n/(n-sum(I)-1))
#   table <- matrix(NA,ncol=4,nrow=kd)
#   rownames(table) <- colnames(d)
#   colnames(table) <- c("coeff.", "se.", "t-value", "p-value")
#   for (i in 1:kd)
#   {
#     if(is.null(I)) {reg2 <- lm(d[,i]~1)}
#     if(!is.null(I)) {reg2 <- lm(d[,i]~x[,-i, drop=FALSE])}
#     v <- reg2$residuals
#     var <- 1/n*1/mean(v^2)*mean(v^2*xi^2)*1/mean(v^2)
#     alpha <- coef(reg1)[1+i]
#     table[i,1] <- alpha
#     table[i,2] <- sqrt(var)
#     table[i,3] <- alpha/sqrt(var)
#     table[i,4] <- 2*pnorm(-abs(alpha/sqrt(var)))
#   }
#   return(table)
# }


################# Methods for rlassoLM

#' Methods for S3 object \code{rlassoLM}
#'
#' Objects of class \code{rlassoLM} are constructed by \code{rlassoLM.formula} or \code{rlassoLM.default}.
#' \code{print.rlassoLM} prints and displays some information about fitted \code{rlassoLM} objects.
#' \code{summary.rlassoLM} summarizes information of a fitted \code{rlassoLM} object.
#' \code{confint.rlassoLM} extracts the confidence intervals.
#' \code{plot.rlassoLM} plots the estimates.
#'
#' @param object An object of class \code{rlassoLM}
#' @param x An object of class \code{rlassoLM}
#' @param digits significant digits in printout
#' @param ... arguments passed to the print function and other methods.
#' @keywords methods rlassoLM
#' @rdname methods.rlassoLM
#' @aliases methods.rlassoLM print.rlassoLM summary.rlassoLM confint.rlassoLM plot.rlassoLM
#' @export

print.rlassoLM <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (length(coef(x))) {
      cat("Coefficients:\n")
      print.default(format(coef(x), digits = digits), print.gap = 2L,
                    quote = FALSE)
    }
  else cat("No coefficients\n")
  cat("\n")
  invisible(coef(x))
}


#' @rdname methods.rlassoLM
#' @export

summary.rlassoLM <- function(object, digits = max(3L, getOption("digits") - 3L), ...) {
  if (length(coef(object))) {
    k <- length(object$coefficients)
    table <- matrix(NA,ncol=4,nrow=k)
    rownames(table) <- names(object$coefficient)
    colnames(table) <- c("coeff.", "se.", "t-value", "p-value")
    table[,1] <- object$coefficient
    table[,2] <- object$se
    table[,3] <- object$t
    table[,4] <- object$pval
    print("Estimation of the effect of selected variables in a high-dimensional regression")
    printCoefmat(table, digits=digits, P.values=TRUE, has.Pvalue=TRUE)
    cat("\n")
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(table)
}

#' @rdname methods.rlassoLM
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level	the confidence level required
#' @param joint logical, if joint confidence intervals should be clalculated.
#' @export

confint.rlassoLM <- function(object, parm, level=0.95, joint=FALSE, ...) {
  B <- 500 # number of bootstrap repitions
  n <- object$samplesize
  k <- p1 <- length(object$coefficient)
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm))
    parm <- pnames
  else if (is.numeric(parm))
    parm <- pnames[parm]
  if (!joint) {
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qt(a, n-k)
  pct <- format.perc(a, 3)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm,
                                                             pct))
  ses <- object$se[parm]
  ci[] <- cf[parm] + ses %o% fac
  }

  if (joint) {
    phi <- object$residuals$e*object$residuals$v
    m <- 1/sqrt(colMeans(phi^2))
    phi <- t(t(phi)/m)
    sigma <- sqrt(colMeans(phi^2))
    sim <- vector("numeric", length=B)
    for (i in 1:B) {
      xi <- rnorm(n)
      phi_temp <- phi*xi
      Nstar <- 1/sqrt(n)*colSums(phi_temp)
      sim[i] <- max(abs(Nstar))
    }
    a <- (1 - level)/2
    ab <- c(a, 1 - a)
    pct <- format.perc(ab, 3)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm,
                                                               pct))
    hatc <- quantile(sim, probs=1-a)
    ci[,1] <- cf[parm] - hatc*1/sqrt(n)*sigma
    ci[,2] <- cf[parm] + hatc*1/sqrt(n)*sigma
  }
  return(ci)
}

#' @rdname methods.rlassoLM
#' @export

plot.rlassoLM <- function(x, ..., var.names=T, col="black", xlim=NULL, xlab="", ylab="", main="", sub="",
                          pch=19, cex=1.3, lwd=2, lty=1){
#plot.rlassoLM <- function(x,  var.names=T, ...){

    coefmatrix <- cbind(summary(x), confint(x))[,c(1,5,6)]
    rownames(coefmatrix) <- names(x$coefficients)
    #colnames(coefmatrixneu) <- c("estimate", "lower CI", "upper CI")
   coefmatrix <- coefmatrix[order(abs(coefmatrix[,1])),]

  p <- nrow(coefmatrix)

  # arranging of plots

  if(p<=20){
    par(mfrow=c(1,1))
    ngr <- p
  }else{
    if(p<=40){
      par(mfrow=c(1,2))
      ngr <- ceiling(p/2)
    }else{
      par(mfrow=c(2,2))
      ngr <- 10
    }
  }

  # Anzahl benötigter Fenster
  nWin <- ceiling(p/ngr)

  # noch nicht geplottete Parameter
  puebrig <- p

  # Festlegung des Wertebereichs
  if(is.null(xlim)){
    low <- min(coefmatrix)
    high <- max(coefmatrix)
  } else{
    low <- xlim[1]
    high <- xlim[2]
  }

  # Erzeugen der Grafikfenster
  for(i in 1:nWin){

    plot(1,xlim=c(low,high),ylim=c(0,ngr),yaxt="n",xlab=xlab,ylab=ylab,type="n",main=main,sub=sub)
    abline(v=0)

    # Achsenbeschriftung der y-Achse mit den Variablennamen? default: ja
    if(var.names==T){
      axis(side=2,at=c(ngr:(ngr-min(puebrig,ngr)+1)),labels=rownames(coefmatrix)[1:min(puebrig,ngr)],las=2)
    }

    # Punkte und Linien entsprechend der Koeffizientenwerte
    for(j in 1:min(puebrig,ngr)){
      points(coefmatrix[j,1],ngr+1-j,pch=pch,cex=cex,col=col)
    }
    for(j in 1:min(puebrig,ngr)){
      lines(x=c(coefmatrix[j,2],coefmatrix[j,3]),y=c(ngr+1-j,ngr+1-j),lwd=lwd,col=col,lty=lty)
      lines(x=c(coefmatrix[j,2],coefmatrix[j,2]),y=c(ngr+1-j+0.2,ngr+1-j-0.2),lwd=lwd,col=col,lty=lty)
      lines(x=c(coefmatrix[j,3],coefmatrix[j,3]),y=c(ngr+1-j+0.2,ngr+1-j-0.2),lwd=lwd,col=col,lty=lty)
    }

    # Matrix und übrige Parameter anpassen
    puebrig <- puebrig-min(puebrig,ngr)
    coefmatrix <- coefmatrix[-(1:ngr),,drop=F]
  }
  par(mfrow=c(1,1))
  return(NULL)
}
