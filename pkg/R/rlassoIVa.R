#' Post-Selection and Post-Regularization Inference in Linear Models with Many
#' Controls and Instruments
#'
#' The function estimates a treatment effect in a setting with very many
#' controls and very many instruments (even larger than the sample size).
#'
#' The implementation follows the procedure described in Chernozhukov et al.
#' (2015) and is built on "triple selection" to achieve an orthogonal moment
#' function. The function returns an object of S3 class \code{rlassoManyIV}.
#' Moreover, it is wrap function for the case that selection should be done only with the instruments Z (\code{rlassoIVselectZ}) or with 
#' the control variables X (\code{rlassoIVselectX}).
#'
#' @aliases rlassoManyIV rlassolassoManyIVmult
#' @param x matrix of exogenous variables
#' @param z matrix of instrumental variables
#' @param y outcome / dependent variable (vector or matrix)
#' @param d endogenous variable
#' @param \dots arguments passed to the function \code{rlasso}
#' @return An object of class \code{rlassoManyIV} containing at least the following
#' components: \item{coefficients}{estimated parameter value}
#' \item{se}{variance-covariance matrix}
#' @references V. Chernozhukov, C. Hansen, M. Spindler (2015). Post-selection
#' and post-regularization inference in linear models with many controls and
#' instruments. American Economic Review: Paper & Proceedings 105(5), 1--7.
#' @keywords Lasso Many controls and instruments Instrumental Variable
#' @rdname rlassoManyIV
#' @export

rlassoManyIV <- function(x,z,y,d,...) {
  d <- as.matrix(d)
  if (is.null(colnames(d))) colnames(d) <- paste("d", 1:ncol(d), sep="")
  if (is.null(colnames(x)) & !is.null(x)) colnames(x) <- paste("x", 1:ncol(x), sep="")
  if (is.null(colnames(z)) & !is.null(z)) colnames(z) <- paste("z", 1:ncol(z), sep="")
  n <- length(y)
  Z <- cbind(z,x)
  lasso.d.zx <- rlasso(Z,d,...)
  lasso.y.x <- rlasso(x,y,...)
  lasso.d.x <- rlasso(x,d,...)
  if (sum(lasso.d.zx$index)==0) return(list(alpha=NA, se=NA))
  ind.dzx <- lasso.d.zx$index
  PZ <- Z[,ind.dzx]%*%MASS::ginv(t(Z[,ind.dzx])%*%Z[,ind.dzx])%*%t(Z[,ind.dzx])%*%d
  lasso.PZ.x <- rlasso(x,PZ,...)
  ind.PZx <- lasso.PZ.x$index
  Dr <- d- x[,ind.PZx]%*%MASS::ginv(t(x[,ind.PZx])%*%x[,ind.PZx])%*%t(x[,ind.PZx])%*%PZ
  Yr <- lasso.y.x$residuals
  Zr <- lasso.PZ.x$residuals
  result <- tsls(Yr,Dr,x=NULL,Zr)
  coef <- as.vector(result$coefficient)
  se <- diag(sqrt(result$vcov))
  names(coef) <- names(se) <- colnames(d)
  res <- list(coefficients=coef, se=se, vcov=vcov, call=match.call(), samplesize=n)
  class(res) <- "rlassoManyIV"
  return(res)
}

#' @rdname rlassoManyIV
#' @export

rlassoManyIVmult <- function(x,z,y,d,...) {
  #browser()
  d <- as.matrix(d)
  if (is.null(colnames(d))) colnames(d) <- paste("d", 1:ncol(d), sep="")
  if (is.null(colnames(x)) & !is.null(x)) colnames(x) <- paste("x", 1:ncol(x), sep="")
  if (is.null(colnames(z)) & !is.null(z)) colnames(z) <- paste("z", 1:ncol(z), sep="")
  d <- as.matrix(d)
  n <- dim(x)[1]
  d <- as.matrix(d)
  kd <- dim(d)[2]
  Z <- cbind(z,x)
  if (is.null(colnames(d))) colnames(d) <- paste("d", 1:kd, sep="")

  lasso.y.x <- rlasso(x,y,...)
  Yr <- lasso.y.x$residuals
  Drhat <- NULL
  Zrhat <- NULL
  for (i in 1:kd) {
  lasso.d.x <- rlasso(x,d[,i],...)
  lasso.d.zx <- rlasso(Z,d[,i],...)
  if (sum(lasso.d.zx$index)==0) {
    Drhat <- cbind(Drhat, d[,i] - mean(d[,i]))
    Zrhat <- cbind(Zrhat, d[,i] - mean(d[,i]))
    next
  }
  ind.dzx <- lasso.d.zx$index
  PZ <- Z[,ind.dzx,drop=FALSE]%*%solve(t(Z[,ind.dzx,drop=FALSE])%*%Z[,ind.dzx,drop=FALSE])%*%t(Z[,ind.dzx,drop=FALSE])%*%d[,i,drop=FALSE]
  lasso.PZ.x <- rlasso(x,PZ,...)
  ind.PZx <- lasso.PZ.x$index
  Dr <- d[,i]- x[,ind.PZx,drop=FALSE]%*%solve(t(x[,ind.PZx,drop=FALSE])%*%x[,ind.PZx,drop=FALSE])%*%t(x[,ind.PZx,drop=FALSE])%*%PZ
  Zr <- lasso.PZ.x$residuals
  Drhat <- cbind(Drhat, Dr)
  Zrhat <- cbind(Zrhat, Zr)
  }
  result <- tsls(Yr,Drhat,x=NULL,Zrhat)
  coef <- as.vector(result$coefficient)
  se <- sqrt(diag(result$vcov))
  names(coef) <- names(se) <- colnames(d)
  res <- list(coefficients = coef, se=se, vcov=result$vcov, call=match.call(), samplesize=n)
  class(res) <- "rlassoManyIV"
  return(res)
}


################# Methods for rlassoIV

#' Methods for S3 object \code{rlassoManyIV}
#'
#' Objects of class \code{rlassoManyIV} are constructed by \code{rlassoManyIV}.
#' \code{print.rlassoManyIV} prints and displays some information about fitted \code{rlassoManyIV} objects.
#' \code{summary.rlassoManyIV} summarizes information of a fitted \code{rlassoManyIV} object.
#' \code{confint.rlassoManyIV} extracts the confidence intervals.
#' @param object An object of class \code{rlassoManyIV}
#' @param x An object of class \code{rlassoManyIV}
#' @param digits significant digits in printout
#' @param ... arguments passed to the print function and other methods
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level	the confidence level required.
#' @keywords methods rlassoManyIV
#' @rdname methods.rlassoManyIV
#' @aliases methods.rlassoManyIV print.rlassoManyIV summary.rlassoManyIV
#' @export

print.rlassoManyIV <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
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

#' @rdname methods.rlassoManyIV
#' @export

summary.rlassoManyIV <- function(object, digits = max(3L, getOption("digits") - 3L), ...) {
  if (length(coef(object))) {
    k <- length(object$coefficient)
    table <- matrix(NA,ncol=4,nrow=k)
    rownames(table) <- names(object$coefficients)
    colnames(table) <- c("coeff.", "se.", "t-value", "p-value")
    table[,1] <- object$coefficients
    table[,2] <- object$se
    table[,3] <- table[,1]/table[,2]
    table[,4] <- 2*pnorm(-abs(table[,3]))
    cat("Estimation of the effect of selected variables in a high-dimensional regression", "\n")
    printCoefmat(table, digits=digits,  P.values=TRUE, has.Pvalue=TRUE)
    cat("\n")
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(table)
}

#' @rdname methods.rlassoManyIV
#' @export

confint.rlassoManyIV <- function(object, parm, level=0.95, ...) {
  n <- object$samplesize
  k <- length(object$coefficients)
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm))
    parm <- pnames
  else if (is.numeric(parm))
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qt(a, n-k)
  pct <- format.perc(a, 3)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm,
                                                             pct))
  ses <- object$se[parm]
  ci[] <- cf[parm] + ses %o% fac
  print(ci)
  invisible(ci)
}




