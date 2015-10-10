# TBD How does selection on X work??

#' Instrumental Variable Estimation with Lasso
#'
#' This function selects the instrumental variables in the first stage by
#' Lasso. First stage predictions are then used in the second stage as optimal
#' instruments to estimate the parameter vector. The function returns an element of class \code{rlassoIVselectX}
#'
#' The implementation follows the procedure described in Belloni et al. (2012).
#' option \code{post=TRUE} conducts post-lasso estimation, i.e. a refit of the
#' model with the selected variables, to estimate the optimal instruments. The
#' parameter vector of the structural equation is then fitted by two-stage
#' least square (tsls) estimation. If variables of the exogenous variables in
#' \code{x} should be used as instruments, they have to be added to the
#' instrument set \code{z} explicitly.
#'
#' @param x exogenous variables in the structural equation (matrix)
#' @param d endogenous variables in the structural equation (vector or matrix)
#' @param y outcome or dependent variable in the structural equation (vector or matrix)
#' @param z set of potential instruments for the endogenous variables.
#' Exogenous variables serve as their own instruments.
#' @param post logical. If \code{TRUE}, post-lasso estimation is conducted.
#' @param \dots arguments passed to the function \code{rlasso}
#' @return An object of class \code{rlassoIVselectX} containing at least the following
#' components: \item{coefficients}{estimated parameter vector}
#' \item{vcov}{variance-covariance matrix} \item{residuals}{
#' residuals} \item{samplesize}{sample size}
#' @references D. Belloni, D. Chen, V. Chernozhukov and C. Hansen (2012).
#' Sparse models and methods for optimal instruments with an application to
#' eminent domain. \emph{Econometrica} 80 (6), 2369--2429.
#' @keywords Instrumental Variables Lasso Hig-dimensional setting
#' @export
#' @rdname rlassoIVselectX
rlassoIVselectX <- function(x,d,y,z, post=TRUE, ...) {

  d <- as.matrix(d)
  if(is.vector(x)) x <-  as.matrix(x)
  n <- length(y)
  kex <- dim(x)[2]
  ke <- dim(d)[2]
  kiv <- dim(z)[2]

  if (is.null(colnames(d))) colnames(d) <- paste("d", 1:ke, sep="")
  if (is.null(colnames(x)) & !is.null(x)) colnames(x) <- paste("x", 1:kex, sep="")
  # first stage regression
  Dhat <- NULL
  for (i in 1:ke) {
    di <- d[,i]
    lasso.fit <- rlasso(z, di, post=post, ...)
    if (sum(lasso.fit$ind)==0) {
      dihat <- rep(mean(di),n) #dihat <- mean(di)
    } else {
      dihat <- z%*%lasso.fit$coefficients
    }
    Dhat <- cbind(Dhat, dihat)
  }
  Dhat <- cbind(Dhat, x)
  d <- cbind(d,x)
  # calculation coefficients
  alpha.hat <- MASS::ginv(t(Dhat)%*%d)%*%(t(Dhat)%*%y)
  # calcualtion of the variance-covariance matrix
  residuals <- y - d%*%alpha.hat
  Omega.hat <- t(Dhat)%*%diag(as.vector(residuals^2))%*%Dhat #  Dhat.e <- Dhat*as.vector(residuals);  Omega.hat <- t(Dhat.e)%*%Dhat.e
  Q.hat.inv <- MASS::ginv(t(d)%*%Dhat) #solve(t(d)%*%Dhat)
  vcov <- Q.hat.inv%*%Omega.hat%*%t(Q.hat.inv)
  rownames(alpha.hat) <- c(colnames(d))
  colnames(vcov) <- rownames(vcov) <- rownames(alpha.hat)
  res <- list(coefficients=alpha.hat[,1], vcov=vcov, residuals=residuals, samplesize=n, call=match.call())
  class(res) <- "rlassoIVselectX"
  return(res)
}


################# Methods for rlassoIVselectX

#' Methods for S3 object \code{rlassoIVselectX}
#'
#' Objects of class \code{rlassoIVselectX} are constructed by \code{rlassoIVselectX.formula} or \code{rlassoIVselectX.default}.
#' \code{print.rlassoIVselectX} prints and displays some information about fitted \code{rlassoIVselectX} objects.
#' \code{summary.rlassoIVselectX} summarizes information of a fitted \code{rlassoIVselectX} object.
#' \code{confint.rlassoIVselectX} extracts the confidence intervals.
#' @param object an object of class \code{rlassoIVselectX}
#' @param x an object of class \code{rlassoIVselectX}
#' @param digits significant digits in printout
#' @param ... arguments passed to the print function and other methods
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level	the confidence level required.
#' @keywords methods rlassoIVselectX
#' @rdname methods.rlassoIVselectX
#' @aliases methods.rlassoIVselectX print.rlassoIVselectX summary.rlassoIVselectX
#' @export

print.rlassoIVselectX <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
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

#' @rdname methods.rlassoIVselectX
#' @export

summary.rlassoIVselectX <- function(object, digits = max(3L, getOption("digits") - 3L), ...) {
  if (length(coef(object))) {
    k <- length(object$coefficient)
    table <- matrix(NA,ncol=4,nrow=k)
    rownames(table) <- names(object$coefficient)
    colnames(table) <- c("coeff.", "se.", "t-value", "p-value")
    table[,1] <- object$coefficient
    table[,2] <- sqrt(diag(as.matrix(object$vcov)))
    table[,3] <- table[,1]/table[,2]
    table[,4] <- 2*pnorm(-abs(table[,3]))
    print("Estimation of the effect of selected variables in a high-dimensional regression")
    printCoefmat(table, digits=digits,  P.values=TRUE, has.Pvalue=TRUE)
    cat("\n")
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(table)
}

#' @rdname methods.rlassoIVselectX
#' @export

confint.rlassoIVselectX <- function(object, parm, level=0.95, ...) {
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
  ses <- sqrt(diag(object$vcov))[parm]
  ci[] <- cf[parm] + ses %o% fac
  print(ci)
  invisible(ci)
}

