#' Two-Stage Least Squares Estimation (TSLS)
#'
#' The function does Two-Stage Least Squares Estimation (TSLS).
#'
#' The function computes tsls estimate (coefficients) and variance-covariance-matrix assuming homoskedasticity
#' for outcome variable \code{y} where \code{d} are endogenous variables in structural equation, \code{x} are exogensous variables in
#' structural equation and z are instruments.  \code{x} should include a constant term.
#'
#' @param y outcome variable
#' @param d endogenous variables
#' @param x exogenous variables
#' @param z instruments
#' @param intercept logical, if intercept is included
#' @return The function returns a list with the following elements \item{coefficients}{coefficients}
#' \item{vcov}{variance-covariance matrix} \item{residuals}{outcome minus predicted values}
#' @keywords Instrumental Variables
#' @keywords Endogeneity
#' @keywords 2SLS
#' @keywords TSLS
#' @export


tsls <- function(y, d, x, z, intercept=TRUE) {
  n <- length(y)
  
  if (intercept==TRUE && is.null(x)) x <- rep(1,n)
  if (intercept==TRUE && !is.null(x)) x <- cbind(1,x)
  a1 <- dim(d)[2]
  a2 <- dim(x)[2]
  if (is.null(x)) {
    a2 <- 0
  }
  if (is.vector(x)) {
    a2 <- 1
  }
  if (is.vector(d)) {
    a1 <- 1
  }
  k <- a1 + a2
  X <- cbind(d, x)
  Z <- cbind(z, x)

  Mxz <- t(X) %*% Z
  Mzz <- MASS::ginv(t(Z) %*% Z)
  M <- MASS::ginv(Mxz %*% Mzz %*% t(Mxz))

  b <- M %*% Mxz %*% Mzz %*% (t(Z) %*% y)
  Dhat <- Z %*% MASS::ginv(t(Z) %*% Z) %*% t(Z) %*% X
  b2 <- MASS::ginv(t(Dhat) %*% X) %*% (t(Dhat) %*% y)
  e <- y - X %*% b
  VC1 <- as.numeric((t(e) %*% e/(n - k))) * M
  return(list(coefficients = b, vcov = VC1, se=sqrt(diag(VC1)), residuals = e, call=match.call(), samplesize=n))
}
