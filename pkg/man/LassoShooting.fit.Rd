\name{LassoShooting.fit}
\alias{LassoShooting.fit}
\title{
Lasso Shooting
}
\description{
Implement of the Shooting Lasso (Fu, 1998) with variable dependent penalization weights.
}
\usage{
LassoShooting.fit(x, y, lambda, control = list(maxIter = 1000, optTol = 10^(-5), zeroThreshold = 10^(-6)), XX = NULL, Xy = NULL, beta.start = NULL)
}
\arguments{
  \item{x}{matrix of regressor variables (\code{n} times \code{p} where \code{n} denotes the number of observations and \code{p} the number of regressors)}
  \item{y}{dependent variable (vector or matrix)}
  \item{lambda}{vector of length \code{p} of penalization parameters for each regressor}
  \item{control}{list wiht control parameters: \code{maxIter} maximal number of iterations, \code{optTol} tolerance for parameter precision, \code{zeroThreshold} threshold applied to the estimated cofficients to clean for numerical issues.}
  \item{XX}{optional, precalculated matrix \eqn{t(X)\%*\%X}}
  \item{Xy}{optional, precalculated matrix \eqn{t(X)\%*\%y}}
  \item{beta.start}{start value for beta}
}
\details{
The function impements the Shooting Lasso (Fu, 1998) with variable dependent penalization. The argruments \code{XX} and \code{Xy} are optional. Precalculated matrices might improve performance.}
\value{
\item{coefficients}{estimated coefficients by the Shooting Lasso Algorithm}
\item{coef.list}{list of coefficient from each iteration}
\item{num.it}{number of iterations run}
}
\references{
Fu, W. (1998). Penalized regressions: the bridge vs the lasso. \emph{Journal of Computational and Graphical Software} 7, 397-416.
}
\author{Martin Spindler}

\keyword{Lasso}
\keyword{Shooting Lasso}
