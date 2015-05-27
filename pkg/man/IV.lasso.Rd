\name{IV.lasso}
\alias{IV.lasso}
\title{
Instrumental Variable Estimation with Lasso
}
\description{
This function selects the instrumental variables in the first stage by Lasso. First stage predictions are then used in the second stage as optimal instruments to estimate the parameter vector.
}
\usage{
IV.lasso(x, d, y, z, post=TRUE, ...)
}
\arguments{
  \item{x}{exogenous variables in the structural equation}
  \item{d}{endogenous variables in the structural equation}
  \item{y}{outcome or depenent variable in the structural equation}
  \item{z}{set of potential instruments for the endogenous variables. Exogenous variables serve as their own instruments.}
  \item{post}{logical. If \code{TRUE}, post-lasso estimation is conducted.}
  \item{\dots}{arguments passed to the function \code{lasso}
}
}
\details{The implementation follows the procedure described in Belloni et al. (2012). In the first stage the selected with lasso with data-driven penalties. The option \code{post=TRUE} conducts post-lasso estimation, i.e. a refit of the model with the selected variables, to estimate the optimal instruments. The parameter vector of the structural equation is then fitted by two-stage least square (tsls) estimation.}
\value{
An object of class list containing at least the following components:
\item{coefficients}{estimated parameter vector}
\item{vcov}{variance-covariance matrix}
\item{residuals}{estimated residuals}
\item{samplesize}{sample size}
}
\references{
D. Belloni, D. Chen, V. Chernozhukov and C. Hansen (2012). Sparse models and methods for optimal instruments with an application to eminent domain. \emph{Econometrica} 80 (6), 2369--2429.
}
%\author{
%Martin Spindler
%}
\keyword{Instrumental Variables}
\keyword{Lasso}
\keyword{Hig-dimensional setting}