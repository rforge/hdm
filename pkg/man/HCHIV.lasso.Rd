\name{HCHIV.lasso}
\alias{HCHIV.lasso}
\alias{HCHIV.lasso.mult}
\title{
Post-Selection and Post-Regularization Inference in Linear Models with Many Controls and Instruments
}
\description{
The function estimates a treatment effect in a setting with very many controls and very many instruments (even larger than the sample size).
}
\usage{
HCHIV.lasso(x, z, y, d, ...)
}
\arguments{
  \item{x}{matrix of exogenous variables}
  \item{z}{matrix of instrumental variables}
  \item{y}{outcome / dependent variable (vector or matrix)}
  \item{d}{endogenous variable}
  \item{\dots}{arguments passed to the function \code{lasso}}
}
\details{
The implementation follows the procedure described in Chernozhukov et al. (2015) and is built on "triple selection" to achieve an orthogonal moment function.
}
\value{
An object of class list containing at least the following components:
\item{coefficient}{estimated parameter value}
\item{se}{variance-covariance matrix}
}
\references{
V. Chernozhukov, C. Hansen, M. Spindler (2015). Post-selection and post-regularization inference in linear models with many controls and instruments. American Economic Review: Paper & Proceedings 105(5), 1--7.
}
%\author{
%Martin Spindler
%}
\keyword{Lasso}
\keyword{Many controls and instruments}
\keyword{Instrumental Variable}
