\name{lasso.significance}
\alias{lasso.significance}

\title{
Joint Confidence Intervals for Lasso
}

\description{
The function calculates joint confidence intervals for coefficients in a Lasso regression
}
\usage{
lasso.significance(x, y, index = c(1:ncol(x)), alpha = 0.05, method = "bootstrap", B = 500, ...)
}
\arguments{
  \item{x}{regressors}
  \item{y}{dependent variable}
  \item{index}{index of the the variables for which confidence intervals should be calculated.}
  \item{alpha}{significane level}
  \item{method}{bootstrap}
  \item{B}{number of bootstrap repititions}
  \item{\dots}{arguments passed to the function \code{HC.lasso}}
}
\details{
The function implements the method proposed in Belloni et al. (2014) for linear Lasso models.
}
\value{
The function returns a list consisting of the following components
\item{alpha}{estimated coefficient}
\item{se}{standard error}
\item{lower CI}{lower bound of the Confidence Interval}
\item{upper CI}{upper bound of the Confidence Interval}
}
\references{
A. Belloni, C. Hansen and V. Chernozhukov (2014). Unifrom Post Selection Inference for Least Absolute Deviation Regression and Other Z-estimation Problems. Biometrika, 1--17.
}
%\author{
%Martin Spindler
%}

\keyword{Lasso}
\keyword{Significance}
\keyword{Confidence Intervals}
