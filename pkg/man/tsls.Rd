\name{tsls}
\alias{tsls}
\title{Two-Stage Least Squares Estimation (TSLS)
}
\description{The function computes tsls estimate (coefficients) and variance covariance matrix assuming homoskedasticity for outcome variable y where d are endogenous variables, in structural equation, x are exogensous variables in structural equation and z are instruments.  x should include the constant term.
}
\usage{
tsls(y, d, x, z)
}
\arguments{
  \item{y}{outcome variable}
  \item{d}{endogenous variables}
  \item{x}{exogenous variables}
  \item{z}{instruments}
}
\details{The function does Two-Stage Least Squares Estimation (TSLS).}
\value{
\item{coefficients}{coefficients}
\item{vcov}{variance-covariance matrix}
\item{residuals}{difference of outcome minus predicted values}
}
%\author{
%Martin Spindler
%}
\keyword{Instrumental Variables}
\keyword{Endogeneity}
\keyword{2SLS}
\keyword{TSLS}
