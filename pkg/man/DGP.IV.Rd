\name{DGP.IV}
\alias{DGP.IV}
\title{
Data Generating Process for IV Estimation in a High-dimensional Setting
}
\description{
The function implements the data generating process used in Belloni et al. (2012)
}
\usage{
DGP.IV(n = 250, p = 100, pnz = 50, Fstat = 180, control = list(s2e = 1, Cev = 0.6, s2z = 1, szz = 0.5, pi1 = 1, alpha = 1))
}
\arguments{
  \item{n}{sample size}
  \item{p}{number of instruments}
  \item{pnz}{number of relevant instruments}
  \item{Fstat}{concentration paramteter}
  \item{control}{list with control parameters}
}
\details{
The function implements the DGP in Belloni et al. (2012), p. 2397-2398. For the cut-off design set \code{pi1=1} and \code{pnz<100} and for the exponential design set \code{pnz=100} and \code{pi<1}.
}
\value{
\item{d}{endogenous variable}
\item{y}{dependent variable}
\item{Z}{instrumental variables}
\item{setup}{\eqn{\alpha} treatment effect (coefficient) and \code{Fstat} the concentration parameter}
}
\references{
A. Belloni, D. Chen, V. Chernozhukov and C. Hansen (2012). Sparse models and methods for optimal instruments with an application to eminent domain. \emph{Econometrica} 80 (6), 2369--2429.
}
%\author{
%Martin Spindler
%}

\keyword{DGP}
\keyword{Instrumental Variables}
