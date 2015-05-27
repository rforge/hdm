\name{late}
\alias{late}
\alias{latt}
\alias{lasf}
\alias{lasft}
\alias{localQuantile}
\title{
Functions for estimation of local average structural / treatment effects
}
\description{The functions estimate local average treatment effects (late), the late of the tretated (latt), local average structural effects (lasf), and the lasf of the treated (lasft). The functions rely on immunized / orthogonal moment conditions. Further details can be found in Belloni et al. (2014).
}
\usage{
late(y, d, z, my_z1x, my_z0x, mz_x, md_z1x, md_z0x, bootstrap = NULL, nRep = 500)
late(y,d,z,my_z1x,my_z0x,mz_x,md_z1x,md_z0x, bootstrap=NULL, nRep=500)
lasf(y,d,z,myd1_z1x,myd1_z0x,myd0_z1x,myd0_z0x,mz_x,md_z1x,md_z0x,d0, weights=NULL, bootstrap=NULL, nRep=500)
lasft(y,d,z,myd1_z0x,myd0_z0x,mz_x,md_z0x,d0, weights=NULL, bootstrap=NULL, nRep=500)
localQuantile(yF,d,z,taus,percs,myd1_z1x,myd0_z1x,myd1_z0x,myd0_z0x,mz_x,md_z1x,md_z0x, bootstrap=NULL, nRep=500, alpha=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{outcome variable / dependent variable}
  \item{yF}{matrix with columns of form 1(y < u) for u approximately covering the support of y}
  \item{d}{treatment variable}
  \item{z}{instrumental variables}
  \item{my_z1x}{E[Y | Z=1, X]}
  \item{my_z0x}{E[Y | Z=0, X]}
  \item{mz_x}{E[Z | X]}
  \item{myd1_z1x}{E[Y*1(D=1) | Z=1, X]}
  \item{myd1_z0x}{E[Y*1(D=1) | Z=0, X]}
  \item{myd0_z1x}{E[Y*1(D=0) | Z=1, X]}
  \item{myd0_z0x}{E[Y*1(D=0) | Z=0, X]}
  \item{md_z1x}{E[D | Z=1, X]}
  \item{md_z0x}{E[D | Z=0, X]}
  \item{d0}{d0=0 or d0=1}
  \item{weights}{weights}
  \item{bootstrap}{boostrapp method which should be employed: "none", "Bayes", "normal", "wild"}
  \item{nRep}{number of replications for the bootstrap}
  \item{alpha}{(1-alpha) confidence intervals will be calculated}
  \item{taus}{quantiles to be calculated}
  \item{percs}{unique quantiles of the data}
}
\details{
Details can be found in Belloni et al. (2014).
}
\value{
Functions return the estimated effects with standard errors and individual effects in the form of a \code{list}.
}
\references{
A. Belloni, V. Chernozhukov, I. Fernandez-Val, and C. Hansen (2014). Program evaluation with high-dimensional data. Working Paper.
}
%\author{
%Martin Spindler
%}

\keyword{local average treatment effect}
\keyword{late}
\keyword{latt}
\keyword{local average structural effect}
\keyword{lasf}
\keyword{lasft}
