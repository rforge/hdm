\name{ProgEval}
\alias{ProgEval}
\title{Function for Program Evaluation in a high-dimensional setting}
\description{
The function calculates local average treatment effects (LATE) and local quantile treatment effects (LQTE), and also the versions for the treated (LATT and LQTT).
}
\usage{
ProgEval(x, y, z, d, post = TRUE, intercept = TRUE, normalize = TRUE, tau=(5:95)/100, alpha=0.05, bootstrap = NULL, nRep = 500)
}
\arguments{
  \item{x}{exogenous / control variables (matrix)}
  \item{y}{outcome variable}
  \item{z}{instrumental variables}
  \item{d}{treatment variable (binary, 0/1)}
  \item{post}{logical. If \code{TRUE}, post-lasso estimation is conducted.}
  \item{intercept}{logical. If \code{TRUE}, intercept is included which is not penalized.}
  \item{normalize}{logical. If \code{TRUE}, design matrix \code{x} is scaled.}
  \item{tau}{quantiles for which the LQTE and LQTT are calculated}
  \item{bootstrap}{bootstrap method: "none", "Bayes", "normal", "wild"}
  \item{nRep}{number of repititions for the bootstrap}
  \item{alpha}{confidence level for the confidence bands of the LQTE and LQTT}
}
\details{Details of the implementation and technical details can be found in Belloni et al. (2014).  }
\value{
The function returns an object of type \code{list} with the following components
  \item{LATE}{list "LATE"" with value "late", standard errors "s.e.", bootstrap standard errors "boot.se" and the individual               effects}
  \item{LATT}{list "LATT"" with value "latt", standard errors "s.e.", bootstrap standard errors "boot.se" and the individual effects}
    \item{LQTE}{list "LQTE"" with values "lqte" for each tau, bootstrap standard errors "se" and t-statistic}
  \item{LQTT}{list "LQTT"" with values "lqtt" for each tau, bootstrap standard errors "se" and t-statistic}
  \item{bootstrap}{chosen bootstrap method}
  \item{tau}{quantiles for wich the quantile effects should be calculated}
}
\references{
A. Belloni, V. Chernozhukov, I. Fernandez-Val, and C. Hansen (2014). Program evaluation with high-dimensional data. Working Paper.
}
%\author{
%%  Martin Spindler
%}

\keyword{Program Evaluation}
\keyword{Treatment Effects}
