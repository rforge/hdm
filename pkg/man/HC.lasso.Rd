\name{HC.lasso}
\alias{HC.lasso}
\alias{HC.lasso.mult}
\title{
Estimation and inference on treatment effects after selection amongst high-dimensional controls
}
\description{
The functions estimates the effect of a treatment effect alpha0 in a setting of high-dimensional controls. The so-called post-double-selection estimator is implemented. The idea is to select variable by regression of the outcome variable on the control variables and the treatment variable on the control variables. The final estimation is done by a regression of the outcome on the treatment effect and the union of the selected variables in the first two steps. The resulting estimator for alpha0 is normal distributed which allows inference on th treatment effect.
}
\usage{
HC.lasso(x, y, d, I3 = NULL, ...)
}

\arguments{
  \item{x}{control variables}
  \item{y}{outcome / dependent variable}
  \item{d}{treatment variable}
  \item{I3}{logical vector with same length as the number of controls; indicates if variables (TRUE) should be included in any case.}
  \item{\dots}{parameters passed to the \code{\link{lasso}} function.}
}
\details{
TBD
}
\value{
The function returns a list with the following components
alpha, se=se, t=tval, pval=pval, no.selected=no.selected))
\item{alpha}{estimated value of the treatment effect}
\item{se}{standard error}
\item{t}{t-statistic}
\item{pval}{p-value}
\item{no.selected}{0/1 variable indicating if no variables were selected (\code{no.selected=1}) or control variables were selected (\code{no.selected=0}). }
\item{cofficients}{parameter vector of the estimated coefficients, including the treatment effect and the control variables.}
}
\references{
A. Belloni, V. Chernozhukov, C. Hansen (2014). Inference on treatment effects after selection among high-dimensional controls. The Review of Economic Studies 81(2), 608--650.
}
%\author{
%Martin Spindler
%}

\keyword{Estimation}
\keyword{Inference}
\keyword{Treatment effect}
\keyword{High-dimensional controls}
