\name{HC.lasso.wrap}
\alias{HC.lasso.wrap}
\title{
Wrap function for \code{\link{HC.lasso}}
}
\description{
The functions treats each of the variables contained in \code{x} or some subset specified by \code{index} as treatment variable (and all other variables as controls) and applies the function \code{\link{HC.lasso}}.
}
\usage{
HC.lasso.wrap(x, y, index = c(1:ncol(x)), I3 = NULL, ...)
}
\arguments{
  \item{x}{matrix of regressor variables serving as controls and potential treatments
}
  \item{y}{outcome variable
}
  \item{index}{vector of integers indicating the position (column) of variables of x which should be used as treatment variables.
}
  \item{I3}{logical vector with same length as the number of controls; indicates if variables (TRUE) should be included in any case.}
  \item{\dots}{parameters passed to the \code{\link{lasso}} function.}
}
\details{TBD
}
\value{
The function returns a table where for each treatment variable (rows) the following quantitites (columns) are given
\item{alpha}{estimated value of the treatment effect}
\item{se}{standard error}
\item{t}{t-statistic}
\item{pval}{p-value}
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
