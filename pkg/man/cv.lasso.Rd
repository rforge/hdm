\name{cv.lasso}
\alias{cv.lasso}
\title{
Cross validation for penalty lambda for Lasso
}
\description{
\code{cv.lasso} computes the K-fold cross-validated mean squared prediction error for lasso.
}
\usage{
cv.lasso(x, y, K = 10, lambda.grid = NULL, post = TRUE, intercept = TRUE, normalize = TRUE)
}

\arguments{
  \item{x}{regressors (matrix); input to \code{lasso}}
  \item{y}{dependent variable; input to \code{lasso}}
  \item{K}{number of folds}
  \item{lambda.grid}{grid for lambda on which the search is conducted.}
  \item{post}{if \code{TRUE}, Post-Lasso estimation is done.}
  \item{intercept}{if \code{TRUE}, intercept is included.}
  \item{normalize}{if \code{TRUE}, regressors are normalized.}
}
\details{
If no grid for lambda is specified, a grid is automatically constructed with size of 100 grid points.
}
\value{
The function returns a list with the following components
\item{lambda.grid}{grid for lambda, as above.}
\item{cv}{CV curve at each value of lambda.grid}
\item{cv.error}{standard error of the CV curve}
\item{lambda.cv}{lambda values which minimize the cross validated mean square error; need not be unique.}
}

%\author{
%Martin Spindler
%}

\keyword{Lasso}
\keyword{Cross Validation}
\keyword{CV}
\keyword{lambda}
\keyword{penalty parameter}
