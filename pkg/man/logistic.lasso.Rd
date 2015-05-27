\name{logistic.lasso}
\alias{logistic.lasso}
\title{
logistic.lasso: Function for logistic Lasso estimation
}
\description{The function estimates the coefficients of a logistic Lasso regression with data-driven penalty. The method of the data-driven penalty can be chosen. The object which is returned is of the S3 class \code{logistic.lasso}
}
\usage{logistic.lasso(x, y, post = TRUE, intercept = TRUE, normalize = TRUE, control = list(c = 1.1, gamma = 0.1, numIter = 15, tol = 10^-5, lambda = "standard", numSim=10000, numfolds=10, lambda.start = NULL, threshold=NULL))
}

\arguments{
  \item{x}{regressors (matrix)}
  \item{y}{dependent variable (vector or matrix)}
  \item{post}{logical. If \code{TRUE}, post-lasso estimation is conducted.}
  \item{intercept}{logical. If \code{TRUE}, intercept is included which is not penalized.}
  \item{normalize}{logical. If \code{TRUE}, design matrix \code{x} is scaled.}
  \item{control}{list with control values. \code{c} and \code{gamma} constants for the penalty, \code{numIter} number of iterations for the algorithm for the estimation of the variance and data-driven penalty, ie. loadings, \code{tol} tolerance for improvement of the estimated variances,  \code{lambda.choice} method for penalty choice ("X-dependent", "X-independent", "standard", "none"), \code{numSim} number of simulations for the "X-dependent" method, \code{lambda.start} initial penalization value, and the \code{threshold} which is applied to the absolute values of the final estimates. 
  %\code{threshold} is applied to the estimated lasso coefficients. Absolute values below the threshold are set to zero.
  }
}
\details{
The function estimates the coefficients of a Logistic Lasso regression with data-driven penalty. The method of the data-driven penalty can be chosen ("X-dependent", "X-independent", "standard", "none"). The description of the methods is given in \code{\link{lasso}}, with the minor difference that here also for the  "X-dependent" and "X-independent" method data-driven loadings are employed (instead of an iteratively estimated variance factor). The option \code{post=TRUE} conducts post-lasso estimation, i.e. a refit of the model with the selected variables.
}
\value{
\code{logistic.lasso} returns an object of class \code{logistic.lasso}
An object of class "logistic.lasso" is a list containing at least the following components:
\item{coefficients}{parameter estimates}
\item{index}{index of selected variables (logicals)}
\item{lambda}{data-driven penalty term for each variable, product of lambda0 (the penalization paramter) and the loadings}
\item{lambda0}{penalty term}
\item{loadings}{loading for each regressor}
\item{residuals}{estimated residuals}
\item{sigma}{root ot the variance of the residuals}
\item{iter}{number of iterations}
\item{call}{function call}
\item{options}{options}
}

\references{
TBD
}
%\author{
%% Martin Spindler
%}


\keyword{logistic lasso}
\keyword{lasso}
\keyword{logistic regression}