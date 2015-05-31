\name{lasso}
\alias{lasso}
\alias{lasso.default}
\alias{lasso.formula}
\alias{glm.lasso}
\title{lasso: Function for Lasso estimation under heterosceastic non-Gaussian disturbances}
\description{
The function estimates the coefficients of a Lasso regression with data-driven penalty under heteroskedasticity and non-Gaussian noise. The method of the data-driven penalty can be chosen. The object which is returned is of the S3 class \code{lasso}.
}

\usage{
\method{lasso}{default}(x, y, post = TRUE, intercept = TRUE, normalize = TRUE, control = list(c = 1.1, gamma = 0.1, numIter = 15, tol = 10^-5, lambda = "standard", numSim=10000, numfolds=10, lambda.start = NULL, threshold=NULL))
\method{lasso}{formula}(formula, data, post = TRUE, intercept = TRUE, normalize = TRUE, control = list(c = 1.1, gamma = 0.1, numIter = 15, tol = 10^-5,lambda = "standard",  numSim=10000, numfolds=10, lambda.start = NULL, threshold=NULL))

}
\arguments{
  \item{x}{regressors (matrix)}
  \item{y}{dependent variable (vector or matrix)}
  \item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted in the form \code{y~x}}
  \item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lasso is called.}
  \item{post}{logical. If \code{TRUE}, post-lasso estimation is conducted.}
  \item{intercept}{logical. If \code{TRUE}, intercept is included which is not penalized.}
  \item{normalize}{logical. If \code{TRUE}, design matrix \code{x} is scaled.}
  \item{control}{list with control values. \code{c} and \code{gamma} constants for the penalty, \code{numIter} number of iterations for the algorithm for the estimation of the variance and data-driven penalty, ie. loadings, \code{tol} tolerance for improvement of the estimated variances,  \code{lambda.choice} method for penalty choice ("X-dependent", "X-independent", "standard", "none", "cv"), \code{numfolds} number of folds for the "cross validation"  method, \code{numSim} number of simulations for the "X-dependent" method, \code{lambda.start} initial penalization value, and the \code{threshold} which is applied to the absolute values of the final estimates. 
  %\code{threshold} is applied to the estimated lasso coefficients. Absolute values below the threshold are set to zero.
  }
}

\details{
The function estimates the coefficients of a Lasso regression with data-driven penalty under heteroskedasticity and non-Gaussian noise. The method of the data-driven penalty can be chosen ("X-dependent", "X-independent", "standard", "none", "cv"). For details of the implementation of the Algorithm for estimation of the data-driven penalty, in particular the regressor-dependent loadings, we refer to Appendix A in Belloni et al.~(2012). When the option "none" is chosen (together with \code{lambda.start}),  lambda is set to \code{lambda.start} and the regressor-dependent loadings are used. The options "X-dependent" and  "X-independent" are described in Belloni et al.~(2013). The method "cv" does cross validation to determine the "optimal" value for lambda (with regressor-dependent loadings). The option \code{nfolds} gives the number of folds used for cross validation. The option \code{post=TRUE} conducts post-lasso estimation, i.e. a refit of the model with the selected variables.
}

\value{
\code{lasso} returns an object of class \code{lasso}
An object of class "lasso" is a list containing at least the following components:
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

\examples{
## DGP
library(hdm)
library(MASS)
n <- 250
p <- 100
px <- 10
X <- matrix(rnorm(n*p), ncol=p)
beta <- c(rep(2,px), rep(0,p-px))
y <- 1 + X\%*\%beta + rnorm(n)
## Lasso estimation
lasso.reg <- lasso(x=X, y=y, post=TRUE, intercept=TRUE)
# Methods for Lasso
print(lasso.reg, all=FALSE)
summary(lasso.reg, all=FALSE)
yhat <- predict(lasso.reg)
Xnew <- matrix(rnorm(n*p), ncol=p)
yhat.new <- predict(lasso.reg, newdata=Xnew)
}

\references{
A. Belloni, D. Chen, V. Chernozhukov and C. Hansen (2012). Sparse models and methods for optimal instruments with an application to eminent domain. \emph{Econometrica} 80 (6), 2369--2429.

A. Belloni, V. Chernozhukov and C. Hansen (2013). Inference for high-dimensional sparse econometric models. In Advances in Economics and Econometrics: 10th World Congress, Vol. 3: Econometrics, Cambirdge University Press: Cambridge, 245-295.
}
%\author{
%Martin Spindler
%}

\keyword{Lasso}
\keyword{data-driven penalty}
\keyword{non-Gaussian}
\keyword{heteroscedasticity}
