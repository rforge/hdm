\name{methods.logistic.lasso}
\alias{methods.logistic.lasso}
\alias{predict.logistic.lasso}
\title{
Methods for objects of class \code{logistic.lasso} 
}
\description{
Objects of class \code{logistic.lasso} are construced by \code{\link{logistic.lasso}}
\code{predict.lasso} predicts values based on a lasso object.
}
\usage{
%\method{print}{lasso}(object, all=TRUE, digits)
%\method{summary}{lasso}(object, all=TRUE, digits)
\method{predict}{logistic.lasso}(object, newdata)
}
\arguments{
  \item{object}{an object of class logistic.lasso}
%  \item{all}{logical, \code{TRUE} or \code{FALSE}. If \code{FALSE} only the coefficients of the selected variables will be displayed.

%  \item{digits}{number of digits displayed
%}
  \item{newdata}{a data frame. If \code{newdata = NULL} then original data will be used.}

%  \item{\dots}{arguments passed to the print function}
}
\keyword{methods}
\keyword{logistic.lasso}