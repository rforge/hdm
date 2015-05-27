\name{methods.lasso}
\alias{methods.lasso}
\alias{print.lasso}
\alias{summary.lasso}
\alias{predict.lasso}
\alias{model.matrix.lasso}
\title{
Methods for displaying information about \code{lasso} objects
}
\description{
Objects of class \code{lasso} are construced by lasso.formula or lasso.default.
\code{print.lasso} prints and displays some information about fitted \code{lasso} objects. 
\code{summary.lasso} summarizes information of a fitted \code{lasso} object.
\code{predict.lasso} predicts values based on a lasso object.
\code{model.matrix.lasso} constructs the model matrix of a lasso object.
}
\usage{
\method{print}{lasso}(object, all=TRUE, digits)
\method{summary}{lasso}(object, all=TRUE, digits)
\method{predict}{lasso}(object, newdata)
\method{model.matrix}{lasso}(object)
}
\arguments{
  \item{object}{
an object of class lasso
}
  \item{all}{logical, \code{TRUE} or \code{FALSE}. If \code{FALSE} only the coefficients of the selected variables will be displayed.
}
  \item{digits}{number of digits displayed
}
  \item{newdata}{
a data frame. If \code{newdata = NULL} then original data will be used.
}

  \item{\dots}{arguments passed to the print function}
}
\keyword{methods}
\keyword{lasso}