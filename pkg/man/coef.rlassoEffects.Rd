% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rlassoEffects.R
\name{coef.rlassoEffects}
\alias{coef.rlassoEffects}
\title{Coefficients from S3 objects \code{rlassoEffects}}
\usage{
\method{coef}{rlassoEffects}(object, complete = TRUE,
  selection.matrix = FALSE, include.targets = FALSE, ...)
}
\arguments{
\item{object}{an object of class \code{rlassoEffects}, usually a result of a call \code{rlassoEffect} or \code{rlassoEffects}.}

\item{complete}{general option of the function \code{coef}.}

\item{selection.matrix}{if TRUE, a selection matrix is returned that indicates the selected variables from each auxiliary regression. 
Default is set to FALSE.}

\item{include.targets}{if FALSE (by default) only the selected control variables are listed in the \code{selection.matrix}. If set to TRUE, 
the selection matrix will also indicate the selection of the target coefficients that are specified in the  \code{rlassoEffects} call.}

\item{...}{further arguments passed to functions coef or print.}
}
\description{
Method to extract coefficients from objects of class \code{rlassoEffects}
}
\details{
Printing coefficients and selection matrix for S3 object \code{rlassoEffects}
}
