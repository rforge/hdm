\name{Unique}
\alias{Unique}
\title{
Extract Unqiue Elements
}
\description{
unique returns a vector, data frame or array like y but with duplicate elements/rows removed. It is comparable to the R function \code{unique} but returns additional information.
}
\usage{
Unique(y)
}
\arguments{
  \item{y}{vector or data frame}
}
\value{
The functions retunrs a list with three components
\item{data}{same type of \code{y} but with duplicated values removed}
\item{index.unique}{index of the unique values}
\item{index.duplicated}{a list which gives for each unique value the index of the emelents with the same value}
}

%\author{
%%%  Martin Spindler
%}
\keyword{unique}
\keyword{index}
