\name{LQplot}
\alias{LQplot}
\title{Function for plotting local quantile treatment effects (LQTE) and LQEs of the treated (LQTE)}
\description{
The functions plots local quantile treatment effects (LQTE) and LQEs of the treated (LQTE) which are extracted from an object returned from the function \code{\link{ProgEval}}. \code{ggplot2} is used for plotting.
}
\usage{
LQplot(object, title = "", ylab = "")
}
\arguments{
  \item{object}{object returned from the function \code{\link{ProgEval}}}
  \item{title}{main title for the plot}
  \item{ylab}{title for the y-axis}
}
\details{
More details are given in the function \code{\link{ProgEval}}.
}
\value{
Plots for the LQTE and the LQTT are returned.
}
%\author{
% Martin Spindler
%}
\keyword{plot}
\keyword{LQTE}
\keyword{LQTT}