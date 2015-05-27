\name{DGP.HCHIV}
\alias{DGP.HCHIV}

\title{
Data Generating Process for Estimation in (Linear) Models with many Controls and Instruments
}
\description{
The function implements the data generating process used in Chernozhukov et al. (2015).
}
\usage{
DGP.HCHIV(n = 250, px = 300, pz = 150, alpha = 1)
}

\arguments{
  \item{n}{sample size}
  \item{px}{number of control variables}
  \item{pz}{number of instrumental variables}
  \item{alpha}{coefficient of the treatment effect}
}

\details{
The functions implements the DGP used in Chernzhukov et al. (2015).
}
\value{
The functions returns a list.
}
\references{
V. Chernozhukov, C. Hansen, M. Spindler (2015). Post-selection and post-regularization inference in linear models with many controls and instruments. American Economic Review: Paper & Proceedings 105(5), 1--7.
}
%\author{
%  Martin Spindler
%}

\keyword{DGP}
\keyword{Many controls and instruments}
\keyword{Instrumental Variable Estimation}
