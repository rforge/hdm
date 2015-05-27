\name{DGP.HC}
\alias{DGP.HC}

\title{
Data Generating Process for inference on treatment effects after selection amongst high-dimensional controls
}
\description{
The function implements the data generating process used in Belloni et al. (2014)
}
\usage{
DGP.HC(n = 250, p = 100, alpha = 0.5, design = 1, R2 = c(0.5, 0.5))
}

\arguments{
  \item{n}{sample size}
  \item{p}{number of variables of exogenous variables in the outcome equation and determing the treatment 
}
  \item{alpha}{coefficient of the treatment variable
}
  \item{design}{design=1 assumes homoscedastic error in both equations; design=2 is the heteroscedastic case.
}
  \item{R2}{vector giving the desired population values for the reduced form R2's
}
}
\details{
The DGP is described in detail in Belloni et al. (2014)
}
\value{
It returns a list with the following components
\item{y}{outcome variable}
\item{d}{treatment variable}
\item{X}{X regressors (matrix)}
}
\references{
A. Belloni, V. Chernozhukov, C. Hansen (2014). Inference on treatment effects after selection among high-dimensional controls. The Review of Economic Studies 81(2), 608--650.}
\keyword{DGP}
\keyword{High-dimensional controls}
