\name{Growth Data}
\alias{Growth Example}
\alias{GrowthData}
\docType{data}
\title{Growth data set}
\description{Data set of growth compiled by Barro Lee.} 
\usage{data(GrowthData)}
\format{
  \describe{Dataframe with the following variables:
    \item{outcome}{dependent variable: national growth rates in GDP per capita for the periods 1965-1975 and 1975-1985}
    \item{x}{covariates which might influence growth}
  }
}
\details{The data set contains growth data of Barro-Lee. The Barro Lee data consists of a panel of 138 countries for the period 1960 to 1985. The depenent variable is national growth rates in GDP per capita for the periods 1965-1975 and 1975-1985. The growth rate in GDP over a period from t1 to t2 is commonly defined as log(GDP_t1/GDP_t2). The number of covariates is p=62. The number of complete observations is 90.

The full data set and further details can be found at http://www.nber.org/pub/barro.lee, http://www.barrolee.com, and, http://www.bristol.ac.uk/Depts/Economics/Growth/barlee.htm. }

\references{R.J. Barro, J.W. Lee (1994). Data set for a oanel of 139 countries. NBER.

R.J. Barro, X. Sala-i-Martin (1995). Economic Growth. McGrwa-Hill, New York.
}

\examples{data(GrwothData)}
\keyword{datasets}
