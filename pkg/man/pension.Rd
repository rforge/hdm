\name{401(k) plans}
\alias{401(k) plans}
\alias{wealth data}
\alias{pension}
\docType{data}
\title{Pension 401(k) data set}
\description{Data set on financial wealth and 401(k) plan participation} 
\usage{data(pension)}
\format{Dataframe with the following variables (amongst others):
  \describe{
      \item{p401}{participation in 401(k)}
      \item{e401}{eligibility for 401(k)}
      \item{a401}{401(k) assets}
      \item{tw}{total wealth (in US $)}
      \item{tfa}{total financial assets (in US $)}
      \item{net_tfa}{net financial assets (in US $)}
      \item{nifa}{non-401k financial assets (in US $)}
      \item{net_nifa}{net non-401k financial assets}
      \item{net_n401}{net non-401(k) assets (in US $)}
      \item{ira}{individual retirement account (IRA)}
      \item{inc}{income (in US $)}
      \item{age}{age}
      \item{fsize}{family size}
      \item{marr}{married}
      \item{pira}{participation in IRA}
      \item{db}{defined benefit pension}
      \item{hown}{home owner}
      \item{educ}{education (in years)}
      \item{male}{male}
      \item{twoearn}{two earners}
      \item{nohs, hs, smcol, col}{dummies for education: no high-school, high-school, some college, college}
      %\item{a1-a5}{}
      %\item{wtw1-wtw7}{}
      %\item{f1-f6}{}
      \item{hmort}{home mortage  (in US $)}
      \item{hequity}{home equity (in US $)}
      \item{hval}{home value (in US $)}

      
     
% "wgt"       "a401"       "tfa_he" "dum91" "nohs"       "hs"         "smcol"     "col"        "a"          "icat"       "ecat"    "f"          "zhat"       "wntfa1"     "wntfa2"     "wntfa3"     "wntfa4"     "wntfa5" "wntfa6"     "wntfa7"     "wnet_nifa1" "wnet_nifa2" "wnet_nifa3" "wnet_nifa4" "wnet_nifa5" "wnet_nifa6" "wnet_nifa7"     "wnet_n4014"
    
  }
}
\details{The sample
is drawn from the 1991 Survey of Income
and Program Participation (SIPP) and consists of 9,915 observations. The observational units are household reference persons aged 25-64 and spouse if present. Households are included in the sample if at least one person is employed and no one is self-employed. The data set was analyzed in Chernozhukov and Hansen (2004) and Belloni et al. (2014) where further details can be found. They examine the effects of 401(k) plans on wealth using data from the Survey of Income
and Program Participation. Using 401(k) eligibility as an instrument for 401(k) participation}

\references{V. Chernohukov, C. Hansen (2004). The impact of 401(k) participation on the wealth distribution: An instrumental quantile regression analysis. The Review of Economic and Statistics 86 (3), 735--751.

A. Belloni, V. Chernozhukov, I. Fernandez-Val, and C. Hansen (2014). Program evaluation with high-dimensional data. Working Paper.}

\examples{data(pension)}
\keyword{datasets}
