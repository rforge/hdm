#' hdm: High-Dimensional Metrics
#'
#' This packages implements high-dimensional methods with a focus on
#' Econometrics.
#'
#' \tabular{ll}{ Package: \tab hdm\cr Type: \tab Package\cr Version: \tab
#' 0.1\cr Date: \tab 2015-05-25\cr License: \tab GPL-3\cr } The package
#' includes functions for fitting a heteroskedastic robust Lasso regression and
#' for instrumental variable (IV) and treatment effect estimation in a
#' high-dimensional setting. Moreover, the methods enable valid post-selection
#' inference.
#'
#' @name hdm-package
#' @aliases hdm-package hdm
#' @docType package
#' @author Victor Chernozhukov, Christian Hansen, Martin Spindler\cr
#'
#' Maintainer: Martin Spindler <spindler@mea.mpisoc.mpg.de>
#' @references TBD
#' @keywords package Lasso Instrumental Variables Endogeneity Microeconometrics
#' Program Evaluation treatment effects
NULL

#' Growth data set
#'
#' Data set of growth compiled by Barro Lee.
#'
#' The data set contains growth data of Barro-Lee. The Barro Lee data consists
#' of a panel of 138 countries for the period 1960 to 1985. The depenent
#' variable is national growth rates in GDP per capita for the periods
#' 1965-1975 and 1975-1985. The growth rate in GDP over a period from t1 to t2
#' is commonly defined as \eqn{\log(GDP_{t_1}/GDP_{t_2})}. The number of covariates is p=62.
#' The number of complete observations is 90.
#'
#' @source The full data set and further details can be found at
#' \url{http://www.nber.org/pub/barro.lee}, \url{http://www.barrolee.com}, and,
#' \url{http://www.bristol.ac.uk//Depts//Economics//Growth//barlee.htm}.
#'
#' @name Growth Data
#' @aliases Growth Example GrowthData GDP
#' @docType data
#' @format \describe{Dataframe with the following variables:
#' \item{outcome}{dependent variable: national growth rates in GDP per capita
#' for the periods 1965-1975 and 1975-1985} \item{x}{covariates which might
#' influence growth} }
#' @references R.J. Barro, J.W. Lee (1994). Data set for a oanel of 139
#' countries. NBER.
#'
#' @references R.J. Barro, X. Sala-i-Martin (1995). Economic Growth. McGrwa-Hill, New York.
#' @keywords datasets Grwoth GDP
#' @examples
#' data(GrwothData)
NULL


#' Pension 401(k) data set
#'
#' Data set on financial wealth and 401(k) plan participation
#'
#' The sample is drawn from the 1991 Survey of Income and Program Participation
#' (SIPP) and consists of 9,915 observations. The observational units are
#' household reference persons aged 25-64 and spouse if present. Households are
#' included in the sample if at least one person is employed and no one is
#' self-employed. The data set was analyzed in Chernozhukov and Hansen (2004)
#' and Belloni et al. (2014) where further details can be found. They examine
#' the effects of 401(k) plans on wealth using data from the Survey of Income
#' and Program Participation. Using 401(k) eligibility as an instrument for
#' 401(k) participation
#'
#' @name pension
#' @aliases 401(k) plans wealth data pension
#' @docType data
#' @format \describe{Dataframe with the following variables (amongst others): 
#' \item{p401}{participation in 401(k)} \item{e401}{eligibility for 401(k)}
#' \item{a401}{401(k) assets}
#' \item{tw}{total wealth (in US $)} \item{tfa}{financial assets (in US $)}
#' \item{net_tfa}{net financial assets (in US $)}
#' \item{nifa}{non-401k financial assets (in US $)}
#' \item{net_nifa}{net non-401k financial assets}
#' \item{net_n401}{net non-401(k) assets (in US $)}
#' \item{ira}{individual retirement account (IRA)}
#' \item{inc}{income (in US $)}
#' \item{age}{age}
#' \item{fsize}{family size}
#' \item{marr}{married}
#' \item{pira}{participation in IRA} \item{db}{defined benefit pension}
#' \item{hown}{home owner} \item{educ}{education (in years)} \item{male}{male}
#' \item{twoearn}{two earners} %\item{i1-i7}{} %\item{a1-a5}{}
#' \item{nohs, hs, smcol, col}{dummies for education: no high-school, high-school, some college, college}
#' \item{hmort}{home mortage (in US $)}
#' \item{hequity}{home equity (in US $)} \item{hval}{home value (in US $)}}
#' @references V. Chernohukov, C. Hansen (2004). The impact of 401(k)
#' participation on the wealth distribution: An instrumental quantile
#' regression analysis. The Review of Economic and Statistics 86 (3), 735--751.
#'
#' @references A. Belloni, V. Chernozhukov, I. Fernandez-Val, and C. Hansen (2014). Program
#' evaluation with high-dimensional data. Working Paper.
#' @keywords datasets 401(k) pension
#' @examples
#' data(pension)
NULL





#' Riboflavin data set
#'
#' Dataset of riboflavin production by Bacillus subtilis containing \eqn{n=71}
#' observations of \eqn{p=4088} predictors (gene expressions) and a
#' one-dimensional response (riboflavin production).
#'
#' Data was made available by DSM (Switzerland). The data set was introduced in
#' B\"uhlmann et al. (2014).
#'
#' @name riboflavin
#' @docType data
#' @format \describe{ \item{y}{Log-transformed riboflavin production rate
#' (original name: q_RIBFLV).} \item{x}{(Co-)variables measuring the logarithm
#' of the expression level of 4088 genes.} }
#' @references P. B\"uhlmann, M. Kalisch, L- Meier (2014). High-dimensional
#' statistics with a view towards applications in biology. Annual Review of
#' Statistics and its Applications Vol. 1, 255--278.
#' @keywords datasets
#' @examples
#' data(riboflavin)
NULL

#' AJR data set
#'
#' Dataset on settler mortality
#'
#' Data set was analyzed in Acemoglu et al. (2001). A detailed decription of the data can be found at \url{http://economics.mit.edu/faculty/acemoglu/data/ajr2001} 
#'
#' @name AJR
#' @docType data
#' @format \describe{ 
#' \item{Mort}{Settler mortality}
#' \item{logMort}{logarithm of Mort}
#' \item{Latitude}{Latitude}
#' \item{Latitude2}{Latitude^2}
#' \item{Africa}{Africa}
#' \item{Asia}{Asia}
#' \item{Namer}{North America}
#' \item{Samer}{South America}
#' \item{Neo}{Neo-Europes}
#' \item{GDP}{GDP}
#' \item{Exprop}{Average protection against expropriation risk} 
#' }
#' @references D. Acemoglu, S. Johnson, J. A. Robinson  (2001). Colonial Origins of Comparative Development: An Empirical Investigation.
#' American Economic Review, 91, 1369--1401.
#' @keywords datasets
#' @examples
#' data(AJR)
NULL

#' Eminent Domain data set
#'
#' Dataset on judicial eminent domain decisions.
#'
#' Data set was analyzed in Belloni et al. (2012).  They estimate the effect of judicial eminent domain decisions on economic outcomes with instrumental variable (IV) in a setting high a large set of potential IVs. 
#' A detailed decription of the data can be found at 
#' \url{https://www.econometricsociety.org/publications/econometrica/2012/11/01/sparse-models-and-methods-optimal-instruments-application} 
#'
#' @name Eminent Domain
#' @docType data
#' @format \describe{ 
#' \item{y}{economic outcome variable, here the CS Home Price Index}
#' \item{x}{set of exogenous variables}
#' \item{d}{eminent domain decisions}
#' \item{z}{set of potential instruments}
#' }
#' @references D. Belloni, D. Chen, V. Chernozhukov and C. Hansen (2012).
#' Sparse models and methods for optimal instruments with an application to
#' eminent domain. \emph{Econometrica} 80 (6), 2369--2429.
#' @keywords datasets
#' @examples
#' data(EminentDomain)
NULL

