# DGP
set.seed(12345)
library(AER)
library(hdm)
#library(hdm2)
# data <- DGP.HCHIV()
#
# y <- data$y
# d <- as.integer(data$d >=0)
# Z <- data$Z[,1]
# z <- Z >= 0
# z <- matrix(as.integer(z), nrow=250)
# x <- data$X

## DGP LATE
n <- 500
eps <- rnorm(n)
z <- sample(c(0,1), n, replace=TRUE, prob=c(0.5,0.5))
u <- runif(n)
ctype <- 1*(u< quantile(u, probs=0.25)) + 2*(u> quantile(u, probs=0.25) & u< quantile(u, probs=0.75)) + 3*(u> quantile(u, probs=0.75))
d <- vector("numeric", length=n)
f <- function(z,type) {
  if (type==2) return(z)
  if (type==3) return(1)
  if (type==1) return(0)
}

d <- mapply(f, z, ctype)
beta <- ctype*2
x <- matrix(rnorm(2*n), ncol=2) #x <- rep(1,n)
y <- x%*%c(1,1) + beta*d +  eps

###
ts <- tsls(y,d,x,z, intercept=TRUE)
iv <- ivreg(y~d+x|z+x)
#

#output1 <- rlassoATE(x,d,y, bootstrap=NULL, nRep=500)
#output2 <- rlassoATET(x,d,y, bootstrap=NULL, nRep=500)
#output4 <- rlassoLATET(x,d,y,z, bootstrap=NULL, nRep=500)

output3a <- rlassoLATE(x,d,y,z, bootstrap=NULL, nRep=500, normalize = FALSE, intercept=FALSE, post=FALSE)
output3b <- rlassoLATE(x,d,y,z, bootstrap=NULL, nRep=500, normalize = FALSE, intercept=FALSE, post=TRUE)
output3c <- rlassoLATE(x,d,y,z, bootstrap=NULL, nRep=500, normalize = FALSE, intercept=TRUE, post=FALSE)
output3d <- rlassoLATE(x,d,y,z, bootstrap=NULL, nRep=500, normalize = FALSE, intercept=TRUE, post=TRUE)
output3e <- rlassoLATE(x,d,y,z, bootstrap=NULL, nRep=500, normalize = TRUE, intercept=FALSE, post=FALSE)
output3f <- rlassoLATE(x,d,y,z, bootstrap=NULL, nRep=500, normalize = TRUE, intercept=FALSE, post=TRUE)
output3g <- rlassoLATE(x,d,y,z, bootstrap=NULL, nRep=500, normalize = TRUE, intercept=TRUE, post=FALSE)
output3h <- rlassoLATE(x,d,y,z, bootstrap=NULL, nRep=500, normalize = TRUE, intercept=TRUE, post=TRUE)

output3i <- rlassoLATE(x,d,y,z, bootstrap="normal", nRep=250, normalize = FALSE, intercept=FALSE, post=FALSE)
output3j <- rlassoLATE(x,d,y,z, bootstrap="normal", nRep=250, normalize = FALSE, intercept=FALSE, post=TRUE)
output3k <- rlassoLATE(x,d,y,z, bootstrap="normal", nRep=250, normalize = FALSE, intercept=TRUE, post=FALSE)
output3l <- rlassoLATE(x,d,y,z, bootstrap="normal", nRep=250, normalize = FALSE, intercept=TRUE, post=TRUE)
output3m <- rlassoLATE(x,d,y,z, bootstrap="wild", nRep=250, normalize = TRUE, intercept=FALSE, post=FALSE)
output3n <- rlassoLATE(x,d,y,z, bootstrap="wild", nRep=250, normalize = TRUE, intercept=FALSE, post=TRUE)
output3o <- rlassoLATE(x,d,y,z, bootstrap="wild", nRep=250, normalize = TRUE, intercept=TRUE, post=FALSE)
output3p <- rlassoLATE(x,d,y,z, bootstrap="wild", nRep=250, normalize = TRUE, intercept=TRUE, post=TRUE)


output4a <- rlassoLATET(x,d,y,z, bootstrap=NULL, nRep=500, normalize = FALSE, intercept=FALSE, post=FALSE)
output4b <- rlassoLATET(x,d,y,z, bootstrap=NULL, nRep=500, normalize = FALSE, intercept=FALSE, post=TRUE)
output4c <- rlassoLATET(x,d,y,z, bootstrap=NULL, nRep=500, normalize = FALSE, intercept=TRUE, post=FALSE)
output4d <- rlassoLATET(x,d,y,z, bootstrap=NULL, nRep=500, normalize = FALSE, intercept=TRUE, post=TRUE)
output4e <- rlassoLATET(x,d,y,z, bootstrap=NULL, nRep=500, normalize = TRUE, intercept=FALSE, post=FALSE)
output4f <- rlassoLATET(x,d,y,z, bootstrap=NULL, nRep=500, normalize = TRUE, intercept=FALSE, post=TRUE)
output4g <- rlassoLATET(x,d,y,z, bootstrap=NULL, nRep=500, normalize = TRUE, intercept=TRUE, post=FALSE)
output4h <- rlassoLATET(x,d,y,z, bootstrap=NULL, nRep=500, normalize = TRUE, intercept=TRUE, post=TRUE)

output4i <- rlassoLATET(x,d,y,z, bootstrap="normal", nRep=250, normalize = FALSE, intercept=FALSE, post=FALSE)
output4j <- rlassoLATET(x,d,y,z, bootstrap="normal", nRep=250, normalize = FALSE, intercept=FALSE, post=TRUE)
output4k <- rlassoLATET(x,d,y,z, bootstrap="normal", nRep=250, normalize = FALSE, intercept=TRUE, post=FALSE)
output4l <- rlassoLATET(x,d,y,z, bootstrap="normal", nRep=250, normalize = FALSE, intercept=TRUE, post=TRUE)
output4m <- rlassoLATET(x,d,y,z, bootstrap="wild", nRep=250, normalize = TRUE, intercept=FALSE, post=FALSE)
output4n <- rlassoLATET(x,d,y,z, bootstrap="wild", nRep=250, normalize = TRUE, intercept=FALSE, post=TRUE)
output4o <- rlassoLATET(x,d,y,z, bootstrap="wild", nRep=250, normalize = TRUE, intercept=TRUE, post=FALSE)
output4p <- rlassoLATET(x,d,y,z, bootstrap="wild", nRep=250, normalize = TRUE, intercept=TRUE, post=TRUE)

output5a <- rlassoATE(x,d,y, bootstrap=NULL, nRep=500, normalize = FALSE, intercept=FALSE, post=FALSE)
output5b <- rlassoATE(x,d,y, bootstrap=NULL, nRep=500, normalize = FALSE, intercept=FALSE, post=TRUE)
output5c <- rlassoATE(x,d,y, bootstrap=NULL, nRep=500, normalize = FALSE, intercept=TRUE, post=FALSE)
output5d <- rlassoATE(x,d,y, bootstrap=NULL, nRep=500, normalize = FALSE, intercept=TRUE, post=TRUE)
output5e <- rlassoATE(x,d,y, bootstrap=NULL, nRep=500, normalize = TRUE, intercept=FALSE, post=FALSE)
output5f <- rlassoATE(x,d,y, bootstrap=NULL, nRep=500, normalize = TRUE, intercept=FALSE, post=TRUE)
output5g <- rlassoATE(x,d,y, bootstrap=NULL, nRep=500, normalize = TRUE, intercept=TRUE, post=FALSE)
output5h <- rlassoATE(x,d,y, bootstrap=NULL, nRep=500, normalize = TRUE, intercept=TRUE, post=TRUE)

output5i <- rlassoATE(x,d,y, bootstrap="normal", nRep=250, normalize = FALSE, intercept=FALSE, post=FALSE)
output5j <- rlassoATE(x,d,y, bootstrap="normal", nRep=250, normalize = FALSE, intercept=FALSE, post=TRUE)
output5k <- rlassoATE(x,d,y, bootstrap="normal", nRep=250, normalize = FALSE, intercept=TRUE, post=FALSE)
output5l <- rlassoATE(x,d,y, bootstrap="normal", nRep=250, normalize = FALSE, intercept=TRUE, post=TRUE)
output5m <- rlassoATE(x,d,y, bootstrap="wild", nRep=250, normalize = TRUE, intercept=FALSE, post=FALSE)
output5n <- rlassoATE(x,d,y, bootstrap="wild", nRep=250, normalize = TRUE, intercept=FALSE, post=TRUE)
output5o <- rlassoATE(x,d,y, bootstrap="wild", nRep=250, normalize = TRUE, intercept=TRUE, post=FALSE)
output5p <- rlassoATE(x,d,y, bootstrap="wild", nRep=250, normalize = TRUE, intercept=TRUE, post=TRUE)

output6a <- rlassoATET(x,d,y, bootstrap=NULL, nRep=600, normalize = FALSE, intercept=FALSE, post=FALSE)
output6b <- rlassoATET(x,d,y, bootstrap=NULL, nRep=600, normalize = FALSE, intercept=FALSE, post=TRUE)
output6c <- rlassoATET(x,d,y, bootstrap=NULL, nRep=600, normalize = FALSE, intercept=TRUE, post=FALSE)
output6d <- rlassoATET(x,d,y, bootstrap=NULL, nRep=600, normalize = FALSE, intercept=TRUE, post=TRUE)
output6e <- rlassoATET(x,d,y, bootstrap=NULL, nRep=600, normalize = TRUE, intercept=FALSE, post=FALSE)
output6f <- rlassoATET(x,d,y, bootstrap=NULL, nRep=600, normalize = TRUE, intercept=FALSE, post=TRUE)
output6g <- rlassoATET(x,d,y, bootstrap=NULL, nRep=600, normalize = TRUE, intercept=TRUE, post=FALSE)
output6h <- rlassoATET(x,d,y, bootstrap=NULL, nRep=600, normalize = TRUE, intercept=TRUE, post=TRUE)

output6i <- rlassoATET(x,d,y, bootstrap="normal", nRep=260, normalize = FALSE, intercept=FALSE, post=FALSE)
output6j <- rlassoATET(x,d,y, bootstrap="normal", nRep=260, normalize = FALSE, intercept=FALSE, post=TRUE)
output6k <- rlassoATET(x,d,y, bootstrap="normal", nRep=260, normalize = FALSE, intercept=TRUE, post=FALSE)
output6l <- rlassoATET(x,d,y, bootstrap="normal", nRep=260, normalize = FALSE, intercept=TRUE, post=TRUE)
output6m <- rlassoATET(x,d,y, bootstrap="wild", nRep=260, normalize = TRUE, intercept=FALSE, post=FALSE)
output6n <- rlassoATET(x,d,y, bootstrap="wild", nRep=260, normalize = TRUE, intercept=FALSE, post=TRUE)
output6o <- rlassoATET(x,d,y, bootstrap="wild", nRep=260, normalize = TRUE, intercept=TRUE, post=FALSE)
output6p <- rlassoATET(x,d,y, bootstrap="wild", nRep=260, normalize = TRUE, intercept=TRUE, post=TRUE)


print(output1)
print(output2)
print(output3)
print(output4)

summary(output1)
summary(output2)
summary(output3)
summary(output4)

confint(output1)
confint(output2)
confint(output3)
confint(output4)


# ############################################################################################################################
# 
# data(pension)
# dat <- pension[,c("age", "inc", "fsize", "educ", "marr","twoearn","db","pira")]
# dat$age <- (dat$age-25)/(64-25)
# dat$inc <- (dat$inc+2652)/(242124+2652)
# dat$fsize <- dat$fsize/13
# dat$educ <- dat$educ/18
# #degree of polynomials for age, eductaion, income, family, size
# nAP <- 4
# nIP <- 8
# nEP <- 2
# nFP <- 2
# ageP <- poly(dat$age, degree=nAP)
# colnames(ageP) <- paste("age", 1:4)
# incP <- poly(dat$inc, degree=nIP)
# colnames(incP) <- paste("inc", 1:8)
# educP <- poly(dat$educ, degree=nEP)
# colnames(educP) <- paste("inc", 1:2)
# fsizeP <- poly(dat$fsize, degree=nFP)
# colnames(fsizeP) <- paste("fsizeP", 1:2)
# x <- as.matrix(cbind(ageP,incP,educP,fsizeP, pension[,c("marr","twoearn","db","pira")]))
# d <- pension$p401 # treatment variable
# z <- pension$e401 # instrument
# y <- pension$tw # outcome variable
# 
# output1a <- rlassoATE(x,d,y, bootstrap="wild", nRep=100)
# output2a <- rlassoATET(x,d,y, bootstrap="wild", nRep=100)
# output3a <- rlassoLATE(x,d,y,z, bootstrap="wild", nRep=100)
# output4a <- rlassoLATET(x,d,y,z, bootstrap="wild", nRep=100)
# 
# 
# 
# ######################################################################################
# library(R.matlab)
# library(hdm)
# data <- readMat("C:\\Users\\Martin\\Dropbox\\HDM\\Matlab Code\\ProgEvalExample\\ProgEvalspec2.mat")
# y <- data$y
# d <- data$d
# z <- data$z
# x <- data$x
# output1 <- rlassoLATE(x,d,y,z, bootstrap="wild", nRep=100)
# xs <- cbind(1,x)
# output2 <- rlassoLATE(xs,d,y,z, bootstrap="wild", nRep=100)
# data <- readMat("C:\\Users\\Martin\\Dropbox\\HDM\\Matlab Code\\ProgEvalExample\\ProgEval2.mat")
# y <- data$y
# d <- data$d
# z <- data$z
# x <- data$x
# 
# 
# data <- readMat("C:\\Users\\Martin\\Dropbox\\HDM\\Matlab Code\\ProgEvalExample\\ProgEvalspec2_intermediateresults.mat")
# data1 <- readMat("C:\\Users\\Martin\\Dropbox\\HDM\\Matlab Code\\ProgEvalExample\\ProgEvalspec2_intermediateresults_x.mat")
# 
# y <- data$y
# d <- data$d
# z <- data$z
# x <- data1$x
# output1 <- rlassoLATE(x,d,y,z, bootstrap="none")
# 
# 
# #############################
# 
# library(R.matlab)
# library(hdm)
# data <- readMat("E:\\R Package hdm\\hdm\\LATE_spec2.mat")
# y <- data$y
# d <- data$d
# z <- data$z
# x <- as.matrix(data$x)
# debug(rlassoLATE)
# #debug(lambdaCalculation)
# output1 <- rlassoLATE(x,d,y,z, bootstrap=NULL, post=TRUE, intercept=TRUE, normalize=FALSE)
# output2 <- rlassoLATE(x,d,y,z, bootstrap=NULL, post=TRUE, intercept=TRUE, normalize=TRUE)
# output3 <- rlassoLATET(x,d,y,z, bootstrap=NULL, post=TRUE, intercept=TRUE, normalize=FALSE)
# output4 <- rlassoLATET(x,d,y,z, bootstrap=NULL, post=TRUE, intercept=TRUE, normalize=TRUE)
# output5 <- rlassoATE(x,d,y,z, bootstrap=NULL, post=TRUE, intercept=TRUE, normalize=FALSE)
# output6 <- rlassoATE(x,d,y,z, bootstrap=NULL, post=TRUE, intercept=TRUE, normalize=TRUE)
# output7 <- rlassoATET(x,d,y,z, bootstrap=NULL, post=TRUE, intercept=TRUE, normalize=FALSE)
# output8 <- rlassoATET(x,d,y,z, bootstrap=NULL, post=TRUE, intercept=TRUE, normalize=TRUE)
# 
# ########################### Analysis simulation
# 
library(R.matlab)
rm(list=ls())
data <- readMat("E:\\R Package hdm\\Misc\\sim.mat")
data_int <- readMat("E:\\R Package hdm\\Misc\\sim_inter.mat")

y <- data$y
d <- data$d
z <- data$z
x <- as.matrix(data$xS)
debug(rlassoLATE)
output1 <- rlassoLATE(x,d,y,z, bootstrap=NULL, post=TRUE, intercept=TRUE, normalize=FALSE)
# output2 <- rlassoLATE(x,d,y,z, bootstrap=NULL, post=TRUE, intercept=TRUE, normalize=TRUE)
# output3 <- rlassoLATE(x,d,y,z, bootstrap=NULL, post=TRUE, intercept=FALSE, normalize=TRUE)
# output4 <- rlassoLATE(x,d,y,z, bootstrap=NULL, post=TRUE, intercept=TRUE, normalize=FALSE)
# 
# debug(rlassoATE)
# output1 <- rlassoATE(x,d,y, bootstrap=NULL, post=TRUE, intercept=FALSE, normalize=FALSE)
