# DGP
set.seed(123456)
library(AER)
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
n <- 10000
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
ts <- hdm::tsls(y,d,cbind(1,x),z)
iv <- ivreg(y~d+x|z+x)
#

output1 <- rlassoATE(x,d,y, bootstrap=NULL, nRep=500)
output2 <- rlassoATET(x,d,y, bootstrap=NULL, nRep=500)
output3 <- rlassoLATE(x,d,y,z, bootstrap=NULL, nRep=500)
output4 <- rlassoLATET(x,d,y,z, bootstrap=NULL, nRep=500)

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


############################################################################################################################

data(pension)
dat <- pension[,c("age", "inc", "fsize", "educ", "marr","twoearn","db","pira")]
dat$age <- (dat$age-25)/(64-25)
dat$inc <- (dat$inc+2652)/(242124+2652)
dat$fsize <- dat$fsize/13
dat$educ <- dat$educ/18
#degree of polynomials for age, eductaion, income, family, size
nAP <- 4
nIP <- 8
nEP <- 2
nFP <- 2
ageP <- poly(dat$age, degree=nAP)
colnames(ageP) <- paste("age", 1:4)
incP <- poly(dat$inc, degree=nIP)
colnames(incP) <- paste("inc", 1:8)
educP <- poly(dat$educ, degree=nEP)
colnames(educP) <- paste("inc", 1:2)
fsizeP <- poly(dat$fsize, degree=nFP)
colnames(fsizeP) <- paste("fsizeP", 1:2)
x <- as.matrix(cbind(ageP,incP,educP,fsizeP, pension[,c("marr","twoearn","db","pira")]))
d <- pension$p401 # treatment variable
z <- pension$e401 # instrument
y <- pension$tw # outcome variable

output1a <- rlassoATE(x,d,y, bootstrap="wild", nRep=100)
output2a <- rlassoATET(x,d,y, bootstrap="wild", nRep=100)
output3a <- rlassoLATE(x,d,y,z, bootstrap="wild", nRep=100)
output4a <- rlassoLATET(x,d,y,z, bootstrap="wild", nRep=100)



######################################################################################
library(R.matlab)
data <- readMat("C:\\Users\\Martin\\Dropbox\\HDM\\Matlab Code\\ProgEvalExample\\ProgEvalspec2.mat")
y <- data$y
d <- data$d
z <- data$z
x <- data$x
output1 <- rlassoLATE(x,d,y,z, bootstrap="wild", nRep=100)
xs <- cbind(1,x)
output2 <- rlassoLATE(xs,d,y,z, bootstrap="wild", nRep=100)
data <- readMat("C:\\Users\\Martin\\Dropbox\\HDM\\Matlab Code\\ProgEvalExample\\ProgEval2.mat")
y <- data$y
d <- data$d
z <- data$z
x <- data$x


data <- readMat("C:\\Users\\Martin\\Dropbox\\HDM\\Matlab Code\\ProgEvalExample\\ProgEvalspec2_intermediateresults.mat")
data1 <- readMat("C:\\Users\\Martin\\Dropbox\\HDM\\Matlab Code\\ProgEvalExample\\ProgEvalspec2_intermediateresults_x.mat")

y <- data$y
d <- data$d
z <- data$z
x <- data1$x
output1 <- rlassoLATE(x,d,y,z, bootstrap="none")
