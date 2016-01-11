context("Test rlassologit")
library(hdm)
set.seed(2)
n <- 250
p <- 100
px <- 10

X <- matrix(rnorm(n*p), ncol=p)
x <- X[,1,drop=F]
beta <- c(rep(2,px), rep(0,p-px))
P <- exp(X %*% beta)/(1+exp(X %*% beta))
y <- numeric(length=250)
for(i in 1:n){
y[i] <- sample(x=c(1,0),size=1,prob=c(P[i],1-P[i]))
}
frame <- as.data.frame(cbind(y,X))
colnames(frame) <- c("Res",paste("V",1:100,sep=""))


test_that("rlassologit - Input check x and y",{
  expect_is(rlassologit(X,y),"rlassologit")
  expect_is(rlassologit(X,y,intercept=F),"rlassologit")
  expect_is(rlassologit(y~X),"rlassologit")
  expect_is(rlassologit(Res~.,data=frame),"rlassologit")
  expect_equivalent(rlassologit(X,y)$coef,rlassologit(y~X)$coef)
  expect_equivalent(rlassologit(X,y)$coef,rlassologit(Res~.,data=frame)$coef)
  expect_is(rlassologit(X,as.vector(y)),"rlassologit")
  expect_error(rlassologit(X,as.data.frame(y)))
  expect_error(rlassologit(X,as.list(y)))
})

test_that("rlassologit - Input check penalty and control",{
  expect_is(rlassologit(X,y,control=list(numIter=15)),"rlassologit")
  expect_is(rlassologit(X,y,control=list(tol=10^-5)),"rlassologit")
  expect_is(rlassologit(X,y,penalty=list(c=1.1,gamma=0.1)),"rlassologit")
  expect_is(rlassologit(X,y,penalty=list(c=1.1)),"rlassologit")
}) 

test_that("rlassologit - check methods",{
  expect_that(summary(rlassologit(X,y)),not(throws_error()))
  expect_that(print(rlassologit(X,y)),not(throws_error()))
  expect_that(model.matrix(rlassologit(X,y)),not(throws_error()))
  expect_that(model.matrix(rlassologit(y~X)),not(throws_error()))
  expect_that(model.matrix(rlassologit(Res~.,data=frame)),not(throws_error()))
  expect_that(predict(rlassologit(X,y)),not(throws_error()))
  expect_that(predict(rlassologit(X,y,intercept=F)),not(throws_error()))
  expect_that(predict(rlassologit(y~X)),not(throws_error()))
  expect_that(predict(rlassologit(Res~.,data=frame)),not(throws_error()))
  expect_that(predict(rlassologit(X,y),as.data.frame(2*X)),not(throws_error()))
  expect_that(predict(rlassologit(y~X),as.data.frame(2*X)),not(throws_error()))
  expect_that(predict(rlassologit(Res~.,data=frame),as.data.frame(2*X)),not(throws_error()))
  expect_error(predict(rlassologit(Res~.,data=frame),as.data.frame(2*X[,1:10])))
  expect_that(predict(rlassologit(X,y),type="link"),not(throws_error()))
})
