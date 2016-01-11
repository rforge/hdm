context("Test rlasso")
library(hdm)
set.seed(2)
n <- 250
p <- 100
px <- 10

X <- matrix(rnorm(n*p), ncol=p)
x <- X[,1,drop=F]
beta <- c(rep(2,px), rep(0,p-px))
y <- X %*% beta + rnorm(n)
frame <- as.data.frame(cbind(y,X))
colnames(frame) <- c("Res",paste("V",1:100,sep=""))

test_that("rlasso - Input check x and y",{
  expect_is(rlasso(X,y),"rlasso")
  expect_is(rlasso(y~X),"rlasso")
  expect_is(rlasso(Res~.,data=frame),"rlasso")
  expect_equivalent(rlasso(X,y)$coef,rlasso(y~X)$coef)
  expect_equivalent(rlasso(X,y)$coef,rlasso(Res~.,data=frame)$coef)
  expect_is(rlasso(as.data.frame(X),y),"rlasso")
  expect_is(rlasso(X,as.vector(y)),"rlasso")
  expect_error(rlasso(X,as.data.frame(y)))
  expect_error(rlasso(X,as.list(y)))
})

test_that("rlasso - Input check x and y: x one-dimensional",{
  expect_is(rlasso(x,y),"rlasso")
  expect_error(rlasso(as.vector(x),y))
  expect_is(rlasso(as.data.frame(x),y),"rlasso")
  expect_is(rlasso(x,as.matrix(y)),"rlasso")
  expect_error(rlasso(x,as.list(y)))
  expect_error(rlasso(as.list(x),y))
})

test_that("rlasso - Input check penalty and control",{
  expect_is(rlasso(X,y,control=list(numIter=15)),"rlasso")
  expect_is(rlasso(X,y,control=list(tol=10^-5)),"rlasso")
  expect_is(rlasso(X,y,penalty=list(method="standard")),"rlasso")
  expect_is(rlasso(X,y,penalty=list(method="standard",c=1.1,gamma=0.1)),"rlasso")
  expect_is(rlasso(X,y,penalty=list(method="X-independent")),"rlasso")
  expect_error(rlasso(X,y,penalty=list(method="none")))
  expect_is(rlasso(X,y, penalty=list(method = "none", lambda.start=1)), "rlasso")
  expect_is(rlasso(X,y,penalty=list(method="X-dependent")),"rlasso")
  expect_is(rlasso(X,y,penalty=list(method="X-dependent",numSim=150)),"rlasso")
})

# test_that("rlasso - check methods",{
#   expect_that(summary(rlasso(X,y)),not(throws_error()))
#   expect_that(print(rlasso(X,y)),not(throws_error()))
#   expect_that(model.matrix(rlasso(X,y)),not(throws_error()))
#   expect_that(model.matrix(rlasso(y~X)),not(throws_error()))
#   expect_that(model.matrix(rlasso(Res~.,data=frame)),not(throws_error()))
#   expect_that(predict(rlasso(X,y)),not(throws_error()))
#   expect_that(predict(rlasso(X,y,intercept=F)),not(throws_error()))
#   expect_that(predict(rlasso(y~X)),not(throws_error()))
#   expect_that(predict(rlasso(Res~.,data=frame)),not(throws_error()))
#   expect_that(predict(rlasso(X,y),as.data.frame(2*X)),not(throws_error()))
#   expect_that(predict(rlasso(y~X),as.data.frame(2*X)),not(throws_error()))
#   expect_that(predict(rlasso(Res~.,data=frame),as.data.frame(2*X)),not(throws_error()))
#   expect_error(predict(rlasso(Res~.,data=frame),as.data.frame(2*X[,1:10])))
# })

