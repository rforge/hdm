context("Test rlasso")

# DGP
n <- 100; p <- 10
X <- matrix(rnorm(n*p), ncol=p)
y <- X%*%rep(c(1,0), c(5,5)) + rnorm(n)

test_that("test inputs type", {
   expect_is(res <- rlasso(X,y), "rlasso")
   expect_is(rlasso(X,as.vector(y)), "rlasso")
   expect_is(rlasso(X,as.matrix(y)), "rlasso")
   #expect_is(rlasso(as.list(X),y), "rlasso")
   expect_is(rlasso(as.data.frame(X),y), "rlasso")
          })

#  test_that("lasso - choice of lambda", {
#    expect_is(lasso(X,y), "lasso")
#    expect_is(lasso(X,y, penalty=list(method = "standard", lambda.start = NULL, c = 1.1, gamma = 0.1)), "lasso")
#    expect_is(lasso(X,y, penalty=list(method = "standard")), "lasso")
#    expect_is(lasso(X,y, penalty=list(method = "X-dependent", numSim=5000)), "lasso")
#    expect_is(lasso(X,y, penalty=list(method = "X-dependent", c=1.2)), "lasso")
#    expect_is(lasso(X,y, penalty=list(method = "X-independent", c=1.2)), "lasso")
#    expect_is(lasso(X,y, penalty=list(method = "none", lambda.start=rep(1,p))), "lasso")
#    expect_is(lasso(X,y, penalty=list(method = "CV")), "lasso")
#  })
