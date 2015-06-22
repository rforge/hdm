# defining functions for the S3 class "lasso"

lasso <- function(x, ...) UseMethod("lasso") # definition generic function

lasso.default <- function(x, y, post=TRUE, intercept=TRUE, normalize=TRUE, 
                          control=list(c=1.1, gamma=.1, numIter=15, tol=10^-5, lambda="standard", numSim=10000, 
                                       numfolds=10, lambda.start=NULL, threshold=NULL)) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  if (is.null(colnames(x))) colnames(x) <- paste("V",1:p, sep="")
  ind.names <- 1:p
  eps <- 10^-9 # needed for scaling
  c <- control$c
  gamma <- control$gamma
  numIter <- control$numIter
  lambda.start <- control$lambda.start
  threshold <- control$threshold
  
  tol <- control$tol
  Ups0 <- Ups1 <- NULL # for the case not used
  # Intercept handling and scaling
  if (intercept) {
    meanx <- colMeans(x)
    x <- scale(x, meanx, FALSE)
    mu <- mean(y)
    y <- y - mu
  } else {
    meanx <- rep(0, p)
    mu <- 0
  }
  
  if (normalize) {
    normx <- sqrt(apply(x,2,var))
    ind <- which(normx < eps)
    if (length(ind)!=0) {
      x <- x[,-ind]
      normx <- normx[-ind]
      ind.names <- ind.names[-ind]
      p <- dim(x)[2]
      lambda.start <- lambda.start[-ind]
    }
    x <- scale(x, FALSE, normx)
  } else {
    normx <- rep(1,p)
  }
  
  XX <- crossprod(x)
  Xy <- crossprod(x,y) 
  
  
  # calculation lambda
  # cross validation
  if (control$lambda=="CV") {
    cv <- cv.lasso(x,y,K=control$numfolds, lambda.grid=NULL, post, intercept, normalize)
    lambda.cv <- max(cv$lambda.cv)
    control <- list(numIter=numIter, tol=10^-5, lambda="none", lambda.start=rep(lambda.cv,p))
    fit <- lasso(x, y, post, intercept, normalize, control=control)
    return(fit)
  }   
  # X-independent
  if (control$lambda=="X-independent") {
    lambda0 <- 2*c*sqrt(n)*qnorm(1-gamma/(2*p))
    lambda <- rep(2*c*sqrt(n)*qnorm(1-gamma/(2*p))*sqrt(var(y)),p)
  }   
  # X-dependent
  if (control$lambda=="X-dependent") {
    R <- control$numSim
    sim <- vector("numeric", length=R)
    for (l in 1:R) {
      g <- matrix(rep(rnorm(n), each=p), ncol=p, byrow=TRUE)
      sim[l] <- n*max(2*colMeans(x*g))
    }
    lambda0 <- c*quantile(sim, probs=1-gamma/2)
    sigma0 <- sqrt(var(y))
    lambda <- rep(c*quantile(sim, probs=1-gamma/2)*sigma0,p)
  }
  # "standard"
  if (control$lambda=="standard") { # choice paper
    lambda0 <- 2*c*sqrt(n)*sqrt(2*log(2*p*log(n)/gamma))
    e0 <- y - mean(y)
    Ups0 <- 1/sqrt(n)*sqrt(t(t(e0^2)%*%(x^2)))
    lambda <- lambda0*Ups0     
  }
  
  if (!is.null(lambda.start)) {
    #lambda0 <- lambda.start
    if (length(lambda.start)==1) {
      lambda.start <- rep(lambda.start,p)
    }
    lambda <-diag(p)%*%as.matrix(lambda.start)
  }
  
  if (control$lambda=="none") {
    lambda0 <- lambda.start
    e0 <- y - mean(y)
    Ups0 <- 1/sqrt(n)*sqrt(t(t(e0^2)%*%(x^2)))
    lambda <- lambda0*Ups0 
  }
  
  
  # calculation first parameters
  coefTemp <- LassoShooting.fit(x, y, lambda, XX=XX, Xy=Xy)$coefficients
  coefTemp[is.na(coefTemp)] <- 0
  ind1 <- (abs(coefTemp) > 0)
  x1 <- as.matrix(x[,ind1, drop=FALSE])
  if (dim(x1)[2]==0) {
    print("No variables selected!")
    #est <- list(coefficients=rep(0,p), index=rep(FALSE,p), lambda=lambda, residuals=y - mean(y), sigma=var(y), 
    #            iter=0, call=match.call(), options=list(post=post, intercept=intercept, normalize=normalize, control=control))
    est <- list(coefficients=rep(0,p), index=rep(FALSE,p), lambda=lambda, lambda0=lambda0, loadings=Ups0, residuals=y - mean(y), sigma=var(y), iter=0, call=match.call(), 
                options=list(post=post, intercept=intercept, normalize=normalize, control=control, mu=mu, meanx=meanx, scalex=normx))
    class(est) <- "lasso"
    return(est)
  }
  # refinement variance estimation
  mm <- 0
  if (post) {
    reg <- lm(y ~ -1 + x1)
    coefT <- coef(reg)
    coefT[is.na(coefT)] <- 0
    e1 <- y-x1%*%coefT
    coefTemp[ind1] <- coefT
    #e1 <- y-x1%*%coef(lm(y ~ -1 + x1))
  }
  if (!post) {
    e1 <- y-x1%*%coefTemp[ind1]
  }
  s1 <- sqrt(var(e1))
  s0 <- sqrt(var(y))
  
  while (mm<numIter) {
    if (dim(x1)[2]==0) {break}
    # X-independent
    if (control$lambda=="X-independent") {
      lambda <- rep(2*c*sqrt(n)*qnorm(1-gamma/(2*p))*s1,p)  
    }   
    # X-dependent
    if (control$lambda=="X-dependent") {
      lambda <- rep(c*quantile(sim, probs=1-gamma/2)*s1,p)
    }
    # standard
    if (control$lambda=="standard") {
      Ups1 <- 1/sqrt(n)*sqrt(t(t(e1^2)%*%(x^2)))
      lambda <- lambda0*Ups1
    }
    # none
    if (control$lambda=="none") {
      Ups1 <- 1/sqrt(n)*sqrt(t(t(e1^2)%*%(x^2)))
      lambda <- lambda0*Ups1
    }
    
    coefTemp <- LassoShooting.fit(x, y, lambda, XX=XX, Xy=Xy)$coefficients
    coefTemp[is.na(coefTemp)] <- 0
    ind1 <- (abs(coefTemp) > 0)
    x1 <- as.matrix(x[,ind1, drop=FALSE])
    
    if (dim(x1)[2]==0) { # added
      est <- list(coefficients=rep(0,p), index=rep(FALSE,p), lambda=lambda, lambda0=lambda0, loadings=Ups0, residuals=y - mean(y), sigma=var(y), iter=0, call=match.call(), 
                  options=list(post=post, intercept=intercept, normalize=normalize, control=control, mu=mu, meanx=meanx, scalex=normx))
      class(est) <- "lasso"
      return(est)
    }
    
    if (post) {
      reg <- lm(y ~ -1 + x1)
      coefT <- coef(reg)
      coefT[is.na(coefT)] <- 0
      e1 <- y-x1%*%coefT
      coefTemp[ind1] <- coefT
    }
    if (!post) {
      #e1 <- y-x1%*%coef(lm(y ~ -1 + x1))
      e1 <- y-x1%*%coefTemp[ind1]
      coefTemp <- coefTemp
    }
    s0 <- s1
    s1 <- sqrt(var(e1))
    mm <- mm+1 
    if(abs(s0-s1)<tol) {break}      
  }
  
  if (dim(x1)[2]==0) {
    coefTemp=NULL
    ind1 <- rep(0,p)
  }
  coefTemp <- scale(t(coefTemp), FALSE, normx)
  coefTemp <- as.vector(coefTemp)
  coefTemp[abs(coefTemp)<threshold] <- 0
  ind1 <- as.vector(ind1)
  names(coefTemp) <- names(ind1) <- colnames(x)
  est <- list(coefficients=coefTemp, index=ind1, lambda=lambda, lambda0=lambda0, loadings=Ups1, residuals=e1, sigma=s1, iter=mm, call=match.call(), options=list(post=post, intercept=intercept, normalize=normalize, control=control, mu=mu, meanx=meanx, scalex=normx))
  class(est) <- "lasso" 
  return(est)
}



lasso.formula <- function(formula, data, post=TRUE, intercept=TRUE, normalize=TRUE, control=list(c=1.1, gamma=.1, numIter=15, tol=10^-5, lambda="standard", numSim=10000, numfolds=10, lambda.start=NULL, threshold=NULL)) {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  attr(mt, "intercept") <- 0
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)
  
  est <- lasso(x,y, post=post, intercept=intercept, normalize=normalize, control=control)
  est$call <- cl
  return(est)
}


print.lasso <- function(object, all=TRUE ,digits = max(3L, getOption("digits") - 3L)) {
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (length(coef(object))) {
    if (all) {
    cat("Coefficients:\n")
    print.default(format(coef(object), digits = digits), print.gap = 2L, 
                  quote = FALSE)
    } else {
      print.default(format(coef(object)[object$index], digits = digits), print.gap = 2L, 
                    quote = FALSE)  
    }
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(object)
  }

summary.lasso <- function(object, all=TRUE, digits = max(3L, getOption("digits") - 3L)) {
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  cat("\nPost-Lasso Estimation: ",  paste(deparse(object$options$post), sep = "\n", collapse = "\n"), "\n", sep = " ")
  coefs <- object$coefficients
  p <- length(coefs)
  num.selected <- sum(abs(object$coefficients)>0) 
  cat("\nTotal number of variables:", p)
  cat("\nNumber of selected variables:", num.selected, "\n", sep=" ")
  resid <- object$residuals
  cat("\nResiduals: \n")
  nam <- c("Min", "1Q", "Median", "3Q", "Max")
  rq <- structure(apply(t(resid), 1L, quantile), dimnames = list(nam, dimnames(resid)[[2L]]))
  print(drop(t(rq)), digits = digits)
  cat("\n") 
  if (all) {
  coefm <- matrix(NA, length(coefs), 1)
  coefm[,1] <- coefs
  colnames(coefm) <- "Estimate"
  rownames(coefm) <- names(coefs)
  printCoefmat(coefm, digits = digits, na.print = "NA")
  } else {
  coefs <- coefs[abs(coefs)>0]
  coefm <- matrix(NA, length(coefs), 1)
  coefm[,1] <- coefs
  colnames(coefm) <- "Estimate"
  rownames(coefm) <- names(coefs)
  printCoefmat(coefm, digits = digits, na.print = "NA")
  }
  cat("\nResidual standard error:", format(signif(object$sigma, digits)))
  cat("\n")
  invisible(object)
}

model.matrix.lasso <- function(object) {
  if (is.call(object$call[[2]])) {
    X <- model.frame(object$call[[2]])
    mt <- attr(X, "terms")
    attr(mt, "intercept") <- 0
    mm <- model.matrix(mt, X)
  } else {
    mm <- eval(object$call[[2]])
  }
  return(mm)
}

# predict.lasso <- function (object, newdata=NULL) {
#   if (missing(newdata) || is.null(newdata)) {
#     X <- model.matrix(object)
#   }
#   else {
#     #mt <- attr(newdata, "terms")
#     #attr(mt, "intercept") <- 0
#     #m <- model.frame(mt, newdata)
#     #X <- model.matrix(mt, m)
#     X <- newdata
#   }
#   n <- length(object$residuals)
#   beta <- object$coefficients
#   if (object$options[["intercept"]]) {
#     yhat <- X%*%beta + object$options$mu - sum(object$options$meanx*beta) # (??)
#   }
#   if (!object$options[["intercept"]]) {
#     yhat <- X%*%beta
#   }
#   return(yhat) 
# }

predict.lasso <- function (object, newdata = NULL) 
{
  if (missing(newdata) || is.null(newdata)) {
    X <- model.matrix(object)
  }
  else {
    X <- newdata
  }
  n <- length(object$residuals)
  beta <- object$coefficients
  if (object$options[["intercept"]]) {
    if (is.null(object$options$mu))  object$options$mu <-0
    if (is.null(object$options$meanx))  object$options$meanx <-0
    yhat <- X %*% beta + object$options$mu - sum(object$options$meanx * 
                                                   beta)
  }
  if (!object$options[["intercept"]]) {
    yhat <- X %*% beta
  }
  return(yhat)
}

###################################### Cross Validation for penalty parameter lambda


cv.lasso <- function (x, y, K = 10, lambda.grid=NULL, post=TRUE, intercept=TRUE, normalize=TRUE)
{
  if (is.null(lambda.grid)) {
    lambda <- lasso(x, y, post, intercept, normalize)$lambda0
    lambda.grid <- seq(0, floor(lambda)*2, length.out=100)
  }
  n <- length(y)
  p <- dim(x)[2]
  folds <- K
  all.folds <- split(sample(1:n), rep(1:folds, length = n))
  residmat <- matrix(NA, length(lambda.grid), K)
  for (j in 1:length(lambda.grid)) {
    for (i in seq(K)) {
      omit <- all.folds[[i]]
      fit <- lasso(x[-omit, , drop = FALSE], y[-omit], post=post, intercept=intercept, normalize=normalize,
                   control=list(numIter=15, tol=10^-5, lambda="none", lambda.start=rep(lambda.grid[j],p)))
      fitvalues <- predict(fit, x[omit, , drop = FALSE])
      if (length(omit) == 1) 
        fitvalues <- matrix(fitvalues, nrow = 1)
      residmat[j, i] <- apply((y[omit] - fitvalues)^2, 2, mean)
    }
  if (sum(fit$index)==0) break
  }
  cv <- apply(residmat, 1, mean)
  cv.error <- sqrt(apply(residmat, 1, var)/K)
  ind <- which(cv==min(cv, na.rm=T))
  lambda.cv <- lambda.grid[ind]
  object <- list(lambda.grid = lambda.grid, cv = cv, cv.error = cv.error, lambda.cv=lambda.cv)
  return(object)
} 