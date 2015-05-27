logistic.lasso <- function(x, y, post=TRUE, intercept=TRUE, normalize=TRUE, 
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
  a0 <- mean(y) # intercept
  # Intercept handling and scaling
  
  #if (normalize) {
  #  normx <- sqrt(apply(x,2,var))
  #  ind <- which(normx < eps)
  #  if (length(ind)!=0) {
  #    x <- x[,-ind]
  #    normx <- normx[-ind]
  #    ind.names <- ind.names[-ind]
  #    p <- dim(x)[2]
  #    lambda.start <- lambda.start[-ind]
  #  }
  #  x <- scale(x, FALSE, normx)
  #} else {
  #  normx <- rep(1,p)
  #}
  
  #if (intercept) {
  #  x <- cbind(1,x)
  #  colnames(x)[1] <- "intercept"
  #}
  
  #XX <- crossprod(x)
  #Xy <- crossprod(x,y) 
  
  
  # calculation lambda
  e0 <- y - mean(y)
  Ups0 <- 1/sqrt(n)*as.vector(sqrt(t(t(e0^2)%*%(x^2))))
  # X-independent
  if (control$lambda=="X-independent") {
    lambda0 <- 2*c*sqrt(n)*qnorm(1-gamma/(2*p))
    
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
  }
  # "standard"
  if (control$lambda=="standard") { # choice paper
    lambda0 <- 2*c*sqrt(n)*sqrt(2*log(2*p*log(n)/gamma))   
  }
  
  if (!is.null(lambda.start)) {
    #lambda0 <- lambda.start
    lambda0 <-diag(p)%*%as.matrix(lambda.start)
  }
  
  if (control$lambda=="none") {
    lambda0 <- lambda.start
  }
  lambda <- lambda0*Ups0 
  lambda0 <- lambda0/(2*n)
  
  # calculation first parameters
  #log.lasso <- glmnet(x, y, family=c("binomial"), alpha = 1, lambda=0.5*lambda0[1],
  #                    standardize = normalize, intercept=intercept, thresh = tol, penalty.factor=Ups0)
  log.lasso <- glmnet(t(t(x)/Ups0), y, family=c("binomial"), alpha = 1, lambda=0.5*lambda0[1],
                       standardize = normalize, intercept=intercept, thresh = tol)
  log.lasso$beta <- log.lasso$beta/Ups0
  coefTemp <- as.vector(log.lasso$beta)
  coefTemp[is.na(coefTemp)] <- 0
  ind1 <- (abs(coefTemp) > 0)
  x1 <- as.matrix(x[,ind1, drop=FALSE])
  if (dim(x1)[2]==0) {
    print("No variables selected!")
    est <- list(coefficients=rep(0,p), index=rep(FALSE,p), a0=a0, lambda=lambda, lambda0=lambda0, loadings=Ups0, residuals=y - mean(y), sigma=var(y), iter=0, call=match.call(), 
                options=list(post=post, intercept=intercept, normalize=normalize, control=control))
    class(est) <- c("logistic.lasso")
    return(est)
  }
  # refinement variance estimation
  mm <- 0
  if (post) {
    #reg <- glm.fit(x1,y, family = binomial(link = "logit"), intercept=intercept)
    reg <- glm(y~x1, family = binomial(link = "logit"))
    coefT <- coef(reg)[-1]
    coefT[is.na(coefT)] <- 0
    e1 <- y-reg$fitted.values
    coefTemp[ind1] <- coefT
  }
  if (!post) {
    e1 <- y- predict(log.lasso, newx=x ,type="response")
  }
  s1 <- sqrt(var(e1))
  s0 <- sqrt(var(y))
  
  while (mm<numIter) {
    if (dim(x1)[2]==0) {break}
    Ups1 <- 1/sqrt(n)*as.vector(sqrt(t(t(e1^2)%*%(x^2))))
    #log.lasso <- glmnet(x, y, family=c("binomial"), alpha = 1, lambda=0.5*lambda0[1],
    #                   standardize = normalize, intercept=intercept, thresh = tol, penalty.factor=Ups1)
    log.lasso <- glmnet(t(t(x)/Ups1), y, family=c("binomial"), alpha = 1, lambda=lambda0[1],
                        standardize = normalize, intercept=intercept, thresh = tol)
    log.lasso$beta <- log.lasso$beta/Ups1
    coefTemp <- as.vector(log.lasso$beta)
    coefTemp[is.na(coefTemp)] <- 0
    ind1 <- (abs(coefTemp) > 0)
    x1 <- as.matrix(x[,ind1, drop=FALSE])
    if (dim(x1)[2]==0) {
      print("No variables selected!")
      est <- list(coefficients=rep(0,p), index=rep(FALSE,p), a0=a0, lambda=lambda, lambda0=lambda0, loadings=Ups0, residuals=y - mean(y), sigma=var(y), iter=0, call=match.call(), 
                  options=list(post=post, intercept=intercept, normalize=normalize, control=control))
      class(est) <- c("logistic.lasso")
      return(est)
    }
    
    if (post) {
      #reg <- glm.fit(x1,y, family = binomial(link = "logit"), intercept=intercept)
      #coefT <- coef(reg)
      #coefT[is.na(coefT)] <- 0
      #e1 <- y-reg$fitted.values
      #coefTemp[ind1] <- coefT
      reg <- glm(y~x1, family = binomial(link = "logit"))
      coefT <- coef(reg)[-1]
      coefT[is.na(coefT)] <- 0
      e1 <- y-reg$fitted.values
      coefTemp[ind1] <- coefT
    }
    if (!post) {
      e1 <- y- predict(log.lasso, newx=x, type="response")
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
  
  if (intercept==TRUE) {
    if (post==TRUE) a0 <- coef(reg)[1]
    if (post==FALSE) a0 <- coef(glmnet)[1]
  }
  coefTemp <- as.vector(coefTemp)
  coefTemp[abs(coefTemp)<threshold] <- 0
  ind1 <- as.vector(ind1)
  names(coefTemp) <- names(ind1) <- colnames(x)
  est <- list(coefficients=coefTemp, a0=a0, index=ind1, lambda=lambda, lambda0=lambda0, loadings=Ups1, residuals=e1, sigma=s1, iter=mm, call=match.call(), options=list(post=post, intercept=intercept, normalize=normalize, control=control))
  class(est) <- c("logistic.lasso") 
  return(est)
}

###############################################################################################################################

################################################################################################

predict.logistic.lasso <- function (object, newdata=NULL) {
  if (missing(newdata) || is.null(newdata)) {
    X <- model.matrix(object)
  }
  else {
    #mt <- attr(newdata, "terms")
    #attr(mt, "intercept") <- 0
    #m <- model.frame(mt, newdata)
    #X <- model.matrix(mt, m)
    X <- newdata
  }
  n <- length(object$residuals)
  beta <- object$coefficients
  if (object$options[["intercept"]]) {
    yp <- object$a0 + X%*%as.vector(beta)
    yhat <- 1/(1+exp(-yp))
  }
  if (!object$options[["intercept"]]) {
    yp <- X%*%as.vector(beta)
    yhat <- 1/(1+exp(-yp))
  }
  return(yhat) 
}

