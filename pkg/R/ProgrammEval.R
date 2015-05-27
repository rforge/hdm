 # glm.lasso: type="logit", "normal"
# to be done: logit, predict logistic

glm.lasso <- function(x, y, post=TRUE, intercept=TRUE, normalize=TRUE, type="normal",control=list(c=1.1, gamma=.1, numIter=15, tol=10^-5, lambda="standard", numSim=10000, 
                                   numfolds=10, lambda.start=NULL, threshold=NULL)) {
  if (type=="normal") {
  est <- lasso(x, y, post=post, intercept=intercept, normalize=normalize, control=control)
  est$type <- type
  }
  
  if (type=="logit") {
    est <- logistic.lasso(x, y,  post=post, intercept=intercept, normalize=normalize, control=control)
                 est$type <- type
  }
  
  return(est)
}

############
###############################################################################################
# Definition

# my_dx = E[Y | D, X];
# md_x = E[D | X];
# mz_x = E[Z | X];
# md_zx = E[D | Z, X];
# my_zx = E[Y | Z, X];
# md = E[D];
# mz = E[Z];
# myd1_zx = E[YD | Z, X];
# myd0_zx = E[Y(1-D) | Z, X];

#######################
# LATE
#######################
late <- function(y,d,z,my_z1x,my_z0x,mz_x,md_z1x,md_z0x, bootstrap=NULL, nRep=500) {
  n <- length(as.vector(y))
  object <- NULL
  eff <- (z*(y-my_z1x)/mz_x -  ((1-z)*(y-my_z0x)/(1-mz_x)) + my_z1x - my_z0x )/
          mean(z*(d-md_z1x)/mz_x - ((1-z)*(d-md_z0x)/(1-mz_x)) + md_z1x - md_z0x)
  
  #eff <- ((y[z==1]-my_z1x)/mz_x[z==1] -  (y[z==0]-my_z0x)/(1-mz_x[z==1]) + my_z1x - my_z0x )/
  #        mean(z*(d-md_z1x)/mz_x - ((1-z)*(d-md_z0x)/(1-mz_x)) + md_z1x - md_z0x)
  
  object$se <- sqrt(var(eff))/sqrt(n)
  object$late <- mean(eff)
  object$individual <- eff
  
  if (!is.null(bootstrap)) {
    boot <- rep(NA,nRep)
    for (i in 1:nRep) {
      if (bootstrap=="Bayes") {
        weights <- rexp(n, rate=1) -1
      }
      if (bootstrap=="normal") {
        weights <- rnorm(n)
      }
      if (bootstrap=="wild") {
        weights <- rnorm(n)/sqrt(2) + (rnorm(n)^2-1)/2
      }
      weights <- weights + 1
      boot[i] <- mean(weights*(z*(y-my_z1x)/mz_x -  ((1-z)*(y-my_z0x)/(1-mz_x)) + my_z1x - my_z0x))/
  mean(weights*(z*(d-md_z1x)/mz_x - ((1-z)*(d-md_z0x)/(1-mz_x)) + md_z1x - md_z0x))
  }
  object$boot.se <- sqrt(var(boot))
  }

  return(object)
}

#######################
# LATT
#######################
latt <- function(y,d,z,my_z0x,mz_x,md_z0x, bootstrap=NULL, nRep=500) {
  n <- length(as.vector(y))
  object <- NULL
  
  eff <- ((y-my_z0x) - (1-z)*(y-my_z0x)/(1-mz_x))/mean((d-md_z0x) - (1-z)*(d-md_z0x)/(1-mz_x))

  object$se <- sqrt(var(eff))/sqrt(n)
  object$latt <- mean(eff)
  object$individual <- eff
  
  if (!is.null(bootstrap)) {
    boot <- rep(NA,nRep)
    for (i in 1:nRep) {
      if (bootstrap=="Bayes") {
        weights <- rexp(n, rate=1) -1
      }
      if (bootstrap=="normal") {
        weights <- rnorm(n)
      }
      if (bootstrap=="wild") {
        weights <- rnorm(n)/sqrt(2) + (rnorm(n)^2-1)/2
      }
      weights <- weights +1
      #boot[i] <- mean(weights*(z*(y-my_z1x)/mz_x -  ((1-z)*(y-my_z0x)/(1-mz_x)) + my_z1x - my_z0x))/
      #  mean(weights*(z*(d-md_z1x)/mz_x - ((1-z)*(d-md_z0x)/(1-mz_x)) + md_z1x - md_z0x))
      boot[i] <- mean(weights*((y-my_z0x) - (1-z)*(y-my_z0x)/(1-mz_x)))/mean(weights*((d-md_z0x) - (1-z)*(d-md_z0x)/(1-mz_x)))
    }
    object$boot.se <- sqrt(var(boot))
  }
  
  return(object)
}

#######################
# LASF
#######################
lasf<- function(y,d,z,myd1_z1x,myd1_z0x,myd0_z1x,myd0_z0x,mz_x,md_z1x,md_z0x,d0, weights=NULL, bootstrap=NULL, nRep=500) {
  n <- length(as.vector(y))
  object <- NULL
  if (is.null(weights)) weights <- rep(1,n)
  eff <-  weights*(d0*(z*(d*y-myd1_z1x)/mz_x - ((1-z)*(d*y-myd1_z0x)/(1-mz_x)) + myd1_z1x - myd1_z0x) + #+/- 
          (d0-1)*(z*((1-d)*y-myd0_z1x)/mz_x - ((1-z)*((1-d)*y-myd0_z0x)/(1-mz_x)) + myd0_z1x - myd0_z0x))/
          mean(weights*(z*(d-md_z1x)/mz_x - ((1-z)*(d-md_z0x)/(1-mz_x)) + md_z1x - md_z0x))
  object$se <- sqrt(var(eff))/sqrt(n)
  object$lasf <- mean(eff)
  object$individual <- eff
  
  if (!is.null(bootstrap)) {
    boot <- rep(NA,nRep)
    for (i in 1:nRep) {
      if (bootstrap=="Bayes") {
        weights <- rexp(n, rate=1) -1
      }
      if (bootstrap=="normal") {
        weights <- rnorm(n)
      }
      if (bootstrap=="wild") {
        weights <- rnorm(n)/sqrt(2) + (rnorm(n)^2-1)/2
      }
      weights <- weights +1
      boot[i] <-mean(weights*(d0*(z*(d*y-myd1_z1x)/mz_x - ((1-z)*(d*y-myd1_z0x)/(1-mz_x)) + myd1_z1x - myd1_z0x) + 
                          (d0-1)*(z*((1-d)*y-myd0_z1x)/mz_x - ((1-z)*((1-d)*y-myd0_z0x)/(1-mz_x)) + myd0_z1x - myd0_z0x)))/mean(weights*(z*(d-md_z1x)/mz_x - ((1-z)*(d-md_z0x)/(1-mz_x)) + md_z1x - md_z0x))
    }
    object$boot.se <- sqrt(var(boot))
  }
  
  return(object)
}

#######################
# LASF-T
#######################

lasft <- function(y,d,z,myd1_z0x,myd0_z0x,mz_x,md_z0x,d0, weights=NULL, bootstrap=NULL, nRep=500) {
  n <- length(as.vector(y))
  object <- NULL
  if (is.null(weights)) weights <- rep(1,n)
  
  eff <- weights*(d0*((y*d-myd1_z0x) - (1-z)*(d*y-myd1_z0x)/(1 - mz_x)) + 
                    (d0-1)*((y*(1 - d)-myd0_z0x) - (1-z)*((1-d)*y - myd0_z0x)/(1-mz_x)))/mean(weights*(d-md_z0x) - (1-z)*(d-md_z0x)/(1-mz_x))
  object$se <- sqrt(var(eff))/sqrt(n)
  object$lasft <- mean(eff)
  object$individual <- eff
  
  if (!is.null(bootstrap)) {
    boot <- rep(NA,nRep)
    for (i in 1:nRep) {
      if (bootstrap=="Bayes") {
        weights <- rexp(n, rate=1) -1
      }
      if (bootstrap=="normal") {
        weights <- rnorm(n)
      }
      if (bootstrap=="wild") {
        weights <- rnorm(n)/sqrt(2) + (rnorm(n)^2-1)/2
      }
      weights <- weights +1
      boot[i] <- mean(weights*(d0*((y*d-myd1_z0x) - (1-z)*(d*y-myd1_z0x)/(1 - mz_x)) 
                          + (d0-1)*((y*(1 - d)-myd0_z0x) - (1-z)*((1-d)*y - myd0_z0x)/(1-mz_x))))/
      mean(weights*((d-md_z0x) - (1-z)*(d-md_z0x)/(1-mz_x)))
    }
    object$boot.se <- sqrt(var(boot))
  }
  
  return(object)
}


#############################
# LQTE and LQTT
#############################


#function [lqte,talpha_lqte,lqtt,talpha_lqtt] = lqteAndlqtt(yF,d,z,taus,percs,myd1_z1x,myd0_z1x,myd1_z0x,myd0_z0x,mz_x,md_z1x,md_z0x,alpha)
localQuantile <- function(yF,d,z,taus,percs,myd1_z1x,myd0_z1x,myd1_z0x,myd0_z0x,mz_x,md_z1x,md_z0x, bootstrap=NULL, nRep=500, alpha=0.05) {

  n <- dim(yF)[1]
  k <- dim(yF)[2]
  object <- NULL
  lasf0 <- matrix(NA, ncol=2, nrow=k)
  lasf1 <- matrix(NA, ncol=2, nrow=k)
  lasft0 <- matrix(NA, ncol=2, nrow=k)
  lasft1 <- matrix(NA, ncol=2, nrow=k)
  for (i in 1:k) {
    temp1 <- lasf(yF[,i],d,z,myd1_z1x[,i], myd1_z0x[,i],myd0_z1x[,i],myd0_z0x[,i],mz_x,md_z1x,md_z0x,0)
    lasf0[i,] <- c(temp1$lasf, temp1$se)
    temp1 <- lasf(yF[,i],d,z,myd1_z1x[,i],myd1_z0x[,i],myd0_z1x[,i],myd0_z0x[,i],mz_x,md_z1x,md_z0x,1)
    lasf1[i,] <-  c(temp1$lasf, temp1$se)
    temp1 <- lasft(yF[,i],d,z,myd1_z0x[,i],myd0_z0x[,i],mz_x,md_z0x,0)
    lasft0[i,] <-  c(temp1$lasft, temp1$se)
    temp1 <- lasft(yF[,i],d,z,myd1_z0x[,i],myd0_z0x[,i],mz_x,md_z0x,1)
    lasft1[i,] <-  c(temp1$lasft, temp1$se)
  }
  if (sum(is.na(lasf0[,1])) == 0) {
    lqsf0 <- approx(lasf0[,1],percs, xout=taus)$y
  } else {
    lqsf0 <- NA*rep(1,length(taus))
  }
  if (sum(is.nan(lasf1[,1])) == 0) {
  lqsf1 <- approx(lasf1[,1],percs,xout=taus)$y
  } else {
    lqsf1 <- NA*rep(1,length(taus))
  }
  if (sum(is.nan(lasft0[,1]))==0) { 
  lqsft0 <- approx(lasft0[,1],percs,xout=taus)$y
  } else {
  lqsft0 <- NA*rep(1,length(taus))
  }
  if (sum(is.nan(lasft1[,1])) == 0) {
  lqsft1 <- approx(lasft1[,1],percs,xout=taus)$y
  } else {
  lqsft1 <- NA*rep(1,length(taus))
  }

  lqte <- lqsf1 - lqsf0
  lqtt <- lqsft1 - lqsft0
  bs_std_lqte <-  bs_std_lqtt <- NULL
  talpha_lqte <- talpha_lqtt <- NULL
    if (!is.null(bootstrap)) {
      
      bs_lqte <- matrix(NA, ncol=k, nrow=nRep)
      bs_lqtt <- matrix(NA, ncol=k, nrow=nRep)
      
      for (i in 1:nRep) {
        if (bootstrap=="Bayes") {
          weights <- rexp(n, rate=1) -1
        }
        if (bootstrap=="normal") {
          weights <- rnorm(n)
        }
        if (bootstrap=="wild") {
          weights <- rnorm(n)/sqrt(2) + (rnorm(n)^2-1)/2
        }
        weights <- weights + 1

        bs_lasf0 <- matrix(NA, ncol=2,nrow=k)
        bs_lasf1 <- matrix(NA, ncol=2,nrow=k)
        bs_lasft0 <- matrix(NA, ncol=2,nrow=k)
        bs_lasft1 <- matrix(NA, ncol=2,nrow=k)
          for (j in 1:k) {
            temp1 <- lasf(yF[,j],d,z,myd1_z1x[,j], myd1_z0x[,j],myd0_z1x[,j],myd0_z0x[,j],mz_x,md_z1x,md_z0x,0, weights=weights)
            bs_lasf0[j,] <- c(temp1$lasf, temp1$se)
            temp1 <- lasf(yF[,j],d,z,myd1_z1x[,j],myd1_z0x[,j],myd0_z1x[,j],myd0_z0x[,j],mz_x,md_z1x,md_z0x,1, weights=weights)
            bs_lasf1[j,] <- c(temp1$lasf, temp1$se)
            temp1 <- lasft(yF[,j],d,z,myd1_z0x[,j],myd0_z0x[,j],mz_x,md_z0x,0, weights=weights)
            bs_lasft0[j,] <- c(temp1$lasft, temp1$se)
            temp1 <- lasft(yF[,j],d,z,myd1_z0x[,j],myd0_z0x[,j],mz_x,md_z0x,1, weights=weights)
            bs_lasft1[j,] <- c(temp1$lasft, temp1$se)
          }
        if (sum(is.na(bs_lasf0[,1])) == 0) {
          bs_lqsf0 <- approx(bs_lasf0[,1],percs, xout=taus)$y
        } else {
          bs_lqsf0 <- NA*rep(1,length(taus))
        }
        if (sum(is.nan(bs_lasf1[,1])) == 0) {
          bs_lqsf1 <- approx(bs_lasf1[,1],percs,xout=taus)$y
        } else {
          bs_lqsf1 <- NA*rep(1,length(taus))
        }
        if (sum(is.nan(bs_lasft0[,1]))==0) { 
          bs_lqsft0 <- approx(bs_lasft0[,1],percs,xout=taus)$y
        } else {
          bs_lqsft0 <- NA*rep(1,length(taus))
        }
        if (sum(is.nan(bs_lasft1[,1])) == 0) {
          bs_lqsft1 <- approx(bs_lasft1[,1],percs,taus)$y
        } else {
          bs_lqsft1 <- NA*rep(1,length(taus))
        }

        bs_lqte[i,] <- bs_lqsf1 - bs_lqsf0
        bs_lqtt[i,] <- bs_lqsft1 - bs_lqsft0

      }
    bs_std_lqte <- sqrt(apply(bs_lqte,2,var))
    bs_std_lqtt <- sqrt(apply(bs_lqtt,2,var))

    bs_t_lqte <- abs(bs_lqte - matrix(rep(lqte, nRep), byrow=T, ncol=k)) / matrix(rep(bs_std_lqte, nRep), byrow=T, ncol=k)
    bs_tmax_lqte <- apply(bs_t_lqte, 2, max)
    talpha_lqte <- quantile(bs_tmax_lqte,probs=1-alpha, na.rm=TRUE)

    bs_t_lqtt <-  abs(bs_lqtt - matrix(rep(lqtt, nRep), byrow=T, ncol=k)) / matrix(rep(bs_std_lqtt, nRep), byrow=T, ncol=k) 
    bs_tmax_lqtt <- apply(bs_t_lqtt, 2,max)
    talpha_lqtt <- quantile(bs_tmax_lqtt, probs=1-alpha, na.rm=TRUE)
    }
    object$lqte <- list(lqte=lqte, se=bs_std_lqte, t=talpha_lqte)
    object$lqtt <- list(lqtt=lqtt, se=bs_std_lqtt, t=talpha_lqtt)
  return(object)
}

################################################################################################################################################
#################################### Main functioin for Programm Evaluation
################################################################################################################################################
ProgEval <- function(x,y,z,d, post=TRUE, intercept=TRUE, normalize=TRUE,  tau = (5:95)/100,  alpha=0.05, bootstrap=NULL, nRep=500) {
  normalize <- FALSE
  post <- TRUE
  n <- dim(x)[1]
  p <- dim(x)[2]
  object <- NULL
  
  # LATE and LATT
  lambda <- 2.2*sqrt(n)*qnorm(1-(1/log(n))/(2*(2*p)))
  control <- list(c = 1.1, gamma = 0.1, numIter = 15, tol = 10^-5, lambda = "none",   lambda.start = rep(lambda, p))
  indz1 <- (z==1)
  indz0 <- (z==0)
  # E[Y|Z = 1,X] = my_z1x
  b_y_z1xL <- lasso(x[indz1,,drop=FALSE], y[indz1], post=post, intercept=intercept, normalize=normalize, control=control)
  my_z1x <- predict(b_y_z1xL, newdata=x)
  # E[Y|Z = 0,X] = my_z0x
  b_y_z0xL <- lasso(x[indz0,], y[indz0],  post=post, intercept=intercept, normalize=normalize, control=control)
  my_z0x <- predict(b_y_z0xL, newdata=x)
  # E[D|Z = 1,X] = md_z1x
  lambda <- 2.2*sqrt(n)*qnorm(1-(1/log(n))/(2*(2*p)))
  control <- list(c = 1.1, gamma = 0.1, numIter = 15, tol = 10^-5, lambda = "none",   lambda.start = rep(lambda, p))
  b_d_z1xL <- logistic.lasso(x[indz1,], d[indz1],  post=post, intercept=intercept, normalize=normalize, control=control)
  md_z1x <- predict(b_d_z1xL, newdata=x)
  # E[D|Z = 0,X] = md_z0x
  md_z0x <- rep(0,n)
  
  # E[Z|X] = mz_x
  lambdaP <- 2.2*sqrt(n)*qnorm(1-(1/log(n))/(2*p))
  control <- list(c = 1.1, gamma = 0.1, numIter = 15, tol = 10^-5, lambda = "none",   lambda.start = rep(lambdaP, p))
  b_z_xL <- logistic.lasso(x, z, post=post, intercept=intercept, normalize=normalize, control=control)
  mz_x <- predict(b_z_xL, newdata=x)
  mz_x <- mz_x*(mz_x > 1e-12 & mz_x < 1-1e-12) + (1-1e-12)*(mz_x > 1-1e-12) + 1e-12*(mz_x < 1e-12)
  
  # LATE 
  object$LATE <- late(y,d,z,my_z1x,my_z0x,mz_x,md_z1x,md_z0x, bootstrap=bootstrap, nRep=nRep)
  # LATT 
  object$LATT <- latt(y,d,z,my_z0x,mz_x,md_z0x, bootstrap=bootstrap, nRep=nRep)
  
  ## LQTE and LQTT 
  
  # Create matrix with columns of form 1(y < u) for u approximately covering
  # the support of y
  
  pt <- quantile(y,tau)
  percs <- Unique(pt)$data
  Itaus <- Unique(pt)$index.unique
  Ipercs <-  Unique(pt)$index.duplicate
  taus <- tau[Itaus]
  # calU <- (taus >= .1) & (taus <= .9)  # Will look at QTE over this interval
  yF <- as.matrix(y)%*%rep(1,length(taus)) <= as.matrix(rep(1,n))%*%percs
  

  ##############################
  # LQTE and LQTT
  
  # Need to do LOTS of selection.  Need myd1_z1x = E[YD | Z = 1, X], 
  # myd0_z1x = E[Y(1-D) | Z = 1, X], myd1_z0x = E[YD | Z = 0, X], 
  # myd0_z0x = E[Y(1-D) | Z = 0, X] for every column in yF
  
  # Note that E[YD|Z = 0,X] = 0
  
  # Note that E[D|Z = 1,X], E[D|Z = 0,X], and E[Z|X] are as above, but 
  # adjusting penalty level to account for functional data.
  lambdaU <- 2.2*sqrt(n)*qnorm(1-(1/log(n))/(2*(2*p)*n))
  control <- list(c = 1.1, gamma = 0.1, numIter = 15, tol = 10^-5, lambda = "none",   lambda.start = rep(lambdaU, p))
  # E[D|Z = 1,X] = md_z1x
  b_d_z1xLU <- logistic.lasso(x[indz1,], d[indz1], post=post, intercept=intercept, normalize=normalize, control=control)
  md_z1xU <- predict(b_d_z1xLU, newdata=x)
  # E[D|Z = 0,X] = md_z0x
  md_z0xU <- rep(0,n)
  # E[Z|X] = mz_x
  lambdaUP <- 2.2*sqrt(n)*qnorm(1-(1/log(n))/(2*p*n))
  control <- list(c = 1.1, gamma = 0.1, numIter = 15, tol = 10^-5, lambda = "none",   lambda.start = rep(lambdaUP, p))
  b_z_xLU <- logistic.lasso(x,z,post=post, intercept=intercept, normalize=normalize, control=control)
  mz_xU  <- predict(b_z_xLU, x)
  mz_xU <- mz_xU*(mz_xU > 1e-12 & mz_xU < 1-1e-12) + (1-1e-12)*(mz_xU > 1-1e-12) + 1e-12*(mz_xU < 1e-12)
  
  # Variable selection using feasible Lasso
  k <- dim(yF)[2]
  b_yd1_z1xL <- matrix(NA, ncol=k, nrow=p)
  b_yd0_z1xL <- matrix(NA, ncol=k, nrow=p)
  b_yd0_z0xL <- matrix(NA, ncol=k, nrow=p)
  use_yd1_z1x <- matrix(FALSE, ncol=k, nrow=p)
  use_yd0_z1x <- matrix(FALSE, ncol=k, nrow=p)
  use_yd0_z0x <- matrix(FALSE, ncol=k, nrow=p)
  myd1_z1xU <- matrix(FALSE, ncol=k, nrow=dim(yF)[1])
  myd0_z1xU <- matrix(FALSE, ncol=k, nrow=dim(yF)[1])
  myd1_z0xU <- matrix(FALSE, ncol=k, nrow=dim(yF)[1])
  myd0_z0xU <- matrix(FALSE, ncol=k, nrow=dim(yF)[1])
  
  for (ii in 1:k) {
  lambdaU <- 2.2*sqrt(n)*qnorm(1-(1/log(n))/(2*(2*p)*n))
  control <- list(c = 1.1, gamma = 0.1, numIter = 15, tol = 10^-5, lambda = "none",   lambda.start = rep(lambdaU, p))
  # Selection for E[YD|Z=1,X]
  yd1_z1xL <- logistic.lasso(x[indz1,], yF[indz1,ii]*d[indz1], post=post, intercept=intercept, normalize=normalize, control=control)
  b_yd1_z1xL[,ii] <- as.vector(yd1_z1xL$coef)
  myd1_z1xU[,ii] <- predict(yd1_z1xL, newdata=x)
  # Selection for E[Y(1-D)|Z=1,X]
  yd0_z1xL <- logistic.lasso(x[indz1,], yF[indz1,ii]*(1-d[indz1]), post=post, intercept=intercept, normalize=normalize, control=control)
  b_yd0_z1xL[,ii] <- yd0_z1xL$coef
  myd0_z1xU[,ii] <- predict(yd0_z1xL, newdata=x) 
  # Selection for E[YD|Z=0,X]
  myd1_z0xU[,ii] <- rep(0,n)
  # Selection for E[Y(1-D)|Z=0,X]
  yd0_z0xL <- logistic.lasso(x[indz0,], yF[indz0,ii]*(1-d[indz0]),post=post, intercept=intercept, normalize=normalize, control=control)
  b_yd0_z0xL[,ii] <- yd0_z0xL$coef
  myd0_z0xU[,ii] <- predict(yd0_z0xL, newdata=x)
  }
  quant <- localQuantile(yF,d,z,taus,percs,myd1_z1xU,myd0_z1xU,myd1_z0xU,myd0_z0xU,mz_xU,md_z1xU,md_z0xU, bootstrap=bootstrap, nRep=nRep, alpha=alpha)
  object$LQTE <- quant$lqte
  object$LQTT <- quant$lqtt
  object$bootstrap <- bootstrap
  object$tau <- Itaus
  return(object)
}

############################################################################################
# help functions
############################################################################################
Unique <- function(y) {
  object <- NULL
  object$data <- unique(y)
  object$index.unique <- seq_along(y)[!duplicated(y)]
  object$index.duplicated <- split(seq_along(y), y)
  return(object)
}

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

############### plot quantile

LQplot <- function(object, title="", ylab="") {
#   #LQTE
#   dat <- data.frame(lqte=object$LQTE$lqte, quantile= object$tau/100)
#   dat$ucl <- object$LQTE$lqte +  object$LQTE$t * object$LQTE$se
#   dat$lcl <- object$LQTE$lqte -  object$LQTE$t * object$LQTE$se
#   dat <- dat[complete.cases(dat),]
#   lqte_plot <- qplot(quantile, lqte, data=dat) + 
#     geom_smooth(aes(ymin = lcl, ymax = ucl), data=dat, stat="identity") +
#     ylab(ylab) + 
#     xlab("Quantile") +
#     ggtitle(paste(title, "(QLTE)"))
#   # LQTT
#   dat <- data.frame(lqtt=object$LQTT$lqtt, quantile= object$tau/100)
#   dat$ucl <- object$LQTT$lqtt +  object$LQTT$t * object$LQTT$se
#   dat$lcl <- object$LQTT$lqtt -  object$LQTT$t * object$LQTT$se
#   dat <- dat[complete.cases(dat),]
#   lqtt_plot <- qplot(quantile, lqtt, data=dat) + 
#     geom_smooth(aes(ymin = lcl, ymax = ucl), data=dat, stat="identity") +
#     ylab(ylab) + 
#     xlab("Quantile") +
#     ggtitle(paste(title, "(QLTE-T)"))
#   # Define grid layout to locate plots and print each graph
#   pushViewport(viewport(layout = grid.layout(1, 2)))
#   print(lqte_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
#   print(lqtt_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  
  ###
  k <- length(object$tau)
  effect <- ucl <- lcl <- NULL
  dat1 <- data.frame(effect= object$LQTE$lqte, quantile= object$tau/100, treatment=rep("LQTE",k))
  dat1$ucl <- object$LQTE$lqte +  object$LQTE$t * object$LQTE$se
  dat1$lcl <- object$LQTE$lqte -  object$LQTE$t * object$LQTE$se
  
  dat2 <- data.frame(effect=object$LQTT$lqtt, quantile= object$tau/100, treatment=rep("LQTT",k))
  dat2$ucl <- object$LQTT$lqtt +  object$LQTT$t * object$LQTT$se
  dat2$lcl <- object$LQTT$lqtt -  object$LQTT$t * object$LQTT$se
  dat <- rbind(dat1, dat2)
  
  return(
  qplot(quantile, effect, data=dat) + 
    geom_smooth(aes(ymin = lcl, ymax = ucl), data=dat, stat="identity") + facet_grid(. ~ treatment) +
    ylab(ylab) + 
    xlab("Quantile") +
    ggtitle(title)
  )
}

#LQplot(res, title="Hello World")
