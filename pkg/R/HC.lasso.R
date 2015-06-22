HC.lasso <- function(x, y, d, I3=NULL,  ...) {
  d <- as.matrix(d, ncol=1)
  y <- as.matrix(y, ncol=1)
  kx <- dim(x)[2]
  if (is.null(colnames(d))) colnames(d) <- "d1"
  if (is.null(colnames(x)) & !is.null(x)) colnames(x) <- paste("x", 1:kx, sep="")
  I1 <- lasso(x, d, ...)$index
  I2 <- lasso(x, y, ...)$index
  n <- dim(x)[1]
  
  if (is.logical(I3)) {
    I <- I1+I2+I3
    I <- as.logical(I)
  } else {
    I <- I1+I2
    I <- as.logical(I)
  }
  if (sum(I)==0) {I <- NULL}
  x <- cbind(d,x[,I, drop=FALSE])
  reg1 <- lm(y~x)
  alpha <- coef(reg1)[2]
  xi<-reg1$residuals*sqrt(n/(n-sum(I)-1))
  if(is.null(I)) {reg2 <- lm(d~1)}
  if(!is.null(I)) {reg2 <- lm(d~x[,-1, drop=FALSE])}
  v <- reg2$residuals
  var <- 1/n*1/mean(v^2)*mean(v^2*xi^2)*1/mean(v^2)
  se <- sqrt(var)
  tval <- alpha/sqrt(var)
  pval <- 2*pnorm(-abs(tval))
  if (is.null(I)) {
    no.selected <- 1
  } else {
    no.selected <- 0
  }
  res <- list(epsilon = xi, v=v)
  return(list(alpha=unname(alpha), se=drop(se), t=unname(tval), pval=unname(pval), no.selected=no.selected, coefficients=coef(reg1), residuals=res))
}



HC.lasso.wrap <- function(x, y, index=c(1:ncol(x)), I3=NULL, ...) {
  k <- length(index)
  n <- dim(x)[1]
  if (is.null(colnames(x))) colnames(x) <- paste("V", 1:dim(x)[2], sep="")
  table <- matrix(NA,ncol=4,nrow=k)
  rownames(table) <- colnames(x)[index]
  colnames(table) <- c("coeff.", "se.", "t-value", "p-value")
  for (i in 1:k)
  {
    d <- x[,index[i], drop=FALSE]
    Xt <- x[,-index[i], drop=FALSE]
    
    result <- try(col <- HC.lasso(Xt,y,d, I3=I3, ...))
    if(class(result) == "try-error") {
      next
    }
    else {
      table[i,1] <- col$alpha
      table[i,2] <- col$se
      table[i,3] <- col$t
      table[i,4] <- col$pval
    }
  }
  return(table)
}

HC.lasso.mult <- function(x, y, d, I3=NULL,  ...) {
  y <- as.matrix(y, ncol=1)
  d <- as.matrix(d)
  if (is.null(colnames(d))) colnames(d) <- paste("d", 1:dim(d)[2], sep="")
  
  kd <- dim(d)[2]
  kx <- dim(x)[2]
  n <- dim(x)[1]
  
  I2 <- lasso(x, y, ...)$index
  I1 <- vector("logical", length=kx)
  for (i in 1:kd) {
  Id <- lasso(x, d[,i], ...)$index
  I1 <- I1 + Id
  }
  
  if (is.logical(I3)) {
    I <- I1+I2+I3
    I <- as.logical(I)
  } else {
    I <- I1+I2
    I <- as.logical(I)
  }
  
  if (sum(I)==0) {I <- NULL}
  x <- cbind(d,x[,I, drop=FALSE])
  reg1 <- lm(y~x)
  xi<-reg1$residuals*sqrt(n/(n-sum(I)-1))
  table <- matrix(NA,ncol=4,nrow=kd)
  rownames(table) <- colnames(d)
  colnames(table) <- c("coeff.", "se.", "t-value", "p-value")
  for (i in 1:kd)
  {
    if(is.null(I)) {reg2 <- lm(d[,i]~1)}
    if(!is.null(I)) {reg2 <- lm(d[,i]~x[,-i, drop=FALSE])}
    v <- reg2$residuals
    var <- 1/n*1/mean(v^2)*mean(v^2*xi^2)*1/mean(v^2)
    alpha <- coef(reg1)[1+i]
    table[i,1] <- alpha
    table[i,2] <- sqrt(var)
    table[i,3] <- alpha/sqrt(var)
    table[i,4] <- 2*pnorm(-abs(alpha/sqrt(var)))     
  }
  return(table)
}


