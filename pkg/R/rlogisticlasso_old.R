# #' rlogisticlasso: Function for logistic Lasso estimation
# #'
# #' The function estimates the coefficients of a logistic Lasso regression with
# #' data-driven penalty. The method of the data-driven penalty can be chosen.
# #' The object which is returned is of the S3 class \code{rlogisticlasso}
# #'
# #' The function estimates the coefficients of a Logistic Lasso regression with
# #' data-driven penalty. The method of the data-driven penalty can be chosen
# #' ("X-dependent", "X-independent", "standard", "none"). The description of the
# #' methods is given in \code{\link{rlasso}}, with the minor difference that here
# #' also for the "X-dependent" and "X-independent" method data-driven loadings
# #' are employed (instead of an iteratively estimated variance factor). The
# #' option \code{post=TRUE} conducts post-lasso estimation, i.e. a refit of the
# #' model with the selected variables.
# #' @param x regressors (matrix)
# #' @param y dependent variable (vector or matrix)
# #' @param post logical. If \code{TRUE}, post-lasso estimation is conducted.
# #' @param intercept logical. If \code{TRUE}, intercept is included which is not
# #' penalized.
# #' @param normalize logical. If \code{TRUE}, design matrix \code{x} is scaled.
# #' @param penalty list with options for the calculation of the penalty.  \code{c} and \code{gamma} constants for the penalty (all methods except "CV"), \code{method} method for penalty choice ("X-dependent",
# #' "X-independent", "standard", "none", "CV"), \code{numfolds} number of folds
# #' for the "cross validation" method, \code{numSim} number of simulations for
# #' the "X-dependent" method, \code{lambda.start} initial penalization value, compulsory for method "none"
# #' @param control list with control values.
# #' \code{numIter} number of iterations for the algorithm for
# #' the estimation of the variance and data-driven penalty, ie. loadings,
# #' \code{tol} tolerance for improvement of the estimated variances. \code{threshold} is applied to the final estimated lasso
# #' coefficients. Absolute values below the threshold are set to zero.
# #' @param ... further parameters passed to glmnet
# #' @return \code{rlogisticlasso} returns an object of class
# #' \code{rlogisticlasso} An object of class \code{rlogisticlasso} is a list
# #' containing at least the following components: \item{coefficients}{parameter
# #' estimates (without intercept)} \item{intercept.value}{value of intercept} \item{index}{index of selected variables (logicals)}
# #' \item{lambda}{data-driven penalty term for each variable, product of lambda0
# #' (the penalization parameter) and the loadings} \item{lambda0}{penalty term}
# #' \item{loadings}{loading for each regressor} \item{residuals}{residuals}
# #' \item{sigma}{root of the variance of the residuals}
# #' \item{iter}{number of iterations} \item{call}{function call}
# #' \item{options}{options}
# #' @references TBD
# #' @keywords logistic lasso lasso logistic regression
# #' @export
# #' @rdname rlogisticlasso
# rlogisticlasso <- function(x, ...)
#   UseMethod("rlogisticlasso") # definition generic function
#
#
# #' @rdname rlogisticlasso
# #' @export
#
# rlogisticlasso.default <- function(x, y, post=TRUE, intercept=TRUE, normalize=TRUE,  penalty = list(method = "standard", lambda.start = NULL, c = 1.1, gamma = 0.1),
#                            control = list(numIter = 15, tol = 10^-5, threshold = NULL),...) {
#   n <- dim(x)[1]
#   p <- dim(x)[2]
#   if (is.null(colnames(x))) colnames(x) <- paste("V",1:p, sep="")
#   ind.names <- 1:p
#
#   # checking input numIter, tol
#   if (!exists("numIter", where = control)) {
#     control$numIter = 15
#     message("numIter in control not provided. Set to default 15")
#   }
#
#   if (!exists("tol", where = control)) {
#     control$tol = 10^-5
#     message("tol in control not provided. Set to default 10^-5")
#   }
#
#   pen <- lambdaCalculation(penalty = penalty, y = y, x = x)
#   lambda <- pen$lambda
#   Ups0 <- Ups1 <- as.vector(pen$Ups0)
#   lambda0 <- pen$lambda0/(2*n)
#
#   ###
#   mm <- 1
#   s0 <- sqrt(var(y))
#   while (mm < control$numIter) {
#     # calculation parameters
#     log.lasso <- glmnet::glmnet(x, y, family=c("binomial"), alpha = 1, lambda = 0.5*lambda0[1],
#                         standardize = normalize, intercept = intercept, thresh = control$tol, penalty.factor = Ups1)
#     #log.lasso$beta <- log.lasso$beta/Ups1
#     coefTemp <- as.vector(log.lasso$beta)
#     coefTemp[is.na(coefTemp)] <- 0
#     ind1 <- (abs(coefTemp) > 0)
#     x1 <- as.matrix(x[,ind1, drop=FALSE])
#     if (dim(x1)[2]==0) {
#       print("No variables selected!")
#       est <- list(coefficients=rep(0,p), a0=mean(y), index=rep(FALSE,p), s0=s0, lambda=lambda, lambda0=pen$lambda0, loadings=Ups1, residuals=y - mean(y), sigma=var(y), iter=0, call=match.call(),
#                   options=list(post=post, intercept=intercept, normalize=normalize, control=control))
#       class(est) <- c("rlogisticlasso")
#       return(est)
#     }
#
#     # refinement variance estimation
#     if (post) {
#       if(intercept){
#         #reg <- glm.fit(x1,y, family = binomial(link = "logit"), intercept=intercept)
#         reg <- glm(y~x1, family = binomial(link = "logit"))
#         coefT <- coef(reg)[-1]
#         coefT[is.na(coefT)] <- 0
#         e1 <- y-reg$fitted.values
#         coefTemp[ind1] <- coefT
#       }
#
#       if(!intercept){
#         reg <- glm(y~-1+x1, family = binomial(link = "logit"))
#         coefT <- coef(reg)
#         coefT[is.na(coefT)] <- 0
#         e1 <- y-reg$fitted.values
#         coefTemp[ind1] <- coefT
#       }
#     }
#     if (!post) {
#       e1 <- y- predict(log.lasso, newx=x ,type="response")
#     }
#
#       s1 <- sqrt(var(e1))
#
#       if(penalty$method=="standard") Ups1 <- 1/sqrt(n)*as.vector(sqrt(t(t(e1^2)%*%(x^2))))
#       if(penalty$method=="none") Ups1 <- 1/sqrt(n)*as.vector(sqrt(t(t(e1^2)%*%(x^2))))
#       if(penalty$method=="X-dependent") Ups1 <- s1
#       if(penalty$method=="X-independent") Ups1 <- s1
#
#       mm <- mm+1
#       if(abs(s0-s1)<control$tol) {break}
#       s0 <- s1
#     }
#     if (intercept==TRUE) {
#       if (post==TRUE) a0 <- coef(reg)[1]
#       if (post==FALSE) a0 <- coef(log.lasso)[1]
#     }
#
#   if (intercept==FALSE) {
#     a0 <- NA
#   }
#
#     coefTemp <- as.vector(coefTemp)
#     coefTemp[abs(coefTemp)<control$threshold] <- 0
#     ind1 <- as.vector(ind1)
#     names(coefTemp) <- names(ind1) <- colnames(x)
#     est <- list(coefficients=coefTemp, a0=a0, index=ind1, lambda=lambda, lambda0=lambda0, loadings=Ups1, residuals=e1, sigma=s1, iter=mm, call=match.call(), options=list(post=post, intercept=intercept, normalize=normalize, control=control))
#     class(est) <- c("rlogisticlasso")
#     return(est)
# }
#
# #' @rdname rlogisticlasso
# #' @param formula an object of class "formula" (or one that can be coerced to
# #' that class): a symbolic description of the model to be fitted in the form
# #' \code{y~x}
# #' @param data an optional data frame, list or environment (or object coercible
# #' @export
# rlogisticlasso.formula <- function(formula, data, post = TRUE, intercept = TRUE,
#                            normalize = TRUE, penalty = list(method = "standard", lambda.start = NULL,
#                                                             c = 1.1, gamma = 0.1), control = list(numIter = 15, tol = 10^-5,
#                                                                                                   threshold = NULL), ...) {
#   cl <- match.call()
#   mf <- match.call(expand.dots = FALSE)
#   m <- match(c("formula", "data"), names(mf), 0L)
#   mf <- mf[c(1L, m)]
#   mf$drop.unused.levels <- TRUE
#   mf[[1L]] <- quote(stats::model.frame)
#   mf <- eval(mf, parent.frame())
#   mt <- attr(mf, "terms")
#   attr(mt, "intercept") <- 0
#   y <- model.response(mf, "numeric")
#   x <- model.matrix(mt, mf)
#
#   est <- rlogisticlasso(x, y, post = post, intercept = intercept, normalize = normalize, penalty=penalty,
#                 control = control)
#   est$call <- cl
#   return(est)
# }
#
# ###############################################################################################################################
#
# ################# Methods for logistic Lasso
#
# #' Methods for S3 object \code{rlogisticlasso}
# #'
# #' Objects of class \code{rlogisticlasso} are construced by \code{rlogisticlasso.formula} or \code{rlogistic.default}.
# #' \code{print.rlogisticlasso} prints and displays some information about fitted \code{rlogisticlasso} objects.
# #' \code{summary.rlogisticlasso} summarizes information of a fitted \code{rlasso} object.
# #' \code{predict.rlogisticlasso} predicts values based on a lasso object.
# #' \code{model.matrix.rlogisticlasso} constructs the model matrix of a lasso object.
# #' @param object an object of class \code{rlogisticlasso}
# #' @param x an object of class \code{rlogisticlasso}
# #' @param all logical, indicates if coefficients of all variables (TRUE) should be displayed or only the non-zero ones (FALSE)
# #' @param digits significant digits in printout
# #' @param type the type of prediction required. The default ("response) is on the scale of the response variable; the alternative "link" is on the scale of the linear predictors.
# #' @param newdata new data set for prediction
# #' @param ... arguments passed to the print function and other methods
# #' @keywords methods rlogisticlasso
# #' @rdname methods.rlogisticlasso
# #' @aliases methods.rlogisticlasso print.rlogisticlasso summary.rlogisticlasso predict.rlogisticlasso model.matrix.rlogisticlasso
# #' @export
#
# predict.rlogisticlasso <- function (object, newdata=NULL, type="response", ...) {
#   if (missing(newdata) || is.null(newdata)) {
#     X <- model.matrix(object)
#   } else {
#     varcoef <- names(object$coefficients)
#       if(all(is.element(varcoef,colnames(newdata)))){
#         X <- as.matrix(newdata[,varcoef])
#         } else {
#           X <- as.matrix(newdata)
#            stopifnot(ncol(X)==length(object$coefficients))
#         }
#   }
#   n <- dim(X)[1] #length(object$residuals)
#   beta <- object$coefficients
#   if (object$options[["intercept"]]) {
#     yp <- object$a0 + X%*%as.vector(beta)
#     if (dim(X)[2]==0) yp <- rep(object$a0,n)
#     if (type=="response") yhat <- 1/(1+exp(-yp))
#     if (type=="link") yhat <- yp
#   }
#   if (!object$options[["intercept"]]) {
#     yp <- X%*%as.vector(beta)
#     if (dim(X)[2]==0) yp <- rep(0,n)
#     yhat <- 1/(1+exp(-yp))
#     if (type=="response") yhat <- 1/(1+exp(-yp))
#     if (type=="link") yhat <- yp
#   }
#   return(yhat)
# }
#
# # predict.rlogisticlasso <- function (object, newdata=NULL, ...) {
# #   if (missing(newdata) || is.null(newdata)) {
# #     X <- model.matrix(object)
# #   }
# #   else {
# #     #mt <- attr(newdata, "terms")
# #     #attr(mt, "intercept") <- 0
# #     #m <- model.frame(mt, newdata)
# #     #X <- model.matrix(mt, m)
# #     X <- newdata
# #   }
# #   n <- length(object$residuals)
# #   beta <- object$coefficients
# #   if (object$options[["intercept"]]) {
# #     yp <- object$a0 + X%*%as.vector(beta)
# #     yhat <- 1/(1+exp(-yp))
# #   }
# #   if (!object$options[["intercept"]]) {
# #     yp <- X%*%as.vector(beta)
# #     yhat <- 1/(1+exp(-yp))
# #   }
# #   return(yhat)
# # }
#
#
# #' @rdname methods.rlogisticlasso
# #' @export
#
# model.matrix.rlogisticlasso <- function(object, ...) {
#
#   # falls formula
#   if (is.call(object$call[[2]])) {
#     # falls kein Datensatz uebergeben
#     if(is.null(object$call$data)){
#       X <- model.frame(object$call[[2]])
#       mt <- attr(X, "terms")
#       attr(mt, "intercept") <- 0
#       mm <- model.matrix(mt, X)
#       # falls Datensatz uebergeben
#     } else {
#       dataev <- eval(object$call$data)
#       mm <- as.matrix(dataev[,names(object$coefficients)])
#     }
#   } else {
#     # falls default
#     mm <- eval(object$call[[2]])
#   }
#   return(mm)
# }
#
#
# #' @rdname methods.rlogisticlasso
# #' @export
#
# print.rlogisticlasso <- function(x, all=TRUE ,digits = max(3L, getOption("digits") - 3L), ...) {
#   cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
#   if (length(coef(x))) {
#     coeffs <- coef(x)
#     if (x$options$intercept) {
#       coeffs <- c(x$a0, coeffs)
#       names(coeffs)[1] <- "(Intercept)"
#       index <- cbind(1,x$index)
#     }
#     if (all) {
#       cat("Coefficients:\n")
#       print.default(format(coeffs , digits = digits), print.gap = 2L,
#                     quote = FALSE)
#     } else {
#       print.default(format(coeffs[index], digits = digits), print.gap = 2L,
#                     quote = FALSE)
#     }
#   }
#   else cat("No coefficients\n")
#   cat("\n")
#   invisible(x)
# }
#
# #' @rdname methods.rlogisticlasso
# #' @export
#
# summary.rlogisticlasso <- function(object, all=TRUE, digits = max(3L, getOption("digits") - 3L), ...) {
#   cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n", sep = "")
#   cat("\nPost-Lasso Estimation: ",  paste(deparse(object$options$post), sep = "\n", collapse = "\n"), "\n", sep = " ")
#   coefs <- object$coefficients
#   p <- length(coefs)
#   num.selected <- sum(abs(object$coefficients)>0)
#   cat("\nTotal number of variables:", p)
#   cat("\nNumber of selected variables:", num.selected, "\n", sep=" ")
#   #resid <- object$residuals
#   #cat("\nResiduals: \n")
#   #nam <- c("Min", "1Q", "Median", "3Q", "Max")
#   #rq <- structure(apply(t(resid), 1L, quantile), dimnames = list(nam, dimnames(resid)[[2L]]))
#   #print(drop(t(rq)), digits = digits)
#   cat("\n")
#   if (all) {
#     coefm <- matrix(NA, length(coefs), 1)
#     coefm[,1] <- coefs
#     colnames(coefm) <- "Estimate"
#     rownames(coefm) <- names(coefs)
#     printCoefmat(coefm, digits = digits, na.print = "NA")
#   } else {
#     coefs <- coefs[abs(coefs)>0]
#     coefm <- matrix(NA, length(coefs), 1)
#     coefm[,1] <- coefs
#     colnames(coefm) <- "Estimate"
#     rownames(coefm) <- names(coefs)
#     printCoefmat(coefm, digits = digits, na.print = "NA",  P.values=TRUE, has.Pvalue=TRUE)
#   }
#   cat("\n")
#   invisible(object)
# }
