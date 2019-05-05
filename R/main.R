# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

##------self-defined function-----------------------

#-- generate sketch matrix
sketchfun <- function(m, n, type='ROS'){
  library(Matrix)
  Praw <- Diagonal(n, x = 1)
  ind <- sample(n, m)
  R <- Diagonal(n, 1/2 - rbinom(n, 1, 0.5))
  H <- Praw
  S <- switch(type,
              subGaussian = matrix(rnorm(n*m),m,n),
              ROS = sqrt(n/m)*Praw[ind,]%*% H %*% R,
              subSampling = sqrt(n/m)*Praw[ind,])
  return(as.matrix(S))
}
#-- assist to evaluate reconstruction error
recerr= function(x,data,alpha,alphaKrow,sumalpha,Ksum ,kern){
  # x <- tsdata
  # data <- trdata
  ns <- nrow(x)
  d <- length(sumalpha)
  if(!inherits(x,'matrix')) x <- as.matrix(x)
  if(!inherits(data,'matrix')) data <- as.matrix(data)
  kmat = kernelMatrix(kern, data,x)
  f <- t(kmat) %*% (alpha) - matrix(apply(kmat, 2,mean)- Ksum, ncol=1)%*% t(sumalpha) -
    matrix(alphaKrow, nrow=ns, ncol=d, byrow=T)
  f <- as.matrix(f)
  res <- diag(kernelMatrix(kern, x,x)) - 2*apply(kmat, 2, mean) + Ksum - apply(f*Conj(f),1, sum)
  return(as.vector(res))
}

#-- main function to conduct sketched KPCA.
sKPCA <- function(trdata, tsdata=NULL,d, m=ceiling(sqrt(nrow(trdata))), Stype='ROS',  kern=rbfdot(sigma=0.5), seed=2){
  require(kernlab)
  n <- nrow(trdata)
  K <- kernelMatrix(kern, as.matrix(trdata))
  Krow <- apply(K,2, mean)
  Ksum <- mean(Krow)
  K <- K - matrix(Krow, n,n, byrow=T) - matrix(Krow, n, n) + Ksum
  set.seed(seed) # reproducible
  S <- sketchfun(m,n, Stype)
  Kt <- S%*%K%*%t(S)

  eigv <- eigen(Kt, symmetric = T)
  alpha_tilde <- eigv$vectors[,1:d]%*% diag(1/sqrt(eigv$values[1:d]))
  alpha <- t(S) %*% alpha_tilde
  # norm(cbind(alpha[,1]),'f')
  alpha_rotate <- K %*% alpha # rotation for principal component scores.
  res <- list()
  res$alpha <- alpha
  res$alpha_rotate <- alpha_rotate
  sumalpha <- apply(alpha, 2, sum)
  alphaKrow <- Krow %*% alpha
  if(is.null(tsdata)){
    err <- recerr(trdata,trdata,alpha,alphaKrow,sumalpha,Ksum,kern)
  }else{
    err <- recerr(tsdata,trdata,alpha,alphaKrow,sumalpha,Ksum,kern)
  }
  res$rec_err <- err
  res$prop_var <- sum(eigv$values[1:d]) / sum(eigv$values)
  class(res) <- 'sKPCA'
  return(res)
}


# original KPCA
oKPCA <- function(trdata, tsdata=NULL,d, kern=rbfdot(sigma=0.5)){
  require(kernlab)
  n <- nrow(trdata)
  K <- kernelMatrix(kern, as.matrix(trdata))
  Krow <- apply(K,2, mean)
  Ksum <- mean(Krow)
  Kt <- K - matrix(Krow, n,n, byrow=T) - matrix(Krow, n, n) + Ksum
  eigv <- eigen(Kt, symmetric = T)
  alpha <- eigv$vectors[,1:d]
  alpha <- alpha %*% diag(1/sqrt(eigv$values[1:d]))

  alpha_rotate <- Kt %*% alpha # rotation for principal component scores.
  res <- list()
  res$alpha <- alpha
  res$alpha_rotate <- alpha_rotate
  sumalpha <- apply(alpha, 2, sum)
  alphaKrow <- Krow %*% alpha
  if(is.null(tsdata)){
    err <- recerr(trdata,trdata,alpha,alphaKrow,sumalpha,Ksum,kern)
  }else{
    err <- recerr(tsdata,trdata,alpha,alphaKrow,sumalpha,Ksum,kern)
  }
  res$rec_err <- err
  res$prop_var <- sum(eigv$values[1:d]) / sum(eigv$values)
  class(res) <- 'oKPCA'
  return(res)
}

#-- evaluate reconstruction error from fKPCA or oKPCA class.
recerrfun <- function(obj){
  if((!inherits(obj, 'sKPCA')) && (!inherits(obj, 'oKPCA')))
    stop('obj must be sKPCA or oKPCA class!')
  recerr <- obj$rec_err
  return(mean(recerr^2))
}

#--
vec2difmat <- function(x){
  nx <- length(x)
  matrix(x, nrow=nx, ncol=nx) - matrix(x, nrow=nx, ncol=nx, byrow=T)
}

sigmaslect <- function(X){
  p <- ncol(X)
  n <- nrow(X)
  res <- array(0, c(n,n, p))
  for(j in 1:p){
    res[,,j] <- vec2difmat(X[,j])
  }
  return(sqrt(sum(res^2))/n^2)
}
