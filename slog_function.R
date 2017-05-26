slog <- function(X, y, group, alpha=0.05) {
  #######################
  #transform group to specify the DAG and prox
  #my working out of these
  gsizes <- rle(group)$lengths
  
  var <- c(as.list(1:p), vector("list", length(gsizes)))
  w <- c(rep(1,p), sqrt(gsizes))
  map <- cbind(1:p, rep(p,p) + group)
  
  prox <- prox.log <- function(grad, lam, alpha) {hsm::hsm(grad, lam, w=c(rep(alpha, p), rep(1-alpha, length(gsizes)))*w, map=map, var=var)$beta} #xiaohan's way to compute LOG prox
  
  ############################
  #PGD
  
  est <- function(lam, alpha, prox, y, X, beta = rep(0, length=ncol(X)), eps=.00001) { #inputting phi allows warmstarts
    s <- base::norm(X, type="2")^(-2)
    
    beta2 <- beta1 <- beta
    Xty <- t(X) %*% y
    for (m in 1:1000) {
      beta <- beta1 + (m-2)/(m+1)*(beta1 - beta2)
      beta2 <- beta1
      
      grad <- beta + s*(Xty - t(X) %*% (X %*% beta))
      beta1 <- prox(grad, lam, alpha)
      if (max(abs(beta - beta1)) <= eps) {
        break
      }
      if (m == 1000) {
        warning("M GOT TO 1000!")
      }
    }
    return(beta)
  }
  
  cv <- function(y, X, fold.num=5) {
    dat <- cbind(y,X)[sample(nrow(X)),]
    
    #Create 10 equally size folds
    folds <- cut(seq(1,nrow(dat)),breaks=fold.num,labels=FALSE)
    
    #Perform 10 fold cross validation
    for(i in 1:fold.num){
      #Segement your data by fold using the which() function 
      inds <- which(folds==i,arr.ind=TRUE)
      test <- dat[inds, ]
      train <- dat[-inds, ]
      #Use the test and train data partitions however you desire...
    }
    
    return(list(train=train, test=test))
  }
  
  
  ooserror <- function(lam, alpha, prox, y, X, beta = rep(0, length=ncol(X)), fold.num=5, traintest=NA) {
    if (any(is.na(traintest))) {
      obj <- cv(y, X, fold.num)
      train <- obj$train
      test <- obj$test
    } else {
      train <- traintest$train
      test <- traintest$test
    }
    
    
    beta <- est(lam, alpha, prox, train[,1], train[,-1], beta = beta)
    return(list(error=sqrt(sum((test[,1] - test[,-1] %*% beta)^2))/sqrt(nrow(test)), beta=beta))
  }
  
  tuning <- function(lamlist, alpha, prox, y, X, fold.num=5, traintest=NA) {
    if (any(is.na(traintest))) {
      obj <- cv(y, X, fold.num)
      train <- obj$train
      test <- obj$test
    } else {
      train <- traintest$train
      test <- traintest$test
    }
    
    
    beta <- rep(0, ncol(X))
    
    err <- c()
    for (i in length(lamlist):1) { #assumes lamlist is increasing
      lam <- lamlist[i]
      res <- ooserror(lam, alpha, prox, y, X, beta=beta, traintest=list(train=train, test=test))
      err[i] <- res$error
      beta <- res$beta
    }
    
    ind <- which.min(err)
    return(ooserror(lamlist[ind], alpha, prox, y, X, beta=beta, traintest=list(train=train, test=test))$beta)
  }
  
  logscale <- function(len=30, lammin=0.001, lammax=10) {lammin * (lammax/lammin)^seq(from=0, to=1, length.out=len)} #used to get deltagrid
  
  return(tuning(logscale(len=2), alpha=alpha, prox, y, X))
}



