


###########################
#3.1

#read in bhat's, bstar 
#bhat.sgl has first col corresponding to 50, second to 150, third to 500
#bstar is a vector
groupsizes <- c(50, 150, 500)

bhat <- cbind(bhat.sgl, bhat.basu, bhat.lasso)
spac <- 1 + (ncol(bhat)/length(groupsizes))*(0:(length(groupsizes)-1))

par(mfrow=c(1,ncol(bhat)/length(groupsizes)))

for (name in c("precision", "recall", "mcc")) {
  f <- match.fun(name)
  
  if(name == "mcc") {
    lims <- c(-1,1)
  } else {
    lims <- c(0,1)
  }
  
  y <- apply(bhat, 2, function(bhat) f(bhat, bstar))
  plot(groupsizes, y[spac], type="l", col=1, ylim=lims, ylab=name, xlab = "group size")
  lines(groupsizes, y[1+spac], col=2)
  lines(groupsizes, y[2+spac], col=3)
  legend("topright", legend=c("SGL", "Basu", "Lasso"), lty=rep(1, 3), col=1:3)
}

###########################
#3.2

#read data like before

ksizes <- c(5, 50, 125)

bhat <- cbind(bhat.sgl, bhat.basu, bhat.lasso)
spac <- 1 + (ncol(bhat)/length(ksizes))*(0:(length(ksizes)-1))

par(mfrow=c(1,ncol(bhat)/length(ksizes)))

for (name in c("precision", "recall", "mcc")) {
  f <- match.fun(name)
  
  if(name == "mcc") {
    lims <- c(-1,1)
  } else {
    lims <- c(0,1)
  }
  
  y <- apply(bhat, 2, function(bhat) f(bhat, bstar))
  plot(ksizes, y[spac], type="l", col=1, ylim=lims, ylab=name, xlab = "within-group sparsity")
  lines(ksizes, y[1+spac], col=2)
  lines(ksizes, y[2+spac], col=3)
  legend("topright", legend=c("SGL", "Basu", "Lasso"), lty=rep(1, 3), col=1:3)
}

#############################
#5.3


#read data like before

index <- c(1, 10, 100, 1000)

bhat <- cbind(bhat.sgl, bhat.basu, bhat.lasso)
spac <- 1 + (ncol(bhat)/length(index))*(0:(length(s)-1))

par(mfrow=c(1,ncol(bhat)/length(index)))

for (name in c("precision", "recall", "mcc")) {
  f <- match.fun(name)
  
  if(name == "mcc") {
    lims <- c(-1,1)
  } else {
    lims <- c(0,1)
  }
  
  y <- apply(bhat, 2, function(bhat) f(bhat, bstar))
  plot(index, y[spac], type="l", col=1, ylim=lims, ylab=name, xlab = "element size")
  lines(index, y[1+spac], col=2)
  lines(index, y[2+spac], col=3)
  legend("topright", legend=c("SGL", "Basu", "Lasso"), lty=rep(1, 3), col=1:3)
}

####################################
#3.4

#read data like before
#each column of bhat.method is for a different setting

mat <- matrix(0, ncol=4, nrow=3) #4 is number of parameter settings, 3 is number of methods
rownames(mat) <- c("sgl", "basu", "lasso")

mat[1,] <- apply(bhat.sgl, 2, function(bhat) precision(bhat, bstar))
mat[2,] <- apply(bhat.basu, 2, function(bhat) precision(bhat, bstar))
mat[3,] <- apply(bhat.lasso, 2, function(bhat) precision(bhat, bstar))

print(mat)





