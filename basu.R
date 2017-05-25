basu <- function(X, y, group) {
  mod <- cv.gglasso(X, y, group, pred.loss="L2")
  mod <- gglasso(X, y, group,lambda=mod$lambda.min) #fits GL via cv
  
  beta <- hbeta <- mod$beta
  
  logscale <- function(len=30, lammin=0.001, lammax=10) {lammin * (lammax/lammin)^seq(from=0, to=1, length.out=len)} #used to get lamlist
  
  #choose optimal delta (threshold) via cv
  delta.depth <- 10
  pred.error <- rep(0, length=delta.depth)
  delta.grid <- logscale(len=delta.depth, lammin=.01, lammax=1)
  for (j in 1:length(delta.grid)) {
    delta <- delta.grid[j]
    for (i in 1:group[length(group)]) {
      g <- (group==i)
      beta[g] <- hbeta[g] * (abs(hbeta[g]) > 2*delta*sqrt(sum((hbeta[g])^2)))*1
    }
    pred.error[j] <- norm(y - X %*% beta)
  }
  delta <- delta.grid[which.min(pred.error)] #now refit the best one
  for (i in 1:group[length(group)]) {
    g <- (group==i)
    beta[g] <- hbeta[g] * (abs(hbeta[g]) > 2*delta*sqrt(sum((hbeta[g])^2)))*1
  }
  
  return(beta)
}