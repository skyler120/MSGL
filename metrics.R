pn <- function(bhat, bstar) {
  zer <- bhat[bstar==0]; nzer <- bhat[bstar != 0]
  return(list(tp = sum(nzer!=0), tn = sum(zer==0), fp = sum(zer!=0), fn = sum(nzer==0)))
}

precision <- function(bhat, bstar) {
  with(pn(bhat, bstar), return(tp/(tp+fp)))
}

recall <- function(bhat, bstar) {
  with(pn(bhat, bstar), return(tp/(tp+fn)))
}

mcc <- function(bhat, bstar) {
  with(pn(bhat, bstar), return((tp*tn - fp*fn)/((tp + fp)^(1/2)*(tp + fn)^(1/2)*(tn + fp)^(1/2)*(tn + fn)^(1/2))))
}
