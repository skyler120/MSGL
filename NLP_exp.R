source("basu.R")

#Read data
X = as.matrix(read.csv("mat.csv", header=F))
W = as.matrix(read.csv("W.csv", header=F))
H = as.matrix(read.csv("H.csv", header=F))
Y = scan(file="ratings.txt")[1:nrow(X)]
P <- nrow(X)
Q <- ncol(X)
G <- ncol(W)

#Get groups as most common words
grps = numeric(Q)
for(i in 1:Q){
  grps[i] = which.max(H[,i])
}


#Use SGL logit
#YY = as.matrix(1*(Y>=3))
#dat = list(x = X, y = YY)
#sgl_res = SGL(dat, grps, type="logit")

#OR SGL linear
YY = as.matrix(Y)
dat = list(x = X, y = YY)
beta_sgl = SGL(dat, grps, type="linear")$beta

#basu linear
YY = as.matrix(Y)
beta_basu = basu(X,YY,grps)


