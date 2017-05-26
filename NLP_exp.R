library(invgamma)
library(MASS)
library(MCMCpack)
library(gglasso)
library(glmnet)
library(hms)
library(SGL)
source("basu_cv.R")
source("sgl_gibbs.R")
source("sgl_beta.R")
source("lasso.R")
source("slog_function.R")

#Read data
X = as.matrix(read.csv("mat.csv", header=F))
W = as.matrix(read.csv("W.csv", header=F))
H = as.matrix(read.csv("H.csv", header=F))
Y = scan(file="ratings.txt")[1:nrow(X)]
YY = as.matrix(Y)
vocab = scan(file="vocabulary.txt", what="character")

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

res_1 = list(v = vocab,
               sgl=sgl_beta(X, YY, grps),
               basu = basu(X, YY, grps),
               bayes = NULL,
               slog = slog(X, YY, grps),
               lasso = lasso(X, YY, grps))
saveRDS(res_1,file="Rating_Exp/betas")



