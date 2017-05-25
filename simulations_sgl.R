#simulations
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


############## Get Results ##################################
get_res_mis <- function(res, file_name){
  res_1 = list(data = res,
               sgl=sgl_beta(res$x, res$y, res$grps_mis),
               basu = basu(res$x, res$y, res$grps_mis),
               bayes = NULL,
               slog = NULL,
               lasso = NULL)
  saveRDS(res_1,file=file_name)
}
get_res <- function(res, file_name){
  res_1 = list(data = res,
               sgl=sgl_beta(res$x, res$y, res$grps),
               basu = basu(res$x, res$y, res$grps),
               bayes = NULL,
               slog = slog(res$x, res$y,res$grps),
               lasso = lasso(res$x, res$y,res$grps))
  saveRDS(res_1,file=file_name)
}
###############Varying Group Sizes############################
gen_varying_grps <- function(g,n,p){
  X = matrix(rnorm(n*p,mean=0, sd=1), n,p)
  beta_ls = c(1:5, rep(0,p/g-5))
  beta = rep(beta_ls,g)
  grps = rep(1:g,each=p/g)
  eps = rnorm(n,mean=0,sd=1)
  Y = numeric(n)
  for(l in 1:g){
    Y = Y + X[,grps==l]%*%beta[grps==l] 
  }
  Y = Y + + 0.5*eps
  return(list(x = X, y = Y, beta = beta, grps = grps))
}

#small group
N = 60; P = 1500; G = 3
res_small = gen_varying_grps(G,N,P)
get_res(res_small, 'Results/res_small')
#medium group
N = 60; P = 1500; G = 10
res_medium = gen_varying_grps(G,N,P)
get_res(res_medium, 'Results/res_medium')
#large group
N = 60; P = 1500; G = 30
res_large = gen_varying_grps(G,N,P)
get_res(res_large, 'Results/res_large')



###############Mostly zero to mostly not############################
gen_varying_empty <- function(k,n,p){
  g = 10
  X = matrix(rnorm(n*p,mean=0, sd=1), n,p)
  beta_ls = c(1:k, rep(0,p/g-k))
  beta = rep(beta_ls,g)
  grps = rep(1:g,each=p/g)
  eps = rnorm(n,mean=0,sd=1)
  Y = numeric(n)
  for(l in 1:g){
    Y = Y + X[,grps==l]%*%beta[grps==l] 
  }
  Y = Y + + 0.5*eps
  return(list(x = X, y = Y, beta = beta, grps = grps))
}

N = 60; P = 1500; K = 5
res_low = gen_varying_empty(K,N,P)
get_res(res_low, 'Results/res_low')

#medium group
N = 60; P = 1500; K = 50
res_middle = gen_varying_empty(K,N,P)
get_res(res_middle, 'Results/res_middle')
#large group
N = 60; P = 1500; K = 125
res_high = gen_varying_empty(K,N,P)
get_res(res_high, 'Results/res_high')


###############Beta Size############################
gen_varying_bsize <- function(k,n,p){
  g = 10
  X = matrix(rnorm(n*p,mean=0, sd=1), n,p)
  beta_ls = c((1:10)*k,rep(0.1,20), rep(0,p/g-30))
  beta = rep(beta_ls,g)
  grps = rep(1:g,each=p/g)
  eps = rnorm(n,mean=0,sd=1)
  Y = numeric(n)
  for(l in 1:g){
    Y = Y + X[,grps==l]%*%beta[grps==l] 
  }
  Y = Y + + 0.5*eps
  return(list(x = X, y = Y, beta = beta, grps = grps))
}

N = 60; P = 1500; K = 1
res_bl = gen_varying_bsize(K,N,P)
get_res(res_bl, 'Results/res_bl')
#medium group
N = 60; P = 1500; K = 10
res_bm = gen_varying_bsize(K,N,P)
get_res(res_bm, 'Results/res_bm')
#large group
N = 60; P = 1500; K = 100
res_bh = gen_varying_bsize(K,N,P)
get_res(res_bh, 'Results/res_bh')
#VLarge group
N = 60; P = 1500; K = 1000
res_bvh = gen_varying_empty(K,N,P)
get_res(res_bvh, 'Results/res_bvh')

###############Relative Group Size############################
#Groups must be larger than 5 and sum to p
gen_varying_grpsize <- function(grp_sizes,n,p){
  g = length(grp_sizes)
  X = matrix(rnorm(n*p,mean=0, sd=1), n,p)
  beta = c()
  grps = c()
  for(m in 1:length(grp_sizes)){
    beta_ls = c(1:5, rep(0,grp_sizes[m]-5))
    beta = c(beta,beta_ls)
    grps = c(grps, rep(m, grp_sizes[m]))
  }
  eps = rnorm(n,mean=0,sd=1)
  Y = numeric(n)
  for(l in 1:g){
    Y = Y + X[,grps==l]%*%beta[grps==l] 
  }
  Y = Y + + 0.5*eps
  return(list(x = X, y = Y, beta = beta, grps = grps))
}

#all small 1 large
N = 60; P = 1500; GRPSIZE = c(rep(10,9),p-90)
res_gl = gen_varying_grpsize(GRPSIZE,N,P)
get_res(res_gl, 'Results/res_gl')
#all equal size
N = 60; P = 1500; GRPSIZE = rep(p/10,10)
res_gm = gen_varying_grpsize(GRPSIZE,N,P)
get_res(res_gm, 'Results/res_gm')
#half small half large
N = 60; P = 1500; GRPSIZE = rep(c(10,290),each=5)
res_hshl = gen_varying_grpsize(GRPSIZE,N,P)
get_res(res_hshl, 'Results/res_hshl')
#Descending groups
N = 60; P = 1500; GRPSIZE = c(500, 400,300,200,50,rep(10,5))
res_desc = gen_varying_grpsize(GRPSIZE,N,P)
get_desc(res_bvh, 'Results/res_desc')

###############Proportion of Misspecified Groups############################
gen_misspec_grps <- function(n,p,expmnt){
  X = matrix(rnorm(n*p,mean=0, sd=1), n,p)
  if(expmnt==0){
    beta = c(1:5, rep(0,p-5))
    grps = c(rep(1,5),rep(2,p-5))
  }else if(expmnt==1){
    beta = c(1:5, rep(0,p-5))
    grps = c(rep(1,2),rep(2,3), rep(3,p-5))
  }else if(expmnt==2){
    beta = c(1:5, rep(0,p-5))
    grps = c(1:5, rep(6,p-5))
  }else if(expmnt==3){
    beta = c(1:3, rep(0,(p/2-3)), 4:5, rep(0,p/2-2))
    grps = c(rep(1,p/2),rep(2,2),rep(3,p/2-2))
  }
  else{
    beta = c(1, rep(0,(p/5-1)), 2, rep(0,p/5-1), 3, rep(0,p/5-1), 4, rep(0,p/5-1), 5, rep(0,p/5-1))
    grps = rep(1:5,each=p/5)
  }
  eps = rnorm(n,mean=0,sd=1)
  Y = numeric(n)
  #for(l in 1:max(grps)){
  #  Y = Y + X[,grps==l]%*%beta[grps==l] 
  #}
  Y = X%*%beta
  Y = Y + + 0.5*eps
  return(list(x = X, y = Y, beta = beta, grps = grps))
}

#two misspecified groups
N = 60; P = 1500;
res_0mg = gen_misspec_grps(N,P,0)
get_res(res_0mg, 'Results/res_nomis')
res_1mg = gen_misspec_grps(N,P,1)
get_res(res_1mg, 'Results/res_2nomis')
res_2mg = gen_misspec_grps(N,P,2)
get_res(res_2mg, 'Results/res_5nomis')
res_3mg = gen_misspec_grps(N,P,3)
get_res(res_3mg, 'Results/res_smallmis')
res_4mg = gen_misspec_grps(N,P,4)
get_res(res_3mg, 'Results/res_lotsmis')
