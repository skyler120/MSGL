library(MCMCpack)

# a toy example
  X = matrix(rnorm(3*4,mean=0, sd=1), 3,4)
  beta = c(0,0,2,2)
  
  eps = rnorm(3,mean=0,sd=1)

  Y = X%*%beta + + 0.5*eps
  gvec = c(1,1,2,2)


#X = res_small$x
#Y = res_small$y
#gvec = res_small$grp




sum_square <- function(x) sum(x^2) # function that calculates the length of a vector

N =100 # Total number of iterations


sgl_gibbs <- function(X,Y,gvec){

  
  # Initialize
  simb <- list()   # a list of list of vectors bg
  simtau <- list() # a list of vectors with entries tau
  simsigmasquare <- vector()
  simssquare <- vector()
  simPi1 <- vector()
  simPi0 <- vector()
  
  simPi0[1] = 0.5
  simPi1[1] = 0.5
  simssquare[1] = 0.1
  simsigmasquare[1] = 1
  simtau[[1]] <- integer(length(gvec))+0.1
  s = list()
  for (t in 1:length(unique(gvec))) {
  s = c(s, list(integer(sum(gvec==t))+0.1))
  }
  simb[[1]] = s
  
  
  # specify values of parameters
  c1=0.1
  c2=0.1
  a1=0.1
  a2=0.1
  t=0.1
  alpha=0.1
  gamma=0.1  
  # Gibbs sampler
  for(i in 2:N){ 
    simb[[i]]= list()
     for (g in 1: length(unique(gvec))){  # loop through number of groups
       simtau[[i]] <-  integer(length(gvec))+0.1
      simb[[i]][[g]] <- numeric(length(unique(gvec)))+0.1
       for (j in 1: sum(gvec==g)){ # loop through the group size
          
         betangj <- diag(simtau[[i-1]][-(g*j)])%*%Reduce(c, simb[[i-1]])[-(g*j)]
        
         vgjsquare <- ginv((1/simssquare[i-1])^2+ (1/simsigmasquare[i-1])^2*t(X[,(g*j)])%*%
           X[,(g*j)]*(Reduce(c, simb[[i-1]])[(g*j)])^2)
         
         ugj <- (1/simssquare[i-1])^2*vgjsquare*t(Y-(X[,-(g*j)]%*%betangj))%*%X[,(g*j)]*
           (Reduce(c, simb[[i-1]])[(g*j)])
         
         qgj <- simPi1[i-1]/(simPi1[i-1]+2*(1-simPi1[i-1])*(simssquare[i-1])^(-1/2)*(vgjsquare^(1/2))*
                               exp(ugj^2/(2*vgjsquare))*pnorm(ugj/sqrt(vgjsquare)))
         
         simtau[[i]][g*j] <- (runif(1)<1-qgj)*abs(rnorm(n = 1, mean=ugj, sd=vgjsquare))
          
       }
       
       #Vg <- diag(simtau[[i]][((sum(sum(gvec == 1:(g-1))))+1): sum(sum(gvec == 1:g))])
       
       if (g==1){
         Vg <- diag(simtau[[i-1]][1: sum(gvec %in% 1:g)])  # this has changed to %in% from ==
         
        
         Vng <- diag(simtau[[i-1]][-(1:sum(gvec %in% 1:g))])
         
         Xng <- X[,-(1: sum(gvec %in% 1:g))]
         Xg <-  X[,(1: sum(gvec %in% 1:g))]
         bng <- Reduce(c,simb[[i-1]][-g])
         
         Sigmag <- ginv(diag(sum(gvec == g)) + (1/(simsigmasquare[i-1]))*Vg%*%t(Xg)%*%Xg%*%Vg)
         rootSigmag <- (svd(Sigmag)$u)%*%diag(sqrt(svd(Sigmag)$d))%*%(svd(Sigmag)$v)
         # rootSigmag <- sqrtm(Sigmag)
         mug <- (1/simsigmasquare[i-1])*rootSigmag%*%Vg%*%t(Xg)%*%(Y-Xng%*%Vng%*%bng)
         
         lg <- simPi0[i-1]/(simPi0[i-1]+(1-simPi0[i-1])*sqrt(abs(det(Sigmag)))*
                              exp((1/(2*simsigmasquare[i-1]^4))*
                                    sum_square(rootSigmag%*%Vg%*%t(Xg)%*%(Y-Xng%*%Vng%*%bng)))) ####problem here
         
         simb[[i]][[g]] <-  (runif(1)<1-lg)*mvrnorm(n = 1, mug, Sigmag) # bg|rest for g=1,...G
         
       } 
       else {
       Vg <- diag(simtau[[i-1]][(sum(gvec %in% 1:(g-1))+1): sum(gvec %in% 1:g)])
       Vng <- diag(simtau[[i-1]][-((sum(gvec %in% 1:(g-1))+1): sum(gvec %in% 1:g))])
       Xng <- X[,-((sum(gvec %in% 1:(g-1))+1): sum(gvec %in% 1:g))]
       Xg <-  X[,((sum(gvec %in% 1:(g-1))+1): sum(gvec %in% 1:g))]
       bng <- Reduce(c,simb[[i-1]][-g])
       
       Sigmag <- ginv(diag(sum(gvec == g)) + (1/(simsigmasquare[i-1])*Vg%*%t(Xg)%*%Xg%*%Vg))
       rootSigmag <- (svd(Sigmag)$u)%*%diag(sqrt(svd(Sigmag)$d))%*%(svd(Sigmag)$v)
       #rootSigmag <- sqrtm(Sigmag)
       mug <- (1/simsigmasquare[i-1])*rootSigmag%*%Vg%*%t(Xg)%*%(Y-Xng%*%Vng%*%bng)
       
       lg <- simPi0[i-1]/(simPi0[i-1]+(1-simPi0[i-1])*sqrt(abs(det(Sigmag)))*
                            exp((1/(2*simsigmasquare[i-1]^4))*
                                  sum_square(rootSigmag%*%Vg%*%t(Xg)%*%(Y-Xng%*%Vng%*%bng))))
       
       simb[[i]][[g]] <-  (runif(1)<1-lg)*mvrnorm(n = 1, mug, Sigmag) # bg|rest for g=1,...G
       }
         
     }
    
     simsigmasquare[i] <- rinvgamma(1, n/2+alpha, (1/2)*sum_square(Y-X%*%Reduce(c,simb[[i]]))+gamma)
     
     for(g in 1:length(unique(gvec))){
       count=0
     if(all(simb[[i]][[g]]==0))
       {
      count = count+1 
     }
     }
     
     simPi0[i] <-  rbeta(1,count+a1, length(unique(gvec))-
                           count+a2)
     
     simPi1[i] <- rbeta(1, sum(simtau[[i]]==0)+c1, length(gvec)-sum(simtau[[i]]==0)+c2)
     
     simssquare[i] <- rinvgamma(1,1+(1/2)*sum(simtau[[i]]==0), t+(1/2)*sum_square(simtau[[i]]))
  }
  
  # burnin and thinning
  burnin = N/2
  thinning = 25
  
  # Gibbs sample
  sim = Reduce(c,simb[[N]])
  #plot(ts(sim))
  return(sim)
}

#####


