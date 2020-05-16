################################################
rm(list=ls())
library(MASS)
library(Matrix)
my.filepath = ""
rowNorms = function(M){
  return(apply(M,1,FUN=function(x){norm(x,"2")}))
}
############################
######    Functions    #####
##### Run code start at ####
#####     Line 240      ####
############################
run_uneven_mod <- function(numberofgroups,seed,kappa){
  ############################
  # Fix sparsity, all dense shuffle, implement lambda based on groupsize 
  # Higher penalty for higher groupsize
  ############################
  nlambda = 50
  set.seed(seed)
  # true coef W and group info are from real data (summstat_final.RData)
  W.true <- BetaBALL
  grp.info.save <- grp.info 
  p = 300
  mykappa = kappa
  grp.select <- 1:numberofgroups
  grp.info = grp.info.save[grp.info.save$g.index%in%c(grp.select),] 
  N = nrow(grp.info)
  ugrp = unique(grp.info$g.index)
  n = length(ugrp)
  
  # Generate Data
  specX.bygrp <- NULL
  X.mean  <-  NULL
  grpcenters = mvrnorm(n=n,mu=rep(0,p),Sigma=diag(rep(1,p)))
  
  ########################################
  # start for loop
  for(i in 1:n){
    grpsize = length(which(grp.info$g.index==ugrp[i]))
    
    # upweight i-th grp in the mixture 
    myalpha <- rep(1,n); 
    myalpha[i] <-  2 
    dat.grp <- rmovMF(n=grpsize,theta=mykappa*grpcenters/apply(grpcenters,1,norm,"2"),alpha=myalpha)
    X.mean <- rbind(X.mean,dat.grp)
    specX.bygrp <- cbind(specX.bygrp,range(svd(dat.grp)$d))
  }
  # end for loop
  ########################################
  
  # specify true Pi: sparse shuffle within group 
  A <- diag(1,N,N)
  d=2
  print(mykappa) 
  ########################################
  # start for loop     
  # Modify to different sparsity (%1-m mapping) for each group 
  # n is the number of groups
  # extreme case: a group only have pure 1-1 or pure 1-m; equal group size
  ##each group has (35*2)% mismatch; either pure 1-1 or pure 1-m according to pure.1to1
  sparsity =rep(0.35,numberofgroups)
  pure.1to1 = sample(c(0,1),size=numberofgroups,replace=TRUE) #rep(1,numberofgroups)#
  ### tmpweight is defined as the proportion of 1-1 mapping in each group
  tmpweight=code1=code2 <- NULL
  for(i in 1:n){
    ind <- which(grp.info$g.index==ugrp[i])
    n.shuffle = floor(length(ind)*sparsity[i]/(d/2))
    if(n.shuffle>0){
      ind.shuffle = sample(ind,size=min(length(ind),n.shuffle*d),replace=F)
      for(j in 1:n.shuffle){
        if(pure.1to1[i]==1){
          A[ind.shuffle[j],ind.shuffle[n.shuffle*d+1-j]]=1
          A[ind.shuffle[n.shuffle*d+1-j],ind.shuffle[j]]=1
          A[ind.shuffle[j],ind.shuffle[j]]=0
          A[ind.shuffle[n.shuffle*d+1-j],ind.shuffle[n.shuffle*d+1-j]]=0
        }else{
          w1=runif(n=length(ind)); w1=w1/norm(w1,type="2")
          A[ind.shuffle[j],ind]=w1
          w2=runif(n=length(ind)); w2=w2/norm(w2,type="2")
          A[ind.shuffle[n.shuffle+j],ind]=w2
        }
      }#end j in 1:n.shuffle
      prop.1to1=length(which(apply(A[ind,ind],1,max)==1))/length(ind)
      tmpweight=c(tmpweight,prop.1to1)
    }else{
      tmpweight=c(tmpweight,1)
    }#end if(n.shuffle>0)
  }# end for loop
  ########################################
  S = length(which(diag(A)!=1)) ##how many mismatch
  N-S;S/N;log(S)/log(N) ## how many match
  X.save=X=X.mean
  Y.mean = A%*%X%*%W.true
  
  # Generate X and Y that follow VMF  
  Y <- t(apply(Y.mean/rowNorms(Y.mean),1,FUN=function(y){
    rmovMF(n=1,theta=mykappa*y)
  }))  
  Y.save <- Y   	
  index = as.vector(A%*%c(1:N)) ## the permutation vector
  
  #################################################
  ####### Step 1: Estimate W 
  #################################################
  Beta.spr = gradient_update_nogrp(X,Y,alpha=1,convergence=1e-10)
  Yhat.spr = X%*%Beta.spr
  (Sbias.spr = norm(Beta.spr - W.true,type="2"))
  (Fbias.spr = norm(Beta.spr - W.true,type="F"))##influenced by kappa and A
  orig = c(Truenorm = norm(W.true,type="F"),estnorm = norm(Beta.spr,type="F"))
  
  
  ## method 2: OLS
  Beta.ols = solve(t(X)%*%X)%*%t(X)%*%Y
  (Sbias.ols = norm(Beta.ols - W.true,type="2"))
  (Fbias.ols = norm(Beta.ols - W.true,type="F"))
  
  #################################################
  ####### Step 2: Estimate Pi
  #################################################
  rslt=ROC.lambda.mod(p,N,Y,Yhat.spr,nlambda,n,grp.info,ugrp,groupsparsity=tmpweight,A)
  
  return(rslt)
}

ROC.lambda.mod = function(p,N,Y,Yhat.spr,nlambda,n,grp.info,ugrp,groupsparsity,A){
  ## info about A: to evaluate Pi
  # A.max.1 = which(apply(A,1,max)==1)
  # M.ind = which(diag(A)!=1)
  # P.ind = intersect(M.ind,A.max.1)
  # C.ind = setdiff(M.ind,A.max.1)
  A.1to1 = which(apply(A,1,max)==1)
  A.1toM = which(apply(A,1,max)!=1)
  ## series of lambda
  lambda.all=c(-1,seq(0,2,length.out=nlambda));count=0
  adaptive.all=old.all=NULL
  for(lambda.cv in lambda.all){
    ### adaptive Pi
    Pi = fitpi_CV_mod(estPi="OLS",Ytrgt=Y,
                      Yhat=Yhat.spr,lambda.cv,n,grp.info,ugrp,
                      ridge.lambda=1e-10,groupsparsity)
    Pi.old = fitpi_CV(estPi="OLS",Ytrgt=Y,
                      Yhat=Yhat.spr,lambda.cv,n,grp.info,ugrp,
                      ridge.lambda=1e-10)
    
    evaluate.Pi=function(Pi){
      Pi.1to1 = which(apply(Pi,1,max)==1)
      Pi.1toM = which(apply(Pi,1,max)!=1)
      C.indmatch = intersect(Pi.1toM,A.1toM)
      (perc.C.indmatch=length(C.indmatch)/length(A.1toM))
      P.indmatch = intersect(Pi.1to1,A.1to1)
      (perc.P.indmatch=length(P.indmatch)/length(A.1to1))
      return(c(perc.C.indmatch=perc.C.indmatch,
               perc.P.indmatch=perc.P.indmatch,
               n.C.indmatch=length(C.indmatch),
               n.C.ind=length(A.1toM),
               n.P.indmatch=length(P.indmatch),
               n.P.ind=length(A.1to1)))
    }
    adaptive.all=rbind(adaptive.all,c(lambda.cv,evaluate.Pi(Pi)))
    old.all=rbind(old.all,c(lambda.cv,evaluate.Pi(Pi.old)))
    print(lambda.cv)
  }
  return(list(adaptive.all=adaptive.all,old.all=old.all))
}

fitpi_CV_mod = function(estPi,Ytrgt,Yhat,lambda.cv,n,grp.info,ugrp,ridge.lambda=1e-10,groupsparsity){
  Pi.all = Pi.all2 = Pi.all3 = list()
  tmpweight <- groupsparsity/max(groupsparsity)
  for(i in 1:n){
    ind = which(grp.info$g.index==ugrp[i])
    if(length(ind)>1){
      X.pi = t(Yhat[ind,])
      Y.pi = t(Ytrgt[ind,])
      Pi.ols = t( solve(t(X.pi)%*%X.pi+diag(ridge.lambda,ncol(X.pi)))%*%t(X.pi)%*%Y.pi  )
      Pi.perm = t(apply(Pi.ols,1,FUN=function(x){
        if((1-max(x/norm(x,type="2"))<(lambda.cv*tmpweight[i]))){
          x[which.max(x)]=1; x[-which.max(x)]=0
        }
        return(x)
      }))
    }else{
      Pi.perm = 1
    }
    Pi.all = append(Pi.all,list(Pi.perm))
  }
  Pi = bdiag(Pi.all)
  return(Pi)
}

fitpi_CV = function(estPi,Ytrgt,Yhat,lambda.cv,n,grp.info,ugrp,ridge.lambda=1e-10){
  Pi.all = Pi.all2 = Pi.all3 = list()
  for(i in 1:n){
    ind = which(grp.info$g.index==ugrp[i])
    if(length(ind)>1){
      X.pi = t(Yhat[ind,])
      Y.pi = t(Ytrgt[ind,])
      Pi.ols = t( solve(t(X.pi)%*%X.pi+diag(ridge.lambda,ncol(X.pi)))%*%t(X.pi)%*%Y.pi  )
      Pi.perm = t(apply(Pi.ols,1,FUN=function(x){
        if((1-max(x/norm(x,type="2"))<(lambda.cv))){
          x[which.max(x)]=1; x[-which.max(x)]=0
        }
        return(x)
      }))
    }else{
      Pi.perm = 1
    }
    Pi.all = append(Pi.all,list(Pi.perm))
  }
  Pi = bdiag(Pi.all)
  return(Pi)
}

gradient_update_nogrp <- function(X,Y,alpha=1,convergence=1e-4){
  p = ncol(X)
  n = nrow(X)
  
  oldW <- W <- matrix(0,p,p)
  gradient = t(as.matrix(X))%*%(as.matrix(Y))
  
  error <- 10000
  
  while(error >=convergence){
    
    # Gradient update
    W <- W + alpha*gradient
    
    # perform SVD approximation
    tmp <- svd(W)
    W <- tmp$u%*%t(tmp$v)
    
    error <- sqrt(sum((oldW-W)^2))
    oldW <- W
  }
  return(W)
  
}

#### args: select some N (numberofgroups) and variation (kappa)
args = commandArgs(TRUE)
seed <- as.numeric(args[[1]])
numberofgroups <- as.numeric(args[[2]]) #100~1500
kappa <- as.numeric(args[[3]])          #150
if(!"seed"%in%ls()){seed=1;numberofgroups=100;kappa=150}
#### load data
load(paste0(my.filepath,"summstat_final.RData"))
#### do one run
rslt <- run_uneven_mod(numberofgroups,seed,kappa)
save(rslt,file=
       paste0(my.filepath,"FIG_G2_",
              paste0("groups",numberofgroups,"kappa",kappa,"mod","seed",seed),".RData"))


#### make plot (after sufficient iteration of run_uneven_mod)
old.all.save = adaptive.all.save = 0
n.rep = 30
for(i in 1:n.rep){
  rslt <- run_uneven_mod(numberofgroups,seed=i,kappa)
  old.all.save = old.all.save + rslt$old.all
  adaptive.all.save = adaptive.all.save + rslt$adaptive.all
}
old.all=old.all.save/n.rep
adaptive.all=adaptive.all.save/n.rep
plot(old.all[,"perc.C.indmatch"],old.all[,"perc.P.indmatch"],
     type="l",xlim=c(0,1),ylim=c(0,1),cex=0.8,
     xlab="% correct selection of\none-to-many mapping",
     ylab="% correct selection of\none-to-one mapping")
lines(adaptive.all[,"perc.C.indmatch"],
      adaptive.all[,"perc.P.indmatch"],col="red",lty=6)
legend("bottomleft",bty="n",cex=0.8,
       lty=c(1,6),
       lwd=c(2,2),
       col=c("black","red"),
       legend=c("Equal threshold","Adaptive threshold")
)
