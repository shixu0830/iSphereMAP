####################################################################################
######        Functions for simulation (Figures 2-3, F.1-F.4)                 ######
######                    Contact: shixu@umich.edu                            ######
####################################################################################

######################################
##        iSphereMAP step I         ##
##       sperical regression        ##
######################################
gradient_update_nogrp <- function(X,Y,alpha=1,convergence=1e-4){
  ## length normalize X and Y just in case
  X = X/apply(X,1,norm,type="2")
  Y = Y/apply(Y,1,norm,type="2")
  # Gradient update with stepsize alpha
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
  }## actually converges after one-step update
  return(W)
  
}

######################################
##       iSphereMAP step II         ##
######################################
fitpi_CV_application = function(estPi,Ytrgt,Yhat,lambda.cv,n,grp.info,ugrp){
  ## this function also incorporates MT method
  ## estPi="OLS" for iSphereMAP; estPi="cosine" for MT method
  Pi.all = Pi.all2 = Pi.all3 = list()
  ugrp=unique(grp.info$g.index)
  n=length(ugrp)
  for(i in 1:n){ ## estimation by group
    ind = which(grp.info$g.index==ugrp[i])
    if(length(ind)>1){
      X.pi = t(Yhat[ind,])
      Y.pi = t(Ytrgt[ind,])
      if(estPi=="OLS"){
        ## step 1: get initial pi
        Pi.initial = try(solve(t(X.pi)%*%X.pi) %*%t(X.pi)%*%Y.pi,silent=T)
        if(inherits(Pi.initial,"try-error")){ Pi.initial = ginv(t(X.pi)%*%X.pi) %*%t(X.pi)%*%Y.pi }
        ## step 2: hard-thresholding
        Pi.perm = t(apply(Pi.initial,1,FUN=function(x){
          xx = x/norm(x,type="2") ## xx = cosine(Pi_i, I_j)
          ind.max = which((1-xx) < lambda.cv)
          if(length(ind.max)>=1){
            #length(ind.max) never >1 if lambda.cv < 1-1/sqrt(2)
            x[which.max(x)]=1; x[-which.max(x)]=0
          }
          return(x)
        }))
      }else if(estPi=="cosine"){
        ## step 1: get initial pi
        Pi.initial = t( 
          ( t(X.pi)/rowNorms(t(X.pi)) )%*%t( t(Y.pi)/rowNorms(t(Y.pi)) )
        )
        ## step 2: hard-thresholding
        Pi.perm = t(apply(Pi.initial,1,FUN=function(x){
          x[which.max(x)]=1; x[-which.max(x)]=0
          return(x)
        }))
      }
      row.names(Pi.perm)=colnames(Y.pi)
      colnames(Pi.perm)=colnames(X.pi)
    }else{
      Pi.perm = 1
      Pi.perm=as.matrix(Pi.perm)
      row.names(Pi.perm)=row.names(Ytrgt)[ind]
      colnames(Pi.perm)=row.names(Yhat)[ind]
    }# end of "if(length(ind)>1)"
    
    if(prod(dim(Pi.perm))!=(length(ind)^2)){print("error");break}
    Pi.all = append(Pi.all,list(Pi.perm))
  }
  Pi = bdiag(Pi.all)
  row.names(Pi)=unlist(lapply(Pi.all,FUN=function(x){row.names(x)}))
  colnames(Pi)=unlist(lapply(Pi.all,FUN=function(x){colnames(x)}))
  return(Pi)
}
fitpi_CV = function(estPi,Ytrgt,Yhat,lambda.cv,n,grp.info,ugrp){
  ## this function also incorporates MT method
  ## estPi="OLS" for iSphereMAP; estPi="cosine" for MT method
  Pi.all = Pi.all2 = Pi.all3 = list()
  ugrp=unique(grp.info$g.index)
  n=length(ugrp)
  for(i in 1:n){ ## estimation by group
    ind = which(grp.info$g.index==ugrp[i])
    if(length(ind)>1){
      X.pi = t(Yhat[ind,])
      Y.pi = t(Ytrgt[ind,])
      if(estPi=="OLS"){
        ## step 1: get initial pi
        Pi.initial = try(solve(t(X.pi)%*%X.pi) %*%t(X.pi)%*%Y.pi,silent=T)
        if(inherits(Pi.initial,"try-error")){ Pi.initial = ginv(t(X.pi)%*%X.pi) %*%t(X.pi)%*%Y.pi }
        ## step 2: hard-thresholding
        Pi.perm = t(apply(Pi.initial,1,FUN=function(x){
          xx = x/norm(x,type="2") ## xx = cosine(Pi_i, I_j)
          ind.max = which((1-xx) < lambda.cv)
          if(length(ind.max)==1){
            #length(ind.max) never >1 if lambda.cv < 1-1/sqrt(2)
            x[which.max(x)]=1; x[-which.max(x)]=0
          }
          return(x)
        }))
      }else if(estPi=="cosine"){
        ## step 1: get initial pi
        Pi.initial = t( 
          ( t(X.pi)/rowNorms(t(X.pi)) )%*%t( t(Y.pi)/rowNorms(t(Y.pi)) )
        )
        ## step 2: hard-thresholding
        Pi.perm = t(apply(Pi.initial,1,FUN=function(x){
          x[which.max(x)]=1; x[-which.max(x)]=0
          return(x)
        }))
      }
      row.names(Pi.perm)=colnames(Y.pi)
      colnames(Pi.perm)=colnames(X.pi)
    }else{
      Pi.perm = 1
      Pi.perm=as.matrix(Pi.perm)
      row.names(Pi.perm)=row.names(Ytrgt)[ind]
      colnames(Pi.perm)=row.names(Yhat)[ind]
    }# end of "if(length(ind)>1)"
    
    if(prod(dim(Pi.perm))!=(length(ind)^2)){print("error");break}
    Pi.all = append(Pi.all,list(Pi.perm))
  }
  Pi = bdiag(Pi.all)
  row.names(Pi)=unlist(lapply(Pi.all,FUN=function(x){row.names(x)}))
  colnames(Pi)=unlist(lapply(Pi.all,FUN=function(x){colnames(x)}))
  return(Pi)
}
#####################################
## Cross-validation to find lambda ##
#####################################
find.lambda.cv = function(p,N,Y,Yhat.spr,nlambda,n,grp.info,ugrp){
  #Create n.folds equally size folds: save and use for each lambda
  nfolds=2; randomorder=sample(1:p,size=p,replace=F)
  folds <- cut(randomorder,breaks=nfolds,labels=FALSE)
  #Perform cross validation
  lambda.all=seq(1e-5,1-1e-5,length.out=nlambda);count=0
  # lambda.all=seq(0,1-1/sqrt(2),length.out=nlambda);count=0
  cv.err.all=cv.perm.rows.all=cv.match.rows.all=
    cv.err.Yfit.all=cv.err.Y.all=
    matrix(0,nrow=length(lambda.all),ncol=nfolds)
  for(lambda.cv in lambda.all){
    #print(lambda.cv)
    count=count+1
    for(i in 1:nfolds){
      #Segement data by fold
      testIndexes <- which(folds==i,arr.ind=TRUE)
      Pi.i = fitpi_CV(estPi="OLS",Ytrgt=Y[,testIndexes],
                      Yhat=Yhat.spr[,testIndexes],lambda.cv,n,grp.info,ugrp)
      Yfit=Pi.i%*%Yhat.spr; Yfit=Yfit/rowNorms(Yfit) 
      ### note: we normalize full row, but later error only on half columns
      
      ### compute prediction error (F norm)
      cv.err.i = norm(Y[,-testIndexes]-Yfit[,-testIndexes],type="F")
      cv.err.all[count,i] = cv.err.i
      cv.err.Yfit.all[count,i] = norm(Yfit[,-testIndexes],type="F")
      cv.err.Y.all[count,i] = norm(Y[,-testIndexes],type="F")
      ### compute number of permutation rows
      ind.onehot = which(apply(Pi.i,1,max)==1)
      ind.perm = which(apply(Pi.i-diag(rep(1,N)),1,max)==1)
      cv.perm.rows.all[count,i] = length(ind.perm)
      cv.match.rows.all[count,i] = length(ind.onehot)
    }
  }
  cv.err=apply(cv.err.all,1,sum)
  cv.err.Yfit=apply(cv.err.Yfit.all,1,sum)
  cv.err.Y=apply(cv.err.Y.all,1,sum)
  cv.perm.rows=apply(cv.perm.rows.all,1,mean)
  cv.match.rows=apply(cv.match.rows.all,1,mean)
  return(list(cv.err=cv.err,
              cv.err.Yfit=cv.err.Yfit,
              cv.err.Y=cv.err.Y,
              cv.perm.rows=cv.perm.rows,
              cv.match.rows=cv.match.rows,
              lambda.all=lambda.all))
}

######################################
##       Simulation one run         ##
######################################
run=function(kk,tt,myseed,MVN,add_mismatch_1to1,perm_only,wrong_grp_info){
  ####################################################################################
  ####################################################################################
  ####### generate data
  ####################################################################################
  ####################################################################################
  set.seed(myseed)
  ##### first simulate X
  if(MVN==T){##generate X and Y that follow MVN then normalize
    X = mvrnorm(n=N,mu=rep(0,p),Sigma=diag(rep(1,p)))
    X = X/apply(X,1,norm,"2")
  }else{##generate X and Y that follow vMF
    grpcenters = mvrnorm(n=n,mu=rep(0,p),Sigma=diag(rep(1,p)))
    X=NULL
    for(i in 1:n){
      grpsize = length(which(grp.info$g.index==ugrp[i]))
      myalpha = rep(1,n); myalpha[i]=2 ## upweight i-th grp in the mixture
      dat.grp = rmovMF(n=grpsize,theta=mykappa*grpcenters/apply(grpcenters,1,norm,"2"),alpha=myalpha)
      X = rbind(X,dat.grp)
    }
  }## end of "if(MVN==T)"
  ##### specify true Pi, called "A": sparse shuffle within group 
  A = diag(1,N,N)
  A = shuffle_A(A, grp.info,sparsity,n) ## edit A to include 1-1 and/or 1-m mismatch
  S = length(which(diag(A)!=1)) ## how many mismatched rows
  ##### second simulate Y
  Y.mean = A%*%X%*%W.true
  if(MVN==T){##generate Y that follow MVN then normalize
    Y = Y.mean*mykappa/p+mvrnorm(n=N,mu=rep(0,p),Sigma=diag(rep(1,p)/p))
    Y=Y/apply(Y,1,norm,"2")
  }else{##generate Y that follow vMF
    Y=t(apply(Y.mean/rowNorms(Y.mean),1,FUN=function(y){
      rmovMF(n=1,theta=mykappa*y)
    }))
  }
  ####################################################################################
  ####################################################################################
  ####### Analysis: steps I to III
  ####################################################################################
  ####################################################################################
  if(wrong_grp_info==T){## misspecified group structure
    dd=5 ## misspecification of group structure for every dd groups
    grp.info$g.index[grp.info$g.index%%dd==0]=
      grp.info$g.index[grp.info$g.index%%dd==0]-1
    ugrp=unique(grp.info$g.index)
  }
  ####################################################################################
  ####################################################################################
  ####### Step I: estimate W-one-hat [note: norm(W.true,type="F")=norm(Beta.spr,type="F")=sqrt(p)]
  ####################################################################################
  ####################################################################################
  Beta.spr = gradient_update_nogrp(X,Y,alpha=1,convergence=1e-10)
  Beta.ols = try(solve(t(X)%*%X)%*%t(X)%*%Y,silent=T)
  if(inherits(Beta.ols,"try-error")){Beta.ols = ginv(t(X)%*%X)%*%t(X)%*%Y}
  Yhat.spr = X%*%Beta.spr
  Yhat.ols = X%*%Beta.ols
  ### save result
  (Fbias.spr = norm(Beta.spr - W.true,type="F"))##influenced by kappa and A
  (Fbias.ols = norm(Beta.ols - W.true,type="F"))
  W.rslt = c(Truenorm = norm(W.true,type="F"),estnorm = norm(Beta.spr,type="F"))
  ####################################################################################
  ####################################################################################
  ####### Step II: est Pi-two-hat
  ####################################################################################
  ####################################################################################
  ##### (1) iSphereMAP
  ## cross validation to find best lambda
  cv.rslt=find.lambda.cv(p,N,Y,Yhat.spr,nlambda,
                         n=length(ugrp),grp.info=grp.info,ugrp=ugrp)
  cv.err=cv.rslt$cv.err;#plot(cv.err)
  lambda.all=cv.rslt$lambda.all
  lambda.cv = lambda.all[which.min(cv.err)]
  Pi = fitpi_CV(estPi="OLS",Ytrgt=Y,Yhat=Yhat.spr,lambda.cv,
                n=length(ugrp),grp.info=grp.info,ugrp=ugrp)
  ##### (2) Machine translation (MT) method by Mikolov
  Pi.MT = fitpi_CV(estPi="cosine",Ytrgt=Y,Yhat=Yhat.ols,lambda.cv=1e10,n=1,grp.info=grp.info2,ugrp=1)
  ### NB: neither the saved A nor Pi are re-normalized
  Pirslt = list(S=S,A=A,Pi=Pi,Pi.MT=Pi.MT,cv.rslt=cv.rslt)
  ####################################################################################
  ####################################################################################
  ####### Step III: est W-hat-two
  ####################################################################################
  ####################################################################################
  #### get clean data for W refinement
  if(add_mismatch_1to1==T){
    ind=which(apply(Pi/rowNorms(Pi),2,max)==1);length(ind)
  }else{ #add_mismatch_1to1==F
    ind=which(diag(Pi/rowNorms(Pi))==1);length(ind)
  }
  n.update=length(ind) #save sample size for W refinement
  X.match=(Pi%*%X)[ind,];Y.match=Y[ind,]
  X.match=X.match/rowNorms(X.match) ## no need to do length normalization for Y.match
  index.matched=grp.info$g.index[ind]
  ugrp.matched=unique(index.matched)
  ###compute W using spherical regress and OLS
  Beta.spr = gradient_update_nogrp(X.match,Y.match,alpha=1,convergence=1e-10)
  Beta.ols = try(solve(t(X.match)%*%X.match)%*%t(X.match)%*%Y.match,silent=T)
  if(inherits(Beta.ols,"try-error")){Beta.ols = ginv(t(X.match)%*%X.match)%*%t(X.match)%*%Y.match}
  Yhat.spr=X%*%Beta.spr
  Yhat.ols=X%*%Beta.ols
  ###compare spr and ols
  (Fbias.spr.update = norm(Beta.spr - W.true,type="F")) ##influenced by kappa and A
  (Fbias.ols.update = norm(Beta.ols - W.true,type="F"))
  W.rslt.update = c(Truenorm.update = norm(W.true,type="F"),estnorm.update = norm(Beta.spr,type="F"))
  ### save all results for W estimation: step I and step III
  W.rslt.all = c(S=S,N=N,n=n,sparsity=sparsity,mykappa=mykappa,
                 Fbias.spr=Fbias.spr,Fbias.ols=Fbias.ols,
                 Fbias.spr.update=Fbias.spr.update,Fbias.ols.update=Fbias.ols.update,
                 W.rslt,W.rslt.update,n.update=n.update)
  
  if(add_mismatch_1to1==F&perm_only==F&wrong_grp_info==F){
    ####################################################################################
    ####################################################################################
    ####### Step IV: est Pi-hat-three (use when theoretical setting is violated)
    ####################################################################################
    ####################################################################################
    ##### (1) iSphereMAP
    ## cross validation to find best lambda
    cv.rslt=find.lambda.cv(p,N,Y,Yhat.spr,nlambda,
                           n=length(ugrp),grp.info=grp.info,ugrp=ugrp)
    cv.err=cv.rslt$cv.err
    lambda.all=cv.rslt$lambda.all
    lambda.cv = lambda.all[which.min(cv.err)]
    Pi = fitpi_CV(estPi="OLS",Ytrgt=Y,Yhat=Yhat.spr,lambda.cv,
                  n=length(ugrp),grp.info=grp.info,ugrp=ugrp)
    ##### (2) Machine translation (MT) method by Mikolov
    Pi.MT = fitpi_CV(estPi="cosine",Ytrgt=Y,Yhat=Yhat.ols,lambda.cv=1e10,n=1,grp.info=grp.info2,ugrp=1)
    ### NB: neither the saved A nor Pi are re-normalized
    Pirslt.update = list(S=S,A=A,Pi=Pi,Pi.MT=Pi.MT,cv.rslt=cv.rslt)
  }else{#if(!(add_mismatch_1to1==F&perm_only==F&wrong_grp_info==F))
    Pirslt.update = NULL
  }
  
  #### all results from simulation
  Pierr=readinPi(Pirslt)
  if(!is.null(Pirslt.update)){Pierr.update=readinPi(Pirslt.update)}else{Pierr.update=NULL}
  all_results = list(Wrslt=W.rslt.all,Pierr=Pierr,Pierr.update=Pierr.update)
  return(all_results)
}


###########################
##     miscellaneous     ##
###########################
##rowNorms() instead of using library(wordspace)
rowNorms = function(M){
  return(apply(M,1,norm,type="2"))
}

### Generate true Pi
shuffle_A = function(A, grp.info,sparsity,n){
  ### edit rows of A to either 1-1 mismatch or 1-m weight vector
  if(perm_only==T){
    ## edit A into a permutation matrix
    d=2
    for(i in 1:n){
      ind = which(grp.info$g.index==ugrp[i])
      n.shuffle = floor(length(ind)*sparsity)
      if(n.shuffle>0){
        ind.shuffle = sample(ind,size=n.shuffle*d,replace=F)
        for(j in 1:n.shuffle){
          A[ind.shuffle[j],ind.shuffle[n.shuffle*d+1-j]]=1
          A[ind.shuffle[n.shuffle*d+1-j],ind.shuffle[j]]=1
          A[ind.shuffle[j],ind.shuffle[j]]=0
          A[ind.shuffle[n.shuffle*d+1-j],ind.shuffle[n.shuffle*d+1-j]]=0
        }
      }
    }
  }else{#perm_only==F
    ## edit A to include both 1-1 and 1-m mismatch
    d=4
    for(i in 1:n){
      ind = which(grp.info$g.index==ugrp[i])
      n.shuffle = floor(length(ind)*sparsity/(d/2))
      if(n.shuffle>0){
        ind.shuffle = sample(ind,size=min(length(ind),n.shuffle*d),replace=F)
        for(j in 1:n.shuffle){
          ### for n.shuffle rows, half of them are 1-1 mismatch
          A[ind.shuffle[j],ind.shuffle[n.shuffle*d+1-j]]=1
          A[ind.shuffle[n.shuffle*d+1-j],ind.shuffle[j]]=1
          A[ind.shuffle[j],ind.shuffle[j]]=0
          A[ind.shuffle[n.shuffle*d+1-j],ind.shuffle[n.shuffle*d+1-j]]=0
          ### the other half are 1-m mismatch
          w1=runif(n=length(ind)); w1=w1/norm(w1,type="2")
          A[ind.shuffle[n.shuffle+j],ind]=w1
          w2=runif(n=length(ind)); w2=w2/norm(w2,type="2")
          A[ind.shuffle[n.shuffle*2+j],ind]=w2
        }
      }
    }
  }## end of "if(perm_only==T)"
  return(A)
}

