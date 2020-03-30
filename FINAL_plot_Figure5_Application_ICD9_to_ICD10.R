###################################################################
###################################################################
### Application (Figures 5)                                     ###
### Mapping From ICD-9 to ICD-10                                ###
### Contact: shixu@umich.edu                                    ###
###################################################################
###################################################################
rm(list=ls())
library(Matrix)
myfilepath=""
source(paste0(myfilepath,"FINAL_Application_functions.R"))
######################################
##             get data             ##
######################################
p=600 ##this is what is used in the paper
#### Load ICD10 embeddings
load(paste0(myfilepath,"vecs_VA18M_PheWASicdFALSE.RData"))
W <- rslt$fit$u[,1:p] %*% diag(sqrt(rslt$fit$d[1:p]))
X <- W[which(rowSums(W) != 0),]
row.names(X)=row.names(rslt$vecs)
#### Load ICD9 embeddings
load(paste0(myfilepath,"vecs_VA18M_PheWASicdTRUE.RData"))
W <- rslt$fit$u[,1:p] %*% diag(sqrt(rslt$fit$d[1:p]))
Z <- W[which(rowSums(W) != 0),]
row.names(Z)=row.names(rslt$vecs)
rm(rslt,W)
### save up just in case
X.save=X;Z.save=Z
### clean up code name
X = preprocess(X) ## ICD10
Z = preprocess(Z) ## ICD9
### remove zero columns
X.ind.zero=which(apply(abs(X),2,mean)<1e-10)
Z.ind.zero=which(apply(abs(Z),2,mean)<1e-10)
if(length(X.ind.zero)>0){X.nonzero=as.matrix(X[,-X.ind.zero])}else{X.nonzero=X}
if(length(Z.ind.zero)>0){Z.nonzero=as.matrix(Z[,-Z.ind.zero])}else{Z.nonzero=Z}
dim(X.nonzero);dim(Z.nonzero)

row.names(X)##icd10
row.names(Z)##icd9
Y.save=Y=Z#icd9
X.save=X=X#icd10
load(paste0(myfilepath,"Group_information.RData"))
grp.info=grp.info[grp.info$icd9cm%in%row.names(Y)&grp.info$icd10cm%in%row.names(X),]
#### select X and Y
X.save=X=as.matrix(X[row.names(X)%in%grp.info$icd10cm,])
Y.save=Y=as.matrix(Y[row.names(Y)%in%grp.info$icd9cm,])
which(!row.names(X.save)%in%grp.info$icd10cm)
which(!row.names(Y)%in%grp.info$icd9cm)
#### now match X and Y (should link by icd910 pair!)
### joint grp.info with (X,Y,grp.info)
ind.x=NULL
for(i in 1:nrow(grp.info)){
  ind.x=c(ind.x,which(row.names(X.save)==grp.info$icd10cm[i]))
}
ind.y=NULL
for(i in 1:nrow(grp.info)){
  ind.y=c(ind.y,which(row.names(Y.save)==grp.info$icd9cm[i]))
}
X=X.save[ind.x,];Y=Y.save[ind.y,]
c(dim(grp.info),dim(X),dim(Y))
### check if grp.info is ordered by grp, and X and Y match with grp.info
all(grp.info$g.index == cummax(grp.info$g.index))##monotonic g.index
which(row.names(X)!=grp.info$icd10cm)
which(row.names(Y)!=grp.info$icd9cm)
summary(1-grp.info$combination)
X.save=X;Y.save=Y
c(nrow(grp.info),nrow(X),nrow(Y))
#### remove data in X (and Y correspondingly) that
#### have grp.size>300
biggrps=names(which(table(grp.info$g.index)>300))
ind=which(grp.info$g.index%in%biggrps)
if(length(ind)>0){X=X[-ind,];Y=Y[-ind,];grp.info=grp.info[-ind,]}
#### save dimensions
N = nrow(grp.info)
ugrp = unique(grp.info$g.index);table(is.na(ugrp))
n = length(ugrp)
p=ncol(X)
c(N,nrow(X),nrow(Y),nrow(grp.info),nrow(X.save),nrow(Y.save),n) ##1463 groups
##### normalize X and Y
X.norms = rowNorms(as.matrix(X))
if(length(which(X.norms!=1))){
  X = X/rowNorms(as.matrix(X))
  Y = Y/rowNorms(as.matrix(Y))
}
summary(grp.info$combination)
mean(1-grp.info$combination)
######################################
##        iSphereMAP step I         ##
##       sperical regression        ##
######################################
### get data for fitting spherical regression
### estimate W using only the 1-1 mapping data
### i.e. do not use 1-m mapping data
ind1to1=which(row.names(Y)%in%grp.info$icd9cm[grp.info$combination==0])
Xfit=X[ind1to1,];Yfit=Y[ind1to1,]
c(nrow(Xfit),nrow(Xfit),nrow(grp.info),nrow(X),nrow(Y))
Beta.spr = gradient_update_nogrp(Xfit,Yfit,alpha=1,convergence=1e-10)
Yhat.spr = X%*%Beta.spr
norm(Yhat.spr-Y,"F")

######################################
##        iSphereMAP step II        ##
##             Mapping              ##
######################################
lambda.cv=0.6 ## this is what is used in the paper
Pi = fitpi_CV(estPi="OLS",Ytrgt=Y,Yhat=Yhat.spr,lambda.cv=lambda.cv,n,grp.info,ugrp)
Pi= Pi/rowNorms(Pi)#%*%X)
Pi.spike=apply(Pi,1,max)
perm.true=which(grp.info$combination==0)
comb.true=which(grp.info$combination==1)
perm.est=which(abs(Pi.spike-1)<=1e-10)
comb.est=which(abs(Pi.spike-1)>1e-10)
true.matched=grp.info[intersect(perm.est,perm.true),"icd9cm"]
matched=row.names(Y)[
  pmax(1,round(as.numeric(Pi[intersect(perm.est,perm.true),]%*%c(1:N)),0))
  ]
print(paste(c(
  lambda.cv,
  ### among all 1-1, how many correctly identified as 1-1
  length(which(perm.est%in%perm.true))/length(perm.true),
  ### among correctly identified as 1-1, how many are correct match
  table(matched==true.matched)["TRUE"]/length(true.matched),
  ### among all 1-m, how many correctly identified as 1-m
  length(which(comb.est%in%comb.true))/length(comb.true)
),collapse=","))

######################################
##        iSphereMAP step II        ##
##        Compare to Mikolov        ##
######################################
W.ols=solve(t(Xfit)%*%Xfit+diag(rep(1e-10,ncol(Xfit))))%*%t(Xfit)%*%Xfit
Yhat.ols = X%*%W.ols
grp.info2 = grp.info; grp.info2$g.index=1
Mikolov_w_GroupInfo = T
if(Mikolov_w_GroupInfo==T){
  Pi = fitpi_CV(estPi="cosine",Ytrgt=Y,Yhat=Yhat.ols,lambda.cv=0,
                n=length(unique(grp.info$g.index)),grp.info=grp.info,
                ugrp=unique(grp.info$g.index))
}else{
  Pi = fitpi_CV(estPi="cosine",Ytrgt=Y,Yhat=Yhat.ols,lambda.cv=0,
                n=1,grp.info=grp.info2,ugrp=1)
}
Pi.spike=apply(Pi,1,max)
perm.true=which(grp.info$combination==0)
comb.true=which(grp.info$combination==1)  
perm.est=which(abs(Pi.spike-1)<=1e-10)
comb.est=which(abs(Pi.spike-1)>1e-10)
true.matched=grp.info[intersect(perm.est,perm.true),"icd10cm"]
matched=colnames(Pi)[c(as.matrix(Pi[intersect(perm.est,perm.true),])%*%c(1:N))]
c(
  ### among all 1-1, how many correctly identified as 1-1
  length(which(perm.est%in%perm.true))/length(perm.true),
  ### among correctly identified as 1-1, how many are correct match
  table(matched==true.matched)["TRUE"]/length(true.matched),
  ### among all 1-m, how many correctly identified as 1-m
  length(which(comb.est%in%comb.true))/length(comb.true)
)


#####################################
#####################################
### Make plot (Figure 5)          ###
#####################################
#####################################
ind1to1=which(row.names(Y)%in%grp.info$icd9cm[grp.info$combination==0])
ind1tom=which(!row.names(Y)%in%grp.info$icd9cm[grp.info$combination==0])
ind=ind1to1 
Xfit=X[ind,];Yfit=Y[ind,]
c(nrow(Xfit),nrow(Xfit),nrow(grp.info),nrow(X),nrow(Y))
Beta.spr = gradient_update_nogrp(Xfit,Yfit,alpha=1,convergence=1e-10)
Yhat.spr = X%*%Beta.spr
mymethod="OLS"
lambda.cv=0.6
Pi = fitpi_CV(estPi=mymethod,Ytrgt=Y,Yhat=Yhat.spr,lambda.cv=lambda.cv,n,grp.info,ugrp)
Pi= Pi/rowNorms(Pi)
kk="714" ## choose between "714" and "E957"
icd9.select=paste0(kk,".")
t=grp.info[grep(icd9.select,grp.info$ICD9),]
ind=t$g.index 
(ind=which(grp.info$g.index%in%ind))
Pi.select=Pi[ind,ind]
Pi.true=diag(rep(1,nrow(Pi.select)));row.names(Pi.true)=row.names(Pi.select);colnames(Pi.true)=colnames(Pi.select)
image(Pi.select)
image(bdiag(Pi.true,Pi.true)[1:ncol(Pi.true),1:ncol(Pi.true)])
dim(Pi.select)

if(kk=="E957"){myh=4.1;myw=10}else{myh=5;myw=10}
pdf(file=paste0(myfilepath,"ICD910_RA2_fineRollup_p",p,"_",kk,mymethod,T,"0830",Sys.time(),".pdf"),h=myh,w=myw)
par(mfrow=c(2,1),mar=c(0.2,0,0,0),mai=c(0.2,0.1,0.1,0.1),omi=c(0.1,0,0,0))
if(kk=="714"){
  plotPi_714(Pi.plot=Pi.select)
  plotPi_714(Pi.plot=Pi.true)
}else if(kk=="E957"){
  plotPi_E957(Pi.plot=Pi.select)
  plotPi_E957(Pi.plot=Pi.true)
}else{
  plotPi(Pi.plot=Pi.select)
  plotPi(Pi.plot=Pi.true)
}
dev.off();
system("say done")
