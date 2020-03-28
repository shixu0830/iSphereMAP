##############################################
## clean ICD9 codes                         ##
##############################################
cleancode = function(code){
  t = unlist(lapply(strsplit(code,"[.]"),FUN=function(x){
    if(nchar(x[1])<3 & length(grep("[0-9]",substr(x[1],start=1,stop=1)))>0){
      ## x[1] is a number, but length is smaller than 3
      x[1] = paste0(paste(rep("0",3-nchar(x[1])),collapse=""), x[1])}
    if(nchar(x[2])<2 & length(grep("[0-9]",substr(x[2],start=1,stop=1)))>0){
      ## x[2] is a number, but length is smaller than 2
      x[2] = paste0(x[2], paste(rep("0",2-nchar(x[2])),collapse=""))}
    x = paste(x,collapse="")
    if(nchar(x)<5 & length(grep("[0-9]",substr(x[1],start=1,stop=1)))>0){
      x = paste0(x,paste(rep("0",5-nchar(x)),collapse=""))}
    return(x)
  }))
  return(t)
}

##############################################
## pre-process code vec: clean up code name ##
##############################################
preprocess = function(dat){
  ### extract icd
  if(class(dat[,1])!="numeric"){icd=as.character(dat[,1]); dat=dat[,-1]
  }else{icd=row.names(dat)}
  ### clean up code to without [.]
  codes=unlist(lapply(strsplit(as.character(icd),"[.]"),FUN=function(x){paste(x,collapse="")}))
  ### remove duplicated code rows
  ind.dup=which(duplicated(codes))
  if(length(ind.dup)>0){
    codes=codes[-ind.dup]
    dat=dat[-ind.dup,]
  }
  ### rename rownames
  row.names(dat) = codes
  return(dat)
}

##############################################
## Step I: sperical regression              ##
##############################################
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
##############################################
## Step II: estimate PI                     ##
##############################################
### define mapping function
fitpi_CV = function(estPi,Ytrgt,Yhat,lambda.cv,n,grp.info,ugrp){
  ## redo row-normalization just in case
  Yhat = Yhat/rowNorms(as.matrix(Yhat))
  Ytrgt = Ytrgt/rowNorms(as.matrix(Ytrgt))
  ## redo compute ugrp and n just in case
  ugrp=unique(grp.info$g.index)
  n=length(ugrp)
  Pi.all = Pi.all2 = Pi.all3 = list()
  ugrp=unique(grp.info$g.index) ##i=994/1004 for Phe714
  n=length(ugrp)
  for(i in 1:n){
    ind = which(grp.info$g.index==ugrp[i])
    if(length(ind)>1){
      X.pi = t(Yhat[ind,])
      Y.pi = t(Ytrgt[ind,])
      if(estPi=="OLS"){
        Pi.ols = t( solve(t(X.pi)%*%X.pi+diag(1e-10,length(ind)))%*%t(X.pi)%*%Y.pi  )
      }else if(estPi=="cosine"){
        Pi.ols = t( 
          ( t(X.pi)/rowNorms(t(X.pi)) )%*%t( t(Y.pi)/rowNorms(t(Y.pi)) )
        )
      }else if(estPi=="spherical"){
        Pi.ols = gradient_update_nogrp(X.pi,Y.pi,alpha=1,convergence=1e-10)
      }
      Pi.perm = t(apply(Pi.ols,1,FUN=function(x){
        if(estPi=="cosine"){
          x[which.max(x)]=1; x[-which.max(x)]=0
        }else{
          xx = x/norm(x,type="2")
          ind.max = which((xx)>lambda.cv)
          if(length(ind.max)>=1){  #>1
            x[which.max(x)]=1; x[-which.max(x)]=0
          }
        }
        return(x)
      }))
      row.names(Pi.perm)=colnames(Y.pi)
      colnames(Pi.perm)=colnames(X.pi)
    }else{#e.g.  i=41
      Pi.perm = 1
      Pi.perm=as.matrix(Pi.perm)
      row.names(Pi.perm)=row.names(Ytrgt)[ind]
      colnames(Pi.perm)=row.names(Yhat)[ind]
    }
    if(prod(dim(Pi.perm))!=(length(ind)^2)){print("error");break}
    Pi.all = append(Pi.all,list(Pi.perm))
  }
  Pi = bdiag(Pi.all)
  # row.names(Pi)=row.names(Ytrgt); colnames(Pi)=row.names(Yhat)
  row.names(Pi)=unlist(lapply(Pi.all,FUN=function(x){row.names(x)}))
  colnames(Pi)=unlist(lapply(Pi.all,FUN=function(x){colnames(x)}))
  return(Pi)
}




rowNorms = function(M){
  return(apply(M,1,FUN=function(x){norm(x,type="2")}))
}

##############################################
## Code to make plot (Figure 5)             ##
##############################################
plotPi=function(Pi.plot=Pi.select){
  ind.cc=which(grp.info$g.index%in%grp.info[grep(icd9.select,grp.info$ICD9),]$g.index)
  grp.info.cc=grp.info[ind.cc,]
  cc=grp.info.cc$icd9cm
  which(cc!=row.names(Pi.plot))##QC
  print(length(cc)==nrow(grp.info.cc))##QC
  mycex1=0.75;mycex2=0.65
  
  dic=grp.info[,c("ICD9","ICD9.String","icd10cm")]
  L = matrix(NA,nrow=length(cc),ncol=5)
  row.names(L) = paste(grp.info.cc[,1],grp.info.cc[,5])
  R=L
  R.loc = data.frame(R); L.loc = data.frame(L);  
  L.loc[,1]=1; R.loc[,1]=3; ## x-axis
  L.loc[,4]=grp.info.cc[,"ICD9"];R.loc[,4]=grp.info.cc[,"icd10cm"]
  L.loc[,3]=grp.info.cc[,"icd9cm"]; R.loc[,3]=grp.info.cc[,"icd10cm"]
  L.loc[,4]=as.character(L.loc[,4]);L.loc[,3]=as.character(L.loc[,3])
  R.loc[,4]=as.character(R.loc[,4]);R.loc[,3]=as.character(R.loc[,3])
  L.loc[,5]=dic[match(L.loc[,4],dic[,1]),2];L.loc[,5]=as.character(L.loc[,5]);
  ###load ICD10 description
  icd10.desc=read.csv("icd10cm_codes_2018.csv",header=F)
  R.loc[,5]=icd10.desc[match(R.loc[,4],icd10.desc[,1]),2];R.loc[,5]=as.character(R.loc[,5]);
  R.loc[,5]=gsub("\\n", "\n", R.loc[,5], fixed=TRUE)
  R.loc[row.names(R)=="",1]=NA
  # L.loc=L.loc[order(as.character(L.loc[,4]),decreasing=F),]
  # R.loc=R.loc[order(as.character(L.loc[,4]),decreasing=F),]
  R.loc[,2]=seq(0,length(R),length.out=nrow(R))
  L.loc[,2]=seq(0,length(R),length.out=nrow(L))
  ### delete duplicated rows in R and L (make 1-m & m-1 mapping prettier)
  R.loc=R.loc[!duplicated(R.loc[,3]),]
  if(icd9.select=="714."){
    L.loc=L.loc[-c(3,5),] 
  }else{
    L.loc=L.loc[!duplicated(L.loc[,3]),]
  }
  
  
  if(icd9.select=="E957."){
    L.loc[,2]=L.loc[,2]
    R.loc[,2]=c(L.loc[1:3,2],16.875-0.65,22.5-5+0.25,28.125-5.625)-5.625
  }
  
  dl=-1.5;dr=-1.5
  # dl=0;dr=0
  plot(
    x=L.loc[,1]+dl,y=L.loc[,2],
    xlim=c(0-10,4+10),ylim=range(c(L.loc[,2],R.loc[,2])),
    type="p",pch=16,
    xlab="",ylab="",axes=F
  )
  text(x=0.5+0.25+dl,y=L.loc[,2],srt=0,adj=1,font=2,
       labels=L.loc[,4],xpd=TRUE,cex=mycex1,lwd=1)
  text(x=0.5-1+dl,y=L.loc[,2],srt=0,adj=1,
       labels=L.loc[,5],
       xpd=TRUE,cex=mycex2,lwd=1)
  
  points(
    x=R.loc[,1]+dr,y=R.loc[,2],
    pch=20
  )
  text(x=3.5-0.25+dr,y=R.loc[,2],srt=0,adj=0,font=2,
       labels=R.loc[,4],xpd=TRUE,cex=mycex1,lwd=1)
  text(x=3.5+1+dr,y=R.loc[,2],srt=0,adj=0,
       labels=R.loc[,5],
       xpd=TRUE,cex=mycex2,lwd=1)
  ###now delete duplicated leftcodes in cc
  ind.dup.cc=which(duplicated(grp.info.cc$icd9cm))
  if(length(ind.dup.cc)>0){
    for(i in unique(cc[ind.dup.cc])){
      mean.weights=apply(Pi.plot[row.names(Pi.plot)==i,],2,mean)
      mean.weights=mean.weights/norm(mean.weights,"2")
      for(j in 1:length(which(row.names(Pi.plot)==i))){
        Pi.plot[row.names(Pi.plot)==i,][j,]=mean.weights
      }
    }
    Pi.plot=Pi.plot[-ind.dup.cc,]
    grp.info.cc=grp.info.cc[-ind.dup.cc,]
    cc=cc[-ind.dup.cc]
  }
  
  ###plot grey lines
  for(i in 1:length(cc)){
    rightweights = Pi.plot[i,]
    rightweights = rightweights[rightweights>0]
    if(anyDuplicated(names(rightweights))>0){
      dup.name=names(rightweights)[duplicated(names(rightweights))]
      dup.ind=which(names(rightweights)%in%dup.name)
      dup.max.ind=which.max(rightweights[dup.ind])
      rightweights = rightweights[-dup.ind[-dup.max.ind]]
    }
    rightcodes = names(rightweights)[rightweights>0]
    rightcodes.ind = sapply(rightcodes,FUN=function(x){which(R.loc[,3]%in%x)})
    leftcodes.ind = which(L.loc[,3]==cc[i])
    for(j in 1:length(rightcodes)){
      segments(x0=L.loc[leftcodes.ind,1]+dl,y0=L.loc[leftcodes.ind,2],
               x1=R.loc[rightcodes.ind[j],1]+dr,y1=R.loc[rightcodes.ind[j],2],
               lwd = as.numeric(rightweights[j])*5,
               col="#00000050")
    }
  }
}

plotPi_E957=function(Pi.plot=Pi.select){
  ind.cc=which(grp.info$g.index%in%grp.info[grep(icd9.select,grp.info$ICD9),]$g.index)
  grp.info.cc=grp.info[ind.cc,]
  cc=grp.info.cc$icd9cm
  which(cc!=row.names(Pi.plot))##QC
  print(length(cc)==nrow(grp.info.cc))##QC
  #mycex1=0.75;mycex2=0.65
  mycex1=1;mycex2=0.9
  dic=matrix(c(
    "E957.9","SSI by jumping from unspecified site",
    "E957.2","SSI by jumping from natural sites",
    "E957.1","SSI by jumping from other man-made structures",
    "E957.0","SSI by jumping from residential premises",
    "Y929","Unspecified place or not applicable",
    "Y92838","Other recreation area",
    "Y92828","Other wilderness area",
    "Y9289", "Other specified places",
    "Y92009","Unspecified place in unspecified\n non-institutional (private) residence",
    "X80XXXA","Intentional self-harm by jumping\n from a high place, initial encounter"
  ),ncol=2,byrow=T)
  L = matrix(NA,nrow=length(cc),ncol=5)
  row.names(L) = paste(grp.info.cc[,1],grp.info.cc[,5])
  R=L
  R.loc = data.frame(R); L.loc = data.frame(L);  
  L.loc[,1]=1; R.loc[,1]=3; ## x-axis
  L.loc[,4]=grp.info.cc[,"ICD9"];R.loc[,4]=grp.info.cc[,"icd10cm"]
  L.loc[,3]=grp.info.cc[,"icd9cm"]; R.loc[,3]=grp.info.cc[,"icd10cm"]
  L.loc[,4]=as.character(L.loc[,4]);L.loc[,3]=as.character(L.loc[,3])
  R.loc[,4]=as.character(R.loc[,4]);R.loc[,3]=as.character(R.loc[,3])
  L.loc[,5]=dic[match(L.loc[,4],dic[,1]),2];L.loc[,5]=as.character(L.loc[,5]);
  R.loc[,5]=dic[match(R.loc[,4],dic[,1]),2];R.loc[,5]=as.character(R.loc[,5]);
  
  R.loc[,2]=seq(0,length(R),length.out=nrow(R))
  L.loc[,2]=seq(0,length(R),length.out=nrow(L))
  ### delete duplicated rows in R and L (make 1-m & m-1 mapping prettier)
  R.loc=R.loc[!duplicated(R.loc[,3]),]
  if(icd9.select=="714."){
    L.loc=L.loc[-c(3,5),] 
  }else{
    L.loc=L.loc[!duplicated(L.loc[,3]),]
  }
  
  
  if(icd9.select=="E957."){
    L.loc[,2]=L.loc[,2]
    R.loc[,2]=c(L.loc[1:3,2],16.875-1.65,22.5-4+0.25,28.125-5.625)-5.625
  }
  
  dl=-1.5;dr=-1.5
  plot(
    x=L.loc[,1]+dl,y=L.loc[,2],
    xlim=c(0-10,4+10),ylim=range(c(L.loc[,2],R.loc[,2])),
    type="p",pch=16,
    xlab="",ylab="",axes=F
  )
  text(x=0.5+0.35+dl,y=L.loc[,2],srt=0,adj=1,font=2,
       labels=L.loc[,4],xpd=TRUE,cex=mycex1,lwd=1)
  text(x=0.5-1.05+dl,y=L.loc[,2],srt=0,adj=1,
       labels=L.loc[,5],
       xpd=TRUE,cex=mycex2,lwd=1)
  
  points(
    x=R.loc[,1]+dr,y=R.loc[,2],
    pch=20
  )
  text(x=3.5-0.35+dr,y=R.loc[,2],srt=0,adj=0,font=2,
       labels=R.loc[,4],xpd=TRUE,cex=mycex1,lwd=1)
  text(x=3.5+1.6+dr,y=R.loc[,2],srt=0,adj=0,
       labels=R.loc[,5],
       xpd=TRUE,cex=mycex2,lwd=1)
  ###now delete duplicated leftcodes in cc
  ind.dup.cc=which(duplicated(grp.info.cc$icd9cm))
  if(length(ind.dup.cc)>0){
    for(i in unique(cc[ind.dup.cc])){
      mean.weights=apply(Pi.plot[row.names(Pi.plot)==i,],2,mean)
      mean.weights=mean.weights/norm(mean.weights,"2")
      for(j in 1:length(which(row.names(Pi.plot)==i))){
        Pi.plot[row.names(Pi.plot)==i,][j,]=mean.weights
      }
    }
    Pi.plot=Pi.plot[-ind.dup.cc,]
    grp.info.cc=grp.info.cc[-ind.dup.cc,]
    cc=cc[-ind.dup.cc]
  }
  
  ###plot grey lines
  for(i in 1:length(cc)){
    rightweights = Pi.plot[i,]
    rightweights = rightweights[rightweights>0]
    if(length(rightweights)>0){
      if(anyDuplicated(names(rightweights))>0){
        dup.name=names(rightweights)[duplicated(names(rightweights))]
        dup.ind=which(names(rightweights)%in%dup.name)
        dup.max.ind=which.max(rightweights[dup.ind])
        rightweights = rightweights[-dup.ind[-dup.max.ind]]
      }
      rightcodes = names(rightweights)[rightweights>0]
      rightcodes.ind = sapply(rightcodes,FUN=function(x){which(R.loc[,3]%in%x)})
      leftcodes.ind = which(L.loc[,3]==cc[i])
      for(j in 1:length(rightcodes)){
        segments(x0=L.loc[leftcodes.ind,1]+dl,y0=L.loc[leftcodes.ind,2],
                 x1=R.loc[rightcodes.ind[j],1]+dr,y1=R.loc[rightcodes.ind[j],2],
                 lwd = as.numeric(rightweights[j])*5,
                 col="#00000050")
      }
    }
  }
}

plotPi_714=function(Pi.plot=Pi.select){
  #image(Pi.plot)
  ind.cc=which(grp.info$g.index%in%grp.info[grep(icd9.select,grp.info$ICD9),]$g.index)
  grp.info.cc=grp.info[ind.cc,]
  cc=grp.info.cc$icd9cm
  which(cc!=row.names(Pi.plot))##QC
  print(length(cc)==nrow(grp.info.cc))##QC
  #mycex1=0.75;mycex2=0.65
  mycex1=1;mycex2=0.9
  
  dic=matrix(c(
    "M069","RA, unspecified",
    "M0500","Felty's syndrome, unspecified site",
    "M0530","Rheumatoid heart disease with RA of unspecified site",
    "M0560","RA of unspecified site with involvement of other organs and systems",
    "M061","Adult-onset Still's disease",
    "M0510","Rheumatoid lung disease with RA of unspecified size",
    "M0800","Unspecified juvenile RA of unspecified site",
    "M083","Juvenile rheumatoid polyarthritis (seronegative)",
    "M0840","Pauciarticular juvenile RA, unspecified site",
    "M1200","Chronic postrheumatic arthropathy, unspecified site",
    "M064","Inflammatory polyarthropathy",
    "714.0","RA",
    "714.1","Felty's syndrome",
    "714.2","Other RA with visceral or systemic involvement",
    "714.81","Rheumatoid lung",
    "714.30","Polyarticular juvenile RA, chronic or unspecified",
    "714.31","Polyarticular juvenile RA, acute",
    "714.32","Pauciarticular juvenile RA",
    "714.33","Monoarticular juvenile RA",
    "714.4","Chronic postrheumatic arthropathy",
    "714.89","Other specified inflammatory polyarthropathies",
    "714.9","Unspecified inflammatory polyarthropathy"
  ),ncol=2,byrow=T)
  L = matrix(NA,nrow=length(cc),ncol=5)
  row.names(L) = paste(grp.info.cc[,1],grp.info.cc[,5])
  R=L
  R.loc = data.frame(R); L.loc = data.frame(L);  
  L.loc[,1]=1; R.loc[,1]=3; ## x-axis
  L.loc[,4]=grp.info.cc[,"ICD9"];R.loc[,4]=grp.info.cc[,"icd10cm"]
  L.loc[,3]=grp.info.cc[,"icd9cm"]; R.loc[,3]=grp.info.cc[,"icd10cm"]
  L.loc[,4]=as.character(L.loc[,4]);L.loc[,3]=as.character(L.loc[,3])
  R.loc[,4]=as.character(R.loc[,4]);R.loc[,3]=as.character(R.loc[,3])
  L.loc[,5]=dic[match(L.loc[,4],dic[,1]),2];L.loc[,5]=as.character(L.loc[,5]);
  R.loc[,5]=dic[match(R.loc[,4],dic[,1]),2]
  R.loc[,5]=as.character(R.loc[,5]);
  R.loc[,5]=gsub("\\n", "\n", R.loc[,5], fixed=TRUE)
  R.loc[row.names(R)=="",1]=NA
  R.loc[,2]=seq(0,length(R),length.out=nrow(R))
  L.loc[,2]=seq(0,length(R),length.out=nrow(L))
  ### delete duplicated rows in R and L (make 1-m & m-1 mapping prettier)
  R.loc=R.loc[!duplicated(R.loc[,3]),]
  if(icd9.select=="714."){
    L.loc=L.loc[-c(3,5),] 
  }else{
    L.loc=L.loc[!duplicated(L.loc[,3]),]
  }
  
  
  if(icd9.select=="E957."){
    L.loc[,2]=L.loc[,2]
    R.loc[,2]=c(L.loc[1:3,2],16.875-0.65,22.5-5+0.25,28.125-5.625)-5.625
  }
  
  dl=-1.5;dr=-1.5
  # dl=0;dr=0
  plot(
    x=L.loc[,1]+dl,y=L.loc[,2],
    xlim=c(0-10,4+10),ylim=range(c(L.loc[,2],R.loc[,2])),
    type="p",pch=16,
    xlab="",ylab="",axes=F
  )
  text(x=0.5+0.3+dl,y=L.loc[,2],srt=0,adj=1,font=2,
       labels=L.loc[,4],xpd=TRUE,cex=mycex1,lwd=1)
  text(x=0.5-1+dl,y=L.loc[,2],srt=0,adj=1,
       labels=L.loc[,5],
       xpd=TRUE,cex=mycex2,lwd=1)
  
  points(
    x=R.loc[,1]+dr,y=R.loc[,2],
    pch=20
  )
  text(x=3.5-0.4+dr,y=R.loc[,2],srt=0,adj=0,font=2,
       labels=R.loc[,4],xpd=TRUE,cex=mycex1,lwd=1)
  text(x=3.5+1+dr,y=R.loc[,2],srt=0,adj=0,
       labels=R.loc[,5],
       xpd=TRUE,cex=mycex2,lwd=1)
  ###now delete duplicated leftcodes in cc
  ind.dup.cc=which(duplicated(grp.info.cc$icd9cm))
  if(length(ind.dup.cc)>0){
    for(i in unique(cc[ind.dup.cc])){
      mean.weights=apply(Pi.plot[row.names(Pi.plot)==i,],2,mean)
      mean.weights=mean.weights/norm(mean.weights,"2")
      for(j in 1:length(which(row.names(Pi.plot)==i))){
        Pi.plot[row.names(Pi.plot)==i,][j,]=mean.weights
      }
    }
    Pi.plot=Pi.plot[-ind.dup.cc,]
    grp.info.cc=grp.info.cc[-ind.dup.cc,]
    cc=cc[-ind.dup.cc]
  }
  
  ###plot grey lines
  for(i in 1:length(cc)){
    rightweights = Pi.plot[i,]
    rightweights = rightweights[rightweights>0]
    if(length(rightweights)>0){
      if(anyDuplicated(names(rightweights))>0){
        dup.name=names(rightweights)[duplicated(names(rightweights))]
        dup.ind=which(names(rightweights)%in%dup.name)
        dup.max.ind=which.max(rightweights[dup.ind])
        rightweights = rightweights[-dup.ind[-dup.max.ind]]
      }
      rightcodes = names(rightweights)[rightweights>0]
      rightcodes.ind = sapply(rightcodes,FUN=function(x){which(R.loc[,3]%in%x)})
      leftcodes.ind = which(L.loc[,3]==cc[i])
      for(j in 1:length(rightcodes)){
        segments(x0=L.loc[leftcodes.ind,1]+dl,y0=L.loc[leftcodes.ind,2],
                 x1=R.loc[rightcodes.ind[j],1]+dr,y1=R.loc[rightcodes.ind[j],2],
                 lwd = as.numeric(rightweights[j])*5,
                 col="#00000050")
      }
    }
  }
}

