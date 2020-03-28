###################################################################
###################################################################
### Application (Figures 4)                                     ###
### Mapping From Partners HealthCare to VA                      ###
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
vec.size = 300; fine_grp=T
load(paste0(myfilepath,"cleaned_dat_nvec",vec.size,
            "finegrp",fine_grp,"VApartner_08202018rollup.RData"))
### only use join rows for training dat
join.names = intersect(row.names(X.nonzero),row.names(Z.nonzero))
X.nonzero = X.nonzero[match(join.names,row.names(X.nonzero)),]
Z.nonzero = Z.nonzero[match(join.names,row.names(Z.nonzero)),]
grp.info = phewas.grp.X$X.phewas
grp.info = grp.info[match(join.names,grp.info$icd),]
X = X.test = X.train = as.matrix(X.nonzero[match(join.names,row.names(X.nonzero)),])
Y = Y.test = Y.train = as.matrix(Z.nonzero[match(join.names,row.names(Z.nonzero)),])
nrow(grp.info);nrow(X);nrow(Y)
plot(match(grp.info$icd,row.names(Y)))
### order X Y and grp.info by g.index
order.g.index = order(grp.info$g.index,grp.info$icd)
X = X.test = X.train = X[order.g.index,] ##VHA
Y = Y.test = Y.train = Y[order.g.index,] ##PHS
grp.info = grp.info[order.g.index,]
plot(match(grp.info$icd,row.names(Y))) ## should be a monotonic increasing line


######################################
##        iSphereMAP step I         ##
##       sperical regression        ##
######################################
X.norms = rowNorms(as.matrix(X))
if(length(which(X.norms!=1))){
  X = X/rowNorms(as.matrix(X))
  Y = Y/rowNorms(as.matrix(Y))
}
Beta.spr = gradient_update_nogrp(X,Y,alpha=1,convergence=1e-10)
Yhat.spr = X%*%Beta.spr
#####################################################################################################
######################################################################################################
## remove any K>p case
(biggrps=names(which(table(grp.info$g.index)>vec.size)))
ind=which(grp.info$g.index%in%biggrps)
if(length(ind)>0){X=X[-ind,];Y=Y[-ind,];grp.info=grp.info[-ind,]}
## update dimensions
(p=ncol(X))
(N = nrow(grp.info))
ugrp = unique(grp.info$g.index);table(is.na(ugrp))
(n = length(ugrp))
######################################
##       iSphereMAP step II         ##
######################################
###
lambda.cv =  1/sqrt(2)
mymethod="OLS" 
########### make final plot
final=T
gg=as.numeric(grp.info$g.index);ugg=unique(gg)
if(final==T){
  pdf(file=paste0(myfilepath,
                  "plot_PHStoVHA_final3groups.pdf"),h=3.7,w=8)
  ggpick.all=ugg
  ind.gg=which(ggpick.all%in%c(172,277,448)) ## selected groups
  Pi = fitpi_CV(estPi=mymethod,Ytrgt=Y,Yhat=Yhat.spr,
                lambda.cv=1/sqrt(2),n,grp.info,ugrp)
  Pi = Pi/rowNorms(Pi)
  gg=as.numeric(grp.info$g.index);ugg=unique(gg)
  plot(match(row.names(Y),grp.info$icd))###must be abline(0,1)otherwise WRONG
  # plot(match(row.names(Y),colnames(Pi)))
  # plot(match(row.names(X),row.names(Pi)))
  ### reorder grp.info,Pi
  grp.info.save=grp.info
  grp.info=grp.info.save[order(grp.info.save$icd),]
  Pi=Pi[order(grp.info$icd),order(grp.info$icd)]
  row.names(Pi)=row.names(Y)
  colnames(Pi)=row.names(Yhat.spr)
  t=grp.info[grp.info$g.index%in%ggpick.all[ind.gg],]
  labels=read.csv(paste0(myfilepath,"phecode_icd9_rolled.csv"))
  labels$newcode = cleancode(as.character(labels$ICD9))
  labels$ICD9=as.character(labels$ICD9)
  labels$ICD9.String=as.character(labels$ICD9.String)
  labels[grep("Abnormality of gait",labels$ICD9.String),]$ICD9.String="Abnormality of gait"
  labels=labels[-which(labels$ICD9.String%in%c("Pain in joint","Symptoms involving nervous and musculoskeletal systems")),]
  t=merge(t,labels,by.x="icd",by.y="newcode")
  t[,4]=as.character(t[,4])
  t[t[,4]=="719.4",4]="719.40"
  t[t[,4]=="519",4]="519.0"
  # par(mfrow=c(3,1))
  for(ggpick in ggpick.all[ind.gg]){
    Pi.select=Pi[which(gg%in%ggpick),which(gg%in%ggpick)]
    cc =row.names(Pi.select)
    grp.info.cc=t[match(cc,t$icd),]
    print(length(cc)==nrow(grp.info.cc))
    cc.match = list()
    for(i in 1:length(cc)){
      pi=Pi.select[row.names(Pi.select)==cc[i],]
      cc.match=append(cc.match,list(cc[pi>0]))
    }
    
    L = matrix(NA,nrow=length(cc),ncol=5)
    row.names(L) = paste(grp.info.cc[,1],grp.info.cc[,5])
    R=L
    R.loc = data.frame(R); R.loc[,1]=3; R.loc[,2]=seq(0,length(R),length.out=nrow(R))
    L.loc = data.frame(L); L.loc[,1]=1; L.loc[,2]=seq(0,length(R),length.out=nrow(L))
    L.loc[,5]=R.loc[,5]=grp.info.cc[,5]
    L.loc[,4]=R.loc[,4]=grp.info.cc[,4]
    L.loc[,3]=R.loc[,3]=grp.info.cc[,1]
    R.loc[row.names(R)=="",1]=NA
    if(ggpick==277){
      ind.dyspnea=grep("Other dyspnea and",L.loc[,5])
      L.loc[,5][ind.dyspnea]=R.loc[,5][ind.dyspnea]=NA
      #"Other dyspnea and\nrespiratory abnormality"
    }
    
    mycex2=1.1;mycex1=1.2
    plot(
      x=L.loc[,1],y=L.loc[,2],
      #xlim=c(0,4),ylim=range(R.loc[,2]),
      xlim=c(0-10,4+10),ylim=range(R.loc[,2]),
      type="p",pch=16,
      xlab="",ylab="",axes=F
    )
    text(x=0.5+0.25,y=L.loc[,2],srt=0,adj=1,font=2,
         labels=L.loc[,4],xpd=TRUE,cex=mycex1,lwd=1)
    text(x=0.5-2.25,y=L.loc[,2],srt=0,adj=1,
         labels=L.loc[,5],
         xpd=TRUE,cex=mycex2,lwd=1)
    if(length(which(is.na(L.loc[,5])))>0&&ggpick==277){
      text(x=0.5-2.25,
           y=c(L.loc[is.na(L.loc[,5]),2]+1.3,L.loc[is.na(L.loc[,5]),2]-1.3),
           srt=0,adj=1,
           labels=c("Other dyspnea and","respiratory abnormality"),
           xpd=TRUE,cex=mycex2,lwd=1)
    }
    
    
    points(
      x=R.loc[,1],y=R.loc[,2],
      pch=16
    )
    text(x=3.5-0.25,y=R.loc[,2],srt=0,adj=0,font=2,
         labels=R.loc[,4],xpd=TRUE,cex=mycex1,lwd=1)
    text(x=3.5+2.25,y=R.loc[,2],srt=0,adj=0,
         labels=R.loc[,5],
         xpd=TRUE,cex=mycex2,lwd=1)
    if(length(which(is.na(R.loc[,5])))>0&&ggpick==277){
      text(x=3.5+2.25,
           y=c(R.loc[is.na(R.loc[,5]),2]+1.3,R.loc[is.na(R.loc[,5]),2]-1.3),
           srt=0,adj=0,
           labels=c("Other dyspnea and","respiratory abnormality"),
           xpd=TRUE,cex=mycex2,lwd=1)
    }
    codes=cc.match
    for(i in 1:length(codes)){
      rightcodes = codes[[i]]
      rightweights = Pi.select[i,]
      rightweights = rightweights[rightweights>0]
      # print(paste0(length(rightcodes),"_",length(rightweights)))
      rightcodes.ind = which(
        R.loc[,3]%in%rightcodes)
      for(j in 1:length(rightcodes)){
        segments(x0=L.loc[i,1],y0=L.loc[i,2],
                 x1=R.loc[rightcodes.ind[j],1],y1=R.loc[rightcodes.ind[j],2],
                 lwd = as.numeric(rightweights[j])*5,
                 col="#00000050")
      }
    }
    # title(paste(ggpick,unique(grp.info.cc$Excl..Phenotypes)))
  }
  
  dev.off()
}



system("say done")
