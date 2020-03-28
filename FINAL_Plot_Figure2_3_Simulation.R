#########################################################################################
#########################################################################################
### Make plot (Figure 2) for simulation results                                       ###
### Contact: shixu@umich.edu                                                          ###
#########################################################################################
#########################################################################################
rm(list=ls())
at1=35
at2=6000
mycex=0.8
plot.W_percmismatch = function(){
  myfilepath=""
  scale.p=300
  perm_only=F
  pdf(paste0("~/Dropbox/percent_bias_W_perm",perm_only,Sys.time(),".pdf"),h=3.2,w=8)
  par(mfrow=c(1,2),mar=c(0,0,0,0),mai=c(0.8,0.6,0.5,0.15),omi=c(0,0,0,0))
  ### %mismatch vs bias
  mykk=17
  
  load(paste0(myfilepath,
              paste0("_",mykk,"_",perm_only,"alltt_Wtable"),".RData"))
  perc.bias.spr = Werr.fixkk[,"Fbias.spr"]^2/scale.p
  perc.bias.ols = Werr.fixkk[,"Fbias.ols"]^2/scale.p
  sparsity=Werr.fixkk[,"S"]/Werr.fixkk[,"N"]*100#Werr.fixkk[,"sparcity"]*2*100
  norder=log(Werr.fixkk[,"S"]+1)/log(unique(Werr.fixkk[,"N"]))
  lo.spr <- loess(perc.bias.spr~sparsity)
  lo.ols <- loess(perc.bias.ols~sparsity)
  ymin=min(predict(lo.spr),predict(lo.ols))
  ymax=max(predict(lo.spr),predict(lo.ols))
  plot(sparsity,predict(lo.spr),
       xlab="",ylab="",axes=F,
       type="l",lty=1,lwd=2,
       ylim=c(ymin,ymax));box()
  ####
  axis(side=3,at=seq(min(sparsity),max(sparsity),by=10),tick=T,line=0,labels=rep("",6))
  axis(side=3,at=seq(min(sparsity),max(sparsity),by=10),
       tick=F,line=-0.5,
       labels=round(seq(min(sparsity),max(sparsity),by=10),0))
  mtext(side=3,line=1.5,text="% Mismatch",cex=mycex)
  ####
  axis(side=1,at=sparsity[c(3,7,9,12,16,17)],tick=T,line=0,labels=rep("",6))
  axis(side=1,at=sparsity[c(3,7,9,12,16,17)],tick=F,line=-0.5,
       labels=c(0.5,0.7,0.8,0.85,0.9,0.92))
  nn=as.numeric(round(Werr.fixkk[1,1]/1000,0)*1000)
  mtext(side=1,line=1.5,text=bquote(alpha))
  mtext(side=1,line=1.5,text=paste("(n = ",nn,")"),at=at1)
  # mtext(side=1,line=3,text=
  #         bquote("Number of Mismatched Rows = "*n^alpha*" (n = "*.(nn)*")"))
  #         # bquote("|D("*symbol(I)*", "*Pi*")| = "*n^alpha*" (n = "*.(nn)*")"))
  ####
  axis(side=2,at=seq(ymin,ymax,length.out=5),tick=T,line=0,labels=rep("",5))
  axis(side=2,at=seq(ymin,ymax,length.out=5),
       tick=F,line=-0.5,
       labels=round(seq(ymin,ymax,length.out=5),1))
  mtext(side=2,line=1.5,text="MSE",cex=mycex)
  lines(sparsity,predict(lo.ols),lty=5,lwd=2)
  ## add updated W rslt
  perc.bias.spr = Werr.fixkk[,"Fbias.spr.update"]^2/scale.p
  perc.bias.ols = Werr.fixkk[,"Fbias.ols.update"]^2/scale.p
  lo.spr.update <- loess(perc.bias.spr~sparsity)
  lo.ols.update <- loess(perc.bias.ols~sparsity)
  lines(sparsity,predict(lo.spr.update),
        lty=6,lwd=1)
  lines(sparsity,predict(lo.ols.update),
        lty=3,lwd=1)
  legend("topleft",bty="n",cex=mycex,
         lty=c(1,5,6,3)[c(1,3,2,4)],
         lwd=c(2,2,1,1)[c(1,3,2,4)],
         legend=c(expression(paste("iSphereMAP ",widehat(W)^"[1]")),
                  expression(paste("MT method ",widehat(W)^"[1]")),
                  expression(paste("iSphereMAP ",widehat(W)^"[2]")),
                  expression(paste("MT method ",widehat(W)^"[2]")))[c(1,3,2,4)])
  
  ### N vs bias
  mytt=9
  if(perm_only==T){tt=mytt}else{tt=mytt}
  perc.bias.spr=perc.bias.ols=N=NULL
  for(kk in 2:mykk){
    load(paste0(myfilepath,
                paste0("_",kk,"_",perm_only,"alltt_Wtable"),".RData"))
    perc.bias.spr = c(perc.bias.spr,Werr.fixkk[tt,"Fbias.spr"]^2/scale.p)
    perc.bias.ols = c(perc.bias.ols,Werr.fixkk[tt,"Fbias.ols"]^2/scale.p)
    N=c(N,Werr.fixkk[tt,"N"])
  }
  lo.spr <- loess(perc.bias.spr~N)
  lo.ols <- loess(perc.bias.ols~N)
  ymin=min(predict(lo.spr),predict(lo.ols))
  ymax=max(predict(lo.spr),predict(lo.ols))
  plot(N,predict(lo.spr),
       xlab="",ylab="",axes=F,
       type="l",lty=1,lwd=2,
       ylim=c(ymin,ymax));box()
  lines(N,predict(lo.ols),lty=5,lwd=2)
  ####
  axis(side=1,at=seq(2000,8000,by=2000),tick=T,line=0,labels=rep("",4))
  axis(side=1,at=seq(2000,8000,by=2000),
       tick=F,line=-0.5,
       labels=seq(2000,8000,by=2000))
  mtext(side=1,line=1.5,text="n")
  mtext(side=1,line=1.65,text=bquote("("*alpha*" = 0.8)"),at=at2)
  ####
  axis(side=2,at=seq(ymin,ymax,length.out=5),tick=T,line=0,labels=rep("",5))
  axis(side=2,at=seq(ymin,ymax,length.out=5),
       tick=F,line=-0.5,
       labels=round(seq(ymin,ymax,length.out=5),1))
  mtext(side=2,line=1.5,text="MSE",cex=mycex)
  ## add updated W rslt
  if(perm_only==T){tt=mytt}else{tt=mytt}
  perc.bias.spr=perc.bias.ols=N=NULL
  for(kk in 2:mykk){
    load(paste0(myfilepath,
                paste0("_",kk,"_",perm_only,"alltt_Wtable"),".RData"))
    perc.bias.spr = c(perc.bias.spr,Werr.fixkk[tt,"Fbias.spr.update"]^2/scale.p)
    perc.bias.ols = c(perc.bias.ols,Werr.fixkk[tt,"Fbias.ols.update"]^2/scale.p)
    N=c(N,Werr.fixkk[tt,"N"])
  }
  lo.spr.update <- loess(perc.bias.spr~N)
  lines(N,predict(lo.spr.update),lty=6,lwd=1)
  lo.ols.update <- loess(perc.bias.ols~N)
  lines(N,predict(lo.ols.update),lty=3,lwd=1)
  
  legend("topright",bty="n",cex=mycex,
         lty=c(1,5,6,3)[c(1,3,2,4)],
         lwd=c(2,2,1,1)[c(1,3,2,4)],
         legend=c(expression(paste("iSphereMAP ",widehat(W)^"[1]")),
                  expression(paste("MT method ",widehat(W)^"[1]")),
                  expression(paste("iSphereMAP ",widehat(W)^"[2]")),
                  expression(paste("MT method ",widehat(W)^"[2]")))[c(1,3,2,4)])
  mtext(side=1,line=-1.3,outer=T,cex=mycex,
        text=bquote("Number of Mismatched Rows = "*n^alpha))
  dev.off()
}
plot.W_percmismatch()

#########################################################################################
#########################################################################################
### Make plot (Figure 2) for simulation results                                       ###
### Contact: shixu@umich.edu                                                          ###
#########################################################################################
#########################################################################################
rm(list=ls())
at1=35
at2=6000
mycex=0.8
myloess=T
plot.Pi_percmismatch = function(){
  my.filepath=""
  mykk=17
  perm_only=F
  pdf(paste0("~/Dropbox/percent_bias_Pi_perm",perm_only,Sys.time(),".pdf"),h=3.2,w=8)
  par(mfrow=c(1,2),mar=c(0,0,0,0),mai=c(0.8,0.6,0.5,0.15),omi=c(0,0,0,0))
  ##########################################################################################
  ##########################################################################################
  ### N vs bias
  mytt=9
  if(perm_only==T){tt=mytt}else{tt=mytt}
  perc.bias.spr=perc.bias.ols=N=NULL
  for(kk in 2:mykk){
    load(paste0(my.filepath,
                paste0("_",kk,"_",perm_only,"alltt_Pitable"),".RData"))
    Werr.fixkk = t(data.frame(lapply(Wrslt.allseed.alltt,FUN=function(x){apply(x,2,mean,na.rm=T)})))
    Pierr.fixkk = t(data.frame(lapply(Pirslt.allseed.alltt,FUN=function(x){apply(x,2,mean,na.rm=T)})))
    perc.bias.spr = c(perc.bias.spr,Pierr.fixkk[tt,"match.rate"])
    perc.bias.ols = c(perc.bias.ols,Pierr.fixkk[tt,"match.rate.nogrp"])
    N=c(N,Werr.fixkk[tt,"N"])
  }
  if(myloess==T){
    lo.spr <- loess(perc.bias.spr~N)
    lo.ols <- loess(perc.bias.ols~N)
    ymin=min(predict(lo.spr),predict(lo.ols))
    ymax=max(predict(lo.spr),predict(lo.ols))
    plot(N,predict(lo.spr),
         xlab="",ylab="",axes=F,
         type="l",lty=1,lwd=2,
         ylim=c(ymin,ymax));box()
    lines(N,predict(lo.ols),lty=5,lwd=2)
  }else{
    ymin=min(perc.bias.spr,perc.bias.ols)
    ymax=max(perc.bias.spr,perc.bias.ols)
    plot(N,perc.bias.spr,
         xlab="",ylab="",axes=F,
         type="l",lty=1,lwd=2,
         ylim=c(ymin,ymax));box()
    lines(N,perc.bias.ols,lty=5,lwd=2)
  }
  ####
  axis(side=1,at=seq(2000,8000,by=2000),tick=T,line=0,labels=rep("",4))
  axis(side=1,at=seq(2000,8000,by=2000),
       tick=F,line=-0.5,
       labels=seq(2000,8000,by=2000))
  mtext(side=1,line=1.5,text="n")
  mtext(side=1,line=1.65,text=bquote("("*alpha*" = 0.8)"),at=at2)
  ####
  axis(side=2,at=seq(ymin,ymax,length.out=5),tick=T,line=0,labels=rep("",5))
  axis(side=2,at=seq(ymin,ymax,length.out=5),
       tick=F,line=-0.5,
       labels=round(seq(ymin,ymax,length.out=5),1))
  mtext(side=2,line=1.5,text="1-to-1 Match Rate",cex=mycex)
  legend("bottomright",bty="n",cex=mycex,
         lty=c(1,5),
         lwd=c(2,2),
         legend=c("iSphereMAP","MT method w/o grp info"))
  
  
  ### N vs bias
  mytt=9
  if(perm_only==T){tt=mytt}else{tt=mytt}
  perc.bias.spr=perc.bias.ols=N=NULL
  for(kk in 2:mykk){
    load(paste0(my.filepath,
                paste0("_",kk,"_",perm_only,"alltt_Pitable"),".RData"))
    Werr.fixkk = t(data.frame(lapply(Wrslt.allseed.alltt,FUN=function(x){apply(x,2,mean,na.rm=T)})))
    Pierr.fixkk = t(data.frame(lapply(Pirslt.allseed.alltt,FUN=function(x){apply(x,2,mean,na.rm=T)})))
    perc.bias.spr = c(perc.bias.spr,Pierr.fixkk[tt,"C.error_F"])
    perc.bias.ols = c(perc.bias.ols,Pierr.fixkk[tt,"C.error_F.nogrp"])
    N=c(N,Werr.fixkk[tt,"N"])
  }
  lo.spr <- loess(perc.bias.spr~N)
  lo.ols <- loess(perc.bias.ols~N)
  ymin=min(predict(lo.spr),predict(lo.ols))
  ymax=max(predict(lo.spr),predict(lo.ols))
  plot(N,predict(lo.spr),
       xlab="",ylab="",axes=F,
       type="l",lty=1,lwd=2,
       ylim=c(ymin,ymax));box()
  lines(N,predict(lo.ols),lty=5,lwd=2)
  ####
  axis(side=1,at=seq(2000,8000,by=2000),tick=T,line=0,labels=rep("",4))
  axis(side=1,at=seq(2000,8000,by=2000),
       tick=F,line=-0.5,
       labels=seq(2000,8000,by=2000))
  mtext(side=1,line=1.5,text="n")
  mtext(side=1,line=1.65,text=bquote("("*alpha*" = 0.8)"),at=at2)
  ####
  axis(side=2,at=seq(ymin,ymax,length.out=5),tick=T,line=0,labels=rep("",5))
  axis(side=2,at=seq(ymin,ymax,length.out=5),
       tick=F,line=-0.5,
       labels=round(seq(ymin,ymax,length.out=5)*1e5,1))
  mtext(side=2,line=1.5,text=bquote("MSE of 1-to-m Weight (x"*10^{-5}*")"),cex=mycex)
  legend("topright",bty="n",cex=mycex,
         lty=c(1,5),
         lwd=c(2,2),
         legend=c("iSphereMAP","MT method w/o grp info"))
  
  dev.off()
  
}
plot.Pi_percmismatch()



