#########################################################################################
#########################################################################################
### Make plot (Figure F.4) for simulation results                                     ###
### Contact: shixu@umich.edu                                                          ###
#########################################################################################
#########################################################################################
rm(list=ls());dev.off()
at1=35
at2=4500
mycex=0.8
scale.p=300
myloess=F
allargs = list(c(0,0,0,0,0),c(1,0,0,0,0),c(0,1,0,0,0),c(0,0,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1))
my.filepath=""
mykk=7;
mytt=(9-1)/2;tt=mytt
read.rslt.fixN=function(allargs.i){
  args=allargs[[allargs.i]]
  wrong_grp_info=(args[1]==1)
  compute_W_k=(args[2]==1)
  add_mismatch_1to1=(args[3]==1)
  perm_only=(args[4]==1)
  MVN=(args[5]==1)
  load(paste0(my.filepath,"allrslt_",mykk,"_",wrong_grp_info,"_",
              compute_W_k,"_",add_mismatch_1to1,"_",perm_only,"_",MVN,"_.RData"))
  Werr.fixkk = t(data.frame(lapply(Wrslt.allseed.alltt,FUN=function(x){apply(x,2,mean,na.rm=T)})))
  sparsity=Werr.fixkk[,"S"]/Werr.fixkk[,"N"]*100#Werr.fixkk[,"sparcity"]*2*100
  norder=log(Werr.fixkk[,"S"]+1)/log(unique(Werr.fixkk[,"N"]))
  #####
  perc.bias.spr = Werr.fixkk[,"Fbias.spr"]^2/scale.p
  perc.bias.ols = Werr.fixkk[,"Fbias.ols"]^2/scale.p
  perc.bias.spr.update = Werr.fixkk[,"Fbias.spr.update"]^2/scale.p
  perc.bias.ols.update = Werr.fixkk[,"Fbias.ols.update"]^2/scale.p
  N = Werr.fixkk[,"N"]
  return(list(perc.bias.spr=perc.bias.spr,perc.bias.ols=perc.bias.ols,
              perc.bias.spr.update=perc.bias.spr.update,perc.bias.ols.update=perc.bias.ols.update,
              sparsity=sparsity,norder=norder,N=N))
}
plot.rslt.fixN=function(rslt,addlabel=T,ylabname){
  perc.bias.spr = rslt$perc.bias.spr
  perc.bias.ols = rslt$perc.bias.ols
  perc.bias.spr.update=rslt$perc.bias.spr.update
  perc.bias.ols.update=rslt$perc.bias.ols.update
  sparsity=rslt$sparsity
  norder=rslt$norder
  N=rslt$N
  #####
  lo.spr <- loess(perc.bias.spr~sparsity)
  lo.ols <- loess(perc.bias.ols~sparsity)
  lo.spr.update <- loess(perc.bias.spr.update~sparsity)
  lo.ols.update <- loess(perc.bias.ols.update~sparsity)
  #####
  ymin=min(na.rm=T,predict(lo.spr),predict(lo.ols),predict(lo.spr.update),predict(lo.ols.update))
  ymax=max(na.rm=T,predict(lo.spr),predict(lo.ols),predict(lo.spr.update),predict(lo.ols.update))
  plot(sparsity,predict(lo.spr),
       xlab="",ylab="",axes=F,
       type="l",lty=1,lwd=2,
       ylim=c(ymin,ymax));box()
  lines(sparsity,predict(lo.ols),lty=5,lwd=2)
  ## add updated W rslt
  lines(sparsity,predict(lo.spr.update),
        lty=6,lwd=1,col="red")
  lines(sparsity,predict(lo.ols.update),
        lty=3,lwd=1,col="red")
  ####
  axis(side=3,at=seq(min(na.rm=T,sparsity),max(na.rm=T,sparsity),by=10),tick=T,line=0,labels=rep("",length(seq(min(na.rm=T,sparsity),max(na.rm=T,sparsity),by=10))))
  axis(side=3,at=seq(min(na.rm=T,sparsity),max(na.rm=T,sparsity),by=10),
       tick=F,line=-0.5,
       labels=round(seq(min(na.rm=T,sparsity),max(na.rm=T,sparsity),by=10),0),cex.axis=mycex)
  mtext(side=3,line=1.5,text="% Mismatch",cex=mycex)
  ####
  axis(side=1,at=sparsity[(c(3,9,13,17,19)-1)/2],tick=T,line=0,labels=rep("",5))
  axis(side=1,at=sparsity[(c(3,9,13,17,19)-1)/2],tick=F,line=-0.5,
       labels=c(0.5,0.8,0.85,0.9,0.92),cex.axis=mycex)
  # axis(side=1,at=sparsity[(c(3,7,9,13,17,19)-1)/2],tick=T,line=0,labels=rep("",6))
  # axis(side=1,at=sparsity[(c(3,7,9,13,17,19)-1)/2],tick=F,line=-0.5,
  #      labels=c(0.5,0.7,0.8,0.85,0.9,0.92),cex.axis=mycex)
  nn=as.numeric(round(mean(N)/1000,0)*1000)
  mtext(side=1,line=1.5,text=bquote(alpha))
  mtext(side=1,line=1.5,text=paste0("(n = ",nn,")"),at=at1)
  ####
  axis(side=2,at=seq(ymin,ymax,length.out=5),cex.axis=mycex,tick=T,line=0,labels=rep("",5))
  if(ymax<1e-1){
    axis(side=2,at=seq(ymin,ymax,length.out=5),cex.axis=mycex,
         tick=F,line=-0.5,
         labels=round(seq(ymin,ymax,length.out=5)*1e2,2))
    if(ylabname=="MSE"){
      mtext(side=2,line=1.5,text=bquote("MSE (x"*10^{-2}*")"),cex=mycex)
    }else{
      mtext(side=2,line=1.5,text=bquote("MSE (x"*10^{-2}*")"),cex=mycex)
    }
    
  }else{
    axis(side=2,at=seq(ymin,ymax,length.out=5),cex.axis=mycex,
         tick=F,line=-0.5,
         labels=round(seq(ymin,ymax,length.out=5),2))
    mtext(side=2,line=1.5,text=ylabname,cex=mycex)
  }
  if(addlabel){
    legend("topleft",bty="n",cex=mycex,
           lty=c(1,5,6,3)[c(1,3,2,4)],
           lwd=c(2,2,1,1)[c(1,3,2,4)],
           col=c("black","black","red","red")[c(1,3,2,4)],
           legend=c(expression(paste("iSphereMAP step I: ",widehat(W)^"[1]")),
                    expression(paste("MT method  step I")),
                    expression(paste("iSphereMAP step III: ",widehat(W)^"[2]")),
                    expression(paste("MT method step III")))[c(1,3,2,4)])
    # legend("topleft",bty="n",cex=mycex,
    #        lty=c(1,5,6,3)[c(1,3,2,4)],
    #        lwd=c(2,2,1,1)[c(1,3,2,4)],
    #        col=c("black","black","red","red")[c(1,3,2,4)],
    #        legend=c(expression(paste("iSphereMAP ",widehat(W)^"[1]")),
    #                 expression(paste("MT method ",widehat(W)^"[1]")),
    #                 expression(paste("iSphereMAP ",widehat(W)^"[2]")),
    #                 expression(paste("MT method ",widehat(W)^"[2]")))[c(1,3,2,4)])
  }
}
compare.rslt.fixN=function(rslt0,rslt1){
  perc.bias.spr0 = rslt0$perc.bias.spr
  perc.bias.ols0 = rslt0$perc.bias.ols
  perc.bias.spr.update0=rslt0$perc.bias.spr.update
  perc.bias.ols.update0=rslt0$perc.bias.ols.update
  sparsity0=rslt0$sparsity
  norder0=rslt0$norder
  N0=rslt0$N
  perc.bias.spr1 = rslt1$perc.bias.spr
  perc.bias.ols1 = rslt1$perc.bias.ols
  perc.bias.spr.update1=rslt1$perc.bias.spr.update
  perc.bias.ols.update1=rslt1$perc.bias.ols.update
  sparsity1=rslt1$sparsity
  norder1=rslt1$norder
  N1=rslt1$N
  
  ####compare
  # perc.bias.spr=perc.bias.spr1-perc.bias.spr0
  # perc.bias.ols=perc.bias.ols1-perc.bias.ols0
  # perc.bias.spr.update=perc.bias.spr.update1-perc.bias.spr.update0
  # perc.bias.ols.update=perc.bias.ols.update1-perc.bias.ols.update0
  # sparsity=(sparsity1+sparsity0)/2
  # norder=(norder1+norder0)/2
  # N=(N1+N0)/2
  perc.bias.spr=perc.bias.spr.update0
  perc.bias.ols=perc.bias.ols.update0 
  perc.bias.spr.update=perc.bias.spr.update1   
  perc.bias.ols.update=perc.bias.ols.update1
  sparsity=(sparsity1+sparsity0)/2
  norder=(norder1+norder0)/2
  N=(N1+N0)/2
  
  rslt=list(perc.bias.spr=perc.bias.spr,perc.bias.ols=perc.bias.ols,
            perc.bias.spr.update=perc.bias.spr.update,perc.bias.ols.update=perc.bias.ols.update,
            sparsity=sparsity,norder=norder,N=N)
  return(rslt)
}
read.rslt=function(allargs.i){
  args=allargs[[allargs.i]]
  wrong_grp_info=(args[1]==1)
  compute_W_k=(args[2]==1)
  add_mismatch_1to1=(args[3]==1)
  perm_only=(args[4]==1)
  MVN=(args[5]==1)
  perc.bias.spr=perc.bias.ols=N=perc.bias.spr.update=perc.bias.ols.update=NULL#sparsity=norder=NULL
  for(kk in 1:mykk){
    readrslt=try(load(paste0(my.filepath,"allrslt_",kk,"_",wrong_grp_info,"_",
                             compute_W_k,"_",add_mismatch_1to1,"_",perm_only,"_",MVN,"_.RData")),silent=T)
    if(!inherits(readrslt, "try-error")){
      Werr.fixkk = t(data.frame(lapply(Wrslt.allseed.alltt,FUN=function(x){apply(x,2,mean,na.rm=T)})))
      perc.bias.spr = c(perc.bias.spr,Werr.fixkk[tt,"Fbias.spr"]^2/scale.p)
      perc.bias.ols = c(perc.bias.ols,Werr.fixkk[tt,"Fbias.ols"]^2/scale.p)
      perc.bias.spr.update=c(perc.bias.spr.update,Werr.fixkk[tt,"Fbias.spr.update"]^2/scale.p)
      perc.bias.ols.update=c(perc.bias.ols.update,Werr.fixkk[tt,"Fbias.ols.update"]^2/scale.p)
      N=c(N,Werr.fixkk[tt,"N"])
      # sparsity = cbind(sparsity,Werr.fixkk[,"S"]/Werr.fixkk[,"N"]*100)#Werr.fixkk[,"sparcity"]*2*100
      # norder = cbind(norder,log(Werr.fixkk[,"S"]+1)/log(unique(Werr.fixkk[,"N"])))
    }else{
      perc.bias.spr = c(perc.bias.spr,NA)
      perc.bias.ols = c(perc.bias.ols,NA)
      perc.bias.spr.update=c(perc.bias.spr.update,NA)
      perc.bias.ols.update=c(perc.bias.ols.update,NA)
      N=c(N,NA)
      # sparsity=cbind(sparsity,rep(NA,9))
      # norder=cbind(norder,rep(NA,9))
    }
  }
  return(list(perc.bias.spr=perc.bias.spr,perc.bias.ols=perc.bias.ols,
              perc.bias.spr.update=perc.bias.spr.update,perc.bias.ols.update=perc.bias.ols.update,
              # sparsity=sparsity,norder=norder,
              N=N))
}
plot.rslt=function(rslt,addlabel=T,ylabname){
  perc.bias.spr = rslt$perc.bias.spr
  perc.bias.ols = rslt$perc.bias.ols
  perc.bias.spr.update=rslt$perc.bias.spr.update
  perc.bias.ols.update=rslt$perc.bias.ols.update
  N=rslt$N
  
  lo.spr <- loess(perc.bias.spr~N)
  lo.ols <- loess(perc.bias.ols~N)
  lo.spr.update <- loess(perc.bias.spr.update~N)
  lo.ols.update <- loess(perc.bias.ols.update~N)
  ymin=min(na.rm=T,predict(lo.spr),predict(lo.ols),predict(lo.spr.update),predict(lo.ols.update))
  ymax=max(na.rm=T,predict(lo.spr),predict(lo.ols),predict(lo.spr.update),predict(lo.ols.update))
  N=N[!is.na(N)]
  plot(N,predict(lo.spr),
       xlab="",ylab="",axes=F,
       type="l",lty=1,lwd=2,
       ylim=c(ymin,ymax));box()
  lines(N,predict(lo.ols),lty=5,lwd=2)
  lines(N,predict(lo.spr.update),lty=6,lwd=1,col="red")
  lines(N,predict(lo.ols.update),lty=3,lwd=1,col="red")
  ####
  axis(side=1,at=seq(2000,8000,by=2000),tick=T,line=0,labels=rep("",length(seq(2000,8000,by=2000))))
  axis(side=1,at=seq(2000,8000,by=2000),
       tick=F,line=-0.5,
       labels=seq(2000,8000,by=2000))
  mtext(side=1,line=1.5,text="n")
  mtext(side=1,line=1.65,text=bquote("("*alpha*" = 0.7)"),at=at2)
  ####
  axis(side=2,at=seq(ymin,ymax,length.out=5),cex.axis=mycex,tick=T,line=0,labels=rep("",5))
  if(ymax<1e-1){
    axis(side=2,at=seq(ymin,ymax,length.out=5),cex.axis=mycex,
         tick=F,line=-0.5,
         labels=round(seq(ymin,ymax,length.out=5)*1e2,2))
    if(ylabname=="MSE"){
      mtext(side=2,line=1.5,text=bquote("MSE (x"*10^{-2}*")"),cex=mycex)
    }else{
      mtext(side=2,line=1.5,text=bquote("MSE (x"*10^{-2}*")"),cex=mycex)
    }
    
  }else{
    axis(side=2,at=seq(ymin,ymax,length.out=5),cex.axis=mycex,
         tick=F,line=-0.5,
         labels=round(seq(ymin,ymax,length.out=5),2))
    mtext(side=2,line=1.5,text=ylabname,cex=mycex)
  }
  mtext(side=1,line=-1.3,outer=T,cex=mycex,
        text=bquote("Number of Mismatched Rows = "*n^alpha))
  if(addlabel){
    legend("topleft",bty="n",cex=mycex,
           lty=c(1,5,6,3)[c(1,3,2,4)],
           lwd=c(2,2,1,1)[c(1,3,2,4)],
           col=c("black","black","red","red")[c(1,3,2,4)],
           legend=c(expression(paste("iSphereMAP step I: ",widehat(W)^"[1]")),
                    expression(paste("MT method  step I")),
                    expression(paste("iSphereMAP step III: ",widehat(W)^"[2]")),
                    expression(paste("MT method step III")))[c(1,3,2,4)])
    # legend("topright",bty="n",cex=mycex,
    #        lty=c(1,5,6,3)[c(1,3,2,4)],
    #        lwd=c(2,2,1,1)[c(1,3,2,4)],
    #        col=c("black","black","red","red")[c(1,3,2,4)],
    #        legend=c(expression(paste("iSphereMAP ",widehat(W)^"[1]")),
    #                 expression(paste("MT method ",widehat(W)^"[1]")),
    #                 expression(paste("iSphereMAP ",widehat(W)^"[2]")),
    #                 expression(paste("MT method ",widehat(W)^"[2]")))[c(1,3,2,4)])
  }
}
compare.rslt=function(rslt0,rslt1){
  perc.bias.spr0 = rslt0$perc.bias.spr
  perc.bias.ols0 = rslt0$perc.bias.ols
  perc.bias.spr.update0=rslt0$perc.bias.spr.update
  perc.bias.ols.update0=rslt0$perc.bias.ols.update
  N0=rslt0$N
  perc.bias.spr1 = rslt1$perc.bias.spr
  perc.bias.ols1 = rslt1$perc.bias.ols
  perc.bias.spr.update1=rslt1$perc.bias.spr.update
  perc.bias.ols.update1=rslt1$perc.bias.ols.update
  N1=rslt1$N
  ####compare
  # perc.bias.spr=perc.bias.spr1-perc.bias.spr0
  # perc.bias.ols=perc.bias.ols1-perc.bias.ols0
  # perc.bias.spr.update=perc.bias.spr.update1-perc.bias.spr.update0
  # perc.bias.ols.update=perc.bias.ols.update1-perc.bias.ols.update0
  # (N=(N1+N0)/2)
  perc.bias.spr=perc.bias.spr.update0
  perc.bias.ols=perc.bias.ols.update0
  perc.bias.spr.update=perc.bias.spr.update1
  perc.bias.ols.update=perc.bias.ols.update1
  (N=(N1+N0)/2)
  rslt=list(perc.bias.spr=perc.bias.spr,perc.bias.ols=perc.bias.ols,
            perc.bias.spr.update=perc.bias.spr.update,perc.bias.ols.update=perc.bias.ols.update,N=N)
  return(rslt)
}
wrong_grp_info=F
compute_W_k=F
add_mismatch_1to1=T
perm_only=F
MVN=F

pdf(paste0("percent_bias_W_add_mismatch_1to1_",Sys.time(),".pdf"),h=3.2*1.25,w=8*1.25)
par(mfrow=c(1,2),mar=c(0,0,0,0),mai=c(0.8,0.6,0.5,0.15),omi=c(0,0,0,0))
# #################################################################
# ### N vs bias
# ##################################################################
# rslt=read.rslt.fixN(allargs.i=4)
# plot.rslt.fixN(rslt,addlabel=T,ylabname="MSE")
# rslt=read.rslt(allargs.i=4)
# plot.rslt(rslt,addlabel=F,ylabname="MSE")
##################################################################
### Compare allargs FFFFF to FTFFF
##################################################################
rslt0=read.rslt.fixN(allargs.i=1)
rslt1=read.rslt.fixN(allargs.i=4)
rslt=compare.rslt.fixN(rslt0,rslt1)
plot.rslt.fixN(rslt,addlabel=F,ylabname="MSE");abline(h=0)
legend(x=0.5,y=0.35,bty="n",cex=mycex,
       lty=c(1,5,6,3)[c(1,3,2,4)],
       lwd=c(2,2,1,1)[c(1,3,2,4)],
       col=c("black","black","red","red")[c(1,3,2,4)],
       title=expression(paste("Step III: refined estimate ",widehat(W)^"[2]")),
       legend=c("iSphereMAP w/ correctly matched data only",
                "MT method w/ correctly matched data only",
                "iSphereMAP w/ 1-1 (mis)matched data",
                "MT method w/ 1-1 (mis)matched data")[c(1,3,2,4)])
rslt0=read.rslt(allargs.i=1)
rslt1=read.rslt(allargs.i=4)
rslt=compare.rslt(rslt0,rslt1)
plot.rslt(rslt,addlabel=F,ylabname="MSE");abline(h=0)
dev.off()



