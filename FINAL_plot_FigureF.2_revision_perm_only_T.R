#########################################################################################
#########################################################################################
### Make plot (Figure F.2) for simulation results                                     ###
### Contact: shixu@umich.edu                                                          ###
#########################################################################################
#########################################################################################
rm(list=ls());dev.off()
at1=35
at2=6600
mycex=0.8
myloess=F
allargs = list(c(0,0,0,0,0),c(1,0,0,0,0),c(0,1,0,0,0),c(0,0,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1))
my.filepath=""
mykk=9;mytt=(9-1)/2;tt=mytt
read.rslt=function(allargs.i,name1="match.rate",name2="match.rate.nogrp"){
  args=allargs[[allargs.i]]
  wrong_grp_info=(args[1]==1)
  compute_W_k=(args[2]==1)
  add_mismatch_1to1=(args[3]==1)
  perm_only=(args[4]==1)
  MVN=(args[5]==1)
  perc.bias.spr=perc.bias.ols=N=perc.bias.spr.update=perc.bias.ols.update=NULL
  for(kk in 2:mykk){
    readrslt=try(load(paste0(my.filepath,"allrslt_",kk,"_",wrong_grp_info,"_",
                             compute_W_k,"_",add_mismatch_1to1,"_",perm_only,"_",MVN,"_.RData")),silent=T)
    if(!inherits(readrslt, "try-error")){
      Werr.fixkk = t(data.frame(lapply(Wrslt.allseed.alltt,FUN=function(x){apply(x,2,mean,na.rm=T)})))
      Pierr.fixkk = t(data.frame(lapply(Pirslt.allseed.alltt,FUN=function(x){apply(x,2,mean,na.rm=T)})))
      Pierr.update.fixkk = t(data.frame(lapply(Pirslt.update.allseed.alltt,FUN=function(x){apply(x,2,mean,na.rm=T)})))
      perc.bias.spr = c(perc.bias.spr,Pierr.fixkk[tt,name1])
      perc.bias.ols = c(perc.bias.ols,Pierr.fixkk[tt,name2])
      perc.bias.spr.update = c(perc.bias.spr.update,Pierr.update.fixkk[tt,name1])
      perc.bias.ols.update = c(perc.bias.ols.update,Pierr.update.fixkk[tt,name2])
      N=c(N,Werr.fixkk[tt,"N"])
    }else{
      perc.bias.spr = c(perc.bias.spr,NA)
      perc.bias.ols = c(perc.bias.ols,NA)
      perc.bias.spr.update=c(perc.bias.spr.update,NA)
      perc.bias.ols.update=c(perc.bias.ols.update,NA)
      N=c(N,NA)
    }
  }
  return(list(perc.bias.spr=perc.bias.spr,perc.bias.ols=perc.bias.ols,
              perc.bias.spr.update=perc.bias.spr.update,perc.bias.ols.update=perc.bias.ols.update,N=N))
}
plot.rslt=function(rslt,ylabname="1-to-1 Match Rate",addlabel=T){
  perc.bias.spr = rslt$perc.bias.spr
  perc.bias.ols = rslt$perc.bias.ols
  perc.bias.spr.update=rslt$perc.bias.spr.update
  perc.bias.ols.update=rslt$perc.bias.ols.update
  N=rslt$N
  if(myloess==T){
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
  }else{
    ymin=min(na.rm=T,perc.bias.spr,perc.bias.ols,perc.bias.spr.update,perc.bias.ols.update)
    ymax=max(na.rm=T,perc.bias.spr,perc.bias.ols,perc.bias.spr.update,perc.bias.ols.update)
    N=N[!is.na(N)]
    perc.bias.spr=perc.bias.spr[!is.na(perc.bias.spr)]
    perc.bias.ols=perc.bias.ols[!is.na(perc.bias.ols)]
    plot(N,perc.bias.spr,
         xlab="",ylab="",axes=F,
         type="l",lty=1,lwd=2,
         ylim=c(ymin,ymax));box()
    lines(N,perc.bias.ols,lty=5,lwd=2)
  }
  #### add update rslt
  if(myloess==T){
    N=N[!is.na(N)]
    lines(N,predict(lo.spr.update),lty=6,lwd=1,col="red")
    lines(N,predict(lo.ols.update),lty=3,lwd=1,col="red")
  }else{
    N=N[!is.na(N)]
    perc.bias.spr.update=perc.bias.spr.update[!is.na(perc.bias.spr.update)]
    perc.bias.ols.update=perc.bias.ols.update[!is.na(perc.bias.ols.update)]
    lines(N,perc.bias.spr.update,lty=6,lwd=1,col="red")
    lines(N,perc.bias.ols.update,lty=3,lwd=1,col="red")
  }
  ######
  axis(side=1,at=seq(2000,8000,by=2000),tick=T,line=0,labels=rep("",4))
  axis(side=1,at=seq(2000,8000,by=2000),
       tick=F,line=-0.5,cex.axis=mycex,
       labels=seq(2000,8000,by=2000))
  mtext(side=1,line=1.5,text="n")
  mtext(side=1,line=1.65,text=bquote("("*alpha*" = 0.7)"),at=at2)
  ####
  axis(side=2,at=seq(ymin,ymax,length.out=5),tick=T,line=0,labels=rep("",5))
  axis(side=2,at=seq(ymin,ymax,length.out=5),
       tick=F,line=-0.5,cex.axis=mycex,
       labels=round(seq(ymin,ymax,length.out=5)*1e2,1))
  mtext(side=2,line=1.5,text=ylabname,cex=mycex)
  if(addlabel){
    legend("bottomright",bty="n",cex=mycex,
           lty=c(1,5,6,3)[c(1,3,2,4)],
           lwd=c(2,2,1,1)[c(1,3,2,4)],
           col=c("black","black","red","red")[c(1,3,2,4)],
           legend=c(expression(paste("iSphereMAP step II: ",widehat(Pi)^"[2]")),
                    expression(paste("MT method")),
                    expression(paste("iSphereMAP step IV: ",widehat(Pi)^"[3]")),
                    expression(paste("MT method step IV")))[c(1,3,2,4)])
    # legend("bottomright",bty="n",cex=mycex,
    #        lty=c(1,5,6,3)[c(1,3,2,4)],
    #        lwd=c(2,2,1,1)[c(1,3,2,4)],
    #        col=c("black","black","red","red")[c(1,3,2,4)],
    #        legend=c(expression(paste("iSphereMAP ",widehat(Pi)^"[2]")),
    #                 expression(paste("MT method ",widehat(Pi)^"[2]")),
    #                 expression(paste("iSphereMAP ",widehat(Pi)^"[3]")),
    #                 expression(paste("MT method ",widehat(Pi)^"[3]")))[c(1,3,2,4)])
  }
  
}
plot.rslt.only1to1=function(rslt,ylabname="1-to-1 Match Rate",addlabel=T){
  perc.bias.spr = rslt$perc.bias.spr
  perc.bias.ols = rslt$perc.bias.ols
  perc.bias.spr.update=rslt$perc.bias.spr#.update
  perc.bias.ols.update=rslt$perc.bias.ols#.update
  N=rslt$N
  if(myloess==T){
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
  }else{
    ymin=min(na.rm=T,perc.bias.spr,perc.bias.ols,perc.bias.spr.update,perc.bias.ols.update)
    ymax=max(na.rm=T,perc.bias.spr,perc.bias.ols,perc.bias.spr.update,perc.bias.ols.update)
    N=N[!is.na(N)]
    perc.bias.spr=perc.bias.spr[!is.na(perc.bias.spr)]
    perc.bias.ols=perc.bias.ols[!is.na(perc.bias.ols)]
    plot(N,perc.bias.spr,
         xlab="",ylab="",axes=F,
         type="l",lty=1,lwd=2,
         ylim=c(ymin,ymax));box()
    lines(N,perc.bias.ols,lty=5,lwd=2)
  }
  # #### add update rslt
  # if(myloess==T){
  #   N=N[!is.na(N)]
  #   lines(N,predict(lo.spr.update),lty=6,lwd=1,col="red")
  #   lines(N,predict(lo.ols.update),lty=3,lwd=1,col="red")
  # }else{
  #   N=N[!is.na(N)]
  #   perc.bias.spr.update=perc.bias.spr.update[!is.na(perc.bias.spr.update)]
  #   perc.bias.ols.update=perc.bias.ols.update[!is.na(perc.bias.ols.update)]
  #   lines(N,perc.bias.spr.update,lty=6,lwd=1,col="red")
  #   lines(N,perc.bias.ols.update,lty=3,lwd=1,col="red")
  # }
  ######
  axis(side=1,at=seq(2000,8000,by=2000),tick=T,line=0,labels=rep("",4))
  axis(side=1,at=seq(2000,8000,by=2000),
       tick=F,line=-0.5,cex.axis=mycex,
       labels=seq(2000,8000,by=2000))
  mtext(side=1,line=1.5,text="n")
  mtext(side=1,line=1.65,text=bquote("("*alpha*" = 0.7)"),at=at2)
  ####
  axis(side=2,at=seq(ymin,ymax,length.out=5),tick=T,line=0,labels=rep("",5))
  axis(side=2,at=seq(ymin,ymax,length.out=5),
       tick=F,line=-0.5,cex.axis=mycex,
       labels=round(seq(ymin,ymax,length.out=5)*1e2,1))
  mtext(side=2,line=1.5,text=ylabname,cex=mycex)
  if(addlabel){
    legend("bottomright",bty="n",cex=mycex,
           lty=c(1,5,6,3)[c(1,2)],
           lwd=c(2,2,1,1)[c(1,2)],
           col=c("black","black","red","red")[c(1,2)],
           legend=c(expression(paste("iSphereMAP step II: ",widehat(Pi)^"[2]")),
                    expression(paste("MT method")))
    )
  }
  
}

compare.rslt=function(rslt0,rslt1,myupdate){
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
  if(myupdate==F){
    perc.bias.spr=perc.bias.spr0
    perc.bias.ols=perc.bias.ols0 
    perc.bias.spr.update=perc.bias.spr1   
    perc.bias.ols.update=perc.bias.ols1
  }else{
    perc.bias.spr=perc.bias.spr.update0
    perc.bias.ols=perc.bias.ols.update0 
    perc.bias.spr.update=perc.bias.spr.update1   
    perc.bias.ols.update=perc.bias.ols.update1
  }
  (N=(N1+N0)/2)
  rslt=list(perc.bias.spr=perc.bias.spr,perc.bias.ols=perc.bias.ols,
            perc.bias.spr.update=perc.bias.spr.update,perc.bias.ols.update=perc.bias.ols.update,N=N)
  return(rslt)
}
wrong_grp_info=F
compute_W_k=F
add_mismatch_1to1=F
perm_only=T
MVN=F
pdf(paste0("percent_bias_Pi_Perm_only_T_",Sys.time(),".pdf"),h=3.2*1.25,w=4*1.25)
par(mfrow=c(1,1),mar=c(0,0,0,0),mai=c(0.8,0.6,0.5,0.15),omi=c(0,0,0,0))
##################################################################
### N vs bias
##################################################################
rslt=read.rslt(allargs.i=5,name1="match.rate",name2="match.rate.nogrp")
# plot.rslt(rslt,ylabname="1-to-1 Match Rate %",addlabel=T)
plot.rslt.only1to1(rslt,ylabname="1-to-1 Match Rate %",addlabel=T)
# rslt=read.rslt(allargs.i=5,name1="C.error_F",name2="C.error_F.nogrp")
# rslt=plot.rslt(rslt,ylabname=bquote("MSE of 1-to-m Weight (x"*10^{-2}*")"),addlabel=F)
# ##################################################################
# ### Compare allargs FFFFF to FTFFF
# ##################################################################
# rslt0=read.rslt(allargs.i=1,name1="match.rate",name2="match.rate.nogrp")
# rslt1=read.rslt(allargs.i=5,name1="match.rate",name2="match.rate.nogrp")
# rslt=compare.rslt(rslt0,rslt1,myupdate=F)
# plot.rslt(rslt,ylabname="1-to-1 Match Rate %",addlabel=F);abline(h=0)
# legend("bottomright",bty="n",cex=mycex,
#        lty=c(1,5,6,3)[c(1,3,2,4)],
#        lwd=c(2,2,1,1)[c(1,3,2,4)],
#        col=c("black","black","red","red")[c(1,3,2,4)],
#        title=expression(paste("Step II: estimation of ",Pi)),
#        legend=c("iSphereMAP w/ both 1-1 and 1-m mappings",
#                 "MT method w/ both 1-1 and 1-m mappings",
#                 "iSphereMAP w/ only 1-1 mapping",
#                 "MT method w/ only 1-1 mapping")[c(1,3,2,4)])
# rslt0=read.rslt(allargs.i=1,name1="C.error_F",name2="C.error_F.nogrp")
# rslt1=read.rslt(allargs.i=5,name1="C.error_F",name2="C.error_F.nogrp")
# rslt=compare.rslt(rslt0,rslt1,myupdate=F)
# plot.rslt(rslt,ylabname=bquote("MSE of 1-to-m Weight (x"*10^{-2}*")"),addlabel=F)
dev.off()



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
  perc.bias.spr=perc.bias.spr1-perc.bias.spr0
  perc.bias.ols=perc.bias.ols1-perc.bias.ols0
  perc.bias.spr.update=perc.bias.spr.update1-perc.bias.spr.update0
  perc.bias.ols.update=perc.bias.ols.update1-perc.bias.ols.update0
  (N=(N1+N0)/2)
  rslt=list(perc.bias.spr=perc.bias.spr,perc.bias.ols=perc.bias.ols,
            perc.bias.spr.update=perc.bias.spr.update,perc.bias.ols.update=perc.bias.ols.update,N=N)
  return(rslt)
}
