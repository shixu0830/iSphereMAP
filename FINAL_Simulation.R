#########################################################################################
#########################################################################################
### Simulation for Figures 2-3,and Figures F.1-F.4                                    ###
### Fig 2-3: wrong_grp_info=F; add_mismatch_1to1=F; perm_only=F; MVN=F; mykappa=150   ###
### Fig F.1: wrong_grp_info=T; add_mismatch_1to1=F; perm_only=F; MVN=F; mykappa=150   ###
### Fig F.2: wrong_grp_info=F; add_mismatch_1to1=F; perm_only=T; MVN=F; mykappa=150   ###
### Fig F.3: wrong_grp_info=F; add_mismatch_1to1=F; perm_only=F; MVN=F; mykappa=3000  ###
### Fig F.4: wrong_grp_info=F; add_mismatch_1to1=T; perm_only=T; MVN=F; mykappa=150   ###
### Contact: shixu@umich.edu                                                          ###
#########################################################################################
#########################################################################################
rm(list=ls())
allseeds=1:10
library(MASS)
library(Matrix)
library(movMF) ##simulate von Mises-Fisher data
my.filepath = ""
p=300
mykappa = 150
nlambda=5                         ## search from how many lambdas in Cross-validation
allk = seq(100,1700,by=200)       ## a range of number of groups
allsparse=seq(0,0.4,by=0.05)      ## a range of sparsity values (alpha), n_mis = n^alpha
#### args: select some N and some sparsity
args = commandArgs(TRUE)
if(length(args)>0){
  tt = as.numeric(args[[1]])      ## index for allsparse (the alpha in n^alpha)
  kk = as.numeric(args[[2]])      ## index for allk (number of groups)
  wrong_grp_info=(args[[3]]==1)   ## whether grp-structure is misspecified
  # compute_W_k=(args[[4]]==1)    ## whether estimate group-wise rotation W_k, always F
  add_mismatch_1to1=(args[[5]]==1)## whether use all 1-1 data (including mismatched) in step III
  perm_only=(args[[6]]==1)        ## whether data only include 1-1 mismatch, no 1-m mismatch
  MVN=(args[[7]]==1)              ## whether generate normalized MVN data or vMF data
  myseed=as.numeric(args[[8]])    ## seed number
}
source(paste0(my.filepath,"FINAL_Simulation_functions.R"))
## read in true W ("BetaBALL") and grouping info ("grp.info") from real data
load(paste0(my.filepath,"FINALsummstat.RData")) 

####################################################################################
####################################################################################
####### simu setting
####################################################################################
####################################################################################
##### choose a number of groups and a sparsity
grp.select = 1:allk[kk]
sparsity = allsparse[tt]
##### W.true is from real data "FINALsummstat.RData"
W.true = BetaBALL
##### "grp.info" is the PheCode group that each code belongs to
## note that the "real data" defined by grp.select may not cover all phewas groups
grp.info = grp.info[grp.info$g.index%in%c(grp.select),]
### "grp.info2" is for the MT method -- no grouping
grp.info2 = grp.info; grp.info2$g.index=1
### save dimensions
N = nrow(grp.info)
ugrp = unique(grp.info$g.index)
n = length(ugrp)
c(N,n,p)

####################################################################################
####################################################################################
####### simu iterations
####################################################################################
####################################################################################
for(myseed in allseeds){
  for(tt in 1:9){## index for allsparse (the alpha in n^alpha)
    rslt=try(run(kk,tt,myseed,MVN,add_mismatch_1to1,
                 perm_only,wrong_grp_info),silent=TRUE)
    if(inherits(rslt, "try-error")){
      Wrslt=Pirslt=Pirslt.update=NA
    }else{
      Wrslt=rslt$Wrslt
      Pierr=rslt$Pierr
      Pierr.update=rslt$Pierr.update
    }
    save(Wrslt,Pierr,Pierr.update,file=
           paste0(my.filepath,"rslt",
                  paste0("_",kk,"_",tt,"_",
                         wrong_grp_info,"_",
                         add_mismatch_1to1,"_",
                         perm_only,"_",
                         MVN,"_",
                         myseed),".RData"))
    print(paste0(tt,"_",myseed))
  }#endtt
}#endmyseed
