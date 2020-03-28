rm(list=ls())
library(Matrix)
library(wordspace)
vec.size = 300; fine_grp=T
load("ICD9_embeddings_from_VHA_PHS.RData")
c(nrow(X.nonzero),nrow(Z.nonzero))
join.names = intersect(row.names(X.nonzero),row.names(Z.nonzero))
X.nonzero = X.nonzero[match(join.names,row.names(X.nonzero)),]
Z.nonzero = Z.nonzero[match(join.names,row.names(Z.nonzero)),]
grp.info = phewas.grp.X$X.phewas
grp.info = grp.info[match(join.names,grp.info$icd),]
length(unique(grp.info$g.index))
X = X.test = X.train = as.matrix(X.nonzero[match(join.names,row.names(X.nonzero)),])
Y = Y.test = Y.train = as.matrix(Z.nonzero[match(join.names,row.names(Z.nonzero)),])
### order X Y and grp.info by g.index
order.g.index = order(grp.info$g.index)
X = X.test = X.train = X[order.g.index,]
Y = Y.test = Y.train = Y[order.g.index,]
grp.info = grp.info[order.g.index,]

dim(X);dim(Y);dim(grp.info)
pcomp=T
normalizeall=T
if(normalizeall==T){
  X.norms = rowNorms(as.matrix(X))
  if(length(which(X.norms!=1))){
    X = X/rowNorms(as.matrix(X))
    Y = Y/rowNorms(as.matrix(Y))
  }
}
######
if(pcomp==T){
  Y = prcomp(Y)$x[,1:3];Y/rowNorms((Y))
  X = prcomp(X)$x[,1:3];X/rowNorms((X))
}


tab=table(grp.info$g.index)
table(tab);
used_in_paper=T
select.grpsize=median(as.numeric(names(table(tab))))
##infact shoud be median(as.numeric((table(tab))))
ntab=c(floor(select.grpsize),ceiling(select.grpsize))
select.grp=as.numeric(names(tab)[which(tab%in%ntab)])
ind=grp.info$g.index%in%select.grp
colors <- rainbow(length(select.grp))
colors <- colors[as.numeric(as.factor(grp.info$g.index[ind]))]
t=grp.info$g.index[grp.info$g.index%in%select.grp]
plot(t,type="p",col=colors)
length(unique(grp.info$g.index))
#####clean for plot
#####choose to plot either X or Y!!!!!!!!!!!!!
dat.save = cbind(Y[,1:3],grp.info$Grp,grp.info)
colnames(dat.save)[1:4]=c("x","y","z","grp")
dat = as.data.frame(dat.save)
dat=dat[ind,]
dat[,1]=x=as.numeric(as.character(dat[,1]))
dat[,2]=y=as.numeric(as.character(dat[,2]))
dat[,3]=z=as.numeric(as.character(dat[,3]))
####plot
library(rgl)
open3d() 
spheres3d(0,0,0,radius=1,lit=F,size=20,
          color="grey",back="lines",front="lines", alpha = 0.2)
spheres3d(dat[,1:3]/rowNorms(as.matrix(dat[,1:3])), 
          col=colors,radius=0.02,lit=T,size=3)
### rotate till good view, then save as pdf
rgl.postscript(paste0("X",Sys.time(),".pdf"),"pdf")
