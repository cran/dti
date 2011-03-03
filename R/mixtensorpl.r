#
#  Model without isotropic part
#
mfunpl0 <- function(par,siq,grad,pen=1e2){
#
#   evaluate rss for Mixtensor-model (without isotropic component)
#
lpar <- length(par)
m <- (lpar-1)/2
ngrad <- dim(grad)[1]
.Fortran("mfunpl0",as.double(par),#par(lpar)
                as.double(siq),#siq(ngrad)
                as.double(t(grad)),#grad(3,ngrad)
                as.integer(m),#number of components
                as.integer(lpar),#number of parameters
                as.integer(ngrad),#number of gradients
                as.double(pen),#penalty for negative weights
                double(ngrad*m),#z(ngrad,m) working array
                double(ngrad),#w(ngrad) working array
                erg = double(1),#residual sum of squares
                PACKAGE="dti")$erg
}
gmfunpl0 <- function(par,siq,grad,pen=1e2){
#
#   evaluate rss for Mixtensor-model
#
lpar <- length(par)
m <- (lpar-1)/2
ngrad <- dim(grad)[1]
.Fortran("mfunpl0g",
         as.double(par),#par(lpar)
         as.double(siq),#s(n)
         as.double(t(grad)),#g(3,n)
         as.integer(m),
         as.integer(lpar),
         as.integer(ngrad),
         double(3*m),# d(3,m)
         double(ngrad*m),# z(n,m)
         double(m*m),# v(m,m)
         double(m),# w(m)
         double(ngrad*m),# dkgj(n,m)
         double(ngrad*m),# dkgj2(n,m)
         double(ngrad*m),# ddkdphig(n,m)
         double(ngrad*m),# ddkdetag(n,m)
         double(3*m),# dddphi(3,m)
         double(3*m),# dddeta(3,m)
         double(m*m),# dvdth(m,m)
         double(m*m*m),# dvdphi(m,m,m)
         double(m*m*m),# dvdeta(m,m,m)
         double(ngrad*m*lpar),# dzdpars(n,m,lpar)
         double(m*lpar),# dwdpars(m,lpar)
         double(ngrad*m),# zs(n,m)
         double(ngrad*m),# work1(n,m)
         double(ngrad*m),# work2(n,m)
         double(ngrad),# scopy(n)
         as.double(pen),# pen
         dfdpar=double(lpar),# dfdpar(lpar)
         PACKAGE="dti")$dfdpar
}
mfunplwghts0 <- function(par,siq,grad,pen=1e2){
#
#   get weights for Mixtensor-model (without isotropic component) and extract parameters 
#
lpar <- length(par)
m <- (lpar-1)/2
ngrad <- dim(grad)[1]
par[1] <- max(0,par[1])
w<-.Fortran("mfunpl0",as.double(par),#par(lpar)
                as.double(siq),#siq(ngrad)
                as.double(t(grad)),#grad(3,ngrad)
                as.integer(m),#number of components
                as.integer(lpar),#number of parameters
                as.integer(ngrad),#number of gradients
                as.double(pen),#penalty for negative weights
                double(ngrad*m),#z(ngrad,m) working array
                w = double(ngrad),#w(ngrad) working array
                double(1),#residual sum of squares
                PACKAGE="dti")$w[1:m]
#                if(any(is.na(w)|abs(w)>1e10)){
#                   cat("par",par,"\n siq",siq,"\n m",m,"\n lpar",lpar,"\n pen",
#                        pen,"\n w",w,"\n")
#                }
erg<-.Fortran("mfunpl0w",as.double(par),#par(lpar)
                as.double(pmax(w,0)),
                as.double(siq),#siq(ngrad)
                as.double(t(grad)),#grad(3,ngrad)
                as.integer(m),#number of components
                as.integer(lpar),#number of parameters
                as.integer(ngrad),#number of gradients
                double(ngrad*m),#z(ngrad,m) working array
                erg = double(1),#residual sum of squares
                PACKAGE="dti")$erg
           o <- order(w,decreasing=TRUE)
           ord <- sum(w>0)
           if(ord<m){
              o <- o[1:ord]
           }
           sw <- sum(w[w>0])
           lev <- c(par[1],-log(sw))
           if(ord>0){
           mix <- w[o]/sw
           } else {
           mix <- NULL
           } 
#           if(par[1]<0) cat("par",par,"\n")
           or <- matrix(par[2:lpar],2,m)[,o,drop=FALSE]
           or[1,or[1,]<0] <- or[1,or[1,]<0]+pi
           or[1,or[1,]>pi] <- or[1,or[1,]>pi]-pi
           or[2,or[2,]<0] <- or[2,or[2,]<0]+2*pi
           or[2,or[2,]>2*pi] <- or[2,or[2,]>2*pi]-2*pi
           par <- c(par[1],or[,1:ord])
list(ord=ord,lev=lev,mix=mix,orient=or,par=par,value=erg)
}
#
#  Model with isotropic part
#
mfunpl1 <- function(par,siq,grad,pen=1e2){
#
#   evaluate rss for Mixtensor-model (without isotropic component)
#
lpar <- length(par)
m <- (lpar-1)/2
mp1 <- m+1
ngrad <- dim(grad)[1]
z<-.Fortran("mfunpl1",as.double(par),#par(lpar)
                as.double(siq),#siq(ngrad)
                as.double(t(grad)),#grad(3,ngrad)
                as.integer(m),#number of components
                as.integer(mp1),#number of components+1
                as.integer(lpar),#number of parameters
                as.integer(ngrad),#number of gradients
                as.double(pen),#penalty for negative weights
                double(ngrad*mp1),#z(ngrad,m+1) working array
                double(ngrad),#w(ngrad) working array
                erg = double(1),#residual sum of squares
                PACKAGE="dti")$erg
#         cat("par:",signif(par,3),"fw:",signif(z,3),"\n")
z
}
gmfunpl1 <- function(par,siq,grad,pen=1e2){
#
#   evaluate rss for Mixtensor-model
#
lpar <- length(par)
m <- (lpar-1)/2
mp1 <- m+1
ngrad <- dim(grad)[1]
z<-.Fortran("mfunpl1g",
         as.double(par),#par(lpar)
         as.double(siq),#s(n)
         as.double(t(grad)),#g(3,n)
         as.integer(m),
         as.integer(mp1),
         as.integer(lpar),
         as.integer(ngrad),
         double(3*m),# d(3,m)
         double(ngrad*mp1),# z(n,mp1)
         double(mp1*mp1),# v(mp1,mp1)
         double(mp1),# w(mp1)
         double(ngrad*m),# dkgj(n,m)
         double(ngrad*m),# dkgj2(n,m)
         double(ngrad*m),# ddkdphig(n,m)
         double(ngrad*m),# ddkdetag(n,m)
         double(3*m),# dddphi(3,m)
         double(3*m),# dddeta(3,m)
         double(mp1*mp1),# dvdth(mp1,mp1)
         double(mp1*mp1*m),# dvdphi(mp1,mp1,m)
         double(mp1*mp1*m),# dvdeta(mp1,mp1,m)
         double(ngrad*mp1*lpar),# dzdpars(n,mp1,lpar)
         double(mp1*lpar),# dwdpars(mp1,lpar)
         double(ngrad*mp1),# zs(n,mp1)
         double(ngrad*mp1),# work1(n,mp1)
         double(ngrad*mp1),# work2(n,mp1)
         double(ngrad),# scopy(n)
         as.double(pen),# pen
         dfdpar=double(lpar),# dfdpar(lpar)
         PACKAGE="dti")$dfdpar
#         cat("par:",signif(par,3),"\n")
#         cat("gradient:",signif(z,3),"\n")
         z
}
mfunplwghts1 <- function(par,siq,grad,pen=1e2){
#
#   get weights for Mixtensor-model (without isotropic component) and extract parameters 
#
lpar <- length(par)
m <- (lpar-1)/2
mp1 <- m+1
ngrad <- dim(grad)[1]
w<-.Fortran("mfunpl1",as.double(par),#par(lpar)
                as.double(siq),#siq(ngrad)
                as.double(t(grad)),#grad(3,ngrad)
                as.integer(m),#number of components
                as.integer(mp1),#number of components+1
                as.integer(lpar),#number of parameters
                as.integer(ngrad),#number of gradients
                as.double(pen),#penalty for negative weights
                double(ngrad*mp1),#z(ngrad,m) working array
                w = double(ngrad),#w(ngrad) working array
                double(1),#residual sum of squares
                PACKAGE="dti")$w[1:mp1]
           erg <- .Fortran("mfunpl1w",as.double(par),#par(lpar)
                as.double(pmax(w,0)),#par(lpar) 
                as.double(siq),#siq(ngrad)
                as.double(t(grad)),#grad(3,ngrad)
                as.integer(m),#number of components
                as.integer(mp1),#number of components+1
                as.integer(lpar),#number of parameters
                as.integer(ngrad),#number of gradients
                double(ngrad*mp1),#z(ngrad,m+1) working array
                erg = double(1),#residual sum of squares
                PACKAGE="dti")$erg
           w0 <- w[1]
           w <- w[-1]
           o <- order(w,decreasing=TRUE)
           ord <- sum(w>0)
           if(ord<m){
              o <- o[1:ord]
           }
           problem <- FALSE
           if(sum(w[w>0])+w0<=0){
#           cat("w0",w0,"w",w,"\n")
#           cat("par",par,"\n")
           problem <- TRUE
           }
           sw <- sum(w[w>0])+max(w0,0)
           lev <- c(par[1],-log(sw))
           if(ord>0){
           mix <- w[o]/sw
           } else {
           mix <- NULL
           } 
           mix0 <- w0/sw
           or <- matrix(par[2:lpar],2,m)[,o,drop=FALSE]
           or[1,or[1,]<0] <- or[1,or[1,]<0]+pi
           or[1,or[1,]>pi] <- or[1,or[1,]>pi]-pi
           or[2,or[2,]<0] <- or[2,or[2,]<0]+2*pi
           or[2,or[2,]>2*pi] <- or[2,or[2,]>2*pi]-2*pi
           par <- c(par[1],or[,1:ord])
           if(problem) cat("ord",ord,"lev",lev,"mix",mix,"mix0",mix0,"or",or,"par",par,"\n")
list(ord=ord,lev=lev,mix=mix,mix0=mix0,orient=or,par=par,value=erg)
}
#
#   Initial estimates
#
selisample <- function(ngrad,maxcomp,nguess,dgrad,maxc){
saved.seed <- .Random.seed
set.seed(1)
isample <- matrix(sample(ngrad,maxcomp*nguess,replace=TRUE),maxcomp,nguess)
ind <- rep(TRUE,nguess)
if(maxcomp>1){
ind <- .Fortran("selisamp",
                as.integer(isample),
                as.integer(nguess),
                as.integer(maxcomp),
                as.double(dgrad),
                as.integer(dim(dgrad)[1]),
                ind = logical(nguess),
                as.double(maxc),
                DUPL=FALSE,
                PACKAGE="dti")$ind 
#for(i in 1:nguess) for(j in 1:(maxcomp-1)) for(k in (j+1):maxcomp){
#if(dgrad[isample[j,i],isample[k,i]]>maxc) ind[i] <- FALSE
#}
.Random.seed <- saved.seed
}
isample[,ind]
}

paroforient <- function(dir){
  theta <- acos(dir[3])
  sth <- sin(theta)
  phi <- 0
  if(sth<1e-8) {
    theta <- 0
  } else {
    z <- dir[1]/sth
    if(abs(z)>=1) {
      phi <- if(z<0) 0 else pi
    } else {
      phi <- acos(z)*sign(dir[2])
    }
    if(phi < 0) phi <- phi+2*pi
  }
  c(theta, phi)
}

getsiind3 <- function(si,mask,sigma2,grad,vico,th,indth,ev,fa,andir,maxcomp=3,maxc=.866,nguess=100){
# assumes dim(grad) == c(ngrad,3)
# assumes dim(si) == c(n1,n2,n3,ngrad)
# SO removed
ngrad <- dim(grad)[1]
nvico <- dim(vico)[1]
ddim <- dim(fa)
nsi <- dim(si)[4]
dgrad <- matrix(abs(grad%*%t(vico)),ngrad,nvico)
dgrad <- dgrad/max(dgrad)
dgradi <- matrix(abs(vico%*%t(vico)),nvico,nvico)
dgradi <- dgradi/max(dgradi)
nth <- length(th)
nvoxel <- prod(ddim)
landir <- fa>.3
landir[is.na(landir)] <- FALSE
if(any(is.na(andir))) {
cat(sum(is.na(andir)),"na's in andir")
andir[is.na(andir)]<-sqrt(1/3)
}
if(any(is.na(landir))) {
cat(sum(is.na(landir)),"na's in landir")
landir[is.na(landir)]<-0
}
if(any(is.na(fa))) {
cat(sum(is.na(fa)),"na's in fa")
fa[is.na(fa)]<-0
}
iandir <- array(.Fortran("iandir",
                   as.double(t(vico)),
                   as.integer(nvico),
                   as.double(andir),
                   as.integer(nvoxel),
                   as.logical(landir),
                   iandir=integer(2*prod(ddim)),
                   DUPL=FALSE,
                   PACKAGE="dti")$iandir,c(2,nvoxel))
isample0 <- selisample(nvico,maxcomp,nguess,dgradi,maxc)
if(maxcomp>1) isample1 <- selisample(nvico,maxcomp-1,nguess,dgradi,maxc)
if(maxcomp==1) isample1 <- sample(ngrad, nguess, replace = TRUE)
#
#  eliminate configurations with close directions 
#
# this provides configurations of initial estimates with minimum angle between 
# directions > acos(maxc)
nvoxel <- prod(dim(si)[-4])
cat("using ",nguess,"guesses for initial estimates\n")
siind <- matrix(as.integer(0),maxcomp+2,nvoxel)
krit <- numeric(nvoxel)
# first voxel with fa<.3
cat(sum(mask&!landir),"voxel with small FA\n")
nguess <- length(isample0)/maxcomp
z <- .Fortran("getsii30",
         as.double(aperm(si,c(4,1:3))),
         as.double(sigma2),
         as.integer(nsi),
         as.integer(nvoxel),
         as.integer(maxcomp),
         as.double(dgrad),
         as.integer(nvico),
         as.double(th),
         as.integer(nth),
         as.integer(indth),
         double(ngrad*nvico),
         as.integer(isample0),
         as.integer(nguess),
         double(nsi),
         double(nsi*(maxcomp+2)),
         siind=integer((maxcomp+2)*nvoxel),
         krit=double(nvoxel),
         as.integer(maxcomp+2),
         as.logical(mask&!landir),
         PACKAGE="dti")[c("siind","krit")]
dim(z$siind) <- c(maxcomp+2,nvoxel)
siind[,!landir] <- z$siind[,!landir]
krit[!landir] <- z$krit[!landir]
# now voxel where first tensor direction seems important
if(maxcomp >0){
cat(sum(mask&landir),"voxel with distinct first eigenvalue \n")
nguess <- if(maxcomp>1) length(isample1)/(maxcomp-1) else length(isample1)
z <- .Fortran("getsii31",
         as.double(aperm(si,c(4,1:3))),
         as.double(sigma2),
         as.integer(nsi),
         as.integer(nvoxel),
         as.integer(maxcomp),
         as.double(dgrad),
         as.integer(nvico),
         as.integer(iandir[1,]),
         as.double(th),
         as.integer(nth),
         as.integer(indth),
         double(ngrad*nvico),
         as.integer(isample1),
         as.integer(nguess),
         double(nsi),
         double(nsi*(maxcomp+2)),
         siind=integer((maxcomp+2)*nvoxel),
         krit=double(nvoxel),
         as.integer(maxcomp+2),
         as.logical(mask&landir),
         as.double(dgradi),
         as.double(maxc),
         PACKAGE="dti")[c("siind","krit")]
dim(z$siind) <- c(maxcomp+2,nvoxel)
siind[,landir] <- z$siind[,landir]
krit[landir] <- z$krit[landir]
}
failed <- (krit^2/ngrad) > (sigma2-1e-10)
if(any(failed[mask])){
print((krit[mask])[failed[mask]])
print(((1:prod(dim(si)[1:3]))[mask])[failed[mask]])
print(sum(failed[mask]))
}
list(siind=array(siind,c(maxcomp+2,dim(si)[-4])),
     krit=array(krit,dim(si)[-4]))
}

getsiind2 <- function(si,mask,sigma2,grad,vico,th,maxcomp=3,maxc=.866,nguess=100){
# assumes dim(grad) == c(ngrad,3)
# assumes dim(si) == c(n1,n2,n3,ngrad)
# SO removed
ngrad <- dim(grad)[1]
nvico <- dim(vico)[1]
nsi <- dim(si)[4]
dgrad <- matrix(abs(grad%*%t(vico)),ngrad,nvico)
dgrad <- dgrad/max(dgrad)
dgradi <- matrix(abs(vico%*%t(vico)),nvico,nvico)
dgradi <- dgradi/max(dgradi)
nth <- length(th)
isample <- selisample(nvico,maxcomp,nth*nguess,dgradi,maxc)
#
#  eliminate configurations with close directions 
#
if(maxcomp>1) { 
nguess <- trunc(dim(isample)[2]/nth) 
isample <- array(isample[,1:(nguess*nth)],c(maxcomp,nguess,nth))
} 
# this provides configurations of initial estimates with minimum angle between 
# directions > acos(maxc)
cat("using ",nguess,"guesses for initial estimates\n")
nvoxel <- prod(dim(si)[-4])
siind <- .Fortran("getsiin2",
         as.double(aperm(si,c(4,1:3))),
         as.double(sigma2),
         as.integer(nsi),
         as.integer(nvoxel),
         as.integer(maxcomp),
         as.double(dgrad),
         as.integer(nvico),
         as.double(th),
         as.integer(nth),
         double(ngrad*nvico),
         as.integer(isample),
         as.integer(nguess),
         double(nsi),
         double(nsi*(maxcomp+2)),
         siind=integer((maxcomp+2)*nvoxel),
         krit=double(prod(dim(si)[1:3])),
         as.integer(maxcomp+2),
         as.logical(mask),
         PACKAGE="dti")[c("siind","krit")]
failed <- (siind$krit^2/ngrad) > (sigma2-1e-10)
if(any(failed[mask])){
print((siind$krit[mask])[failed[mask]])
print(((1:prod(dim(si)[1:3]))[mask])[failed[mask]])
print(sum(failed[mask]))
}
list(siind=array(siind$siind,c(maxcomp+2,dim(si)[-4])),
     krit=array(siind$krit,dim(si)[-4]))
}

dwiMixtensor <- function(object, ...) cat("No dwiMixtensor calculation defined for this class:",class(object),"\n")

setGeneric("dwiMixtensor", function(object,  ...) standardGeneric("dwiMixtensor"))

setMethod("dwiMixtensor","dtiData",function(object, maxcomp=3,  p=40, method="mixtensor", reltol=1e-6, maxit=5000,ngc=1000, optmethod="BFGS", nguess=100*maxcomp^2,msc="BIC",pen=NULL){
#
#  uses  S(g)/s_0 = w_0 exp(-l_1) +\sum_{i} w_i exp(-l_2-(l_1-l_2)(g^T d_i)^2)
#
#  
#
  set.seed(1)
  if(is.null(pen)) pen <- 1e2
  theta <- .5
  maxc <- .866
  args <- sys.call(-1)
  args <- c(object@call,args)
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  ngrad0 <- ngrad - ns0
  if(5*(1+3*maxcomp)>ngrad0){
#     maxcomp <- max(1,trunc((ngrad0-5)/15))
     cat("Maximal number of components reduced to", maxcomp,"due to insufficient
          number of gradient directions\n")
  }
#
#  First tensor estimates to generate eigenvalues and -vectors
#
  prta <- Sys.time()
  cat("Start tensor estimation at",format(prta),"\n")
  tensorobj <- dtiTensor(object)
  cat("Start evaluation of eigenstructure at",format(Sys.time()),"\n")
  z <- .Fortran("dtieigen",
                as.double(tensorobj@D),
                as.integer(ddim[1]),
                as.integer(ddim[2]),
                as.integer(ddim[3]),
                as.logical(tensorobj@mask),
                fa=double(prod(ddim)),
                ev=double(3*prod(ddim)),
                andir=double(6*prod(ddim)),
                DUP=FALSE,
                PACKAGE="dti")[c("fa","ev","andir")]
  rm(tensorobj)
  gc()
  fa <- array(z$fa,ddim)
  ev <- array(z$ev,c(3,ddim))
  andir <- array(z$andir,c(3,2,ddim))
  rm(z)
  gc()
  nth <- 11
  th <- ev[1,,,] - (ev[2,,,]+ev[3,,,])/2
  falevel <- min(quantile(fa[fa>0],.75),.3)
  rth <- quantile(th[fa>=falevel],c(.1,.99))
  th[th<rth[1]] <- rth[1]
  th[th>rth[2]] <- rth[2]
  if(diff(rth)>0){
     indth <- trunc((nth-1)*(th-rth[1])/diff(rth)+1)
     th <- seq(rth[1],rth[2],length=nth)
  } else {
     th <- rep(max(th),nth)
     indth <- rep(1,length(th))
  }
  cat("Start search outlier detection at",format(Sys.time()),"\n")
#
#  replace physically meaningless S_i by mena S_0 values
#
  z <- .Fortran("outlier",#misc.f 
                as.double(object@si),
                as.integer(prod(ddim)),
                as.integer(ngrad),
                as.logical((1:ngrad)%in%s0ind),
                as.integer(ns0),
                si=integer(prod(ddim)*ngrad),
                index=integer(prod(ddim)),
                lindex=integer(1),
                DUP=FALSE,
                PACKAGE="dti")[c("si","index","lindex")]
  cat("End search outlier detection at",format(Sys.time()),"\n")
  si <- array(as.integer(z$si),c(ddim,ngrad))
  index <- if(z$lindex>0) z$index[1:z$lindex] else numeric(0)
  rm(z)
  gc()
  cat("Start generating auxiliary objects",format(Sys.time()),"\n")
#
#  compute mean S_0, s_i/S_0 (siq), var(siq) and mask
#
  z <- .Fortran("sweeps0",# mixtens.f
                as.integer(si[,,,-s0ind,drop=FALSE]),
                as.integer(si[,,,s0ind,drop=FALSE]),
                as.integer(ddim[1]),
                as.integer(ddim[2]),
                as.integer(ddim[3]),
                as.integer(ns0),
                as.integer(ngrad0),
                as.integer(object@level),
                siq=double(prod(ddim[1:3])*ngrad0),
                s0=double(prod(ddim[1:3])),
                vsi=double(prod(ddim[1:3])),
                mask=logical(prod(ddim[1:3])),
                DUPL=FALSE,
                PACKAGE="dti")[c("siq","s0","vsi","mask")]
  rm(si)
  s0 <- array(z$s0,ddim[1:3])
  siq <- array(z$siq,c(ddim[1:3],ngrad0))
  sigma2 <- array(z$vsi,ddim[1:3])
  mask <- array(z$mask,ddim[1:3])
  rm(z)
  gc()
  npar <- if(method=="mixtensor") 1+3*(0:maxcomp) else c(1,2+3*(1:maxcomp))
#
#   compute penalty for model selection, default BIC
#
  penIC <- switch(msc,"AIC"=2*npar/ngrad0,"BIC"=log(ngrad0)*npar/ngrad0,
                  "AICC"=(1+npar/ngrad0)/(1-(npar+2)/ngrad0),
                  log(ngrad0)*npar/ngrad0)
  cat("End generating auxiliary objects",format(Sys.time()),"\n")
#
#  avoid situations where si's are larger than s0
#
  grad <- t(object@gradient[,-s0ind])
#
#   determine initial estimates for orientations 
#
  cat("Start search for initial directions at",format(Sys.time()),"\n")
  data("polyeders")
  polyeder <- icosa3
  vert <- polyeder$vertices
# remove redundant directions
  vind <- rep(TRUE,dim(vert)[2])
  vind[vert[1,]<0] <- FALSE
  vind[vert[1,]==0 & vert[2,] <0] <- FALSE
  vind[vert[1,]==0 & vert[2,] == 0 &vert[3,]<0] <- FALSE
  vert <- vert[,vind]
#
#  compute initial estimates (EV from grid and orientations from icosa3$vertices)
#
  siind <- getsiind3(siq,mask,sigma2,grad,t(vert),th,indth,ev,fa,andir,maxcomp,maxc=maxc,nguess=nguess)
  krit <- siind$krit # sqrt(sum of squared residuals) for initial estimates
  siind <- siind$siind # components 1: model order 2: 
                       # grid index for EV 2+(1:m) index of orientations
  cat("Model orders for initial estimates")
  print(table(siind[1,,,]))
  cat("End search for initial values at",format(Sys.time()),"\n")
  order <- array(0,ddim)
#  logarithmic eigen values
  mix <- array(0,c(maxcomp,ddim))
  orient <- array(0,c(2,maxcomp,ddim))
  lev <- array(0,c(2,ddim))
  n1 <- ddim[1]
  n2 <- ddim[2]
  n3 <- ddim[3]
  igc <- 0
  ingc <- 0
  prt0 <- Sys.time()
#
#   loop over voxel in volume
#
  for(i1 in 1:n1) for(i2 in 1:n2) for(i3 in 1:n3){ # begin loop
     if(mask[i1,i2,i3]){ # begin mask
#   only analyze voxel within mask
     mc0 <- maxcomp
     ord <- mc0+1
     for(j in 1:mc0) { 
          iv <- siind[j+2,i1,i2,i3]
          if(iv==0) iv <- j # this should never happen
          orient[,j,i1,i2,i3] <- paroforient(vert[,iv])
     }
#
#   these are the gradient vectors corresponding to minima in spherical coordinates
#
     if(method=="mixtensor"||method=="mixtensoriso"){
     par <- numeric(2*mc0+1)
#  initialize EV-parameter
     if(siind[2,i1,i2,i3]>0){
     par[1] <- th[siind[2,i1,i2,i3]]
     } else {
     par[1] <- .001
     }
#   initialize orientations
     par[rep(2*(1:mc0),rep(2,mc0))+c(0,1)] <- orient[,1:mc0,i1,i2,i3]
     } 
     sigmai <- sigma2[i1,i2,i3]
     krit <- log(sigmai)+penIC[1]
#
#  use AIC/ngrad0, BIC/ngrad0 or AICC/ngrad0 respectively
#
     for(k in mc0:1){ # begin model order
        if(k<ord) {
#
#  otherwise we would reanalyze a model
#
        if(method=="mixtensor"){
           lpar <- 2*k+1
#
           if(optmethod=="BFGS"){
                 z <- optim(par[1:(2*k+1)],mfunpl0,gmfunpl0,siq=siq[i1,i2,i3,],grad=grad,pen=pen,
                         method="BFGS",control=list(maxit=maxit,reltol=reltol))
           } else {
              z <- optim(par[1:(2*k+1)],mfunpl0,siq=siq[i1,i2,i3,],grad=grad,pen=pen,
                         method=optmethod,control=list(maxit=maxit,reltol=reltol))
           }
        } else if (method=="mixtensoriso"){
           lpar <- 2*k+1
#
           if(optmethod=="BFGS"){
                 z <- optim(par[1:(2*k+1)],mfunpl1,gmfunpl1,siq=siq[i1,i2,i3,],grad=grad,pen=pen,
                         method="BFGS",control=list(maxit=maxit,reltol=reltol))
           } else {
              z <- optim(par[1:(2*k+1)],mfunpl1,siq=siq[i1,i2,i3,],grad=grad,pen=pen,
                         method=optmethod,control=list(maxit=maxit,reltol=reltol))
           }
        }        
#
#   estimate of sigma from the best fitting model
#
        if(method=="mixtensor"){
            zz <- mfunplwghts0(z$par[1:lpar],siq[i1,i2,i3,],grad,pen)
        } else if (method=="mixtensoriso"){
            zz <- mfunplwghts1(z$par[1:lpar],siq[i1,i2,i3,],grad,pen)
        }
        ord <- zz$ord
#  replace sigmai by best variance estimate from currently best model
        value <- zz$value 
# thats sum of squared residuals for the restricted model (w>0)
        if(any(zz$lev<0)||ord<k){
           ttt <- krit
#   parameters not interpretable reduce order
        } else {
           si2new <- value/(ngrad0-3*ord-1)
           ttt <- log(si2new)+penIC[1+ord]
           par <- zz$par
        }
#
#     use directions corresponding to largest weights as initial directions
#
        if(ttt < krit) {
           krit <- ttt
           order[i1,i2,i3] <- ord
           lev[,i1,i2,i3] <- zz$lev
           mix[,i1,i2,i3] <- if(ord==maxcomp) zz$mix else c(zz$mix,rep(0,maxcomp-ord))
           orient[,1:ord,i1,i2,i3] <- zz$orient
           sigma2[i1,i2,i3] <- si2new
       }
     }
   } # end model order
    if(igc<ngc){
       igc <- igc+1
    } else {
       igc <- 1
       ingc <- ingc+1
       prt1 <- Sys.time()
       gc()
       cat("Nr. of voxel",ingc*ngc,"time elapsed:",format(difftime(prt1,prta),digits=3),"remaining time:",
            format(difftime(prt1,prt0)/(ingc*ngc)*(sum(mask)-ingc*ngc),digits=3),"\n")
    }
  }# end mask
  }# end loop
  invisible(new("dwiMixtensor",
                model = "homogeneous_prolate",
                call   = args,
                ev     = lev,
                mix    = mix,
                orient = orient,
                order  = order,
                p      = p,
                th0    = s0,
                sigma  = sigma2,
                scorr  = array(1,c(1,1,1)), 
                bw     = c(0,0,0), 
                mask   = mask,
                hmax   = 1,
                gradient = object@gradient,
                btb    = object@btb,
                ngrad  = object@ngrad, # = dim(btb)[2]
                s0ind  = object@s0ind,
                replind = object@replind,
                ddim   = object@ddim,
                ddim0  = object@ddim0,
                xind   = object@xind,
                yind   = object@yind,
                zind   = object@zind,
                voxelext = object@voxelext,
                level = object@level,
                orientation = object@orientation,
                rotation = object@rotation,
                source = object@source,
                outlier = index,
                scale = 1,
                method = method)
            )
   }
)
