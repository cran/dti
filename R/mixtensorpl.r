#
#  Model without isotropic part and explicit weights
#
mfunpl <- function(par,siq,grad){
#
#   evaluate rss for Mixtensor-model (without isotropic component)
#
lpar <- length(par)
m <- (lpar-1)/3
lpar <- 2*m+1
ngrad <- dim(grad)[1]
erg <- .Fortran("mfunpl",as.double(par[1:lpar]),#par(lpar)
                as.double(par[-(1:lpar)]),#w(m)
                as.double(siq),#siq(ngrad)
                as.double(t(grad)),#grad(3,ngrad)
                as.integer(m),#number of components
                as.integer(lpar),#number of parameters (without weights)
                as.integer(ngrad),#number of gradients
                double(ngrad*m),#z(ngrad,m) working array
                erg = double(1),#residual sum of squares
                PACKAGE="dti")$erg
erg
}
mfunpli <- function(par,siq,grad){
#
#   evaluate rss for Mixtensor-model (with isotropic component)
#
lpar <- length(par)
m <- (lpar-2)/3
lpar <- 2*m+1
ngrad <- dim(grad)[1]
erg <- .Fortran("mfunpli",as.double(par[1:lpar]),#par(lpar)
                as.double(par[-(1:lpar)]),#w(m)
                as.double(siq),#siq(ngrad)
                as.double(t(grad)),#grad(3,ngrad)
                as.integer(m),#number of components
                as.integer(lpar),#number of parameters
                as.integer(ngrad),#number of gradients
                double(ngrad*m),#z(ngrad,m) working array
                erg = double(1),#residual sum of squares
                PACKAGE="dti")$erg
erg
}
gmfunpl <- function(par,siq,grad){
#
#   evaluate rss for Mixtensor-model (without isotropic component)
#
lparw <- length(par)
m <- (lparw-1)/3
lpar <- 2*m+1
ngrad <- dim(grad)[1]
erg <- .Fortran("gmfunpl",as.double(par[1:lpar]),#par(lpar)
                as.double(par[-(1:lpar)]),#w(m)
                as.double(siq),#siq(ngrad)
                as.double(t(grad)),#grad(3,ngrad)
                as.integer(m),#number of components
                as.integer(lpar),#number of parameters (without weights)
                as.integer(ngrad),#number of gradients
                double(ngrad*m),#z(ngrad,m) working array
                double(ngrad),#res
                double(ngrad),#resd
                double(ngrad*m),#dkgj
                double(ngrad*m),#dkgj2
                double(ngrad*m),#ddkdphig
                double(ngrad*m),#ddkdetag
                double(ngrad*m*3),#dzdpars
                double(ngrad*m),#work1
                double(ngrad*m),#work2
                dfdparw = double(lparw),#gradient vector
                PACKAGE="dti")$dfdparw
erg
}
gmfunpli <- function(par,siq,grad){
#
#   evaluate rss for Mixtensor-model (without isotropic component)
#
lparw <- length(par)
m <- (lparw-2)/3
lpar <- 2*m+1
ngrad <- dim(grad)[1]
erg <- .Fortran("gmfunpli",as.double(par[1:lpar]),#par(lpar)
                as.double(par[-(1:lpar)]),#w(m)
                as.double(siq),#siq(ngrad)
                as.double(t(grad)),#grad(3,ngrad)
                as.integer(m),#number of components
                as.integer(lpar),#number of parameters (without weights)
                as.integer(ngrad),#number of gradients
                double(ngrad*m),#z(ngrad,m) working array
                double(ngrad),#res
                double(ngrad),#resd
                double(ngrad*m),#dkgj
                double(ngrad*m),#dkgj2
                double(ngrad*m),#ddkdphig
                double(ngrad*m),#ddkdetag
                double(ngrad*m*3),#dzdpars
                double(ngrad*m),#work1
                double(ngrad*m),#work2
                dfdparw = double(lparw),#gradient vector
                PACKAGE="dti")$dfdparw
erg
}

mfunwghts <- function(par,siq,grad){
lparw <- length(par)
m <- (lparw-1)/3
lpar <- 2*m+1
ngrad <- dim(grad)[1]
w <- par[-(1:lpar)]           
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

or <- matrix(par[2:lpar],2,m)[,o,drop=FALSE]

while(length(or[1,or[1,]<0])!=0 && or[1,or[1,]<0]){
   or[1,or[1,]<0] <- or[1,or[1,]<0]+pi
}


while(length(or[1,or[1,]>pi])!=0 && or[1,or[1,]>pi]){
   or[1,or[1,]>pi] <- or[1,or[1,]>pi]-pi
}


while(length(or[2,or[2,]<0])!=0 && or[2,or[2,]<0]){
   or[2,or[2,]<0] <- or[2,or[2,]<0]+2*pi
}


while(length(or[2,or[2,]>2*pi])!=0 && or[2,or[2,]>2*pi]){
   or[2,or[2,]>2*pi] <- or[2,or[2,]>2*pi]-2*pi
}

par <- c(par[1],or[,1:ord])
list(ord=ord,lev=lev,mix=mix,orient=or,par=par,w=w[o])
}

mfunwghtsi <- function(par,siq,grad){
lparw <- length(par)
m <- (lparw-2)/3
if(m>0){
   lpar <- 2*m+1
   ngrad <- dim(grad)[1]
   w <- par[-(1:lpar)]           
   w0 <- w[1]
   w <- w[-1]
   o <- order(w,decreasing=TRUE)
   ord <- sum(w>0)
   if(ord<m){
      o <- o[1:ord]
   }
   sw <- sum(w[w>0])+max(w0,0)
   lev <- c(par[1],-log(sw))
   if(ord>0){
      mix <- w[o]/sw
   } else {
      mix <- NULL
   } 
   
   or <- matrix(par[2:lpar],2,m)[,o,drop=FALSE]

   while(length(or[1,or[1,]<0])!=0 && or[1,or[1,]<0]){
      or[1,or[1,]<0] <- or[1,or[1,]<0]+pi
   }
   
   
   while(length(or[1,or[1,]>pi])!=0 && or[1,or[1,]>pi]){
      or[1,or[1,]>pi] <- or[1,or[1,]>pi]-pi
   }
   
   
   while(length(or[2,or[2,]<0])!=0 && or[2,or[2,]<0]){
      or[2,or[2,]<0] <- or[2,or[2,]<0]+2*pi
   }
   
   
   while(length(or[2,or[2,]>2*pi])!=0 && or[2,or[2,]>2*pi]){
      or[2,or[2,]>2*pi] <- or[2,or[2,]>2*pi]-2*pi
   }

   par <- c(par[1],or[,1:ord])
   w <- w[o]
} else {
   ord <- 0
   lev <- c(par[1],0)
   or <- NULL
   w0 <- 1
   w <- NULL
   par <- par[1]
   mix <- NULL
}
list(ord=ord,lev=lev,mix=mix,orient=or,par=par,w=c(w0,w))
}

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
erg <- .Fortran("mfunpl0",as.double(par),#par(lpar)
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
erg
}


mfunpl0h <- function(par,siq,grad){
#
#   evaluate rss for Mixtensor-model (without isotropic component)
#   uses LawsonHanson-nnls code
#
lpar <- length(par)
m <- (lpar-1)/2
ngrad <- dim(grad)[1]
erg <- .Fortran("mfunpl0h",as.double(par),#par(lpar)
                as.double(siq),#siq(ngrad)
                as.double(t(grad)),#grad(3,ngrad)
                as.integer(m),#number of components
                as.integer(lpar),#number of parameters
                as.integer(ngrad),#number of gradient
                double(ngrad*m),#z(ngrad,m) working array
                double(ngrad),#w(ngrad) working array
                double(ngrad),# b(ngrad) working array
                double(ngrad),# work1(ngrad) working array                
                erg = double(1),#residual sum of squares
                PACKAGE="dti")$erg
erg
}

gmfunpl0 <- function(par,siq,grad,pen=1e2){
#
#   evaluate rss for Mixtensor-model
#
lpar <- length(par)
m <- (lpar-1)/2
ngrad <- dim(grad)[1]
dfdpar<-.Fortran("mfunpl0g",
         as.double(par),#par(lpar)
         as.double(siq),#s(n)
         as.double(t(grad)),#g(3,n)
         as.integer(m),
         as.integer(lpar),
         as.integer(ngrad),
         double(ngrad*m),# z(n,m)
         double(m*m),# v(m,m)
         double(ngrad),# w(n) need w((m+1):n) for solver dgelsy
         double(ngrad*m),# dkgj(n,m)
         double(ngrad*m),# dkgj2(n,m)
         double(ngrad*m),# ddkdphig(n,m)
         double(ngrad*m),# ddkdetag(n,m)
         double(m*m),# dvdth(m,m)
         double(m*m*m),# dvdphi(m,m,m)
         double(m*m*m),# dvdeta(m,m,m)
         double(ngrad*m*3),# dzdpars(n,m,3)
         double(m*lpar),# dwdpars(m,lpar)
         double(m*lpar),# dwdpars2(m,lpar)
         double(ngrad*m),# zs(n,m)
         double(ngrad*m),# work1(n,m)
         double(ngrad*m),# work2(n,m)
         double(ngrad),# scopy(n)
         as.double(pen),# pen
         dfdpar=double(lpar),# dfdpar(lpar)
         PACKAGE="dti")$dfdpar
dfdpar
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
                w <- pmax(0,w)

                if(all(is.finite(w))) {
erg<-.Fortran("mfunpl0w",as.double(par),#par(lpar)
                as.double(w),
                as.double(siq),#siq(ngrad)
                as.double(t(grad)),#grad(3,ngrad)
                as.integer(m),#number of components
                as.integer(lpar),#number of parameters
                as.integer(ngrad),#number of gradients
                double(ngrad*m),#z(ngrad,m) working array
                erg = double(1),#residual sum of squares
                PACKAGE="dti")$erg
           } else {
             cat("got w=",w,"\n")
             erg <- 1e20
           }

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

           or <- matrix(par[2:lpar],2,m)[,o,drop=FALSE]

      while(length(or[1,or[1,]<0])!=0 && or[1,or[1,]<0]){
         or[1,or[1,]<0] <- or[1,or[1,]<0]+pi
      }


      while(length(or[1,or[1,]>pi])!=0 && or[1,or[1,]>pi]){
         or[1,or[1,]>pi] <- or[1,or[1,]>pi]-pi
      }


      while(length(or[2,or[2,]<0])!=0 && or[2,or[2,]<0]){
         or[2,or[2,]<0] <- or[2,or[2,]<0]+2*pi
      }


      while(length(or[2,or[2,]>2*pi])!=0 && or[2,or[2,]>2*pi]){
         or[2,or[2,]>2*pi] <- or[2,or[2,]>2*pi]-2*pi
      }

           par <- c(par[1],or[,1:ord])

  list(ord=ord,lev=lev,mix=mix,orient=or,par=par,value=erg)
}

mfunplwghts0h <- function(par,siq,grad){
#
#   get weights for Mixtensor-model (without isotropic component) and extract parameters 
#
par[1] <- max(0,par[1])
#   uses LawsonHanson-nnls code
#
lpar <- length(par)
m <- (lpar-1)/2
ngrad <- dim(grad)[1]
z <- .Fortran("mfunpl0h",as.double(par),#par(lpar)
                as.double(siq),#siq(ngrad)
                as.double(t(grad)),#grad(3,ngrad)
                as.integer(m),#number of components
                as.integer(lpar),#number of parameters
                as.integer(ngrad),#number of gradients
                double(ngrad*m),#z(ngrad,m) working array
                w=double(ngrad),#w(ngrad) working array
                double(ngrad),# b(ngrad) working array
                double(ngrad),# work1(ngrad) working array                
                erg = double(1),#residual sum of squares
                PACKAGE="dti")[c("erg","w")]
           erg <- z$erg^2
           w <- z$w[1:m]

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

           or <- matrix(par[2:lpar],2,m)[,o,drop=FALSE]

           while(length(or[1,or[1,]<0])!=0 && or[1,or[1,]<0]){
              or[1,or[1,]<0] <- or[1,or[1,]<0]+pi
              }


      while(length(or[1,or[1,]>pi])!=0 && or[1,or[1,]>pi]){
         or[1,or[1,]>pi] <- or[1,or[1,]>pi]-pi
         }


      while(length(or[2,or[2,]<0])!=0 && or[2,or[2,]<0]){
         or[2,or[2,]<0] <- or[2,or[2,]<0]+2*pi
         }


      while(length(or[2,or[2,]>2*pi])!=0 && or[2,or[2,]>2*pi]){
         or[2,or[2,]>2*pi] <- or[2,or[2,]>2*pi]-2*pi
         }

           par <- c(par[1],or[,1:ord])

  list(ord=ord,lev=lev,mix=mix,orient=or,par=par,value=erg)
}
#
#  Model with isotropic part
#
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

getsiind3 <- function(si,mask,sigma2,grad,vico,th,indth,ev,fa,andir,maxcomp=3,
maxc=.866,nguess=100,mc.cores = getOption("mc.cores", 2L)){
# assumes dim(grad) == c(ngrad,3)
# assumes dim(si) == c(ngrad,n1,n2,n3)
# SO removed
ngrad <- dim(grad)[1]
nvico <- dim(vico)[1]
ddim <- dim(fa)
nsi <- dim(si)[1]
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
iandir <- .Fortran("iandir",
                   as.double(t(vico)),
                   as.integer(nvico),
                   as.double(andir),
                   as.integer(nvoxel),
                   as.logical(landir),
                   iandir=integer(prod(ddim)),
                   DUPL=FALSE,
                   PACKAGE="dti")$iandir
isample0 <- selisample(nvico,maxcomp,nguess,dgradi,maxc)
if(maxcomp>1) isample1 <- selisample(nvico,maxcomp-1,nguess,dgradi,maxc)
if(maxcomp==1) isample1 <- sample(ngrad, nguess, replace = TRUE)
#
#  eliminate configurations with close directions 
#
# this provides configurations of initial estimates with minimum angle between 
# directions > acos(maxc)
nvoxel <- prod(dim(si)[-1])
cat("using ",nguess,"guesses for initial estimates\n")
siind <- matrix(as.integer(0),maxcomp+2,nvoxel)
krit <- numeric(nvoxel)
# first voxel with fa<.3
if(sum(mask&!landir)>0){
cat(sum(mask&!landir),"voxel with small FA\n")
nguess <- length(isample0)/maxcomp
if(mc.cores<=1){
z <- .Fortran("getsii30",
         as.double(si),
         as.double(sigma2),
         as.integer(nsi),
         as.integer(nvoxel),
         as.integer(maxcomp),
         as.double(dgrad),
         as.integer(nvico),
         as.double(th),
         as.integer(nth),
         as.integer(indth),
         double(nsi*nvico),
         as.integer(isample0),
         as.integer(nguess),
         double(nsi),
         double(nsi*(maxcomp+2)),
         siind=integer((maxcomp+2)*nvoxel),
         krit=double(nvoxel),
         as.integer(maxcomp+2),
         as.logical(mask&!landir),
         DUP=FALSE,
         PACKAGE="dti")[c("siind","krit")]
dim(z$siind) <- c(maxcomp+2,nvoxel)
siind[,!landir] <- z$siind[,!landir]
krit[!landir] <- z$krit[!landir]
} else {
   x <- matrix(0,nsi+2,sum(mask&!landir))
   x[1:nsi,] <- matrix(si,nsi,nvoxel)[,mask&!landir]
   x[nsi+1,] <- sigma2[mask&!landir]
   x[nsi+2,] <- indth[mask&!landir]
   z <- plmatrix(x,pgetsii30,maxcomp=maxcomp,dgrad=dgrad,th=th,isample0=isample0,
                nsi=nsi,nth=length(th),nvico=nvico,nguess=nguess,
                mc.cores=mc.cores)
   dim(z) <- c(maxcomp+3,sum(mask&!landir))
   siind[,mask&!landir] <- z[-1,]
   krit[mask&!landir] <- z[1,]
}
}
# now voxel where first tensor direction seems important
if(sum(mask&landir)>0){
if(maxcomp >0){
cat(sum(mask&landir),"voxel with distinct first eigenvalue \n")
nguess <- if(maxcomp>1) length(isample1)/(maxcomp-1) else length(isample1)
if(mc.cores<=1){
z <- .Fortran("getsii31",
         as.double(si),
         as.double(sigma2),
         as.integer(nsi),
         as.integer(nvoxel),
         as.integer(maxcomp),
         as.double(dgrad),
         as.integer(nvico),
         as.integer(iandir),
         as.double(th),
         as.integer(nth),
         as.integer(indth),
         double(nsi*nvico),
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
         DUP=FALSE,
         PACKAGE="dti")[c("siind","krit")]
dim(z$siind) <- c(maxcomp+2,nvoxel)
siind[,landir] <- z$siind[,landir]
krit[landir] <- z$krit[landir]
} else {
   x <- matrix(0,nsi+3,sum(mask&landir))
   x[1:nsi,] <- matrix(si,nsi,nvoxel)[,mask&landir]
   x[nsi+1,] <- sigma2[mask&landir]
   x[nsi+2,] <- indth[mask&landir]
   x[nsi+3,] <- iandir[mask&landir]
   z <- plmatrix(x,pgetsii31,maxcomp=maxcomp,dgrad=dgrad,th=th,isample1=isample1,
                nsi=nsi,nth=length(th),nvico=nvico,nguess=nguess,
                dgradi=dgradi,maxc=maxc,mc.cores=mc.cores)
   dim(z) <- c(maxcomp+3,sum(mask&landir))
   siind[,mask&landir] <- z[-1,]
   krit[mask&landir] <- z[1,]
}
}
}
failed <- (krit^2/ngrad) > (sigma2-1e-10)
if(any(failed[mask])){
}
list(siind=array(siind,c(maxcomp+2,dim(si)[-1])),
     krit=array(krit,dim(si)[-1]))
}
getsiind3iso <- function(si,mask,sigma2,grad,vico,th,indth,ev,fa,andir,maxcomp=3,maxc=.866,nguess=100){
# assumes dim(grad) == c(ngrad,3)
# assumes dim(si) == c(ngrad,n1,n2,n3)
# SO removed
ngrad <- dim(grad)[1]
nvico <- dim(vico)[1]
ddim <- dim(fa)
nsi <- dim(si)[1]
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
iandir <- .Fortran("iandir",
                   as.double(t(vico)),
                   as.integer(nvico),
                   as.double(andir),
                   as.integer(nvoxel),
                   as.logical(landir),
                   iandir=integer(prod(ddim)),
                   DUPL=FALSE,
                   PACKAGE="dti")$iandir
isample0 <- selisample(nvico,maxcomp,nguess,dgradi,maxc)
if(maxcomp>1) isample1 <- selisample(nvico,maxcomp-1,nguess,dgradi,maxc)
if(maxcomp==1) isample1 <- sample(ngrad, nguess, replace = TRUE)
#
#  eliminate configurations with close directions 
#
# this provides configurations of initial estimates with minimum angle between 
# directions > acos(maxc)
nvoxel <- prod(dim(si)[-1])
cat("using ",nguess,"guesses for initial estimates\n")
siind <- matrix(as.integer(0),maxcomp+2,nvoxel)
krit <- numeric(nvoxel)
# first voxel with fa<.3
cat(sum(mask&!landir),"voxel with small FA\n")
nguess <- length(isample0)/maxcomp
z <- .Fortran("getsi30i",
         as.double(si),
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
         double(nsi*(maxcomp+3)),
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
z <- .Fortran("getsi31i",
         as.double(si),
         as.double(sigma2),
         as.integer(nsi),
         as.integer(nvoxel),
         as.integer(maxcomp),
         as.double(dgrad),
         as.integer(nvico),
         as.integer(iandir),
         as.double(th),
         as.integer(nth),
         as.integer(indth),
         double(ngrad*nvico),
         as.integer(isample1),
         as.integer(nguess),
         double(nsi),
         double(nsi*(maxcomp+3)),
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
#print((krit[mask])[failed[mask]])
#print(((1:prod(dim(si)[1:3]))[mask])[failed[mask]])
#print(sum(failed[mask]))
}
list(siind=array(siind,c(maxcomp+2,dim(si)[-1])),
     krit=array(krit,dim(si)[-4]))
}


dwiMixtensor <- function(object, ...) cat("No dwiMixtensor calculation defined for this class:",class(object),"\n")

setGeneric("dwiMixtensor", function(object,  ...) standardGeneric("dwiMixtensor"))

setMethod("dwiMixtensor","dtiData",function(object, maxcomp=3, method="mixtensor", reltol=1e-6, maxit=5000,ngc=1000, optmethod="BFGS", nguess=100*maxcomp^2,msc="BIC",pen=NULL,code="C",thinit=NULL, 
    mc.cores = getOption("mc.cores", 1L)){
#
#  uses  S(g)/s_0 = w_0 exp(-l_1) +\sum_{i} w_i exp(-l_2-(l_1-l_2)(g^T d_i)^2)
#
#  choices for optmethod:
#  BFGS  -  BFGS with analytic gradients and penalization
#  CG - Conjugate gradients with analytic gradients and penalization
#  L-BFGS-B  -  constrained BFGS with analytic gradients 
#  Nelder-Mead - using LawsonHanson-nnls code
#
#  Defaults: 
#     BFGS for tensor mixture models without isotropic compartment
#     L-BFGS-B for tensor mixture models with isotropic compartment
#
  if(method!="mixtensor"||optmethod!="BFGS") {
     cat("Using R-code")
     code <- "R"
  }
  set.seed(1)
  bvalue <- object@bvalue
  bvalue <- bvalue[bvalue>.1*median(bvalue)]
  if(sd(bvalue)>.1*median(bvalue)){
     warning("b-values indicate measurements on multiple shells,
         tensor mixtures not yet implemented\n returning original object")
     return(object)
  }
  if(is.null(pen)) pen <- 100
  if(method=="mixtensoriso") optmethod <- "L-BFGS-B"
  theta <- .5
  maxc <- .866
  args <- sys.call(-1)
  args <- c(object@call,args)
  ngrad <- object@ngrad
  ddim <- object@ddim
  nvox <- prod(ddim)
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
  tensorobj <- dtiTensor(object, mc.cores = mc.cores)
  cat("Start evaluation of eigenstructure at",format(Sys.time()),"\n")
  z <- dtieigen(tensorobj@D, tensorobj@mask, mc.cores = mc.cores)
  rm(tensorobj)
  gc()
  fa <- array(z$fa,ddim)
  ev <- array(z$ev,c(3,ddim))*median(bvalue)
#
#  rescale by bvalue to go to implemented scale
#
  andir <- array(z$andir,c(3,2,ddim))
  rm(z)
  gc()
#  nth <- 11
  if(is.null(thinit)){
  nth <- 5
  th <- ev[1,,,] - (ev[2,,,]+ev[3,,,])/2
  falevel <- min(quantile(fa[fa>0],.75),.4)
  cat("falevel",falevel,"\n")
  qth <- unique(quantile(th[fa>=falevel&fa<.95],seq(.8,.99,length=nth)))
  nth <- length(qth)
  if(nth>1){
     indth <- cut(th,qth,labels=FALSE)
     indth[th<=qth[1]] <- 1
     indth[th>=qth[nth]] <- nth
     th <- qth
  } else {
    indth <- rep(1,nvox)
    th <- qth
  }
  } else {
    nth <- 1
    indth <- rep(1,nvox)
    th <- thinit
  }
cat("using th:::",th,"\n")
  cat("Start search outlier detection at",format(Sys.time()),"\n")
#
#  replace physically meaningless S_i by mena S_0 values
#
  z <- sioutlier(object@si,s0ind,mc.cores=mc.cores)
#
#  this does not scale well with openMP
#
  cat("End search outlier detection at",format(Sys.time()),"\n")
  si <- array(as.integer(z$si),c(ngrad,ddim))
  index <- z$index
  rm(z)
  gc()
  cat("Start generating auxiliary objects",format(Sys.time()),"\n")
#
#  compute mean S_0, s_i/S_0 (siq), var(siq) and mask
#
  nvox <- prod(ddim[1:3])
  cat("sweeps0:")
  t1 <- Sys.time()
  if(mc.cores==1||ngrad0>250){
  z <- .Fortran("sweeps0",# mixtens.f
                as.integer(si[-s0ind,,,,drop=FALSE]),
                as.integer(si[s0ind,,,,drop=FALSE]),
                as.integer(nvox),
                as.integer(ns0),
                as.integer(ngrad0),
                as.integer(object@level),
                siq=double(nvox*ngrad0),
                s0=double(nvox),
                vsi=double(nvox),
                mask=logical(nvox),
                DUPL=FALSE,
                PACKAGE="dti")[c("siq","s0","vsi","mask")]
  t2 <- Sys.time()
  cat(difftime(t2,t1),"for",nvox,"voxel\n")
  s0 <- array(z$s0,ddim[1:3])
  siq <- array(z$siq,c(ngrad0,ddim[1:3]))
#
#  siq is permutated c(4,1:3)
#
  sigma2 <- array(z$vsi,ddim[1:3])
  mask <- array(z$mask,ddim[1:3])
  } else {
  mc.cores.old <- setCores(,reprt=FALSE)
  setCores(mc.cores)
  z <- matrix(.Fortran("sweeps0p",# mixtens.f
                as.integer(si[-s0ind,,,,drop=FALSE]),
                as.integer(si[s0ind,,,,drop=FALSE]),
                as.integer(nvox),
                as.integer(ns0),
                as.integer(ngrad0),
                as.integer(object@level),
                siq=double(nvox*(ngrad0+3)),
                as.integer(ngrad0+3),
                DUPL=FALSE,
                PACKAGE="dti")$siq,ngrad0+3,nvox)
  t2 <- Sys.time()
  cat(difftime(t2,t1),"for",nvox,"voxel\n")
  setCores(mc.cores.old,reprt=FALSE)
  s0 <- array(z[ngrad0+1,],ddim[1:3])
  siq <- array(z[1:ngrad0,],c(ngrad0,ddim[1:3]))
#
#  siq is permutated c(4,1:3)
#
  sigma2 <- array(z[ngrad0+2,],ddim[1:3])
  mask <- array(as.logical(z[ngrad0+3,]),ddim[1:3])
  }
  rm(si)
  rm(z)
  gc()
  npar <- if(method=="mixtensor") 1+3*(0:maxcomp) else c(1,2+3*(1:maxcomp))
#
#   compute penalty for model selection, default BIC
#
  penIC <- switch(msc,"AIC"=2*npar/ngrad0,"BIC"=log(ngrad0)*npar/ngrad0,
                  "AICC"=(1+npar/ngrad0)/(1-(npar+2)/ngrad0),
                  "None"=log(ngrad0)-log(ngrad0-npar),
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
  siind <- if(method=="mixtensor")  getsiind3(siq,mask,sigma2,grad,t(vert),th,indth,ev,fa,andir,maxcomp,maxc=maxc,nguess=nguess,mc.cores=mc.cores) else getsiind3iso(siq,mask,sigma2,grad,t(vert),th,indth,ev,fa,andir,maxcomp,maxc=maxc,nguess=nguess)
  krit <- siind$krit # sqrt(sum of squared residuals) for initial estimates
  siind <- siind$siind # components 1: model order 2: 
                       # grid index for EV 2+(1:m) index of orientations
  cat("Model orders for initial estimates")
  print(table(siind[1,,,]))
  cat("End search for initial values at",format(Sys.time()),"\n")
#  logarithmic eigen values
  orient <- array(0,c(2,maxcomp,ddim))
  n1 <- ddim[1]
  n2 <- ddim[2]
  n3 <- ddim[3]
  igc <- 0
  ingc <- 0
  prt0 <- Sys.time()
#
#   loop over voxel in volume
#
  if(code=="C"){  
#
#     C-Code
#
  if(method=="mixtensor") meth = 1 else meth = 2
  optmeth <- switch(optmethod, "BFGS" = 1,
                    "CG" = 2, "Nelder-Mead" = 3, "L-BFGS-B" = 4)

if(mc.cores<=1){
  cat("Starting parameter estimation and model selection (C-code)",format(Sys.time()),"\n")
  dim(siq) <- c(ngrad0,nvox)
  dim(siind) <- c(2+maxcomp,nvox)
  nvoxm <- sum(mask)
  z <- .C("mixture2", 
          as.integer(meth),
          as.integer(optmeth), 
          as.integer(1), 
          as.integer(1), 
          as.integer(nvoxm),
          as.integer(rep(1L,nvoxm)), 
          as.integer(siind[,mask]), 
          as.integer(ngrad0),
          as.integer(maxcomp),
          as.integer(maxit),
          as.double(pen),
          as.double(t(grad)),
          as.double(reltol),
          as.double(th),
          as.double(penIC),
          as.double(sigma2[mask]),
          as.double(vert),
          as.double(t(siq[,mask])),
          sigma2  = double(nvoxm),# error variance 
          orient  = double(2*maxcomp*nvoxm), # phi/theta for all mixture tensors
          order   = integer(nvoxm),   # selected order of mixture
          lev     = double(2*nvoxm),         # logarithmic eigenvalues
          mix     = double(maxcomp*nvoxm),   # mixture weights
          DUPL=FALSE, PACKAGE="dti")[c("sigma2","orient","order","lev","mix")]
  cat("End parameter estimation and model selection (C-code)",format(Sys.time()),"\n")
  sigma2 <-  array(0,ddim)
  sigma2[mask] <- z$sigma2
  orient <- matrix(0,2*maxcomp,nvox)
  orient[,mask] <- z$orient
  dim(orient) <- c(2, maxcomp, ddim)
  order <- array(0, ddim)
  order[mask] <- z$order
  lev <- matrix(0,2,nvox)
  lev[,mask] <- z$lev
  dim(lev) <- c(2,ddim)
  mix <- matrix(0,maxcomp,nvox)
  mix[,mask] <- z$mix
  dim(mix) <- c(maxcomp, ddim)
} else {
  cat("Starting parameter estimation and model selection (C-code) on",mc.cores," cores",format(Sys.time()),"\n")
  x <- matrix(0,ngrad0+3+maxcomp,sum(mask))
  dim(siq) <- c(ngrad0,nvox)
  x[1:ngrad0,] <- siq[,mask]
  x[ngrad0+1,] <- sigma2[mask]
  dim(siind) <- c(2+maxcomp,nvox)
  x[ngrad0+2:(3+maxcomp),] <- siind[,mask] 
  res <- matrix(0,4+3*maxcomp,nvox)
  res[,mask] <- plmatrix(x,pmixtens,
                      meth=meth,optmeth=optmeth,
                      ngrad0=ngrad0,maxcomp=maxcomp,maxit=maxit,
                      pen=pen,grad=grad,reltol=reltol,th=th,
                      penIC=penIC,vert=vert,
                      mc.cores=mc.cores)
  cat("End parameter estimation and model selection (C-code)",format(Sys.time()),"\n")
  rm(x)
  gc()
  sigma2 <-  array(res[2,],ddim)
  orient <- array(res[maxcomp+4+1:(2*maxcomp),], c(2, maxcomp, ddim))
  order <- array(as.integer(res[1,]), ddim)
  lev <- array(res[3:4,], c(2,ddim))
  mix <- array(res[4+(1:maxcomp),], c(maxcomp, ddim))
  }
  method <- "mixtensor"
  } else {
  order <- array(0,ddim)
  mix <- array(0,c(maxcomp,ddim))
  lev <- array(0,c(2,ddim))
  cat("Starting parameter estimation and model selection (R-code)",format(Sys.time()),"\n")
#
#     R/Fortran-Code 
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
     if(optmethod != "L-BFGS-B"){
        param <- numeric(2*mc0+1) 
     } else {
        if(method == "mixtensor"){
           param <- numeric(3*mc0+1) # additional m weights
        } else {
           param <- numeric(3*mc0+2) # additional isotropic weight
        }     
     }
#  initialize EV-parameter
     if(siind[2,i1,i2,i3]>0){
       param[1] <- th[siind[2,i1,i2,i3]]
     } else {
       param[1] <- .001
     }
#   initialize orientations
     param[2:(2*mc0+1)] <- orient[,1:mc0,i1,i2,i3]
     sigmai <- sigma2[i1,i2,i3]
     krit <- log(sigmai)+penIC[1]

#      cat("param from orient:\n", param, "\n")
#
#  use AIC/ngrad0, BIC/ngrad0 or AICC/ngrad0 respectively
#
     for(k in mc0:1){ # begin model order
        if(k<ord) {
#
#  otherwise we would reanalyze a model
#
           if(optmethod=="L-BFGS-B"){
      param <- if(method=="mixtensor") param[1:(3*k+1)] else param[1:(3*k+2)]
      if(k==mc0){
         param[-(1:(2*k+1))] <- if(method=="mixtensor") rep(1/k,k) else rep(1/k,k+1)
      } else {
         param[-(1:(2*k+1))] <- zz$w[1:(if(method=="mixtensor") k else (k+1))]
      }
          }

          if(method=="mixtensor"){

            z <- switch(optmethod,
                    "BFGS"=optim(param[1:(2*k+1)],mfunpl0,gmfunpl0,
                           siq=siq[,i1,i2,i3],grad=grad,pen=pen,
                           method="BFGS",control=list(maxit=maxit,reltol=reltol)),
                    "CG"=optim(param[1:(2*k+1)],mfunpl0,gmfunpl0,
                           siq=siq[,i1,i2,i3],grad=grad,pen=pen,
                           method="CG",control=list(maxit=maxit,reltol=reltol)),
                    "Nelder-Mead"=optim(param[1:(2*k+1)],mfunpl0h,
                           siq=siq[,i1,i2,i3],grad=grad,method="Nelder-Mead",
                           control=list(maxit=maxit,reltol=reltol)),
                    "L-BFGS-B"=optim(param[1:(3*k+1)],mfunpl,gmfunpl,
                           siq=siq[,i1,i2,i3],grad=grad,method="L-BFGS-B",
                           lower=c(0,rep(-Inf,2*k),rep(0,k)),
                           control=list(maxit=maxit,factr=reltol/.Machine$double.eps)))
          } else { # method=="mixtensoriso"
            z <- optim(param[1:(3*k+2)],fn=mfunpli,gr=gmfunpli,
                           siq=siq[,i1,i2,i3],grad=grad,method="L-BFGS-B",
                           lower=c(0,rep(-Inf,2*k),rep(0,k+1)),
                           control=list(maxit=maxit,factr=reltol/.Machine$double.eps))
#
#  other methods seem less numerically stable in this situation
#
          }        
#
#   estimate of sigma from the best fitting model
#
          if(method=="mixtensor"){
            zz <- switch(optmethod,
                    "BFGS"=mfunplwghts0(z$par[1:(2*k+1)],siq[,i1,i2,i3],grad,pen),
                    "CG"=mfunplwghts0(z$par[1:(2*k+1)],siq[,i1,i2,i3],grad,pen),
                    "Nelder-Mead"=mfunplwghts0h(z$par[1:(2*k+1)],siq[,i1,i2,i3],grad),
                    "L-BFGS-B"=mfunwghts(z$par[1:(3*k+1)],siq[,i1,i2,i3],grad))
          } else {
            zz <- mfunwghtsi(z$par[1:(3*k+2)],siq[,i1,i2,i3],grad)
          }
          value <- if(optmethod=="L-BFGS-B") z$value else zz$value
          ord <- zz$ord

#  replace sigmai by best variance estimate from currently best model
# thats sum of squared residuals for the restricted model (w>0)
          if(any(zz$lev<0)||ord<k){
            ttt <- krit
          } else {
#           si2new <- value/(ngrad0-3*ord-1) 
            si2new <- value/ngrad0 # changed 2011/05/13 
            ttt <- log(si2new)+penIC[1+ord]
            param[1:(2*k+1)] <- zz$par
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
  cat("End parameter estimation and model selection (R-code)",format(Sys.time()),"\n")
  } # end Code-verzweigung 
  invisible(new("dwiMixtensor",
                model = "homogeneous_prolate",
                call   = args,
                ev     = lev,
                mix    = mix,
                orient = orient,
                order  = order,
                p      = 0,
                th0    = s0,
                sigma  = sigma2,
                scorr  = array(1,c(1,1,1)), 
                bw     = c(0,0,0), 
                mask   = mask,
                hmax   = 1,
                gradient = object@gradient,
                bvalue = object@bvalue,
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
