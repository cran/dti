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
           o <- order(w,decreasing=TRUE)
           ord <- sum(w>0)
           if(ord<m-1){
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
           or[1,or[1,]<0] <- or[1,or[1,]<0]+pi
           or[1,or[1,]>pi] <- or[1,or[1,]>pi]-pi
           or[2,or[2,]<0] <- or[2,or[2,]<0]+2*pi
           or[2,or[2,]>2*pi] <- or[2,or[2,]>2*pi]-2*pi
           par <- c(par[1],or[,1:ord])
list(ord=ord,lev=lev,mix=mix,orient=or,par=par)
}
selisample <- function(ngrad,maxcomp,nguess,dgrad,maxc){
saved.seed <- .Random.seed
set.seed(1)
isample <- matrix(sample(ngrad,maxcomp*nguess,replace=TRUE),maxcomp,nguess)
ind <- rep(TRUE,nguess)
for(i in 1:nguess) for(j in 1:(maxcomp-1)) for(k in (j+1):maxcomp){
if(dgrad[isample[j,i],isample[k,i]]>maxc) ind[i] <- FALSE
}
.Random.seed <- saved.seed
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

getsiind2 <- function(si,mask,grad,th,maxcomp=3,maxc=.866,nguess=100){
# assumes dim(grad) == c(ngrad,3)
# assumes dim(si) == c(n1,n2,n3,ngrad)
# SO removed
ngrad <- dim(grad)[1]
nsi <- dim(si)[4]
dgrad <- matrix(abs(grad%*%t(grad)),ngrad,ngrad)
dgrad <- dgrad/max(dgrad)
isample <- selisample(ngrad,maxcomp,nguess,dgrad,maxc)
nth <- length(th)
#
#  eliminate configurations with close directions 
#
nguess <- dim(isample)[2]
# this provides configurations of initial estimates with minimum angle between 
# directions > acos(maxc)
cat("using ",nguess,"guesses for initial estimates\n")
siind <- .Fortran("getsiin2",
         as.double(aperm(si,c(4,1:3))),
         as.integer(nsi),
         as.integer(dim(si)[1]),
         as.integer(dim(si)[2]),
         as.integer(dim(si)[3]),
         as.integer(maxcomp),
         as.double(dgrad),
         as.double(th),
         as.integer(nth),
         double(ngrad*ngrad),
         as.integer(isample),
         as.integer(nguess),
         double(nsi),
         double(nsi*(maxcomp+1)),
         siind=integer((maxcomp+2)*prod(dim(si)[-4])),
         krit=double(prod(dim(si)[1:3])),
         as.integer(maxcomp+2),
         as.logical(mask),
         PACKAGE="dti")[c("siind","krit")]
array(siind$siind,c(maxcomp+2,dim(si)[-4]))
}

dwiMixtensor <- function(object, ...) cat("No dwiMixtensor calculation defined for this class:",class(object),"\n")

setGeneric("dwiMixtensor", function(object,  ...) standardGeneric("dwiMixtensor"))

setMethod("dwiMixtensor","dtiData",function(object, maxcomp=3,  p=40, method="mixtensor", reltol=1e-6, maxit=5000,ngc=1000, optmethod="BFGS", nguess=50*maxcomp^2,penalty="BIC"){
#
#  uses  S(g)/s_0 = w_0 exp(-l_1) +\sum_{i} w_i exp(-l_2-(l_1-l_2)(g^T d_i)^2)
#
  pen <- 1e2
  theta <- .5
  maxc <- .866
  args <- sys.call(-1)
  args <- c(object@call,args)
  prta <- proc.time()[1]
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  cat("Start search outlier detection at",date(),"\n")
  z <- .Fortran("outlier",
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
  cat("End search outlier detection at",date(),"\n")
  si <- array(z$si,c(ddim,ngrad))
  index <- if(z$lindex>0) z$index[1:z$lindex] else numeric(0)
  rm(z)
  ngrad0 <- ngrad - ns0
  s0 <- si[,,,s0ind,drop=FALSE]
  if(length(s0ind)>1) s0 <- apply(s0,1:3,mean) else dim(s0) <- dim(s0)[1:3]
  mask <- s0 > object@level
  siq <- si[,,,-s0ind,drop=FALSE]
  dim(siq) <- c(prod(ddim),ngrad0)
  siq[mask,] <- sweep(siq[mask,],1,s0[mask],"/")
  dim(siq) <- c(ddim,ngrad0)
  siqmed <- apply(siq,1:3,median)
  siqmed[siqmed<.9] <- .9
  siqmed[siqmed>.99] <- .99
  siq <- sweep(siq,1:3,siqmed,pmin)
  if(penalty=="BIC") penIC <- log(ngrad0) else penIC <- 2
#
#  avoid situations where si's are larger than s0
#
  grad <- t(object@gradient[,-s0ind])
#
# initial estimates for eigenvalues
#
  cat("Start search for initial estimates of eigenvalues at",date(),"\n")
  lev <- array(.Fortran("getev0",
               as.double(aperm(siq,c(4,1:3))),
               as.integer(ngrad0),
               as.integer(ddim[1]),
               as.integer(ddim[2]),
               as.integer(ddim[3]),
               lev=double(2*prod(ddim)),
               DUPL=FALSE,
               PACKAGE="dti")$lev,c(2,ddim))
  mlev <- median(lev[2,,,][mask])
  theta <- theta*mlev
  minth <- .2*theta
  maxth <- 1.4*theta
  nth <- 7
  th <- seq(minth,maxth,length=nth)
  cat("using theta=",th,"\n")
#
#   determine initial estimates for orientations 
#
  cat("Start search for initial directions at",date(),"\n")
  data("polyeders")
#  siind <- getsiind2(siq,mask,grad,theta,maxcomp,maxc=maxc,nguess=nguess)
  siind <- getsiind2(siq,mask,t(icosa3$vertices),th,maxcomp,maxc=maxc,nguess=nguess)
  print(table(siind[2,,,]))
 cat("End search for initial values at",date(),"\n")
  order <- array(0,ddim)
#  logarithmic eigen values
  mix <- array(0,c(maxcomp,ddim))
  orient <- array(0,c(2,maxcomp,ddim))
  sigma2 <- apply(siq,1:3,var)
  n1 <- ddim[1]
  n2 <- ddim[2]
  n3 <- ddim[3]
  igc <- 0
  ingc <- 0
  prt0 <- proc.time()[1]
  prta <- prt0-prta
  for(i1 in 1:n1) for(i2 in 1:n2) for(i3 in 1:n3){
     if(mask[i1,i2,i3]){
     mc0 <- maxcomp
     ord <- mc0+1
#     for(j in 1:mc0) orient[,j,i1,i2,i3] <- paroforient(grad[siind[j+1,i1,i2,i3],])
     for(j in 1:mc0) orient[,j,i1,i2,i3] <- paroforient(icosa3$vertices[,siind[j+2,i1,i2,i3]])
#
#   these are the gradient vectors corresponding to minima in spherical coordinates
#
     if(method=="mixtensor"){
     par <- numeric(2*mc0+1)
     par[1] <- th[siind[2,i1,i2,i3]]#*lev[2,i1,i2,i3]/mlev
#
#  these is an initial estimate for the eigen-value parameter
#
     par[rep(2*(1:mc0),rep(2,mc0))+c(0,1)] <- orient[,1:mc0,i1,i2,i3]
     } 
     sigmai <- sigma2[i1,i2,i3]
     krit <- sigmai*(ngrad0-1)
     maxcomp0 <- maxcomp
     for(k in mc0:1){
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
        }         
        value <- z$value
#
#   estimate of sigma from the best fitting model
#
        if(method=="mixtensor"){
            zz <- mfunplwghts0(z$par[1:lpar],siq[i1,i2,i3,],grad)
        } 
        ord <- zz$ord
        sigmai <- min(value/(ngrad0-3*ord-1),sigmai)
        if(any(zz$lev<0)||ord<k){
           ttt <- krit
           maxcomp0 <- k-1
#   parameters not interpretable reduce order
        } else {
           penalty <- penIC*(3*ord+1)
           ttt <- value+penalty*sigmai
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
           sigma2[i1,i2,i3] <- sigmai
       }
     }
   }
    if(igc<ngc){
       igc <- igc+1
    } else {
       igc <- 1
       ingc <- ingc+1
       prt1 <- proc.time()[1]
       gc()
       cat("Nr. of voxel",ingc*ngc,"time elapsed:",prta+prt1-prt0,"remaining time:",
            (prt1-prt0)/(ingc*ngc)*(sum(mask)-ingc*ngc),"\n")
    }
  }
  }
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
