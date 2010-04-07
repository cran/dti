

#
#  Misc. functions
#

neighbors <- function(grad, m=7) {
  ngrad <- dim(grad)[1]
  z <- apply(as.matrix(dist(rbind(grad,-grad)))[1:ngrad,],1,order)[1:m,]
  # also consider opposite directions
  z[z>ngrad] <- z[z>ngrad]-ngrad
  # identify opposite directions
  z
}

locmin <- function(si,sneighbors,grad){
  ngrad <- length(si)
  z <- si[sneighbors]
  dim(z) <- dim(sneighbors)
  ind <- (1:ngrad)[apply(t(z)>=si,1,all)]
  ind <- ind[order(si[ind])]
  list(ind = ind, si = si[ind], grad = grad[ind,,drop=FALSE])
}

locmin2 <- function(si, sneighbors){
  ddim <- dim(si)[1:3]
  ngrad <- dim(si)[4]
  si <- .C("smoothsi",
           as.integer(ddim[1]),
           as.integer(ddim[2]),
           as.integer(ddim[3]),
           as.double(si),
           as.integer(ngrad),
           as.integer(sneighbors),
           as.integer(dim(sneighbors)[1]),
           si = double(length(si)),
           PACKAGE = "dti")$si
  dim(si) <- c(ddim,ngrad)
  minsi <- function(x, sneighbors) {
    z <- x[sneighbors]
    dim(z) <- dim(sneighbors)
    ind <- (1:ngrad)[apply(t(z) >= x, 1, all)]
    si2 <- rep(1, dim(sneighbors)[2])
    si2[ind] <- x[ind]
    si2
  }
  s <- apply(si, 1:3, minsi, sneighbors)
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

getw <- function(cw){
   m <- length(cw)
   ew <- exp(cw)
   z<-ew/(1+ew)
   if(m>1){
   for(i in 2:m) z[i] <- (1-sum(z[1:(i-1)]))*z[i]
   }
   c(z,1-sum(z))
}



mfun <- function(par,siq,grad){
lpar <- length(par)
m <- (lpar-1)/3
ngrad <- dim(grad)[1]
.Fortran("mfun",as.double(par),
                as.double(siq),
                as.double(t(grad)),
                as.integer(m),
                as.integer(lpar),
                as.integer(ngrad),
                double(ngrad),
                erg = double(1),
                PACKAGE="dti")$erg
}

mfun2 <- function(par,siq,grad,ep=50){
lpar <- length(par)
m <- (lpar-1)/3
ngrad <- dim(grad)[1]
.Fortran("mfun2",as.double(par),
                as.double(siq),
                as.double(t(grad)),
                as.integer(m),
                as.double(ep),
                as.integer(lpar),
                as.integer(ngrad),
                double(ngrad),
                erg = double(1),
                PACKAGE="dti")$erg
}
mfun3 <- function(par,siq,grad){
lpar <- length(par)
m <- (lpar-2)/3
ngrad <- dim(grad)[1]
ep <- max(1,min(par[lpar],50))
.Fortran("mfun2",as.double(par[-lpar]),
                as.double(siq),
                as.double(t(grad)),
                as.integer(m),
                as.double(ep),
                as.integer(lpar-1),
                as.integer(ngrad),
                double(ngrad),
                erg = double(1),
                PACKAGE="dti")$erg
}


dwiMixtensor.old <- function(object, ...) cat("No dwiMixtensor.old calculation defined for this class:",class(object),"\n")

setGeneric("dwiMixtensor.old", function(object,  ...) standardGeneric("dwiMixtensor.old"))

setMethod("dwiMixtensor.old","dtiData",function(object, maxcomp=2, p=40, maxneighb=7, method="mixtensor", reltol=1e-8, maxit=5000){
  args <- sys.call(-1)
  args <- c(object@call,args)
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
#  sdcoef <- object@sdcoef
#  if(all(sdcoef==0)) {
#    cat("No parameters for model of error standard deviation found\n estimating these parameters\n You may prefer to run sdpar before calling dtiTensor")
#    sdcoef <- sdpar(object,interactive=FALSE)@sdcoef
#  }
  z <- .Fortran("outlier",
                as.integer(object@si),
                as.integer(prod(ddim)),
                as.integer(ngrad),
                as.logical((1:ngrad)%in%s0ind),
                as.integer(ns0),
                si=integer(prod(ddim)*ngrad),
                index=integer(prod(ddim)),
                lindex=integer(1),
                DUP=FALSE,
                PACKAGE="dti")[c("si","index","lindex")]
  si <- array(z$si,c(ddim,ngrad))
  index <- if(z$lindex>0) z$index[1:z$lindex] else numeric(0)
  rm(z)
  ngrad0 <- ngrad - length(s0ind)
  s0 <- si[,,,s0ind]
  if(length(s0ind)>1) s0 <- apply(s0,1:3,mean)
  mask <- s0>0
  si <- si[,,,-s0ind]
  siq <- sweep(si,1:3,s0,"/")
  siqmed <- apply(siq,1:3,median)
  siqmed[siqmed<.9] <- .9
  siqmed[siqmed>.99] <- .99
  siq <- sweep(siq,1:3,siqmed,pmin)
#
#  avoid situations where si's are larger than s0
#
  grad <- t(object@gradient[,-s0ind])
  sneighbors <- neighbors(grad,maxneighb)
  order <- array(maxcomp,ddim)
  lev <- array(as.integer(0),c(2,ddim))
#  logarithmic eigen values
  mix <- array(0,c(maxcomp,ddim))
  orient <- array(0,c(2,maxcomp,ddim))
  sigma2 <- array(0,ddim)
  n1 <- ddim[1]
  n2 <- ddim[2]
  n3 <- ddim[3]
  if(method=="Jian2")  p <- array(50,ddim) 
  swghts <- 1/(1:maxneighb)
  for(i1 in 1:n1) for(i2 in 1:n2) for(i3 in 1:n3){
     if(mask[i1,i2,i3]){
 #  cat("voxel",i1,i2,i3,"\n")
# smooth the si's by weighted nearest neighbors 
# to gain more stable locations of minima
     smsi <- .Fortran("smsi",
                      as.double(si[i1,i2,i3,]),
                      as.integer(ngrad0),
                      as.integer(sneighbors),
                      as.integer(maxneighb),
                      as.double(swghts),
                      smsi=double(ngrad0),
                      DUP=FALSE,
                      PACKAGE="dti")$smsi
#     z <- locmin(si[i1,i2,i3,],sneighbors,grad)
     z <- locmin(smsi,sneighbors,grad)
     mc0 <- min(maxcomp,dim(z$grad)[1])
#     cat("mc0",mc0,"\n")
     for(j in 1:mc0) orient[,j,i1,i2,i3] <- paroforient(z$grad[j,])
     par <- numeric(3*mc0+1)
     lev[,i1,i2,i3] <- par[1:2] <- log(-log(siq[i1,i2,i3,z$ind[1]])*c(.8,.2))
     par[rep(3*(1:mc0),rep(2,mc0))+c(0,1)] <- orient[,1:mc0,i1,i2,i3] 
     rss <- Inf
     krit <- Inf
     for(k in mc0:1){
        if(k>1) par[3*(2:k)-1] <- -log((k-1):1)
        lpar <- 3*k+1;
        z <-switch(method,
                  "mixtensor" =  optim(par[1:(3*k+1)],mfun,siq=siq[i1,i2,i3,],grad=grad,control=list(maxit=maxit,reltol=reltol)),
#                  "mixtensor" = .C("mixtensor",
#                                   as.integer(lpar),
#                                   as.double(par[1:(3*k+1)]),
#                                   par = double(lpar),
#                                   as.integer(ngrad0),
#                                   as.double(siq[i1,i2,i3,]),
#                                   as.double(as.matrix(grad)),
#                                   as.integer(maxit),
#                                   as.double(reltol),
#                                   value = as.double(1))[c("par", "value")],
                 "Jian"       = optim(par[1:(3*k+1)],mfun2,siq=siq[i1,i2,i3,],grad=grad,ep=p,control=list(maxit=maxit,reltol=reltol)),
                 "Jian2"      = optim(c(par[1:(3*k+1)],25),mfun3,siq=siq[i1,i2,i3,],grad=grad,control=list(maxit=maxit,reltol=reltol)))
        rss <- min(z$value,rss)
        ttt <- z$value+2*(3*k+1)/(ngrad-3*maxcomp-1)*rss
#        cat("risk",z$value,ttt,"\n")
        if(ttt < krit) {
           krit <- ttt
           order[i1,i2,i3] <- as.integer(k)
           lev[,i1,i2,i3] <- z$par[1:2]
           mix[1:k,i1,i2,i3] <- if(k>1) getw(z$par[3*(2:k)-1]) else 1
           or <- matrix(z$par[rep(3*(1:k),rep(2,k))+c(0,1)],2,k)
           or[1,or[1,]<0] <- or[1,or[1,]<0]+pi
           or[1,or[1,]>pi] <- or[1,or[1,]>pi]-pi
           or[2,or[2,]<0] <- or[2,or[2,]<0]+2*pi
           or[2,or[2,]>2*pi] <- or[2,or[2,]>2*pi]-2*pi
           orient[,1:k,i1,i2,i3] <- or
           if(method=="Jian2") p[i1,i2,i3] <- z$par[length(z$par)]
       }
     }
   sigma2[i1,i2,i3] <- rss/(ngrad0-3*mc0-1)
#   cat("order",order[i1,i2,i3],"\n")
#   cat("error variance",sigma2[i1,i2,i3]*s0[i1,i2,i3]^2,"\n")
#   cat("ev",c(exp(lev[1,i1,i2,i3]),0,0)+exp(lev[2,i1,i2,i3]),"\n")
#   cat("mix",mix[,i1,i2,i3],"\n")
#   cat("orient",orient[,,i1,i2,i3],"\n")
#   if(method=="Jian2") cat("p",p[i1,i2,i3],"\n") 
    gc()
  }
  }
  invisible(new("dwiMixtensor",
                call   = args,
                ev     = exp(lev),
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
                source = object@source,
                outlier = index,
                scale = 1,
                method = method)
            )
}
)

dwiMixtensor <- function(object, ...) cat("No dwiMixtensor calculation defined for this class:",class(object),"\n")

setGeneric("dwiMixtensor", function(object,  ...) standardGeneric("dwiMixtensor"))

setMethod("dwiMixtensor", "dtiData", function(object, maxcomp=2, p=2, maxneighb=7, method="mixtensor", reltol=1e-8, maxit=5000) {

  cat("entering method dwiMixtensor\n")

  if (maxcomp > 3) {
    stop("maximum number of tensor components is 3\n");
  }
  args <- sys.call(-1)
  args <- c(object@call, args)
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  ns0 <- length(s0ind)

  cat("determine outliers ... ")
  # determine outlier ????
  z <- .Fortran("outlier",
                as.integer(object@si),
                as.integer(prod(ddim)),
                as.integer(ngrad),
                as.logical((1:ngrad)%in%s0ind),
                as.integer(ns0),
                si      = integer(prod(ddim)*ngrad),
                index   = integer(prod(ddim)),
                lindex  = integer(1),
                DUP   = FALSE,
                PACKAGE = "dti")[c("si","index","lindex")]
  si <- array(z$si, c(ddim, ngrad))
  index <- if (z$lindex>0) z$index[1:z$lindex] else numeric(0)
  rm(z)
  cat("done\n")

  cat("prepare initial estimates ... ")
  # prepare data for optim
  s0 <- si[,,,s0ind]
  si <- si[,,,-s0ind]
  if (length(s0ind)>1) s0 <- apply(s0, 1:3, mean)
  # normalized DW data
  siq <- sweep(si,1:3,s0,"/")
  # heuristics to avoid DWI that are larger than s0.
  siqmed <- apply(siq,1:3,median)
  siqmed[siqmed<.9] <- .9
  siqmed[siqmed>.99] <- .99
  siq <- sweep(siq,1:3,siqmed,pmin)
  # mask for calculation
  mask <- s0 > 0
  grad <- t(object@gradient[,-s0ind])
  sneighbors <- neighbors(grad, maxneighb)
  s1 <- locmin2(siq, sneighbors) # NOTE! dim(s1) == c(ngrad-length(s0ind), ddim)
  ordersi <- function(x, ncomp) order(x)[1:ncomp]
  siind <- numeric( (maxcomp + 1) * prod(ddim))
  dim(siind) <- c( maxcomp + 1, ddim)
  siind[-1,,,] <- apply(s1, 2:4, ordersi, maxcomp)
  countmin <- function(x, ncomp) min(sum(x < 1), ncomp)
  siind[1,,,] <- apply(s1, 2:4, countmin, maxcomp)
  cat("done\n")

  cat("optimizing ... ")
  mm <- switch(method, "mixtensor" = 1,
                       "Jian"      = 2,
#                       "Jian2"     = 3,
                       1)
  pl <- if (method == "Jian2") prod(ddim) else 1;

  # perform voxelwise optimization and order selection of tensor mixture model
  a <- .C("mixture",
          as.integer(mm),                         # select mixture method
          as.integer(ddim[1]),                    # voxel dim x
          as.integer(ddim[2]),                    # voxel dim y
          as.integer(ddim[3]),                    # voxel dim z
          as.integer(mask),                       # calculation mask
          as.double(siq),                         # DWI without s0
          as.integer(siind),                      # DWI indices of local minima
          as.integer(ngrad - length(s0ind)),      # number of DWI
          as.double(grad),                        # gradient directions
          as.integer(sneighbors),                 # (maxneighb x ngrad) matrix of nearest gradient neighbors
          as.integer(maxcomp),                  # max number of gradient neighbors
          as.integer(pl),                         # exp for Jian model
          as.integer(maxit),                      # max number of iterations for optim
          as.double(reltol),                      # reltol crit for optim
          order   = double(prod(ddim)),   # selected order of mixture
          lev     = double(2*prod(ddim)),         # logarithmic eigenvalues
          mix     = double(maxcomp*prod(ddim)),   # mixture weights
          orient  = double(2*maxcomp*prod(ddim)), # phi/theta for all mixture tensors
          p       = double(pl),                   # decay const in Jian
          sigma2  = double(prod(ddim)),           # parameter variance ???
          PACKAGE = "dti")[c("order", "lev", "mix", "orient", "p", "sigma2")]

  # set dimension attr
  dim(a$order) <- ddim;
  dim(a$lev) <- c(2, ddim);
  dim(a$mix) <- c(maxcomp, ddim);
  dim(a$orient) <- c(2, maxcomp, ddim);
  if (method == "Jian2") dim(a$p) <- ddim;
  dim(a$sigma2) <- ddim;
  cat("done\n")

  # create and return new object
  invisible(new("dwiMixtensor",
                call        = args,
                ev          = exp(a$lev),
                mix         = a$mix,
                orient      = a$orient,
                order       = a$order,
                p           = if (method == "Jian2") a$p else p,
                th0         = s0,
                sigma       = a$sigma2,
                scorr       = array(1, c(1,1,1)), # ???
                bw          = c(0,0,0),           # ???
                mask        = mask,
                hmax        = 1,                  # ???
                gradient    = object@gradient,
                btb         = object@btb,
                ngrad       = object@ngrad,
                s0ind       = object@s0ind,
                replind     = object@replind,
                ddim        = object@ddim,
                ddim0       = object@ddim0,
                xind        = object@xind,
                yind        = object@yind,
                zind        = object@zind,
                voxelext    = object@voxelext,
                level       = object@level,
                orientation = object@orientation,
                source      = object@source,
                outlier     = index,
                scale       = 1,
                method      = method)
            )
}
)
