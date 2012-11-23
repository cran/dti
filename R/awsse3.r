# This file contains the implementation of dti.smooth() for 
# "dtiData" lines 12--
# and 
# "dtiTensor" lines 225--
# objects. The latter can be done with 
# "Euclidian Metric" lines 237--
# "Riemann Metric" lines 414--

dwi.smooth <- function(object, ...) cat("No DTI smoothing defined for this class:",class(object),"\n")

setGeneric("dwi.smooth", function(object, ...) standardGeneric("dwi.smooth"))

setMethod("dwi.smooth", "dtiData", function(object,kstar,lambda=6,kappa0=NULL,ncoils=1,sigma=NULL,sigma2=NULL,level=NULL,minsb=5,vred=4,xind=NULL,yind=NULL,zind=NULL,verbose=FALSE,dist=1,model="Chi2"){
  args <- sys.call(-1)
  args <- c(object@call,args)
  sdcoef <- object@sdcoef
  vext <- object@voxelext[2:3]/object@voxelext[1]
  if(length(sigma)==1) {
     cat("using supplied sigma",sigma,"\n")
  } else {
  mask <- getmask(object,level)$mask
  sigma <- numeric(object@ngrad)
  for(i in 1:object@ngrad){
  sigma[i] <- awssigmc(object@si[,,,i],12,mask,ncoils,vext,h0=1.25,verbose=verbose)$sigma
  cat("image ",i," estimated sigma",sigma[i],"\n")
  }
 cat("quantiles of estimated sigma values",quantile(sigma),"\n")
  sigma <- median(sigma)
 cat("using median estimated sigma",sigma,"\n")
  }
  model <- if(model=="Chi2") 1 else 0 
  if(!(is.null(xind)&is.null(yind)&is.null(zind))){
  if(is.null(xind)) xind <- 1:object@ddim[1]
  if(is.null(yind)) yind <- 1:object@ddim[2]
  if(is.null(zind)) zind <- 1:object@ddim[3]
  object <- object[xind,yind,zind]
  }
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  ngrad <- ngrad - ns0
  grad <- object@gradient[,-s0ind]
  sb <- object@si[,,,-s0ind]
  s0 <- object@si[,,,s0ind]
  if(is.null(kappa0)){
#  select kappa based on variance reduction on the sphere
   if(is.null(vred)) {
     warning("You need to specify either kappa0 or vred\n returning unsmoothed object")
     return(object)
   }
   if(!is.numeric(vred)|vred<1){
     warning("vred needs to be >= 1\n returning unsmoothed object")
     return(object)
   }
   kappa0 <- suggestkappa(grad,vred,dist)$kappa
  }
  if(model==1){
#
#   use squared values for Chi^2
#
     sb <- sb^2
     s0 <- s0^2
     sigma <- sigma^2
  }
  if(ns0>1){
     dim(s0) <- c(prod(ddim),ns0)
     s0 <- s0%*%rep(1/ns0,ns0)
     dim(s0) <- ddim
  }
     th0 <- s0
     ni0 <- array(1,ddim)
  mask <- apply(sb,1:3,mean) > minsb
  masksb <- array(mask,c(ddim,ngrad))
  gradstats <- getkappas(grad,dist=dist)
  hseq <- gethseqfullse3(kstar,gradstats,kappa0,vext=vext)
  kappa <- hseq$kappa
  nind <- hseq$n
  hseq <- hseq$h
# make it nonrestrictive for the first step
  ni <- array(1,dim(sb))
  z <- list(th=sb, ni = ni)
  prt0 <- Sys.time()
  cat("adaptive smoothing in SE3, kstar=",kstar,if(verbose)"\n" else " ")
  kinit <- if(lambda<1e10) 1 else kstar
  mc.cores <- setCores(,reprt=FALSE)
  for(k in kinit:kstar){
    gc()
    hakt <- hseq[,k]
    param <- lkfullse3(hakt,kappa/hakt,gradstats,vext,nind) 
    if(length(sigma)==1) {
    z <- .Fortran("adsmse3p",
                as.single(sb),
                as.single(z$th),
                ni=as.single(z$ni),
                as.logical(mask),
                as.integer(ddim[1]),
                as.integer(ddim[2]),
                as.integer(ddim[3]),
                as.integer(ngrad),
                as.double(lambda),
                as.integer(ncoils),
                as.integer(mc.cores),
                as.integer(param$ind),
                as.double(param$w),
                as.integer(param$n),
                th=single(prod(ddim)*ngrad),
                single(prod(ddim)*ngrad),#ldf (to precompute lgamma)
                as.double(sigma),
                double(ngrad*mc.cores),
                double(ngrad*mc.cores),
                as.integer(model),
                DUPL=FALSE,
                PACKAGE="dti")[c("ni","th")]
    gc()
    } else {
    warning("not yet implemented for heterogenious variances\n
             returning original object")
    return(object)
    }
if(verbose){
   dim(z$ni) <- c(prod(ddim),ngrad)
   cat("k:",k,"h_k:",signif(max(hakt),3)," quartiles of ni",signif(quantile(z$ni[mask,]),3),
  "mean of ni",signif(mean(z$ni[mask,]),3),
  " time elapsed:",format(difftime(Sys.time(),prt0),digits=3),"\n")
  } else {
      cat(".")
    }
      param <- reduceparam(param)
     z0 <- .Fortran("asmse30p",
                    as.double(s0),
                    as.double(th0),
                    ni=as.double(ni0),
                    as.logical(mask),
                    as.integer(ddim[1]),
                    as.integer(ddim[2]),
                    as.integer(ddim[3]),
                    as.integer(ns0),
                    as.double(lambda),
                    as.integer(ncoils),
                    as.integer(param$ind),
                    as.double(param$w),
                    as.integer(param$n),
                    as.integer(param$starts),
                    as.integer(param$nstarts),
                    th0=double(prod(ddim)),
                    double(prod(ddim)),#ldf (to precompute lgamma)
                    as.double(sigma),
                    double(param$nstarts),#swi
                    as.integer(model),
                    DUPL=FALSE,
                    PACKAGE="dti")[c("ni","th0")]
      th0 <- z0$th0
      ni0 <- z0$ni
      rm(z0)
      gc()
if(verbose){
   cat("End smoothing s0: quartiles of ni",signif(quantile(ni0[mask]),3),
  "mean of ni",signif(mean(ni0[mask]),3),
  " time elapsed:",format(difftime(Sys.time(),prt0),digits=3),"\n")
  }
  }
  ngrad <- ngrad+1
  si <- array(0,c(ddim,ngrad))
  si[,,,1] <-  th0 
  si[,,,-1] <- z$th
  object@si <- if(model==1) sqrt(si) else si
  object@gradient <- grad <- cbind(c(0,0,0),grad)
  object@bvalue <- c(0,object@bvalue[-object@s0ind])
  object@btb <- create.designmatrix.dti(grad)
  object@s0ind <- as.integer(1)
  object@replind <- as.integer(1:ngrad)
  object@ngrad <- as.integer(ngrad)
  object@call <- args
  object
}
)

koayinv <- function(r,th0,eps=1e-6){
eps <- max(max(1e-8,eps))
if(r<400){
r <- max(sqrt(pi/(4-pi))+eps,r)
thsq <- th0*th0
db <-(2+thsq)*besselI(thsq/4,0,TRUE)+thsq*besselI(thsq/4,1,TRUE)
ksi <- 2+thsq-pi/8*db*db
th1 <- sqrt(ksi*(1+r*r)-2)
while(abs(th0-th1)>eps){
th0 <- th1
thsq <- th0*th0
db <-(2+thsq)*besselI(thsq/4,0,TRUE)+thsq*besselI(thsq/4,1,TRUE)
ksi <- 2+thsq-pi/8*db*db
th1 <- sqrt(ksi*(1+r*r)-2)
}
} else {
th1 <- r*0.9999953
}
th1
}
lkfullse3 <- function(h,kappa,gradstats,vext,n){
      ngrad <- dim(gradstats$bghat)[2]
      if(length(h)<ngrad) h <- rep(h[1],ngrad)
      z <- .Fortran("lkfulse3",
                    as.double(h),
                    as.double(kappa),
                    as.double(gradstats$k456),
                    as.integer(ngrad),
                    as.double(vext),
                    ind=integer(5*n),
                    w=double(n),
                    n=as.integer(n),
                    as.integer(gradstats$dist),
                    DUPL=FALSE,
                    PACKAGE="dti")[c("ind","w","n")]
      dim(z$ind) <- c(5,n)
list(h=h,kappa=kappa,ind=z$ind[,1:z$n],w=z$w[1:z$n],nind=z$n)
}


