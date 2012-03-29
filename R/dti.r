################################################################
#                                                              #
# Section for DTI functions                                    #
#                                                              #
################################################################

dtiTensor <- function(object,  ...) cat("No DTI tensor calculation defined for this class:",class(object),"\n")

setGeneric("dtiTensor", function(object,  ...) standardGeneric("dtiTensor"))

setMethod("dtiTensor","dtiData",function(object, method="nonlinear",varmethod="replicates",varmodel="local",mc.cores=getOption("mc.cores", 2L)) {
#  available methods are 
#  "linear" - use linearized model (log-transformed)
#  "nonlinear" - use nonlinear model with parametrization according to Koay et.al. (2006)
  if(is.null(mc.cores)) mc.cores <- 1
  mc.cores <- min(mc.cores,detectCores())
  args <- sys.call(-1)
  args <- c(object@call,args)
  ngrad <- object@ngrad
  grad <- object@gradient
  ddim <- object@ddim
  nvox <- prod(ddim)
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  sdcoef <- object@sdcoef
  require(parallel)
  if(all(sdcoef[1:4]==0)) {
    cat("No parameters for model of error standard deviation found\n estimating these parameters\n You may prefer to run sdpar before calling dtiTensor")
    sdcoef <- sdpar(object,interactive=FALSE)@sdcoef
  }
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
  si <- array(z$si,c(ddim,ngrad))
  index <- if(z$lindex>0) z$index[1:z$lindex] else numeric(0)
  rm(z)
  gc()
  if(method=="linear"){
     ngrad0 <- ngrad - length(s0ind)
     s0 <- si[,,,s0ind]
     si <- si[,,,-s0ind]
     if(ns0>1) {
         dim(s0) <- c(prod(ddim),ns0)
         s0 <- s0 %*% rep(1/ns0,ns0)
         dim(s0) <- ddim
     }
     mask <- s0 > object@level
     mask <- connect.mask(mask)
     dim(s0) <- dim(si) <- NULL
     ttt <- -log(si/s0)
     ttt[is.na(ttt)] <- 0
     ttt[(ttt == Inf)] <- 0
     ttt[(ttt == -Inf)] <- 0
     dim(ttt) <- c(prod(ddim),ngrad0)
     ttt <- t(ttt)
     cat("Data transformation completed ",format(Sys.time()),"\n")

     btbsvd <- svd(object@btb[,-s0ind])
     solvebtb <- btbsvd$u %*% diag(1/btbsvd$d) %*% t(btbsvd$v)
     D <- solvebtb%*% ttt
     cat("Diffusion tensors generated ",format(Sys.time()),"\n")

     res <- ttt - t(object@btb[,-s0ind]) %*% D
     rss <- res[1,]^2
     for(i in 2:ngrad0) rss <- rss + res[i,]^2
     dim(rss) <- ddim
     sigma2 <- rss/(ngrad0-6)
     D[c(1,4,6),!mask] <- 1e-6
     D[c(2,3,5),!mask] <- 0
#  replace non-tensors (with negative eigenvalues) by a small isotropic tensor 
#      ind <- array(.Fortran("dti3Dev",
#                           as.double(D),
#                           as.integer(nvox),
#                           as.logical(mask),
#                           ev=double(3*prod(ddim)),
#                           DUP=FALSE,
#                           PACKAGE="dti")$ev,c(3,ddim))[1,,,]<1e-6
#       if(sum(ind&mask)>0){
#           D[c(1,4,6),ind&mask] <- 1e-6
#           D[c(2,3,5),ind&mask] <- 0
#       }
#  Replace tensors with negative eigenvalues
     D <- .Fortran("dti3Dreg",
                   D=as.double(D),
                   as.integer(prod(ddim)),
                   DUP=FALSE,
                   PACKAGE="dti")$D                   
     dim(D) <- c(6,ddim)
     dim(res) <- c(ngrad0,ddim)
     cat("Variance estimates generated ",format(Sys.time()),"\n")
     th0 <- array(s0,object@ddim)
     th0[!mask] <- 0
     gc()
  } else {
#  method == "nonlinear" 
     ngrad0 <- ngrad
     si <- aperm(si,c(4,1:3))
     s0 <- si[s0ind,,,]
     if(ns0>1) {
         dim(s0) <- c(ns0,nvox)
         s0 <- rep(1/ns0,ns0)%*%s0
         dim(s0) <- ddim
     }
     mask <- s0 > object@level
     mask <- connect.mask(mask)
     df <- sum(table(object@replind)-1)
     cat("start nonlinear regression",format(Sys.time()),"\n")
     if(mc.cores==1){
     z <- .Fortran("nlrdtirg",
                as.integer(si),
                as.integer(ngrad),
                as.integer(nvox),
                as.logical(mask),
                as.double(object@btb),
                as.double(sdcoef),
                th0=as.double(s0),
                D=double(6*nvox),
                as.integer(200),
                as.double(1e-6),
                res=double(ngrad*nvox),
                rss=double(nvox),
                double(ngrad),
                PACKAGE="dti",DUP=FALSE)[c("th0","D","res","rss")]
     res <- z$res
     D <- z$D
     rss <- z$rss
     th0 <- z$th0
     } else {
        dim(si) <- c(dim(si)[1],nvox)
        z <- matrix(0,8+ngrad,nvox)
        z[,mask] <- pmatrix(si[,mask],pnlrdtirg,btb=object@btb,sdcoef=sdcoef,s0ind=s0ind,
                     mc.cores=mc.cores,mc.silent = TRUE)
        th0 <- z[1,]
        D <- z[2:7,]
        rss <- z[8,]
        res <- z[8+(1:ngrad),]
     }
     dim(th0) <- ddim
     dim(D) <- c(6,ddim)
     dim(res) <- c(ngrad,ddim)
     dim(rss) <- ddim
#  handle points where estimation failed
     n <- prod(ddim)
     dim(mask) <- c(1, ddim)
     indD <- (1:n)[D[2,,,, drop=FALSE]==0&D[3,,,, drop=FALSE]==0&D[5,,,, drop=FALSE]==0&mask]
     dim(mask) <- ddim
# this does not work in case of 2D data: 
#     indD <- (1:n)[D[2,,,]==0&D[3,,,]==0&D[5,,,]==0&mask]
     if(length(indD)>0){
     dim(si) <- c(ngrad,n)
     dim(D) <- c(6,n)
     dim(res) <- c(ngrad,n)
#     cat("length of IndD",length(indD),"\n")
     for(i in indD){
        zz <- optim(c(1,0,0,1,0,1),opttensR,method="BFGS",si=si[-s0ind,i],s0=s0[i],grad=grad[,-s0ind],sdcoef=sdcoef)
        D[,i] <- rho2D(zz$par)
        th0[i] <- s0[i]
        rss[i] <- zz$value
        res[s0ind,i] <- 0
        res[-s0ind,i] <- tensRres(zz$par,si[-s0ind,i],s0[i],grad[,-s0ind])
     }
     dim(D) <- c(6,ddim)
     dim(res) <- c(ngrad,ddim)
     }
     cat("successfully completed nonlinear regression ",format(Sys.time()),"\n")
     sigma2 <- array(0,c(1,1,1))
     rm(z)
     gc()
  }
  lags <- c(5,5,3)
  scorr <- .Fortran("mcorr",as.double(res),
                   as.logical(mask),
                   as.integer(ddim[1]),
                   as.integer(ddim[2]),
                   as.integer(ddim[3]),
                   as.integer(ngrad0),
                   double(prod(ddim)),
                   double(prod(ddim)),
                   scorr = double(prod(lags)),
                   as.integer(lags[1]),
                   as.integer(lags[2]),
                   as.integer(lags[3]),
                   PACKAGE="dti",DUP=FALSE)$scorr
  dim(scorr) <- lags
  scorr[is.na(scorr)] <- 0
  cat("estimated spatial correlations",format(Sys.time()),"\n")
  cat("first order  correlation in x-direction",signif(scorr[2,1,1],3),"\n")
  cat("first order  correlation in y-direction",signif(scorr[1,2,1],3),"\n")
  cat("first order  correlation in z-direction",signif(scorr[1,1,2],3),"\n")

  scorr[is.na(scorr)] <- 0
  bw <- optim(c(2,2,2),corrrisk,method="L-BFGS-B",lower=c(.2,.2,.2),
  upper=c(3,3,3),lag=lags,data=scorr)$par
  bw[bw <= .25] <- 0
  cat("estimated corresponding bandwidths",format(Sys.time()),"\n")
  if(mc.cores<=1) ev <- .Fortran("dti3Dev",
                       as.double(D),
                       as.integer(nvox),
                       as.logical(mask),
                       ev=double(3*nvox),
                       DUP=FALSE,
                       PACKAGE="dti")$ev
   else {
      ev <- matrix(0,3,prod(ddim))
      dim(D) <- c(6,prod(ddim))
      ev[,mask] <- pmatrix(D[,mask],pdti3Dev,mc.cores=mc.cores,mc.silent = TRUE)
   }
  dim(ev) <- c(3,ddim)   
  dim(D) <- c(6,ddim)   
  scale <- quantile(ev[3,,,][mask],.95)
  cat("estimated scale information",format(Sys.time()),"\n")  
  invisible(new("dtiTensor",
                call  = args,
                D     = D,
                th0   = th0,
                sigma = sigma2,
                scorr = scorr, 
                bw = bw, 
                mask = mask,
                hmax = 1,
                gradient = object@gradient,
                btb   = object@btb,
                ngrad = object@ngrad, # = dim(btb)[2]
                s0ind = object@s0ind,
                replind = object@replind,
                ddim  = object@ddim,
                ddim0 = object@ddim0,
                xind  = object@xind,
                yind  = object@yind,
                zind  = object@zind,
                voxelext = object@voxelext,
                level = object@level,
                orientation = object@orientation,
                rotation = object@rotation,
                source = object@source,
                outlier = index,
                scale = scale,
                method = method)
            )
})

opttensR <- function(param,si,s0,grad,sdcoef){
      .Fortran("opttensR",
               as.double(param),
               as.double(si),
               as.double(s0),
               as.double(grad),
               as.integer(length(si)),
               as.double(sdcoef),
               erg=double(1),
               DUP=FALSE,
               PACKAGE="dti")$erg
}
tensRres <- function(param,si,s0,grad){
      .Fortran("tensRres",
               as.double(param),
               as.double(si),
               as.double(s0),
               as.double(grad),
               as.integer(length(si)),
               res=double(length(si)),
               DUP=FALSE,
               PACKAGE="dti")$res
}
rho2D <- function(param){
      .Fortran("rho2D0",
               as.double(param),
               D=double(6),
               DUP=FALSE,
               PACKAGE="dti")$D
}
D2rho <- function(D){
      .Fortran("D2rho0",
               as.double(D),
               rho=double(6),
               DUP=FALSE,
               PACKAGE="dti")$rho
}
#############

dtiIndices <- function(object, ...) cat("No DTI indices calculation defined for this class:",class(object),"\n")

setGeneric("dtiIndices", function(object, ...) standardGeneric("dtiIndices"))

setMethod("dtiIndices","dtiTensor",
function(object, which, mc.cores=getOption("mc.cores", 2L)) {
  args <- sys.call(-1)
  args <- c(object@call,args)
  ddim <- object@ddim
  n <- prod(ddim)
  z <- if(mc.cores==1) .Fortran("dtiind3D",
                as.double(object@D),
                as.integer(n),
                as.logical(object@mask),
                fa=double(n),
                ga=double(n),
                md=double(n),
                andir=double(3*n),
                bary=double(3*n),
                DUP=FALSE,
                PACKAGE="dti")[c("fa","ga","md","andir","bary")] else {
           D <- matrix(object@D,6,prod(ddim))[,object@mask]
           res <- matrix(0,9,prod(ddim))
           res[,object@mask] <- pmatrix(D,pdtiind3D,mc.cores=mc.cores,mc.silent = TRUE)
           list(andir=res[1:3,],
                   fa=res[4,],
                   ga=res[5,],
                   md=res[6,],
                 bary=res[7:9,])
                }

  invisible(new("dtiIndices",
                call = args,
                fa = array(z$fa,ddim),
                ga = array(z$ga,ddim),
                md = array(z$md,ddim),
                andir = array(z$andir,c(3,ddim)),
                bary = array(z$bary,c(3,ddim)),
                gradient = object@gradient,
                btb   = object@btb,
                ngrad = object@ngrad, # = dim(btb)[2]
                s0ind = object@s0ind,
                ddim  = ddim,
                ddim0 = object@ddim0,
                voxelext = object@voxelext,
                orientation = object@orientation,
                rotation = object@rotation,
                xind  = object@xind,
                yind  = object@yind,
                zind  = object@zind,
                method = object@method,
                level = object@level,
                source= object@source)
            )
})

