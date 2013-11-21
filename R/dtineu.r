################################################################
#                                                              #
# Section for DTI functions                                    #
#                                                              #
################################################################

dtiTensor <- function(object,  ...) cat("No DTI tensor calculation defined for this class:",class(object),"\n")

setGeneric("dtiTensor", function(object,  ...) standardGeneric("dtiTensor"))

setMethod( "dtiTensor", "dtiData",
           function( object, 
                     method = c( "nonlinear", "linear"),
                     mc.cores = setCores( , reprt = FALSE)) {

## check method! available are: 
##   "linear" - use linearized model (log-transformed)
##   "nonlinear" - use nonlinear model with parametrization according to Koay et.al. (2006)
             method <- match.arg( method)
                          
   if(is.null(mc.cores)) mc.cores <- 1
   mc.cores <- min(mc.cores,detectCores())
   args <- sys.call(-1)
   args <- c(object@call,args)
   ngrad <- object@ngrad
   grad <- object@gradient
   ddim <- object@ddim
   ntotal <- prod(ddim)
   s0ind <- object@s0ind
   ns0 <- length(s0ind)
   sdcoef <- object@sdcoef
   #if(mc.cores>1&method=="nonlinear") require(parallel)
   if(all(sdcoef[1:4]==0)) {
      cat("No parameters for model of error standard deviation found\n
          Using constant weights \n You may prefer to run sdpar before calling dtiTensor\n")
#      sdcoef <- sdpar(object,interactive=FALSE)@sdcoef
       sdcoef <- c(1,0,1,1)
   }
   z <- sioutlier(object@si,s0ind,mc.cores=mc.cores)
#
#  this does not scale well with openMP
#
cat("sioutlier completed\n")
   si <- array(z$si,c(ngrad,ddim))
   index <- z$index
   ngrad0 <- ngrad - length(s0ind)
   s0 <- si[s0ind,,,]
   si <- si[-s0ind,,,]
   if(ns0>1) {
      dim(s0) <- c(ns0,prod(ddim))
      s0 <- rep(1/ns0,ns0)%*%s0
      dim(s0) <- ddim
   }
   mask <- s0 > object@level
   mask <- connect.mask(mask)
   dim(si) <- c(ngrad0,prod(ddim))
   ttt <- array(0,dim(si))
   ttt[,mask] <- -log1p(sweep(si[,mask],2,as.vector(s0[mask]),"/")-1)
#  suggestion by B. Ripley
#   idsi <- 1:length(dim(si))
#   ttt <- -log(sweep(si,idsi[-1],s0,"/"))
   ttt[is.na(ttt)] <- 0
   ttt[(ttt == Inf)] <- 0
   ttt[(ttt == -Inf)] <- 0
   dim(ttt) <- c(ngrad0,prod(ddim))
   cat("Data transformation completed ",format(Sys.time()),"\n")
               
   btbsvd <- svd(object@btb[,-s0ind])
   solvebtb <- btbsvd$u %*% diag(1/btbsvd$d) %*% t(btbsvd$v)
   D <- solvebtb%*% ttt
   cat("Diffusion tensors (linearized model) generated ",format(Sys.time()),"\n")
   D[c(1,4,6),!mask] <- 1e-6
   D[c(2,3,5),!mask] <- 0
   D <- dti3Dreg(D,mc.cores=mc.cores)
   dim(D) <- c(6,ntotal)
   th0 <- array(s0,ntotal)
   th0[!mask] <- 0
   nvox <- sum(mask)
   if(method== "nonlinear"){
#  method == "nonlinear" uses linearized model for initial estimates
       ngrad0 <- ngrad
       df <- sum(table(object@replind)-1)
       param <- matrix(0,7,nvox)
       ms0 <- mean(s0[mask])
       mbv <- mean(object@bvalue[-object@s0ind])
##  use ms0 and mbv to rescale parameters such that they are of comparable magnitude
       param[1,] <- s0[mask]/ms0
       param[-1,] <- D2Rall(D[,mask]*mbv)
## use reparametrization D = R^T R
       sdcoef[-2] <- sdcoef[-2]/ms0# effect of rescaling of signal
       cat("start nonlinear regression",format(Sys.time()),"\n")
       if(mc.cores==1){
          param <- matrix(.C("dtens",
                         as.integer(nvox),
                         param=as.double(param),
                         as.double(matrix(z$si,c(ngrad,ntotal))[,mask]/ms0),
                         as.integer(ngrad),
                         as.double(object@btb/mbv),
                         as.double(sdcoef),
                         as.double(rep(0,ngrad)),#si
                         as.double(rep(1,ngrad)),#var                         
                         as.integer(1000),#maxit
                         as.double(1e-7),#reltol
                         PACKAGE="dti",DUP=TRUE)$param,7,nvox)
       } else {
          x <- matrix(0,ngrad+7,nvox)
          x[1:7,] <- param
          x[-(1:7),] <- matrix(z$si,c(ngrad,ntotal))[,mask]/ms0
          param <- plmatrix(x,ptensnl,ngrad=ngrad,btb=object@btb/mbv,
                            sdcoef=sdcoef,maxit=1000,reltol=1e-7)
       }
       th0[mask] <- param[1,]*ms0
       D[,mask] <- R2Dall(param[-1,])/mbv
       cat("successfully completed nonlinear regression ",format(Sys.time()),"\n")
   }
   res <- matrix(0,ngrad,ntotal)
   z <- .Fortran("tensres",
                 as.double(th0[mask]),
                 as.double(D[,mask]),
                 as.double(matrix(z$si,c(ngrad,ntotal))[,mask]),
                 as.integer(nvox),
                 as.integer(ngrad),
                 as.double(object@btb),
                 res = double(ngrad*nvox),
                 rss = double(nvox),
                 PACKAGE="dti",DUP=TRUE)[c("res","rss")]
    res <- matrix(0,ngrad,ntotal)
    res[,mask] <- z$res
    rss <- numeric(ntotal)
    rss[mask] <- z$rss/(ngrad0-6)
    rm(z)
    dim(th0) <- ddim
    dim(D) <- c(6,ddim)
    dim(res) <- c(ngrad,ddim)
    dim(rss) <- ddim
    dim(mask) <-  ddim
    gc()
#
#   get spatial correlation
#
    if(any(is.na(res))){
       dim(res) <- c(ngrad,ntotal)
       indr <- (1:ntotal)[apply(is.na(res),2,any)]
       cat("NA's in res in voxel",indr,"\n")
       res[,indr] <- 0
    }
    if(any(is.na(D))|any(abs(D)>1e10)){
       dim(D) <- c(6,total)
       indD <- (1:total)[apply(is.na(D),2,any)]
       cat("NA's in D in ", length(indD),"voxel:",indD,"\n")
       D[,indD] <- c(1,0,0,1,0,1)
       mask[indD] <- FALSE
       indD <- (1:ntotal)[apply(abs(D)>1e10,2,any)]
       cat("Inf's in D in", length(indD)," voxel:",indD,"\n")
       D[,indD] <- c(1,0,0,1,0,1)
       mask[indD] <- FALSE
    }
    scorr <- mcorr(res,mask,ddim,ngrad0,lags=c(5,5,3),mc.cores=mc.cores)
    ev <- dti3Dev(D,mask,mc.cores=mc.cores)
    dim(ev) <- c(3,ddim)   
    dim(D) <- c(6,ddim)   
    scale <- quantile(ev[3,,,][mask],.95,na.rm=TRUE)
    cat("estimated scale information",format(Sys.time()),"\n")  
    invisible(new("dtiTensor",
                  call  = args,
                  D     = D,
                  th0   = th0,
                  sigma = rss,
                  scorr = scorr$scorr, 
                  bw = scorr$bw, 
                  mask = mask,
                  hmax = 1,
                  gradient = object@gradient,
                  bvalue = object@bvalue,
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

D2Rall <- function(D){
#
#  used for initial values
# in case of negative eigenvalues this uses c(1,0,0,1,0,1) for initial rho
#
      nvox <- dim(D)[2]
      matrix(.Fortran("D2Rall",
               as.double(D),
               rho=double(6*nvox),
               as.integer(nvox),
               DUP=FALSE,
               PACKAGE="dti")$rho,6,nvox)
}
R2Dall <- function(R){
#
#  used for initial values
# in case of negative eigenvalues this uses c(1,0,0,1,0,1) for initial rho
#
      nvox <- dim(R)[2]
      matrix(.Fortran("R2Dall",
               as.double(R),
               D=double(6*nvox),
               as.integer(nvox),
               DUP=FALSE,
               PACKAGE="dti")$D,6,nvox)
}
#############
#############

dtiIndices <- function(object, ...) cat("No DTI indices calculation defined for this class:",class(object),"\n")

setGeneric("dtiIndices", function(object, ...) standardGeneric("dtiIndices"))

setMethod("dtiIndices","dtiTensor",
function(object, mc.cores=setCores(,reprt=FALSE)) {
  args <- sys.call(-1)
  args <- c(object@call,args)
  ddim <- object@ddim
  n <- prod(ddim)
  z <- dtiind3D(object@D,object@mask,mc.cores=mc.cores)
  invisible(new("dtiIndices",
                call = args,
                fa = array(z$fa,ddim),
                ga = array(z$ga,ddim),
                md = array(z$md,ddim),
                andir = array(z$andir,c(3,ddim)),
                bary = array(z$bary,c(3,ddim)),
                gradient = object@gradient,
                bvalue = object@bvalue,
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


