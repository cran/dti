#
#  Misc. functions
#

getsiind <- function(si,mask,grad,maxcomp=3,maxc=.866){
# assumes dim(grad) == c(ngrad,3)
# assumes dim(si) == c(n1,n2,n3,ngrad)
# SO removed
ngrad <- dim(grad)[1]
dgrad <- matrix(abs(grad%*%t(grad)),ngrad,ngrad)
dgrad <- dgrad/max(dgrad)
# this provides a matrix of 1-cos(alpha) where alpha
# is the angle between gradient directions
array(.Fortran("getsiind",
         as.double(aperm(si,c(4,1:3))),
         as.integer(ngrad),
         as.integer(dim(si)[1]),
         as.integer(dim(si)[2]),
         as.integer(dim(si)[3]),
         as.integer(maxcomp),
         as.double(maxc),
         as.double(dgrad),
         double(ngrad),
         siind=integer((maxcomp+1)*prod(dim(si)[-4])),
         as.integer(maxcomp+1),
         as.logical(mask),
         PACKAGE="dti")$siind,c(maxcomp+1,dim(si)[-4]))
}


orientofpar <- function(par){
c(sin(par[1])*cos(par[2]),sin(par[1])*sin(par[2]),cos(par[1]))
}

dwiMixtensor.old <- function(object, ...) cat("No dwiMixtensor calculation defined for this class:",class(object),"\n")

setGeneric("dwiMixtensor.old", function(object,  ...) standardGeneric("dwiMixtensor.old"))

setMethod("dwiMixtensor.old", "dtiData", function(object, maxcomp=2, p=2, maxneighb=7, method="mixtensor", reltol=1e-8, maxit=5000) {

  cat("entering method dwiMixtensor (old C-based version)\n")

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
                as.double(object@si),
                as.integer(prod(ddim)),
                as.integer(ngrad),
                as.logical((1:ngrad)%in%s0ind),
                as.integer(ns0),
                si      = integer(prod(ddim)*ngrad),
                index   = integer(prod(ddim)),
                lindex  = integer(1),
                DUP     = FALSE,
                PACKAGE = "dti")[c("si","index","lindex")]
  si <- array(z$si, c(ddim, ngrad))
  index <- if (z$lindex>0) z$index[1:z$lindex] else numeric(0)
  rm(z)
  cat("done\n")

  cat("prepare initial estimates ... ")
  # prepare data for optim
  s0 <- si[,,,s0ind]
  if (length(s0ind)>1) s0 <- apply(s0, 1:3, mean)
  # normalized DW data
  mask <- s0 > object@level
  siq <- si[,,,-s0ind]
  dim(siq) <- c(prod(ddim),ngrad-ns0)
  siq[mask,] <- sweep(siq[mask,],1,s0[mask],"/")
  dim(siq) <- c(ddim,ngrad-ns0)
  # heuristics to avoid DWI that are larger than s0.
  siqmed <- apply(siq,1:3,median)
  siqmed[siqmed<.9] <- .9
  siqmed[siqmed>.99] <- .99
  siq <- sweep(siq,1:3,siqmed,pmin)
  # mask for calculation
  grad <- t(object@gradient[,-s0ind])
#
#   determine initial estimates for orientations 
#
  siind <- getsiind(siq,mask,grad,maxcomp,maxc=.866)
  cat("done\n")

  cat("optimizing ... ")
  mm <- switch(method, "mixtensor" = 1,
                       "Jian"      = 2,
                       1)
  pl <- if (method == "Jian2") prod(ddim) else 1;

  # perform voxelwise optimization and order selection of tensor mixture model
  a <- .C("mixture",
          as.integer(mm),                         # select mixture method
#          as.integer(ddim[1]),                    # voxel dim x
#          as.integer(ddim[2]),                    # voxel dim y
#          as.integer(ddim[3]),                    # voxel dim z
          as.integer(prod(ddim)),                 # number of voxels
          as.integer(mask),                       # calculation mask
          as.double(siq),                         # DWI without s0
          as.integer(siind),                      # DWI indices of local minima
          as.integer(ngrad - ns0),      # number of DWI
          as.double(grad),                        # gradient directions
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
          DUP     = FALSE,
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
                scorr       = pl, # ???
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
                rotation    = object@rotation,
                source      = object@source,
                outlier     = index,
                scale       = 1,
                method      = method)
            )
}
)
