################################################################
#                                                              #
# Section for Utility functions                                #
#                                                              #
################################################################

sdpar <- function(object,  ...) cat("No method defined for class:",class(object),"\n")

setGeneric("sdpar", function(object,  ...) standardGeneric("sdpar"))

setMethod("sdpar","dtiData",function(object,level=NULL,sdmethod="sd",interactive=TRUE,threshfactor=1){
  # determine interval of linearity
  if(!(sdmethod%in%c("sd","mad"))){
    warning("sdmethod needs to be either 'sd' or 'mad'")
    return(object)
  }
  if(is.null(level)) level <- object@level
  s0ind<-object@s0ind
  s0 <- object@si[,,,s0ind]
  ls0ind <- length(s0ind)
  A0 <- level
  if(ls0ind>1) {
    dim(s0) <- c(prod(object@ddim),ls0ind)
    s0mean <- s0%*%rep(1/ls0ind,ls0ind)
    A1 <- quantile(s0mean[s0mean>0],.98)
    dim(s0mean) <- object@ddim
  } else {
    A1 <- quantile(s0[s0>0],.98)
  }
  if(interactive) {
    accept <- FALSE
    ddim <- object@ddim
    bw <- min(bw.nrd(if(ls0ind>1) s0mean[s0mean>0] else s0[s0>0]),diff(range(if(ls0ind>1) s0mean else s0))/256)
    z <- density(if(ls0ind>1) s0mean[s0mean>0&s0mean<A1] else s0[s0>0&s0<A1],bw = bw,n=1024)
    indx1 <- trunc(0.05*ddim[1]):trunc(0.95*ddim[1])
    indx2 <- trunc(0.1*ddim[1]):trunc(0.9*ddim[1])
    indx3 <- trunc(0.15*ddim[1]):trunc(0.85*ddim[1])
    indy1 <- trunc(0.05*ddim[2]):trunc(0.95*ddim[2])
    indy2 <- trunc(0.1*ddim[2]):trunc(0.9*ddim[2])
    indy3 <- trunc(0.15*ddim[2]):trunc(0.85*ddim[2])
    indz1 <- trunc(0.05*ddim[3]):trunc(0.95*ddim[3])
    indz2 <- trunc(0.1*ddim[3]):trunc(0.9*ddim[3])
    indz3 <- trunc(0.15*ddim[3]):trunc(0.85*ddim[3])
    z1 <- density(if(ls0ind>1) s0mean[indx1,indy1,indz1][s0mean[indx1,indy1,indz1]>0] else s0[indx1,indy1,indz1][s0[indx1,indy1,indz1]>0],bw=bw,n=1024)
    z2 <- density(if(ls0ind>1) s0mean[indx2,indy2,indz2][s0mean[indx2,indy2,indz2]>0] else s0[indx2,indy2,indz2][s0[indx2,indy2,indz2]>0],bw=bw,n=1024)
    z3 <- density(if(ls0ind>1) s0mean[indx3,indy3,indz3][s0mean[indx3,indy3,indz3]>0] else s0[indx3,indy3,indz3][s0[indx3,indy3,indz3]>0],bw=bw,n=1024)
    n <- prod(ddim)
    n1 <- length(indx1)*length(indy1)*length(indz1)
    n2 <- length(indx2)*length(indy2)*length(indz2)
    n3 <- length(indx3)*length(indy3)*length(indz3)
    ylim <- range(z$y,z1$y*n1/n,z2$y*n2/n,z3$y*n3/n)
    while(!accept){
      plot(z,type="l",main="Density of S0 values and cut off point",ylim=ylim)
      lines(z1$x,z1$y*n1/n,col=2)
      lines(z2$x,z2$y*n2/n,col=3)
      lines(z3$x,z3$y*n3/n,col=4)
      lines(c(A0,A0),c(0,max(z$y)/2),col=2,lwd=2)
      legend(min(A0,0.25*max(z$x)),ylim[2],c("Full cube",paste("Central",(n1*100)%/%n,"%"),
      paste("Central",(n2*100)%/%n,"%"),paste("Central",(n3*100)%/%n,"%")),col=1:4,lwd=rep(1,4))
      cat("A good cut off point should be left of support of the density of grayvalues within the head\n")
      a <- readline(paste("Accept current cut off point",A0," (Y/N):"))
      if (toupper(a) == "N") {
        cutpoint <-  readline("Provide value for cut off point:")
        cutpoint <- if(!is.null(cutpoint)) as.numeric(cutpoint) else A0
        if(!is.na(cutpoint)) {
          level <-A0 <- cutpoint
        }
      } else {
        accept <- TRUE
      }
    }
  } else {
    ddim <- object@ddim
    indx1 <- trunc(0.4*ddim[1]):trunc(0.6*ddim[1])
    indy1 <- trunc(0.4*ddim[2]):trunc(0.6*ddim[2])
    indz1 <- trunc(0.7*ddim[3]):trunc(0.7*ddim[3])
    A0a <- quantile(if(ls0ind>1) s0mean[indx1,indy1,indz1][s0mean[indx1,indy1,indz1]>1] else s0[indx1,indy1,indz1][s0[indx1,indy1,indz1]>1],.01)/(1+1/length(object@s0ind))
#  A0a provides a guess for a threshold based on lower quantiles of intensities
#  in a central cube (probably contained within the head)
#  the last factor adjusts for increased accuracy with replicated s0-values
    indx1 <- c(1:trunc(0.15*ddim[1]),trunc(0.85*ddim[1]):ddim[1])
    indy1 <- c(1:trunc(0.15*ddim[2]),trunc(0.85*ddim[2]):ddim[2])
    indz1 <- c(1:trunc(0.15*ddim[3]),trunc(0.85*ddim[3]):ddim[3])
    A0b <- quantile(if(ls0ind>1) s0mean[indx1,indy1,indz1] else s0[indx1,indy1,indz1],.99)
#  A0a provides a guess for a threshold based on upper quantiles of intensities
#  in cubes located at the edges (probably only containing noise
    level <- A0 <- min(A0a,A0b)*threshfactor
  }
  # determine parameters for linear relation between standard deviation and mean
  if(ls0ind>1) {
    s0sd <- apply(s0,1,sdmethod)
    ind <- s0mean>A0&s0mean<A1
    sdcoef <- coefficients(lm(s0sd[ind]~s0mean[ind]))
  } else {
    sdcoef <- awslinsd(s0,hmax=5,mask=NULL,A0=A0,A1=A1)$vcoef
  }
  object@level <- level
  object@sdcoef <- c(sdcoef,A0,A1)
  cat("Estimated parameters:",signif(sdcoef[1:2],3),"Interval of linearity",signif(A0,3),"-",signif(A1,3),"\n")
  object
})

############### [

setMethod("[","dtiData",function(x, i, j, k, drop=FALSE){
  args <- sys.call(-1)
  args <- c(x@call,args)
  if (missing(i)) i <- TRUE
  if (missing(j)) j <- TRUE
  if (missing(k)) k <- TRUE
  if (is.logical(i)) ddimi <- x@ddim[1] else ddimi <- length(i)
  if (is.logical(j)) ddimj <- x@ddim[2] else ddimj <- length(j)
  if (is.logical(k)) ddimk <- x@ddim[3] else ddimk <- length(k)

  invisible(new("dtiData",
                call   = args,
                si     = x@si[i,j,k,,drop=FALSE],
                gradient = x@gradient,
                btb    = x@btb,
                ngrad  = x@ngrad,
                s0ind  = x@s0ind,
                replind = x@replind,
                ddim   = c(ddimi,ddimj,ddimk),
                ddim0  = x@ddim0,
                xind   = x@xind[i],
                yind   = x@yind[j],
                zind   = x@zind[k],
                sdcoef = x@sdcoef,
                level  = x@level,
                voxelext = x@voxelext,
                orientation = x@orientation,
                source = x@source)
            )
})

##############

setMethod("[","dtiTensor",function(x, i, j, k, drop=FALSE){
  args <- sys.call(-1)
  args <- c(x@call,args)
  if (missing(i)) i <- TRUE
  if (missing(j)) j <- TRUE
  if (missing(k)) k <- TRUE
  if (is.logical(i)) ddimi <- x@ddim[1] else ddimi <- length(i)
  if (is.logical(j)) ddimj <- x@ddim[2] else ddimj <- length(j)
  if (is.logical(k)) ddimk <- x@ddim[3] else ddimk <- length(k)

  ind <- 1:prod(x@ddim)
  if(length(x@outlier)>0){
    ind <- rep(FALSE,prod(x@ddim))
    ind[x@outlier] <- TRUE
    dim(ind) <- x@ddim
    ind <- ind[i,j,k]
    outlier <- (1:length(ind))[ind]
  } else {
    outlier <- numeric(0)
  }

  invisible(new("dtiTensor",
                call  = args, 
                D     = x@D[,i,j,k,drop=FALSE],
                th0   = x@th0[i,j,k,drop=FALSE],
                sigma = if(x@method=="linear") x@sigma[i,j,k,drop=FALSE] else array(1,c(1,1,1)),
                scorr = x@scorr, 
                bw = x@bw,
                mask = x@mask[i,j,k,drop=FALSE],
                hmax = x@hmax,
                gradient = x@gradient,
                btb   = x@btb,
                ngrad = x@ngrad,
                s0ind = x@s0ind,
                replind = x@replind,
                ddim  = c(ddimi,ddimj,ddimk),
                ddim0 = x@ddim0,
                xind  = x@xind[i],
                yind  = x@yind[j],
                zind  = x@zind[k],
                voxelext = x@voxelext,
                level = x@level,
                orientation = x@orientation,
                outlier = outlier,
                scale = x@scale,
                source = x@source,
                method = x@method)
            )
})
#############
setMethod("[","dwiMixtensor",function(x, i, j, k, drop=FALSE){
  args <- sys.call(-1)
  args <- c(x@call,args)
  if (missing(i)) i <- TRUE
  if (missing(j)) j <- TRUE
  if (missing(k)) k <- TRUE
  if (is.logical(i)) ddimi <- x@ddim[1] else ddimi <- length(i)
  if (is.logical(j)) ddimj <- x@ddim[2] else ddimj <- length(j)
  if (is.logical(k)) ddimk <- x@ddim[3] else ddimk <- length(k)

  ind <- 1:prod(x@ddim)
  if(length(x@outlier)>0){
    ind <- rep(FALSE,prod(x@ddim))
    ind[x@outlier] <- TRUE
    dim(ind) <- x@ddim
    ind <- ind[i,j,k]
    outlier <- (1:length(ind))[ind]
  } else {
    outlier <- numeric(0)
  }

  invisible(new("dwiMixtensor",
                call  = args, 
                ev     = x@ev[,i,j,k,drop=FALSE],
                mix    = x@mix[,i,j,k,drop=FALSE],
                orient = x@orient[,,i,j,k,drop=FALSE],
                order  = x@order[i,j,k,drop=FALSE],
                p      = x@p,
                th0   = x@th0[i,j,k,drop=FALSE],
                sigma = x@sigma[i,j,k,drop=FALSE],
                scorr = x@scorr, 
                bw = x@bw,
                mask = x@mask[i,j,k,drop=FALSE],
                hmax = x@hmax,
                gradient = x@gradient,
                btb   = x@btb,
                ngrad = x@ngrad,
                s0ind = x@s0ind,
                replind = x@replind,
                ddim  = c(ddimi,ddimj,ddimk),
                ddim0 = x@ddim0,
                xind  = x@xind[i],
                yind  = x@yind[j],
                zind  = x@zind[k],
                voxelext = x@voxelext,
                level = x@level,
                orientation = x@orientation,
                outlier = outlier,
                scale = x@scale,
                source = x@source,
                method = x@method)
            )
})

#############

setMethod("[","dtiIndices",function(x, i, j, k, drop=FALSE){
  args <- sys.call(-1)
  args <- c(x@call,args)
  if (missing(i)) i <- TRUE
  if (missing(j)) j <- TRUE
  if (missing(k)) k <- TRUE
  if (is.logical(i)) ddimi <- x@ddim[1] else ddimi <- length(i)
  if (is.logical(j)) ddimj <- x@ddim[2] else ddimj <- length(j)
  if (is.logical(k)) ddimk <- x@ddim[3] else ddimk <- length(k)

  invisible(new("dtiIndices",
                call = args,
                fa = x@fa[i,j,k,drop=FALSE],
                ga = x@ga[i,j,k,drop=FALSE],
                md = x@md[i,j,k,drop=FALSE],
                andir = x@andir[,i,j,k,drop=FALSE],
                bary = x@bary[,i,j,k,drop=FALSE],
                gradient = x@gradient,
                btb   = x@btb,
                ngrad = x@ngrad,
                s0ind = x@s0ind,
                ddim  = c(ddimi,ddimj,ddimk),
                ddim0 = x@ddim0,
                voxelext = x@voxelext,
                orientation = x@orientation,
                xind  = x@xind[i],
                yind  = x@yind[j],
                zind  = x@zind[k],
                method = x@method,
                level = x@level,
                source= x@source)
            )
})

###########

setMethod("[","dwiQball",function(x, i, j, k, drop=FALSE){
  args <- sys.call(-1)
  args <- c(x@call,args)
  if (missing(i)) i <- TRUE
  if (missing(j)) j <- TRUE
  if (missing(k)) k <- TRUE
  if (is.logical(i)) ddimi <- x@ddim[1] else ddimi <- length(i)
  if (is.logical(j)) ddimj <- x@ddim[2] else ddimj <- length(j)
  if (is.logical(k)) ddimk <- x@ddim[3] else ddimk <- length(k)

  ind <- 1:prod(x@ddim)
  if(length(x@outlier)>0){
    ind <- rep(FALSE,prod(x@ddim))
    ind[x@outlier] <- TRUE
    dim(ind) <- x@ddim
    ind <- ind[i,j,k]
    outlier <- (1:length(ind))[ind]
  } else {
    outlier <- numeric(0)
  }

  invisible(new("dwiQball",
                call  = args, 
                order = x@order,
                lambda = x@lambda,
                sphcoef = x@sphcoef[,i,j,k,drop=FALSE],
                th0   = x@th0[i,j,k,drop=FALSE],
                sigma = x@sigma[i,j,k,drop=FALSE],
                scorr = x@scorr, 
                bw = x@bw,
                mask = x@mask[i,j,k,drop=FALSE],
                hmax = x@hmax,
                gradient = x@gradient,
                btb   = x@btb,
                ngrad = x@ngrad,
                s0ind = x@s0ind,
                replind = x@replind,
                ddim  = c(ddimi,ddimj,ddimk),
                ddim0 = x@ddim0,
                xind  = x@xind[i],
                yind  = x@yind[j],
                zind  = x@zind[k],
                voxelext = x@voxelext,
                level = x@level,
                orientation = x@orientation,
                outlier = outlier,
                scale = x@scale,
                source = x@source,
                what = x@what)
            )
})


########## extract()

extract <- function(x, ...) cat("Data extraction not defined for this class:",class(x),"\n")

setGeneric("extract", function(x, ...) standardGeneric("extract"))

setMethod("extract","dtiData",function(x, what="data", xind=TRUE, yind=TRUE, zind=TRUE){
  what <- tolower(what) 
  x <- x[xind,yind,zind]

  z <- list(NULL)
  if("gradient" %in% what) z$gradient <- x@gradient
  if("btb" %in% what) z$btb <- x@btb
  if("s0" %in% what) z$S0 <- x@si[,,,x@s0ind,drop=FALSE]
  if("sb" %in% what) z$Si <- x@si[,,,-x@s0ind,drop=FALSE]
  if("data" %in% what) z$data <- x@si
  invisible(z)
})

#############

setMethod("extract","dwiMixtensor",function(x, what="andir", xind=TRUE, yind=TRUE, zind=TRUE){
  what <- tolower(what) 
  x <- x[xind,yind,zind]
  n1 <- x@ddim[1]
  n2 <- x@ddim[2]
  n3 <- x@ddim[3]
  z <- list(NULL)
  if("order" %in% what) z$order <- x@order
  if("ev" %in% what) { 
     ev <- array(0,c(3,dim(x@ev)[-1]))
     ev[1,,,] <- x@ev[1,,,] + x@ev[2,,,]
     ev[2,,,] <- x@ev[2,,,]
     ev[3,,,] <- x@ev[2,,,]
     z$ev <- ev
     }
  if("mix" %in% what) z$mix <- x@mix
  if("andir" %in% what) {
     orient <- x@orient
     andir <- array(0,c(3,prod(dim(orient))/2))
     dim(orient) <- c(2,prod(dim(orient))/2)
     sth <- sin(orient[1,])
     andir[1,] <- sth*cos(orient[2,])
     andir[2,] <- sth*sin(orient[2,])
     andir[3,] <- cos(orient[1,])
     z$andir <- array(andir,c(3,dim(x@orient)[-1]))
     }
  if("s0" %in% what) z$S0 <- x@S0
  if("mask" %in% what) z$mask <- x@mask
  if("gfa" %in% what) z$gfa <- x@ev[1,,,]/sqrt((x@ev[1,,,]+x@ev[2,,,])^2+2*x@ev[2,,,]^2)
  invisible(z)
})

setMethod("extract","dtiTensor",function(x, what="tensor", xind=TRUE, yind=TRUE, zind=TRUE){
  what <- tolower(what) 

  x <- x[xind,yind,zind]
  n1 <- x@ddim[1]
  n2 <- x@ddim[2]
  n3 <- x@ddim[3]
  needev <- ("fa" %in% what) || ("ga" %in% what) || ("md" %in% what) || ("evalues" %in% what)
  needall <- needev && ("andir" %in% what)

  z <- list(NULL)
  if(needall){
    erg <- .Fortran("dti3Dall",
                    as.double(x@D),
                    as.integer(n1),
                    as.integer(n2),
                    as.integer(n3),
                    as.logical(x@mask),
                    fa=double(n1*n2*n3),
                    ga=double(n1*n2*n3),
                    md=double(n1*n2*n3),
                    andir=double(3*n1*n2*n3),
                    ev=double(3*n1*n2*n3),
                    DUPL=FALSE,
                    PACKAGE="dti")[c("fa","ga","md","andir","ev")]
    if("fa" %in% what) z$fa <- array(erg$fa,c(n1,n2,n3))
    if("ga" %in% what) z$ga <- array(erg$ga,c(n1,n2,n3))
    if("md" %in% what) z$md <- array(erg$md,c(n1,n2,n3))
    if("evalues" %in% what) z$evalues <- array(erg$ev,c(3,n1,n2,n3))
    if("andir" %in% what) z$andir <- array(erg$andir,c(3,n1,n2,n3))
  } else {
    if(needev){
      ev <- array(.Fortran("dti3Dev",
                           as.double(x@D),
                           as.integer(n1),
                           as.integer(n2),
                           as.integer(n3),
                           as.logical(x@mask),
                           ev=double(3*n1*n2*n3),
                           DUPL=FALSE,
                           PACKAGE="dti")$ev,c(3,n1,n2,n3))
      if("fa" %in% what) {
        dd <- apply(ev^2,2:4,sum)
        md <- (ev[1,,,]+ev[2,,,]+ev[3,,,])/3
        sev <- sweep(ev,2:4,md)
        z$fa <- sqrt(1.5*apply(sev^2,2:4,sum)/dd)
      }
      if("ga" %in% what) {
        sev <- log(ev)
        md <- (sev[1,,,]+sev[2,,,]+sev[3,,,])/3
        sev <- sweep(sev,2:4,md)
        ga <- sqrt(apply(sev^2,2:4,sum))
        ga[is.na(ga)] <- 0
        z$ga <- ga 
      }
      if("md" %in% what) z$md <- (ev[1,,,]+ev[2,,,]+ev[3,,,])/3
      if("evalues" %in% what) z$evalues <- ev
    }
    if("andir" %in% what){
      z$andir <- array(.Fortran("dti3Dand",
                                as.double(x@D),
                                as.integer(n1),
                                as.integer(n2),
                                as.integer(n3),
                                as.logical(x@mask),
                                andir=double(3*n1*n2*n3),
                                DUPL=FALSE,
                                PACKAGE="dti")$andir,c(3,n1,n2,n3))
    }
  }
  if("tensor" %in% what) z$tensor <- x@D
  if("s0" %in% what) z$s0 <- x@th0
  if("mask" %in% what) z$mask <- x@mask
  if("outlier" %in% what) {
    ind <- 1:prod(x@ddim)
    ind <- rep(FALSE,prod(x@ddim))
    if(length(x@outlier)>0) ind[x@outlier] <- TRUE
    dim(ind) <- x@ddim
  }
  invisible(z)
})

##############

setMethod("extract","dtiIndices",function(x, what=c("fa","andir"), xind=TRUE, yind=TRUE, zind=TRUE){
  what <- tolower(what) 

  x <- x[xind,yind,zind]
  n1 <- x@ddim[1]
  n2 <- x@ddim[2]
  n3 <- x@ddim[3]

  z <- list(NULL)
  if("fa" %in% what) z$fa <- x@fa
  if("ga" %in% what) z$ga <- x@ga
  if("md" %in% what) z$md <- x@md
  if("andir" %in% what) z$andir <- x@andir
  if("bary" %in% what) z$bary <- x@bary
  invisible(z)
})

##############

setMethod("extract","dwiQball",function(x, what="sphcoef", xind=TRUE, yind=TRUE, zind=TRUE){
  what <- tolower(what) 

  x <- x[xind,yind,zind]
  n1 <- x@ddim[1]
  n2 <- x@ddim[2]
  n3 <- x@ddim[3]

  z <- list(NULL)
  if("sphcoef" %in% what) z$sphcoef <- x@sphcoef
  if("s0" %in% what) z$s0 <- x@th0
  if("mask" %in% what) z$mask <- x@mask
  if("outlier" %in% what) {
    ind <- 1:prod(x@ddim)
    ind <- rep(FALSE,prod(x@ddim))
    if(length(x@outlier)>0) ind[x@outlier] <- TRUE
    dim(ind) <- x@ddim
  }
  invisible(z)
})

