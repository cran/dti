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
  if(prod(object@ddim)==1){
    warning("you need more than one voxel to model variances")
    return(object)
  }
  level0 <- if(is.null(level)) object@level else max(0,level)
  s0ind<-object@s0ind
  s0 <- object@si[,,,s0ind]
  ls0ind <- length(s0ind)
  A0 <- level0
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
    z <- density(if(ls0ind>1) s0mean[s0mean>0&s0mean<A1] else s0[s0>0&s0<A1],bw = max(bw,.01),,n=1024)
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
    if(is.null(level)){
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
  }
  # determine parameters for linear relation between standard deviation and mean
  if(ls0ind>1) {
    s0sd <- apply(s0,1,sdmethod)
    ind <- s0mean>A0&s0mean<A1
    if(length(ind)<2){
         warning("you need more than one voxel to model variances choice of A0/A1 to restrictive")
         return(object)
         }
    sdcoef <- coefficients(lm(s0sd[ind]~s0mean[ind]))
    if(sdcoef[1]<0){
       sdcoef <- numeric(2)
       sdcoef[1] <- .25  # this is an arbitrary (small) value to avaoid zero variances
       sdcoef[2] <- coefficients(lm(s0sd[ind]~s0mean[ind]-1))
       }
    if(sdcoef[2]<0){
       sdcoef <- numeric(2)
       sdcoef[1] <- max(0.25,mean(s0sd[ind]))
       sdcoef[2] <- 0
       }
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
  swap <- rep(FALSE,3)
  if (!is.logical(i)) swap[1] <- i[1] > i[length(i)]
  if (!is.logical(j)) swap[2] <- j[1] > j[length(j)]
  if (!is.logical(k)) swap[3] <- k[1] > k[length(k)]
  orientation <- x@orientation
  gradient <- x@gradient
  if(swap[1]) {
     orientation[1] <- (orientation[1]+1)%%2
     gradient[1,] <- -gradient[1,]
  }
  if(swap[2]) {
     orientation[2] <- (orientation[2]+1)%%2+2
     gradient[2,] <- -gradient[2,]
  }
  if(swap[3]) {
     orientation[3] <- (orientation[3]+1)%%2+4
     gradient[3,] <- -gradient[3,]
  }
  invisible(new("dtiData",
                call   = args,
                si     = x@si[i,j,k,,drop=FALSE],
                gradient = gradient,
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
                orientation = as.integer(orientation),
                rotation = x@rotation,
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
  swap <- rep(FALSE,3)
  if (!is.logical(i)) swap[1] <- i[1] > i[length(i)]
  if (!is.logical(j)) swap[2] <- j[1] > j[length(j)]
  if (!is.logical(k)) swap[3] <- k[1] > k[length(k)]
  orientation <- x@orientation
  gradient <- x@gradient
  btb <- x@btb
  D <- x@D
  if(swap[1]) {
     orientation[1] <- (orientation[1]+1)%%2
     gradient[1,] <- -gradient[1,]
     btb[2:3,] <- - btb[2:3,]
     D[2:3,,,] <- - D[2:3,,,]
  }
  if(swap[2]) {
     orientation[2] <- (orientation[2]+1)%%2+2
     gradient[2,] <- -gradient[2,]
     btb[c(2,5),] <- - btb[c(2,5),]
     D[c(2,5),,,] <- - D[c(2,5),,,]
  }
  if(swap[3]) {
     orientation[3] <- (orientation[3]+1)%%2+4
     gradient[3,] <- -gradient[3,]
     btb[c(3,5),] <- - btb[c(3,5),]
     D[c(3,5),,,] <- - D[c(3,5),,,]
  }
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
                D     = D[,i,j,k,drop=FALSE],
                th0   = x@th0[i,j,k,drop=FALSE],
                sigma = if(x@method=="linear") x@sigma[i,j,k,drop=FALSE] else array(1,c(1,1,1)),
                scorr = x@scorr, 
                bw = x@bw,
                mask = x@mask[i,j,k,drop=FALSE],
                hmax = x@hmax,
                gradient = gradient,
                btb   = btb,
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
                orientation = as.integer(orientation),
                rotation = x@rotation,
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
#  cat("indix i",i,"\n")
#  cat("indix j",j,"\n")
#  cat("indix k",k,"\n")
  swap <- rep(FALSE,3)
  if (!is.logical(i)) swap[1] <- i[1] > i[length(i)]
  if (!is.logical(j)) swap[2] <- j[1] > j[length(j)]
  if (!is.logical(k)) swap[3] <- k[1] > k[length(k)]
#  cat("swap",swap,"\n")
  orientation <- x@orientation
  gradient <- x@gradient
  btb <- x@btb
  orient <- x@orient
  if(swap[1]) {
     orientation[1] <- (orientation[1]+1)%%2
     gradient[1,] <- -gradient[1,]
     btb[2:3,] <- - btb[2:3,]
     orient <- - orient
  }
  if(swap[2]) {
     orientation[2] <- (orientation[2]+1)%%2+2
     gradient[2,] <- -gradient[2,]
     btb[c(2,5),] <- - btb[c(2,5),]
     orient[2,,,,] <- - orient[2,,,,]
 }
  if(swap[3]) {
     orientation[3] <- (orientation[3]+1)%%2+4
     gradient[3,] <- -gradient[3,]
     btb[c(3,5),] <- - btb[c(3,5),]
     orient[1,,,,] <- pi - orient[1,,,,]
  }
#  cat("new orientation",orientation,"\n")
#  cat("indix i",i,"\n")
#  cat("indix j",j,"\n")
#  cat("indix k",k,"\n")
  invisible(new("dwiMixtensor",
                call  = args, 
                ev     = x@ev[,i,j,k,drop=FALSE],
                mix    = x@mix[,i,j,k,drop=FALSE],
                orient = orient[,,i,j,k,drop=FALSE],
                order  = x@order[i,j,k,drop=FALSE],
                p      = x@p,
                th0   = x@th0[i,j,k,drop=FALSE],
                sigma = x@sigma[i,j,k,drop=FALSE],
                scorr = x@scorr, 
                bw = x@bw,
                mask = x@mask[i,j,k,drop=FALSE],
                hmax = x@hmax,
                gradient = gradient,
                btb   = btb,
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
                orientation = as.integer(orientation),
                rotation = x@rotation,
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
  swap <- rep(FALSE,3)
  if (!is.logical(i)) swap[1] <- i[1] > i[length(i)]
  if (!is.logical(j)) swap[2] <- j[1] > j[length(j)]
  if (!is.logical(k)) swap[3] <- k[1] > k[length(k)]
  orientation <- x@orientation
  gradient <- x@gradient
  btb <- x@btb
  andir <- x@andir
  if(swap[1]) {
     orientation[1] <- (orientation[1]+1)%%2
     gradient[1,] <- -gradient[1,]
     btb[2:3,] <- - btb[2:3,]
     andir[1,,,] <- - andir[1,,,]
  }
  if(swap[2]) {
     orientation[2] <- (orientation[2]+1)%%2+2
     gradient[2,] <- -gradient[2,]
     btb[c(2,5),] <- - btb[c(2,5),]
     andir[2,,,] <- - andir[2,,,]
  }
  if(swap[3]) {
     orientation[3] <- (orientation[3]+1)%%2+4
     gradient[3,] <- -gradient[3,]
     btb[c(3,5),] <- - btb[c(3,5),]
     andir[3,,,] <- - andir[3,,,]
  }

  invisible(new("dtiIndices",
                call = args,
                fa = x@fa[i,j,k,drop=FALSE],
                ga = x@ga[i,j,k,drop=FALSE],
                md = x@md[i,j,k,drop=FALSE],
                andir = andir[,i,j,k,drop=FALSE],
                bary = x@bary[,i,j,k,drop=FALSE],
                gradient = gradient,
                btb   = btb,
                ngrad = x@ngrad,
                s0ind = x@s0ind,
                ddim  = c(ddimi,ddimj,ddimk),
                ddim0 = x@ddim0,
                voxelext = x@voxelext,
                orientation = as.integer(orientation),
                rotation = x@rotation,
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
    swap <- rep(FALSE,3)
  if (!is.logical(i)) swap[1] <- i[1] > i[length(i)]
  if (!is.logical(j)) swap[2] <- j[1] > j[length(j)]
  if (!is.logical(k)) swap[3] <- k[1] > k[length(k)]
  if(any(swap)) {
     warning("can't reverse order of indices")
     return(invisible(x))
  }
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
                rotation = x@rotation,
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
  swap <- rep(FALSE,3)
  if(is.numeric(xind)) swap[1] <- xind[1]>xind[length(xind)]
  if(is.numeric(yind)) swap[2] <- yind[1]>yind[length(yind)]
  if(is.numeric(zind)) swap[3] <- zind[1]>zind[length(zind)]
  if(any(swap)) {
     warning("can't reverse order of indices ")
     return(NULL)
  }
  x <- x[xind,yind,zind]
  z <- list(NULL)
  if("gradient" %in% what) z$gradient <- x@gradient
  if("btb" %in% what) z$btb <- x@btb
  if("s0" %in% what) z$s0 <- x@si[,,,x@s0ind,drop=FALSE]
  if("sb" %in% what) z$sb <- x@si[,,,-x@s0ind,drop=FALSE]
  if("siq" %in% what) {
     S0 <- x@si[,,,x@s0ind,drop=FALSE]
     Si <- x@si[,,,-x@s0ind,drop=FALSE]
     z$siq <- sweep(Si,1:3,apply(S0,1:3,mean),"/")
     z$siq[is.na(z$siq)] <- 0
  }
  if("data" %in% what) z$data <- x@si
  invisible(z)
})

#############

setMethod("extract","dwiMixtensor",function(x, what="andir", xind=TRUE, yind=TRUE, zind=TRUE){
  what <- tolower(what) 
  swap <- rep(FALSE,3)
  if(is.numeric(xind)) swap[1] <- xind[1]>xind[length(xind)]
  if(is.numeric(yind)) swap[2] <- yind[1]>yind[length(yind)]
  if(is.numeric(zind)) swap[3] <- zind[1]>zind[length(zind)]
  if(any(swap)){
     warning("can't reverse order of indices ")
     return(NULL)
  }
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
  if("s0" %in% what) z$s0 <- x@S0
  if("mask" %in% what) z$mask <- x@mask
  if("fa" %in% what){
      fa <- x@ev[1,,,]/sqrt((x@ev[1,,,]+x@ev[2,,,])^2+2*x@ev[2,,,]^2)
      fa[x@order==0] <- 0
      dim(fa) <- x@ddim
      z$fa <- fa
    }
  if("eorder" %in% what) {
     maxorder <- dim(x@mix)[1]
     mix <- x@mix
     dim(mix) <- c(maxorder,n1*n2*n3)     
     z$eorder <- array((2*(1:maxorder)-1)%*%mix,x@ddim)
  }
  if("bic" %in% what) {
      ngrad <- x@ngrad      
      ns0 <- length(x@s0ind)
      iso <- apply(x@mix,-1,sum)
      iso <- iso>0&&iso<1e0-1e-8
      penBIC <- log(ngrad-ns0)/(ngrad-ns0)*(iso+1+2*x@order)
      z$bic <- array(log(pmax(1e-10,x@sigma))+penBIC,dim(x@sigma))
  }
  if("aic" %in% what) {
      ngrad <- x@ngrad      
      ns0 <- length(x@s0ind)
      iso <- apply(x@mix,-1,sum)
      iso <- iso>0&&iso<1e0-1e-8
      penAIC <- 2/(ngrad-ns0)*(iso+1+2*x@order)
      z$aic <- array(log(pmax(1e-10,x@sigma))+penAIC,dim(x@sigma))
  }
  invisible(z)
})



setMethod("extract","dtiTensor",function(x, what="tensor", xind=TRUE, yind=TRUE, zind=TRUE){
  what <- tolower(what) 
  swap <- rep(FALSE,3)
  if(is.numeric(xind)) swap[1] <- xind[1]>xind[length(xind)]
  if(is.numeric(yind)) swap[2] <- yind[1]>yind[length(yind)]
  if(is.numeric(zind)) swap[3] <- zind[1]>zind[length(zind)]
  if(any(swap)) {
     warning("can't reverse order of indices ")
     return(NULL)
  }

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
                    DUP=FALSE,
                    PACKAGE="dti")[c("fa","ga","md","andir","ev")]
    if("fa" %in% what) z$fa <- array(erg$fa,x@ddim)
    if("ga" %in% what) z$ga <- array(erg$ga,x@ddim)
    if("md" %in% what) z$md <- array(erg$md,x@ddim)
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
                           DUP=FALSE,
                           PACKAGE="dti")$ev,c(3,n1,n2,n3))
      if("fa" %in% what) {
        dd <- apply(ev^2,2:4,sum)
        md <- (ev[1,,,]+ev[2,,,]+ev[3,,,])/3
        sev <- sweep(ev,2:4,md)
        z$fa <- array(sqrt(1.5*apply(sev^2,2:4,sum)/dd),x@ddim)
      }
      if("ga" %in% what) {
        sev <- log(ev)
        md <- (sev[1,,,]+sev[2,,,]+sev[3,,,])/3
        sev <- sweep(sev,2:4,md)
        ga <- sqrt(apply(sev^2,2:4,sum))
        ga[is.na(ga)] <- 0
        z$ga <- array(ga,x@ddim)
      }
      if("md" %in% what) z$md <- array((ev[1,,,]+ev[2,,,]+ev[3,,,])/3,x@ddim)
      if("evalues" %in% what) z$evalues <- array(ev,c(3,x@ddim))
    }
    if("andir" %in% what){
      z$andir <- array(.Fortran("dti3Dand",
                                as.double(x@D),
                                as.integer(n1),
                                as.integer(n2),
                                as.integer(n3),
                                as.logical(x@mask),
                                andir=double(3*n1*n2*n3),
                                DUP=FALSE,
                                PACKAGE="dti")$andir,c(3,n1,n2,n3))
    }
  }
  if("tensor" %in% what) z$tensor <- array(x@D,c(6,x@ddim))
  if("s0" %in% what) z$s0 <- array(x@th0,x@ddim)
  if("mask" %in% what) z$mask <- x@mask
  if("bic" %in% what) {
      ngrad <- x@ngrad      
      ns0 <- length(x@s0ind)
      penBIC <- log(ngrad-ns0)/(ngrad-ns0)*6
      z$bic <- array(log(pmax(1e-10,x@sigma))+penBIC,dim(x@sigma))
  }
  if("aic" %in% what) {
      ngrad <- x@ngrad      
      ns0 <- length(x@s0ind)
      penAIC <- 12/(ngrad-ns0)
      z$aic <- array(log(pmax(1e-10,x@sigma))+penAIC,dim(x@sigma))
  }
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
  swap <- rep(FALSE,3)
  if(is.numeric(xind)) swap[1] <- xind[1]>xind[length(xind)]
  if(is.numeric(yind)) swap[2] <- yind[1]>yind[length(yind)]
  if(is.numeric(zind)) swap[3] <- zind[1]>zind[length(zind)]
  if(any(swap)) {
     warning("can't reverse order of indices ")
     return(NULL)
  }

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
  swap <- rep(FALSE,3)
  if(is.numeric(xind)) swap[1] <- xind[1]>xind[length(xind)]
  if(is.numeric(yind)) swap[2] <- yind[1]>yind[length(yind)]
  if(is.numeric(zind)) swap[3] <- zind[1]>zind[length(zind)]
  if(any(swap)) {
     warning("can't reverse order of indices ")
     return(NULL)
  }

  x <- x[xind,yind,zind]
  n1 <- x@ddim[1]
  n2 <- x@ddim[2]
  n3 <- x@ddim[3]

  z <- list(NULL)
  if("sphcoef" %in% what) z$sphcoef <- x@sphcoef
  if("s0" %in% what) z$s0 <- x@th0
  if("mask" %in% what) z$mask <- x@mask
  if("bic" %in% what) {
      ngrad <- x@ngrad      
      ns0 <- length(x@s0ind)
      ord <- x@order
      penBIC <- log(ngrad-ns0)/(ngrad-ns0)*(ord+1)*(ord+2)/2
      z$bic <- array(log(pmax(1e-10,x@sigma))+penBIC,dim(x@sigma))
  }
  if("aic" %in% what) {
      ngrad <- x@ngrad      
      ns0 <- length(x@s0ind)
      ord <- x@order
      penAIC <- (ord+1)*(ord+2)/(ngrad-ns0)
      z$aic <- array(log(pmax(1e-10,x@sigma))+penAIC,dim(x@sigma))
  }
  if("outlier" %in% what) {
    ind <- 1:prod(x@ddim)
    ind <- rep(FALSE,prod(x@ddim))
    if(length(x@outlier)>0) ind[x@outlier] <- TRUE
    dim(ind) <- x@ddim
  }
  invisible(z)
})

