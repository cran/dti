################################################################
#                                                              #
# Section for HARDI functions                                  #
#                                                              #
################################################################

dwiQball <- function(object,  ...) cat("No DWI Q-ball calculation defined for this class:",class(object),"\n")

setGeneric("dwiQball", function(object,  ...) standardGeneric("dwiQball"))

setMethod("dwiQball","dtiData",function(object,what="wODF",order=4,lambda=0){
  args <- sys.call(-1)
  args <- c(object@call,args)
  if (!(what %in% c("ODF","wODF","aODF","ADC"))) {
      stop("what should specify either ODF, wODF, aODF, or ADC\n")
          }
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  sdcoef <- object@sdcoef
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

  # prepare data including mask
  ngrad0 <- ngrad - length(s0ind)
  s0 <- si[,,,s0ind]
  si <- si[,,,-s0ind]
  if (ns0>1) {
    dim(s0) <- c(prod(ddim),ns0)
    s0 <- s0 %*% rep(1/ns0,ns0)
    dim(s0) <- ddim
  }
  mask <- s0 > object@level
  mask <- connect.mask(mask)

  lord <- rep(seq(0,order,2),2*seq(0,order,2)+1)
  while(length(lord)>=ngrad0){
     order <- order-2
     lord <- rep(seq(0,order,2),2*seq(0,order,2)+1)
     cat("Reduced order of spherical harmonics to",order,"\n")
  }
  cat("Using",length(lord),"sperical harmonics\n")
  L <- -diag(lord*(lord+1))
  # now switch for different cases
  if (what=="ODF") {
     cat("Data transformation started ",format(Sys.time()),"\n")
     dim(s0) <- dim(si) <- NULL
     si[is.na(si)] <- 0
     si[(si == Inf)] <- 0
     si[(si == -Inf)] <- 0
     dim(si) <- c(prod(ddim),ngrad0)
     si <- t(si)
     cat("Data transformation completed ",format(Sys.time()),"\n")

     # get SH design matrix ...
     z <- design.spheven(order,object@gradient[,-s0ind],lambda)
     # ... and estimate coefficients of SH series of ODF
     # see Descoteaux et al. (2007)
     # include FRT(SH) -> P_l(0)
     sicoef <- z$matrix%*% si
     sphcoef <- plzero(order)%*%sicoef
     cat("Estimated coefficients for ODF (order=",order,") ",format(Sys.time()),"\n")
  } else if (what=="wODF") {
     cat("Data transformation started ",format(Sys.time()),"\n")
     dim(s0) <- dim(si) <- NULL
     si <- si/s0
     si[is.na(si)] <- 0
# Regularization following Aganj et al. (2010) delta=1e-3
     ind1 <- si<0
     ind2 <- (si<1e-3)&(!ind1)
     ind3 <- si>1 
     ind4 <- (si>1-1e-3)&(!ind3)
     si[ind1] <- 5e-4
     si[ind2] <- 5e-4+5e2*si[ind2]^2
     si[ind3] <- 1 - 5e-4
     si[ind4] <- 1 - 5e-4 - 5e2*(1-si[ind4])^2
#     si[si>=1] <- 1-.Machine$double.neg.eps
     si <- log( -log(si))
     si[is.na(si)] <- 0
     si[(si == Inf)] <- 0
     si[(si == -Inf)] <- 0
     dim(si) <- c(prod(ddim),ngrad0)
     si <- t(si)
     cat("Data transformation completed ",format(Sys.time()),"\n")

     # get SH design matrix ...
     z <- design.spheven(order,object@gradient[,-s0ind],lambda)
     # ... and estimate coefficients of SH series of ODF
     # see Aganj et al. (2009)
     # include FRT(SH) -> P_l(0)
     sicoef <- z$matrix%*% si
     plz <- plzero(order)/2/pi
     sphcoef <- plz%*%L%*%sicoef
     coef0 <- sphcoef[1,]
     sphcoef[1,] <- 1/2/sqrt(pi)
     sphcoef[-1,] <- sphcoef[-1,]/8/pi
     cat("Estimated coefficients for wODF (order=",order,") ",format(Sys.time()),"\n")
  } else if (what=="aODF") {
     cat("Data transformation started ",format(Sys.time()),"\n")
     dim(s0) <- dim(si) <- NULL
     si <- si/s0
     si[is.na(si)] <- 0
     si[si>=1] <- 1-.Machine$double.neg.eps
     si <- 1/(-log(si))
     si[is.na(si)] <- 0
     si[(si == Inf)] <- 0
     si[(si == -Inf)] <- 0
     dim(si) <- c(prod(ddim),ngrad0)
     si <- t(si)
     cat("Data transformation completed ",format(Sys.time()),"\n")

     # get SH design matrix ...
     z <- design.spheven(order,object@gradient[,-s0ind],lambda)
     # ... and estimate coefficients of SH series of ODF
     # see Descoteaux et al. (2007)
     # include FRT(SH) -> P_l(0)
     sicoef <- z$matrix%*% si
     sphcoef <- plzero(order)%*%sicoef
     cat("Estimated coefficients for aODF (order=",order,") ",format(Sys.time()),"\n")
  } else { #  what == "ADC" 
     cat("Data transformation started ",format(Sys.time()),"\n")
     dim(s0) <- dim(si) <- NULL
     si <- si/s0
     si[is.na(si)] <- 0
     si[si>=1] <- 1-.Machine$double.neg.eps
     si <- -log(si)
#     si <- -log(si)
     si[is.na(si)] <- 0
     si[(si == Inf)] <- 0
     si[(si == -Inf)] <- 0
     dim(si) <- c(prod(ddim),ngrad0)
     si <- t(si)
     cat("Data transformation completed ",format(Sys.time()),"\n")

     # get SH design matrix ...
     z <- design.spheven(order,object@gradient[,-s0ind],lambda)
     # ... and estimate coefficients of SH series of ADC
     sphcoef <- sicoef <- z$matrix%*% si
     cat("Estimated coefficients for ADC expansion in spherical harmonics (order=",order,") ",format(Sys.time()),"\n")
  }
  res <- si - t(z$design) %*% sicoef
  rss <- res[1,]^2
  for(i in 2:ngrad0) rss <- rss + res[i,]^2
  sigma2 <- rss/(ngrad0-length(lord))
  if(what %in% c("ODF","aODF")){
     varcoef <- outer(diag(plzero(order))^2*diag(z$matrix%*%t(z$matrix)),sigma2,"*")
  } else if(what=="wODF"){
     varcoef <- outer(diag(plzero(order))^2*diag(L)^2*diag(z$matrix%*%t(z$matrix)),sigma2,"*")
     varcoef[-1,] <- varcoef[-1,]/256/pi^4
     varcoef[1,] <- 0
  } else {
     varcoef <- outer(diag(z$matrix%*%t(z$matrix)),sigma2,"*")
  }
  dim(sigma2) <- ddim
  sphcoef[,!mask] <- 0
  dim(sphcoef) <- dim(varcoef) <- c((order+1)*(order+2)/2,ddim)
  dim(res) <- c(ngrad0,ddim)
  cat("Variance estimates generated ",format(Sys.time()),"\n")
  th0 <- array(s0,object@ddim)
  th0[!mask] <- 0
  gc()
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
  invisible(new("dwiQball",
                call  = args,
                order = as.integer(order),
                lambda = lambda,
                sphcoef = sphcoef,
                varsphcoef = varcoef,
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
                scale = 0.5,
                what = what)
            )
})

design.spheven <- function(order,gradients,lambda){
#
#  compute design matrix for Q-ball
#
  order <- as.integer(max(0,order))
  if(order%%2==1){
    warning("maximum order needs to be even, increase order by one")
    order <- order+1
  } 

  # calculate spherical angles theta and phi corresponding to the gradients
  n <- dim(gradients)[2]
  theta <- phi <- numeric(n)
  for( i in 1:n){
    angles <- sphcoord(gradients[,i])
    theta[i] <- angles[1]
    phi[i] <-  angles[2]
  }

  # values of SH on specified spherical angles
  sphharmonics <- getsphericalharmonicseven(order,theta,phi)
  # Laplace-Beltrami-Regularization term
  lord <- rep(seq(0,order,2),2*seq(0,order,2)+1)
  L <- lambda*diag(lord^2*(lord+1)^2)
  # transformation matrix for SH coefficients
  ttt <- solve(sphharmonics%*%t(sphharmonics)+L)%*%sphharmonics
  # results
  list(design = sphharmonics,
       matrix = ttt,
       theta = theta,
       phi = phi)
}

plzero <- function(order){
  l <- seq(2,order,2)
  pl <- l
  for(i in 1:length(l)) pl[i] <- (-1)^(l[i]/2)*prod(seq(1,(l[i]-1),2))/prod(seq(2,l[i],2))
  2*pi*diag(rep(c(1,pl),2*seq(0,order,2)+1))
}

