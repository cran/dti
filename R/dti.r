#
#
#

setMethod("print", "dti",
function(x){
    cat("  DTI object of class", class(x),"\n")
    cat("  Dimension            :", paste(x@ddim, collapse="x"), "\n")
    cat("  Number of Gradients  :", paste(x@ngrad, collapse="x"), "\n")
    cat("  Filename             :", x@source, "\n")
    cat("  Slots                :\n")
    print(slotNames(x))
    invisible(NULL)
})
setMethod("summary", "dti",
function(object){
    cat("  DTI object of class", class(object),"\n")
    cat("  Filename             :", object@source, "\n")
    cat("  Dimension            :", paste(object@ddim, collapse="x"), "\n")
    cat("  Number of Gradients  :", paste(object@ngrad, collapse="x"), "\n")
    cat("  Voxel extensions     :", paste(object@voxelext, collapse="x"), "\n")
    cat("  Ind. S0-Images       :", paste(object@s0ind, collapse="x"), "\n")
    if(class(object)=="dtiTensor"){
    cat("  Voxel in mask        :", paste(sum(object@mask), collapse="x"), "\n")
    cat("  Spatial smoothness   :", paste(signif(object@bw,3), collapse="x"), "\n")
    cat("  mean variance        :", paste(signif(mean(object@sigma[object@mask]),3), collapse="x"), "\n")
    cat("  hmax                 :", paste(object@hmax, collapse="x"), "\n")
}
    cat("\n")
    invisible(NULL)
})

setMethod("plot", "dtiTensor", function(x, y, slice=1, view="axial", quant=0, minanindex=NULL, contrast.enh=1, qrange=c(.01,.99),xind=NULL,yind=NULL,zind=NULL, mar=c(2,2,2,.2),mgp=c(2,1,0),...) {
  if(is.null(x@D)) cat("No diffusion tensor yet")
  adimpro <- require(adimpro)
  if(is.null(xind)) xind<-(1:x@ddim[1])
  if(is.null(yind)) yind<-(1:x@ddim[2])
  if(is.null(zind)) zind<-(1:x@ddim[3])
  if (view == "sagittal") {
    D <- x@D[,slice,yind,zind]
    mask <- x@mask[slice,yind,zind]
    n1 <- length(yind)
    n2 <- length(zind)
  } else if (view == "coronal") {
    D <- x@D[,xind,slice,zind]
    mask <- x@mask[xind,slice,zind]
    n1 <- length(xind)
    n2 <- length(zind)
  } else {
    D <- x@D[,xind,yind,slice]
    mask <- x@mask[xind,yind,slice]
    n1 <- length(xind)
    n2 <- length(yind)
  }
  z <- .Fortran("dtiind2D",
                as.double(D),
                as.integer(n1),
                as.integer(n2),
                as.logical(mask),
                fa=double(n1*n2),
                md=double(n1*n2),
                andir=double(3*n1*n2),
                DUPL=FALSE,
                PACKAGE="dti")[c("fa","md","andir")]
   oldpar <- par(mfrow=c(3,3),mar=mar,mgp=mgp,...)
#  now draw information to graphical device
   on.exit(par(oldpar))
   img<-D[1,,]
   rg<-quantile(img,qrange)
   img[img>rg[2]]<-rg[2]
   show.image(make.image(65535*img/max(img)))
   title(paste("Dxx: mean",signif(mean(D[mask]),3),"max",signif(max(D[1,,][mask]),3)))
   img<-D[2,,]
   rg<-quantile(img,qrange)
   img[img>rg[2]]<-rg[2]
   img[img<rg[1]]<-rg[1]
   show.image(make.image(img))
   title(paste("Dxy: min",signif(min(D[2,,][mask]),3),"max",signif(max(D[2,,][mask]),3)))
   img<-D[3,,]
   rg<-quantile(img,qrange)
   img[img>rg[2]]<-rg[2]
   img[img<rg[1]]<-rg[1]
   show.image(make.image(img))
   title(paste("Dxz: min",signif(min(D[3,,][mask]),3),"max",signif(max(D[3,,][mask]),3)))
   show.image(make.image(matrix(z$fa,n1,n2)))
   title(paste("Anisotropy index (FA)  range:",signif(min(z$fa[mask]),3),"-",
                signif(max(z$fa[mask]),3)))
   img<-D[4,,]
   rg<-quantile(img,qrange)
   img[img>rg[2]]<-rg[2]
   img[img<rg[1]]<-rg[1]
   show.image(make.image(65535*img/max(img)))
   title(paste("Dyy: min",signif(min(D[4,,][mask]),3),"max",signif(max(D[4,,][mask]),3)))
   img<-D[5,,]
   rg<-quantile(img,qrange)
   img[img>rg[2]]<-rg[2]
   img[img<rg[1]]<-rg[1]
   show.image(make.image(img))
   title(paste("Dyz: min",signif(min(D[5,,][mask]),3),"max",signif(max(D[5,,][mask]),3)))
   andir.image(matrix(z$fa,n1,n2),array(z$andir,c(3,n1,n2)),quant=quant,minanindex=minanindex)
   title(paste("Anisotropy directions"))
   img <- matrix(z$md,n1,n2)
   show.image(make.image(65535*img/max(img)))
   title(paste("Mean diffusivity   range:",signif(min(z$md[mask]),3),"-",
                signif(max(z$md[mask]),3)))
   img<-D[6,,]
   rg<-quantile(img,qrange)
   img[img>rg[2]]<-rg[2]
   img[img<rg[1]]<-rg[1]
   show.image(make.image(65535*img/max(img)))
   title(paste("Dzz: min",signif(min(D[6,,][mask]),3),"max",signif(max(D[6,,][mask]),3)))
   invisible(NULL)
})

setMethod("plot", "dtiData", function(x, y,slice=1, gradient=NULL, view= "axial", show=TRUE,xind=NULL,yind=NULL,zind=NULL, mar=c(3,3,3,.3),mgp=c(2,1,0), ...) {
if(is.null(x@si)) cat("No dwi data yet")
maxsi <- max(x@si)
  if(is.null(xind)) xind<-(1:x@ddim[1])
  if(is.null(yind)) yind<-(1:x@ddim[2])
  if(is.null(zind)) zind<-(1:x@ddim[3])
if(is.null(gradient)) gradient <- x@s0ind[1]
if(gradient<1||gradient>x@ngrad) {
   warning("gradient number out of range, show s0 image")
   gradient <- x@s0ind[1]
}
adimpro <- require(adimpro)
if (view == "sagittal") {
   if(slice<1||slice>x@ddim[1]) {
      warning("slice number out of range, show central slice")
      slice <- x@ddim[1]%/%2
   }
   img <- x@si[slice,yind,zind,gradient]
  } else if (view == "coronal") {
   if(slice<1||slice>x@ddim[2]) {
      warning("slice number out of range, show central slice")
      slice <- x@ddim[2]%/%2
   }
   img <- x@si[xind,slice,zind,gradient]
  } else {
   if(slice<1||slice>x@ddim[3]) {
      warning("slice number out of range, show central slice")
      slice <- x@ddim[3]%/%2
   }
   img <- x@si[xind,yind,slice,gradient]
  }
  oldpar <- par(mar=mar,mgp=mgp, ...)
  if(adimpro) {
     img <- make.image(65535*img/maxsi)
     if(show) show.image(img,...)
    } else if(show) {
      image(img,...)
    }
    par(oldpar)
    invisible(img)
}
)
setMethod("plot", "dti", function(x, y, ...) cat("No implementation for class dti\n"))

setMethod("plot", "dtiIndices", 
function(x, y, slice=1, view= "axial", method=1, quant=0, minanindex=NULL, show=TRUE, contrast.enh=1,xind=NULL,yind=NULL,zind=NULL, mar=c(3,3,3,.3),mgp=c(2,1,0), ...) {
  if(is.null(x@fa)) cat("No anisotropy index yet")
  if(!(method %in% 1:3)) {
      warning("method out of range, reset to 1")
      method <- 1
  }
  if(is.null(xind)) xind<-(1:x@ddim[1])
  if(is.null(yind)) yind<-(1:x@ddim[2])
  if(is.null(zind)) zind<-(1:x@ddim[3])
  adimpro <- require(adimpro)
  oldpar <- par(mar=mar,mgp=mgp, ...)
  if (view == "sagittal") {
    anindex <- x@fa[slice,yind,zind]
    andirection <- x@andir[,slice,yind,zind]
  } else if (view == "coronal") {
    anindex <- x@fa[xind,slice,zind]
    andirection <- x@andir[,xind,slice,zind]
  } else {
    anindex <- x@fa[xind,yind,slice]
    andirection <- x@andir[,xind,yind,slice]
  }
    anindex[anindex>1]<-0
    anindex[anindex<0]<-0
  if ((method==1) || (method==2)) {
    if(contrast.enh<1&&fa.contrast.enh>0) anindex <- pmin(anindex/contrast.enh,1)
    if(is.null(minanindex)) minanindex <- quantile(anindex,quant,na.rm=TRUE)
    if (diff(range(anindex,na.rm=TRUE)) == 0) minanindex <- 0
    if(method==1) {
      andirection[1,,] <- abs(andirection[1,,])
      andirection[2,,] <- abs(andirection[2,,])
      andirection[3,,] <- abs(andirection[3,,])
    } else {
      ind<-andirection[1,,]<0
      dim(andirection) <- c(3,prod(dim(ind)))
      andirection[,ind] <- - andirection[,ind]
      andirection[2,] <- (1+andirection[2,])/2
      andirection[3,] <- (1+andirection[3,])/2
      dim(andirection) <- c(3,dim(ind))
    }
    andirection <- aperm(andirection,c(2,3,1))
    andirection <- andirection*as.vector(anindex)*as.numeric(anindex>minanindex)
    if(adimpro) {
      andirection[is.na(andirection)] <- 0
      andirection <- make.image(andirection,gamma=TRUE)
      if(show) show.image(andirection,...)
    } else if(show) {
      dim(anindex) <- dim(andirection)[2:3]
      image(anindex,...)
    }
    par(oldpar)
    invisible(andirection)
  } else if (method==3) {
    if(adimpro) {
      bary[is.na(bary)] <- 0
      bary <- make.image(aperm(bary,c(2,3,1)))
      if(show) show.image(bary,...)
    } else if(show) {
      image(bary[1,,],...)
    }
    par(oldpar)
    invisible(bary)
  } 
})

#
#
#

dtiData <- function(gradient,imagefile,ddim,xind=NULL,yind=NULL,zind=NULL,level=0,mins0value=0,maxvalue=10000,voxelext=c(1,1,1),orientation=c(1,3,5)) {
  if (any(sort((orientation)%/%2) != 0:2)) stop("invalid orientation \n")
  if (dim(gradient)[2]==3) gradient <- t(gradient)
  if (dim(gradient)[1]!=3) stop("Not a valid gradient matrix")
  ngrad <- dim(gradient)[2]
  s0ind <- (1:ngrad)[apply(abs(gradient),2,max)==0] 
  if (!(file.exists(imagefile))) stop("Image file does not exist")
  zz <- file(imagefile,"rb")
#  si now contains all images (S_0 and S_I), ngrad includes 
#  number of zero gradients
  si <- readBin(zz,"integer",prod(ddim)*ngrad,2,FALSE)
  close(zz)
  cat("Data successfully read",date(),proc.time(), "\n")

  if (is.null(xind)) xind <- 1:ddim[1]
  if (is.null(yind)) yind <- 1:ddim[2]
  if (is.null(zind)) zind <- 1:ddim[3]
  dim(si) <- c(ddim,ngrad)
  si <- si[xind,yind,zind,] 
  dimsi <- dim(si)
#
#   set correct orientation
#
  xyz <- (orientation)%/%2+1
  swap <- orientation-2*(orientation%/%2)
  if(any(xyz!=1:3)) {
      abc <- 1:3
      abc[xyz] <- abc
      si <- aperm(si,c(abc,4))
      swap[xyz] <- swap
      voxelext[xyz] <- voxelext
      dimsi[xyz] <- dimsi[1:3]
      ddim[xyz] <- ddim[1:3]
      gradient[xyz,] <- gradient
  }
  if(swap[1]==1) {
      si <- si[dimsi[1]:1,,,] 
      gradient[1,] <- -gradient[1,]
      }
  if(swap[2]==1) {
      si <- si[,dimsi[2]:1,,]  
      gradient[2,] <- -gradient[2,]
      }
  if(swap[3]==0) {
      si <- si[,,dimsi[3]:1,]    
      gradient[3,] <- -gradient[3,]
      }
#
#   orientation set to radiological convention
#
  si <- .Fortran("initdata",
                 si=as.integer(si),
                 as.integer(dimsi[1]),
                 as.integer(dimsi[2]),
                 as.integer(dimsi[3]),
                 as.integer(dimsi[4]),
                 as.integer(maxvalue),
                 PACKAGE="dti")$si
#  this replaces the content off all voxel with elements <=0 or >maxvalue by 0
     dim(si) <- dimsi
  level <- max(mins0value,level*mean(si[,,,s0ind][si[,,,s0ind]>0])) # set level to level*mean  of positive s_0 values
  ddim0 <- as.integer(ddim)
  ddim <- as.integer(dim(si)[1:3])

  cat("Create auxiliary statistics",date(),proc.time(), " \n")
  btb <- create.designmatrix.dti(gradient)
  rind <- replind(gradient)
  
  invisible(new("dtiData",
                si     = si,
                btb    = btb,
                ngrad  = ngrad, # = dim(btb)[2]
                s0ind  = s0ind, # indices of S_0 images
                replind = rind,
                ddim   = ddim,
                ddim0  = ddim0,
                xind   = xind,
                yind   = yind,
                zind   = zind,
                level  = level,
                voxelext = voxelext,
                orientation = as.integer(c(0,2,5)),
                source = imagefile)
            )
}

#
#
#


dti <- function(object,  ...) cat("This object has class",class(object),"\n")
setGeneric("dti", function(object,  ...) 
standardGeneric("dti"))


dtiTensor <- function(object,  ...) cat("No DTI tensor calculation defined for this class:",class(object),"\n")


setGeneric("dtiTensor", function(object,  ...) standardGeneric("dtiTensor"))

setMethod("dtiTensor","dtiData",
function(object, method="nonlinear",varmethod="replicates",varmodel="local") {
#  available methods are 
#  "linear" - use linearized model (log-transformed)
#  "nonlinear" - use nonlinear model with parametrization according to Koay et.al. (2006)
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  if(method=="linear"){
     ngrad0 <- ngrad - length(s0ind)
     s0 <- object@si[,,,s0ind]
     si <- object@si[,,,-s0ind]
     if(length(s0ind)>1) s0 <- apply(s0,1:3,mean) 
     mask <- s0 > object@level
     dim(s0) <- dim(si) <- NULL
     ttt <- -log(si/s0)
     ttt[is.na(ttt)] <- 0
     ttt[(ttt == Inf)] <- 0
     ttt[(ttt == -Inf)] <- 0
     dim(ttt) <- c(prod(ddim),ngrad0)
     ttt <- t(ttt)
     cat("Data transformation completed ",date(),proc.time(),"\n")

     btbsvd <- svd(object@btb[,-s0ind])
     solvebtb <- btbsvd$u %*% diag(1/btbsvd$d) %*% t(btbsvd$v)
     D <- solvebtb%*% ttt
     cat("Diffusion tensors generated ",date(),proc.time(),"\n")

     res <- ttt - t(object@btb[,-s0ind]) %*% D
     rss <- res[1,]^2
     for(i in 2:ngrad0) rss <- rss + res[i,]^2
     dim(rss) <- ddim
     sigma2 <- rss/(ngrad0-6)
     D[c(1,4,6),!mask] <- 1e-6
     D[c(2,3,5),!mask] <- 0
     dim(D) <- c(6,ddim)
     dim(res) <- c(ngrad0,ddim)
     cat("Variance estimates generated ",date(),proc.time(),"\n")
     th0 <- array(s0,object@ddim)
     th0[!mask] <- 0
     gc()
  } else {
#  method == "nonlinear" 
     ngrad0 <- ngrad
     si <- aperm(object@si,c(4,1:3))
     s0 <- si[s0ind,,,]
     if(length(s0ind)>1) s0 <- apply(s0,2:4,mean)
     dim(s0) <- ddim
     mask <- s0 > object@level
     cat("start nonlinear regression",date(),proc.time(),"\n")
     z <- .Fortran("nlrdtirg",
                as.integer(si),
                as.integer(ngrad),
                as.integer(ddim[1]),
                as.integer(ddim[2]),
                as.integer(ddim[3]),
                as.logical(mask),
                as.double(object@btb),
                th0=as.double(s0),
                D=double(6*prod(ddim)),
                as.integer(200),
                as.double(1e-6),
                res=double(ngrad*prod(ddim)),
                rss=double(prod(ddim)),
                PACKAGE="dti",DUP=FALSE)[c("th0","D","res","rss")]
     cat("successfully completed nonlinear regression ",date(),proc.time(),"\n")
     dim(z$th0) <- ddim
     dim(z$D) <- c(6,ddim)
     dim(z$res) <- c(ngrad,ddim)
     dim(z$rss) <- ddim
     df <- sum(table(object@replind)-1)
     res <- z$res
     D <- z$D
     rss <- z$rss
     th0 <- z$th0
     rm(z)
     gc()
     cat("Start variance estimation ",date(),proc.time(),"\n")
     if(df<1||varmethod!="replicates"){
        sigma2 <- rss/(ngrad-7)
     } else {
#
#  We may want something more sophisticated here in case of
#  replicated designs !!!
#
        df <- sum(table(object@replind)-1)
        hmax <- max(1,(125/df)^(1/3))
        z <- replvar(si,object@replind)
#
#   need to correct for underestimtion of variances due to 
#   truncation of si to integer
#   standard deviation is underestimated by about 0.385 for sd > 2
#
        z <- (sqrt(z)+0.385)^2
        dim(z) <- ddim
#  adaptive bw to achive approx. 200 degrees of freedom
           sigma2 <- awsvar(z,shape=df,hmax=pmax(1,(125/df)^(1/3)),mask=mask)
     }
  }
  if(varmodel=="global") sigma2 <- array(median(sigma2[sigma2>0]),dim(sigma2))
     cat("successfully completed variance estimation ",date(),proc.time(),"\n")
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
  cat("estimated spatial correlations",date(),proc.time(),"\n")
  cat("first order  correlation in x-direction",signif(scorr[2,1,1],3),"\n")
  cat("first order  correlation in y-direction",signif(scorr[1,2,1],3),"\n")
  cat("first order  correlation in z-direction",signif(scorr[1,1,2],3),"\n")

  scorr[is.na(scorr)] <- 0
  bw <- optim(c(2,2,2),corrrisk,method="L-BFGS-B",lower=c(.2,.2,.2),
  upper=c(3,3,3),lag=lags,data=scorr)$par
  bw[bw <= .25] <- 0
  cat("estimated corresponding bandwidths",date(),proc.time(),"\n")

  invisible(new("dtiTensor",
                D     = D,
                th0   = th0,
                sigma = sigma2,
                scorr = scorr, 
                bw = bw, 
                mask = mask,
                hmax = 1,
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
                source = object@source,
                method = method)
            )
})

#
#
#

create.designmatrix.dti <- function(gradient, bvalue=1) {
  dgrad <- dim(gradient)
  if (dgrad[2]==3) gradient <- t(gradient)
  if (dgrad[1]!=3) stop("Not a valid gradient matrix")

  btb <- matrix(0,6,dgrad[2])
  btb[1,] <- gradient[1,]*gradient[1,]
  btb[4,] <- gradient[2,]*gradient[2,]
  btb[6,] <- gradient[3,]*gradient[3,]
  btb[2,] <- 2*gradient[1,]*gradient[2,]
  btb[3,] <- 2*gradient[1,]*gradient[3,]
  btb[5,] <- 2*gradient[2,]*gradient[3,]

  btb * bvalue
}


#
#
#

dtiIndices <- function(object, ...) cat("No DTI indices calculation defined for this class:",class(object),"\n")

setGeneric("dtiIndices", function(object, ...) standardGeneric("dtiIndices"))

setMethod("dtiIndices","dtiTensor",
function(object, which) {
  ddim <- object@ddim

  z <- .Fortran("dtiind3D",
                as.double(object@D),
                as.integer(object@ddim[1]),
                as.integer(object@ddim[2]),
                as.integer(object@ddim[3]),
                as.logical(object@mask),
                fa=double(prod(object@ddim)),
                md=double(prod(object@ddim)),
                andir=double(3*prod(object@ddim)),
                bary=double(3*prod(object@ddim)),
                DUPL=FALSE,
                PACKAGE="dti")[c("fa","md","andir","bary")]

  invisible(new("dtiIndices",
                fa = array(z$fa,object@ddim),
                md = array(z$md,object@ddim),
                andir = array(z$andir,c(3,object@ddim)),
                bary = array(z$bary,c(3,object@ddim)),
                btb   = object@btb,
                ngrad = object@ngrad, # = dim(btb)[2]
                s0ind = object@s0ind,
                ddim  = object@ddim,
                ddim0 = object@ddim0,
                voxelext = object@voxelext,
                orientation = object@orientation,
                xind  = object@xind,
                yind  = object@yind,
                zind  = object@zind,
                method = object@method,
                level = object@level,
                source= object@source)
            )
})

