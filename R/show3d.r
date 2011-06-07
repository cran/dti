################################################################
#                                                              #
# Section for show3d() functions (public)                      #
#                                                              #
################################################################

show3d <- function(obj,  ...) cat("3D Visualization not implemented for this class:",class(obj),"\n")

setGeneric("show3d", function(obj,  ...) standardGeneric("show3d"))

setMethod("show3d","dtiData", function(obj,nx=NULL,ny=NULL,nz=NULL,center=NULL,quant=.8,scale=.4,bgcolor="black",add=FALSE,maxobjects=729,what="ADC",minalpha=1,nn=1,normalize=FALSE,box=FALSE,title=FALSE,...){
  if(!require(rgl)) stop("Package rgl needs to be installed for 3D visualization")
  if(is.null(nx)) nx <- obj@ddim[1]
  if(is.null(ny)) ny <- obj@ddim[2]
  if(is.null(nz)) nz <- obj@ddim[3]
  n <- nx*ny*nz
  if(is.null(center)) center <- floor((obj@ddim+1)/2)
  if(nx*ny*nz>maxobjects) {
  cat("size of data cube",n," exceeds maximum of",maxobjects,"\n")
  if(nz > maxobjects^(1/3)) n3 <- 1 else n3 <- nz
    n1 <- n2 <- floor(sqrt(maxobjects/n3))
  } else {
    n1 <- nx
    n2 <- ny
    n3 <- nz
  }
  xind <- (center[1]-(n1%/%2)):(center[1]+(n1%/%2))
  yind <- (center[2]-(n2%/%2)):(center[2]+(n2%/%2))
  zind <- (center[3]-(n3%/%2)):(center[3]+(n3%/%2))
  xind <- xind[xind>0&xind<=obj@ddim[1]]
  yind <- yind[yind>0&yind<=obj@ddim[2]]
  zind <- zind[zind>0&zind<=obj@ddim[3]]
  if (obj@orientation[1]==1) xind <- min(xind)+max(xind)-xind
  if (obj@orientation[2]==3) yind <- min(yind)+max(yind)-yind
  if (obj@orientation[3]==4) zind <- min(zind)+max(zind)-zind
  n1 <- length(xind)
  n2 <- length(yind)
  n3 <- length(zind)
  n <- n1*n2*n3
  if(n==0) stop("Empty cube specified")
  cat(" selected cube specified by \n xind=",min(xind),":",max(xind),
      "\n yind=",min(yind),":",max(yind),
      "\n zind=",min(zind),":",max(zind),"\n")
  obj <- obj[xind,yind,zind,drop=FALSE]
  vext <- obj@voxelext
  tmean <- array(0,c(3,n1,n2,n3))
  tmean[1,,,] <- xind*vext[1]
  tmean[2,,,] <- outer(rep(1,n1),yind)*vext[2]
  tmean[3,,,] <- outer(rep(1,n1),outer(rep(1,n2),zind))*vext[3]
  dim(tmean) <- c(3,n)
  radii <- extract(obj,"sb")$sb
  s0 <- extract(obj,"s0")$s0
  if(length(dim(s0))==4) s0 <- apply(s0,1:3,mean)
  radii <- sweep(radii,1:3,s0,"/")
  if(what=="ADC") radii <- array(pmax(0,-log(radii)),dim(radii))
# avoid using negative ADC's caused by random effects 
  ngrad <- dim(radii)[length(dim(radii))]
  dim(radii) <- c(length(radii)/ngrad,ngrad)
  radii <- t(radii)
  sscale <- scale
#  if(what=="colorcoded") sscale <- 1
  if(normalize){
     minradii <- apply(radii,2,min)
     maxradii <- apply(radii,2,max)
     radii <- sweep(radii,2,minradii,"-")
     radii <- sweep(radii,2,maxradii-minradii,"/")*sscale
  } else {
     radii <- radii/quantile(apply(radii,2,max),quant)*sscale
  }
  gradient <- obj@gradient[,-obj@s0ind]
  if(!add) {
     open3d()
     par3d(...)
     rgl.bg(color=bgcolor)
     }
#  if(what=="colorcoded") {
#     polyeder <- switch(subdivide+1,icosa0,icosa1,icosa2,icosa3,icosa4)
#     ngrad <- dim(gradient)[2]
#     n <- dim(radii)[2]
#     cat("radii",dim(radii))
#     radii <- matrix(.Fortran("datinter",
#                        as.double(radii),
#                        as.integer(n),
#                        as.double(gradient),
#                        as.integer(ngrad),
#                        as.double(polyeder$vertices),
#                        as.integer(polyeder$nv),
#                        as.integer(nn),#number of nearest neighbors
#                        double(nn),#auxiliary 
#                        integer(nn),#auxiliary 
#                        polyradii=double(polyeder$nv*n),
#                        DUP=FALSE,
#                        PACKAGE="dti")$polyradii,polyeder$nv,n)
#     cat("newradii",dim(radii))
#     show3dCdata(radii,polyeder,centers=tmean,minalpha=minalpha,scale=scale,...)
#  } else {
     show3dData(radii,gradient,centers=tmean,minalpha=minalpha,...)
#  }
  if(box) bbox3d()
  if(is.character(title)) {
     title3d(title,color="white",cex=1.5)
  } else {
     if(title) title3d(switch(tolower(what),"data"="observed DWI data","adc"="observed ADC"),color="white",cex=1.5)
  }
  cat("\n rgl-device",rgl.cur(),switch(tolower(what),"data"="observed diffusion weighted data","adc"="apparent diffusion coefficients from data"),"\n",
  if(normalize) "normalized","\n")
  invisible(rgl.cur())
})
##############

setMethod("show3d","dtiTensor", function(obj,nx=NULL,ny=NULL,nz=NULL,center=NULL,method=1,minfa=.3,mask=NULL,fibers=FALSE,maxangle=30,level=0,quant=.8,scale=.4,bgcolor="black",add=FALSE,subdivide=2,maxobjects=729,what="tensor",odfscale=3,minalpha=.25,normalize=NULL,box=FALSE,title=FALSE,...){
  if(!require(rgl)) stop("Package rgl needs to be installed for 3D visualization")
  if(!exists("icosa0")) data("polyeders")
  if(subdivide<0||subdivide>4) subdivide <- 3
  if(is.null(nx)) nx <- obj@ddim[1]
  if(is.null(ny)) ny <- obj@ddim[2]
  if(is.null(nz)) nz <- obj@ddim[3]
  n <- nx*ny*nz
  if(is.null(center)) center <- floor((obj@ddim+1)/2)
  if(nx*ny*nz>maxobjects) {
  cat("size of data cube",n," exceeds maximum of",maxobjects,"\n")
  if(nz > maxobjects^(1/3)) n3 <- 1 else n3 <- nz
    n1 <- n2 <- floor(sqrt(maxobjects/n3))
  } else {
    n1 <- nx
    n2 <- ny
    n3 <- nz
  }
  xind <- (center[1]-(n1%/%2)):(center[1]+(n1%/%2))
  yind <- (center[2]-(n2%/%2)):(center[2]+(n2%/%2))
  zind <- (center[3]-(n3%/%2)):(center[3]+(n3%/%2))
  xind <- xind[xind>0&xind<=obj@ddim[1]]
  yind <- yind[yind>0&yind<=obj@ddim[2]]
  zind <- zind[zind>0&zind<=obj@ddim[3]]
  if (obj@orientation[1]==1) xind <- min(xind)+max(xind)-xind
  if (obj@orientation[2]==3) yind <- min(yind)+max(yind)-yind
  if (obj@orientation[3]==4) zind <- min(zind)+max(zind)-zind
  n1 <- length(xind)
  n2 <- length(yind)
  n3 <- length(zind)
  n <- n1*n2*n3
  if(n==0) stop("Empty cube specified")
  cat(" selected cube specified by \n xind=",min(xind),":",max(xind),
      "\n yind=",min(yind),":",max(yind),
      "\n zind=",min(zind),":",max(zind),"\n")
  mask <- if(is.null(mask)) obj@mask else obj@mask&mask
  mask[-xind,,] <- FALSE
  mask[,-yind,] <- FALSE
  mask[,,-zind] <- FALSE
#  obj <- obj[xind,yind,zind]
  vext <- obj@voxelext
  center <- center*vext
  n <- prod(obj@ddim) 
  D <- obj@D
  D <- D/max(D)
  dim(D) <- c(6,n)
  mask <- mask & (D[1,]*D[4,]*D[6,]>0)
  tmean <- array(0,c(3,obj@ddim))
  tmean[1,,,] <- (1:obj@ddim[1])*vext[1]
  tmean[2,,,] <- outer(rep(1,obj@ddim[1]),1:obj@ddim[2])*vext[2]
  tmean[3,,,] <- outer(rep(1,obj@ddim[1]),outer(rep(1,obj@ddim[2]),1:obj@ddim[3]))*vext[3]
  dim(tmean) <- c(3,n)
  z <- extract(obj,what=c("andir","fa"))
  if(minfa>0) mask <- mask&(z$fa>=minfa)
  maxev <- extract(obj,what="evalues")$evalues[3,,,,drop=FALSE][mask]
  dim(mask) <- NULL
  andir <- matrix(z$andir,3,n)[,mask]
  fa <- z$fa[mask]
  if(method==1) {
    andir <- abs(andir)
  } else {
    ind<-andir[1,]<0
    andir[,ind] <- - andir[,ind]
    andir[2,] <- (1+andir[2,])/2
    andir[3,] <- (1+andir[3,])/2
  }
  colorvalues <- rgb(andir[1,],andir[2,],andir[3,])
  D <- D[,mask]
  tmean <- tmean[,mask]
  n <- sum(mask)
  if(is.null(normalize)) normalize <- switch(tolower(what),"tensor"=FALSE,"adc"=TRUE,"odf"=FALSE)
  polyeder <- switch(subdivide+1,icosa0,icosa1,icosa2,icosa3,icosa4)
  radii <- .Fortran(switch(tolower(what),tensor="ellradii",adc="adcradii",odf="odfradii"),
                    as.double(polyeder$vertices),
                    as.integer(polyeder$nv),
                    as.double(D),
                    as.integer(n),
                    radii=double(n*polyeder$nv),
                    DUP=FALSE,
                    PACKAGE="dti")$radii
  dim(radii) <- c(polyeder$nv,n)
  if(tolower(what)=="odf") normalize <- FALSE
  if(normalize){
     minradii <- apply(radii,2,min)
     maxradii <- apply(radii,2,max)
     radii <- sweep(radii,2,minradii,"-")
     radii <- sweep(radii,2,maxradii-minradii,"/")*scale
  } else {
    if (tolower(what)=="odf"){
#
#   use a sphere of radius level as baseline for the ODF
#
        radii <- radii+level
#
#   to display results in a form that the volumes are comparable,
#   need volume elements around vertices that have volume proportional to radii,
#   i.e. a radial extension of radii^(1/3)
#   odfscale = 1 corresponds to using radii directly for ODF-values
#   values inbetween are possible
#
        radii <- radii^(1/odfscale)
        radii <- radii/quantile(apply(radii,2,max),quant)*scale
     } else {
     radii <- (radii+level)/(quantile(apply(radii,2,max),quant)+level)*scale
     }
  }
  if(!add) {
     open3d()
     par3d(...)
     rgl.bg(color=bgcolor)
     }
  if(tolower(what)=="odf"){
  show3dODF(radii,polyeder,centers=tmean,minalpha=minalpha,...)
     } else {
  show3dTens(radii,polyeder,centers=tmean,colors=colorvalues,alpha=minalpha+(1-minalpha)*fa)
    }
  if(fibers){
     tracks <- tracking(obj,mask=mask,minfa=minfa,maxangle=maxangle)
     dd <- tracks@fibers
     startind <- tracks@startind
     dd <- expandFibers(dd,startind)$fibers
     rgl.lines(dd[,1]+vext[1]/2,dd[,2]+vext[2]/2,dd[,3]+vext[3]/2,
            color=rgb(abs(dd[,4]),abs(dd[,5]),abs(dd[,6])),
            size=3)
  }  
  if(box) bbox3d()
  if(is.character(title)) {
     title3d(title,color="white",cex=1.5)
  } else {
     if(title) title3d(switch(tolower(what),"tensor"="estimated tensors","adc"="estimated ADC (tensor)"),color="white",cex=1.5)
  }
  cat("\n rgl-device",rgl.cur(),switch(tolower(what),"tensor"="estimated tensors","adc"="apparent diffusion coefficients from estimated tensors"),"\n",
  if(obj@hmax>1) paste("smoothed with hmax=",obj@hmax),if(normalize) "normalized","\n")
  invisible(rgl.cur())
})

setMethod("show3d","dwiMixtensor", function(obj,nx=NULL,ny=NULL,nz=NULL,center=NULL,minfa=.3,minorder=1,mineo=1,fibers=FALSE,maxangle=30,level=0,quant=.8,scale=.4,bgcolor="black",add=FALSE,subdivide=3,maxobjects=729,what="ODF",odfscale=3,minalpha=1,lwd=3,box=FALSE,title=FALSE,...){
  if(!require(rgl)) stop("Package rgl needs to be installed for 3D visualization")
  if(!exists("icosa0")) data("polyeders")
  if(subdivide<0||subdivide>4) subdivide <- 3
  if(is.null(nx)) nx <- obj@ddim[1]
  if(is.null(ny)) ny <- obj@ddim[2]
  if(is.null(nz)) nz <- obj@ddim[3]
  n <- nx*ny*nz
  if(is.null(center)) center <- floor((obj@ddim+1)/2)
  if(nx*ny*nz>maxobjects) {
  cat("size of data cube",n," exceeds maximum of",maxobjects,"\n")
  if(nz > maxobjects^(1/3)) n3 <- 1 else n3 <- nz
    n1 <- n2 <- floor(sqrt(maxobjects/n3))
  } else {
    n1 <- nx
    n2 <- ny
    n3 <- nz
  }
  xind <- (center[1]-(n1%/%2)):(center[1]+(n1%/%2))
  yind <- (center[2]-(n2%/%2)):(center[2]+(n2%/%2))
  zind <- (center[3]-(n3%/%2)):(center[3]+(n3%/%2))
  xind <- xind[xind>0&xind<=obj@ddim[1]]
  yind <- yind[yind>0&yind<=obj@ddim[2]]
  zind <- zind[zind>0&zind<=obj@ddim[3]]
  if (obj@orientation[1]==1) xind <- min(xind)+max(xind)-xind
  if (obj@orientation[2]==3) yind <- min(yind)+max(yind)-yind
  if (obj@orientation[3]==4) zind <- min(zind)+max(zind)-zind
  n1 <- length(xind)
  n2 <- length(yind)
  n3 <- length(zind)
  n <- n1*n2*n3
  if(n==0) stop("Empty cube specified")
  cat(" selected cube specified by \n xind=",min(xind),":",max(xind),
      "\n yind=",min(yind),":",max(yind),
      "\n zind=",min(zind),":",max(zind),"\n")
  mask <- obj@mask
  mask[-xind,,] <- FALSE
  mask[,-yind,] <- FALSE
  mask[,,-zind] <- FALSE
#  obj <- obj[xind,yind,zind]
  vext <- obj@voxelext
  scale <- scale*min(vext)
  center <- center*vext
  order <- obj@order
  ev <- obj@ev
  mix <- obj@mix
  maxorder <- dim(mix)[1]
  orient <- obj@orient
  n <- prod(obj@ddim) 
  tmean <- array(0,c(3,obj@ddim))
  tmean[1,,,] <- (1:obj@ddim[1])*vext[1]
  tmean[2,,,] <- outer(rep(1,obj@ddim[1]),1:obj@ddim[2])*vext[2]
  tmean[3,,,] <- outer(rep(1,obj@ddim[1]),outer(rep(1,obj@ddim[2]),1:obj@ddim[3]))*vext[3]
  if(minfa > 0){
     fa <- extract(obj,"fa")$fa
     mask <- (fa>=minfa)&mask
  }
  if(mineo > 1) mask <- (extract(obj,"eorder")$eorder>=mineo)&mask
  if(minorder > 1) mask <- (order>=minorder)&mask
  dim(ev) <- c(2,n)
  dim(tmean) <- c(3,n)
  dim(orient) <- c(dim(orient)[1:2],n)
  dim(mix) <- c(dim(mix)[1],n)
  if(minfa>0||mineo>1||minorder>1){
    indpos <- (1:n)[mask]
    ev <- ev[,indpos]
    tmean <- tmean[,indpos]
    orient <- orient[,,indpos]
    mix <- mix[,indpos]
    order <- order[indpos]
    fa <- fa[indpos]
    n <- length(indpos)
  }
  gc()
  polyeder <- switch(subdivide+1,icosa0,icosa1,icosa2,icosa3,icosa4)
  if(toupper(what) %in% c("ODF","BOTH")){
  radii <- .Fortran("mixtradi",
                    as.double(polyeder$vertices),
                    as.integer(polyeder$nv),
                    as.double(ev),
                    as.double(orient),
                    as.double(mix),
                    as.integer(order),
                    as.integer(maxorder),
                    as.integer(n),
                    radii=double(n*polyeder$nv),
                    DUP=FALSE,
                    PACKAGE="dti")$radii
  dim(radii) <- c(polyeder$nv,n)
#
#   use a sphere of radius level as baseline for the ODF
#
  radii <- radii+level
#
#   to display results in a form that the volumes are comparable,
#   need volume elements around vertices that have volume proportional to radii,
#   i.e. a radial extension of radii^(1/3)
#   odfscale = 1 corresponds to using radii directly for ODF-values
#   values inbetween are possible
#
  radii <- radii^(1/odfscale)
  radii <- radii/quantile(apply(radii,2,max),quant)*scale
  }
  if(toupper(what) %in% c("AXIS","BOTH")){
   colors <- rainbow(1024,end=2/3)
   ranger <- range(fa)
   ind <- 1024-(fa-ranger[1])/(ranger[2]-ranger[1])*1023
   colorvalues <- rep(colors[ind],rep(2*dim(mix)[1],length(ind)))
  andir <- array(.Fortran("mixandir",
                    as.double(orient),
                    as.double(mix),
                    as.integer(order),
                    as.integer(maxorder),
                    as.integer(n),
                    andir=double(3*n*dim(mix)[1]),
                    DUP=FALSE,
                    PACKAGE="dti")$andir,c(3,dim(mix)[1],n1,n2,n3))
  lcoord <- array(0,c(3,2*dim(mix)[1],n))
  for(i in 1:dim(mix)[1]){
       lcoord[,2*i-1,] <-  andir[,i,]*scale+tmean
       lcoord[,2*i,] <-  -andir[,i,]*scale+tmean
     }          
     dim(lcoord) <- c(3,2*dim(mix)[1]*n)
  }
  dim(tmean) <- c(3,n)
  if(!add) {
     open3d()
     par3d(...)
     rgl.bg(color=bgcolor)
     }
  if(toupper(what) %in% c("ODF","BOTH")) show3dODF(radii,polyeder,centers=tmean,minalpha=minalpha,...)
  if(toupper(what) %in% c("AXIS","BOTH"))  rgl.lines(lcoord[1,],lcoord[2,],lcoord[3,],color=colorvalues,size=lwd)
  if(fibers){
     tracks <- tracking(obj,mask=mask,minfa=minfa,maxangle=maxangle)
     dd <- tracks@fibers
     startind <- tracks@startind
     dd <- expandFibers(dd,startind)$fibers
     rgl.lines(dd[,1]+vext[1]/2,dd[,2]+vext[2]/2,dd[,3]+vext[3]/2,
            color=rgb(abs(dd[,4]),abs(dd[,5]),abs(dd[,6])),
            size=lwd)
  }
  if(box) bbox3d()
  if(is.character(title)) {
     title3d(title,color="white",cex=1.5)
  } else {
     if(title) title3d(switch(tolower(what),"ODF"="estimated ODF"),color="white",cex=1.5)
  }
  cat("\n rgl-device",rgl.cur(),switch(tolower(what),"ODF"="estimated ODF"),"\n")
  if(obj@hmax>1) paste("smoothed with hmax=",obj@hmax,"\n")
  invisible(rgl.cur())
})
##############

setMethod("show3d","dtiIndices",function(obj, index="FA", nx=NULL, ny=NULL, nz=NULL, center=NULL, method=1, minfa=0, bgcolor="black", add=FALSE, lwd=1,box=FALSE,title=FALSE,...){
  if(!require(rgl)) stop("Package rgl needs to be installed for 3D visualization")
  index <- tolower(index) 
  if(!(index%in%c("fa","ga"))) stop("index should be either 'FA' or 'GA'\n")
  if(is.null(center)) center <- floor((obj@ddim+1)/2)
  if(is.null(nx)) nx <- obj@ddim[1]
  if(is.null(ny)) ny <- obj@ddim[2]
  if(is.null(nz)) nz <- obj@ddim[3]
  xind <- (center[1]-(nx%/%2)):(center[1]+(nx%/%2))
  yind <- (center[2]-(ny%/%2)):(center[2]+(ny%/%2))
  zind <- (center[3]-(nz%/%2)):(center[3]+(nz%/%2))
  xind <- xind[xind>0&xind<=obj@ddim[1]]
  yind <- yind[yind>0&yind<=obj@ddim[2]]
  zind <- zind[zind>0&zind<=obj@ddim[3]]
  if (obj@orientation[1]==1) xind <- min(xind)+max(xind)-xind
  if (obj@orientation[2]==3) yind <- min(yind)+max(yind)-yind
  if (obj@orientation[3]==4) zind <- min(zind)+max(zind)-zind
  n1 <- length(xind)
  n2 <- length(yind)
  n3 <- length(zind)
  n <- n1*n2*n3
  vext <- obj@voxelext
  ind <- switch(index,"fa"=obj@fa[xind,yind,zind], "ga"=obj@ga[xind,yind,zind])
  ind[ind<minfa] <- 0
  ind <- ind*min(vext)
  tmean <- array(0,c(3,n1,n2,n3))
  tmean[1,,,] <- xind*vext[1]
  tmean[2,,,] <- outer(rep(1,n1),yind)*vext[2]
  tmean[3,,,] <- outer(rep(1,n1),outer(rep(1,n2),zind))*vext[3]
  andir <- obj@andir[,xind,yind,zind]
  if(method==1) {
    andir <- abs(andir)
    dim(andir) <- c(3,n1*n2*n3)
  } else {
    ind1 <- andir[1,,]<0
    dim(andir) <- c(3,n1*n2*n3)
    andir[,ind1] <- - andir[,ind1]
    andir[2,] <- (1+andir[2,])/2
    andir[3,] <- (1+andir[3,])/2
  }
  colorvalues <- rgb(andir[1,],andir[2,],andir[3,])
  dim(andir) <- c(3,n1,n2,n3)
  andir <- sweep(obj@andir[,xind,yind,zind],2:4,ind,"*")
  lcoord <- array(0,c(3,2,n1,n2,n3))
  lcoord[,1,,,] <-  andir/2+tmean[,,,,drop=FALSE]
  lcoord[,2,,,] <-  -andir/2+tmean[,,,,drop=FALSE]
  dim(lcoord) <- c(3,2*n1*n2*n3)
  colorvalues <- c(rbind(colorvalues,colorvalues))
  if(!add) {
    open3d()
    par3d(...)
    rgl.bg(color=bgcolor)
  }
  rgl.lines(lcoord[1,],lcoord[2,],lcoord[3,],color=colorvalues,size=lwd)
  if(box) bbox3d()
  if(is.character(title)) {
     title3d(title,color="white",cex=1.5)
  } else {
     if(title) title3d("Main directions",color="white",cex=1.5)
  }
  cat("\n rgl-device",rgl.cur(),"Main directions of diffusion estimated from the tensor model\n\n")
  if(box) bbox3d()
  invisible(rgl.cur())
})

##############

setMethod("show3d","dwiQball", function(obj,nx=NULL,ny=NULL,nz=NULL,center=NULL,level=0,quant=.8,scale=.4,odfscale=3,bgcolor="black",add=FALSE,subdivide=3,maxobjects=729,minalpha=1,box=FALSE,title=FALSE,...){
  if(!require(rgl)) stop("Package rgl needs to be installed for 3D visualization")
  if(!exists("icosa0")) data("polyeders")
  if(subdivide<0||subdivide>4) subdivide <- 3
  if(is.null(nx)) nx <- obj@ddim[1]
  if(is.null(ny)) ny <- obj@ddim[2]
  if(is.null(nz)) nz <- obj@ddim[3]
  n <- nx*ny*nz
  if(is.null(center)) center <- floor((obj@ddim+1)/2)
  if(nx*ny*nz>maxobjects) {
  cat("size of data cube",n," exceeds maximum of",maxobjects,"\n")
  if(nz > maxobjects^(1/3)) n3 <- 1 else n3 <- nz
    n1 <- n2 <- floor(sqrt(maxobjects/n3))
  } else {
    n1 <- nx
    n2 <- ny
    n3 <- nz
  }
  xind <- (center[1]-(n1%/%2)):(center[1]+(n1%/%2))
  yind <- (center[2]-(n2%/%2)):(center[2]+(n2%/%2))
  zind <- (center[3]-(n3%/%2)):(center[3]+(n3%/%2))
  xind <- xind[xind>0&xind<=obj@ddim[1]]
  yind <- yind[yind>0&yind<=obj@ddim[2]]
  zind <- zind[zind>0&zind<=obj@ddim[3]]
  if (obj@orientation[1]==1) xind <- min(xind)+max(xind)-xind
  if (obj@orientation[2]==3) yind <- min(yind)+max(yind)-yind
  if (obj@orientation[3]==4) zind <- min(zind)+max(zind)-zind
  n1 <- length(xind)
  n2 <- length(yind)
  n3 <- length(zind)
  n <- n1*n2*n3
  if(n==0) stop("Empty cube specified")
  cat(" selected cube specified by \n xind=",min(xind),":",max(xind),
      "\n yind=",min(yind),":",max(yind),
      "\n zind=",min(zind),":",max(zind),"\n")
  obj <- obj[xind,yind,zind]
  vext <- obj@voxelext
  center <- center*vext
  sphcoef <- obj@sphcoef
  dim(sphcoef) <- c(dim(sphcoef)[1],prod(dim(sphcoef)[-1]))
  tmean <- array(0,c(3,n1,n2,n3))
  tmean[1,,,] <- xind*vext[1]
  tmean[2,,,] <- outer(rep(1,n1),yind)*vext[2]
  tmean[3,,,] <- outer(rep(1,n1),outer(rep(1,n2),zind))*vext[3]
  dim(tmean) <- c(3,n)
  polyeder <- switch(subdivide+1,icosa0,icosa1,icosa2,icosa3,icosa4)
  sphdesign <- design.spheven(obj@order,polyeder$vertices,obj@lambda)$design
  radii <- t(sphdesign)%*%sphcoef
  radii <- array(pmax(0,radii),dim(radii))
#  avoid negative ODF's, otherwise scaling by volume produces
#  strange results
  mradii <- apply(radii,2,mean)
#
#   rescale and use a sphere of radius level as baseline for the ODF
#
  radii <- sweep(radii,2,mradii,"/")+level
  radii[is.na(radii)] <- 0
#
#   to display results in a form that the volumes are comparable,
#   need volume elements around vertices that have volume proportional to radii,
#   i.e. a radial extension of radii^(1/3)
#   odfscale = 1 corresponds to using radii directly for ODF-values
#   values inbetween are possible
#
  radii <- radii^(1/odfscale)
  radii <- radii/quantile(apply(radii,2,max),quant)*scale
  if(!add) {
     open3d()
     par3d(...)
     rgl.bg(color=bgcolor)
  }
  show3dODF(radii,polyeder,centers=tmean,minalpha=minalpha,...)
  if(box) bbox3d()
  if(is.character(title)) {
     title3d(title,color="white",cex=1.5)
  } else {
     if(title) title3d(switch(tolower(obj@what),"ODF"="ODF","wODF"="Weighted ODF","aODF"="alternative ODF","adc"="ADC (Sph. Harmonics)"),color="white",cex=1.5)
  }
  cat("\n rgl-device",rgl.cur(),switch(tolower(obj@what),"ODF"="Estimated orientation density function (Qball)","aODF"="Estimated orientation density function (Qball)","adc"="estimated apparent diffusion coefficients (sperical harmonics","wODF"="Estimated orientation density function (Aganji et.al. 2009)"),"\n")
  invisible(rgl.cur())
})

setMethod("show3d","dwiFiber", function(obj,add=FALSE,bgcolor="black",box=FALSE,title=FALSE,lwd=1,...){
  if(!require(rgl)) stop("Package rgl needs to be installed for 3D visualization")
  if(!add) {
     open3d()
     par3d(...)
     rgl.bg(color=bgcolor)
  }
  dd <- obj@fibers
  startind <- obj@startind
  dd <- expandFibers(dd,startind)$fibers
  rgl.lines(dd[,1],dd[,2],dd[,3],
            color=rgb(abs(dd[,4]),abs(dd[,5]),abs(dd[,6])),
            size=lwd)
  if(box) bbox3d()
  if(is.character(title)) {
     title3d(title,color="white",cex=1.5)
  } else {
     if(title) title3d("Fiber tracks",color="white",cex=1.5)
  }

  invisible(rgl.cur())
})

################################################################
#                                                              #
# Section for show3d() functions (misc)                        #
#                                                              #
################################################################

show3dTens <- function(radii,polyeder,centers=NULL,colors=NULL,alpha=1,...){
   if(is.null(centers)){
      centers <- matrix(0,3,1)
      n <- 1
   } else {
      dcenters <- dim(centers)
      if(length(dcenters)!=2 || dcenters[1]!=3) stop("centers needs to be NULL or a matrix 
      with dimension (3,n)")
      n <- dcenters[2]
   }
   if(is.null(colors)){
      colors <- heat.colors(1)
   } 
   if(length(colors)!=n){
      nc <- length(colors)
      nnc <- n%/%nc+1
      colors <- rep(colors,nnc)[1:n]
   }
   if(is.null(alpha)){
      alpha <- 1
   } 
   if(length(alpha)!=n){
      nc <- length(alpha)
      nnc <- n%/%nc+1
      alpha <- rep(alpha,nnc)[1:n]
   }
   nv <- polyeder$nv
   ni <- polyeder$ni*3
   colors <- t(matrix(colors,n,ni))
   alpha <- t(matrix(alpha,n,ni))
   vertices <- array(polyeder$vertices,c(3,nv,n))
   indices <- matrix(polyeder$indices,c(ni,n))
   if(length(radii)!=nv*n) stop("wrong length of radii, needs to be 
             dim(polyeder$vertices)[2]*dim(centers)[2]")
   vertices <- sweep(vertices,2:3,radii,"*")
   vertices <- sweep(vertices,c(1,3),centers,"+")
   dim(vertices) <- c(3,nv*n)
   indices <- sweep(matrix(indices,ni,n),2,((1:n)-1)*nv,"+")
   rgl.triangles(vertices[1,indices],vertices[2,indices],vertices[3,indices],
                 color=colors,alpha=alpha,...)
}

#############

show3dData <- function(radii,vertices,centers=NULL,minalpha=1,...){
#
#   use gradients directly
#
   if(is.null(centers)){
      centers <- matrix(0,3,1)
      n <- 1
   } else {
      dcenters <- dim(centers)
      if(length(dcenters)!=2 || dcenters[1]!=3) stop("centers needs to be NULL or a matrix 
      with dimension (3,n)")
      n <- dcenters[2]
   }
   maxradii <- apply(radii,2,max)
   alpha <- minalpha+(1-minalpha)*sweep(radii,2,maxradii,"/")
   colors <- rgb(abs(vertices[1,]),abs(vertices[2,]),abs(vertices[3,]))
   if(length(alpha)!=n){
      nc <- length(alpha)
      nnc <- n%/%nc+1
      alpha <- rep(alpha,nnc)[1:n]
   }
   nv <- dim(vertices)[2]
   lines <- array(0,c(2,3,n,nv))
   vertices <- array(vertices,c(3,nv,n))
   colors <- matrix(colors,nv,n)
   ind1 <- rep(1:(n*nv),rep(2,n*nv))
   if(length(radii)!=nv*n) stop("wrong length of radii, needs to be 
             dim(vertices)[2]*dim(centers)[2]")
   vertices0 <- sweep(vertices,2:3,radii,"*")
   vertices <- sweep(.9*vertices0,c(1,3),centers,"+")
   lines[1,,,] <- aperm(vertices,c(1,3,2))
   vertices <- sweep(vertices0,c(1,3),centers,"+")
   lines[2,,,] <- aperm(vertices,c(1,3,2))
   colors <- array(colors,c(nv,n))
#   rgl.lines(lines[,1,,],lines[,2,,],lines[,3,,],color=t(colors)[ind1],lwd=1)
   rgl.lines(lines[,1,,],lines[,2,,],lines[,3,,],lwd=1)
   rgl.points(vertices[1,,],vertices[2,,],vertices[3,,],color=colors,size=4)
   vertices <- sweep(-.9*vertices0,c(1,3),centers,"+")
   lines[1,,,] <- aperm(vertices,c(1,3,2))
   vertices <- sweep(-vertices0,c(1,3),centers,"+")
   lines[2,,,] <- aperm(vertices,c(1,3,2))
#   rgl.lines(lines[,1,,],lines[,2,,],lines[,3,,],color=t(colors)[ind1],lwd=1)
   rgl.lines(lines[,1,,],lines[,2,,],lines[,3,,],lwd=1)
   rgl.points(vertices[1,,],vertices[2,,],vertices[3,,],color=colors,size=4)
}

#############

show3dData0 <- function(radii,vertices,centers=NULL,minalpha=1,...){
   if(is.null(centers)){
      centers <- matrix(0,3,1)
      n <- 1
   } else {
      dcenters <- dim(centers)
      if(length(dcenters)!=2 || dcenters[1]!=3) stop("centers needs to be NULL or a matrix 
      with dimension (3,n)")
      n <- dcenters[2]
   }
   maxradii <- apply(radii,2,max)
   alpha <- minalpha+(1-minalpha)*sweep(radii,2,maxradii,"/")
   colors <- rgb(abs(vertices[1,]),abs(vertices[2,]),abs(vertices[3,]))
   if(length(alpha)!=n){
      nc <- length(alpha)
      nnc <- n%/%nc+1
      alpha <- rep(alpha,nnc)[1:n]
   }
   nv <- dim(vertices)[2]
   indices <- gettriangles(vertices)$triangles
   ni <- length(indices)
   indices <- sweep(matrix(indices,ni,n),2,((1:n)-1)*nv,"+")
   colors <- matrix(colors,nv,n)
   alpha <- matrix(alpha,nv,n)
   vertices <- array(vertices,c(3,nv,n))
   if(length(radii)!=nv*n) stop("wrong length of radii, needs to be 
             dim(vertices)[2]*dim(centers)[2]")
   vertices <- sweep(vertices,2:3,radii,"*")
   vertices <- sweep(vertices,c(1,3),centers,"+")
   dim(vertices) <- c(3,nv*n)
   rgl.triangles(vertices[1,indices],vertices[2,indices],vertices[3,indices],
                 color=colors[indices],alpha=alpha[indices],...)
}
#############

show3dCdata <- function(radii,polyeder,centers=NULL,minalpha=1,scale=.5,...){
   if(is.null(centers)){
      centers <- matrix(0,3,1)
      n <- 1
   } else {
      dcenters <- dim(centers)
      if(length(dcenters)!=2 || dcenters[1]!=3) stop("centers needs to be NULL or a matrix 
      with dimension (3,n)")
      n <- dcenters[2]
   }
   vertices <- polyeder$vertices
   colors <- rainbow(1000,start=0,end=.7)[pmax(1,pmin(1000,as.integer(1000*radii)))]
   alpha <- minalpha+(1-minalpha)*radii
   nv <- polyeder$nv
   ni <- polyeder$ni*3
   colors <- matrix(colors,nv,n)
   alpha <- matrix(alpha,nv,n)
   vertices <- array(polyeder$vertices,c(3,nv,n))
   indices <- matrix(polyeder$indices,c(ni,n))
   if(length(radii)!=nv*n) stop("wrong length of radii, needs to be 
             dim(polyeder$vertices)[2]*dim(centers)[2]")
   vertices <- sweep(scale*vertices,c(1,3),centers,"+")
   dim(vertices) <- c(3,nv*n)
   indices <- sweep(matrix(indices,ni,n),2,((1:n)-1)*nv,"+")
   rgl.triangles(vertices[1,indices],vertices[2,indices],vertices[3,indices],
                 color=colors[indices],alpha=alpha[indices],...)
}

#############

show3dODF <- function(radii,polyeder,centers=NULL,minalpha=1,...){
   if(is.null(centers)||length(centers)==3){
      centers <- matrix(0,3,1)
      n <- 1
   } else {
      dcenters <- dim(centers)
      if(length(dcenters)!=2 || dcenters[1]!=3) stop("centers needs to be NULL or a matrix 
      with dimension (3,n)")
      n <- dcenters[2]
   }
   vertices <- polyeder$vertices
   alpha <- rep(minalpha,length(radii))
   nv <- polyeder$nv
   ni <- polyeder$ni*3
   colors <- rainbow(1024,end=2/3)
   rradii <- apply(radii,2,range)
   sradii <- sweep(sweep(radii,2,rradii[1,],"-"),2,rradii[2,]-rradii[1,],"/")
   ind <- 1024-sradii*1023
   alpha <- matrix(alpha,nv,n)
   vertices <- array(polyeder$vertices,c(3,nv,n))
   indices <- array(polyeder$indices,c(ni,n))
   if(length(radii)!=nv*n) stop("wrong length of radii, needs to be 
             dim(polyeder$vertices)[2]*dim(centers)[2]")
   vertices <- sweep(vertices,2:3,radii,"*")
   vertices <- sweep(vertices,c(1,3),centers,"+")
   dim(vertices) <- c(3,nv*n)
   indices <- sweep(matrix(indices,ni,n),2,((1:n)-1)*nv,"+")
   rgl.triangles(vertices[1,indices],vertices[2,indices],vertices[3,indices],
                 color=colors[ind][indices],alpha=alpha[indices],...)
}
