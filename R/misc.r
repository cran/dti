kldist <- function(L,eta1,eta2){
#
#  L -number of coils
#  eta1 - m1_1/sigma
#  eta2 - m1_2/sigma
n1 <- length(eta1)
n2 <- length(eta2)
f1 <- (2*L+eta1^2)^2/(2*L+2*eta1^2)
c1 <- (2*L+2*eta1^2)/(2*L+eta1^2)
f2 <- (2*L+eta2^2)^2/(2*L+2*eta2^2)
c2 <- (2*L+2*eta2^2)/(2*L+eta2^2)
lGf1 <- lgamma(f1/2)
lGf2 <- lgamma(f2/2)
flc1 <- f1/2*log(c1)
flc2 <- f2/2*log(c2)
psif1 <- digamma(f1/2)
-outer(lGf1,lGf2,"-")-outer(flc1,flc2,"-")+outer(f1,f2,"-")/2*outer(log(c1)+psif1,rep(1,n2),"*")+
(outer(c1,c2,"/")-1)*outer(f1,rep(1/2,n2),"*")
}


fwhm2bw <- function(hfwhm) hfwhm/sqrt(8*log(2))

replind <- function(gradient){
#
#  determine replications in the design that may be used for 
#  variance estimates
#
  if (dim(gradient)[1]!=3) stop("Not a valid gradient matrix")
  ngrad <- dim(gradient)[2]
  replind <- numeric(ngrad)
  while(any(replind==0)){
     i <- (1:ngrad)[replind==0][1]
     ind <- (1:ngrad)[apply(abs(gradient-gradient[,i]),2,max)==0]
     replind[ind] <- i
  }
  as.integer(replind)
}

Spatialvar.gauss<-function(h,h0,d,interv=1){
#
#   Calculates the factor of variance reduction obtained for Gaussian Kernel and bandwidth h in 
#
#   case of colored noise that was produced by smoothing with Gaussian kernel and bandwidth h0
#
#   Spatialvar.gauss(lkern,h,h0,d)/Spatialvar.gauss(lkern,h,1e-5,d) gives the 
#   a factor for lambda to be used with bandwidth h 
#
#
#  interv allows for further discretization of the Gaussian Kernel, result depends on
#  interv for small bandwidths. interv=1  is correct for kernel smoothing, 
#  interv>>1 should be used to handle intrinsic correlation (smoothing preceeding 
#  discretisation into voxel) 
#
  h0 <- pmax(h0,1e-5)
  h <- pmax(h,1e-5)
  h<-h/2.3548*interv
  if(length(h)==1) h<-rep(h,d)
  ih<-trunc(4*h)
  ih<-pmax(1,ih)
  dx<-2*ih+1
  penl<-dnorm(((-ih[1]):ih[1])/h[1])
  if(d==2) penl<-outer(dnorm(((-ih[1]):ih[1])/h[1]),dnorm(((-ih[2]):ih[2])/h[2]),"*")
  if(d==3) penl<-outer(dnorm(((-ih[1]):ih[1])/h[1]),outer(dnorm(((-ih[2]):ih[2])/h[2]),dnorm(((-ih[3]):ih[3])/h[3]),"*"),"*")
  dim(penl)<-dx
  h0<-h0/2.3548*interv
  if(length(h0)==1) h0<-rep(h0,d)
  ih<-trunc(4*h0)
  ih<-pmax(1,ih)
  dx0<-2*ih+1
  x<- ((-ih[1]):ih[1])/h0[1]
  penl0<-dnorm(((-ih[1]):ih[1])/h0[1])
  if(d==2) penl0<-outer(dnorm(((-ih[1]):ih[1])/h0[1]),dnorm(((-ih[2]):ih[2])/h0[2]),"*")
  if(d==3) penl0<-outer(dnorm(((-ih[1]):ih[1])/h0[1]),outer(dnorm(((-ih[2]):ih[2])/h0[2]),dnorm(((-ih[3]):ih[3])/h0[3]),"*"),"*")
  dim(penl0)<-dx0
  penl0<-penl0/sum(penl0)
  dz<-dx+dx0-1
  z<-array(0,dz)
  if(d==1){
    for(i1 in 1:dx0) {
      ind1<-c(0:(i1-1),(dz-dx0+i1):dz+1)
      ind1<-ind1[ind1<=dz][-1]
      z[-ind1]<-z[-ind1]+penl*penl0[i1]
    }
  } else if(d==2){
    for(i1 in 1:dx0[1]) for(i2 in 1:dx0[2]){
      ind1<-c(0:(i1-1),(dz[1]-dx0[1]+i1):dz[1]+1)
      ind1<-ind1[ind1<=dz[1]][-1]
      ind2<-c(0:(i2-1),(dz[2]-dx0[2]+i2):dz[2]+1)
      ind2<-ind2[ind2<=dz[2]][-1]
      z[-ind1,-ind2]<-z[-ind1,-ind2]+penl*penl0[i1,i2]
    }
  } else if(d==3){
    for(i1 in 1:dx0[1]) for(i2 in 1:dx0[2]) for(i3 in 1:dx0[3]){
      ind1<-c(0:(i1-1),(dz[1]-dx0[1]+i1):dz[1]+1)
      ind1<-ind1[ind1<=dz[1]][-1]
      ind2<-c(0:(i2-1),(dz[2]-dx0[2]+i2):dz[2]+1)
      ind2<-ind2[ind2<=dz[2]][-1]
      ind3<-c(0:(i3-1),(dz[3]-dx0[3]+i3):dz[3]+1)
      ind3<-ind3[ind3<=dz[3]][-1]
      z[-ind1,-ind2,-ind3]<-z[-ind1,-ind2,-ind3]+penl*penl0[i1,i2,i3]
    }
  }
  sum(z^2)/sum(z)^2*interv^d
}

Varcor.gauss<-function(h){
#
#   Calculates a correction for the variance estimate obtained by (IQRdiff(y)/1.908)^2
#
#   in case of colored noise that was produced by smoothing with lkern and bandwidth h
#
  h<-pmax(h/2.3548,1e-5)
  ih<-trunc(4*h)+1
  dx<-2*ih+1
  d<-length(h)
  penl <- dnorm(((-ih[1]):ih[1])/h[1])
  if(d==2) penl <- outer(penl,dnorm(((-ih[2]):ih[2])/h[2]),"*")
  if(d==3) penl <- outer(penl,outer(dnorm(((-ih[2]):ih[2])/h[2]),dnorm(((-ih[3]):ih[3])/h[3]),"*"),"*")
  2*sum(penl)^2/sum(diff(penl)^2)
}


corrrisk <- function(bw,lag,data){
  z <- thcorr3D(bw,lag)
  mean((data-z)^2/outer(outer(1:lag[1],1:lag[2],"*"),1:lag(3),"*"))
}

thcorr3D <- function(bw,lag=rep(5,3)){
  g <- trunc(fwhm2bw(bw)*4)
  gw1 <- dnorm(-(g[1]):g[1],0,fwhm2bw(bw[1]))
  gw2 <- dnorm(-(g[2]):g[2],0,fwhm2bw(bw[2]))
  gw3 <- dnorm(-(g[3]):g[3],0,fwhm2bw(bw[3]))
  gwght <- outer(gw1,outer(gw2,gw3,"*"),"*")
  gwght <- gwght/sum(gwght)
  dgw <- dim(gwght)
  scorr <- .Fortran("thcorr",as.double(gwght),
                    as.integer(dgw[1]),
                    as.integer(dgw[2]),
                    as.integer(dgw[3]),
                    scorr=double(prod(lag)),
                    as.integer(lag[1]),
                    as.integer(lag[2]),
                    as.integer(lag[3]),
                    PACKAGE="dti",DUP=TRUE)$scorr
  # bandwidth in FWHM in voxel units
  dim(scorr) <- lag
  scorr
}

andir2.image <- function(dtobject,slice=1,method=1,quant=0,minfa=NULL,show=TRUE,xind=NULL,yind=NULL,...){
  if(!("dti" %in% class(dtobject))) stop("Not an dti-object")
  if(is.null(dtobject$anindex)) stop("No anisotropy index yet")
  adimpro <- require(adimpro)
  anindex <- dtobject$anindex
  dimg <- dim(anindex)[1:2]
  if(is.null(xind)) xind <- 1:dimg[1]
  if(is.null(yind)) yind <- 1:dimg[2]
  if(is.null(slice)) slice <- 1
  anindex <- anindex[xind,yind,slice]
  dimg <- dim(anindex)[1:2]
  andirection <- dtobject$andirection[,xind,yind,slice]
  anindex[anindex>1]<-0
  anindex[anindex<0]<-0
  dim(andirection)<-c(3,prod(dimg))
  if(is.null(minfa)) minfa <- quantile(anindex,quant)
  if(method==1) {
    andirection[1,] <- abs(andirection[1,])
    andirection[2,] <- abs(andirection[2,])
    andirection[3,] <- abs(andirection[3,])
  } else {
    ind<-andirection[1,]<0
    andirection[,ind] <- - andirection[,ind]
    andirection[2,] <- (1+andirection[2,])/2
    andirection[3,] <- (1+andirection[3,])/2
  }
  andirection <- t(andirection)
  andirection <- andirection*as.vector(anindex)*as.numeric(anindex>minfa)
  dim(andirection)<-c(dimg,3)
  if(adimpro) {
    andirection <- make.image(andirection)
    if(show) show.image(andirection,...)
  } else if(show) {
    dim(anindex) <- dimg
    image(anindex,...)
  }
  invisible(andirection)
} 

andir.image <- function(anindex,andirection,quant=0,minfa=NULL){
  dimg <- dim(anindex)
  anindex[anindex>1]<-0
  anindex[anindex<0]<-0
  dim(andirection)<-c(3,prod(dimg))
  if(is.null(minfa)) minfa <- quantile(anindex,quant)
  andirection[1,] <- abs(andirection[1,])
  andirection[2,] <- abs(andirection[2,])
  andirection[3,] <- abs(andirection[3,])
  andirection <- t(andirection)*as.vector(anindex)*as.numeric(anindex>minfa)
  dim(andirection)<-c(dimg,3)
  show.image(make.image(andirection))
  invisible(NULL)
} 

connect.mask <- function(mask){
  dm <- dim(mask)
  n1 <- dm[1]
  n2 <- dm[2]
  n3 <- dm[3]
  n <- n1*n2*n3
  mask1 <- .Fortran("lconnect",
                 as.logical(mask),
                 as.integer(n1),
                 as.integer(n2),
                 as.integer(n3),
                 as.integer((n1+1)/2),
                 as.integer((n2+1)/2),
                 as.integer((n3+1)/2),
                 integer(n),
                 integer(n),
                 integer(n),
                 mask=logical(n),
                 DUP=FALSE,
                 PACKAGE="dti")$mask
  dim(mask1) <- dm
  mask1
}

sphcoord <- function(ccoord){
#
#  transform cartesian into sherical coordinates
#
  ccoord <- ccoord/sqrt(sum(ccoord^2))
  phi <- atan2(ccoord[2],ccoord[1])
  theta <- atan2(sqrt(ccoord[2]^2+ccoord[1]^2),ccoord[3])
  c(theta,phi)
}


gettriangles <- function(gradients){
  dgrad <- dim(gradients)
  if(dgrad[2]==3) gradients <- t(gradients)
  ngrad <- dim(gradients)[2]
  ndist <- ngrad*(ngrad-1)/2
  z <- .Fortran("distvert",
                 as.double(gradients),
                 as.integer(ngrad),
                 ab=integer(2*ndist),
                 distab=double(ndist),
                 as.integer(ndist),
                 DUPL=FALSE,
                 PACKAGE="dti")[c("ab","distab")]
  o <- order(z$distab)
  distab <- z$distab[o]
  ab <- matrix(z$ab,2,ndist)[,o]
  z <- .Fortran("triedges",
                 as.integer(ab),
                 as.double(distab),
                 iab=integer(ndist),
                 as.integer(ndist),
                 triangles=integer(3*5*ngrad),
                 ntriangles=as.integer(5*ngrad),
                 DUPL=FALSE,
                 PACKAGE="dti")[c("iab","triangles","ntriangles")]
  list(triangles=matrix(z$triangles,3,5*ngrad)[,1:z$ntriangle], edges=ab[,z$iab==2])
}

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

identify.fa <- function(view,slice,xind,yind,zind){
n1 <- switch(view,"sagittal"=length(yind),length(xind))
n2 <- switch(view,"axial"=length(yind),length(zind))
x <- as.vector(outer(1:n1,rep(1,n2),"*"))
y <- as.vector(outer(rep(1,n1),1:n2,"*"))
cat("Please use left mouse click to identify a voxel,\n terminate selection process by right mouse click\n")
z <- identify(x,y,plot=FALSE)
coord <- matrix(0,3,length(z))
if(view=="sagittal"){
coord[1,] <- slice
coord[2,] <- yind[x[z]]
coord[3,] <- zind[y[z]]
} else if(view=="coronal") {
coord[1,] <- xind[x[z]]
coord[2,] <- slice
coord[3,] <- zind[y[z]]
} else {
coord[1,] <- xind[x[z]]
coord[2,] <- yind[y[z]]
coord[3,] <- slice
}
coord
}

mthomogen <- function(object,minw=.1,maxangle=30){
#
#  homogenize Mt-objects (used to construct pseudo-realistic examples)
#
andir <- extract(object,"andir")$andir
mix <- extract(object,"mix")$mix
order <- extract(object,"order")$order
mask <- extract(object,"mask")$mask
ddim <- object@ddim
z <- .Fortran("mthomog",
              as.double(andir),
              mix=as.double(mix),
              order=as.integer(order),
              as.integer(ddim[1]),
              as.integer(ddim[2]),
              as.integer(ddim[3]),
              as.integer(dim(mix)[1]),
              as.logical(mask),
              as.double(minw),
              as.double(maxangle/180*pi),
              as.double(object@voxelext),
              andir=as.double(andir),
              DUPL=TRUE,
              PACKAGE="dti")[c("andir","order","mix")]
object@orient <- array(.Fortran("parofor",
                                as.double(z$andir),
                                as.integer(prod(dim(mix))),
                                orient=double(2*prod(dim(mix))),
                                DUPL=FALSE,
                                PACKAGE="dti")$orient,c(2,dim(mix)))
object@mix <- array(z$mix,dim(mix))
object@order <- array(z$order,ddim)
object
}

vcrossp <- function(a, b) {
   c(a[2] * b[3] - a[3] * b[2],
     a[3] * b[1] - a[1] * b[3],
     a[1] * b[2] - a[2] * b[1])
}
