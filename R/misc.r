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

andir2.image <- function(dtobject,slice=1,method=1,quant=0,minanindex=NULL,show=TRUE,xind=NULL,yind=NULL,...){
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
if(is.null(minanindex)) minanindex <- quantile(anindex,quant)
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
andirection <- andirection*as.vector(anindex)*as.numeric(anindex>minanindex)
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

andir.image <- function(anindex,andirection,quant=0,minanindex=NULL){
dimg <- dim(anindex)
anindex[anindex>1]<-0
anindex[anindex<0]<-0
dim(andirection)<-c(3,prod(dimg))
if(is.null(minanindex)) minanindex <- quantile(anindex,quant)
andirection[1,] <- abs(andirection[1,])
andirection[2,] <- abs(andirection[2,])
andirection[3,] <- abs(andirection[3,])
andirection <- t(andirection)*as.vector(anindex)*as.numeric(anindex>minanindex)
dim(andirection)<-c(dimg,3)
show.image(make.image(andirection))
invisible(NULL)
} 

tensor2medinria <- function(obj, filename, xind=NULL, yind=NULL, zind=NULL) {
  if (!require(fmri)) stop("cannot execute function without package fmri, because of missing write.NIFTI() function")

  if (is.null(xind)) xind <- 1:obj@ddim[1]
  if (is.null(yind)) yind <- 1:obj@ddim[2]
  if (is.null(zind)) zind <- 1:obj@ddim[3]

  header <- list()
  header$dimension <- c(5,length(xind),length(yind),length(zind),1,6,0,0)
  header$pixdim <- c(-1, obj@voxelext[1:3], 0, 0, 0, 0)
  header$intentcode <- 1007
  header$datatype <- 16
  header$bitpix <- 192
  header$sclslope <- 1
  header$xyztunits <- "\002" # ???
  header$qform <- 1
  header$sform <- 1
  header$quaternd <- 1
  header$srowx <- c(-2,0,0,0)
  header$srowy <- c(0,2,0,0)
  header$srowz <- c(0,0,2,0)

  write.NIFTI(aperm(obj@D,c(2:4,1))[xind,yind,zind,c(1,2,4,3,5,6)],header,filename)
  return(NULL)
}

medinria2tensor <- function(filename) {
  if (!require(fmri)) stop("cannot execute function without package fmri, because of missing read.NIFTI() function")
  args <- sys.call() 
  data <- read.NIFTI(filename)
 
  invisible(new("dtiTensor",
                call  = list(args),
                D     = aperm(extract.data(data),c(4,1:3))[c(1,2,4,3,5,6),,,],
                sigma = array(0,data$dim[1:3]),
                scorr = array(0,c(1,1,1)),
                bw    = rep(0,3),
                mask  = array(TRUE,data$dim[1:3]),
                method = "unknown",
                hmax  = 1,
                th0   = array(0,dim=data$dim[1:3]),
                btb   = matrix(0,1,1),
                ngrad = as.integer(0), # = dim(btb)[2]
                s0ind = as.integer(0),
                ddim  = data$dim[1:3],
                ddim0 = data$dim[1:3],
                xind  = 1:data$dim[1],
                yind  = 1:data$dim[2],
                zind  = 1:data$dim[3],
                voxelext = data$delta,
                scale = 1,
                source= "unknown")
            )

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

