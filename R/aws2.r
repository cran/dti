#
#   Nonlinear regression; regularized version according to Koay et. al. (2006)
#   this is also based on an a statistical penalty defined using log-likelihood difference
#
dtireg.smooth <- function(object,hmax=5,hinit=1,lambda=20,rho=1,graph=FALSE,slice=NULL,quant=.8,
                         minanindex=NULL,hsig=2.5,lseq=NULL,varmethod="residuals",rician=TRUE,niter=5,varmodel="local",wlse=TRUE,volseq=TRUE){
#
#     lambda and lseq adjusted for alpha=0.2
#
  eps <- 1e-6
  if (graph) {
    adimpro <- require(adimpro)
    if (!adimpro) cat("No graphical output! Install package adimpro from CRAN!\n")
    graph <- graph & adimpro
  }
  args <- match.call()
  s0ind <- object@s0ind
  si <- aperm(object@si,c(4,1:3))
  ngrad <- object@ngrad
  ddim0 <- object@ddim0
  ddim <- object@ddim
  xind <- object@xind
  yind <- object@yind
  zind <- object@zind
  source <- object@source
  btb <- object@btb
  voxelext <- object@voxelext
  if(is.null(voxelext)) vext <- c(1,1,1) else vext <- voxelext/min(voxelext)
  dtobject <- dtiTensor(object,method="nonlinear",varmethod=varmethod,varmodel=varmodel)
  mask <- dtobject@mask
  th0 <- dtobject@th0
  D <- dtobject@D
  sigma2 <- dtobject@sigma
  scorr <- dtobject@scorr
  h0 <- dtobject@bw
  cat("Corresponding bandwiths for specified correlation:",h0,"\n")
  gc()
  dimy <- dim(D)
  if(length(dimy)!=4||dimy[1]!=6) stop("D does not contain 3D diffusion tensor image")
  n1<-dimy[2]
  n2<-dimy[3]
  n3<-dimy[4]
  n<-n1*n2*n3
  sigma2[sigma2<=mean(sigma2)*1e-5]<- mean(sigma2)*1e-5
  z <- .Fortran("projdt2",
                as.double(D),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                D=double(6*n),
                anindex=double(n),
                andirection=double(3*n),
                det=double(n),
                as.double(eps),
                DUP=FALSE,
                PACKAGE="dti")[c("D","anindex","andirection","det")]
  dim(z$D) <- dimy
  z$rss <- array(ngrad*sigma2,dimy[-1])
  z$th0 <- th0
  dim(z$anindex) <-dim(z$det) <- dimy[-1]
  dim(z$andirection) <- c(3,dimy[-1]) 
  z$sigma2hat <- sigma2
  z$sihat <- si
   z$bi <- array(1,dimy[-1])#/z$sigma2hat
   dim(z$sigma2hat) <- dimy[-1]
#
#  initial state for h=1
#
  if(graph){
     oldpar <- par(mfrow=c(3,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
     on.exit(par(oldpar))
     if(is.null(slice)) slice<-n3%/%2
     class(z) <- "dti"
     img<-z$D[1,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dxx: mean",signif(mean(z$D[1,,,][mask]),3),"max",signif(max(z$D[1,,,][mask]),3)))
     img<-z$D[2,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     img[img<rg[1]]<-rg[1]
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dxy: min",signif(min(z$D[2,,,][mask]),3),"max",signif(max(z$D[2,,,][mask]),3)))
     img<-z$D[3,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     img[img<rg[1]]<-rg[1]
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dxz: min",signif(min(z$D[3,,,][mask]),3),"max",signif(max(z$D[3,,,][mask]),3)))
     show.image(make.image(z$anindex[,,slice]),xaxt="n",yaxt="n")
     title(paste("Anisotropy index  range:",signif(min(z$anindex[mask]),3),"-",
                  signif(max(z$anindex[mask]),3)))
     img<-z$D[4,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dyy: mean",signif(mean(z$D[4,,,][mask]),3),"max",signif(max(z$D[4,,,][mask]),3)))
     img<-z$D[5,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     img[img<rg[1]]<-rg[1]
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dyz: min",signif(min(z$D[5,,,][mask]),3),"max",signif(max(z$D[5,,,][mask]),3)))
     andir2.image(z,slice,quant=quant,minanindex=minanindex,xaxt="n",yaxt="n")
     title(paste("Directions (h=1), slice",slice))
     ni <- array(1,dimy[-1])*as.integer(mask)
     show.image(make.image((65535*ni/max(ni))[,,slice]),xaxt="n",yaxt="n")
     title(paste("sum of weights  mean=",signif(1,3)))
     img<-z$D[6,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dzz: mean",signif(mean(z$D[6,,,][mask]),3),"max",signif(max(z$D[6,,,][mask]),3)))
  }
  hincr <- 1.25^(1/3)
  maxvol <- getvofh(hmax,c(1,0,0,1,0,1),vext)
  if(volseq) kstar <- as.integer(log(maxvol)/log(1.25)) else kstar <- as.integer(log(hmax)/log(hincr))
  if(is.null(hinit)){
  hakt0 <- 1
  hakt <- hincr
  hinit <- 1
  } else {
  hakt0 <- max(1,hinit/hincr)
  hakt <- hinit
  }
  steps <- as.integer(log(hmax/hinit)/log(hincr)+1)

  # define lseq
  if (is.null(lseq)) {
# this is optimized for lkern="Gaussian" such that alpha approx 0.04 -- 0.1 and probability of separated points is approx. 1e-4
    lseq <- c(1.5,.9,.8,.85,.8,0.85,.85,0.95,1.25,1.15)# alpha=0.2    abs. deviation in S0
  }
  if (length(lseq)<steps) lseq <- c(lseq,rep(1,steps-length(lseq)))
  lseq <- lseq[1:steps]
  k <- 1
  if(volseq) lambda0 <- lambda else lambda0 <- lambda*lseq[k]
  while(k <= kstar) {
    if(volseq) {
      hakt0 <- gethani(1,10,1.25^(k-1),c(1,0,0,1,0,1),vext,1e-4)
      hakt <- gethani(1,10,1.25^k,c(1,0,0,1,0,1),vext,1e-4)
    } else {
      hakt0 <- hincr^(k-1)
      hakt <- hincr^k
    }
    if (any(h0 >= 0.25)) {
       corrfactor <- Spatialvar.gauss(hakt0/0.42445/4,h0,3) /
       Spatialvar.gauss(h0,1e-5,3) /
       Spatialvar.gauss(hakt0/0.42445/4,1e-5,3)
       lambda0 <- lambda0 * corrfactor
       cat("Correction factor for spatial correlation",signif(corrfactor,3),"\n")
    }
    if(volseq) z <- .Fortran("awsrgdt2",
                    as.integer(si),
                    sihat=as.integer(z$sihat), # needed for statistical penalty
                    double(ngrad*n),# array for predicted Si's from the tensor model 
                    as.integer(ngrad),
                    as.integer(n1),
                    as.integer(n2),
                    as.integer(n3),
                    as.logical(mask),
                    as.double(btb),
                    as.double(sigma2),
                    as.logical(wlse),
                    as.double(z$th0),
                    th0=double(n),
                    as.double(z$D), 
                    D=double(6*n),
                    rss=as.double(z$rss),
                    bi=as.double(z$bi),
                    anindex=as.double(z$anindex),
                    andirection=as.double(z$andirection),
                    det=as.double(z$det),
                    as.double(z$sigma2hat),
                    sigma2hat=double(n),
                    sigma2r=double(n),
                    as.double(1.25^k),
                    as.integer(niter),
                    as.double(vext),
                    as.double(rho),
                    as.double(lambda0),
                    double(ngrad),#swsi
                    double(ngrad),#swsi2
                    double(ngrad),#swsi4
                    double(ngrad),#F
                    as.double(eps),
                    as.logical(rician), # based on x <- seq(0,100,.1) !!!
                    DUP=FALSE,PACKAGE="dti")[c("th0","D","rss","bi","anindex","andirection","det","sigma2hat","sigma2r","sihat")] else z <- .Fortran("awsrgdti",
                    as.integer(si),
                    sihat=as.integer(z$sihat), # needed for statistical penalty
                    double(ngrad*n),# array for predicted Si's from the tensor model 
                    as.integer(ngrad),
                    as.integer(n1),
                    as.integer(n2),
                    as.integer(n3),
                    as.logical(mask),
                    as.double(btb),
                    as.double(sigma2),
                    as.logical(wlse),
                    as.double(z$th0),
                    th0=double(n),
                    as.double(z$D), 
                    D=double(6*n),
                    rss=as.double(z$rss),
                    bi=as.double(z$bi),
                    anindex=as.double(z$anindex),
                    andirection=as.double(z$andirection),
                    det=as.double(z$det),
                    as.double(z$sigma2hat),
                    sigma2hat=double(n),
                    sigma2r=double(n),
                    as.double(hakt),
                    as.integer(niter),
                    as.double(vext),
                    as.double(rho),
                    as.double(lambda0),
                    double(ngrad),#swsi
                    double(ngrad),#swsi2
                    double(ngrad),#swsi4
                    double(ngrad),#F
                    as.double(eps),
                    as.logical(rician), # based on x <- seq(0,100,.1) !!!
                    DUP=FALSE,PACKAGE="dti")[c("th0","D","rss","bi","anindex","andirection","det","sigma2hat","sigma2r","sihat")]
     if(hakt<hsig){
        eta <- (hsig^3 - hakt^3)/hsig^3
        z$sigma2hat <- eta*sigma2+(1-eta)*z$sigma2hat
     }
     dim(z$th0) <- dim(z$rss) <- dim(z$bi) <- dim(z$anindex) <- dim(z$det) <- dim(z$sigma2hat) <- dim(z$sigma2r) <- dimy[-1]
     dim(z$D) <- dimy
     dim(z$andirection) <- c(3,dimy[-1]) 
     if(graph){
     class(z) <- "dti"
     img<-z$D[1,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dxx: mean",signif(mean(z$D[1,,,][mask]),3),"max",signif(max(z$D[1,,,][mask]),3)))
     img<-z$D[2,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     img[img<rg[1]]<-rg[1]
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dxy: min",signif(min(z$D[2,,,][mask]),3),"max",signif(max(z$D[2,,,][mask]),3)))
     img<-z$D[3,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     img[img<rg[1]]<-rg[1]
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dxz: min",signif(min(z$D[3,,,][mask]),3),"max",signif(max(z$D[3,,,][mask]),3)))
     show.image(make.image(z$anindex[,,slice]),xaxt="n",yaxt="n")
     title(paste("Anisotropy index  range:",signif(min(z$anindex[mask]),3),"-",
                  signif(max(z$anindex[mask]),3)))
     img<-z$D[4,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dyy: mean",signif(mean(z$D[4,,,][mask]),3),"max",signif(max(z$D[4,,,][mask]),3)))
     img<-z$D[5,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     img[img<rg[1]]<-rg[1]
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dyz: min",signif(min(z$D[5,,,][mask]),3),"max",signif(max(z$D[5,,,][mask]),3)))
     andir2.image(z,slice,quant=quant,minanindex=minanindex,xaxt="n",yaxt="n")
     title(paste("Directions (h=",signif(hakt,3),"), slice",slice))
     ni<-z$bi[,,slice]*if(wlse) z$sigma2hat[,,slice] else 1
     show.image(make.image(65535*ni/max(ni)),xaxt="n",yaxt="n")
     title(paste("sum of weights  mean=",signif(mean((z$bi*if(wlse) z$sigma2hat else 1)[mask]),3)))
     img<-z$D[6,,,slice]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dzz: mean",signif(mean(z$D[6,,,][mask]),3),"max",signif(max(z$D[6,,,][mask]),3)))
     }
     cat("h=",signif(hakt,3),"Quantiles (.5, .75, .9, .95, 1) of anisotropy index",signif(quantile(z$anindex[mask],c(.5, .75, .9, .95, 1)),3),"\n")
#    c1 <- (prod(h0+1))^(1/3)
#    c1 <- 2.7214286 - 3.9476190*c1 + 1.6928571*c1*c1 - 0.1666667*c1*c1*c1
#    x <- (prod(1.25^(k-1)/c(1,1,1)))^(1/3)
#    scorrfactor <- (c1+x)/(c1*prod(h0+1)+x)
    k <- k+1
#    if(volseq) lambda0 <- lambda*scorrfactor else lambda*lseq[k]*scorrfactor    
     if(volseq) lambda0 <- lambda else lambda*lseq[k]
  }
  invisible(new("dtiTensor",
                list(s2rician=if(rician) z$sigma2r else NULL, ni=z$bi*if(wlse) z$sigma2hat else 1),
                D = z$D,
                th0 = z$th0,
                sigma = z$sigma2hat,
                scorr = dtobject@scorr,
                bw = dtobject@bw,
                mask = dtobject@mask,
                btb   = dtobject@btb,
                hmax  = hmax,
                ngrad = ngrad, # = dim(btb)[2]
                s0ind = object@s0ind,
                ddim  = as.integer(ddim),
                ddim0 = as.integer(ddim0),
                xind  = xind,
                yind  = yind,
                zind  = zind,
                voxelext = object@voxelext,
                source= object@source,
                method= dtobject@method)
            )
}

gethani <- function(x,y,value,a,vext,eps=1e-2){
.Fortran("gethani",
         as.double(x),
         as.double(y),
         as.double(value),
         as.double(a),
         as.double(vext),
         as.double(eps),
         bw=double(1),
         PACKAGE="dti")$bw
}
getvofh <- function(bw,a,vext){
.Fortran("getvofh",
         as.double(a),
         as.double(bw),
         as.double(vext),
         vol=double(1),
         PACKAGE="dti")$vol
}

