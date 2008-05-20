#
#
#   check  propagation  
#
#
dti.propagation <- function(object,hmax=5,hinit=NULL,lambda=47,
                                            rho=1,graph=FALSE,slice=NULL,quant=.8,
                                            minanindex=NULL,zext=1,eps=1e-6,hsig=2.5,lseq=NULL){
  if (graph) {
    adimpro <- require(adimpro)
    if (!adimpro) cat("No graphical output! Install package adimpro from CRAN!\n")
    graph <- graph & adimpro
  }
  args <- match.call()
  si <- object$si
  s0 <- object$s0
  ngrad <- object@ngrad
  ddim0 <- object@ddim0
  ddim <- object@ddim
  xind <- object@xind
  yind <- object@yind
  zind <- object@zind
  source <- object@source
  btb <- object@btb
  Bcov <- btb%*%t(btb)
  btbsvd <- svd(btb)
  solvebtb <- btbsvd$u %*% diag(1/btbsvd$d) %*% t(btbsvd$v)
  solvebtb0 <- btbsvd$u %*% diag(1/btbsvd$d) %*% t(btbsvd$u)
  
  dtobject <- as(object,"dtiTensor")
  y <- dtobject$theta
  sigma2 <- dtobject$sigma
  scorr <- dtobject$scorr
  h0 <- dtobject$bw
  cat("Corresponding bandwiths for specified correlation:",h0,"\n")

  rm(object,dtobject)
  gc()
  dimy <- dim(y)
  if(length(dimy)!=4||dimy[1]!=6) stop("y does not contain 3D diffusion tensor image")
  n1<-dimy[2]
  n2<-dimy[3]
  n3<-dimy[4]
  n<-n1*n2*n3
  sigma2[sigma2<=mean(sigma2)*1e-5]<- mean(sigma2)*1e-5
  z <- .Fortran("projdt2",
                as.double(y),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                theta=double(6*n),
                anindex=double(n),
                andirection=double(3*n),
                det=double(n),
                as.double(eps),
                DUP=FALSE,
                PACKAGE="dti")[c("theta","anindex","andirection","det")]
  y <- array(z$theta,dimy)
  z$bi <- 1/sigma2
  dim(z$theta) <- dimy
  dim(z$anindex) <-dim(z$det) <- dimy[-1]
  dim(z$andirection) <- c(3,dimy[-1]) 
  z$s0hat <- s0
  sigma2hat <- .Fortran("smsigma",
                       as.double(sigma2),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       as.double(hsig),
                       as.double(zext),
                       sigma2hat=double(n1*n2*n3),
                       DUP=FALSE,
                       PACKAGE="dti")$sigma2hat
   z$sigma2hat <- sigma2hat
#
#  initial state for h=1
#
  if(graph){
     oldpar <- par(mfrow=c(3,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
     on.exit(par(oldpar))
     if(is.null(slice)) slice<-n3%/%2
     class(z) <- "dti"
     img<-z$theta[1,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dxx: mean",signif(mean(img),3),"max",signif(max(img),3)))
     img<-z$theta[2,,,slice]
     show.image(make.image(img))
     title(paste("Dxy: min",signif(min(img),3),"max",signif(max(img),3)))
     img<-z$theta[3,,,slice]
     show.image(make.image(img))
     title(paste("Dxz: min",signif(min(img),3),"max",signif(max(img),3)))
     show.image(make.image(z$anindex[,,slice]))
     title(paste("Anisotropy index  range:",signif(min(z$anindex),3),"-",
                  signif(max(z$anindex),3)))
     img<-z$theta[4,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dyy: mean",signif(mean(img),3),"max",signif(max(img),3)))
     img<-z$theta[5,,,slice]
     show.image(make.image(img))
     title(paste("Dyz: min",signif(min(img),3),"max",signif(max(img),3)))
     andir2.image(z,slice,quant=quant,minanindex=minanindex)
     title(paste("Directions (h=1), slice",slice))
     ni<-z$bi[,,slice]*sigma2[,,slice]
     show.image(make.image(65535*ni/max(ni)))
     title(paste("sum of weights  mean=",signif(mean(z$bi*sigma2),3)))
     img<-z$theta[6,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dzz: mean",signif(mean(img),3),"max",signif(max(img),3)))
  }
  hincr <- 1.25^(1/3)
  if(is.null(hinit)){
  hakt0 <- 1
  hakt <- hincr
  } else {
  hakt0 <- max(1,hinit/hincr)
  hakt <- hinit
  }
  steps <- as.integer(log(hmax)/log(hincr)+1)

  # define lseq
  if (is.null(lseq)) {
# this is optimized for lkern="Gaussian" such that alpha approx 0.04 -- 0.1 and probability of separated points is approx. 1e-4
    lseq <- c(1.5,.9,.8,.85,.8,0.85,.85,0.95,1.25,1.15)# alpha=0.2  (absolute deviation in S0)
  }
  if (length(lseq)<steps) lseq <- c(lseq,rep(1,steps-length(lseq)))
  lseq <- lseq[1:steps]
  k <- 1
  lambda0 <- lambda*lseq[k]
  z0 <- z
  while(hakt <= hmax) {
    if (any(h0 >= 0.25)) {
       corrfactor <- Spatialvar.gauss(hakt0/0.42445/4,h0,3) /
       Spatialvar.gauss(h0,1e-5,3) /
       Spatialvar.gauss(hakt0/0.42445/4,1e-5,3)
       lambda0 <- lambda0 * corrfactor
       cat("Correction factor for spatial correlation",signif(corrfactor,3),"\n")
    }
     z <- .Fortran("awssidti",
                as.double(s0),
                as.double(si),
                as.double(z$theta),
                bi=as.double(z$bi),
                anindex=as.double(z$anindex),
                andirection=as.double(z$andirection),
                det=as.double(z$det),
                as.double(Bcov),
                as.double(solvebtb),
                as.double(sigma2),
                as.double(z$sigma2hat),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                as.integer(ngrad),
                as.double(hakt),
                as.double(zext),
                as.double(rho),
                as.double(lambda0),
                theta=double(6*n),
                s0hat=double(n),
                sigma2hat=double(n),
                double(ngrad),
                as.double(eps),
                DUP=FALSE,
                PACKAGE="dti")[c("theta","bi","anindex","andirection","det","s0hat","sigma2hat")]
     z0 <- .Fortran("awssidti",
                as.double(s0),
                as.double(si),
                as.double(z0$theta),
                bi=as.double(z$bi),
                anindex=as.double(z0$anindex),
                andirection=as.double(z0$andirection),
                det=as.double(z0$det),
                as.double(Bcov),
                as.double(solvebtb),
                as.double(sigma2),
                as.double(z0$sigma2hat),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                as.integer(ngrad),
                as.double(hakt),
                as.double(zext),
                as.double(rho),
                as.double(1e50),
                theta=double(6*n),
                s0hat=double(n),
                sigma2hat=double(n),
                double(ngrad),
                as.double(eps),
                DUP=FALSE,
                PACKAGE="dti")[c("theta","bi","anindex","andirection","det","s0hat","sigma2hat")]
     if(hakt<hsig){
        eta <- (hsig^3 - hakt^3)/hsig^3
        z$sigma2hat <- eta*sigma2hat+(1-eta)*z$sigma2hat
     }
     dim(z$s0hat) <- dim(z$bi) <- dim(z$anindex) <- dim(z$det) <- dim(z$sigma2hat) <- dimy[-1]
     dim(z$theta) <- dimy
     dim(z$andirection) <- c(3,dimy[-1]) 
     dim(z0$s0hat) <- dim(z0$bi) <- dim(z0$anindex) <- dim(z0$det) <- dim(z0$sigma2hat) <- dimy[-1]
     dim(z0$theta) <- dimy
     dim(z0$andirection) <- c(3,dimy[-1]) 
     cat("hakt=",hakt,"lseq[k]",lseq[k],"alpha=",risks(z0,z,solvebtb0),mean(z0$bi)/mean(z$bi)-1,risks2(z0,z),risks3(z0,z),"\n")
     if(graph){
     class(z) <- "dti"
     img<-z$theta[1,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dxx: mean",signif(mean(img),3),"max",signif(max(img),3)))
     img<-z$theta[2,,,slice]
     show.image(make.image(img))
     title(paste("Dxy: min",signif(min(img),3),"max",signif(max(img),3)))
     img<-z$theta[3,,,slice]
     show.image(make.image(img))
     title(paste("Dxz: min",signif(min(img),3),"max",signif(max(img),3)))
     show.image(make.image(z$anindex[,,slice]))
     title(paste("Anisotropy index  range:",signif(min(z$anindex),3),"-",
                  signif(max(z$anindex),3)))
     img<-z$theta[4,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dyy: mean",signif(mean(img),3),"max",signif(max(img),3)))
     img<-z$theta[5,,,slice]
     show.image(make.image(img))
     title(paste("Dyz: min",signif(min(img),3),"max",signif(max(img),3)))
     andir2.image(z,slice,quant=quant,minanindex=minanindex)
     title(paste("Directions (h=",signif(hakt,3),"), slice",slice))
     ni<-z$bi[,,slice]*z$sigma2hat[,,slice]
     show.image(make.image(65535*ni/max(ni)))
     title(paste("sum of weights  mean=",signif(mean(z$bi*z$sigma2hat),3)))
     img<-z$theta[6,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dyy: mean",signif(mean(img),3),"max",signif(max(img),3)))
     }
     cat("h=",signif(hakt,3),"Quantiles (.5, .75, .9, .95, 1) of anisotropy index",signif(quantile(z$anindex,c(.5, .75, .9, .95, 1)),3),"\n")
     hakt0<-hakt
     hakt <- hakt*hincr
    c1 <- (prod(h0+1))^(1/3)
    c1 <- 2.7214286 - 3.9476190*c1 + 1.6928571*c1*c1 - 0.1666667*c1*c1*c1
    x <- (prod(1.25^(k-1)))^(1/3)
    scorrfactor <- (c1+x)/(c1*prod(h0+1)+x)
    k <- k+1
    lambda0 <- lambda*lseq[k]*scorrfactor
  }
  dim(z$s0hat) <- c(n1,n2,n3)
#  z <- list(theta=z$theta,bi=z$bi,anindex=z$anindex,andirection=z$andirection,
#            ddim0=ddim0,xind=xind,yind=yind,zind=zind,InvCov=Bcov,s0hat=z$s0hat,call=args)

  invisible(new("dtiTensor",
                list(theta = z$theta, sigma = z$sigma2hat, scorr = scorr),
                btb   = btb,
                ngrad = ngrad, # = dim(btb)[2]
                ddim  = as.integer(ddim),
                ddim0 = as.integer(ddim0),
                xind  = xind,
                yind  = yind,
                zind  = zind,
                source= source)
            )

}
risks <- function(z0,z,sbtb){
theta <- z$theta
theta0 <- z0$theta
theta <- sweep(theta,1,apply(theta0,1,mean),"-")
theta0 <- sweep(theta0,1,apply(theta0,1,mean),"-")
dim(theta) <- dim(theta0) <- c(6,prod(dim(theta)[-1]))
theta <- sbtb%*%theta
theta0 <- sbtb%*%theta0
cat(mean(apply(theta^2,1,sum)),mean(apply(theta0^2,1,sum)),"\n")
alpha <- mean(theta^2)/mean(theta0^2)-1
alpha
}
risks2 <- function(z0,z){
theta <- z$s0hat
theta0 <- z0$s0hat
mean((theta-mean(theta0))^2)/mean((theta0-mean(theta0))^2)-1
}
risks3 <- function(z0,z){
theta <- z$s0hat
theta0 <- z0$s0hat
mean(abs(theta-mean(theta0)))/mean(abs(theta0-mean(theta0)))-1
}
