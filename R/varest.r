#
#
#      estimate variance parameter in a multicoil system
#
#
awssigmc <- function(y,                 # data
                     steps,             # number of iteration steps for PS
                     mask = NULL,       # data mask, where to do estimation
                     ncoils = 1,        # number of coils for parallel MR image acquisition
                     vext = c( 1, 1),   # voxel extensions
                     lambda = 20,       # adaptation parameter for PS
                     h0 = 2,            # initial bandwidth for first step in PS
                     verbose = FALSE, 
                     sequence = FALSE,  # return estimated sigma for intermediate steps of PS?
                     hadj = 1,          # adjust parameter for density() call for mode estimation
                     q = .25,  # for IQR
                     qni = .8,
                     method=c("VAR","MAD")  # for variance, alternative "MAD" for mean absolute deviation
                     ) {
  method <- match.arg(method)
  ## some functions for pilot estimates
  IQQ <- function (x, q = .25, na.rm = FALSE, type = 7) 
    diff(quantile(as.numeric(x), c(q, 1-q), na.rm = na.rm, names = FALSE, type = type))

  IQQdiff <- function(y, mask, q = .25, verbose = FALSE) {
    cq <- qnorm(1-q)*sqrt(2)*2
    sx <- IQQ( diff(y[mask]), q)/cq
    sy <- IQQ( diff(aperm(y,c(2,1,3))[aperm(mask,c(2,1,3))]), q)/cq
    sz <- IQQ( diff(aperm(y,c(3,1,2))[aperm(mask,c(3,1,2))]), q)/cq
    if(verbose) cat( "Pilot estimates of sigma", sx, sy, sz, "\n")
    min( sx, sy, sz)
  }
  
  estsigma <- function(y, mask, q, L, sigma, verbose=FALSE){
    meany <- mean(y[mask]^2)
    eta <- sqrt(max( 0, meany/sigma^2-2*L))
    m <- sqrt(pi/2)*gamma(L+1/2)/gamma(L)/gamma(3/2)*hyperg_1F1(-1/2,L,-eta^2/2)
    v <- max( .01, 2*L+eta^2-m^2)
    if(verbose) cat( eta, m, v, "\n")
    IQQdiff( y, mask, q)/sqrt(v)
  }
  
  varstats <- sofmchi(ncoils)
  if(length(vext)==3) vext <- vext[2:3]/vext[1]
  ## dimension and size of cubus
  ddim <- dim(y)
  n <- prod(ddim)

  ## test dimension
  if (length(ddim) != 3) stop("first argument should be a 3-dimentional array")

  ## check mask
  if (is.null(mask)) mask <- array(TRUE, ddim)
  if(length(mask) != n) stop("dimensions of data array and mask should coincide")

  ## initial value for sigma_0 
  # sigma <- sqrt( mean( y[mask]^2) / 2 / ncoils)
  sigma <- IQQdiff( y, mask, q, verbose=verbose)
#  cat( "sigmahat1", sigma, "\n")
  sigma <- estsigma( y, mask, q, ncoils, sigma)
#  cat( "sigmahat2", sigma, "\n")
  sigma <- estsigma( y, mask, q, ncoils, sigma)
#  cat( "sigmahat3", sigma, "\n")
  sigma <- estsigma( y, mask, q, ncoils, sigma, verbose=verbose)
#  cat( "sigmahat4", sigma,"\n")

  ## define initial arrays for parameter estimates and sum of weights (see PS)
  th <- array( 1, ddim)
  ni <- array( 1, ddim)
#  y <- y/sigma # rescale to avoid passing sigma to awsvchi
  minlev <- sqrt(2)*gamma(ncoils+.5)/gamma(ncoils)
  if (sequence) sigmas <- lind <- minni <- numeric(steps)
  mc.cores <- setCores(,reprt=FALSE)
  ## iterate PS starting with bandwidth h0
  for (i in 1:steps) {

    h <- h0 * 1.25^((i-1)/3)
    nw <- prod(2*as.integer(h/c(1,vext))+1)
#    cat("nw=",nw,"h/vext=",h/c(1,1/vext),"ih=",as.integer(h/c(1,1/vext)),"\n")
    param <- .Fortran("paramw3",
                      as.double(h),
                      as.double(vext),
                      ind=integer(3*nw),
                      w=double(nw),
                      n=as.integer(nw),
                      DUPL = FALSE,
                      PACKAGE = "dti")[c("ind","w","n")]
    nw <- param$n
    param$ind <- param$ind[1:(3*nw)]
    dim(param$ind) <- c(3,nw)
    param$w   <- param$w[1:nw]    
    fncchi <- fncchiv(th/sigma,varstats)
## correction factor for variance of NC Chi distribution
    ## perform one step PS with bandwidth h
    if(method=="VAR"){
    z <- .Fortran("awsvchi",
                  as.double(y),        # data
                  as.double(th),       # previous estimates
                  ni = as.double(ni/fncchi),
                  as.logical(mask),
                  as.integer(ddim[1]),
                  as.integer(ddim[2]),
                  as.integer(ddim[3]),
                  as.integer(param$ind),
                  as.double(param$w),
                  as.integer(nw),
                  as.double(lambda),
                  as.double(sigma),
                  th = double(n),
                  sy = double(n),
                  DUPL = FALSE,
                  PACKAGE = "dti")[c("ni","th","sy")]
    } else {
    z <- .Fortran("awsadchi",
                  as.double(y),        # data
                  as.double(th),       # previous estimates
                  ni = as.double(ni/fncchi),
                  as.logical(mask),
                  as.integer(ddim[1]),
                  as.integer(ddim[2]),
                  as.integer(ddim[3]),
                  as.integer(param$ind),
                  as.double(param$w),
                  as.integer(nw),
                  as.double(lambda),
                  as.double(sigma),
                  double(nw*mc.cores),
                  as.integer(mc.cores),
                  th = double(n),
                  sy = double(n),
                  DUPL = FALSE,
                  PACKAGE = "dti")[c("ni","th","sy")]       
    }
    ## extract sum of weigths (see PS) and consider only voxels with ni larger then mean
    th <- z$th
    ni <- z$ni
    ni[!mask]<-1
    ind <- (ni > .9999*quantile(ni[ni>1],qni))#&(z$th>sigma*minlev)
## use correction factor for sd of NC Chi distribution
    sy <- z$sy[ind]/fncchis(th[ind]/sigma,varstats)
    ## use the maximal mode of estimated local sd parameters, exclude largest values for better precision
    dsigma <- density( sy[sy>0], n = 4092, adjust = hadj, to = min( max(sy[sy>0]), median(sy[sy>0])*5) )
    sigma <- dsigma$x[dsigma$y == max(dsigma$y)][1]

    if (sequence) {
         sigmas[i] <- sigma
         lind[i] <- sum(ind)
         minni[i] <- min(ni[ind])
         }

    if (verbose) {
      plot(dsigma, main = paste( "estimated sigmas step", i, "h=", signif(h,3)))
      cat( "step", i, "h=", signif( h, 3), "quantiles of ni", signif( quantile(ni[ind]), 3), "mean", signif( mean(ni[ind]), 3), "\n")
      cat( "quantiles of sigma", signif( quantile(sy[sy>0]), 3), "mode", signif( sigma, 3), "\n")
    }
  }
  ## END PS iteration


  ## this is the result (th is expectation, not the non-centrality parameter !!!)
  invisible(list(sigma = if(sequence) sigmas else sigma,
                 theta = th, 
                 lind  = if(sequence) lind else sum(ind), 
                 minni  = if(sequence) minni else min(ni[ind])))
}
awssigmc0 <- function(y,                 # data
                     steps,             # number of iteration steps for PS
                     mask = NULL,       # data mask, where to do estimation
                     ncoils = 1,        # number of coils for parallel MR image acquisition
                     vext = c( 1, 1),   # voxel extensions
                     lambda = 10,       # adaptation parameter for PS
                     h0 = 2,            # initial bandwidth for first step in PS
                     verbose = FALSE, 
                     sequence = FALSE,  # return estimated sigma for intermediate steps of PS?
                     eps = 1e-5,        # accuracy for fixpoint iteration 
                     hadj = 1,          # adjust parameter for density() call for mode estimation
                     q = .25            # for IQR
                     ) {

  ## some functions for pilot estimates
  IQQ <- function (x, q = .25, na.rm = FALSE, type = 7) 
    diff(quantile(as.numeric(x), c(q, 1-q), na.rm = na.rm, names = FALSE, type = type))

  IQQdiff <- function(y, mask, q = .25) {
    cq <- qnorm(1-q)*sqrt(2)*2
    sx <- IQQ( diff(y[mask]), q)/cq
    sy <- IQQ( diff(aperm(y,c(2,1,3))[aperm(mask,c(2,1,3))]), q)/cq
    sz <- IQQ( diff(aperm(y,c(3,1,2))[aperm(mask,c(3,1,2))]), q)/cq
    cat( "Pilot estimates of sigma", sx, sy, sz, "\n")
    min( sx, sy, sz)
  }
  
  estsigma <- function(y, mask, q, L, sigma){
    meany <- mean(y[mask]^2)
    eta <- sqrt(max( 0, meany/sigma^2-2*L))
    m <- sqrt(pi/2)*gamma(L+1/2)/gamma(L)/gamma(3/2)*hyperg_1F1(-1/2,L,-eta^2/2)
    v <- max( .01, 2*L+eta^2-m^2)
    cat( eta, m, v, "\n")
    IQQdiff( y, mask, q)/sqrt(v)
  }

  ## dimension and size of cubus
  ddim <- dim(y)
  n <- prod(ddim)

  ## test dimension
  if (length(ddim) != 3) stop("first argument should be a 3-dimentional array")

  ## check mask
  if (is.null(mask)) mask <- array(TRUE, ddim)
  if(length(mask) != n) stop("dimensions of data array and mask should coincide")

  ## initial value for sigma_0 
  # sigma <- sqrt( mean( y[mask]^2) / 2 / ncoils)
  sigma <- IQQdiff( y, mask, q)
#  cat( "sigmahat1", sigma, "\n")
  sigma <- estsigma( y, mask, q, ncoils, sigma)
#  cat( "sigmahat2", sigma, "\n")
  sigma <- estsigma( y, mask, q, ncoils, sigma)
#  cat( "sigmahat3", sigma, "\n")
  sigma <- estsigma( y, mask, q, ncoils, sigma)
#  cat( "sigmahat4", sigma,"\n")

  ## define initial arrays for parameter estimates and sum of weights (see PS)
  th <- array( 1, ddim)
  ni <- array( 1, ddim)
  y <- y/sigma # rescale to avoid passing sigma to awsvchi2
  ## use squared quantities for chisq method
  y <- y^2
    
  if (sequence) sigmas <- numeric(steps)

  ## iterate PS starting with bandwidth h0
  for (i in 1:steps) {

    h <- h0 * 1.25^((i-1)/3)

    ## perform one step PS with bandwidth h
    z <- .Fortran("awsvchi2",
                  as.double(y),        # data
                  as.double(th),       # previous estimates
                  ni = as.double(ni),
                  as.logical(mask),
                  as.integer(ddim[1]),
                  as.integer(ddim[2]),
                  as.integer(ddim[3]),
                  as.double(lambda),
                  as.integer(ncoils),
                  th = double(n),
                  th2 = double(n),
                  ni2 = double(n),
                  double(n),           # array to precompute lgamma
#                  as.double(sigma^2),
                  as.double(h),
                  as.double(vext),
                  DUPL = FALSE,
                  PACKAGE = "dti")[c("ni","th","th2","ni2")]

    ## extract sum of weigths (see PS) and consider only voxels with ni larger then mean
    ni <- z$ni
    ind <- ni > mean(ni)
    cw <- z$ni2[ind]/ni[ind]^2

      th <- z$th
      m1 <- th[ind]
      mu <- pmax( 1/(1-cw)*(z$th2[ind]-m1^2), 0)
      p <- 2*ncoils
      indt <- m1^2 - mu*ncoils < 0
      s2 <- sqrt((m1-sqrt(pmax( 0, m1^2-mu*ncoils)))/p)

    ## use the maximal mode of estimated local variance parameters, exclude largest values for better precision
    dsigma <- density( s2[s2>0], n = 4092, adjust = hadj, to = min( max(s2[s2>0]), median(s2[s2>0])*5) )
    sigmaf <- dsigma$x[dsigma$y == max(dsigma$y)][1]
      th <- th/sigmaf
      y <- y/sigmaf 
      sigma <- sigma*sqrt(sigmaf)

    if (sequence) sigmas[i] <- sigma

    if (verbose) {
      plot(dsigma, main = paste( "estimated sigmas step", i, "h=", signif(h,3)))
      cat( "step", i, "h=", signif( h, 3), "quantiles of ni", signif( quantile(ni), 3), "mean", signif( mean(ni), 3), "\n")
      cat( "quantiles of sigma", signif( quantile(s2[s2>0]), 3), "mode", signif( sigma, 3), "\n")
    }
  }
  ## END PS iteration

  ## estimate parameter
  eta <- sqrt(pmax( 0, th-2*ncoils)) 
  dim(eta) <- ddim

  ## this is the result
  invisible(list(sigma = if(sequence) sigmas else sigma,
                 theta = eta*sigma))
}


afsigmc <- function(y,                 # data
                    mask = NULL,       # data mask, where to do estimation
                    ncoils = 1,        # number of coils for parallel MR image acquisition
                    vext = c( 1, 1),   # voxel extensions
                    h = 2,             # initial bandwidth for first step in PS
                    verbose = FALSE,
                    hadj = 1           # adjust parameter for density() call for mode estimation
                    ) {

  ## dimension and size of cubus
  ddim <- dim(y)
  n <- prod(ddim)

  ## test dimension
  if (length(ddim) != 3) stop("first argument should be a 3-dimentional array")

  ## check mask
  if (is.null(mask)) mask <- array(TRUE, ddim)
  if(length(mask) != n) stop("dimensions of data array and mask should coincide")

  ## let FORTRAN do the calculation
  sigma <- .Fortran("afvarest",
                    as.double(y),
                    as.integer(ddim[1]),
                    as.integer(ddim[2]),
                    as.integer(ddim[3]),
                    as.logical(mask),
                    as.double(h),
                    as.double(vext),
                    sigma = double(n),
                    DUPL = FALSE,
                    PACKAGE = "dti")$sigma
  sigma <- array( sqrt(sigma), ddim)

  ##  use the maximal mode of estimated local variance parameters, exclude largest values for better precision
  dsigma <- density( sigma[sigma>0], n = 4092, adjust = hadj, to = min( max(sigma[sigma>0]), median(sigma[sigma>0])*5) )
  sigmag <- dsigma$x[dsigma$y == max(dsigma$y)][1]

  if(verbose){
    plot(dsigma, main = paste( "estimated sigmas h=", signif( h, 3)))
    cat("quantiles of sigma", signif( quantile(sigma[sigma>0]), 3), "mode", signif( sigmag, 3), "\n")
  }

  ## this is the estimate
  sigmag
}




#
#    R - function  aws  for likelihood  based  Adaptive Weights Smoothing (AWS)
#    for local constant Gaussian, Bernoulli, Exponential, Poisson, Weibull and  
#    Volatility models                                                         
#
#    emaphazises on the propagation-separation approach 
#
#    Copyright (C) 2006 Weierstrass-Institut fuer
#                       Angewandte Analysis und Stochastik (WIAS)
#
#    Author:  Joerg Polzehl
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
#  USA.
#
#     default parameters:  see function setawsdefaults
#       
awslinsd <- function(y,hmax=NULL,hpre=NULL,h0=NULL,mask=NULL,
                ladjust=1,wghts=NULL,varprop=.1,A0,A1)
{
#
#    first check arguments and initialize
#
homogen <- TRUE
wghts <- NULL

args <- match.call()
dy<-dim(y)
if(length(dy)!=3) stop("Image should be 3D")
#
#   set appropriate defaults
#
if(is.null(wghts)) wghts <- c(1,1,1)
wghts <- wghts[1]/wghts[2:3]
cpar<-setawsdefaults(mean(y),ladjust,hmax,wghts)
if(is.null(mask)) {
    if(length(dy)==0) mask <- rep(TRUE,length(y)) else mask <- array(TRUE,dy)
}
lambda <- cpar$lambda
maxvol <- cpar$maxvol
k <- cpar$k
kstar <- cpar$kstar
hmax <- cpar$hmax
n<-length(y)
# 
#   family dependent transformations 
#
zfamily <- awsgfamily(y,h0,3)
sigma2 <- zfamily$sigma2
h0 <- zfamily$h0
rm(zfamily)
# now check which procedure is appropriate
##  this is the version on a grid
n <- length(y)
n1 <- dy[1]
n2 <- dy[2]
n3 <- dy[3]
#
#    Initialize  for the iteration
#  
zobj<-list(ai=y, bi= rep(1,n), theta= y, fix=!mask)
mae<-NULL
lambda0<-1e50 # that removes the stochstic term for the first step, initialization by kernel estimates
#
#   produce a presmoothed estimate to stabilze variance estimates
#
if(is.null(hpre)) hpre<-20^(1/3)
dlw<-(2*trunc(hpre/c(1,wghts))+1)[1:3]
hobj <- .Fortran("caws03d",as.double(y),
                       as.logical(mask),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       as.double(hpre),
                       theta=as.double(zobj$theta),
                       bi=as.double(zobj$bi),
		       double(prod(dlw)),
		       as.double(wghts),
		       PACKAGE="dti",DUP=FALSE)[c("bi","theta")]
dim(hobj$theta) <- dim(hobj$bi) <- dy
#
#   iteratate until maximal bandwidth is reached
#
cat("Progress:")
total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
while (k<=kstar) {
      hakt0 <- gethani(1,10,1.25^(k-1),c(1,0,0,1,0,1),c(1,wghts),1e-4)
      hakt <- gethani(1,10,1.25^k,c(1,0,0,1,0,1),c(1,wghts),1e-4)
dlw<-(2*trunc(hakt/c(1,wghts))+1)[1:3]
if(any(h0>0)) lambda0<-lambda0*Spatialvar.gauss(hakt0/0.42445/4,h0,3)/Spatialvar.gauss(hakt0/0.42445/4,1e-5,3)
# Correction for spatial correlation depends on h^{(k)} 
hakt0<-hakt
# heteroskedastic Gaussian case
zobj <- .Fortran("cgaws",as.double(y),
                       as.logical(mask),
                       as.double(sigma2),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       hhom=as.double(rep(1,n)),
                       as.double(lambda0),
                       as.double(zobj$theta),
                       bi=as.double(zobj$bi),
		       gi=double(n),
		       gi2=double(n),
                       theta=double(n),
		       double(prod(dlw)),
		       as.double(wghts),
		       PACKAGE="dti",DUP=FALSE)[c("bi","hhom","theta","gi","gi2","hakt")]
dim(zobj$theta)<-dim(zobj$gi)<-dim(zobj$gi2)<-dim(zobj$bi)<-dy
hhom <- zobj$hhom
#
#    Calculate MAE and MSE if true parameters are given in u 
#    this is for demonstration and testing for propagation (parameter adjustments) 
#    only.
#
#   Prepare for next iteration
#
#
#   Create new variance estimate
#
vobj <- awsgsigma2(y,mask,hobj,zobj,varprop,h0,A0,A1)
sigma2 <- vobj$sigma2inv
coef <- vobj$coef
rm(vobj)
lambda0<-lambda
if (max(total) >0) {
      cat(signif(total[k],2)*100,"% . ",sep="")
     }
k <- k+1
gc()
}
cat("\n")
###                                                                       
###            end iterations now prepare results                                                  
###                                 
list(theta = zobj$theta, vcoef=coef, mask=mask)
}
###########################################################################
#
#   Auxialiary functions
#
############################################################################
#
#   transformations for Gaussian case with variance modelling
#
############################################################################
awsgfamily <- function(y,h0,d){
  sigma2 <- max(1,IQRdiff(as.vector(y))^2)
  if(any(h0)>0) sigma2<-sigma2*Varcor.gauss(h0)
#  cat("Estimated variance: ", signif(sigma2,4),"\n")
  sigma2 <- rep(sigma2, length(y))
  dim(sigma2) <- dim(y)
  sigma2 <- 1/sigma2 #  taking the invers yields simpler formulaes 
  list(sigma2=sigma2,h0=h0)
}
############################################################################
#
#  estimate inverse of variances
#
############################################################################
awsgsigma2 <- function(y,mask,hobj,tobj,varprop,h0,thmin,thmax){
## specify sd to be linear to mean
    thrange <- range(y[mask])
#    thmax <- thrange[2]-.2*diff(thrange)
#    thmin <- thrange[1]+.1*diff(thrange)
    ind <- tobj$gi>1.5&mask&tobj$theta>thmin&tobj$theta<thmax
    absresid <- abs((y-tobj$theta)[ind]*tobj$gi[ind]/sqrt(tobj$gi[ind]^2-tobj$gi2[ind]))/.8
#     absresid <- abs((y-tobj$theta)[ind])/.8
    theta <- tobj$theta[ind]
    wght <- (tobj$gi[ind]^2-tobj$gi2[ind])/tobj$gi[ind]^2
    coef <- coefficients(lm(absresid~theta,weights=wght^2))
# force positive variance for positive mean by increasing variance estimate
  if(coef[2] < 0){
     coef[1] <- coefficients(lm(absresid~1,weights=wght^2))
     coef[2] <- 0
  }
  if(coef[1] <0.5){
     coef[2] <- coefficients(lm(absresid~theta-1,weights=wght^2))
     coef[1] <- 0
  }
  gamma <- pmin(tobj$gi/hobj$bi,1)
  theta <- gamma*tobj$theta+(1-gamma)*hobj$theta
  sigma2 <- (coef[1]+coef[2]*theta)^2
  varquantile <- varprop*mean(sigma2[mask])
  sigma2 <- pmax(sigma2,varquantile)
#  cat("Estimated mean variance",signif(mean(sigma2),3)," Variance parameters:",signif(coef,3),"\n")
  list(sigma2inv=1/sigma2,coef=coef)
}
#######################################################################################
#
#        Auxilary functions
#
#######################################################################################
#
#        Set default values
#
# default values for lambda and tau are chosen by propagation condition
# (strong version) with alpha=0.05 (Gaussian) and alpha=0.01 (other family models)
# see script aws_propagation.r
#
#######################################################################################
setawsdefaults <- function(meany,ladjust,hmax,wghts){
qlambda <- .98
if(is.null(hmax)) hmax <- 5
lambda <- ladjust*qchisq(qlambda,1)*2
maxvol <- getvofh(hmax,c(1,0,0,1,0,1),c(1,wghts))
kstar <- as.integer(log(maxvol)/log(1.25))
k <- 6
cat("Estimating variance model using PS with with lambda=",signif(lambda,3)," hmax=",hmax,"number of iterations:",kstar-k+1,"\n")
list(lambda=lambda,hmax=hmax,kstar=kstar,maxvol=maxvol,k=k,wghts=wghts)
}
#######################################################################################
#
#    IQRdiff (for robust variance estimates
#
#######################################################################################
IQRdiff <- function(y) IQR(diff(y))/1.908
