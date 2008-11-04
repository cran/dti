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
  gamma <- pmin(tobj$gi/hobj$bi,1)
  theta <- gamma*tobj$theta+(1-gamma)*hobj$theta
# force positive variance for positive mean by increasing variance estimate
  if(coef[2] <0){
     coef[1] <- coef[1]+(thmax+thmin)/2*coef[2]
     coef[2] <- 0
  }
  if(coef[1] <0.5){
     coef[2] <- coef[2]+coef[1]/thmax
     coef[1] <- 0.5
  }
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
