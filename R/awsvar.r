#
#    R - function  awsvar  for likelihood  based  Adaptive Weights Smoothing (AWS)
#    for 3D local constant variance models                                                         
#
#    emaphazises on the propagation-separation approach 
#
#    Copyright (C) 2007 Weierstrass-Institut fuer
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
awsvar <- function(y,hmax=NULL,qlambda=NULL,
                shape=NULL,wghts=NULL,
		spmin=0.25,lseq=NULL,u=NULL,testprop=FALSE,mask=NULL)
{
#
#    first check arguments and initialize
#
args <- match.call()
dy<-dim(y)
if(length(dy)>3) stop("AWS for more than 3 dimensional grids is not implemented")
#
#   set appropriate defaults
#
lkern<-2
cpar<-setawsdefaults(dy,mean(y),qlambda,lseq,hmax,shape)
lambda <- cpar$lambda
hmax <- cpar$hmax
lseq <- cpar$lseq
shape <- cpar$shape
d <- cpar$d
hinit <- cpar$hinit
hincr <- cpar$hincr
n<-length(y)
# 
#   family dependent transformations that depend on the value of family
#
lambda <- 2*lambda/shape 
# this accounts for the additional 1/2 in Q(\hat{theta},theta) and the degrees of freedom in chisq

# now check which procedure is appropriate
##  this is the version on a grid
n <- length(y)
n1 <- switch(d,n,dy[1],dy[1])
n2 <- switch(d,1,dy[2],dy[2])
n3 <- switch(d,1,1,dy[3])
#
#    Initialize  for the iteration
#
if(is.null(mask)||length(mask)!=n) mask <- rep(TRUE,n)
if(is.null(wghts)) wghts<-c(1,1,1)
hinit<-hinit/wghts[1]
hmax<-hmax/wghts[1]
wghts<-(wghts[2:3]/wghts[1])
tobj<-list(bi= rep(1,n),  thn= y/shape)
hakt <- hinit
hakt0 <- hinit
lambda0<-1e50 # that removes the stochstic term for the first step, Initialization by kernel estimates
if(testprop) {
#
#  prepare to  check for alpha in propagation condition (to adjust qlambda and lseq)
#
       if(is.null(u)) u <- 0
       cpar <- c(cpar, list(n1=n1,n2=n2,n3=n3,n=n1*n2*n3,lkern=lkern,wghts=wghts,spmin=spmin,u=u))
       propagation <- NULL
    } 
#
#   iteratate until maximal bandwidth is reached
#
steps <- as.integer(log(hmax/hinit)/log(hincr))
cat("Progress:")
for(k in 0:steps){
dlw<-(2*trunc(hakt/c(1,wghts))+1)[1:d]
hakt0<-hakt
{
# all other cases
tobj <- .Fortran("cawsvar",as.double(y),
                       as.logical(mask),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(tobj$thn),
                       double(n1*n2*n3),
                       bi=as.double(tobj$bi),
                       thn=double(n),
                       as.double(spmin),
		       double(prod(dlw)),
		       as.double(wghts),
		       PACKAGE="dti",DUP=FALSE)[c("bi","thn","hakt")]
}
dim(tobj$thn)<-dy
dim(tobj$bi)<-dy
#
#  if testprop == TRUE
#  check alpha in propagation condition (to adjust qlambda and lseq)
#  
if(testprop) propagation <- awstestprop(y,tobj,mask,hakt,cpar,u,propagation)
#
#   Prepare for next iteration
#
hakt <- hakt*hincr
lambda0<-lambda*lseq[k+1]
cat(paste(signif(sum(hincr^(2*(0:k)))/sum(hincr^(2*(0:steps)))*100,2),"% ",sep=""))
gc()
}
cat("\n")
###                                                                       
###            end iterations now prepare results                                                  
###                                 
array(tobj$thn,dim(y))
}
#######################################################################################
#
#        Auxilary functions
#
#######################################################################################
#
#        Set default values
#
# default values for qlambda, lseq  are chosen by propagation condition
# (strong version) with alpha=0.1
# see script aws_propagation.r
#
#######################################################################################
setawsdefaults <- function(dy,meany,qlambda,lseq,hmax,shape){
hinit <- 1
if(is.null(dy)||length(dy)==1){ 
      d<-1
      if(is.null(qlambda)) qlambda <- .96
      if(is.null(lseq)) lseq<-9
     }
if(length(dy)==2) {
     d<-2
     if(is.null(qlambda)) qlambda <-  .96
     if(is.null(lseq)) lseq<-11
	}
if(length(dy)==3){
     d<-3
     if(is.null(qlambda)) qlambda <- .96
     if(is.null(lseq)) lseq<-11
}
if(is.null(hmax)) hmax <- switch(d,250,12,5)
if(qlambda<1) lambda <- qchisq(qlambda,1) else lambda <- 1e50
hincr <- 1.25^(1/d)
#
#   set maximal bandwidth
#
# uses a maximum of about 500, 450 and 520  points, respectively.
lambda<-lambda*1.8
if(is.null(shape)) shape<-1
steps <- as.integer(log(hmax/hinit)/log(hincr))
lseq <- c(lseq,rep(1,steps))[1:steps]
if(qlambda<1) 
cat("Running PS with lambda=",signif(lambda,3)," hmax=",hmax," number of iterations=",steps,"\n")
list(lambda=lambda,lseq=lseq,hmax=hmax,d=d,shape=shape,hinit=hinit,hincr=hincr)
}
##################################################################################
#
#   AWS local constant Test for propagation condition 
#
##################################################################################
awstestprop <- function(y,tobj,mask,hakt,cpar,u,propagation){
dlw<-(2*trunc(hakt/c(1,cpar$wghts))+1)[1:cpar$d]
# all other cases
n <- cpar$n1*cpar$n2*cpar$n3
pobj <- .Fortran("cawsvar",as.double(y),
                       as.logical(mask),
                       as.integer(cpar$n1),
                       as.integer(cpar$n2),
                       as.integer(cpar$n3),
                       hakt=as.double(hakt),
                       as.double(1e40),
                       as.double(tobj$thn),
                       double(n),
                       bi=as.double(tobj$bi),
                       thn=double(n),
                       as.double(cpar$spmin),
		       double(prod(dlw)),
		       as.double(cpar$wghts),
		       PACKAGE="dti",DUP=FALSE)[c("bi","thn","hakt")]
ptheta <- pobj$thn
dim(ptheta) <- dim(y) 
narisk <- sum(abs(ptheta-u))
if(narisk==0) narisk<-1e10
propagation <- c(propagation,sum(abs(tobj$thn-ptheta))/narisk)
cat("Propagation with alpha=",max(propagation),"\n")
cat("alpha values:","\n")
print(rbind(cpar$lseq[1:length(propagation[-1])],signif(propagation[-1],3)))
propagation
}
