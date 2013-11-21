# This file contains the implementation of dti.smooth() for 
# "dtiData" Adaptive smoothing in SE(3) considering b=0 as an individual shell

dwi.smooth.ms.wghts <- function(object,kstar,lambda=12,kappa0=.5,ncoils=1,sigma=NULL,ws0=1,level=NULL,vx=NULL,vy=NULL,vz=NULL,verbose=FALSE,usemaxni=TRUE){
#
#   vx,vy,vz  voxel-coordinates
#
  args <- sys.call(-1)
  args <- c(object@call,args)
  sdcoef <- object@sdcoef
  level <- object@level
  vext <- object@voxelext[2:3]/object@voxelext[1]
  varstats <- sofmchi(ncoils)
  if(length(sigma)==1) {
     cat("using supplied sigma",sigma,"\n")
  } else {
  mask <- getmask(object,level)$mask
  sigma <- numeric(object@ngrad)
  for(i in 1:object@ngrad){
  sigma[i] <- awssigmc(object@si[,,,i],12,mask,ncoils,vext,h0=1.25,verbose=verbose)$sigma
  cat("image ",i," estimated sigma",sigma[i],"\n")
  }
 cat("quantiles of estimated sigma values",quantile(sigma),"\n")
  sigma <- median(sigma)
 cat("using median estimated sigma",sigma,"\n")
  }
#
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  ngrad <- ngrad - ns0
  grad <- object@gradient[,-s0ind]
  bvalues <- object@bvalue[-s0ind]
  msstructure <- getnext3g(grad,bvalues)
  nshell <- msstructure$nbv
  sb <- object@si[,,,-s0ind]
  s0 <- object@si[,,,s0ind]
  if(is.null(kappa0)){
#  select kappa based on variance reduction on the sphere
     warning("You need to specify  kappa0  returning unsmoothed object")
     return(object)
   }
#
#  rescale so that we have Chi-distributed values
#
  sb <- sb/sigma
  s0 <- s0/sigma
  if(ns0>1){
     dim(s0) <- c(prod(ddim),ns0)
     s0 <- s0%*%rep(1/sqrt(ns0),ns0)
#  make sure linear combination of s0 has same variance as original 
     dim(s0) <- ddim
  }
  mask <- s0>(sqrt(ns0)*level/sigma)
##
##  check memory size needed for largest vector
##
  vsize <- (nshell+1)*ngrad*prod(ddim)
  if(vsize>2^31-1){
     cat("region specified is to large, a vector of size",vsize,"needs to be passed through
      the .Fortran, this is limited to 2^31-1 in R\n reducing the zind")
     maxzind <- (2^31-1)/(nshell+1)/ngrad/prod(ddim[1:2])
     nzind <- apply(mask,3,sum)
     
  }
  gc()
  ddim0 <- dim(s0)
  ddimb <- dim(sb)
  gradstats <- getkappasmsh3(grad, msstructure)
  hseq <- gethseqfullse3msh(kstar,gradstats,kappa0,vext=vext)
  nind <- as.integer(hseq$n*1.25)
# make it nonrestrictive for the first step
  z <- list(th=array(1,ddimb), th0=array(1,ddim0), ni = array(1,ddimb), ni0 = array(1,ddim0))
  if(usemaxni){
     ni <- array(1,ddimb)
     ni0 <- array(1,ddim0)
     }
  prt0 <- Sys.time()
  cat("adaptive smoothing in SE3, kstar=",kstar,if(verbose)"\n" else " ")
  kinit <- if(lambda<1e10) 0 else kstar
  mc.cores <- setCores(,reprt=FALSE)
  gc()
  pwghts <- list(NULL)
  k1 <- 1
  for(k in kinit:kstar){
     hakt <- hseq$h[,k+1]
     t0 <- Sys.time()
     thnimsh <- interpolatesphere1(z$th,z$th0,z$ni,z$ni0,msstructure,mask)
     gc()
     t1 <- Sys.time()
     param <- lkfullse3msh(hakt,kappa0/hakt,gradstats,vext,nind) 
     hakt0 <- mean(hakt)
     param0 <- lkfulls0(hakt0,vext,nind) 
     vs2 <- varstats$s2[findInterval(thnimsh$mstheta, varstats$mu, all.inside = TRUE)]/2
     vs02 <- varstats$s2[findInterval(thnimsh$msth0, varstats$mu, all.inside = TRUE)]/2
     pind1 <- cbind(param$ind,-param$ind[,param$ind[1,]>0])
     pind1[4:5,] <- abs(pind1[4:5,])
     n1 <- dim(pind1)[2]
     pind01 <- cbind(param0$ind,-param$ind0[,param$ind0[1,]>0])
     n01 <- dim(pind01)[2]
     t2 <- Sys.time()
     z <- .Fortran("adsmse3w",
                as.double(sb),#y
                as.double(s0),#y0
                as.double(thnimsh$mstheta),#th
                as.double(thnimsh$msni),#ni/si^2
                as.double(thnimsh$msth0),#th0
                as.double(thnimsh$msni0),#ni0/si^2
                as.double(vs2),#var/2
                as.double(vs02),#var/2 for s0
                as.logical(mask),#mask
                as.integer(nshell+1),#ns number of shells
                as.integer(ddim0[1]),#n1
                as.integer(ddim0[2]),#n2
                as.integer(ddim0[3]),#n3
                as.integer(ngrad),#ngrad
                as.double(lambda),#lambda
                as.double(ws0),# wghts0 rel. weight for s0 image
                as.integer(mc.cores),#ncores
                as.integer(param$ind),#ind
                as.double(param$w),#w
                as.integer(param$n),#n
                as.integer(param0$ind),#ind0
                as.double(param0$w),#w0
                as.integer(param0$n),#n0
                th=double(prod(ddimb)),#thn
                ni=double(prod(ddimb)),#nin
                th0=double(prod(ddim0)),#th0n
                ni0=double(prod(ddim0)),#ni0n
                double(ngrad*mc.cores),#sw
                double(ngrad*mc.cores),#swy
                double((nshell+1)*mc.cores),#thi
                double((nshell+1)*mc.cores),#nii
                double((nshell+1)*mc.cores),#fsi2  
                as.integer(vx),
                as,integer(vy),
                as.integer(vz),
                as.integer(n1),
                as.integer(n01),
                w=double(n1),
                w0=double(n01),
                DUPL=FALSE,
                PACKAGE="dti")[c("ni","th","ni0","th0","w","w0")]
    t3 <- Sys.time()
    pwghts[[k1]] <- list(ind=pind1,w=z$w,n=n1,ind0=pind01,w0=z$w0,n0=n01)
    k1 <- k1+1
    gc()
    if(usemaxni){
       ni <- z$ni <- if(usemaxni) pmax(ni,z$ni)
       ni0 <- z$ni0 <- if(usemaxni) pmax(ni0,z$ni0)
    }
    dim(z$th) <- ddimb
    dim(z$th0) <- ddim0
    dim(z$ni) <- ddimb
    dim(z$ni0) <- ddim0
    if(verbose){
      dim(z$ni) <- c(prod(ddim0),ngrad)
      cat("k:",k,"h_k:",signif(max(hakt),3)," quartiles of ni",signif(quantile(z$ni[mask,]),3),
          "mean of ni",signif(mean(z$ni[mask,]),3),"\n              quartiles of ni0",signif(quantile(z$ni0[mask]),3),
          "mean of ni0",signif(mean(z$ni0[mask]),3),
          " time elapsed:",format(difftime(Sys.time(),prt0),digits=3),"\n")
      cat("interpolation:",format(difftime(t1,t0),digits=3),
          "param:",format(difftime(t2,t1),digits=3),
          "smoothing:",format(difftime(t3,t2),digits=3),"\n")
    } else {
      cat(".")
    }
  }
  ngrad <- ngrad+1
#
#  one s0 image only
#
  si <- array(object@si,c(ddim,ngrad))
#
#  back to original scale
#
  si[,,,1] <-  z$th0/sqrt(ns0)*sigma
#  go back to original s0 scale
  si[,,,-1] <- z$th*sigma
  lobject <- list(si=si,gradient=grad,bvalues=bvalues,pwghts=pwghts)
  lobject
}



