dwiMtImprove <- function( mtobj,dwiobj, ...) cat("No dwiMixtensor calculation defined for this class:",class(mtobj),class(dwiobj),"\n")

setGeneric("dwiMtImprove", function( mtobj,dwiobj, ...) standardGeneric("dwiMtImprove"))

setMethod("dwiMtImprove",c("dwiMixtensor","dtiData"), function(mtobj, dwiobj, maxcomp=3, method="mixtensor", reltol=1e-6, maxit=5000,ngc=1000, optmethod="BFGS", nguess=100*maxcomp^2,msc="BIC",pen=NULL,where=NULL,new=FALSE){
#
#  uses  S(g)/s_0 = w_0 exp(-l_1) +\sum_{i} w_i exp(-l_2-(l_1-l_2)(g^T d_i)^2)
#
#  
#
  set.seed(1)
  theta <- .5
  maxc <- .866
  args <- sys.call(-1)
  args <- c(mtobj@call,args)
  if(is.null(pen)) pen <- 100
  if(method=="mixtensoriso") optmethod <- "L-BFGS-B"
  ngrad <- mtobj@ngrad
  ddim <- mtobj@ddim
  mask <- mtobj@mask
  if(is.null(where)||any(dim(where)!=ddim[1:3])) where <- mask
  where <- where & mask
  s0ind <- mtobj@s0ind
  vext <- mtobj@voxelext
  if(any(ddim != dwiobj@ddim)){
     warning("incompatible Mixtensor- and dwiData objects, returning original dwiMixtensor-object")
     return(mtobj)
  }
  ns0 <- length(s0ind)
  ngrad0 <- ngrad - ns0
  if(5*(1+3*maxcomp)>ngrad0){
#     maxcomp <- max(1,trunc((ngrad0-5)/15))
     cat("Maximal number of components reduced to", maxcomp,"due to insufficient
          number of gradient directions\n")
  }
#
#  First tensor estimates to generate eigenvalues and -vectors
#
  prta <- Sys.time()
  gc()
  z <- extract(mtobj,c("order","mix","andir"))
  nmix <- mix <- z$mix
  norder <- order <- z$order
  nandir <- andir <- z$andir
  maxorder <- dim(andir)[2]
  maxcomp <- max(min(maxcomp,maxorder),max(mtobj@order))
  rm(z)
  gc()
  nsigma2 <- sigma2 <- mtobj@sigma
  norient <- orient <- mtobj@orient
  nlev <- lev <-mtobj@ev
  mask <- mtobj@mask
  cat("Start generating auxiliary objects",format(Sys.time()),"\n")
#
#  compute mean S_0, s_i/S_0 (siq), var(siq) and mask
#
  sii <- matrix(dwiobj@si,prod(ddim[1:3]),ngrad)[where,]      
  z <- .Fortran("sweepimp",# mixtens.f
                as.integer(sii[,-s0ind]),
                as.integer(sii[,s0ind]),
                as.integer(sum(where)),
                as.integer(ns0),
                as.integer(ngrad0),
                siq=double(sum(where)*ngrad0),
                s0=double(sum(where)),
                DUPL=FALSE,
                PACKAGE="dti")[c("siq","s0")]
  rm(sii)
  s0 <- array(0,ddim[1:3])
  s0[where] <- z$s0
  siq <- matrix(0,prod(ddim[1:3]),ngrad0)
  siq[where,] <- z$siq
  dim(siq) <- c(ddim[1:3],ngrad0)
  rm(z)
  gc()
  gc()
  npar0 <- 1+3*(0:maxcomp)
  npar <- if(method=="mixtensor") npar0 else c(1,2+3*(1:maxcomp))
#
#   compute penalty for model selection, default BIC
#
  penIC0 <- switch(msc,"AIC"=2*npar0/ngrad0,"BIC"=log(ngrad0)*npar0/ngrad0,
                  "AICC"=(1+npar0/ngrad0)/(1-(npar+2)/ngrad0),
                  "None"=log(ngrad0)-log(ngrad0-npar0),
                  log(ngrad0)*npar0/ngrad0)
# used to destinguish between initial estimates with qnd without isotropic part
  penIC <- switch(msc,"AIC"=2*npar/ngrad0,"BIC"=log(ngrad0)*npar/ngrad0,
                  "AICC"=(1+npar/ngrad0)/(1-(npar+2)/ngrad0),
                  "None"=log(ngrad0)-log(ngrad0-npar),
                  log(ngrad0)*npar/ngrad0)
  cat("End generating auxiliary objects",format(Sys.time()),"\n")
#
#  avoid situations where si's are larger than s0
#
  grad <- t(dwiobj@gradient[,-s0ind])
#
#   determine initial estimates for orientations 
#
  n1 <- ddim[1]
  n2 <- ddim[2]
  n3 <- ddim[3]
  igc <- 0
  ingc <- 0
  if(method=="mixtensor") where <- where & order<maxcomp
  prt0 <- Sys.time()
#
#   loop over voxel in volume
#  exclude indices at the sides of the cube (just to save index operations,
#  usually they are outside the mask anyway
#
  for(i3 in 1:n3) for(i2 in 1:n2) for(i1 in 1:n1){ # begin loop
     ordi <- order[i1,i2,i3]
     smix <- sum(mix[,i1,i2,i3])
# used to destinguish between initial estimates with qnd without isotropic part
     krit <- log(sigma2[i1,i2,i3])+ if(smix>1-1e-8) penIC0[ordi+1] else penIC[ordi+1]
     if(new) krit<-1e20
     if(where[i1,i2,i3]){ 
#   only analyze voxel within mask
     mc0 <- maxcomp
     ord <- mc0+1
     z <- .Fortran("imprparb",
                     as.integer(maxcomp),
                     as.integer(maxorder),
                     as.integer(getcube3(order,i1,i2,i3)),
                     as.double(getcube3(mix,i1,i2,i3)),
                     as.double(getcube3(andir,i1,i2,i3)),
                     as.double(getcube3(orient,i1,i2,i3)),
                     as.double(getcube3(lev,i1,i2,i3)),
                     as.double(vext),
                     param=numeric(2*maxcomp+1),
                     as.integer(2*maxcomp+1),
                     npar=integer(1),
                     DUPL=TRUE,
                     PACKAGE="dti")[c("param","npar")]   
     npar <- z$npar
     if(any(is.na(z$param))||any(is.infinite(z$param))){
         cat("voxel",i1,i2,i3,"param",z$param,"npar",npar,"\n")
         npar <- length(!is.na(cumsum(z$param)))
     }
     param <- z$param[1:npar]
     if(optmethod=="L-BFGS-B") param <- c(param,switch(method,
                                        "mixtensor"=rep(1/npar,npar),
                                        "mixtensoriso"=rep(1/npar,(npar+1))))
#
#   these are the gradient vectors corresponding to minima in spherical coordinates
#
#
#  use AIC/ngrad0, BIC/ngrad0 or AICC/ngrad0 respectively
#
#        cat("i",i1,i2,i3,"par",par[1:z$npar],"npar",z$npar,"\n")
     mc0 <- (npar-1)/2

     for(k in mc0:1){ # begin model order
        if(k<ord) {
#
#  otherwise we would reanalyze a model
#
        if(optmethod=="L-BFGS-B"){
        param <- if(method=="mixtensor") param[1:(3*k+1)] else param[1:(3*k+2)]
        param[-(1:(2*k+1))] <- if(method=="mixtensor") rep(1/k,k) else rep(1/(k+1),k+1)
        }
        if(method=="mixtensor"){
#
           z <- switch(optmethod,
                    "BFGS"=optim(param[1:(2*k+1)],mfunpl0,gmfunpl0,
                           siq=siq[i1,i2,i3,],grad=grad,pen=pen,
                           method="BFGS",control=list(maxit=maxit,reltol=reltol)),
                    "CG"=optim(param[1:(2*k+1)],mfunpl0,gmfunpl0,
                           siq=siq[i1,i2,i3,],grad=grad,pen=pen,
                           method="CG",control=list(maxit=maxit,reltol=reltol)),
                    "Nelder-Mead"=optim(param[1:(2*k+1)],mfunpl0h,
                           siq=siq[i1,i2,i3,],grad=grad,method="Nelder-Mead",
                           control=list(maxit=maxit,reltol=reltol)),
                    "L-BFGS-B"=optim(param[1:(3*k+1)],mfunpl,gmfunpl,
                           siq=siq[i1,i2,i3,],grad=grad,method="L-BFGS-B",
                           lower=c(0,rep(-Inf,2*k),rep(0,k)),
                           control=list(maxit=maxit,factr=reltol/.Machine$double.eps)))
        } else { # method=="mixtensoriso"
#           cat("voxel",i1,i2,i3,"param",param[1:(3*k+2)],"\n")
           z <- optim(param[1:(3*k+2)],mfunpli,gmfunpli,
                           siq=siq[i1,i2,i3,],grad=grad,method="L-BFGS-B",
                           lower=c(0,rep(-Inf,2*k),rep(0,k+1)),
                           control=list(maxit=maxit,factr=reltol/.Machine$double.eps))
#
#  other methods seem less numerically stable in this situation
#
        }        
        if(method=="mixtensor"){
            zz <- switch(optmethod,
                    "BFGS"=mfunplwghts0(z$par[1:(2*k+1)],siq[i1,i2,i3,],grad,pen),
                    "CG"=mfunplwghts0(z$par[1:(2*k+1)],siq[i1,i2,i3,],grad,pen),
                    "Nelder-Mead"=mfunplwghts0h(z$par[1:(2*k+1)],siq[i1,i2,i3,],grad),
                    "L-BFGS-B"=mfunwghts(z$par[1:(3*k+1)],siq[i1,i2,i3,],grad))
        } else if (method=="mixtensoriso"){
            zz <- mfunwghtsi(z$par[1:(3*k+2)],siq[i1,i2,i3,],grad)
        }
        value <- if(optmethod=="L-BFGS-B") z$value else zz$value
        ord <- zz$ord
#  replace sigmai by best variance estimate from currently best model
        if(any(zz$lev<0)||ord<k){
           ttt <- krit
#   parameters not interpretable reduce order
        } else {
           si2new <- value/ngrad0 # changed 2011/05/13 
           if(si2new<1e-15){
               cat(i1,i2,i3,ord,si2new,"\n")
               si2new <- 1e-15
               }
           ttt <- log(si2new)+penIC[1+ord]
           par <- zz$par
        }
#        cat("par",par,"value",value,"ord",ord,"w",zz$mix,"lev",zz$lev,"\n")
#
#     use directions corresponding to largest weights as initial directions
#
        if(ttt < krit & ord>0) {
           krit <- ttt
           norder[i1,i2,i3] <- ord
           nlev[,i1,i2,i3] <- zz$lev
           nmix[,i1,i2,i3] <- if(ord==maxorder) zz$mix else c(zz$mix,rep(0,maxorder-ord))
           norient[,1:ord,i1,i2,i3] <- zz$orient
           nsigma2[i1,i2,i3] <- si2new
       }
     }
   } # end model order#
    if(igc<ngc){
       igc <- igc+1
    } else {
       igc <- 1
       ingc <- ingc+1
       prt1 <- Sys.time()
       gc()
       cat("Nr. of voxel",ingc*ngc,"time elapsed:",format(difftime(prt1,prta),digits=3),"remaining time:",
            format(difftime(prt1,prt0)/(ingc*ngc)*(sum(where)-ingc*ngc),digits=3),"\n")
    }
  }# end mask
  }# end loop
  invisible(new("dwiMixtensor",
                model = "homogeneous_prolate",
                call   = args,
                ev     = nlev,
                mix    = nmix,
                orient = norient,
                order  = norder,
                p      = mtobj@p,
                th0    = mtobj@th0,
                sigma  = nsigma2,
                scorr  = mtobj@scorr, 
                bw     = mtobj@bw, 
                mask   = mtobj@mask,
                hmax   = mtobj@hmax,
                gradient = mtobj@gradient,
                btb    = mtobj@btb,
                ngrad  = mtobj@ngrad, # = dim(btb)[2]
                s0ind  = mtobj@s0ind,
                replind = mtobj@replind,
                ddim   = mtobj@ddim,
                ddim0  = mtobj@ddim0,
                xind   = mtobj@xind,
                yind   = mtobj@yind,
                zind   = mtobj@zind,
                voxelext = mtobj@voxelext,
                level = mtobj@level,
                orientation = mtobj@orientation,
                rotation = mtobj@rotation,
                source = mtobj@source,
                outlier = mtobj@outlier,
                scale = mtobj@scale,
                method = mtobj@method)
            )
   }
)
dwiMtCombine <- function(mtobj1, mtobj2, ...) cat("No dwiMixtensor calculation defined for this class:",class(mtobj1),class(mtobj2),"\n")

setGeneric("dwiMtCombine", function(mtobj1, mtobj2, ...) standardGeneric("dwiMtCombine"))

setMethod("dwiMtCombine",c("dwiMixtensor","dwiMixtensor"), function(mtobj1,mtobj2, msc="BIC", where=NULL){
#
#  combine results from two dwiMixtensor objects
#
  set.seed(1)
  args <- sys.call(-1)
  if(class(mtobj1)!="dwiMixtensor"||class(mtobj2)!="dwiMixtensor"){
     warning("First two arguments need to specify dwiMixtensor objects \n returning
     NULL")
     return(invisible(NULL))
  }
  args <- c(mtobj1@call,args)
  ngrad <- mtobj1@ngrad
  ngrad0 <- ngrad - length(mtobj1@s0ind)
  ddim <- mtobj1@ddim
  mask <- mtobj1@mask
  ngrad2 <- mtobj2@ngrad
  ddim2 <- mtobj2@ddim
  mask2 <- mtobj2@mask
  if(any(ddim!=ddim2)||any(ngrad!=ngrad2)||any(mask!=mask2)){
   warning("incompatible objects \n returning first dwiMixtensor object\n")
   return(mtobj1)
  }
  ncomp1 <- dim(mtobj1@mix)[1]
  ncomp2 <- dim(mtobj2@mix)[1]
  if(ncomp1<ncomp2){
   warning("first object should have larger maximum number of components \n 
            switching order\n")
   return(combine(mtobj2,mtobj1,msc,where))
  }
  if(is.null(where)||any(dim(where)!=ddim[1:3])) where <- mask
  where <- where & mask
  gc()
  z1 <- extract(mtobj1,c("order","mix"))
  z2 <- extract(mtobj2,c("order","mix"))
  ev1 <- mtobj1@ev
  ev2 <- mtobj2@ev
  orient1 <- mtobj1@orient
  orient2 <- mtobj2@orient
  sigma1 <- mtobj1@sigma  
  sigma2 <- mtobj2@sigma  
  npar1 <- if(mtobj1@method=="mixtensor") 1+3*(0:ncomp1) else c(1,2+3*(1:ncomp1))
  npar2 <- if(mtobj2@method=="mixtensor") 1+3*(0:ncomp1) else c(1,2+3*(1:ncomp1))
#
#   compute penalty for model selection, default BIC
#
  penIC1 <- switch(msc,"AIC"=2*npar1/ngrad0,"BIC"=log(ngrad0)*npar1/ngrad0,
                  "AICC"=(1+npar1/ngrad0)/(1-(npar1+2)/ngrad0),
                  "None"=log(ngrad0)-log(ngrad0-npar1),
                  log(ngrad0)*npar1/ngrad0)
  penIC2 <- switch(msc,"AIC"=2*npar2/ngrad0,"BIC"=log(ngrad0)*npar2/ngrad0,
                  "AICC"=(1+npar2/ngrad0)/(1-(npar2+2)/ngrad0),
                  "None"=log(ngrad0)-log(ngrad0-npar2),
                  log(ngrad0)*npar2/ngrad0)
     krit1 <- log(sigma1[where])+penIC1[z1$order[where]+1]
     krit2 <- log(sigma2[where])+penIC2[z2$order[where]+1]
     n <- prod(ddim)
     ind <- rep(FALSE,n)
     ind[where][krit1>krit2] <- TRUE
     z1$order[ind] <- z2$order[ind]
     dim(ev1) <- dim(ev2) <- c(2,n)
     ev1[,ind] <- ev2[,ind]
     dim(ev1) <- c(2,ddim) 
     sigma1[ind] <- sigma2[ind]
     dim(z1$mix) <- c(ncomp1,n)
     dim(z2$mix) <- c(ncomp2,n)
     z1$mix[1:ncomp2,ind] <- z2$mix[,ind]
     if(ncomp2<ncomp1) z1$mix[-(1:ncomp2),ind] <- 0
     dim(z1$mix) <- c(ncomp1,ddim)
     dim(orient1) <- c(2,ncomp1,n)
     dim(orient2) <- c(2,ncomp2,n)
     orient1[,1:ncomp2,ind] <- orient2[,,ind]
     dim(orient1) <- c(2,ncomp1,ddim)
     if(sum(ind)>0) cat("Improvements in ",sum(ind)," voxel \n maximal:",
     max(krit1-krit2)," median:",median((krit1-krit2)[krit1>krit2]),"\n")
  invisible(new("dwiMixtensor",
                model = "homogeneous_prolate",
                call   = args,
                ev     = ev1,
                mix    = z1$mix,
                orient = orient1,
                order  = z1$order,
                p      = mtobj1@p,
                th0    = mtobj1@th0,
                sigma  = sigma1,
                scorr  = mtobj1@scorr, 
                bw     = mtobj1@bw, 
                mask   = mtobj1@mask,
                hmax   = mtobj1@hmax,
                gradient = mtobj1@gradient,
                btb    = mtobj1@btb,
                ngrad  = mtobj1@ngrad, # = dim(btb)[2]
                s0ind  = mtobj1@s0ind,
                replind = mtobj1@replind,
                ddim   = mtobj1@ddim,
                ddim0  = mtobj1@ddim0,
                xind   = mtobj1@xind,
                yind   = mtobj1@yind,
                zind   = mtobj1@zind,
                voxelext = mtobj1@voxelext,
                level = mtobj1@level,
                orientation = mtobj1@orientation,
                rotation = mtobj1@rotation,
                source = mtobj1@source,
                outlier = mtobj1@outlier,
                scale = mtobj1@scale,
                method = mtobj1@method)
            )
   }
)
getcube3 <- function(aobj,i1,i2,i3){
ddim <- dim(aobj)
ndim <- length(ddim)
subcube <- array(0,c(switch(ndim-2,NULL,ddim[1],ddim[1:2]),3,3,3))
ind3 <- (-1):1
indz <- ind3[(i3+ind3 > 0) & (i3+ind3 <= ddim[ndim])]
indy <- ind3[(i2+ind3 > 0) & (i2+ind3 <= ddim[ndim-1])]
indx <- ind3[(i1+ind3 > 0) & (i1+ind3 <= ddim[ndim-2])]
switch(ndim-2, subcube[indx+2,indy+2,indz+2] <- aobj[i1+indx,i2+indy,i3+indz], 
               subcube[,indx+2,indy+2,indz+2] <- aobj[,i1+indx,i2+indy,i3+indz], 
               subcube[,,indx+2,indy+2,indz+2] <- aobj[,,i1+indx,i2+indy,i3+indz])
subcube
}
