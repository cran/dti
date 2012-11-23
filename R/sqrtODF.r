
dwiSqrtODF <- function(object,  ...) cat("No DWI Q-ball calculation defined for this class:",class(object),"\n")

setGeneric("dwiSqrtODF", function(object,  ...) standardGeneric("dwiSqrtODF"))

setMethod("dwiSqrtODF","dtiData",function(object,what="sqrtODF",order=4,forder=1,lambda=0,D0=1.4e-3){
#  D_0 set to 1.4e-3 mm^2/s
  cat("Nonnegative Definite EAP and ODF Estimation via a Unified Multi-Shell HARDI Reconstruction (Cheng et al, MICCAI 2012)\n This code is still experimental and not sufficiently testet\n")
  if(is.null(D0)) D0 <- 1.4e-3
  args <- sys.call(-1)
  args <- c(object@call,args)
  if (!(what %in% c("sqrtODF"))) {
      stop("what should specify sqrtODF\n")
          }
  ng <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  grad <- object@gradient[,-s0ind]
  sdcoef <- object@sdcoef
  z <- sioutlier(object@si,s0ind)
  si <- array(z$si,c(ng,ddim))
  index <- z$index
  rm(z)
  gc()

  # prepare data including mask
  ng <- ng - length(s0ind)
  s0 <- si[s0ind,,,]
  si <- si[-s0ind,,,]
  if (ns0>1) {
    dim(s0) <- c(ns0,prod(ddim))
    s0 <- rep(1/ns0,ns0)%*%s0
    dim(s0) <- ddim
  }
  mask <- s0 > object@level
  mask <- connect.mask(mask)
  bvalue <- object@bvalue[-s0ind]
  si <- sweep(si,2:4,s0,"/")
  nbv <- length(unique(bvalue)) # thats number of different b-balues
#  forder <- min(forder,nbv-1)
  cat("Using",forder,"radial fourier basis functions\n")
  while((order+1)*(order+2)/2*(forder+1) >=ng) order <- order-2
  cat("Using",(order+1)*(order+2)/2,"sperical harmonics\n")
  L <- Lnlm(lambda,forder,order)
  if(is.null(D0)) D0 <- getD0(siq,bvalue)  ## check ##
  n <- prod(ddim)
  dim(si) <- c(ng,n)
  nmask <- sum(mask)
  for (i in 1:length(D0)){
  t0 <- Sys.time()
  Kqzeta <- Kofgrad(grad,D0[i],forder,order,bvalue)
#
#   now we have the ingredients
#
  nk <- dim(Kqzeta)[1]
#
#  order of coefficients as in  c((order+1)*(order+2)/2,forder+1,ng)
#
  t1 <- Sys.time()
  mc.cores <- setCores(,reprt=FALSE)
  cat("Kernel computed in",format(difftime(t1,t0)),"\n")
  coefs <- matrix(.Fortran("sqrteap",
                as.double(si[,mask]),
                as.double(Kqzeta),
                as.integer(nk),
                as.integer(ng),
                as.integer(nmask),
                as.double(L),
                double(nk*mc.cores),#w
                double(nk*mc.cores),#nablam
                double(nk*mc.cores),#ck
                double(nk*mc.cores),#ck1
                coefs=double(nk*nmask),
                DUPL=TRUE,
                PACKAGE="dti")$coefs,c(nk,nmask))
  t2 <- Sys.time()
  cat("Obtained parameter estimates in",format(difftime(t2,t1)),"\n")
   z <- .Fortran("Mofcall",
                as.double(coefs),
                as.double(Kqzeta),
                as.integer(nk),
                as.integer(ng),
                as.integer(nmask),
                as.double(si[,mask]),
                as.double(L),
                fvmofc=double(nmask),
                res=double(nmask*ng),
                DUPL=TRUE,
                PACKAGE="dti")
  sigma2 <- z$fvmofc/(ng-nk)*2
  res <- array(z$res,c(ng,nmask))
  rm(z)
  cat("mean sigma2",mean(sigma2),"\n")
  if(length(D0)>1) {
     if(i==1) {
          bcoefs <- coefs
          bres <- res
          bsigma2 <- sigma2
     } else {
          ind <- rep(1,ng)%*%(res^2) < rep(1,ng)%*%(bres^2)
          bcoefs[,ind] <- coefs[,ind]
          bres[,ind] <- res[,ind]
          bsigma2[ind] <- sigma2[ind]
     }
  }
  }
  if(length(D0)>1) {
     coefs <- bcoefs
     res <- bres
     sigma2 <- bsigma2
  }  
  cat("mean sum of squared residuals",mean(rep(1,ng)%*%(res^2)),"\n")
  resa <- res
  res <- matrix(0,ng,n)
  res[,mask] <- resa
  rm(resa)
  gc()
  dim(res) <- c(ng,ddim)
  sigma2a <- sigma2
  sigma2 <- numeric(n)
  sigma2[mask] <- sigma2a
  rm(sigma2a)
  coefsa <- coefs
  coefs  <- matrix(0,nk,n)
  coefs[,mask] <- coefsa
  rm(coefsa)
  gc()
  coefs[1,!mask] <- 1
  dim(coefs) <- c((order+1)*(order+2)/2,forder+1,ddim)
  dim(sigma2) <- ddim
  cat("Variance estimates generated ",format(Sys.time()),"\n")
  th0 <- array(s0,object@ddim)
  th0[!mask] <- 0
  gc()
#
#   get spatial correlation
#
  scorr <- mcorr(res,mask,ddim,ng,lags=c(5,5,3),mc.cores=1)  
  t3 <- Sys.time()
  cat("Obtained statistics in",format(difftime(t3,t2)),"\n")
  invisible(new("dwiQball",
                call  = args,
                order = as.integer(order),
                forder = as.integer(forder),
                D0 = D0,
                lambda = lambda,
                sphcoef = coefs,
                varsphcoef = array(1,dim(coefs)),
                th0   = s0,
                sigma = sigma2,
                scorr = scorr$scorr, 
                bw = scorr$bw, 
                mask = mask,
                hmax = 1,
                gradient = object@gradient,
                bvalue = object@bvalue,
                btb   = object@btb,
                ngrad = object@ngrad, # = dim(btb)[2]
                s0ind = object@s0ind,
                replind = object@replind,
                ddim  = object@ddim,
                ddim0 = object@ddim0,
                xind  = object@xind,
                yind  = object@yind,
                zind  = object@zind,
                voxelext = object@voxelext,
                level = object@level,
                orientation = object@orientation,
                rotation = object@rotation,
                source = object@source,
                outlier = index,
                scale = 0.5,
                what = "sqrtODF")
            )
})
Qlm <- function(l1,l2,a,m1,m2,b){
    sqrt((2*l1+1)*(2*l2+1)*(2*a+1)/4/pi)*
         coupling_3j(2*l1,2*l2,2*a,0,0,0)*
         coupling_3j(2*l1,2*l2,2*a,2*m1,2*m2,2*b)
}
kappan <- function(n){
sqrt(gamma(1:(n+1))/gamma(0:n+3/2))
}
hnni <- function(n){
lni <- function(n,i) (-1)^i*gamma(n+1.5)/gamma(n-i+1)/gamma(i+1.5)/gamma(i+1)
hi <- array(0,c(n+1,n+1,2*n+1))
for(n1 in 0:n){ 
   ln1i <- lni(n1,0:n1)
   for(n2 in 0:n){
      ln2i <- lni(n2,0:n2)
      for(i in 0:n1) for(j in 0:n2){
         k <- i+j
         hi[n1+1,n2+1,k+1] <- hi[n1+1,n2+1,k+1] + ln1i[i+1]*ln2i[j+1]
         }
      }
   }
hi
}
g1f1 <- function(n,L,bvD0){
ng <- length(bvD0)
g1 <- array(0,c(L+1,2*n+1,ng))
for(i in 0:(2*n)) for(j in (0:L)) {
   g1[j+1,i+1,] <- gamma(j+i+1.5)*hyperg_1F1( i+j+1.5, 2*j+1.5, - bvD0)
}
g1
}
Inna <- function(n,L,bvD0){
ng <- length(bvD0)
Inn <- array(0,c(L+1,n+1,n+1,ng))
hi <- hnni(n)
g1 <- g1f1(n,L,bvD0)
# first component: alpha
# second component: i
# third component: bvalue*D0
kn <- kappan(n)
bvl <- t(outer(bvD0,2*(0:L),"^"))/gamma(2*(0:L)+1.5)
# first component: alpha
# second component: bvalue*D0
for(i in 0:n) for(j in 0:n){
   for(il in 1:(L+1)){
      for(k in 0:(i+j)) {
         Inn[il,i+1,j+1,] <- Inn[il,i+1,j+1,]+ hi[i+1,j+1,k+1]*g1[il,k+1,]
      }
   }
   Inn[,i+1,j+1,] <- Inn[,i+1,j+1,]*kn[i+1]*kn[j+1]*bvl
}
Inn
}
Yab <- function(order,gradients){
  order <- as.integer(max(0,order))
  if(order%%2==1){
    warning("maximum order needs to be even, increase order by one")
    order <- order+1
  } 

  # calculate spherical angles theta and phi corresponding to the gradients
  n <- dim(gradients)[2]
  theta <- phi <- numeric(n)
  for( i in 1:n){
    angles <- sphcoord(gradients[,i])
    theta[i] <- angles[1]
    phi[i] <-  angles[2]
  }

  # values of SH on specified spherical angles 
  getsphericalharmonicseven(order,theta,phi)# first index spharmon (a,b)-index 
  #second index - gradient
}
Qlmmat <- function(order){
# l in (0:(L%/%2))*2 
order <- as.integer(max(0,order))
if(order%%2==1){
warning("maximum order needs to be even, increase order by one")
order <- order+1
} 
kseq <- seq(0,order,2)
nind <- (order+1)*(order+2)/2
kseq2 <- seq(0,2*order,2)
nind2 <- (2*order+1)*(order+1)
Qmat <- array(0,c(nind,nind,nind2))
i1 <- 1
for(l1 in kseq){
   m1seq <- seq(-l1,l1,1)
   for(m1 in m1seq){
      i2 <- 1
      for(l2 in kseq){
         m2seq <- seq(-l2,l2,1)      
         for(m2 in m2seq){
            i3 <- 1
            for(a in kseq2){
               if(a<=l1+l2&l1<=l2+a&l2<=a+l1){
                  bseq <- seq(-a,a,1)
                  for(b in bseq){
                     if(m1+m2+b==0) Qmat[i1,i2,i3] <- Qlm(l1,l2,a,m1,m2,b)
#  otherwise the result is zero
                     i3 <- i3+1
                     }
                  } else {
                    i3 <- i3+2*a+1
                  }
               }
               i2 <- i2+1
            }
         }
         i1 <- i1+1
      }
   }
Qmat
}
QlmmatR <- function(order){
# warning
# Spherical integral over triplets of real valued SH basis functions
Qmat <- Qlmmat(order)
# thats the array of corresponding integrals for the complex valued basis
indlm <- function(l,m) (l^2+l+2)/2+m
indl <- function(i,horder){
rep(2*(0:horder),4*(0:horder)+1)[i]
}
indm <- function(i,horder){
il <- indl(i,horder)
i - (il^2+il+2)/2
}
dqm <- dim(Qmat)
QmatR <- array(0,dqm)
for(i1 in 1:dqm[1]){
   l1 <- indl(i1,order/2); m1 <- indm(i1,order/2)
   for(i2 in 1:dqm[2]){
      l2 <- indl(i2,order/2); m2 <- indm(i2,order/2)
      for(i3 in 1:dqm[3]){
         l3 <- indl(i3,order); m3 <- indm(i3,order)
         result <- 0
         if (m1<0 && m2<0 && m3<0){
            if (m1+m2==m3) result <- (-1)^m3*Qmat[indlm(l1,m1),indlm(l2,m2),indlm(l3,-m3)]/sqrt(2)#
            else if (m1+m3==m2) result <- (-1)^m2*Qmat[indlm(l1,m1),indlm(l2,-m2),indlm(l3,m3)]/sqrt(2)#
            else if (m2+m3==m1) result <- (-1)^m1*Qmat[indlm(l1,-m1),indlm(l2,m2),indlm(l3,m3)]/sqrt(2)#
         }
#         if (m1>0 && m2<0 && m3<0 || m1<0 && m2>0 && m3<0 || m1<0 && m2<0 && m3>0)      result = 0
         if (m1>0 && m2>0 && m3<0){ 
            if (m1+m2+m3==0) result <- (-1)^m3*Qmat[indlm(l1,m1),indlm(l2,m2),indlm(l3,m3)]/sqrt(2)#
            else if (m1==m2+m3) result <- (-1)^m2*Qmat[indlm(l1,m1),indlm(l2,-m2),indlm(l3,-m3)]/sqrt(2)#
            else if (m2==m1+m3) result <- (-1)^m1*Qmat[indlm(l1,-m1),indlm(l2,m2),indlm(l3,-m3)]/sqrt(2)#
         }
         if (m1>0 && m2<0 && m3>0){
            if (m1+m2+m3==0) result <- (-1)^m2*Qmat[indlm(l1,m1),indlm(l2,m2),indlm(l3,m3)]/sqrt(2)#
            else if (m1==m2+m3) result <- (-1)^m3*Qmat[indlm(l1,m1),indlm(l2,-m2),indlm(l3,-m3)]/sqrt(2)#
            else if (m3==m1+m2) result <- (-1)^m1*Qmat[indlm(l1,-m1),indlm(l2,-m2),indlm(l3,m3)]/sqrt(2)#
         }
         if (m1<0 && m2>0 && m3>0){
            if (m1+m2+m3==0) result <- (-1)^m1*Qmat[indlm(l1,m1),indlm(l2,m2),indlm(l3,m3)]/sqrt(2)#
            else if (m3==m2+m1) result <- (-1)^m2*Qmat[indlm(l1,-m1),indlm(l2,-m2),indlm(l3,m3)]/sqrt(2)#
            else if (m2==m1+m3) result <- (-1)^m3*Qmat[indlm(l1,-m1),indlm(l2,m2),indlm(l3,-m3)]/sqrt(2)#
         }
#         if (m1>0 && m2>0 && m3>0)     result <- 0
         if (m1==0 && m2==0 && m3==0) result <- Qmat[indlm(l1,m1),indlm(l2,m2),indlm(l3,m3)]
         if (m1==0 && m2==m3) result <-(-1)^m3*Qmat[indlm(l1,m1),indlm(l2,m2),indlm(l3,-m3)]
         if (m2==0 && m3==m1) result <-(-1)^m1*Qmat[indlm(l1,m1),indlm(l2,m2),indlm(l3,-m3)]
         if (m3==0 && m1==m2) result <-(-1)^m2*Qmat[indlm(l1,m1),indlm(l2,-m2),indlm(l3,m3)]
         QmatR[i1,i2,i3] <- result
         }
      }  
   }
QmatR
}

Kofgrad <- function(grad,D0,N=1,L=2,bvalue){
#
#   compute the kernel K(q|\zeta) for all gradient directions
#   according to Jian Cheng "Nonnegative Definite EAP and ODF Estimation ...
# 
#  this differs from Jian Chengs paper using b D_0 instead of \pi^2 \zeta q^2 
#  which simpifies the expression for I_{nn'\alpha}
#
require(gsl)
ngrad <- dim(grad)[2]
Np1 <- N+1
Lind <- (L+1)*(L+2)/2
Lind2 <- (2*L+1)*(L+1)
a <- NULL
for(i in 0:(L)) a <- c(a,rep(2*i,4*i+1))
Kofg <- array(0,c(Lind,N+1,Lind,N+1,ngrad))
Yabm <- Yab(2*L,grad)# dim Lind2, ngrad
Qlmm <- QlmmatR(L)# dim Lind,Lind,Lind2
Inn <- Inna(N,L,bvalue*D0)# dim L+1 (even values only) , N+1, N+1
for(i1 in 1:Lind) for(i2 in 1:i1)
for(n1 in 1:Np1) for(n2 in 1:n1){
z <- rep(0,ngrad)
ind <- (1:Lind2)[Qlmm[i1,i2,]!=0]
for(ab in ind) {
z <- z+4*pi*(-1)^(a[ab]/2)*Inn[a[ab]/2+1,n1,n2,]*Qlmm[i1,i2,ab]*Yabm[ab,]
}
Kofg[i1,n1,i2,n2,] <- z
# use symmetry
Kofg[i2,n1,i1,n2,] <- z
Kofg[i2,n2,i1,n1,] <- z
Kofg[i1,n2,i2,n1,] <- z
}
dim(Kofg) <- c(Lind*(N+1),Lind*(N+1),ngrad)
#save(bvalue,D0,Lind,Np1,ngrad,Lind2,a,Yabm,Qlmm,Inn,Kofg,file="YQI.rsc")
Kofg
}
Lnlm <- function(lambda,N,L){
if(length(lambda)!=2) lambda <- rep(lambda[1],2)
lord <- rep(seq(0,L,2),2*seq(0,L,2)+1)
lord <- lambda[2]*(lord*(lord+1))^2
nord <- lambda[1]*((0:N)*(1:(N+1)))^2     
Lnlm <- array(rep(lord,N+1),c((L+1)*(L+2)/2,N+1))
for(i in 1:(N+1)) Lnlm[,i] <- Lnlm[,i]+nord[i]
Lnlm
}
Mofcall <- function(coef,kern,si,lambda=0){
n <- prod(dim(coef)[-(1:2)])
nk <- prod(dim(coef)[1:2])
forder <- dim(coef)[2]-1
if(any(dim(kern)[1:2]!=c(nk,nk))) return(0)
ng <- dim(kern)[3]
dim(si) <- c(n,ng)
order <- switch(as.character(nk/(forder+1)),"1"=0,"6"=2,"15"=4,"28"=6,"45"=8,"66"=10)
L <- Lnlm(lambda,forder,order)
z<-.Fortran("Mofcall",as.double(coef),
                   as.double(kern),
                   as.integer(nk),
                   as.integer(ng),
                   as.integer(n),
                   as.double(t(si)),
                   as.double(L),
                   krit=double(n),
                   res=double(ng*n),
                   DUPL=TRUE,
                   PACKAGE="dti")[c("krit","res")]
list(Mofc=array(z$krit,dim(coef)[-(1:2)]),residuals=array(z$res,c(ng,dim(coef)[-(1:2)])),
     fitted.values=array(z$res+t(si),c(ng,dim(coef)[-(1:2)])))
}
                   