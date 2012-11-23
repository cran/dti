betagamma <- function(g){
  dg <- dim(g)
  ngrad <- if(!is.null(dg)) dg[2] else 1
  z <- .Fortran("bgstats",
                as.double(g),
                as.integer(ngrad),
                double(2*ngrad),
                bghat = double(2*ngrad*ngrad),
                DUPL = FALSE,
                PACKAGE = "dti")[c("bghat")]
  dim(z$bghat) <- c(2, ngrad, ngrad)
  ## sphaerische Coordinaten fuer Gradienten-Paare
  z
}

matrm <- function(b, g){
  matrix(c(cos(b), 0, sin(b), sin(b)*sin(g), 
           cos(g), -cos(b)*sin(g), -cos(g)*sin(b),
           sin(g), cos(b)*cos(g)),
         3, 3)
}

matrn4 <- function(b){
#  matrix(c(0, tan(b), 0,
#           -tan(b), 0, -1/cos(b),
#           0, 1/cos(b), 0),
#         3, 3)
  matrix(c(0, sin(b), 0,
           -sin(b), 0, -1,
           0, 1, 0),
         3, 3)
}
#
#
#

getkappas <- function(grad, trace = 0, dist = 1){
#
#  dist = 1: k4^2+k5^2+|k6|
#  dist = 2: k4^2+k5^2+k6^2
#  dist = 3: sin^2(g_i,g_j)
# new version with analytic expm(par[1]*mx)
#
  krit <- function(par, matm, beta){
    ## sum((matm-expm(par[1]*m4)%*%expm(par[2]*m5)%*%expm(par[3]*m6))^2)
    .Fortran("k456krb",
             as.double(par),
             as.double(beta),
             as.double(matm),
             erg = double(1),
             DUPL = FALSE,
             PACKAGE = "dti")$erg
  }
  krit5 <- function(x,p,pk4,matm,beta){
# for line search with respect to k5 to get second solution
     p1 <- p+c(pk4/2,x,pi)
     krit(p1,matm,beta)
  }
  ngrad <- dim(grad)[2]
  if(dist<3){
  prta <- Sys.time()
  cat("Start computing spherical distances", format(Sys.time()), "\n")
  kappa456 <- kappa456a <- array(0, c(3, ngrad, ngrad))
  zbg <- betagamma(grad)
  for (i in 1:ngrad) for (j in 1:ngrad) {
    bg <- zbg$bghat[, i, j]
#   fix for discontinuity
    if(abs(cos(bg[1])) < 1.e-6) bg[1] = pi/2 - 1e-6*sign(cos(bg[1]))
    matm <- matrm(bg[1], bg[2])
    cbg1 <- cos(bg[1])
    k456 <- runif(3, -.1, .1)
    maxit <- 10000
    fnscale <- 1
    z <- optim(k456, krit, method = "BFGS", matm = matm, beta=bg[1],
               control = list(trace = trace, fnscale=fnscale, maxit=maxit, reltol = 1e-12, abstol = 1e-16))
    count <- 5
    while (z$value > 1e-14&count>0) {
      ## cat("i",i,"j",j,"value",z$value,"par",z$par,"\n")
      maxit <- maxit*1.5
      fnscale <- fnscale/2
      k456 <- runif(3, -1, 1)
      z <- optim(k456, krit, method = "BFGS", matm = matm, beta=bg[1],
                 control = list(trace = trace, fnscale=fnscale, maxit=maxit, reltol = 1e-12, abstol = 1e-16))
      ## cat(" new value",z$value,"par",z$par,"\n")
      count <- count - 1
      if(count==0)
      cat("failed for bg:",bg,"value",z$value,"par",z$par,"\n")
    }
    z$par[1] <- z$par[1]/cbg1
    kappa456[, i, j] <- z$par
    pk4 <- abs(2*pi*cbg1/(2-cbg1^2)^.5)
    if(kappa456[1,i,j] < -pk4/2) kappa456[1,i,j] <- kappa456[1,i,j] - trunc(kappa456[1, i, j]/pk4-1)*pk4
    if(kappa456[1,i,j] >= pk4/2) kappa456[1,i,j] <- kappa456[1,i,j] - trunc(kappa456[1, i, j]/pk4+1)*pk4
    if(kappa456[2,i,j] < -pi) kappa456[2,i,j] <- kappa456[2,i,j] - trunc(kappa456[2, i, j]/2/pi-1)*2*pi
    if(kappa456[2,i,j] > pi) kappa456[2,i,j] <- kappa456[2,i,j] - trunc(kappa456[2, i, j]/2/pi+1)*2*pi
    if(kappa456[3,i,j] < -pi) kappa456[3,i,j] <- kappa456[3,i,j] - trunc(kappa456[3, i, j]/2/pi-1)*2*pi
    if(kappa456[3,i,j] > pi) kappa456[3,i,j] <- kappa456[3,i,j] - trunc(kappa456[3, i, j]/2/pi+1)*2*pi
    kappa456a[,i,j] <- kappa456[,i,j]+c(pk4/2,0,pi)
    kpar <- kappa456[,i,j]*c(cbg1,1,1)
    kappa456a[2,i,j] <- optimize(krit5,c(-pi,pi),p=kpar,pk4=pk4,
                          matm=matm,beta=bg[1],maximum=FALSE)$minimum
    if(kappa456a[1,i,j] < -pk4/2) kappa456a[1,i,j] <- kappa456a[1,i,j] - trunc(kappa456a[1, i, j]/pk4-1)*pk4
    if(kappa456a[1,i,j] >= pk4/2) kappa456a[1,i,j] <- kappa456a[1,i,j] - trunc(kappa456a[1, i, j]/pk4+1)*pk4
    if(kappa456a[2,i,j] < -pi) kappa456a[2,i,j] <- kappa456a[2,i,j] - trunc(kappa456a[2, i, j]/2/pi-1)*2*pi
    if(kappa456a[2,i,j] > pi) kappa456a[2,i,j] <- kappa456a[2,i,j] - trunc(kappa456a[2, i, j]/2/pi+1)*2*pi
    if(kappa456a[3,i,j] < -pi) kappa456a[3,i,j] <- kappa456a[3,i,j] - trunc(kappa456a[3, i, j]/2/pi-1)*2*pi
    if(kappa456a[3,i,j] > pi) kappa456a[3,i,j] <- kappa456a[3,i,j] - trunc(kappa456a[3, i, j]/2/pi+1)*2*pi
    }
  dka <- switch(dist,kappa456[1,,]^2+kappa456[2,,]^2+abs(kappa456[3,,]),
                     kappa456[1,,]^2+kappa456[2,,]^2+kappa456[3,,]^2)
  dkb <- switch(dist,kappa456a[1,,]^2+kappa456a[2,,]^2+abs(kappa456a[3,,]),
                     kappa456a[1,,]^2+kappa456a[2,,]^2+kappa456a[3,,]^2)
  dim(kappa456) <- dim(kappa456a) <- c(3,ngrad*ngrad)
  kappa456[,dkb<dka] <- kappa456a[,dkb<dka]
  dim(kappa456) <- c(3,ngrad,ngrad)
  prtb <- Sys.time()
  cat("End computing spherical distances", format(Sys.time()), "\n")
  } else {
      kappa456 <- array(0, c(3, ngrad, ngrad))
      zbg <- betagamma(grad)
      for(i in 1:ngrad) kappa456[1,i,] <- 1-(grad[,i]%*%grad)^2
  }
  list(k456 = kappa456, bghat = zbg$bghat, dist=dist)
}
##
##  Correction for Rician Bias in variance estimates 
##  see file ~polzehl/R/dti/SE3smooth/adjust_lambda.r 
##  for parameter selection
##
Lhalf <- function(x) (1-x)*besselI(-x/2,0,TRUE) - x*besselI(-x/2,1,TRUE)

ksiofx <- function(x) 2 + x^2 - pi/2*Lhalf(-x^2/2)^2

improvex <- function(x, E, V, par){
  vcorr <- par[4] + par[3]*(pi/2 + atan(pmax(0, par[1]*x + par[2])^par[5]))/pi
  sqrt(pmax(0, ksiofx(x)*(1 + E^2/(V*vcorr)) - 2))
}

sigmaRicecorrected <- function(E, S, par = c(2.3547413, -2.2819829, -0.5873854, 1.5964876, 3.2652531)){
  xhat <- E/S
  xold <- 0
  count <- 100
  while(abs(xold-xhat) > 1e-10){
    xold <- xhat
    xhat <- improvex(xhat, E, S^2, par)
    count <- count - 1
    if(count < 0) break
  }
  vcorr <- par[4] + par[3]*(pi/2 + atan(pmax(0, par[1]*xhat + par[2])^par[5]))/pi
  S*sqrt(vcorr/ksiofx(xhat))
}

lvar <- function(a,mask,L,q=.9){
#
#   variance estimates for central \chi^2 with 2L df
#  Aja-Fernandez (2009) eqn. (32)
#
da <- dim(a)
ind1 <- 2:(da[1]-1)
ind2 <- 2:(da[2]-1)
ind3 <- 2:(da[3]-1)
ma <- ma2 <- array(0,da-2)
for(i in -1:1) for(j in -1:1) for(k in -1:1){
ma <- ma + a[ind1+i,ind2+j,ind3+k, drop = FALSE]
ma2 <- ma2 + a[ind1+i,ind2+j,ind3+k, drop = FALSE]^2
}
v <- 27/26*(ma2/27 - (ma/27)^2)
to <- quantile(v[mask[ind1,ind2,ind3]],q)
z <- density(v[mask[ind1,ind2,ind3]],to=to,n=1024)
s <- z$x[z$y==max(z$y)]
s2 <- 1/(2*L-2*gamma(L+.5)^2/gamma(L)^2)*s
s2
}

suggestkappa <- function(grad,vred=1,dist=1){
#
#  get a kappa value from variance reduction on the sphere
#
gstats <- getkappas(grad,dist=dist)
ngrad <- dim(grad)[2]
vred <- min(vred,ngrad-1)
d <- switch(dist,apply(gstats$k456[1:2,,]^2,2:3,sum)+abs(gstats$k456[3,,]),
                 apply(gstats$k456^2,2:3,sum),
                 apply(gstats$k456^2,2:3,sum))
kmin <- sqrt(min(d[d>1e-8]))# just to prevent from taking zero
kappa <- kmin
vredk <- 1
while(vredk < vred){
kappa <- kappa*1.005
w <- matrix(pmax(1-d/kappa^2,0),ngrad,ngrad)
vredk <- mean(apply(w,1,sum)^2/apply(w^2,1,sum))
}
list(kappa=kappa,vred=vredk)
}
vredsphere <- function(grad,kappa,dist=1){
#
#  compute initial variance reduction on the sphere 
#  for given kappa
#
gstats <- getkappas(grad,dist=dist)
ngrad <- dim(grad)[2]
d <- switch(dist,apply(gstats$k456[1:2,,]^2,2:3,sum)+abs(gstats$k456[3,,]),
                 apply(gstats$k456^2,2:3,sum),
                 apply(gstats$k456^2,2:3,sum))
w <- matrix(pmax(1-d/kappa^2,0),ngrad,ngrad)
mean(apply(w,1,sum)^2/apply(w^2,1,sum))
}


