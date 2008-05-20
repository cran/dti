# has to be re-implemented!!!!!!!!!!!!!!!!!!!!!!!!!!1
createdata.dti <- function(file,dtensor,btb,s0,sigma,level=250){
#  btb should include zero gradients !!!
  ngrad <- dim(btb)[2]
  ddim <- dim(s0)
  dtensor[1,,,][s0<level] <- 1e-5
  dtensor[4,,,][s0<level] <- 1e-5
  dtensor[6,,,][s0<level] <- 1e-5
  dtensor[2,,,][s0<level] <- 0
  dtensor[3,,,][s0<level] <- 0
  dtensor[5,,,][s0<level] <- 0
  dim(dtensor)<-c(6,prod(ddim))
  dtensor <- t(dtensor)
  si <- exp(-dtensor%*%btb)*as.vector(s0)
  rsi <- pmax(0,rnorm(si,si,pmin(s0/2.5,sigma)))
  zz <- file(file,"wb")
  writeBin(as.integer(rsi),zz,2)
  close(zz)
  dim(si)<-c(ddim,ngrad)
  dtensor <- t(dtensor)
  dim(dtensor)<-c(6,ddim)
  list(s0=s0,si=si,dtensor=dtensor,sigma=sigma,level=level,btb=btb)
}
# END implementation needed!

#
#
#
determine.eigenvalue <- function(y, reduced=FALSE) {
  cat("\nNOTE: This code is still experimental!\n") 

  dy <- dim(y)
  n <- prod(dy[-4]) # number of voxel

  if (reduced) {
    ll <- array(0,dy[1:3])
    th <- array(0,c(dy[1:3],3))
  } else {
    ll <- array(0,c(dy[1:3],3))
    th <- array(0,c(dy[1:3],9))
  }
  ierr <- array(0,dy[1:3])

  for (i in 1:dy[1]) {
    cat(".")
    for (j in 1:dy[2]) {
      for (k in 1:dy[3]) {
        if (reduced) {
          z <- .Fortran("eigen3r",
                        as.double(y[i,j,k,]),
                        lambda = double(1),
                        theta = double(3),
                        ierr = integer(1),
                        PACKAGE="dti")[c("lambda","theta","ierr")]
          ll[i,j,k] <- z$lambda
        } else {
          z <- .Fortran("eigen3",
                        as.double(y[i,j,k,]),
                        lambda = double(3),
                        theta = double(3*3),
                        ierr = integer(1),
                        PACKAGE="dti")[c("lambda","theta","ierr")]
          ll[i,j,k,] <- z$lambda
        }
        th[i,j,k,] <- z$theta
        ierr[i,j,k] <- z$ierr
      }
    }
  }

  if (!reduced) dim(th) <- c(dy[1:3],3,3)
  
  list(lambda = ll, theta = th, ierr = ierr)
}


anisotropy <- function(eigen) {
  cat("\nNOTE: This code is still experimental!\n") 
  dimdt <- dim(eigen)
  dim(eigen) <- c(prod(dimdt[1:3]),3)

  trc <- as.vector(eigen %*% c(1,1,1))/3
  fa <- sqrt(1.5*((sweep(eigen,1,trc)^2)%*% c(1,1,1))/((eigen^2)%*% c(1,1,1)))
  ra <- sqrt(((sweep(eigen,1,trc)^2)%*% c(1,1,1))/(3*trc))

  dim(trc) <- dim(fa) <- dim(ra) <- dimdt[1:3]
  
  list(fa=fa, ra=ra, trace=trc)
}

barycentric <- function(eigen) {
  cat("\nNOTE: This code is still experimental!\n") 
  dimdt <- dim(eigen)
  dim(eigen) <- c(prod(dimdt[1:3]),3)

  trc <- as.vector(eigen %*% c(1,1,1))
  cl <- (eigen[,1] - eigen[,2]) / trc
  cp <- 2*(eigen[,2] - eigen[,3]) / trc
  cs <- 3*eigen[,3] / trc
  dim(bary$cl) <- dimdt[1:3]
  dim(bary$cp) <- dimdt[1:3]
  dim(bary$cs) <- dimdt[1:3]
  
  invisible(list(cl=cl,cp=cp,cs=cs))
}

create.dti <- function(gradient,imagefile,ddim,residuals=TRUE,level=0,part=0){
if(dim(gradient)[2]==3)  gradient<-t(gradient)
if(dim(gradient)[1]!=3)  stop("Not a valid gradient matrix")
ngrad <- dim(gradient)[2]
if(!(file.exists(imagefile))) stop("Image file does not exist")
zz<-file(imagefile,"rb")
s0 <- readBin(zz,"integer",prod(ddim),2,FALSE)
ttt <- readBin(zz,"integer",prod(ddim)*ngrad,2,FALSE)
close(zz)
dim(s0) <- ddim
mask <- s0>level
xind <- apply(mask,1,mean)
xind <- (1:ddim[1])[xind>part]
xind <- min(xind):max(xind)
cat("range in x",range(xind),"\n")
yind <- apply(mask,2,mean)
yind <- (1:ddim[2])[yind>part]
yind <- min(yind):max(yind)
cat("range in y",range(yind),"\n")
zind <- apply(mask,3,mean)
zind <- (1:ddim[3])[zind>part]
zind <- min(zind):max(zind)
cat("range in z",range(zind),"\n")
cat("Data successfully read \n")
s0 <- s0[xind,yind,zind]
mask <- mask[xind,yind,zind]
ddim0 <- ddim
ddim <- dim(s0)
dim(ttt)<-c(ddim0,ngrad)
ttt <- ttt[xind,yind,zind,]
dim(s0)<-dim(ttt)<-NULL
ttt <- -log(ttt/s0)
ttt[is.na(ttt)] <- 0
ttt[(ttt==Inf)] <- 0
ttt[(ttt==-Inf)] <- 0
n <- prod(ddim)
dim(ttt) <- c(n,ngrad)
ttt<-t(ttt)
cat("Data transformation completed \n")
btb <- matrix(0,6,ngrad)
btb[1,]<-gradient[1,]*gradient[1,]
btb[4,]<-gradient[2,]*gradient[2,]
btb[6,]<-gradient[3,]*gradient[3,]
btb[2,]<-2*gradient[1,]*gradient[2,]
btb[3,]<-2*gradient[1,]*gradient[3,]
btb[5,]<-2*gradient[2,]*gradient[3,]
btbsvd <- svd(btb)
theta <- btbsvd$u %*% diag(1/btbsvd$d) %*% t(btbsvd$v)%*% ttt
cat("Diffusion tensors generated \n")
res <- ttt - t(btb) %*% theta
rm(ttt)
gc()
mres2 <- res[1,]^2
for(i in 2:ngrad) mres2 <- mres2 + res[i,]^2
sigma2 <- array(mres2/(ngrad-6),ddim)
cat("Variance estimates generated \n")
rm(mres2)
gc()
z<-list(theta=array(theta,c(6,ddim)),sigma2=sigma2,btb=btb,scorr=c(0,0),s0=array(s0,ddim),
        ddim=ddim,ddim0=ddim0,xind=xind,yind=yind,zind=zind,mask=mask,level=level,
        ngrad=ngrad,file=imagefile,res=if(residuals) res else NULL)
class(z) <- "dti"
invisible(z)
}

getscorr <- function(dtobject){
  if(!("dti" %in% class(dtobject))) stop("Not an dti-object")
  ddim <- dtobject$ddim
  mask <- dtobject$mask
  n <- prod(ddim)
  ngrad <- dtobject$ngrad
  if(is.null(dtobject$res)) {
    imagefile <- dtobject$file
    level <- dtobject$level
    ddim0 <- dtobject$ddim0
    if(!(file.exists(imagefile))) stop("Image file does not exist")
    zz<-file(imagefile,"rb")
    s0 <- readBin(zz,"integer",prod(ddim0),2,FALSE)
    ttt <- readBin(zz,"integer",prod(ddim0)*ngrad,2,FALSE)
    close(zz)
    dim(s0) <- ddim0
    s0 <- s0[dtobject$xind,dtobject$yind,dtobject$zind]
    dim(ttt) <- c(ddim0,ngrad)
    ttt <- ttt[dtobject$xind,dtobject$yind,dtobject$zind,]
    dim(s0) <- dim(ttt) <- NULL
    ttt[is.na(ttt)] <- 0
    ttt[(ttt==Inf)] <- 0
    ttt[(ttt==-Inf)] <- 0
    n <- prod(ddim)
    dim(ttt) <- c(n,ngrad)
    ttt<-t(ttt)
    btb <- dtobject$btb
    res <- ttt - t(btb) %*% dtobject$theta
  } else {
    res <- dtobject$res
  }
  scorr <- dtobject$scorr
  dim(res) <- c(ngrad,ddim)
  res <- aperm(res,c(2:4,1))
  dim(res) <- c(n,ngrad)
  res1 <- as.vector(res[as.vector(mask),])
  scorr[1] <- mean(res1[-1]*res1[-length(res1)])/var(res1)
  cat("correlation in x-direction",signif(scorr[1],3),"\n")
  dim(res) <- c(ddim,ngrad)
  res <- aperm(res,c(2,1,3,4))
  dim(res) <- c(n,ngrad)
  res1 <- as.vector(res[as.vector(aperm(mask,c(2,1,3))),])
  scorr[2] <- mean(res1[-1]*res1[-length(res1)])/var(res1)
  cat("correlation in y-direction",signif(scorr[2],3),"\n")
  dtobject$scorr <- scorr
  dtobject$res <- NULL
  invisible(dtobject)
}

andir.image <- function(dtobject,slice=1,method=1,quant=0,minanindex=NULL,show=TRUE,...){
if(!("dti" %in% class(dtobject))) stop("Not an dti-object")
if(is.null(dtobject$anindex)) stop("No anisotropy index yet")
adimpro <- require(adimpro)
anindex <- dtobject$anindex
dimg <- dim(anindex)[1:2]
if(is.null(slice)) slice <- 1
anindex <- anindex[,,slice]
andirection <- dtobject$andirection[,,,slice]
mask <- dtobject$mask[,,slice]
anindex[anindex>1]<-0
anindex[anindex<0]<-0
dim(andirection)<-c(3,prod(dimg))
if(is.null(minanindex)) minanindex <- quantile(anindex[mask],quant)
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
andirection <- andirection*as.vector(anindex)*as.numeric(mask)*as.numeric(anindex>minanindex)
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

test <- function(mat1,mat2){
rlm<-.Fortran("rlogmap",
             as.double(mat1),
             as.double(mat2),
             as.integer(1),
             rlm=double(9),
             DUP=FALSE,
             PACKAGE="dti")$rlm
rem <- .Fortran("rexpmap",
             as.double(mat1),
             as.double(rlm),
             rem=double(9),
             integer(1),
             DUP=FALSE,
             PACKAGE="dti")$rem
print(mat1)
print(mat2)
print(matrix(rem,3,3))
}


create2.dti <- function(gradient,imagefile,ddim,xind=NULL,yind=NULL,zind=NULL){
  if (dim(gradient)[2]==3) gradient <- t(gradient)
  if (dim(gradient)[1]!=3) stop("Not a valid gradient matrix")
  ngrad <- dim(gradient)[2]

  if (!(file.exists(imagefile))) stop("Image file does not exist")
  zz <- file(imagefile,"rb")
  s0 <- readBin(zz,"integer",prod(ddim),2,FALSE)
  si <- readBin(zz,"integer",prod(ddim)*ngrad,2,FALSE)
  close(zz)
  cat("Data successfully read \n")

  if (is.null(xind)) xind <- 1:ddim[1]
  if (is.null(yind)) yind <- 1:ddim[2]
  if (is.null(zind)) zind <- 1:ddim[3]
  dim(s0) <- ddim
  s0 <- s0[xind,yind,zind]
  dim(si) <- c(ddim,ngrad)
  si <- si[xind,yind,zind,]
  ddim0 <- ddim
  ddim <- dim(s0)
  dim(s0) <- dim(si) <- NULL
  ttt <- -log(si/s0)
  ttt[is.na(ttt)] <- 0
  ttt[(ttt==Inf)] <- 0
  ttt[(ttt==-Inf)] <- 0
  n <- prod(ddim)
  dim(ttt) <- c(n,ngrad)
  ttt <- t(ttt)
  cat("Data transformation completed \n")

  btb <- matrix(0,6,ngrad)
  btb[1,] <- gradient[1,]*gradient[1,]
  btb[4,] <- gradient[2,]*gradient[2,]
  btb[6,] <- gradient[3,]*gradient[3,]
  btb[2,] <- 2*gradient[1,]*gradient[2,]
  btb[3,] <- 2*gradient[1,]*gradient[3,]
  btb[5,] <- 2*gradient[2,]*gradient[3,]
  btbsvd <- svd(btb)
  solvebtb <- btbsvd$u %*% diag(1/btbsvd$d) %*% t(btbsvd$v)
  theta <- solvebtb%*% ttt
  cat("Diffusion tensors generated \n")

  res <- ttt - t(btb) %*% theta
  mres2 <- res[1,]^2
  for(i in 2:ngrad) mres2 <- mres2 + res[i,]^2
  sigma2 <- array(mres2/(ngrad-6),ddim)
  cat("Variance estimates generated \n")

  rm(mres2)
  gc()
  z <- list(theta=array(theta,c(6,ddim)),sigma2=sigma2,btb=btb,
            solvebtb=solvebtb,scorr=c(0,0),s0=array(s0,ddim),
            ddim=ddim,ddim0=ddim0,xind=xind,yind=yind,zind=zind,
            ngrad=ngrad,file=imagefile,si=si,res=res)
  class(z) <- "dti"
  invisible(z)
}

getscorr2 <- function(dtobject,level=0){
  if(!("dti" %in% class(dtobject))) stop("Not an dti-object")

  ddim <- dtobject$ddim
  mask <- dtobject$s0>level
  n <- prod(ddim)
  ngrad <- dtobject$ngrad
  scorr <- c(0,0)
  res <- dtobject$res
  dim(res) <- c(ngrad,ddim)
  res <- aperm(res,c(2:4,1))
  dim(res) <- c(n,ngrad)
  res1 <- as.vector(res[as.vector(mask),])
  scorr[1] <- mean(res1[-1]*res1[-length(res1)])/var(res1)
  cat("correlation in x-direction",signif(scorr[1],3),"\n")
  dim(res) <- c(ddim,ngrad)
  res <- aperm(res,c(2,1,3,4))
  dim(res) <- c(n,ngrad)
  res1 <- as.vector(res[as.vector(aperm(mask,c(2,1,3))),])
  scorr[2] <- mean(res1[-1]*res1[-length(res1)])/var(res1)
  cat("correlation in y-direction",signif(scorr[2],3),"\n")
  dtobject$scorr <- scorr
  dtobject$res <- NULL
  invisible(dtobject)
}









