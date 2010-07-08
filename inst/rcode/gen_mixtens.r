tdatamixtens <- function(w,th,alpha,beta,grad,sigma){
#  
#  w  - matrix of wghts (m+1,n1,n2,n3)
#  th - ev-parameters (2,n1,n2,n3) lambda1 = th[1]+th[2]
#                                  lambda2 = th[2]
#  alpha, beta (m,n1,n2,n3)  spherical coordinates of directions
#
args <- list(sys.call())
ddim <- dim(w)[2:4]
m <- dim(w)[1]-1
ngrad <- dim(grad)[2]
#  check array dimensions
if(any(dim(th)!=c(2,ddim))) return("wrong dimension of th")
if(any(dim(alpha)!=c(m,ddim))) return("wrong dimension of alpha")
if(any(dim(beta)!=c(m,ddim))) return("wrong dimension of beta")
if(any(dim(grad)!=c(3,ngrad))) return("wrong dimension of grad")
sb <- array(0,c(ddim,ngrad+1))
d <- array(0,c(3,m,ddim))
sal <- sin(alpha)
cal <- cos(alpha)
sbe <- sin(beta)
cbe <- cos(beta)
d[1,,,,] <- sal*cbe
d[2,,,,] <- sal*sbe
d[3,,,,] <- cal
el1 <- exp(-th[1,,,]-th[2,,,])
sb[,,,1] <- rnorm(prod(ddim),10000,10)
for(i in 1:ngrad){
   sb[,,,i+1] <- w[1,,,]*el1
   for(j in 1:m){
      dg <- d[1,j,,,]*grad[1,i]+d[2,j,,,]*grad[2,i]+d[3,j,,,]*grad[3,i]
      sb[,,,i+1] <- sb[,,,i+1]+w[j+1,,,]*exp(-th[2,,,]-th[1,,,]*dg^2)
   }
   sb[,,,i+1] <- sb[,,,i+1]*sb[,,,1]
}
sb <- sqrt(rnorm(sb,sb,sigma*sb[,,,1])^2+rnorm(sb,0,sigma*sb[,,,1])^2)
dim(sb) <- c(ddim,ngrad+1)
grad <- cbind(rep(0,3),grad)
ddim0 <- ddim
invisible(new("dtiData",
                call = args,
                si     = sb,
                gradient = grad,
                btb    = create.designmatrix.dti(grad),
                ngrad  = as.integer(ngrad+1), # = dim(btb)[2]
                s0ind  = as.integer(1), # indices of S_0 images
                replind = 1:(ngrad+1),
                ddim   = ddim,
                ddim0  = ddim0,
                xind   = 1:ddim[1],
                yind   = 1:ddim[2],
                zind   = 1:ddim[3],
                level  = 1,
                sdcoef = c(5.e-01, 0, 1, 2.e+04),
                voxelext = c(1,1,1),
                orientation = as.integer(c(0,2,5)),
                rotation = diag(3),
                source = "artificial")
            )
}
truemixtens <- function(w,th,alpha,beta,grad,sigma){
#  
#  w  - matrix of wghts (m+1,n1,n2,n3)
#  th - ev-parameters (2,n1,n2,n3) lambda1 = th[1]+th[2]
#                                  lambda2 = th[2]
#  alpha, beta (m,n1,n2,n3)  spherical coordinates of directions
#
args <- list(sys.call())
ddim <- dim(w)[2:4]
m <- dim(w)[1]-1
ngrad <- dim(grad)[2]
#  check array dimensions
if(any(dim(th)!=c(2,ddim))) return("wrong dimension of th")
if(any(dim(alpha)!=c(m,ddim))) return("wrong dimension of alpha")
if(any(dim(beta)!=c(m,ddim))) return("wrong dimension of beta")
if(any(dim(grad)!=c(3,ngrad))) return("wrong dimension of grad")
grad <- cbind(rep(0,3),grad)
ddim0 <- ddim
orient <- array(0,c(2,m,ddim))
orient[1,,,,] <- alpha
orient[2,,,,] <- beta

  invisible(new("dwiMixtensor",
                model  = "homogeneous_prolate",
                call   = args,
                ev     = th,
                mix    = w[-1,,,],
                orient = orient,
                order  = array(m,ddim),
                p      = 1,
                th0    = array(1000,ddim),
                sigma  = array(1000*sigma,ddim),
                scorr  = array(0,c(1,1,1)), 
                bw     = c(0,0,0), 
                mask   = array(TRUE,ddim),
                hmax   = 1,
                gradient = grad,
                btb    = create.designmatrix.dti(grad),
                ngrad  = as.integer(ngrad+1), # = dim(btb)[2]
                s0ind  = as.integer(1),
                replind = 1:(ngrad+1),
                ddim   = ddim,
                ddim0  = ddim0,
                xind   = 1:ddim[1],
                yind   = 1:ddim[2],
                zind   = 1:ddim[3],
                voxelext = c(1,1,1),
                level = 1,
                orientation = as.integer(c(0,2,5)),
                rotation = diag(3),
                source = "artificial",
                outlier = numeric(0),
                scale = 1,
                method = "mixtensor")
            )
}
create.designmatrix.dti <- function(gradient, bvalue=1) {
  dgrad <- dim(gradient)
  if (dgrad[2]==3) gradient <- t(gradient)
  if (dgrad[1]!=3) stop("Not a valid gradient matrix")

  btb <- matrix(0,6,dgrad[2])
  btb[1,] <- gradient[1,]*gradient[1,]
  btb[4,] <- gradient[2,]*gradient[2,]
  btb[6,] <- gradient[3,]*gradient[3,]
  btb[2,] <- 2*gradient[1,]*gradient[2,]
  btb[3,] <- 2*gradient[1,]*gradient[3,]
  btb[5,] <- 2*gradient[2,]*gradient[3,]

  btb * bvalue
}
odfdist <- function(obj1,obj2,poly=4){
  if(any(obj1@ddim!=obj2@ddim)) return(warning("incompatible dimensions"))
  if(!exists("icosa0")) data("polyeders")
  polyeder <- switch(poly+1,icosa0,icosa1,icosa2,icosa3,icosa4)
  n <- prod(obj1@ddim)
  radii1 <- .Fortran("mixtradi",
                    as.double(polyeder$vertices),
                    as.integer(polyeder$nv),
                    as.double(obj1@ev),
                    as.double(obj1@orient),
                    as.double(obj1@mix),
                    as.integer(obj1@order),
                    as.integer(dim(obj1@mix)[1]),
                    as.integer(n),
                    radii=double(n*polyeder$nv),
                    DUP=FALSE,
                    PACKAGE="dti")$radii
  dim(radii1) <- c(polyeder$nv,obj1@ddim)
  radii2 <- .Fortran("mixtradi",
                    as.double(polyeder$vertices),
                    as.integer(polyeder$nv),
                    as.double(obj2@ev),
                    as.double(obj2@orient),
                    as.double(obj2@mix),
                    as.integer(obj2@order),
                    as.integer(dim(obj2@mix)[1]),
                    as.integer(n),
                    radii=double(n*polyeder$nv),
                    DUP=FALSE,
                    PACKAGE="dti")$radii
  dim(radii2) <- c(polyeder$nv,obj1@ddim)
  apply((radii2-radii1)^2,2:4,mean)
  }
