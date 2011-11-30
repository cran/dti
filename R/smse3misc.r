betagamma <- function(g){
dg <- dim(g)
ngrad <- if(!is.null(dg)) dg[2] else 1
z <- .Fortran("bgstats",as.double(g),
                   as.integer(ngrad),
                   double(2*ngrad),
                   bghat=double(2*ngrad*ngrad),
                   nbg=double(9*ngrad),
                   nbghat=double(9*ngrad*ngrad),
                   DUPL=FALSE,
                   PACKAGE="dti")[c("bghat","nbg","nbghat")]
dim(z$bghat) <- c(2,ngrad,ngrad)
# sphaerische Coordinaten fuer Gradienten-Paare
dim(z$nbg) <- c(3,3,ngrad)
# normalen-vektoren n1,n2,n3 Gradienten
dim(z$nbghat) <- c(3,3,ngrad,ngrad)
# normalen-vektoren n1,n2,n3 Gradienten-Paare
z
}
matrm <- function(b,g){
matrix(c(cos(b), 0, sin(b), sin(b)*sin(g), 
  cos(g), -cos(b) *sin(g), -cos(g)*sin(b), sin(g), cos(b)*cos(g)),3,3)
}
matrn4 <- function(b){
matrix(c(0,tan(b),0,-tan(b),0,-1/cos(b),0,1/cos(b),0),3,3)
}
#
#
#
getkappas <- function(grad,trace=0){
krit <- function(par,matm,m4,m5,m6){
#sum((matm-expm(par[1]*m4)%*%expm(par[2]*m5)%*%expm(par[3]*m6))^2)
.Fortran("k456krit",as.double(par),
                    as.double(matm),
                    as.double(m4),
                    as.double(m5),
                    as.double(m6),
                    erg=double(1),
                    DUPL=FALSE,
                    PACKAGE="dti")$erg
}
prta <- Sys.time()
cat("Start computing spherical distances",format(Sys.time()),"\n")
ngrad <- dim(grad)[2]
m5 <- matrix(c(0,0,1,0,0,0,-1,0,0),3,3)
m6 <- matrix(c(0,-1,0,1,0,0,0,0,0),3,3)
kappa456 <- array(0,c(3,ngrad,ngrad))
zbg <- betagamma(grad)
for(i in 1:ngrad) for(j in 1:ngrad) {
bg <- zbg$bghat[,i,j]
m4 <- matrn4(bg[1])
matm <- matrm(bg[1],bg[2])
k456 <- runif(3,-.01,.01)
z <-  optim(k456,krit,method="BFGS",matm=matm,m4=m4,m5=m5,m6=m6,control=list(trace=trace,reltol=1e-12,abstol=1e-12))
while(z$value>1e-8) {
#cat("i",i,"j",j,"value",z$value,"par",z$par,"\n")
k456 <- runif(3,-.01,.01)
z <- optim(k456,krit,method="BFGS",matm=matm,m4=m4,m5=m5,m6=m6,control=list(trace=trace,reltol=1e-12,abstol=1e-12))
#cat(" new value",z$value,"par",z$par,"\n")
}
kappa456[,i,j] <- z$par
}
while(any(abs(kappa456[2:3,,])>pi)){
kappa456[2:3,,][kappa456[2:3,,]< -pi] <- kappa456[2:3,,][kappa456[2:3,,]< -pi]+2*pi
kappa456[2:3,,][kappa456[2:3,,]> pi] <- kappa456[2:3,,][kappa456[2:3,,]> pi]-2*pi
}
prtb <- Sys.time()
cat("End computing spherical distances",format(Sys.time()),"\n")
list(k456=kappa456,bghat=zbg$bghat,nbg=zbg$nbg,nbghat=zbg$nbghat)
}





