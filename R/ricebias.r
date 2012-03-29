dwiRiceBias <- function(object,  ...) cat("No Rice Bias correction defined for this class:",class(object),"\n")

setGeneric("dwiRiceBias", function(object,  ...) standardGeneric("dwiRiceBias"))

setMethod("dwiRiceBias","dtiData",function(object, 
sigma=NULL, method="1stMoment") {
args <- object@call
if(is.null(sigma)||sigma<1){
cat("please provide a value for sigma \n Returning the original object\n")
return(object)
}
corrected <- FALSE
for(i in 1:length(args)){
if(length(grep("dwiRiceBias",args[i][[1]]))>0){
corrected <- TRUE
cat("Rice bais correction already performed by\n")
print(args[i][[1]])
cat("\n Returning the original object\n")
}
}
if(!corrected){
object@si <- array(ricebiascorr(object@si,sigma),dim(object@si))
object@call <- c(args,sys.call(-1))
}
invisible(object)
})



rrician <- function(x,m=0,s=1){
x1 <- rnorm(x,m,s)
x2 <- rnorm(x,0,s)
sqrt(x1*x1+x2*x2)
}

ricebiascorr <- function(x,s=1,method=1){
Lhalf <- function(x){
(1-x)*besselI(-x/2,0,TRUE) - x*besselI(-x/2,1,TRUE)
}
fofx <- function(x) {
      ind <- x<50
      x[ind] <- sqrt(pi/2)*Lhalf(-x[ind]^2/2)
      x
}
xt <- x/s
if(method==1){
xt <- .Fortran("rcstep",xt=as.double(xt),
               as.integer(length(xt)),
               DUPL=TRUE,
               PACKAGE="dti")$xt
} else {
y <- xt
xt <- sqrt(pmax(0,y^2-2))
}
xt*s
}

betaL <- function(L,eta){
# L-number of coils, eta - noncentrality parameter of\chi_L
require(gsl)
sqrt(pi/2)*poch(L,.5)/factorial(.5)
}
xiofetaL <- function(L,eta){
require(gsl)
2*L+eta^2-betaL(L,eta)^2*hyperg_1F1(-.5,L,-eta^2/2)^2
}
m1chiL <- function(L,eta){
require(gsl)
betaL(L,eta)*hyperg_1F1(-.5,L,-eta^2/2)
}
fixpetaL <- function(L,eta,m1,mu2,eps=1e-8,maxcount=1000){
n <- length(eta)
converged <- rep(FALSE,n)
mquot <- m1^2/mu2+1
count <- 0
while(!all(converged)&count<maxcount){
ind <- (1:n)[!converged]
etanew <- sqrt(pmax(0,xiofetaL(L,eta[ind])*mquot[ind]-2*L))
converged[ind] <- abs(etanew-eta[ind])<=eps
count <- count+1
#cat(count,"lind",length(ind),"rind",range(ind),"krit",
#     mean(abs(etanew-eta[ind])),"\n")
eta[ind] <- etanew
}
eta
}
etasolve <- function(L,m1,sigma,eps=.01){
m1s <- m1/sigma
maxeta <- max(m1s+1e-10)
x <- seq(0,maxeta,eps)
lx <- length(x)
fx <- m1chiL(L,x)
m1s[m1s<min(fx)] <- min(fx)
ind <- cut(m1s+1e-10, breaks = fx, labels=FALSE)
ind1 <- pmin(ind+1,length(fx)-1)
(x[ind1]*(m1s-fx[ind])+ x[ind]*(fx[ind1]-m1s))/(fx[ind1]-fx[ind])
}
