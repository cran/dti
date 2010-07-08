cat("Demo for mix-tensor models\n")
set.seed(1)
source(system.file("rcode/gen_mixtens.r",package="dti"))
grad <- t(read.table(system.file("dat/hardi-grad.txt",package="dti")))
n <- readline("size of data cube (n x n x n) (default n=6)") 

if (!is.null(n)) n <- as.numeric(n) else n <- 6
if(is.na(n)) n <- 6 else n <- as.integer(min(10,max(1,n)))

snr <- readline("signal-noise-ratio (default 50)") 

if (!is.null(snr)) snr <- as.numeric(snr) else snr <- 50
if(is.na(snr)) snr <- 50 else snr <- max(10,min(1000,snr))
sigma <- 1/snr

scenario <- readline("select scenario \n 1: varying directions (default)
                                      \n 2: varying mixtures
                                      \n 3: varying gfa
                                      \n 4: homogeneous")

if (is.null(scenario)) scenario <- 1 else scenario <- as.numeric(scenario)
if(is.na(scenario)) scenario <- 1 else scenario <- as.integer(max(1,min(4,scenario)))

scenario <- max(1,min(4,as.integer(scenario)))
mix <- array(0,c(4,n,n,n))
th <- array(0,c(2,n,n,n))
alpha <- array(0,c(3,n,n,n))
beta <- array(0,c(3,n,n,n))
mix[1,,,] <- 0
if(scenario==1){
th[1,,,] <- 4
th[2,,,] <- .8
alpha[1,,,] <- 0
beta[1,,,] <- 0
#for(i in 1:n) alpha[1,i,,] <- pi*(i-1)/(n-1)
#for(i in 1:n) beta[1,i,,] <- pi*(i-1)/(n-1)
#for(i in 1:n) alpha[2,,i,] <- pi/2/(n-1)*(i-1)
#for(i in 1:n) beta[2,,i,] <- 0
#for(i in 1:n) alpha[3,,,i] <- pi/2/(n-1)*(i-1)
#for(i in 1:n) beta[3,,,i] <- pi/2
# first component in x -direction
for(i in 1:n) alpha[1,i,,] <- pi/2
for(i in 1:n) beta[1,i,,] <- 0
# second component in x-y -plane
for(i in 1:n) alpha[2,,,i] <- pi/2
for(i in 1:n) beta[2,,,i] <- pi/2*(i-1)/(n-1)
# third component in x-z -plane
for(i in 1:n) alpha[3,,i,] <- pi/2*(1-(i-1)/(n-1))
for(i in 1:n) beta[3,,i,] <- 0
} else { 
angle <- readline("angle between directions (in radiant) default (and maximum) pi/2\n ") 

if (!is.null(angle)) angle <- as.numeric(angle) else angle <- pi/2
if(is.na(angle)) angle <- pi/2 else angle <- min(pi/2,max(0,angle))
# first component in x -direction
for(i in 1:n) alpha[1,i,,] <- pi/2  
for(i in 1:n) beta[1,i,,] <- 0
# second component in x-y -plane
for(i in 1:n) alpha[2,,i,] <- pi/2
for(i in 1:n) beta[2,,i,] <- angle  
# third component in x-z -plane
for(i in 1:n) alpha[3,,,i] <- pi/2-angle
for(i in 1:n) beta[3,,,i] <- 0
}
if(scenario!=2){
mix1 <- readline("first mixture coefficient (default 1/3)\n
                  0 corresponds to mixtures of order 2, 1 to order 1") 

if (!is.null(mix1)) mix1 <- as.numeric(mix1) else mix1 <- 1/3
if(is.na(mix1)) mix1 <- 1/3 else mix1 <- min(1,max(0,mix1))

mix[2,,,] <- mix1
mix[3:4,,,] <- (1-mix1)/2
} else {
for(i in 1:n) mix[2,i,,] <- i/2/(n+1) 
for(i in 1:n) mix[3,,i,] <- i/2/(n+1) 
mix[4,,,] <- 1-mix[2,,,]-mix[3,,,]
}
th[2,,,] <- .8
if(scenario!=3){
gfa <- readline("Generalized FA (default gfa = .8)") 

if (!is.null(gfa)) gfa <- as.numeric(gfa) else gfa <- .8
if(is.na(gfa)) gfa <- .8 else gfa <- min(.9,max(.2,gfa))
} else {
gfa <- seq(.5,.9,length=n^3) 
}
th[1,,,] <- .8*(gfa^2+sqrt(3*gfa^2-2*gfa^4))/(1-gfa^2)


maxcomp <- readline("maximal order of mix-tensor model (default 3)") 

if (is.null(maxcomp)) maxcomp <- 3 else maxcomp <- as.numeric(maxcomp)
if(is.na(maxcomp)) maxcomp <- 3 else maxcomp <- min(5,max(1,maxcomp))

z0 <- truemixtens(mix,th,alpha,beta,grad,sigma)
z <- tdatamixtens(mix,th,alpha,beta,grad,sigma)
zt <- dtiTensor(z)
zmix <- dwiMixtensor(z,optmethod="BFGS",maxcomp=maxcomp,reltol=1e-6)
gslexists <- "gsl" %in% .packages(TRUE)
if(gslexists) zqball <- dwiQball(z,order=8,lambda=2e-2)

size <- as.integer(min(.adimpro$xsize/3.2,.adimpro$ysize/2.4))

vodf <- readline("Visualize and compare estimated ODF's (Y/N) :")

if(toupper(vodf)!="N"){
source(system.file("rcode/mousecallbacks.r",package="dti"))
size <- 500
w1 <- show3d(z0,scale=.5,maxobjects=n^3,FOV=1,windowRect = c(1, 1, size, size))
w2 <- show3d(zt,what="odf",minalpha=1,falevel=0,scale=.5,maxobjects=n^3,FOV=1,windowRect = c(size+11, 1, 2*size+10, size))
w3 <- show3d(zmix,scale=.5,maxobjects=n^3,FOV=1,windowRect = c(1, size+11 , size, 2*size+10))
if(gslexists){ 
w4 <- show3d(zqball,scale=.5,maxobjects=n^3,FOV=1,windowRect = c(size+11, size+11, 2*size+10, 2*size+10))
mouseTrackball(dev=c(w1,w2,w3,w4))
mouseZoom(2,dev=c(w1,w2,w3,w4))
mouseFOV(3,dev=c(w1,w2,w3,w4))
cat("True ODF in device",w1,"\n")
cat("Estimated tensor ODF in device",w2,"\n")
cat("Estimated mixtensor ODF in ",w3,"\n")
cat("Estimated q-Ball ODF in ",w4,"\n")
} else {
mouseTrackball(dev=c(w1,w2,w3))
mouseZoom(2,dev=c(w1,w2,w3))
mouseFOV(3,dev=c(w1,w2,w3))
cat("True ODF in device",w1,"\n")
cat("Estimated tensor ODF in device",w2,"\n")
cat("Estimated mixtensor ODF in ",w3,"\n")
}
}

vgfa <- readline("Visualize and compare estimated FA and GFA's (Y/N) :")

if(toupper(vgfa)!="N"){

X11()
par(mfrow=c(3,n),mar=c(1,1,2,.1),mgp=c(2,1,0))
gfa0 <- extract(z0,"gfa")$gfa
fa <- extract(zt,"fa")$fa
gfa <- extract(zmix,"gfa")$gfa
for(i in 1:n) image(gfa0[,,i],col=grey((0:255)/255),zlim=c(0,1),main=paste("True gfa sl.,",i),xaxt="n",yaxt="n")
for(i in 1:n) image(fa[,,i],col=grey((0:255)/255),zlim=c(0,1),main=paste("Est. fa sl.,",i),xaxt="n",yaxt="n")
for(i in 1:n) image(gfa[,,i],col=grey((0:255)/255),zlim=c(0,1),main=paste("Est gfa sl.,",i),xaxt="n",yaxt="n")
}


