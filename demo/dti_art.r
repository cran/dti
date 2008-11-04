cat("K. Tabelow, J. Polzehl, V. Spokoiny, and H.U. Voss,\n Diffusion Tensor Imaging: Structural Adaptive Smoothing,\n Neuroimage, 39(4), 1763--1773 (2008)\n used linear tensor estimation for their examples.\n The package also contains non-linear tensor estimation which is the new default.")
a <- readline("Use non-linear (Y, default) or linear (N) tensor estimates ?")

if (toupper(a) == "N") {
  method <- "linear"
  cat("Note: Contrary to the paper above, due to numeric issues in this demo,\n there will be some additional non-positive definite tensors in the phantoms!\n")
lambda <- 47 
# this value was used in the Neuroimage paper
} else {
  method <- "nonlinear"
lambda <- 30
}
cat("---> using",method,"tensor estimation!\n")

a <- readline("Mask small non-diffusion weighted values ? (y/n)?")

if (toupper(a) == "N") {
  mins0value <- 0
} else {
  mins0value <- 100
}

a <- readline("Provide standard deviation in K-space (default: 2400, example in Neuroimage paper: 1600):")
sigma <- if(!is.null(a)) as.numeric(a) else 2400
if( is.na(sigma)) sigma <- 2400






# define some constants
rho <- 1
ddim <- c(64,64,26)
ngrad <- 25
factor <- 2.5

# read the gradient data, these are 25 gradient directions + one non-zero weighted
bvec <- read.table(system.file("dat/b-directions.txt",package="dti"))

#
#  generate files containing the phantom- and noisy diffusion weighted images
#

source(system.file("rcode/generatedata.r",package="dti"))

# Read Phantom data 

dt0obj <- dtiData(bvec,tmpfile1,mins0value=mins0value,ddim)
dt0obj <- sdpar(dt0obj,interactive=FALSE)

# Compute phantom tensors

dt0 <- dtiTensor(dt0obj, method=method)

# Compute indices of phantom 

dt0aniso <- dtiIndices(dt0)

# Read noisy data 

dtobj <- dtiData(bvec,tmpfile2,mins0value=mins0value,ddim)
dtobj <- sdpar(dtobj,interactive=FALSE)

# Estimate tensors

dthat1 <- dtiTensor(dtobj, method=method)

# Compute indices of estimated tensors 

dthat1aniso <- dtiIndices(dthat1)

# adaptive smoothing bandwidth=4

dthat4 <- dti.smooth(dtobj,hmax=4,graph=TRUE,lambda=lambda,minanindex=0,slice=15,rho=rho,lseq=NULL,method=method)

# Compute indices of estimated smoothed tensors 

dthat4aniso <- dtiIndices(dthat4)

# plot the color-coded directional maps, phantom, noisy, smoothed
par(mfrow=c(1,3))
plot(dt0aniso,slice=15)
plot(dthat1aniso,slice=15)
plot(dthat4aniso,slice=15)

# illustrate what print() does
print(dtobj)
print(dthat1)
print(dthat4)
print(dthat4aniso)

# illustrate what summary() does
summary(dtobj)
summary(dthat1)
summary(dthat4)
summary(dthat4aniso)

# write tensor to a NIFTY-file
tmpfile3 <- tempfile("dti_art")
tensor2medinria(dthat4, tmpfile3)
# read tensor from  NIFTY-file
dthat4b <- medinria2tensor(tmpfile3)
# plot the resulting object
plot(dthat4b,slice=15)

z <- readline("Visualize and compare estimated tensors (Y/N) :")

size <- as.integer(min(.adimpro$xsize/3.2,.adimpro$ysize/2.4))
if(toupper(z)!="N"){
dthat1@scale <- dt0@scale
dthat4@scale <- dt0@scale 
#  use same scale in all plots
w1<-show3d(dt0,level=.3,nz=5,center=c(20,20,13),maxobjects=2000,FOV=1,windowRect = c(1, 1, size, size))
w2<-show3d(dthat1,level=.3,nz=5,center=c(20,20,13),maxobjects=2000,FOV=1,windowRect = c(size+11, 1, 2*size+10, size))
w3<-show3d(dthat4,level=.3,nz=5,center=c(20,20,13),maxobjects=2000,FOV=1,windowRect = c(2*size+21, 1, 3*size+20, size))
source(system.file("rcode/mousecallbacks.r",package="dti"))
#
#  from package rgl::demo(mouseCallbacks)
#

mouseTrackball(dev=c(w1,w2,w3))
mouseZoom(2,dev=c(w1,w2,w3))
mouseFOV(3,dev=c(w1,w2,w3))
cat("True tensor in device",w1,"\n")
cat("Estimated tensor in device",w2,"\n")
cat("Estimated smoothed tensor in device",w3,"\n")
}
z <- readline("Visualize smoothed estimated dtiIndex (Y/N) :")

if(toupper(z)!="N"){
w4<-show3d(dt0aniso,level=.3,center=c(20,20,13),lwd=2,FOV=1,windowRect = c(1, size+21, size, 2*size+20))
w5<-show3d(dthat1aniso,level=.3,center=c(20,20,13),lwd=2,FOV=1,windowRect = c(size+11, size+21, 2*size+10, 2*size+20))
w6<-show3d(dthat4aniso,level=.3,center=c(20,20,13),lwd=2,FOV=1,windowRect = c(2*size+21, size+21, 3*size+20, 2*size+20))
mouseTrackball(dev=c(w4,w5,w6))
mouseZoom(2,dev=c(w4,w5,w6))
mouseFOV(3,dev=c(w4,w5,w6))
cat("True FA in device",w4,"\n")
cat("Estimated FA in device",w5,"\n")
cat("Estimated smoothed FA in device",w6,"\n")
}

z <- readline("End of demo, remove created objects (Y/N) :")

graphics.off()
if(toupper(z)!="N"){
file.remove(tmpfile1)
file.remove(tmpfile2)
file.remove(paste(tmpfile3,".nii",sep=""))
rm(a,btb,bvec,cphi,createdata.dti,ddim,dt0,dt0aniso,dt0obj,dthat1,dthat1aniso,
dthat4,dthat4aniso,dthat4b,dtiso,dtobj,eta,etai,etas,factor,i,ind,j,lambda,method,
mins0value,ngrad,phi,project.cylinder,rad,rad1,rad2,rho,s0,s0offa,sigma,sphi,x,y,z,tmpfile1,tmpfile2,tmpfile3,w1,w2,w3)

}


