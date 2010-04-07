selectFibers <- function(obj,  ...) cat("Selection of fibers is not implemented for this class:",class(obj),"\n")

setGeneric("selectFibers", function(obj,  ...) standardGeneric("selectFibers"))

setMethod("selectFibers","dwiFiber", function(obj, roix=NULL, roiy=NULL, roiz=NULL, nroimask=NULL, minlength=1)
{
#)
#     extract fiber information and descriptions
#
  args <- sys.call(-1)
  args <- c(obj@call,args)
  fibers <- obj@fibers
  if(minlength>1){
#
#   eliminate fibers shorter than minlength
#
  fiberstart <- obj@startind
  fiberlength <- diff(c(fiberstart,dim(fibers)[1]+1))/2
  remove <- (1:length(fiberstart))[fiberlength<minlength]
  if(length(remove)>0){
     inda <- fiberstart[remove]
     inde <- c(fiberstart,dim(fibers)[1]+1)[remove+1]-1
     fibers <- fibers[-c(mapply(":",inda,inde),recursive=TRUE),]
  }
  }
  fiberstart <- ident.fibers(fibers)
  fiberlength <- diff(c(fiberstart,dim(fibers)[1]+1))/2
  roimask <- as.integer(obj@roimask)
  mroimask <- max(roimask)
  if(mroimask>127){
     warning("Recursive use of regions of interest is limited to 7")
     roix <- NULL
     roimask <- NULL
  }
#  if((!(is.null(roix)||is.null(roiy)||is.null(roiz))&&is.null(nroimask))){
if(!(is.null(roix)&&is.null(roiy)&&is.null(roiz)&&is.null(nroimask))){
#
#    no region of interest specified otherwise
#
  if(is.null(nroimask)){ 
     nroimask <- array(0,obj@ddim)
     if(is.null(roix)) roix <- 1:obj@ddim[1]
     if(is.null(roiy)) roiy <- 1:obj@ddim[2]
     if(is.null(roiz)) roiz <- 1:obj@ddim[3]
     nroimask[roix,roiy,roiz] <- 1
  }
  lnewmask <- as.integer(log(mroimask,2))+1
  roimask[nroimask>0] <- roimask[nroimask>0]+2^lnewmask
  cat("dimension of roimask",dim(roimask),"\n")
  z <- .Fortran("roifiber",
                as.double(fibers),
                newfibers=double(prod(dim(fibers))),#new fibers
                as.integer(dim(fibers)[1]),
                integer(3*max(fiberlength)),#array for fiber in vcoord
                as.integer(max(fiberlength)),# maximum fiberlength
                as.integer(fiberstart),
                as.integer(fiberlength),
                as.integer(length(fiberstart)),#number of fibers
                as.logical(nroimask>0),#roi
                as.integer(obj@ddim[1]),
                as.integer(obj@ddim[2]),
                as.integer(obj@ddim[3]),
                as.double(obj@voxelext),
                sizenf=integer(1),
                DUP=FALSE,
                PACKAGE="dti")[c("newfibers","sizenf")]
  if(z$sizenf>1) {
     fibers <- array(z$newfibers,dim(fibers))[1:z$sizenf,]
  }  else {
     warning("No fibers found, return original")
  }
  fiberstart <- ident.fibers(fibers)
  }
  invisible(new("dwiFiber",
                call  = args,
                fibers = fibers,
                startind = as.integer(fiberstart),
                roimask = as.raw(roimask),
                gradient = obj@gradient,
                btb   = obj@btb,
                ngrad = obj@ngrad, # = dim(btb)[2]
                s0ind = obj@s0ind,
                replind = obj@replind,
                ddim  = obj@ddim,
                ddim0 = obj@ddim0,
                xind  = obj@xind,
                yind  = obj@yind,
                zind  = obj@zind,
                voxelext = obj@voxelext,
                level = obj@level,
                orientation = obj@orientation,
                source = obj@source,
                method = obj@method,
                minanindex = obj@minanindex,
                maxangle = obj@maxangle)
            )
})





ident.fibers <- function(mat){
#
#  Identify indices in mat where a new fiber starts
#
   dd <- dim(mat)
   if(dd[2]!=6){
      warning("Incorrect dimensions for fiber array")
   }
   dd <- dd[1]
   mat <- mat[,1:3]
   dim(mat) <- c(2,dd/2,3)
   dmat <- mat[2,-(dd/2),]-mat[1,-1,]
   fiberends <- (1:(dd/2-1))[apply(dmat^2,1,sum)>0]
   c(0,fiberends)*2+1
}

reduce.fibers <- function(mat){
   dd <- dim(mat)
#
#  clean up fiber description in dd
#  removes instances in dd that are either redundant or would not show up on display
#
   if(dd[2]!=6){
      warning("Incorrect dimensions for fiber array")
   }
   dd <- dd[1]
   dim(mat) <- c(2,dd/2,6)
   dmat <- mat[1,,1:3]-mat[2,,1:3]
   remove <- apply(dmat^2,1,sum)==0
   if(sum(remove)>0) warning(paste("Found ",sum(remove)," instances, where begin and end of a line segment coincide"))
   mat <- mat[,!remove,]
   dim(mat) <- c(2*dim(mat)[2],6)
   mat
}
