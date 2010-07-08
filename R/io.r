################################################################
#                                                              #
# Section for data reading functions                           #
#                                                              #
################################################################

dtiData <- function(gradient,imagefile,ddim,xind=NULL,yind=NULL,zind=NULL,level=0,mins0value=0,maxvalue=10000,voxelext=c(1,1,1),orientation=c(0,2,5),rotation=diag(3)) {
  args <- list(sys.call())
  if (any(sort((orientation)%/%2) != 0:2)) stop("invalid orientation \n")
  if (dim(gradient)[2]==3) gradient <- t(gradient)
  if (dim(gradient)[1]!=3) stop("Not a valid gradient matrix")
  ngrad <- dim(gradient)[2]
  s0ind <- (1:ngrad)[apply(abs(gradient),2,max)==0] 
  if (!(file.exists(imagefile))) stop("Image file does not exist")
  cat("Start Data reading",date(), "\n")
  zz <- file(imagefile,"rb")
#  si now contains all images (S_0 and S_I), ngrad includes 
#  number of zero gradients

  if (is.null(xind)) xind <- 1:ddim[1]
  if (is.null(yind)) yind <- 1:ddim[2]
  if (is.null(zind)) zind <- 1:ddim[3]
  si <- numeric()
  for (grad in 1:ngrad) {
    sitemp <- readBin(zz,"integer",prod(ddim),2,FALSE)
    dim(sitemp) <- ddim
    si <- c(si,sitemp[xind,yind,zind])
    cat(".")
  }
  close(zz)
  dim(si) <- c(length(xind),length(yind),length(zind),ngrad)
  dimsi <- dim(si)

  cat("Data successfully read",date(), "\n")

#
#   set correct orientation
#
  xyz <- (orientation)%/%2+1
  swap <- orientation%%2
  if(any(xyz!=1:3)) {
      abc <- 1:3
      abc[xyz] <- abc
      si <- aperm(si,c(abc,4))
      swap[xyz] <- swap
      voxelext[xyz] <- voxelext
      dimsi[xyz] <- dimsi[1:3]
      ddim[xyz] <- ddim[1:3]
      gradient[xyz,] <- gradient
  }
  if(swap[1]==1) {
      si <- si[dimsi[1]:1,,,] 
      gradient[1,] <- -gradient[1,]
      }
  if(swap[2]==1) {
      si <- si[,dimsi[2]:1,,]  
      gradient[2,] <- -gradient[2,]
      }
  if(swap[3]==0) {
      si <- si[,,dimsi[3]:1,]    
      gradient[3,] <- -gradient[3,]
      }
#
#   orientation set to radiological convention
#
  si <- .Fortran("initdata",
                 si=as.integer(si),
                 as.integer(dimsi[1]),
                 as.integer(dimsi[2]),
                 as.integer(dimsi[3]),
                 as.integer(dimsi[4]),
                 as.integer(maxvalue),
                 PACKAGE="dti")$si
#  this replaces the content off all voxel with elements <=0 or >maxvalue by 0
     dim(si) <- dimsi
  level <- max(mins0value,level*mean(si[,,,s0ind][si[,,,s0ind]>0])) # set level to level*mean  of positive s_0 values
  ddim0 <- as.integer(ddim)
  ddim <- as.integer(dim(si)[1:3])

  cat("Create auxiliary statistics",date(), " \n")
  rind <- replind(gradient)
  
  invisible(new("dtiData",
                call = args,
                si     = si,
                gradient = gradient,
                btb    = create.designmatrix.dti(gradient),
                ngrad  = ngrad, # = dim(btb)[2]
                s0ind  = s0ind, # indices of S_0 images
                replind = rind,
                ddim   = ddim,
                ddim0  = ddim0,
                xind   = xind,
                yind   = yind,
                zind   = zind,
                level  = level,
                sdcoef = rep(0,4),
                voxelext = voxelext,
                orientation = as.integer(c(0,2,5)), #   orientation set to radiological convention
                rotation = rotation,
                source = imagefile)
            )
}

############

readDWIdata <- function(gradient, dirlist, format, nslice = NULL, order = NULL,
                        xind=NULL, yind=NULL, zind=NULL,
                        level=0, mins0value=0, maxvalue=10000,
                        voxelext=NULL, orientation=c(0,2,5), rotation=diag(3)) {
  # basic consistency checks
  args <- list(sys.call())
  if (!(format %in% c("DICOM","NIFTI","ANALYZE","AFNI")))
    stop("Cannot handle other formats then DICOM|NIFTI|ANALYZE|AFNI, found:",format)
  if ((format == "DICOM") & is.null(nslice))
    stop("Cannot handle DICOM folders without specifying number of slices nslice!")
#  if (any(sort((orientation)%/%2) != 0:2)) stop("invalid orientation \n")
  if (dim(gradient)[2]==3) gradient <- t(gradient)
  if (dim(gradient)[1]!=3) stop("Not a valid gradient matrix")
  ngrad <- dim(gradient)[2]
  s0ind <- (1:ngrad)[apply(abs(gradient),2,max)==0] 

  # generate file list in specified order
  filelist <- NULL
  for (dd in dirlist) filelist <- c(filelist, paste(dd,list.files(dd),sep=.Platform$file.sep))
  if (format == "DICOM") {
    if (is.null(zind)) zind <- 1:nslice
    if (length(filelist) != ngrad * nslice)
      stop("Number of found files does not match ngrad*nslice",length(filelist))
    if (is.null(order)) {
      order <- 1:(ngrad*nslice)
    } else {
      if (length(order) != ngrad*nslice)
        stop("Length of order vector does not match ngrad*nslice")
    }
    dim(order) <- c(nslice,ngrad)
    order <- order[zind,]
    dim(order) <- NULL
    filelist <- filelist[order]
  } else {
    if (format =="ANALYZE") filelist <- unlist(strsplit(filelist[regexpr("\\.hdr$", filelist) != -1],"\\.hdr"))
    if (format =="AFNI") filelist <- filelist[regexpr("\\.HEAD$", filelist) != -1]
    if (length(filelist) != ngrad)
      stop("Number of found files does not match ngrad",length(filelist),"\nPlease provide each gradient cube in a separate file.")
    if (is.null(order)) {
      order <- 1:ngrad
    } else {
      if (length(order) != ngrad)
        stop("Length of order vector does not match ngrad")
    }
    filelist <- filelist[order]
  }
  # read all DICOM files
  cat("Start reading data",date(), "\n")
  si <- numeric()
  cat("\n")
  ddim <- NULL
  first <- TRUE
  i <- 0
  for (ff in filelist) {
    i <- i+1
    cat(".")
    if (format == "DICOM") {
      data <- read.DICOM(ff)
    } else if (format == "NIFTI") {
      data <- fmri::read.NIFTI(ff,setmask=FALSE)
      nslice <- data$dim[3]
      if (is.null(zind)) zind <- 1:nslice
    } else if (format == "ANALYZE") {
      data <- fmri::read.ANALYZE(ff,setmask=FALSE)
      nslice <- data$dim[3]
      if (is.null(zind)) zind <- 1:nslice
    } else if (format == "AFNI") {
      data <- fmri::read.AFNI(ff,setmask=FALSE)
      nslice <- data$dim[3]
      if (is.null(zind)) zind <- 1:nslice
    } 
    if (is.null(ddim)) ddim <- c(data$dim[1:2],nslice,ngrad)
    if (is.null(voxelext)) {
      if (!is.null(data$delta)) {
        if (!prod(voxelext == data$delta))
          warning("Voxel extension",voxelext,"is not found in data:",data$delta)
        voxelext <- data$delta
      } else {
        warning("Voxel extension neither found nor given!")
      }
    }
    if (is.null(xind)) xind <- 1:data$dim[1]
    if (is.null(yind)) yind <- 1:data$dim[2]
    if (format == "DICOM") {
      if(first){ 
         ttt <- extract.data(data)[xind,yind]
         nttt <- dim(ttt)
         n <- length(filelist)
         si <- numeric(n*prod(nttt))
         dim(si) <- c(nttt,n)
         si[,,1]<- ttt
         first <- FALSE
      } else {
         si[,,i] <- extract.data(data)[xind,yind]
      }
    } else {
      if(first){ 
         ttt <- extract.data(data)[xind,yind,zind,]
         nttt <- dim(ttt)
         n <- length(filelist)
         si <- numeric(n*prod(nttt))
         dim(si) <- c(nttt,n)
         if(length(nttt)==4) si[,,,,1]<- ttt else si[,,,1] <- ttt
         first <- FALSE
     } else {
      ttt <- extract.data(data)[xind,yind,zind,]
      if(length(nttt)==4) si[,,,,i] <- ttt else si[,,,i] <- ttt
      }
    }
  }
  cat("\n")
  dim(si) <- c(length(xind),length(yind),length(zind),ngrad)
  dimsi <- dim(si)
  cat("Data successfully read",date(), "\n")

  # redefine orientation
  xyz <- (orientation)%/%2+1
  swap <- orientation%%2
  if(any(xyz!=1:3)) {
      abc <- 1:3
      abc[xyz] <- abc
      si <- aperm(si,c(abc,4))
      swap[xyz] <- swap
      voxelext[xyz] <- voxelext
      dimsi[xyz] <- dimsi[1:3]
      ddim[xyz] <- ddim[1:3]
      gradient[xyz,] <- gradient
  }
  if(swap[1]==1) {
      si <- si[dimsi[1]:1,,,] 
      gradient[1,] <- -gradient[1,]
      }
  if(swap[2]==1) {
      si <- si[,dimsi[2]:1,,]  
      gradient[2,] <- -gradient[2,]
      }
  if(swap[3]==0) {
      si <- si[,,dimsi[3]:1,]    
      gradient[3,] <- -gradient[3,]
      }
  # orientation set to radiological convention
  si <- .Fortran("initdata",
                 si=as.integer(si),
                 as.integer(dimsi[1]),
                 as.integer(dimsi[2]),
                 as.integer(dimsi[3]),
                 as.integer(dimsi[4]),
                 as.integer(maxvalue),
                 PACKAGE="dti")$si
  # this replaces the content off all voxel with elements <=0 or >maxvalue by 0
  dim(si) <- dimsi
  level <- max(mins0value,level*mean(si[,,,s0ind][si[,,,s0ind]>0])) # set level to level*mean  of positive s_0 values
  ddim0 <- as.integer(ddim)
  ddim <- as.integer(dim(si)[1:3])

  cat("Create auxiliary statistics",date(), " \n")
  rind <- replind(gradient)

  invisible(new("dtiData",
                call = args,
                si     = si,
                gradient = gradient,
                btb    = create.designmatrix.dti(gradient),
                ngrad  = ngrad, # = dim(btb)[2]
                s0ind  = s0ind, # indices of S_0 images
                replind = rind,
                ddim   = ddim,
                ddim0  = ddim0,
                xind   = xind,
                yind   = yind,
                zind   = zind,
                level  = level,
                sdcoef = rep(0,4),
                voxelext = voxelext,
                orientation = as.integer(c(0,2,5)),
                rotation = rotation,
                source = paste(dirlist,collapse="|"))
            )
}

################################################################
#                                                              #
# Section for summary(), print(), show() functions (generic)   #
#                                                              #
################################################################

setMethod("print", "dtiData",
function(x){
    cat("  Object of class", class(x),"\n")
    cat("  Generated by calls    :\n")
    print(x@call)
    cat("  Dimension            :", paste(x@ddim, collapse="x"), "\n")
    cat("  Number of Gradients  :", paste(x@ngrad, collapse="x"), "\n")
    cat("  Source-Filename      :", x@source, "\n")
    cat("  Slots                :\n")
    print(slotNames(x))
    invisible(NULL)
})
setMethod("print", "dtiTensor",
function(x){
    cat("  Object of class", class(x),"\n")
    cat("  Generated by calls    :\n")
    print(x@call)
    cat("  Dimension            :", paste(x@ddim, collapse="x"), "\n")
    cat("  Number of Gradients  :", paste(x@ngrad, collapse="x"), "\n")
    cat("  Source-Filename      :", x@source, "\n")
    cat("  Slots                :\n")
    print(slotNames(x))
    invisible(NULL)
})
setMethod("print", "dwiMixtensor",
function(x){
    cat("  Object of class", class(x),"\n")
    cat("  Generated by calls    :\n")
    print(x@call)
    cat("  Dimension            :", paste(x@ddim, collapse="x"), "\n")
    cat("  Number of Gradients  :", paste(x@ngrad, collapse="x"), "\n")
    cat("  Source-Filename      :", x@source, "\n")
    cat("  Naximal number od mixture components:",max(x@order),"\n") 
    cat("  Slots                :\n")
    print(slotNames(x))
    invisible(NULL)
})
setMethod("print", "dwiQball",
function(x){
    cat("  Object of class", class(x),"\n")
    cat("  Generated by calls    :\n")
    print(x@call)
    cat("  Dimension            :", paste(x@ddim, collapse="x"), "\n")
    cat("  Number of Gradients  :", paste(x@ngrad, collapse="x"), "\n")
    cat("  Source-Filename      :", x@source, "\n")
    cat("  Kind                 :", paste(x@what, collapse="x"), "\n")
    cat("  Order                :", paste(x@order, collapse="x"), "\n")
    cat("  Slots                :\n")
    print(slotNames(x))
    invisible(NULL)
})
setMethod("print","dtiIndices",
function(x){
    cat("  Object of class", class(x),"\n")
    cat("  Generated by calls    :\n")
    print(x@call)
    cat("  Dimension            :", paste(x@ddim, collapse="x"), "\n")
    cat("  Number of Gradients  :", paste(x@ngrad, collapse="x"), "\n")
    cat("  Source-Filename      :", x@source, "\n")
    cat("  Slots                :\n")
    print(slotNames(x))
    invisible(NULL)
})

setMethod("print","dwiFiber",
function(x){
    cat("  Object of class", class(x),"\n")
    cat("  Generated by calls    :\n")
    print(x@call)
    cat("  Dimension            :", paste(x@ddim, collapse="x"), "\n")
    cat("  Number of Gradients  :", paste(x@ngrad, collapse="x"), "\n")
    cat("  Source-Filename      :", x@source, "\n")
    cat("  Minimum FA  :", x@minanindex, "\n")
    cat("  Maximum angle :", x@maxangle , "\n")
    cat("  Slots                :\n")
    print(slotNames(x))
    invisible(NULL)
})

setMethod("show", "dtiData",
function(object){
    cat("  Object of class", class(object),"\n")
    cat("  Generated by calls    :\n")
    print(object@call)
    cat("  Dimension            :", paste(object@ddim, collapse="x"), "\n")
    cat("  Number of Gradients  :", paste(object@ngrad, collapse="x"), "\n")
    cat("  Source-Filename      :", object@source, "\n")
    cat("  Slots                :\n")
    print(slotNames(object))
    invisible(NULL)
})
setMethod("show", "dtiTensor",
function(object){
    cat("  Object of class", class(object),"\n")
    cat("  Generated by calls    :\n")
    print(object@call)
    cat("  Dimension            :", paste(object@ddim, collapse="x"), "\n")
    cat("  Number of Gradients  :", paste(object@ngrad, collapse="x"), "\n")
    cat("  Source-Filename      :", object@source, "\n")
    cat("  Slots                :\n")
    print(slotNames(object))
    invisible(NULL)
})
setMethod("show", "dwiMixtensor",
function(object){
    cat("  Object of class", class(object),"\n")
    cat("  Generated by calls    :\n")
    print(object@call)
    cat("  Dimension            :", paste(object@ddim, collapse="x"), "\n")
    cat("  Number of Gradients  :", paste(object@ngrad, collapse="x"), "\n")
    cat("  Source-Filename      :", object@source, "\n")
    cat("  Naximal number od mixture components:",max(object@order),"\n") 
    cat("  Slots                :\n")
    print(slotNames(object))
    invisible(NULL)
})
setMethod("show", "dtiIndices",
function(object){
    cat("  Object of class", class(object),"\n")
    cat("  Generated by calls    :\n")
    print(object@call)
    cat("  Dimension            :", paste(object@ddim, collapse="x"), "\n")
    cat("  Number of Gradients  :", paste(object@ngrad, collapse="x"), "\n")
    cat("  Source-Filename      :", object@source, "\n")
    cat("  Slots                :\n")
    print(slotNames(object))
    invisible(NULL)
})
setMethod("show","dwiFiber",
function(object){
    cat("  Object of class", class(object),"\n")
    cat("  Generated by calls    :\n")
    print(object@call)
    cat("  Dimension            :", paste(object@ddim, collapse="x"), "\n")
    cat("  Number of Gradients  :", paste(object@ngrad, collapse="x"), "\n")
    cat("  Source-Filename      :", object@source, "\n")
    cat("  Minimum FA  :", object@minanindex, "\n")
    cat("  Maximum angle :", object@maxangle , "\n")
    cat("  Slots                :\n")
    print(slotNames(object))
    invisible(NULL)
})

setMethod("summary", "dtiData",
function(object, ...){
    cat("  Object of class", class(object),"\n")
    cat("  Generated by calls    :\n")
    print(object@call)
    cat("  Source-Filename       :", object@source, "\n")
    cat("  Dimension             :", paste(object@ddim, collapse="x"), "\n")
    cat("  Number of Gradients   :", paste(object@ngrad, collapse="x"), "\n")
    cat("  Voxel extensions      :", paste(object@voxelext, collapse="x"), "\n")
    cat("  Index of S0-Images    :", paste(object@s0ind, collapse="x"), "\n")
    cat("  Quantiles of S0-values:","\n")
    print(signif(quantile(object@si[,,,object@s0ind],...),3))
    cat("  Mean S0-value         :", paste(z <- signif(mean(object@si[,,,object@s0ind]),3),collapse="x"), "\n")
    cat("  Threshold for mask    :", paste(signif(object@level,3),collapse="x"), "\n")
    cat("\n")
    invisible(NULL)
})
setMethod("summary", "dtiTensor",
function(object, ...){
    cat("  Object of class", class(object),"\n")
    cat("  Generated by calls    :\n")
    print(object@call)
    cat("  Source-Filename       :", object@source, "\n")
    cat("  Dimension             :", paste(object@ddim, collapse="x"), "\n")
    cat("  Number of Gradients   :", paste(object@ngrad, collapse="x"), "\n")
    cat("  Voxel extensions      :", paste(object@voxelext, collapse="x"), "\n")
    cat("  Quantiles of S0-values:","\n")
    print(signif(quantile(object@th0,...),3))
    cat("  Mean S0-value         :", paste(z <- signif(mean(object@th0),3),collapse="x"), "\n")
    cat("  Voxel in mask         :", paste(sum(object@mask), collapse="x"), "\n")
    cat("  Spatial smoothness    :", paste(signif(object@bw,3), collapse="x"), "\n")
    cat("  mean variance         :", paste(signif(mean(object@sigma[object@mask]),3), collapse="x"), "\n")
    cat("  hmax                  :", paste(object@hmax, collapse="x"), "\n")
    if(length(object@outlier)>0) cat("  Number of outliers    :", paste(length(object@outlier), collapse="x"), "\n")
    cat("\n")
    invisible(NULL)
})
setMethod("summary", "dwiMixtensor",
function(object, ...){
    cat("  Object of class", class(object),"\n")
    cat("  Generated by calls    :\n")
    print(object@call)
    cat("  Source-Filename       :", object@source, "\n")
    cat("  Dimension             :", paste(object@ddim, collapse="x"), "\n")
    cat("  Number of Gradients   :", paste(object@ngrad, collapse="x"), "\n")
    cat("  Voxel extensions      :", paste(object@voxelext, collapse="x"), "\n")
    cat("  Quantiles of S0-values:","\n")
    print(signif(quantile(object@th0,...),3))
    cat("  Mean S0-value         :", paste(z <- signif(mean(object@th0),3),collapse="x"), "\n")
    cat("  Voxel in mask         :", paste(sum(object@mask), collapse="x"), "\n")
    cat("  Spatial smoothness    :", paste(signif(object@bw,3), collapse="x"), "\n")
    cat("  mean variance         :", paste(signif(mean(object@sigma[object@mask]),3), collapse="x"), "\n")
    cat("  hmax                  :", paste(object@hmax, collapse="x"), "\n")
    if(length(object@outlier)>0) cat("  Number of outliers    :", paste(length(object@outlier), collapse="x"), "\n")
    cat("  Numbers od mixture components:",table(object@order[object@mask]),"\n") 
    cat("\n")
    invisible(NULL)
})
setMethod("summary", "dwiQball",
function(object, ...){
    cat("  Object of class", class(object),"\n")
    cat("  Generated by calls    :\n")
    print(object@call)
    cat("  Source-Filename       :", object@source, "\n")
    cat("  Dimension             :", paste(object@ddim, collapse="x"), "\n")
    cat("  Number of Gradients   :", paste(object@ngrad, collapse="x"), "\n")
    cat("  Voxel extensions      :", paste(object@voxelext, collapse="x"), "\n")
    cat("  Kind                 :", paste(object@what, collapse="x"), "\n")
    cat("  Order                :", paste(object@order, collapse="x"), "\n")
    cat("  Quantiles of S0-values:","\n")
    print(signif(quantile(object@th0,...),3))
    cat("  Mean S0-value         :", paste(z <- signif(mean(object@th0),3),collapse="x"), "\n")
    cat("  Voxel in mask         :", paste(sum(object@mask), collapse="x"), "\n")
    cat("  Spatial smoothness    :", paste(signif(object@bw,3), collapse="x"), "\n")
    cat("  mean variance         :", paste(signif(mean(object@sigma[object@mask]),3), collapse="x"), "\n")
    cat("  hmax                  :", paste(object@hmax, collapse="x"), "\n")
    if(length(object@outlier)>0) cat("  Number of outliers    :", paste(length(object@outlier), collapse="x"), "\n")
    cat("\n")
    invisible(NULL)
})
setMethod("summary", "dtiIndices",
function(object, ...){
    cat("  Object of class", class(object),"\n")
    cat("  Generated by calls    :\n")
    print(object@call)
    cat("  Source-Filename       :", object@source, "\n")
    cat("  Dimension             :", paste(object@ddim, collapse="x"), "\n")
    cat("  Number of Gradients   :", paste(object@ngrad, collapse="x"), "\n")
    cat("  Voxel extensions      :", paste(object@voxelext, collapse="x"), "\n")
    cat("  Percentage of zero values      :",paste(signif(mean(object@fa==0)*100,3), "%",collapse="x"), "\n")
    cat("  Quantiles of positive FA-values:","\n")
    print(signif(quantile(object@fa[object@fa>0],...),3))
    cat("  Quantiles of positive GA-values:","\n")
    print(signif(quantile(object@ga[object@ga>0],...),3))
    cat("  Quantiles of positive MD-values:","\n")
    print(signif(quantile(object@md[object@md>0],...),3))
    cat("\n")
    invisible(NULL)
})
setMethod("summary","dwiFiber",
function(object){
    cat("  Object of class", class(object),"\n")
    cat("  Generated by calls    :\n")
    print(object@call)
    cat("  Source-Filename       :", object@source, "\n")
    cat("  Dimension             :", paste(object@ddim, collapse="x"), "\n")
    cat("  Number of Gradients   :", paste(object@ngrad, collapse="x"), "\n")
    cat("  Voxel extensions      :", paste(object@voxelext, collapse="x"), "\n")
    cat("  Minimum FA  :", object@minanindex, "\n")
    cat("  Maximum angle :", object@maxangle , "\n")
    cat("  Number of fibers :", length(object@startind), "\n")
    cat("  Quantiles of fiber lengths:\n")
    print(quantile(diff(c(object@startind,dim(object@fibers)[1]+1))/2))
    cat("  Total number of line segments :", dim(object@fibers)[1]/2,"\n")
    cat("\n")
    invisible(NULL)
})

################################################################
#                                                              #
# Section for interface functions                              #
# like: tensor2medinria()                                      #
#                                                              #
################################################################

tensor2medinria <- function(obj, filename, xind=NULL, yind=NULL, zind=NULL) {
  if (!require(fmri)) stop("cannot execute function without package fmri, because of missing write.NIFTI() function")

  if (is.null(xind)) xind <- 1:obj@ddim[1]
  if (is.null(yind)) yind <- 1:obj@ddim[2]
  if (is.null(zind)) zind <- 1:obj@ddim[3]
  if (obj@orientation[1]==1) xind <- min(xind)+max(xind)-xind
  if (obj@orientation[2]==3) yind <- min(yind)+max(yind)-yind
  if (obj@orientation[3]==4) zind <- min(zind)+max(zind)-zind
  header <- list()
  header$dimension <- c(5,length(xind),length(yind),length(zind),1,6,0,0)
  header$pixdim <- c(-1, obj@voxelext[1:3], 0, 0, 0, 0)
  header$intentcode <- 1007
  header$datatype <- 16
  header$bitpix <- 192
  header$sclslope <- 1
  header$xyztunits <- "\002" # ???
  header$qform <- 1
  header$sform <- 1
  header$quaternd <- 1
  header$srowx <- c(-2,0,0,0)
  header$srowy <- c(0,2,0,0)
  header$srowz <- c(0,0,2,0)

  write.NIFTI(aperm(obj@D,c(2:4,1))[xind,yind,zind,c(1,2,4,3,5,6)],header,filename)
  return(NULL)
}

medinria2tensor <- function(filename) {
  if (!require(fmri)) stop("cannot execute function without package fmri, because of missing read.NIFTI() function")
  args <- sys.call() 
  data <- read.NIFTI(filename)
 
  invisible(new("dtiTensor",
                call  = list(args),
                D     = aperm(extract.data(data),c(4,1:3))[c(1,2,4,3,5,6),,,],
                sigma = array(0,data$dim[1:3]),
                scorr = array(0,c(1,1,1)),
                bw    = rep(0,3),
                mask  = array(TRUE,data$dim[1:3]),
                method = "unknown",
                hmax  = 1,
                th0   = array(0,dim=data$dim[1:3]),
                gradient = matrix(0,1,1),
                btb   = matrix(0,1,1),
                ngrad = as.integer(0), # = dim(btb)[2]
                s0ind = as.integer(0),
                ddim  = data$dim[1:3],
                ddim0 = data$dim[1:3],
                xind  = 1:data$dim[1],
                yind  = 1:data$dim[2],
                zind  = 1:data$dim[3],
                voxelext = data$delta,
                orientation = as.integer(c(0,2,5)),
                rotation = diag(3),
                scale = 1,
                source= "unknown")
            )

}
