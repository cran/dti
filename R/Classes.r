setClass("dti",
         representation(.Data = "list",
                        call = "list",
                        btb    = "matrix",
                        ngrad  = "integer", # = dim(btb)[2]
                        s0ind  = "integer", # indices of s0 images
                        replind = "integer", # replications in gradient design
                        ddim   = "integer",
                        ddim0  = "integer",
                        xind   = "integer",
                        yind   = "integer",
                        zind   = "integer",
                        voxelext = "numeric",
                        level  = "numeric",
                        orientation = "integer",
                        source = "character"),
         )

setClass("dtiData",
         representation(si   = "array",
                        sdcoef = "numeric"),
         contains=c("list","dti"),
         validity=function(object){
          if (any(dim(object@si)!=c(object@ddim,object@ngrad))) {
            cat("incorrect dimension of image data\n")
            return(invisible(FALSE))
          }
          if (length(object@s0ind)<1) {
            cat("no S_0 images, parameters not identifiable \n")
             return(invisible(FALSE))
         }
          if (length(object@orientation)!=3) {
            cat("invalid orientation \n")
             return(invisible(FALSE))
         }
          if (any(sort((object@orientation)%/%2) != 0:2)) {
            cat("invalid orientation \n")
             return(invisible(FALSE))
         }
          if (object@level < 0) {
            cat("invalid level \n")
             return(invisible(FALSE))
         }
          if (length(object@sdcoef) != 4) {
             cat("invalid model for error standard deviation \n")
             return(invisible(FALSE))
         }
          if (any(object@sdcoef<0)||object@sdcoef[3]>object@sdcoef[4]) {
             cat("illegal interval of linearity in model for  error standard deviation \n")
             return(invisible(FALSE))
         }
         }
         )
setClass("dtiTensor",
         representation(method = "character",
                        D      = "array",
                        th0    = "array",
                        sigma  = "array",
                        scorr  = "array",
                        bw     = "numeric",
                        mask   = "array",
                        hmax   = "numeric",
                        outlier = "numeric",
                        scale  = "numeric"),
         contains=c("list","dti"),
         validity=function(object){
          if (any(dim(object@D)!=c(6,object@ddim))) {
            cat("invalid dimension of tensor array D \n")
            return(invisible(FALSE))
          }
          if (any(dim(object@th0)!=object@ddim)) {
            cat("invalid dimension of array th0\n")
            return(invisible(FALSE))
          }
          if (object@method=="linear"&any(dim(object@sigma)!=object@ddim)) {
            cat("invalid dimension of array sigma\n")
            return(invisible(FALSE))
          }
          if (any(dim(object@mask)!=object@ddim)) {
            cat("dimension of mask:",dim(object@mask),"\n")
            cat("should be:",object@ddim,"\n")
            cat("invalid dimension of array mask\n")
            return(invisible(FALSE))
          }
          if (!is.logical(object@mask)) {
            cat("invalid type of array mask, should be logical\n")
            return(invisible(FALSE))
          }
          if (length(dim(object@scorr))!=3) {
            cat("invalid dimension of scorr\n")
            return(invisible(FALSE))
          }
          if (length(object@bw)!=3) {
            cat("invalid length of bw\n")
            return(invisible(FALSE))
          }
          if (!(object@method %in% c("linear","nonlinear","unknown"))) {
            cat("method should specify either linear or nonlinear or unknown\n")
            return(invisible(FALSE))
          }
         }
         )

setClass("dtiIndices",
         representation(method = "character",
                        fa     = "array",
                        ga     = "array",
                        md     = "array",
                        andir  = "array",
                        bary   = "array"),
         contains=c("list","dti"),
          validity=function(object){
          if (any(dim(object@fa)!=object@ddim)) {
            cat("invalid dimension of array fa\n")
            return(invisible(FALSE))
          }
          if (any(dim(object@ga)!=object@ddim)) {
            cat("invalid dimension of array ga\n")
            return(invisible(FALSE))
          }
          if (any(dim(object@md)!=object@ddim)) {
            cat("invalid dimension of array ra\n")
            return(invisible(FALSE))
          }
          if (any(dim(object@andir)!=c(3,object@ddim))) {
            cat("invalid dimension of array andir\n")
            return(invisible(FALSE))
          }
          if (any(dim(object@bary)!=c(3,object@ddim))) {
            cat("invalid dimension of array bary\n")
            return(invisible(FALSE))
          }
         }
        )
