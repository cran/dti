\name{dtiData-class}
\docType{class}
\alias{dtiData-class}
\alias{dti.smooth,dtiData-method}
\alias{dtiTensor,dtiData-method}
\alias{plot,dtiData-method}

\title{Class "dtiData"  }
\description{ Diffusion Weighted Image (DWI) Data }
\section{Objects from the Class}{
Objects can be created by calls of function \code{\link{dtiData}}.
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} with component \item{\code{si}}{4D integer array of image values. Dimension \code{c(ddim,ngrad)}} }
    \item{\code{level}:}{Object of class \code{"numeric"} minimal valid s0-level. 
     No evaluation will be performed for voxel with s0-values less than \code{level}.  }
    \item{\code{btb}:}{Object of class \code{"matrix"} matrix of dimension \code{c(6,ngrad)} obtained from gradient directions.}
    \item{\code{ngrad}:}{Object of class \code{"integer"} number of gradients (including zero gradients) }
    \item{\code{s0ind}:}{Object of class \code{"integer"} index of zero gradients within sequence \code{1:ngrad} }
    \item{\code{replind}:}{Object of class \code{"integer"} index (identifier) of unique 
           gradient directions. Used to charactreize replications in the gradient design by identical indices. length is \code{ngrad}}
    \item{\code{ddim}:}{Object of class \code{"integer"} dimension of original image cubes. Integer vector of length 3. }
    \item{\code{ddim0}:}{Object of class \code{"integer"} dimension of subcube defined by \code{xind}, \code{yind} and \code{zind}. }
    \item{\code{xind}:}{Object of class \code{"integer"} index for subcube definition in x-direction  }
    \item{\code{yind}:}{Object of class \code{"integer"} index for subcube definition in y-direction }
    \item{\code{zind}:}{Object of class \code{"integer"} index for subcube definition in z-direction }
    \item{\code{voxelext}:}{Object of class \code{"numeric"} voxel extensions in x-, y- and z-direction. vector of length 3. }
    \item{\code{orientation}:}{Object of class \code{"integer"} vector of length 3. Orientation of data according to AFNI convention.}
    \item{\code{source}:}{Object of class \code{"character"} name of the imgfile used to create the data. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{dti}"}, directly.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{dti.smooth}{\code{signature(object = "dtiData")}: Create estimates of diffusion tensors in each voxel using structural adaptive spatial smoothing. }
    \item{dtiTensor}{\code{signature(object = "dtiData")}: Create estimates of diffusion tensors in each voxel. }
    \item{plot}{\code{signature(x = "dtiData")}: not yet implemented }
	 }
}
\references{ K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, NeuroImage 39(4), 1763-1773 (2008).

\url{http://www.wias-berlin.de/projects/matheon_a3/} 
 }
\author{ Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de  }
\seealso{\code{\link{dtiData}},\code{\link{dtiTensor}}, \code{\link{dti.smooth}}, \code{\linkS4class{dti}}, \code{\linkS4class{dtiTensor}}, \code{\linkS4class{dtiIndices}}
}
\examples{
showClass("dtiData")
\dontrun{demo(dti_art)}
}
\keyword{classes}
