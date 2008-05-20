\name{dtiTensor-class}
\docType{class}
\alias{dtiTensor-class}
\alias{coerce,dtiTensor,dtiIndices-method}
\alias{dtiIndices,dtiTensor-method}
\alias{dti.smooth,dtiTensor-method}
\alias{plot,dtiTensor-method}

\title{Class "dtiTensor"}
\description{ Diffusion Tensor object. Containes the results of estimating a diffusion tensor from diffusion weighted image (DWI) data.}
\section{Objects from the Class}{
Objects can be created by calls to functions \code{\link{dtiTensor}},  \code{\link{dti.smooth}} or \code{\link{medinria2tensor}}.
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} with optional components
    \item{\code{s2rician}}{estimated variance parameter in the rice distribution
    (if \code{rician=TRUE} in function \code{dti.smooth}} 
    \item{\code{ni}}{voxelwise sum of weights in adaptive smoothing (from function \code{dti.smooth})} 
}
    \item{\code{D}:}{Object of class \code{"array"} contains estimated tensors, dimension \code{c(6,ddim)}.
          Tensors are stored as upper diagonal matrices.} 
    \item{\code{th0}:}{Object of class \code{"array"} contains estimated intensities in S0 images, dimension \code{ddim}} 
    \item{\code{sigma}:}{Object of class \code{"array"} containing estimated error variances} 
    \item{\code{scorr}:}{Object of class \code{"numeric"} containing  estimated spatial correlations in coordinate directions} 
    \item{\code{bw}:}{Object of class \code{"numeric"} containing bandwidth for a Gaussian kernel that approximately creates
         the estimated spatial correlations. Needed for adjustments of critical values in the adaptive smoothing algorithm used in function \code{dti.smooth}} 
    \item{\code{mask}:}{Object of class \code{"array"}, logical array indicating the voxel where the tensor was estimated.} 
    \item{\code{level}:}{Object of class \code{"numeric"} minimal valid s0-level. 
     No evaluation was be performed for voxel with s0-values less than \code{level}. Used to determine  \code{mask}.}
    \item{\code{hmax}:}{Object of class \code{"numeric"} maximal bandwidth in case of adaptive smoothing. contains 1 otherwise.} 
    \item{\code{method}:}{Object of class \code{"character"} either \code{"linear"} or \code{"nonlinear"} or \code{"unknown"}. Indicates the regression model used for estimating the tensors.}     
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
    \item{dtiIndices}{\code{signature(object = "dtiTensor")}: Calculate fractal anisotropy (FA) and main directions of anisotropy from diffusion tensors. }
    \item{dti.smooth}{\code{signature(object = "dtiTensor")}: Smooth diffusion tensor. For exploration only. We strictly recommend using function \code{dti.smooth} on a \code{dtiData}-object. }
    \item{plot}{\code{signature(x = "dtiTensor")}: not jet implemented. }
	 }
}
\references{ K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, NeuroImage 39(4), 1763-1773 (2008).

\url{http://www.wias-berlin.de/projects/matheon_a3/} 
 }
\author{ Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de  }
\seealso{\code{\link{dtiData}},\code{\link{dtiTensor}},\code{\link{dtiIndices}},\code{\link{dti.smooth}}, \code{\linkS4class{dti}}, \code{\linkS4class{dtiData}}, \code{\linkS4class{dtiIndices}}
}
\examples{
showClass("dtiTensor")
\dontrun{demo(dti_art)}
}
\keyword{classes}
