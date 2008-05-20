\name{plot}
\alias{plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Visualizing of dtiIndices}
\description{
Generic plot function for visualization of objects of class
"dtiData", "dtiIndices" and class "dtiTensor"
}
\usage{
plot(x, y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ object of class \code{dtiIndices}, \code{dtiData} or \code{dtiTensor}}
  \item{y}{ not used }
  \item{\dots}{ additional arguments such as
   \item{slice }{Slice number}
   \item{view}{ choose sagittal, coronal, or axial view here}
   \item{gradient}{index of data cube to chose. (for class
\code{dtiData}) only. Defaults to the first S0 image.}
   \item{method }{\code{method==1} and \code{method==2} specify two different color schemes for 
    directional FA-maps. \code{method==3} specifies visualization of dtiIndices using color coded bary coordinates. (for class
\code{dtiIndices} only)}
   \item{quant }{if \code{is.null(minanindex)} specify \code{minanindex} as corresponding quantile of the fractal anisotropy (FA) index.}
   \item{minanindex }{display only information for voxel with \code{FA>minanindex}}
   \item{show }{Visualize information in a graphics device (for classes
\code{dtiData} and \code{dtiIndices} only)}
   \item{contrast.enh}{ enhance image contrast using \code{min(1,x$anindex/contrast.enh} instead
    of the anisotropy index itself. Effective values are within the interval (0,1).}
   \code{settings for function\code{par} (for class
\code{dtiTensor} only)}
   \code{settings for function\code{show.image} (for class
\code{dtiData} and \code{dtiIndices} only)}
   \code{xind}{If provided restrict display to indices specified in \code{xind} for x-direction.}
   \code{yind}{If provided restrict display to indices specified in \code{xind} for y-direction.}
   \code{zind}{If provided restrict display to indices specified in \code{xind} for z-direction.}
}
}
\details{
This function allows to select a slide (2D) by arguments 
\code{slice} and \code{view}. For objects of class \code{"dtiData"} \code{gradient} can be used to specify a specific data cube
associated with the index of a gradient direction.
For objects of class \code{"dtiData"} images are produces that are
scaled by the maximal observed image value. This guarantees 
that subsequently produced images are on a comparable grey scale.
The resulting image is returned.

For objects of class \code{"dtiTensor"} the tensor itself, fractional anisotropy (FA), mean diffusivity (MD) and a color coded anisotropy map are provided. \code{NULL} is returned.
 
For objects of class \code{"dtiIndices"} color coded anisotropy
maps are produced depending on the specification in \code{method}.
The resulting image is returned.


 
}
\value{
an image of class \code{adimpro} for objects of class \code{"dtiData"} or \code{"dtiIndices"}, \code{invisible(NULL)} otherwise.
}
\references{ K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, NeuroImage 39(4), 1763-1773 (2008).

\url{http://www.wias-berlin.de/projects/matheon_a3/} }
\author{ Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de }
\seealso{ \code{\link{dtiIndices}}, \code{\link{dtiData}}, \code{\link{dtiTensor}}}
\examples{
\dontrun{demo(dti_art)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
