\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,ANY-method}
\alias{plot,dwi-method}
\alias{plot,dtiData-method}
\alias{plot,dtiIndices-method}
\alias{plot,dtiTensor-method}
\title{ Methods for Function `plot' in Package `dti' } 
\description{
  Visualization of objects of class \code{"dtiData"}, \code{"dtiIndices"} and class \code{"dtiTensor"}
}
\usage{
  \S4method{plot}{dtiData}(x, y, slice=1, gradient=NULL, view= "axial", show=TRUE, density=FALSE, xind=NULL, yind=NULL, zind=NULL, mar=c(3,3,3,.3), mgp=c(2,1,0), \dots)
  \S4method{plot}{dtiTensor}(x, y, slice=1, view="axial", quant=0, minanindex=NULL, contrast.enh=1, what="FA", qrange=c(.01,.99), xind=NULL, yind=NULL, zind=NULL, mar=c(2,2,2,.2), mgp=c(2,1,0), \dots)
  \S4method{plot}{dtiIndices}(x, y, slice=1, view= "axial", method=1, quant=0, minanindex=NULL, show=TRUE, density=FALSE, contrast.enh=1, what="FA", xind=NULL, yind=NULL, zind=NULL, mar=c(3,3,3,.3), mgp=c(2,1,0), \dots)
}
\arguments{
  \item{x}{ Object of class \code{"dtiIndices"}, \code{"dtiData"} or \code{"dtiTensor"}}
  \item{y}{ Not used }
  \item{slice}{ Slice number }
  \item{view}{ Choose \code{"sagittal"}, \code{"coronal"}, or \code{"axial"} view here }
  \item{gradient}{ Index of data cube to plot. Defaults to the first S0 image. }
  \item{method}{ Method for color coding tensor indices. }
  \item{quant}{ If \code{is.null(minanindex)} specify \code{minanindex} as corresponding quantile of the fractal anisotropy (FA) index.}
  \item{minanindex}{ Display only information for voxel with \code{FA>minanindex}}
  \item{show}{Visualize information in a graphics device (for classes \code{"dtiData"} and \code{"dtiIndices"} only).  }
  \item{density}{Show density of S0(Sb)-values (for class \code{"dtiData"}) or densities of fractal anisotropy (FA) or geodesic anisotropy (GA) ( for class \code{"dtiIndices"}).  }
  \item{contrast.enh}{ Enhance image contrast using \code{min(1,x$anindex/contrast.enh} instead of the anisotropy index itself. Effective values are within the interval (0,1).}
  \item{what}{ If \code{what="GA"} use geodesic anisotropy (GA) instead of fractional anisotropy (FA).  }
  \item{mar}{ Graphical parameter for \code{par}. }
  \item{mgp}{ Graphical parameter for \code{par}. }
  \item{qrange}{ Cut image intensity to these quantiles to avoid that outliers determine the dynamic range of the image. }
  \item{xind}{If provided restrict display to indices specified in \code{xind} for x-direction.}
  \item{yind}{If provided restrict display to indices specified in \code{yind} for y-direction.}
  \item{zind}{If provided restrict display to indices specified in \code{zind} for z-direction.} 
  \item{\dots}{ currently not used }
}
\section{Methods}{
\describe{
  \item{x = "ANY"}{ Generic function: see \code{\link[graphics]{plot}}. }

  \item{x =  "dwi"}{ Returns a warning. }

  \item{x = "dtiData"}{ \code{gradient} can be used to specify a specific data cube associated with the index of a gradient direction. For objects of class \code{"dtiData"} images are produces that are scaled by the maximal observed image value. This guarantees that subsequently produced images are on a comparable grey scale. The resulting image of class \code{"adimpro"} from package \pkg{adimpro} is returned. }

  \item{x = "dtiIndices"}{ Color coded anisotropy maps are produced depending on the specification in \code{method}. \code{method==1}, \code{method==2}, \code{method==4}  and \code{method==5} specify three different color schemes for directional FA-maps. \code{method==3} specifies visualization of dtiIndices using color coded shape parameters. The resulting image of class \code{"adimpro"} from package \pkg{adimpro} is returned.}

  \item{x = "dtiTensor"}{ The tensor itself, fractional anisotropy (FA), mean diffusivity (MD) and a color coded anisotropy map are provided. \code{NULL} is returned. }
  
}}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\linkS4class{dtiIndices}}, 
  \code{\linkS4class{dtiData}}, 
  \code{\linkS4class{dtiTensor}}
}
\examples{
  \dontrun{demo(dti_art)}
}
\keyword{ methods }
\keyword{ hplot }
