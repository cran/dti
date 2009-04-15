\name{show3d-methods}
\docType{methods}
\alias{show3d}
\alias{show3d-methods}
\alias{show3d,ANY-method}
\alias{show3d,dtiData-method}
\alias{show3d,dtiIndices-method}
\alias{show3d,dtiTensor-method}
\alias{show3d,dwiQball-method}
\title{ Methods for Function `show3d' in Package `dti' }
\description{
  The function provides 3D visualization of \code{"dtiData"}, \code{"dtiTensor"}, \code{"dwiQball"} and \code{"dtiIndices"} objects using the \code{"rgl"}-package. Functionality of the rgl-package allows to rotate and zoom the displayed object.
}
\usage{
  \S4method{show3d}{dtiData}(obj,nx=NULL,ny=NULL,nz=NULL,center=NULL,scale=.5,bgcolor="black",add=FALSE,maxobjects=729,what="ADC",minalpha=1,power=1,nn=1,normalize=FALSE,box=FALSE,title=FALSE,\dots)
  \S4method{show3d}{dtiTensor}(obj,nx=NULL,ny=NULL,nz=NULL,center=NULL,method=1,level=0,scale=.5,bgcolor="black",add=FALSE,subdivide=2,maxobjects=729,what="ADC",minalpha=.25,normalize=NULL,box=FALSE,title=FALSE,\dots)
  \S4method{show3d}{dtiIndices}(obj, index="FA", nx=NULL, ny=NULL, nz=NULL, center=NULL, method=1, level=0, bgcolor="black", add=FALSE, lwd=1,box=FALSE,title=FALSE, \dots)
  \S4method{show3d}{dwiQball}(obj,nx=NULL,ny=NULL,nz=NULL,center=NULL,scale=0.5,bgcolor="black",add=FALSE,subdivide=3,maxobjects=729,minalpha=1,power=1,normalize=TRUE,box=FALSE,title=FALSE, \dots)
}
\arguments{
\item{obj}{ An object of class \code{dtiTensor} or \code{dtiIndices} }
\item{nx}{Number of voxel in x-direction}
\item{ny}{Number of voxel in y-direction}
\item{nz}{Number of slices}
\item{center}{Vector of length 3 specifying the center of the data cube (class \code{dtiData} or \code{dtiTensor}) or center of display (class \code{dtiIndices})}
\item{scale}{Scale factor for the size of polyeders \code{dtiData} and ellipsoids (\code{dtiTensor})}
\item{bgcolor}{Backgroundcolor for rgl-display}
\item{add}{If true information is added to the current device, otherwise a new device is opened.}
\item{maxobjects}{Maximal size of data cube (in voxel) to display}
\item{minalpha}{Minimum value for transparency.}
\item{power}{Exponent in visualization, defaults to 1. heighter values may be used for increased contrast.}
\item{nn}{Number of nearest neighbors used for interpolation onto a regular polyeder.}
\item{normalize}{If \code{TRUE} normalize values (project to interval (0,1) within each voxel). For tensor objects \code{normalize=NULL} specifies a default depending on the content of argument what (\code{normalize <- switch(tolower(what),"tensor"=FALSE,"adc"=TRUE)}).}
\item{box}{Logical, add a bounding box.}
\item{title}{Either a character string specifying a title or a logical. If \code{title==TRUE} a default title characterizing the type of plot is generated. }
\item{method}{\code{method==1} and \code{method==2} specify two different color schemes for directional FA-maps.}
\item{level}{Minimal FA value of tensors.}
\item{subdivide}{Level of subdivisions for meshing, level \code{0:4} correspond to use of
\code{c(12,42,162,642,2562)} vertices per tensor, respectively.}
\item{what}{For dtiTensor-objects either \code{"tensor"} for visualization using ellipsoids or \code{"ADC"} for Apparent Diffusion Coefficients. For dtiData-objects choices are either \code{"data"}  or \code{"ADC"}.}
\item{lwd}{Linewidth for visualization of dtiIndices objects.}
\item{index}{Eiter \code{"FA"} for fractional anisotropy index or 
\code{"GA"} for geodesic anisotropy index.}
\item{\dots}{Additional parameters passed to function \code{rgl.par} from the rgl-package.}
}
\section{Methods}{
\describe{

  \item{obj = "ANY"}{ Returns a warning }

  \item{obj = "dtiData"}{Empirical ADC's are visualized at the voxel centers. Color 
  is determined by gradient directions, ADC values are reflected by both radial extend and
  transparancy. The value of \code{maxobjects} limits the size of datacube and may be increased
 on hardware with suitable graphics capabilities. }

  \item{obj = "dtiIndices"}{ Objects are visualized as a collection of line segments with location given by the voxel center, orientation and color determined by the main direction of inisotropy and length corresponding to either fractional or geodesic anisotropy as specified in \code{index}. 

  Displayed objects are restricted to voxel with an fractional (geodesic) anisotropy larger than \code{level}.}

  \item{obj = "dtiTensor"}{Ellipsoids/ADC's are visualized at the voxel centers. Orientation and size correspond to the tensor
  values, color is determined by the main direction of anisotropy using the colorsceme specified with \code{method}. The fractional anisotropy value is coded as transparency. The value of \code{maxobjects} limits the size of datacube and may be increased
 on hardware with suitable graphics capabilities. }
 
  \item{obj = "dwiQball"}{Estimated ODF/ADC's are visualized at the voxel centers. Color 
  is determined by directions, ODF/ADC values are reflected by both radial extend and
  transparancy. The value of \code{maxobjects} limits the size of datacube and may be increased
 on hardware with suitable graphics capabilities. }

  Displays can be closed using function \code{rgl.close}
}
}
\value{
  The function returns the number of the current rgl-device.
}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\link{dtiIndices-methods}},
  \code{\link{dti.smooth}}, 
  \code{\linkS4class{dtiTensor}},
  \code{\linkS4class{dtiIndices}}
}
\examples{
  \dontrun{demo(dti_art)}
}
\keyword{ methods }
\keyword{ iplot }
\keyword{ hplot }
