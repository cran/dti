\name{show3d-methods}
\docType{methods}
\alias{show3d}
\alias{show3d-methods}
\alias{show3d,ANY-method}
\alias{show3d,dtiData-method}
\alias{show3d,dtiIndices-method}
\alias{show3d,dtiTensor-method}
\alias{show3d,dkiTensor-method}
\alias{show3d,dwiMixtensor-method}
\alias{show3d,dwiQball-method}
\alias{show3d,dwiFiber-method}
\title{ Methods for Function `show3d' in Package `dti' }
\description{
  The function provides 3D visualization of \code{"dtiData"}, \code{"dtiTensor"}, \code{"dwiQball"} and \code{"dtiIndices"} objects using the \code{"rgl"}-package. Functionality of the rgl-package allows to rotate and zoom the displayed object.
}
\usage{
\S4method{show3d}{dtiData}(obj, xind=NULL, yind=NULL, zind=NULL, quant=.8, 
         scale=.4,
         bgcolor="black", add=FALSE, maxobjects=729, what=c("adc","data"), 
         minalpha=1, nn=1, normalize=FALSE, box=FALSE, title=FALSE, \dots)
\S4method{show3d}{dtiTensor}(obj, xind=NULL, yind=NULL, zind=NULL, method=1,
         minfa=.3, mask=NULL, fibers=FALSE, 
         maxangle = 30,level=0, quant=.8, scale=.4, bgcolor="black", add=FALSE, 
         subdivide=2, maxobjects=729, what=c("tensor","adc","odf"), odfscale = 1, 
         minalpha=.25, normalize=NULL, box=FALSE, title=FALSE,\dots)
\S4method{show3d}{dkiTensor}(obj, xind=NULL, yind=NULL, zind=NULL, method=1, 
         minfa=.3, mask=NULL, level=0, quant=.8, scale=.4, bgcolor="black", 
         add=FALSE, subdivide=2, maxobjects=729, what=c("KT", "DT"),  
         minalpha=.25, normalize=NULL, box=FALSE, title=FALSE,\dots)
\S4method{show3d}{dtiIndices}(obj, index=c("fa","ga"), xind=NULL, yind=NULL, 
         zind=NULL, method=1, 
         minfa=0, bgcolor="black", add=FALSE, lwd=1, box=FALSE, 
         title=FALSE, \dots)
\S4method{show3d}{dwiMixtensor}(obj, xind=NULL, yind=NULL, zind=NULL, minfa=.3, 
         minorder = 1, mineo=1, fibers=FALSE, maxangle=30, level=0,           
         quant=.8, scale=.4, bgcolor="black", add=FALSE, 
         subdivide=3, maxobjects=729, what=c("odf","axis","both"), odfscale=1, 
         minalpha=1, lwd=3, box=FALSE, title=FALSE, \dots) 
\S4method{show3d}{dwiQball}(obj, xind=NULL, yind=NULL, zind=NULL, level=0, 
         quant=.8, scale=0.4, odfscale=1, bgcolor="black", add=FALSE, 
         subdivide=3, maxobjects=729, minalpha=1, box=FALSE, 
         title=FALSE, \dots)
\S4method{show3d}{dwiFiber}(obj, add=FALSE, bgcolor="black", box=FALSE, 
         title=FALSE, lwd=1, \dots)
}
\arguments{
\item{obj}{ An object of class \code{dtiData}, \code{dtiTensor}, \code{dtiIndices},\code{dwiMixTensor} or \code{dwiQball}}
  \item{xind}{ vector of x-coordinates, defaults to whole range. }
  \item{yind}{ vector of y-coordinates, defaults to whole range. }
  \item{zind}{ vector of z-coordinates, defaults to whole range. }
\item{quant}{Quantile of maximal radii of objects used for scaling.}
\item{scale}{Scale factor for the size of objects}
\item{bgcolor}{Backgroundcolor for rgl-display}
\item{add}{If true information is added to the current device, otherwise a new device is opened.}
\item{maxobjects}{Maximal size of data cube (in voxel) to display}
\item{minalpha}{Minimum value for transparency.}
\item{nn}{Number of nearest neighbors used for interpolation onto a regular polyeder.}
\item{normalize}{If \code{TRUE} normalize values (project to interval (0,1) within each voxel). For tensor objects \code{normalize=NULL} specifies a default depending on the content of argument what 

(\code{normalize <- switch(what,"tensor"=FALSE,"adc"=TRUE)}).}
\item{box}{Logical, add a bounding box.}
\item{title}{Either a character string specifying a title or a logical. If \code{title==TRUE} a default title characterizing the type of plot is generated. }
\item{method}{\code{method==1} and \code{method==2} specify two different color schemes for directional FA-maps.}
\item{minfa}{Minimal FA value for dtiTensor objects and for dwiMixtensor objects.}
\item{mask}{additional mask for dtiTensor objects.}
\item{minorder}{Minimal order for dwiMixtensor objects.}
\item{mineo}{Minimal effective order for dwiMixtensor objects.}
\item{fibers}{If \code{TRUE} show fibers starting in voxel with \code{fa>=minfa}, \code{order>=minorder} and \code{eorder>=mineo}, the last two effective for dwiMixtensor objects only.}
\item{maxangle}{argument for fibertracking}
\item{level}{Radius of sphere used as support for ODF visualisation}
\item{subdivide}{Level of subdivisions for meshing, level \code{0:4} correspond to use of
\code{c(12,42,162,642,2562)} vertices per tensor, respectively.}
\item{what}{For dtiTensor-objects either \code{"tensor"} for visualization using ellipsoids, \code{"adc"} for Apparent Diffusion Coefficients or \code{"odf"} for the Orientation Density Function. 
For dwiMixtensor-objects possible specifications are \code{"odf"}, \code{"axis"}
and \code{"both"}, with the latter superposing the estimated main directions on the estimated ODF. For \code{"axis"}(and \code{"both"}) the length of the axis corresponds to the mixture weights.
For dtiData-objects choices are either \code{"data"}  or \code{"adc"}.}
\item{odfscale}{Determines visualisation of the Orientation density function (ODF). For \code{odfscale=3} the
  ODF values are rescaled such that the volume of the displayes objects is constant. \code{odfscale=1} uses
  the values of the ODF as radii in the corrsponding vertice direction of the specified polyhedron. This
  can lead to extremely large volumes in case of one mixture component with high excentricity. 
  values of \code{odfscale} inbetween 1 and 3 are possible and allow to balance between volume based
  visualization and emphazising highly structured ODF's.  }
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
 
  \item{obj = "dkiTensor"}{Prelimenary show3d method for diffusion kurtosis tensors.}
 
  \item{obj = "dwiQball"}{Estimated ODF/ADC's are visualized at the voxel centers. Color 
  is determined by directions, ODF/ADC values are reflected by both radial extend and
  transparancy. The value of \code{maxobjects} limits the size of datacube and may be increased
 on hardware with suitable graphics capabilities. }

  \item{obj = "dwiFiber"}{Display and combine fibres generated by function tracking. }

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
