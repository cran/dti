\name{dwi.smooth-methods}
\docType{methods}
\alias{dwi.smooth}
\alias{dwi.smooth-methods}
\alias{dwi.smooth,ANY-method}
\alias{dwi.smooth,dtiData-method}
\title{Smooth DWI data}
\description{
 Adaptive smoothing of DWI data. Smoothing is performed both in space and on the sphere (e.g. between images  obtained for different gradient directions) employing a natural
geometrical distance ( in SE(3)). Structural adaptation is used in space only.}
\usage{
  \S4method{dwi.smooth}{dtiData}(object,kstar,lambda=6,kappa0=NULL,ncoils=1,sigma=NULL,sigma2=NULL,level=NULL,minsb=5,vred=4,xind=NULL,yind=NULL,zind=NULL,verbose=FALSE,dist=1,model="Chi2")
}
\arguments{
  \item{object}{ Object of class \code{"dtiData"}}
  \item{kstar}{ Number of steps in structural adaptation}
  \item{lambda}{Scale parameter in adaptation }
  \item{kappa0}{determines amount of smoothing on the sphere. Larger values correspond to stronger smoothing on the sphere. If \code{kappa0=NULL} a value is that corresponds to a variace reduction with factor \code{vred} on the sphere.}
  \item{ncoils}{Number of coils in MR system }
  \item{sigma}{Error standard deviation. Assumed to be known and homogeneous in the current implementation. A reasonable estimate may be defined 
  as the modal value of standard deviations obtained using method \code{getsdofsb}. }
  \item{sigma2}{estimated variances in each voxel (the corresponding code is not fully tested)}
  \item{level}{Threshold for image intensities when setting mask.}
  \item{minsb}{lower threshold used to define a mask. Voxel with mean intensity (over all diffusion weighted images) less then \code{minsb} are not used. }
  \item{vred}{Used if \code{kappa0=NULL} to specify the variance reduction on the sphere when suggesting a value of kappa0}
  \item{xind}{index for x-coordinate}
  \item{yind}{index for y-coordinate}
  \item{zind}{index for z-coordinate}
  \item{verbose}{If \code{verbose=TRUE} additional reports are given.}
  \item{dist}{Distance in SE3. Reasonable values are 1 (default, see Becker et.al. 2012), 2 ( a slight modification of 1: with k6^2 instead of abs(k6)) and 3 (using a 'naive' distance on the sphere)}
  \item{model}{Determines which quantities are smoothed. Possible values are
  \code{"Chi"} for observed values (assumed to be distributed as noncentral
  Chi with \code{2*ncoils} degrees of freedom) and \code{"Chi2"} for squares of
  observed values (assumed to be distributed as noncentral
  Chi-squared with \code{2*ncoils} degrees of freedom)}
  }
\value{ 
  An object of class \code{"dtiData"} with smoothed diffusion weighted images.
}
\section{Methods}{
\describe{

\item{\code{signature(object) = "ANY"}}{
Returns a warning.
}

\item{\code{signature(object) = "dtiData"}}{
Smoothing of DWI data 
}
}}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\link{dtiData}}, 
  \code{\linkS4class{dtiData}}, 
} 
\keyword{methods}
\keyword{smooth }
