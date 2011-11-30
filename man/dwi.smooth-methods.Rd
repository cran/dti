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
  \S4method{dwi.smooth}{dtiData}(object,kstar,lambda=10,kappa0=0.4,sigma=NULL,sigma2=NULL,minsb=5,smooths0=TRUE,xind=NULL,yind=NULL,zind=NULL,verbose=FALSE)
}
\arguments{
  \item{object}{ Object of class \code{"dtiData"}}
  \item{kstar}{ Number of steps in structural adaptation}
  \item{lambda}{Scale parameter in adaptation }
  \item{kappa0}{determines amaount of smoothing on the sphere. Larger values correspond to stronger smoothing on the sphere. }
  \item{sigma}{Error standard deviation. Assumed to be known and homogeneous in the current implementation. A reasonable estimate may be defined 
  as the modal value of standard deviations obtained using method \code{getsdofsb}. }
  \item{sigma2}{estimated variances in each voxel (the corresponding code is not fully tested)}
  \item{minsb}{lower threshold used to define a mask. Voxel with mean intensity (over all diffusion weighted images) less then \code{minsb} are not used. }
  \item{smooths0}{ also smooth unweighted images. }
  \item{xind}{index for x-coordinate}
  \item{yind}{index for y-coordinate}
  \item{zind}{index for z-coordinate}
  \item{verbose}{If \code{verbose=TRUE} additional reports are given.}
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
