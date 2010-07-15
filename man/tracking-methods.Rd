\name{tracking-methods}
\docType{methods}
\alias{tracking}
\alias{selectFibers}
\alias{tracking-methods}
\alias{selectFibers-methods}
\alias{tracking,ANY-method}
\alias{tracking,dtiTensor-method}
\alias{tracking,dtiIndices-method}
\alias{tracking,dwiMixtensor-method}
\alias{selectFibers,dwiFiber-method}
\title{ Methods for Function `tracking' in Package `dti' }
\description{
  The function provides fiber tracking of \code{"dtiTensor"}, \code{"dtiIndices"}, and \code{"dwiMixtensor"} objects.
}
\usage{
\S4method{tracking}{dtiTensor}(obj, roix=NULL, roiy=NULL, roiz=NULL, method="LINEPROP",
                    minanindex=0.3, maxangle=30, subsample = 1)
\S4method{tracking}{dtiIndices}(obj, roix=NULL, roiy=NULL, roiz=NULL, method="LINEPROP",
                    minanindex=0.3, maxangle=30, subsample = 1)
\S4method{tracking}{dwiMixtensor}(obj, roix=NULL, roiy=NULL, roiz=NULL, method="LINEPROP",
                    minanindex=0.3, maxangle=30, subsample = 1)
\S4method{selectFibers}{dwiFiber}(obj, roix=NULL, roiy=NULL, roiz=NULL, nroimask=NULL, 
                    minlength=1)
}
\arguments{
\item{obj}{ An object of class \code{"dtiTensor"}, \code{"dtiIndices"}, or \code{"dwiMixtensor"} for tracking() and \code{"dwiFiber"} for selectFiber().}
\item{roix}{Indices defining the ROI in x direction. Currently min/max is used to define ROIx}
\item{roiy}{Indices defining the ROI in y direction. Currently min/max is used to define ROIy}
\item{roiz}{Indices defining the ROI in z direction. Currently min/max is used to define ROIz}
\item{method}{Method for fibre tracking. "LINEPROP" is simple line propagation algorithm which is the default.}
\item{minanindex}{Minimal FA to follow the tracks. default 0.3}
\item{maxangle}{Maximal angle between fiber in adjacent voxels. default 30 degree.}
\item{subsample}{Subsampling order of the data to get more dense fibre tracks. Note, that objects become very(!) large.}
\item{nroimask}{Array of image size with 1 at voxel within the ROI to be choosen, and 0 otherwise.}
\item{minlength}{Minimal length of fibers to be selected.}
}
\section{Methods}{
\describe{

  \item{obj = "dtiTensor"}{ Fiber tracking is performed on the estimated vector field of principal diffusion direction using the method \code{method}. Currently only line propagation is implemented. The resulting tracks can be visualized using function \code{\link{show3d}}}. 
 
  \item{obj = "dtiIndices"}{ Fiber tracking is performed on the estimated vector field of principal diffusion direction using the method \code{method}. Currently only line propagation is implemented. The resulting tracks can be visualized using function \code{\link{show3d}}}. 
  
  \item{obj = "dwiMixtensor"}{ Fiber tracking is performed on the estimated vector fields of diffusion direction in the mixed tensor model using the method \code{method}. Currently only line propagation is implemented. The resulting tracks can be visualized using function \code{\link{show3d}}}. 
  Displays can be closed using function \code{rgl.close}
}
}
\value{
  The function returns an object of class \code{\linkS4class{dwiFiber}}.
}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}, Joerg Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\linkS4class{dtiTensor}},\code{\linkS4class{dtiIndices}},\code{\linkS4class{dwiFiber}},
  \code{\link{show3d}}, \code{\link{summary}}, \code{\link{print}}
}
\keyword{ methods }
