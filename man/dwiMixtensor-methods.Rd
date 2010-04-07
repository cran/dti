\name{dwiMixtensor-methods}
\docType{methods}
\alias{dwiMixtensor}
\alias{dwiMixtensor-methods}
\alias{dwiMixtensor,ANY-method}
\alias{dwiMixtensor,dtiData-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Methods for Function `dwiMixtensor' in Package `dti' }
\description{
 The method estimates, in each voxel, a mixture of radial symmetric tensors from the DWI data contained in an object of class \code{"dtiData"}.
}
\usage{
  \S4method{dwiMixtensor}{dtiData}(object, maxcomp=2, p=2, maxneighb=7, method="mixtensor", reltol=1e-8, maxit=5000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ Object of class \code{"dtiData"}}
  \item{maxcomp}{ Maximal number of mixture components. }
  \item{p}{ Exponent in Jian-Model (only effective if \code{method="Jian"}. }
  \item{maxneighb}{ Size of neighborhood used to search for minima of
   the ADC. These are then used to determine initial parameters for the 
   mixture model} 
  \item{method}{ Specifies the mixture model used. \code{method="mixtensor"}   specifies a mixture of tensor models (default), \code{method="Jian"} refers to the model defined in Jian et al. (2007) with fixed \code{p} while with \code{method="Jian2"} also \code{p} is estimated. } 
  \item{reltol}{ Relative tolerance for Rs nmmin() function. }
  \item{maxit}{ Maximal number of iterations in Rs nmmin() function. }
}
\details{
  The method estimates, in each voxel, a mixture of radial symmetric tensors from the DWI data contained in an object of class \code{"dtiData"}. The number of mixture components is selected depending on the data, with a maximum number of components specified by \code{maxcomp}. 

  If \code{method="Jian"} the model parameters from Jian et al. (2007) with fixed \code{p} are estimated. With \code{method="Jian2"} also \code{p} is estimated. 

  In all models tensors are restricted to be rotational sysmmetric with common excentricity and destinct largest eigenvalue. 
}
\value{
  An object of class \code{"dwiMixtensor"}.
}
\references{ 
  Jian et al. (2007), A novel tensor distribution model for the 
  diffusion-weighted MR signal, \emph{NeuroImage} \bold{37}, 164--176.
}
\author{   Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{   \code{\link{dtiData}}, 
  \code{\link{readDWIdata}}, 
  \code{\link{medinria}}, 
  \code{\linkS4class{dtiData}}, 
  \code{\linkS4class{dwiMixtensor}}
}
\examples{
  \dontrun{demo(dti_art)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ manip }
\keyword{ models }
