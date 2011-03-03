\name{dwiMixtensor-methods}
\docType{methods}
\alias{dwiMixtensor}
\alias{dwiMixtensor-methods}
\alias{dwiMixtensor,ANY-method}
\alias{dwiMixtensor,dtiData-method}
\alias{dwiMtImprove}
\alias{dwiMtImprove-methods}
\alias{dwiMtImprove,ANY-method}
\alias{dwiMtImprove,dwiMixtensor,dtiData-method}
\alias{dwiMtCombine}
\alias{dwiMtCombine-methods}
\alias{dwiMtCombine,ANY-method}
\alias{dwiMtCombine,dwiMixtensor,dwiMixtensor-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Methods for Function `dwiMixtensor' in Package `dti' }
\description{
 The method estimates, in each voxel, a mixture of radial symmetric tensors from the DWI data contained in an object of class \code{"dtiData"}.
}
\usage{
\S4method{dwiMixtensor}{dtiData}(object, maxcomp=3,  p=40, method="mixtensor", reltol=1e-6,
          maxit=5000, ngc=1000, optmethod="BFGS",nguess=100*maxcomp^2,
          msc="BIC",pen=NULL)
\S4method{dwiMtImprove}{dwiMixtensor,dtiData}(mtobj, dwiobj, maxcomp=3,  p=40, method="mixtensor", reltol=1e-6, maxit=5000,ngc=1000, optmethod="BFGS", nguess=100*maxcomp^2,msc="BIC",pen=1,where=NULL)
\S4method{dwiMtCombine}{dwiMixtensor,dwiMixtensor}(mtobj1,mtobj2, msc="BIC", where=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ Object of class \code{"dtiData"}}
  \item{maxcomp}{ Maximal number of mixture components. }
  \item{p}{ Exponent in Jian-Model (only effective if \code{method="Jian"}(Not yet implemented) }
  \item{method}{ Specifies the mixture model used. \code{method="mixtensor"}   specifies a mixture of tensor models (default), \code{method="Jian"} refers to the model defined in Jian et al. (2007) with fixed \code{p}. For method \code{dwiMtImprove} also \code{method="mixtensoriso"} may be specified to include an isotropic component in the model.} 
  \item{reltol}{ Relative tolerance for R's nmmin() function. }
  \item{maxit}{ Maximal number of iterations in Rs nmmin() function. }
  \item{ngc}{ provide information on number of voxel processed, elapsed time and estimated remaining time after \code{ngc} voxel.}
  \item{optmethod}{Optimization method used, currently available is 
   \code{optmethod="Nelder-Mead"}, \code{optmethod="BFGS"} using analytic gradients
   will be added}
   \item{nguess}{number of guesses in search for initial estimates}
   \item{msc}{Criterion used to select the order of the mixture model, either
   \code{BIC} (Bayes Information Criterion) \code{AIC} (Akaike Information Criterion) or \code{AICC} ((Bias-)Corrected Akaike Information Criterion)}
   \item{pen}{Penalty used in optimization criterion for negative mixture components.}
   \item{mtobj}{For method \code{"dwiMtImprove"} an initial \code{"dwiMixtensor"}-object.}
   \item{dwiobj}{For method \code{"dwiMtImprove"} the \code{"dwiData"} object corresponding to \code{mtobj}}
   \item{mtobj1}{For method \code{"dwiMtCombine"} an \code{"dwiMixtensor"}-object.}
   \item{where}{Mask of voxel for which  \code{"dwiMtImprove"} or \code{"dwiMtCombine"} 
   should be performed.}
   \item{mtobj2}{For method \code{"dwiMtCombine"} an \code{"dwiMixtensor"}-object obtained from the same \code{"dwiData"} object. The maximum number of components in \code{mtobj2} should preferably  be less or equal to the maximum number of components in \code{mtobj1}.}
   
   }
\details{
  The method \code{"dwiMixtensor"} estimates, in each voxel, a mixture of radial symmetric tensors from the DWI data contained in an object of class \code{"dtiData"}. The number of mixture components is selected depending on the data, with a maximum number of components specified by \code{maxcomp}. 

  If \code{method="Jian"} the model parameters from Jian et al. (2007) with fixed \code{p} are estimated. With \code{method="Jian2"} also \code{p} is estimated. 

  In a voxel tensors are restricted to be rotational sysmmetric with common  excentricity and destinct largest eigenvalue. 
  
  The method \code{"dwiMtImprove"} evaluates the results in \code{mtobj}, including directions identified in neighboring voxel, to obtain an alternative set of initial values for estimating the parameters of the tensor mixture model. The resulting object contains, in a voxel, either the results from \code{mtobj} or the estimated parameters 
  obtained by optimisation starting with the new initial values, depending on a comparison
  of their respective estimated MSEP values as specified by \code{msc}. 
  The specification of \code{method="mixtensoriso"} includes an isotropic component 
  into the model. If the isotropic term leads to an improbement with respect to the estimated MSEP the sum of weights for such a voxel will be less than 1, with 
  the discrepancy to 1 corresponding to the partial volume associated to the isotropic
  compartment.
  
  The method \code{"dwiMtCombine"} enables to combine results obtained for the same
  dwi data set with different specifications, e.g. for maximum number of components
  \code{mcomp} and settings that influence initial estimates. The combined result 
  contains in each voxel the best result from both reconstructions with respect to
  the specified model selection criterion \code{msc}. 
  
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
  \dontrun{demo(mixtens_art)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ manip }
\keyword{ models }
