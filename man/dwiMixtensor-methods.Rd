\name{dwiMixtensor-methods}
\docType{methods}
\alias{dwiMixtensor}
\alias{dwiMixtensor-methods}
\alias{dwiMixtensor,ANY-method}
\alias{dwiMixtensor,dtiData-method}
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
\S4method{dwiMixtensor}{dtiData}(object, maxcomp=3,        
          model=c("MT","MTiso","MTisoFA","MTisoEV"), fa=NULL, 
          lambda=NULL, reltol=1e-6, maxit=5000, ngc=1000, 
          nguess=100*maxcomp^2, msc=c("BIC","AIC","AICC","none"), 
          mc.cores = setCores(,reprt=FALSE))
\S4method{dwiMtCombine}{dwiMixtensor,dwiMixtensor}(mtobj1,mtobj2, msc="BIC", where=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ Object of class \code{"dtiData"}}
  \item{maxcomp}{ Maximal number of mixture components. }
  \item{model}{ Specifies the mixture model used. \code{"MT"} corresponds to a mixture 
   of prolate tensors, \code{"MTiso"} includes an isotropic compartment, \code{"MTisoFA"} additionally fixes FA to the value given in argument \code{fa} and \code{"MTisoEV"} 
   uses eigenvalues specified by \code{fa} and \code{lambda}. } 
  \item{fa}{Value for FA in case of \code{model="MTisoFA"} or \code{model="MTisoEV"}}
  \item{lambda}{Value for first eigenvalue in case of \code{model="MTisoEV"}}
  \item{reltol}{ Relative tolerance for R's optim() function. }
  \item{maxit}{ Maximal number of iterations in R's optim() function. }
  \item{ngc}{ provide information on number of voxel processed, elapsed time and estimated remaining time after \code{ngc} voxel.}
   \item{nguess}{number of guesses in search for initial estimates}
   \item{msc}{Criterion used to select the order of the mixture model, either
   \code{BIC} (Bayes Information Criterion) \code{AIC} (Akaike Information Criterion) or \code{AICC} ((Bias-)Corrected Akaike Information Criterion).
   \code{None} may be specified to only correct for under-estimation of variances.}
   \item{mtobj1}{For method \code{"dwiMtCombine"} an \code{"dwiMixtensor"}-object.}
   \item{where}{Mask of voxel for which  \code{"dwiMtImprove"} or \code{"dwiMtCombine"} 
   should be performed.}
   \item{mtobj2}{For method \code{"dwiMtCombine"} an \code{"dwiMixtensor"}-object obtained from the same \code{"dwiData"} object. The maximum number of components in \code{mtobj2} should preferably  be less or equal to the maximum number of components in \code{mtobj1}.}
   \item{mc.cores}{Number of cores to use. Defaults to number of threads specified for openMP, see documentation of package \pkg{awsMethods}. 
   Our experience suggests to use 4-6 cores if available.}
   }
\details{
  For \code{model=="MT"} the function estimates, in each voxel, a mixture of radial symmetric (prolate) tensors from the DWI data contained in an object of class \code{"dtiData"}. The number of mixture components is selected depending on the data, with a maximum number of components specified by \code{maxcomp}. Optimization is performed usin R's internal BFGS code with mixture weights (volumes of compartments
  corresponding to a tensor component) computed using the Lawson-Hannson NNLS code. 
  In case of \code{model=="MTiso"} the model additionally contains an isotropic compartment. Optimization uses the internal L-BFGS-B code. 
  \code{model=="MTisoFA"} and \code{model=="MTisoEV"} fix FA and eigenvalues 
  of the prolate tensors, respectively, in the tensor mixture model with isotropic compartment.

   
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
