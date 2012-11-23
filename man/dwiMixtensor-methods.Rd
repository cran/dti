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
\S4method{dwiMixtensor}{dtiData}(object, maxcomp=3, method="mixtensor", reltol=1e-6,
          maxit=5000, ngc=1000, optmethod="BFGS", nguess=100*maxcomp^2,
          msc="BIC", pen=NULL, code = "C", thinit = NULL, mc.cores =
          getOption("mc.cores", 1L))
\S4method{dwiMtImprove}{dwiMixtensor,dtiData}(mtobj, dwiobj, maxcomp=3, method="mixtensor", reltol=1e-6, maxit=5000, ngc=1000, optmethod="BFGS", nguess= 100*maxcomp^2, msc="BIC", pen=NULL, where=NULL, new=FALSE)
\S4method{dwiMtCombine}{dwiMixtensor,dwiMixtensor}(mtobj1,mtobj2, msc="BIC", where=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ Object of class \code{"dtiData"}}
  \item{maxcomp}{ Maximal number of mixture components. }
  \item{method}{ Specifies the mixture model used. \code{method="mixtensor"}   (default) and \code{method="mixtensoriso"} specify a mixture of tensor models without and including an isotropic component. } 
  \item{reltol}{ Relative tolerance for R's optim() function. }
  \item{maxit}{ Maximal number of iterations in R's optim() function. }
  \item{ngc}{ provide information on number of voxel processed, elapsed time and estimated remaining time after \code{ngc} voxel.}
  \item{optmethod}{Optimization method used, currently available are
    \code{optmethod="BFGS"} (default for \code{method="mixtensor"}), \code{optmethod="CG"}, 
    \code{optmethod="Nelder-Mead"} and \code{optmethod="L-BFGS-B"} 
    (only method if \code{method="mixtensoriso"}). For \code{optmethod="BFGS"} and \code{optmethod="CG"} weights are obtained as solutions of a linear
    subproblem, i.e. implicit parameters. Positivity of weights is forced by 
    penalization, see \code{pen}. Analytic gradients are used. For \code{optmethod="Nelder-Mead"} weights are obtained as solutions of a linear problem with positivity constraints using the NNLS code by Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
    1973 JUN 15, and published in the book
    SOLVING LEAST SQUARES PROBLEMS, Prentice-Hall, 1974. For \code{optmethod="L-BFGS-B"} weights are handled as regular parameters.}
   \item{nguess}{number of guesses in search for initial estimates}
   \item{msc}{Criterion used to select the order of the mixture model, either
   \code{BIC} (Bayes Information Criterion) \code{AIC} (Akaike Information Criterion) or \code{AICC} ((Bias-)Corrected Akaike Information Criterion).
   \code{None} may be specified to only correct for under-estimation of variances.}
   \item{pen}{Penalty used in optimization criterion for negative mixture components (weights).}
   \item{code}{If \code{code="C"}, \code{method="mixtensor"} and
               \code{optmethod="BFGS"} a much faster C/FORTRAN-code 
               is used}
   \item{thinit}{optional initial estimate for difference of eigenvalues. If not provided (default) an initial estimate
                 is generated from results in using the tensor model}
   \item{mtobj}{For method \code{"dwiMtImprove"} an initial \code{"dwiMixtensor"}-object.}
   \item{dwiobj}{For method \code{"dwiMtImprove"} the \code{"dwiData"} object corresponding to \code{mtobj}}
   \item{mtobj1}{For method \code{"dwiMtCombine"} an \code{"dwiMixtensor"}-object.}
   \item{where}{Mask of voxel for which  \code{"dwiMtImprove"} or \code{"dwiMtCombine"} 
   should be performed.}
   \item{mtobj2}{For method \code{"dwiMtCombine"} an \code{"dwiMixtensor"}-object obtained from the same \code{"dwiData"} object. The maximum number of components in \code{mtobj2} should preferably  be less or equal to the maximum number of components in \code{mtobj1}.}
   \item{new}{if \code{new=TRUE} the original fit supplied to function \code{dwiMtImprove}
             in argument \code{mtobj} will be discarded. Defaults to \code{new=FALSE}, i.e. 
             the original fit will be kept if it performs better with respect to the specified
             criterion in \code{msc}.}
   \item{mc.cores}{Number of cores to use. Defaults to getOption("mc.cores", 2L). 
   This may be larger than the optimal value. Our experience suggests to use 4-6 cores if available.}
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
