\name{dtiTensor-methods}
\docType{methods}
\alias{dtiTensor}
\alias{dtiTensor-methods}
\alias{dtiTensor,ANY-method}
\alias{dtiTensor,dtiData-method}
\title{ Methods for Function `dtiTensor' in Package `dti' }
\description{ 
  The method estimates, in each voxel, the diffusion tensor from the DWI data contained in an object of class \code{"dtiData"}.
}
\usage{
  \S4method{dtiTensor}{dtiData}(object, method="nonlinear", varmethod="replicates", 
                                varmodel="local", mc.cores = getOption("mc.cores", 1L))
}
\arguments{
  \item{object}{ Object of class \code{"dtiData"}}
  \item{method}{ Method for tensor estimation. May be \code{"linear"}, or \code{"nonlinear"}. }
  \item{varmethod}{ Specifies the method for estimating the error variance. May be \code{"replicates"}. }
  \item{varmodel}{ Specifies the model for the variance. May be \code{"global"}, or \code{"local"}. }
   \item{mc.cores}{Number of cores to use. Defaults to getOption("mc.cores", 2L). 
   This may be larger than the optimal value. Our experience suggests to use 4-6 cores if available.}
}
\value{ 
  An object of class \code{"dtiTensor"}.
}
\section{Methods}{
\describe{

  \item{obj = "ANY"}{ Returns a warning. }

  \item{obj = "dtiData"}{ Estimate diffusion tensor from data in each voxel with the different options for the regression type and model for variance estimation.  If \code{method=="linear"} estimates are obtained using a linearization of the tensor model. This was the estimate used in Tabelow et.al. (2008).  \code{method=="nonlinear"} uses a nonlinear regression model with reparametrization that ensures the tensor to be  positive semidefinite, see  Koay et.al. (2006).  If \code{varmethod=="replicates"} the error variance is estimated from replicated gradient directions if possible, otherwise an estimate is obtained from the residual sum of squares.  If \code{varmodel=="global"} a homogeneous variance is assumed and estimated as the median of the local variance estimates. }

}
}
\references{ J. Polzehl and K. Tabelow, 
             \emph{Beyond the diffusion tensor model: The package \pkg{dti}},
             Journal of Statistical Software, to appear.

             K. Tabelow, H.U. Voss and J. Polzehl, 
             \emph{Modeling the orientation distribution function by 
             mixtures of angular central Gaussian distributions},
             Journal of Neuroscience Methods, to appear.

             J. Polzehl and K. Tabelow, 
             \emph{Structural adaptive smoothing in diffusion tensor imaging: 
             The R package dti}, Journal of Statistical Software, 31 (2009) pp. 1--24.
             
              K. Tabelow, J. Polzehl, V. Spokoiny and H.U. Voss. 
              \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, 
              NeuroImage 39(4), 1763-1773 (2008).

              C.G. Koay, J.D. Carew, A.L. Alexander, P.J. Basser and M.E. Meyerand.
              \emph{ Investigation of Anomalous Estimates of Tensor-Derived Quantities in Diffusion Tensor Imaging}, 
              Magnetic Resonance in Medicine, 2006, 55, 930-936.

              \url{ http://www.wias-berlin.de/projects/matheon_a3/} 
}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\link{dtiData}}, 
  \code{\link{readDWIdata}}, 
  \code{\link{dtiIndices-methods}},  
  \code{\link{medinria}}, 
  \code{\linkS4class{dtiData}}, 
  \code{\linkS4class{dtiTensor}}
  \code{\linkS4class{dwiMixtensor}}
}
\examples{
  \dontrun{demo(dti_art)}
}
\keyword{ methods }
\keyword{ manip }
\keyword{ models }
