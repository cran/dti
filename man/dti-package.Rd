\name{dti-package}
\alias{dti-package}
\alias{dti}
\docType{package}
\title{ DTI Analysis }
\description{
  Diffusion Weighted Imaging is a Magnetic Resonance Imaging modality, that measures diffusion of water in tissues like the human brain.  The package contains R-functions to process diffusion-weighted data in the context of the diffusion tensor model (DTI).  This includes the calculation of anisotropy measures and, most important, the implementation of our structural adaptive smoothing algorithm as described in Tabelow et al. (2008).
}
\details{
\tabular{ll}{
Package: \tab dti\cr
Version: \tab 0.5-7\cr
Date: \tab 2008-07-24\cr
Depends: \tab R (>= 2.5.0), adimpro, fmri, rgl\cr
License: \tab GPL (>=2)\cr
Copyright: \tab 2008 Weierstrass Institute for\cr
\tab Applied Analysis and Stochastics.\cr
URL: \tab \url{http://www.wias-berlin.de/projects/matheon_a3}\cr
}

The package is based on S4 classes and methods.  For help on a specific topic use \code{class ? <class-name>} for classes, \code{methods ? <method-name>} for methods and \code{?<function-name>} for all other functions. 

Index:
\preformatted{
dti-class               Classes "dti", "dtiData", "dtiTensor",
                        "dtiIndices"

dtiData                 Read Diffusion Weighted Data
                        from Image File
readDWIdata             Read Diffusion Weighted Data
                        from Directory

dti.smooth              Structural Adaptive Smoothing

dtiTensor-methods       Methods for Function 'dtiTensor'
dtiIndices-methods      Methods for Function 'dtiIndices'

extract-methods         Methods for Function 'extract'
plot-methods            Methods for Function 'plot'
show3d-methods          Methods for Function 'show3d'
print-methods           Methods for Function 'print'
summary-methods         Methods for Function 'summary'

tensor2medinria         Write Tensor Data as NIFTI File
medinria2tensor         Read Tensor Data from NIFTI File
}

}
\author{
Karsten Tabelow \email{tabelow@wias-berlin.de}\cr
Joerg Polzehl \email{polzehl@wias-berlin.de}

Maintainer: Karsten Tabelow \email{tabelow@wias-berlin.de}
}
\references{  Tabelow, K. and Polzehl, J. and Spokoiny, V. and Voss, H. U. (2008)
              Diffusion Tensor Imaging: Structural adaptive smoothing,
              \emph{NeuroImage} \bold{39(4)}, 1763--1773.

              Polzehl, J. and Tabelow, K. (2008)
              Structural Adaptive Smoothing in Diffusion Tensor Imaging:  the R package dti,
              WIAS-Preprint 1382.

              \url{ http://www.wias-berlin.de/projects/matheon_a3/} 
}
\examples{
  \dontrun{demo(dti_art)}
}
\keyword{ package }
