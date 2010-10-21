\name{dti-package}
\alias{dti-package}
\alias{dti}
\docType{package}
\title{ DTI Analysis }
\description{
Diffusion Weighted Imaging (DWI) is a Magnetic Resonance Imaging
             modality, that measures diffusion of water in tissues like the human brain.
             The package contains R-functions to process diffusion-weighted data. The functionality 
             includes diffusion tensor imaging (DTI), structural adaptive smoothing in 
             in case of (DTI) (K. Tabelow, J. Polzehl, V. Spokoiny,
             and H.U. Voss, Diffusion Tensor Imaging: Structural Adaptive Smoothing, 
             Neuroimage 39(4), 1763-1773 (2008)), modeling for high angular resolution 
             diffusion weighted imaging (HARDI) using Q-ball-reconstruction and tensor mixture
             models and a stremaline fiber tracking for tensor and tensor mixture models.
             The package provides functionality to manipulate and visualize results in 2D and 3D.}
\details{
\tabular{ll}{
Package: \tab dti\cr
Version: \tab 0.9-3\cr
Date: \tab 2010-10-21\cr
Depends: \tab R (>= 2.5.0), adimpro, fmri, rgl\cr
License: \tab GPL (>=2)\cr
Copyright: \tab 2008-2010 Weierstrass Institute for\cr
\tab Applied Analysis and Stochastics.\cr
URL: \tab \url{http://www.wias-berlin.de/projects/matheon_a3}\cr
}

The package is based on S4 classes and methods.  For help on a specific topic use \code{class ? <class-name>} for classes, \code{methods ? <method-name>} for methods and \code{?<function-name>} for all other functions. 

Index:
\preformatted{
dti-class               Classes "dti", "dtiData", "dtiTensor",
                        "dtiIndices", "dwiQball", "dwiMixtensor"
                        and "dwiFiber"

dtiData                 Read Diffusion Weighted Data
                        from Image File
readDWIdata             Read Diffusion Weighted Data
                        from Directory

dti.smooth              Structural Adaptive Smoothing

dtiTensor-methods       Methods for Function 'dtiTensor'
dtiIndices-methods      Methods for Function 'dtiIndices'

dwiQball-methods        Methods for Function 'dwiQball'
dwiMixtensor-methods    Methods for Function 'dwiMixtensor'

extract-methods         Methods for Functions 'extract' and "["
plot-methods            Methods for Function 'plot'
show-methods            Methods for Function 'show'
show3d-methods          Methods for Function 'show3d'
print-methods           Methods for Function 'print'
summary-methods         Methods for Function 'summary'

tracking-methods        Methods for Functions 'tracking', 'selectFibers' and 'reduceFibers'
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
