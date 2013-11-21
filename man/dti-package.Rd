\name{dti-package}
\alias{dti-package}
\alias{dti}
\docType{package}
\title{
Analysis of diffusion weighted imaging (DWI) data
}
\description{
Diffusion Weighted Imaging (DWI) is a Magnetic Resonance Imaging
             modality, that measures diffusion of water in tissues like the human 
             brain. The package contains R-functions to process diffusion-weighted 
             data. The functionality includes diffusion tensor imaging (DTI),
             diffusion kurtosis imaging (DKI), modeling for high angular resolution 
             diffusion weighted imaging (HARDI) using Q-ball-reconstruction and 
             tensor mixture models, several methods for structural adaptive 
             smoothing including POAS and msPOAS, and a streamline fiber tracking 
             for tensor and tensor mixture models.
             The package provides functionality to manipulate and visualize results 
             in 2D and 3D.
}
\details{
\tabular{ll}{
Package: \tab dti\cr
Version: \tab 1.1-7\cr
Date: \tab 2013-11-21\cr
Authors@R: \tab c(person("Karsten", "Tabelow", role = c("aut", "cre"),
                    email = "karsten.tabelow@wias-berlin.de"), 
             person("Joerg", "Polzehl", role = c("aut"),
                    email = "joerg.polzehl@wias-berlin.de"),
             person("Felix", "Anker", role = c("ctb")))\cr
Depends: \tab R (>= 2.14.0), awsMethods (>= 1.0-1)\cr
Imports: \tab methods, parallel, oro.nifti (>= 0.3.9), oro.dicom, adimpro, gsl, rgl\cr
Suggests: \tab quadprog\cr
LazyData: \tab TRUE\cr
License: \tab GPL (>= 2)\cr
Copyright: \tab This package is
Copyright (C) 2005-2013 Weierstrass Institute for
Applied Analysis and Stochastics.\cr
URL: \tab http://www.wias-berlin.de/projects/matheon_a3\cr
Built: \tab R 3.0.2; x86_64-unknown-linux-gnu; 2013-11-25 15:38:06 UTC; unix\cr
}

Index:
\preformatted{
awssigmc                Estimate noise variance for multicoil MR
                        systems
colqFA                  FA map color scheme
combineDWIdata          Combine two objects of class "dtiData")
dkiTensor-methods       Diffusion Kurtosis Imaging (DKI)
dti-package             DTI Analysis
dti.smooth-methods      Methods for Function 'dti.smooth' in Package
                        'dti'
dtiIndices-methods      Methods for Function 'dtiIndices' in Package
                        'dti'
dtiTensor-methods       Methods for Function 'dtiTensor' in Package
                        'dti'
dwi-class               Class "dwi"
dwi.smooth-methods      Smooth DWI data
dwiMixtensor-methods    Methods for Function 'dwiMixtensor' in Package
                        'dti'
dwiQball-methods        Methods for Function 'dwiQball' in Package
                        'dti'
dwiRiceBias-methods     Correction for Rician Bias
dwiSqrtODF-methods      Methods for positive definite EAP and ODF
                        estimation in Package 'dti'
extract-methods         Methods for Function 'extract' and '[' in
                        Package 'dti'
getmask-methods         Methods for Function 'getmask' in Package 'dti'
getsdofsb-methods       Estimate the noise standard deviation
medinria                Read/Write Diffusion Tensor Data from/to NIFTI
                        File
optgrad                 Optimal gradient directions
optgradients            Optimal gradient directions for number of
                        gradients between 6 and 162
plot-methods            Methods for Function 'plot' in Package 'dti'
pmatrix                 Parallelize columnwise computations on a matrix
                        using forking (generalizes function 'pvec' from
                        package 'parallel')
polyeder                Polyeders derived from the Icosahedron (icosa0)
                        by sequential triangulation of surface
                        triangles
print-methods           Methods for Function 'print' in Package 'dti'
readDWIdata             Read Diffusion Weighted Data
sdpar-methods           Methods for Function 'sdpar' in Package 'dti'
show-methods            Methods for Function 'show' in Package 'dti'
show3d-methods          Methods for Function 'show3d' in Package 'dti'
showFAColorScale        Writes an image with the colqFA colorscale to
                        disk.
subsetg                 Create an objects of class "dtiData" containing
                        only a subset of gradient directions.
summary-methods         Methods for Function 'summary' in Package 'dti'
tracking-methods        Methods for Function 'tracking' in Package
                        'dti'
}
}
\author{
Karsten Tabelow and Joerg Polzehl, with contributions from F. Anker 

Maintainer: Karsten Tabelow <karsten.tabelow@wias-berlin.de>
}
\references{S. Becker, K. Tabelow, H.U. Voss, A. Anwander, R.M. Heidemann and
            J. Polzehl, \emph{Position-orientation adaptive smoothing of 
            diffusion weighted magnetic resonance data (POAS)}, Medical Image
            Analysis, 16 (2012), pp. 1142-1155.

            J. Polzehl and K. Tabelow, 
             \emph{Beyond the diffusion tensor model: The package \pkg{dti}},
             Journal of Statistical Software, 44 no. 12 (2011) pp. 1-26.

             K. Tabelow, H.U. Voss and J. Polzehl, 
             \emph{Modeling the orientation distribution function by 
             mixtures of angular central Gaussian distributions},
             Journal of Neuroscience Methods, 203 (2012), pp. 200-211.

             J. Polzehl and K. Tabelow, 
             \emph{Structural adaptive smoothing in diffusion tensor imaging: 
             The R package dti}, Journal of Statistical Software, 31 (2009) pp. 1--24.
             
              K. Tabelow, J. Polzehl, V. Spokoiny and H.U. Voss. 
              \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, 
              NeuroImage 39(4), 1763-1773 (2008).

              \url{ http://www.wias-berlin.de/projects/matheon_a3/} 
}
\keyword{ package }
\examples{
  \dontrun{demo(dti_art)}
  \dontrun{demo(mixtens_art)}
}  
