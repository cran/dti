\name{tensor2medinria}
\alias{tensor2medinria}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Write Tensor-Data as NIFTI File }
\description{
  Write Tensor-Data as NIFTI File. Interface function to MedINRIA
}
\usage{
tensor2medinria(obj, filename, xind = NULL, yind = NULL, zind = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ object of class \code{dtiTensor} }
  \item{filename}{ Filename, without extension 'nii' }
  \item{xind}{ index to define a subcube in x-direction. If \code{is.null(xind} all voxel indices are used. }
  \item{yind}{ index to define a subcube in y-direction. If \code{is.null(yind} all voxel indices are used. }
  \item{zind}{ index to define a subcube in z-direction. If \code{is.null(zind} all voxel indices are used. }
}
\details{
The function allows to write tensor data as NIFTY-Files readable by MedINRIA. 
}
\value{\code{NULL}
}
\references{ K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, NeuroImage 39(4), 1763-1773 (2008).

\url{http://www.wias-berlin.de/projects/matheon_a3/} 

P. Fillard, J. Souplet and N. Toussaint \emph{ Medical Image Navigation and Research
Tool by INRIA (MedINRIA)}, INRIA Sophia Antipolis - Research Project ASCLEPIOS 2007}
\author{ Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de }
\seealso{ \code{\link{dtiTensor}} }
\examples{
\dontrun{demo(dti_art)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
