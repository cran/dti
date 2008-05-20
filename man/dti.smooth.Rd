\name{dti.smooth}
\alias{dti.smooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Smoothing of DWI data (Diffusion tensor model)}
\description{
The function provides structural adaptive smoothing for diffusion weighted image data within the context of an diffusion tensor (DTI)
model. It implements smoothing of DWI data using a structural assumption of a local (anisotropic) homogeneous  diffusion tensor model (in case an \code{dtiData}-object is provided). It also 
implements adaptive smoothing of a diffusion tensor using a 
Rimannian metric (in case an \code{dtiTensor}-object is given),
althoug we strictly recommend to use the first variant due to methodological reasons.

}
\usage{
dti.smooth(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ either an object of class \code{dtiData} 
or an object of class \code{dtiTensor}}
  \item{\dots}{ additional parameters
\item{hmax}{Maximal bandwidth}
\item{hinit}{Initial bandwidth (default 1)}
\item{lambda}{Critical parameter (default 20)}
\item{rho}{Regularization parameter for anisotropic vicinities (default 1)}
\item{graph}{Logical: Visualize intermediate results (default FALSE)}
\item{slice}{slice number, determines the slice used in visualization}
\item{quant}{determines \code{minanindex} as corresponding quantile of FA if \code{is.null(minanindex)} }
\item{minanindex}{minimal anisotropy index to use in visualization}
\item{hsig}{bandwidth for presmoothing of variance estimates}
\item{lseq}{ sequence of correction factors for \code{lambda}}
\item{method}{ if \code{model=="linear"} estimates are obtained using a linearization of the tensor model. This was the estimate used in Tabelow et.al. (2008). \code{model=="nonlinear"} uses a 
nonlinear regression model with reparametrization that ensures the tensor to be  positive semidefinite, see  Koay et.al. (2006).}
\item{varmethod}{specifies the method for estimating the error  variance. If \code{varmethod=="replicates"} the error variance is estimated from replicated 
gradient directions if possible. Otherwise an estimate is obtained from the residual sum of squares.}
\item{varmodel}{if \code{varmodel=="global"} a homogeneous variance estimate is assumed and estimated as the median of the local variance estimates.}
\item{rician}{logical: apply a correction for Rician bias. This is still experimental and depends on spatial independence of errors.}
\item{niter}{Maximum number of iterations for tensor estimates using the nonlinear model.}
\item{volseq}{If \code{volseq==TRUE} the sum of location weights is fixed to \code{1.25^k} within iteration \code{k} (does not depend on the actual tensor). Otherwise the ellipsoid of positive location weights is determined by a
bandwidth \code{h_k = 1.25^(k/3)}.}
 }
}
\details{
Effective parameters depend on the class of the supplied object.
We highly recommend to use function \code{dti.smooth} on 
DWI data directly, i.e. on an object of class \code{dtiData},
due to methodological reasons. 
}
\value{
An object of class \code{dtiTensor}.
}
\references{ K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, NeuroImage 39(4), 1763-1773 (2008).

C.G. Koay, J.D. Carew, A.L. Alexander, P.J. Basser and M.E. Meyerand.
\emph{ Investigation of Anomalous Estimates of Tensor-Derived Quantities in Diffusion Tensor Imaging}, Magnetic Resonance in medicine, 2006, 55, 930-936.

\url{http://www.wias-berlin.de/projects/matheon_a3/} }

\author{Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de }
\seealso{ \code{\link{dtiData}}, \code{\link{dtiTensor}}, \code{\link{tensor2medinria}} , \code{\linkS4class{dtiData}}, \code{\linkS4class{dtiIndices}}, \code{\linkS4class{dtiTensor}}}
\examples{
\dontrun{demo(dti_art)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }
