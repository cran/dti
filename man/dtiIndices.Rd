\name{dtiIndices}
\alias{dtiIndices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create estimates of anisotropy indices from diffusion tensors  }
\description{
  The function creates estimates of the fractional anisotropy (FA) and relative anisotropy (RA) indices, the main directions of anisotropy and several statistics used for visualization.
}
\usage{
dtiIndices(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ object of class \code{dtiTensor} }
  \item{\dots}{ Additional arguments specifying the model:
   currently none.
}
}
\value{an object of class \code{dtiIndices}.
}
\references{ K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, NeuroImage 39(4), 1763-1773 (2008).

\url{http://www.wias-berlin.de/projects/matheon_a3/}}
\author{Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de }
\seealso{ \code{\link{dtiTensor}}, \code{\linkS4class{dtiTensor}}, \code{\linkS4class{dtiIndices}}}
\examples{
\dontrun{demo(dti_art)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
\keyword{ models }% __ONLY ONE__ keyword per line
