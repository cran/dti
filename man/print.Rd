\name{print}
\alias{print}
\alias{print,dti-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ generic function for objects of class "dti", "dtiData",
"dtiTensor" and "dtiIndices". }
\description{
  The function provides information on data dimensions, data source and existing slot-names for objects of class \code{dti}, \code{dtiData},
\code{dtiTensor} and \code{dtiIndices}.
}
\usage{
print(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class \code{dti}, \code{dtiData},
\code{dtiTensor} or \code{dtiIndices}}
  \item{\dots}{ without effect }
}
\value{\code{invisible(NULL)}
}
\references{ K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, NeuroImage 39(4), 1763-1773 (2008).

\url{http://www.wias-berlin.de/projects/matheon_a3/} }
\author{ Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de }
\seealso{ \code{\link{dtiIndices}}, \code{\link{dtiData}}, \code{\link{dtiTensor}}}
\examples{
\dontrun{demo(dti_art)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
