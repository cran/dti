\name{combineDWIdata}
\alias{combineDWIdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Combine two objects of class "dtiData")
}
\description{
This function creates a dtiData-object from two compatible 
dtiData-objects. Compatible means that the spatial dimensions
coincide, but gradients and b-values may be different. 
}
\usage{
combineDWIdata(x1, x2, s0strategy = "first")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{Object of class \code{"dtiData"}
}
  \item{x2}{Object of class \code{"dtiData"}
}
  \item{s0strategy}{Character, determines
  how the unweighted S0 images are handled. Six strategies are implemented.
  \code{s0strategy="first"} copies the S0 images from object \code{x1},
  \code{s0strategy="second"} copies the S0 images from object \code{x2},
  \code{s0strategy="both"} used the S0 images from both objects.
  \code{s0strategy="rfirst"} creates one average S0 image from object \code{x1},
  \code{s0strategy="rsecond"} creates one average S0 image from object \code{x2},
  \code{s0strategy="rboth"} creates one average S0 image from the S0 images in both objects.
}
}
\details{
  The function can be used to merge two objects of class "dtiData" under the condition that
  the information in slot \code{ddim} in both objects is identical. Also slots \code{voxelext},
  \code{orientation} and \code{rotation} should be indentical. 
}
\value{
  An object of class "dtiData".
}
\author{
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{dtiData}}, 
  \code{\link{readDWIdata}}, 
  \code{\linkS4class{dtiData}}, 
  \code{\link{subsetg}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
