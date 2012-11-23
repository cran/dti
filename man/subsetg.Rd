\name{subsetg}
\alias{subsetg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create an objects of class "dtiData" containing only a subset of gradient directions.
}
\description{
This function creates an object of class "dtiData" that containes
only a subset, defined by an index vector, of the S0 and diffusion weighted images. 
This function may e.g. be used to separate information measured on different
shells.
}
\usage{
subsetg(x, ind)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class \code{"dtiData"}
}
  \item{ind}{
Indexvector containing values between \code{1} and \code{x@ngrad}.
}
}
\value{
  An object of class "dtiData".
}
\author{
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{
  \code{\link{dtiData}}, 
  \code{\link{readDWIdata}}, 
  \code{\linkS4class{dtiData}}, 
  \code{\link{combineDWIdata}}
}
\keyword{ manip }
