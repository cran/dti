%- Also NEED an '\alias' for EACH other topic documented here.
\name{getmask-methods}
\docType{methods}
\alias{getmask}
\alias{getmask-methods}
\alias{getmask,ANY-method}
\alias{getmask,dtiData-method}
\alias{getmask,array-method}
\title{ Methods for Function `getmask' in Package `dti' }
\description{
Create a mask containing voxel inside the head
}
\usage{
\S4method{getmask}{dtiData}(object, level = NULL, prop = 0.4, size = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{"dtiData"}}
  \item{level}{S0 intensity value to be used to discriminate between
       voxel inside and outside the brain. A good value of level may be determined
       using method \code{sdpar} in advance.}
  \item{prop}{proportion of voxel in test area with s0 value larger than level
    needed to decide for a voxel inside the brain}
  \item{size}{size of a cube defining a test area}
}
\section{Methods}{
\describe{

  \item{obj = "ANY"}{ Returns a warning }

  \item{obj = "dtiData"}{Create a mask containing voxel inside the head }

  \item{obj = "array"}{Create a mask containing voxel inside the head }
}
}
\value{
  The function returns an object of class \code{dtiData}.
}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\link{dtiData}},
  \code{\link{readDWIdata}},
  \code{\link{dti.smooth}}, 
  \code{\link{sdpar}}
}
\keyword{ methods }
