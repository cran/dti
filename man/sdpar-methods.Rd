\name{sdpar-methods}
\docType{methods}
\alias{sdpar}
\alias{sdpar-methods}
\alias{sdpar,ANY-method}
\alias{sdpar,dtiData-method}
\title{ Methods for Function `sdpar' in Package `dti' }
\description{
This function estimates the parameters of a piecewise linear model for 
the dependence between error standard deviation and mean. 
  }
\usage{
  \S4method{sdpar}{dtiData}(object,level=NULL,sdmethod="sd",interactive=TRUE,threshfactor=1)
}
\arguments{
\item{object}{ An object of class \code{dtiData} }
\item{level}{ Suggested value for slot \code{level}. As a default the value in \code{object@level} is used. The value determines the lower endpoint of the linear section in the model for error standard deviation as a function of the mean.}
\item{sdmethod}{Method for estimating voxelwise standard deviations if replicates of zero weighted images are available, either \code{"sd"} or \code{"mad"}.}
\item{interactive}{If \code{TRUE} a density of values in zero weighted images is plotted together with 
the specification of the lower endpoint of the interval of linearity. A good choice of this point should correspond, if present, to the minimum between the first two modes of the density estimate. The value can be changed or accepted. If changed a new value for slot \code{lambda} is set. }
\item{threshfactor}{Factor for threshold-value selected if function is run in interactive mode. May be used to correct results if automatic threshold selection fails.  }
}
\section{Methods}{
\describe{

  \item{obj = "ANY"}{ Returns a warning }

  \item{obj = "dtiData"}{ Estimate parameters of a model for the dependence between error standard deviation and mean.}

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
  \code{\linkS4class{dtiTensor}},
}
\examples{
  \dontrun{demo(dti_art)}
}
\keyword{ methods }
