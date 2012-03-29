\name{awssigmc}
\alias{awssigmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate noise variance for multicoil MR systems
}
\description{
    The distribution of image intensity values \eqn{S_i} divided by the noise standard deviation in \eqn{K}-space \eqn{\sigma} 
    in dMRI experiments is assumed
    to follow a non-central chi-distribution with \eqn{2L} degrees of freedom and noncentrality parameter \eqn{\eta}, where \eqn{L} refers to the number of receiver 
    coils in the system and \eqn{\sigma \eta} is the signal of interest. This is an idealization in the sense that 
    each coil is assumed to have the same contribution at each location. For realistic modeling \eqn{L} should 
    be a locally smooth function in voxel space that reflects the varying local influence of the receiver coils in the
    the reconstruction algorithm used. 
  
    The function assumes \eqn{L} to be known and estimates a global \eqn{\sigma} employing an assumption of local homogeneity for
     the noncentrality parameter \eqn{\eta}.
}
\usage{
awssigmc(y, steps, mask = NULL, ncoils = 1, vext = c(1, 1), lambda = 10, h0 = 2, method = "median", verbose = FALSE, model = "chisq", sequence = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
    3D array, usually obtained from an object of class \code{dwi} as
    \code{obj@si[,,,i]} for some \code{i}, i.e. one 3D image from an dMRI experiment.
}
  \item{steps}{
   number of steps in adapive weights smoothing, used to reveal the unerlying
   mean structure.
}
  \item{mask}{
   restrict computations to voxel in mask, if \code{is.null(mask)} all voxel are used.
}
  \item{ncoils}{
   number of coils, or equivalently number of effective degrees of freedom of non-central chi distribution
   divided by 2.
}
  \item{vext}{
   voxel extentions
}
  \item{lambda}{
   scale parameter in adaptive weights smoothing
}
  \item{h0}{
   initial bandwidth
}
  \item{method}{
   either \code{"median"} or \code{"mode"}.
}
  \item{verbose}{
   if \code{verbose==TRUE} density plots 
   and quantiles of local estimates of \code{sigma} are provided.
}
  \item{model}{
   either \code{"chi"} or \code{"chisq"}. In the latter case 
   smoothing and variance estimation are performed for \code{y^2}
   instead of \code{y} which is considerably faster. 
}
  \item{sequence}{
   if \code{sequence=TRUE} a vector of estimates for the noise
   standard deviation \code{sigma} for the individual steps is returned
   instead of the final value only.  
}
}
\value{
   either a scalar or a vector of estimated noise standard deviations.
}
\references{
K. Tabelow and J. Polzehl (2013). Estimating the noise level in 
MRI using structural adaptive smoothing. Manuscript in preparation.
}
\author{
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\keyword{ smooth }
