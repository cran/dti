\name{dwi-class}
\docType{class}
\alias{dwi-class}
\alias{dtiData-class}
\alias{dtiTensor-class}
\alias{dwiMixtensor-class}
\alias{dtiIndices-class}
\alias{dwiQball-class}
\alias{dwiFiber-class}
\title{Class "dwi" }
\description{ The family of \code{"dwi"} classes is used for Diffusion Weighted Imaging (DWI) data and, within the Diffusion Tensor Model (DTI), diffusion tenors and its indices. }
\section{Objects from the Class}{
  \code{"dwi"} is only a superclass, no instances should be created. However, objects can be created by calls of the form \code{new("dwi", ...)}. \code{"dtiData"}, \code{"dtiTensor"}, and \code{"dtiIndices"} can be created from their correspondingly named functions and methods.
}
\section{Slots}{
\describe{
    \item{ \code{.Data}: }{ Object of class \code{"list"}, usually empty. }
    \item{ \code{gradient}: }{ Object of class \code{"matrix"}, matrix of dimension \code{c(3,ngrad)} containing gradient directions. }
    \item{ \code{btb}: }{ Object of class \code{"matrix"}, matrix of dimension \code{c(6,ngrad)} obtained from gradient directions. }
    \item{ \code{ngrad}: }{ Object of class \code{"integer"}, number of gradients (including zero gradients). }
    \item{ \code{s0ind}: }{ Object of class \code{"integer"}, index of zero gradients within the sequence \code{1:ngrad}. }
    \item{ \code{replind}: }{ Object of class \code{"integer"}, index (identifier) of unique gradient directions. Used to characterize replications in the gradient design by identical indices. length is \code{ngrad}. }
    \item{ \code{ddim}: }{ Object of class \code{"integer"}, dimension of subcube defined by \code{xind}, \code{yind} and \code{zind}. }
    \item{ \code{ddim0}: }{ Object of class \code{"integer"}, dimension of original image cubes. Vector of length 3. }
    \item{ \code{xind}, \code{yind}, \code{zind}:}{ Objects of class \code{"integer"}, index for subcube definition in x-, y- and z-direction. }
    \item{ \code{voxelext}: }{ Object of class \code{"numeric"}, voxel extensions in x-, y- and z-direction. Vector of length 3. }
    \item{ \code{orientation}: }{ Object of class \code{"integer"}, orientation of data according to AFNI convention. Vector of length 3. }
    \item{ \code{rotation}: }{ Object of class \code{"matrix"}, optional rotation matrix for gradient directions. }
    \item{ \code{level}: }{ Object of class \code{"numeric"}, minimal valid S0-level. No evaluation will be performed for voxels with S0-values less than \code{level}. } 
    \item{ \code{source}: }{ Object of class \code{"character"}, name of the source imgage file or source directory. }
    \item{ \code{call}: }{ Object of class \code{"call"}, call that created the object. }
}

For class \code{"dtiData"}:

\describe{
    \item{ \code{si}: }{ Object of class \code{"array"}, Diffusion Weighted Data. }
    \item{ \code{sdcoef}: }{ Object of class \code{"numeric"}, Parameters of the model for error 
      standard deviation as a function of the mean. First two entries refer to intercept and slope of a linear function,
      third and fourth value are the endpoints of the interval of linearity. Contains \code{rep(0,4)} if not set. If the function }
}

For class \code{"dtiTensor"}:

\describe{
    \item{ \code{D}: }{ Object of class \code{"array"}, estimated tensors, dimension \code{c(6,ddim)}.
          Tensors are stored as upper diagonal matrices.} 
    \item{ \code{th0}: }{ Object of class \code{"array"}, estimated intensities in S0 images, dimension \code{ddim}} 
    \item{ \code{sigma}: }{ Object of class \code{"array"}, estimated error variances if \code{method=="linear"}, zero otherwise.} 
    \item{ \code{scorr}: }{ Object of class \code{"numeric"}, estimated spatial correlations in coordinate directions} 
    \item{ \code{bw}: }{ Object of class \code{"numeric"}, bandwidth for a Gaussian kernel that approximately creates the estimated spatial correlations. Needed for adjustments of critical values in the adaptive smoothing algorithm used in function \code{dti.smooth}} 
    \item{ \code{mask}: }{ Object of class \code{"array"}, logical indicating the voxel where the tensor was estimated.} 
    \item{ \code{hmax}: }{ Object of class \code{"numeric"}, maximal bandwidth in case of adaptive smoothing, 1 otherwise.} 
    \item{ \code{outlier}: }{ Object of class \code{"numeric"}, index of voxel where physical constraints are not met, i.e. where the observed values in gradient images Si were larger than the corresponding S0 values. These are probably motion effects or registration errors. Values are replaced by the corresponding (mean) S0 values.}
    \item{\code{scale}:}{Numerical value corresponding to the 95\% quantile of the maximal eigenvalues of estimated tensors within the mask.
          Used for scaling in function \code{show3d.dtiTensor}}
    \item{ \code{method}: }{ Object of class \code{"character"}, either \code{"linear"} or \code{"nonlinear"} or \code{"unknown"}. Indicates the regression model used for estimating the tensors.}
}

For class \code{"dtiIndices"}:

\describe{
    \item{ \code{fa}: }{ Object of class \code{"array"}, Fractional anisotropy values (FA)}
    \item{ \code{ga}: }{ Object of class \code{"array"}, Geodetic anisotropy values (GA)}
    \item{ \code{md}: }{ Object of class \code{"array"}, Mean diffusivity values (MD) }
    \item{ \code{andir}: }{ Object of class \code{"array"}, Main directions of anisotropy }
    \item{ \code{bary}: }{ Object of class \code{"array"}, Shape parameters }
    \item{ \code{method}: }{ Object of class \code{"character"} either \code{"linear"} or \code{"nonlinear"} or \code{"unknown"}. Indicates the regression model used for estimating the tensors.}
}

For class \code{"dwiQball"}:

\describe{
    \item{ \code{order}: }{ Object of class \code{"integer"}, maximal order of Spherical Harmonics to use, needs to be even.} 
    \item{ \code{lambda}: }{ Object of class \code{"numeric"}, nonnegative regularization parameter.} 
    \item{ \code{sphcoef}: }{ Object of class \code{"array"}, estimated coefficients for spherical harmonics, dimension \code{c((order+1)*(order+2)/2,ddim)}.} 
    \item{ \code{sigma}: }{ Object of class \code{"array"}, estimated error variances if \code{method=="linear"}, zero otherwise.} 
    \item{ \code{scorr}: }{ Object of class \code{"numeric"}, estimated spatial correlations in coordinate directions} 
    \item{ \code{bw}: }{ Object of class \code{"numeric"}, bandwidth for a Gaussian kernel that approximately creates the estimated spatial correlations. Needed for adjustments of critical values in the adaptive smoothing algorithm used in function \code{dti.smooth}} 
    \item{ \code{mask}: }{ Object of class \code{"array"}, logical indicating the voxel where the tensor was estimated.} 
    \item{ \code{hmax}: }{ Object of class \code{"numeric"}, maximal bandwidth in case of adaptive smoothing, 1 otherwise.} 
    \item{ \code{outlier}: }{ Object of class \code{"numeric"}, index of voxel where physical constraints are not met, i.e. where the observed values in gradient images Si were larger than the corresponding S0 values. These are probably motion effects or registration errors. Values are replaced by the corresponding (mean) S0 values.}
    \item{\code{scale}:}{Numerical value corresponding to the 95\% quantile of the maximal eigenvalues of estimated tensors within the mask.
          Used for scaling in function \code{show3d.dwiQball}}
    \item{ \code{what}: }{ Object of class \code{"character"}, \code{"ODF"}, \code{"wODF"}, \code{"aODF"} or \code{"ADC"}. Indicates if the object contains coefficients of the orientation density function (ODF (Descoteaux 2007), wODF (Sapiro(2009) or aODF) or the apparent diffusion coefficient (ADC). Coefficients are computed with respect to spherical harmonics of the specified order.}
}
For class \code{"dwiFiber"}:

\describe{
    \item{ \code{fibers}: }{ Object of class \code{"matrix"}, Matrix of fibers. The first three columns contain the coordinates of the track points, the last three columns the direction vectors for each of these points.} 
    \item{ \code{startind}: }{Object of class \code{"integer"}, indices for the first dimension of \code{fibers} where
                              coordinates for a new fiber start.}
    \item{ \code{roix}: }{ Object of class \code{"integer"}, coordinate range of region of interest in x-direction } 
    \item{ \code{roiy}: }{ Object of class \code{"integer"}, coordinate range of region of interest in x-direction } 
    \item{ \code{roiz}: }{ Object of class \code{"integer"}, coordinate range of region of interest in x-direction } 
    \item{ \code{method}: }{ Object of class \code{"character"}, fiber tracking method.} 
    \item{ \code{minanindex}: }{ Object of class \code{"numeric"}, minimal fractional anisotropy index} 
    \item{ \code{maxangle}: }{ Object of class \code{"numeric"}, maximal angle between fibres.} 
}
For class \code{"dwiMixtensor"}:

\describe{
    \item{ \code{model}: }{ Object of class \code{"character"}, characterizes the type of the 
            mixed tensor model. Currently the only implemented model is \code{model="homogeneous_prolate"}.
           } 
    \item{ \code{ev}: }{ Object of class \code{"array"}, estimated eigenvalues, dimension \code{c(2,ddim)}} 
    \item{ \code{mix}: }{ Object of class \code{"array"}, estimated mixture coefficients, dimension \code{c(nmix,ddim)}. \code{nmix} is the number of mixture components specified.} 
    \item{ \code{orient}: }{ Object of class \code{"array"}, estimated tensor orientations, dimension \code{c(2,nmix,ddim)}} 
    \item{ \code{th0}: }{ Object of class \code{"array"}, estimated intensities in S0 images, dimension \code{ddim}} 
    \item{ \code{sigma}: }{ Object of class \code{"array"}, estimated error variances if \code{method=="linear"}, zero otherwise.} 
    \item{ \code{scorr}: }{ Object of class \code{"numeric"}, estimated spatial correlations in coordinate directions} 
    \item{ \code{bw}: }{ Object of class \code{"numeric"}, bandwidth for a Gaussian kernel that approximately creates the estimated spatial correlations. Needed for adjustments of critical values in the adaptive smoothing algorithm used in function \code{dti.smooth}} 
    \item{ \code{mask}: }{ Object of class \code{"array"}, logical indicating the voxel where the tensor was estimated.} 
    \item{ \code{hmax}: }{ Object of class \code{"numeric"}, maximal bandwidth in case of adaptive smoothing, 1 otherwise.} 
    \item{ \code{outlier}: }{ Object of class \code{"numeric"}, index of voxel where physical constraints are not met, i.e. where the observed values in gradient images Si were larger than the corresponding S0 values. These are probably motion effects or registration errors. Values are replaced by the corresponding (mean) S0 values.}
    \item{\code{scale}:}{Numerical value corresponding to the 95\% quantile of the maximal eigenvalues of estimated tensors within the mask.
          Used for scaling in function \code{show3d.dtiTensor}}
    \item{ \code{method}: }{ Object of class \code{"character"}, either \code{"mixtensor"} or \code{"Jian"}. Indicates the regression model used for estimating the tensors.}
}

}
\section{Methods}{
  Methods only operate on subclasses \code{"dtiData"}, \code{"dtiTensor"},  \code{"dtiIndices"}, \code{"dwiQball"} and \code{"dwiFiber"}.
  \describe{
    \item{dti.smooth}{ Create estimates of diffusion tensors in each voxel using structural adaptive spatial smoothing. }
    \item{dtiTensor}{\code{signature(object = "dtiData")}: Create estimates of diffusion tensors in each voxel. }
    \item{dtiIndices}{\code{signature(object = "dtiTensor")}: Create estimates of diffusion tensors indices in each voxel. }
    \item{tracking}{\code{signature(object = "dtiTensor")} or \code{signature(object = "dtiIndices")}: Fiber tracking. }
    \item{dtiQball}{\code{signature(object = "dtiData")}: Create estimates of ADC-parameters with respect to a sherical harmonics ortho-normal system. }
    \item{show3d}{ Method for Function `show3d' in Package `dti'. }
    \item{plot}{ Method for Function `plot' in Package `dti'. }
    \item{print}{ Method for Function `print' in Package `dti'. }
    \item{summary}{ Method for Function `summary' in Package `dti'. }
  }
}
\author{
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{
  \code{\link{dtiData}},
  \code{\link{readDWIdata}},
  \code{\link{sdpar-methods}},
  \code{\link{dtiTensor-methods}}, 
  \code{\link{dwiMixtensor-methods}}, 
  \code{\link{dti.smooth-methods}},
  \code{\link{dtiIndices-methods}},
  \code{\link{dwiQball-methods}}, 
  \code{\link{tracking-methods}},
  \code{\link{show3d-methods}},
  \code{\link{plot-methods}},
  \code{\link{print-methods}},
  \code{\link{summary-methods}},
  \code{\link{extract-methods}}
}
\keyword{ classes }
