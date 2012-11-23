\name{readDWIdata}
\alias{dtiData}
\alias{readDWIdata}
\title{ Read Diffusion Weighted Data }
\description{
  The functions create a \code{"dtiData"} object from Diffusion Weighted Data from medicial imaging files in a list of directories or from an imagefile, where the diffusion weighted data is given as 2-byte integer.
}
\usage{
dtiData(gradient, imagefile, ddim, bvalue = NULL, xind = NULL, yind = NULL, zind = NULL, 
        level = 0, mins0value = 0, maxvalue = 32000, voxelext = c(1, 1, 1), 
        orientation = c(0, 2, 5), rotation = diag(3))
readDWIdata(gradient, dirlist, format = c("DICOM", "NIFTI", "ANALYZE", "AFNI"), nslice = NULL, order = NULL, bvalue = NULL,
        xind = NULL, yind = NULL, zind = NULL, level = 0, mins0value = 0, 
        maxvalue = 32000, voxelext = NULL, orientation = c(0, 2, 5), 
        rotation = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gradient}{ matrix of diffusion gradients (including zero gradients for S0 images) }
  \item{imagefile}{ name of data image file (binary 2Byte integers)}
  \item{ddim}{ dimension of image cube (3D) }
  \item{dirlist}{ list of directories containing the data files }
  \item{format}{ string specifying the medical imaging format, one of ''DICOM'', ''NIFTI'', ''ANALYZE'', or ''AFNI'' }
  \item{nslice}{ number of slices (usually z-direction) }
  \item{order}{ vector, specifying a different order of the data files, i.e. other than alphabetic order in the directories given by \code{dirlist}. If not given, 1:n is used for n data files (no order change). }
  \item{bvalue}{ vector of b-values (default 0 for S0 and 1 for Si)}
  \item{xind}{ subindex for x-direction }
  \item{yind}{ subindex for y-direction }
  \item{zind}{ subindex for z-direction }
  \item{level}{ determine \code{mins0value} as quantile of positive S0-values }
  \item{mins0value}{ set voxel in S0-images with values less than \code{level} \dQuote{inactive} }
  \item{maxvalue}{ set voxel with values larger than \code{maxvalue} inactive }
  \item{voxelext}{ voxel extensions in coordinate directions }
  \item{orientation}{ orientations of data as coded in AFNI }
  \item{rotation}{ optional rotation matrix for the coordinate system.}
  \item{verbose}{ some progress reports if TRUE }
}
\details{
  The function \code{dtiData} creates an object of class \code{"dtiData"} from an image file, where the diffusion weighted data is given as 2-byte integer.  This image file has to be prepared by the user.  Use \code{\link{writeBin}} to write out first all S0 images and than all Si images. The \code{gradient} should be created according to this order. Run the demo in order to have an example, how to do this!
  
  The function \code{readDWIdata} reads the data files given in the directories in \code{dirlist} in alphabetic order. The order can be changed using the \code{order} argument: If \code{filelist} is the vector of files in alphabetic order, they are read in the order \code{filelist[order]}.  If order is not given \code{order <- 1:n} is used (no change!).  The medical imaging format is given by \code{format} and can be one of ''DICOM'', ''NIFTI'', ''ANALYZE'', or ''AFNI''.  The number of slices of the three dimensional data cube is given by \code{nslice}.  The diffusion gradients are provided as matrix \code{gradient}.

  \code{xind}, \code{yind}, and \code{zind} define a region of interest as indices.  If not given \code{1:dim[i]} is used. \code{level} determine \code{mins0value} as quantile of positive S0-values.  \code{mins0value} sets voxel in S0-images with values less than \code{level} \dQuote{inactive}.  \code{maxvalue} sets voxel with values larger than \code{maxvalue} inactive.

  \code{voxelext} defines the voxel extension, overwrites the values found in the imaging files.  \code{orientation} codes the data orientation in AFNI notation.
}
\value{
  An object of class \code{"dtiData"}.
}
\references{  
  \url{http://afni.nimh.nih.gov/pub/dist/src/README.attributes}
}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{
  \code{\link{dti.smooth}}, 
  \code{\link{dtiTensor-methods}}, 
  \code{\linkS4class{dtiData}} 
}
\examples{
  \dontrun{demo(dti_art)}
}
\keyword{ IO }
