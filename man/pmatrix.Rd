\name{pmatrix}
\alias{pmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parallelize columnwise computations on a matrix using forking 
(generalizes function 'pvec' from package 'parallel')
}
\description{
     parellelizes the execution of a function on matrix elements
     by splitting the matrix by columns and submitting each part to one core. 
     The function must take a matrix argument, with computations depending on
     row entries. It creates a matrix with same number of columns as on input,
     but possibly different number of rows. 
     The function relies on forking and hence is not available on 
     Windows unless 'mc.cores = 1'.
     The function is designed following the lines of function \code{pvec} from
     \pkg{parallel}.
}
\usage{
pmatrix(v, FUN, ..., mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L), mc.cleanup = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{ matrix to operate on.
   
}
  \item{FUN}{
      function depending on a matrix argument, and providing a matrix valued 
      result. 
}
  \item{\dots}{
   see documentation for function \code{pvec}
}
  \item{mc.set.seed}{ See  'mcparallel'
   see documentation for function \code{pvec}
}
  \item{mc.silent}{
}
  \item{mc.cores}{
   see documentation for function \code{pvec}
}
  \item{mc.cleanup}{
   see documentation for function \code{pvec}
}
}
\details{
   see documentation for function \code{pvec}
}
\value{
Returns a matrix containing the results for the individual colums
  of \code{v} as columns.
}
\author{
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\note{
This function generalizes function \code{pvec} from package \pkg{parallel}
}

\seealso{
\code{\link{pvec}}
}
