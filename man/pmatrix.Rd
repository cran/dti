\name{pmatrix}
\alias{pmatrix}
\alias{plmatrix}
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
}
\usage{
pmatrix(x, FUN, ..., mc.cores = getOption("mc.cores", 2L))
plmatrix(x, FUN, ..., mc.cores = getOption("mc.cores", 2L))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ matrix to operate on.
   
}
  \item{FUN}{
      function depending on a matrix argument, and providing a matrix valued 
      result. 
}
  \item{\dots}{
   see documentation for function \code{pvec}
}
  \item{mc.cores}{
   see documentation for function \code{pvec}
}
}
\details{
   this functions are wrappers to functions  \code{parCapply} and \code{parLapply} from \pkg{parallel}.
}
\value{
Returns a matrix containing the results for the individual colums
  of \code{x} as columns.
}
\author{
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\note{
   These functions are wrappers to functions  \code{parCapply} and \code{parLapply} from \pkg{parallel}.
}

