\name{polyeder}
\docType{data}
\alias{polyeder}
\alias{icosa0}
\alias{icosa1}
\alias{icosa2}
\alias{icosa3}
\alias{icosa4}
\title{Polyeders derived from the Icosahedron (icosa0) by 
       sequential triangulation of surface triangles}
\description{
       icosa0 - icosa4 provide a description of regular polyeders
       derived from the Icosahedron (icosa0) by 
       sequential triangulation of surface triangles
}
\usage{icosax}
\format{a list with components 
    \enumerate{
       \item vertices - array of dimension \code{c(3,nv)}. containing cartesian coordinate of the \code{nv} vertices.
       \item indices - Indices of vertices that define surface triangles of the polyeder.
       \item edges - Indices of vertices that define edges of the polyeder.
       \item nv - number of vertices
       \item ne - number of edges
       \item ni - number of triangles
    }
}
\keyword{datasets}
