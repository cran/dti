show3d.tens <- function(radii,polyeder,centers=NULL,colors=NULL,alpha=1,...){
   if(is.null(centers)){
      centers <- matrix(0,3,1)
      n <- 1
   } else {
      dcenters <- dim(centers)
      if(length(dcenters)!=2 || dcenters[1]!=3) stop("centers needs to be NULL or a matrix 
      with dimension (3,n)")
      n <- dcenters[2]
   }
   if(is.null(colors)){
      colors <- heat.colors(1)
   } 
   if(length(colors)!=n){
      nc <- length(colors)
      nnc <- n%/%nc+1
      colors <- rep(colors,nnc)[1:n]
   }
   if(is.null(alpha)){
      alpha <- 1
   } 
   if(length(alpha)!=n){
      nc <- length(alpha)
      nnc <- n%/%nc+1
      alpha <- rep(alpha,nnc)[1:n]
   }
   nv <- polyeder$nv
   ni <- polyeder$ni*3
   colors <- t(matrix(colors,n,ni))
   alpha <- t(matrix(alpha,n,ni))
   vertices <- array(polyeder$vertices,c(3,nv,n))
   indices <- matrix(polyeder$indices,c(ni,n))
   if(length(radii)!=nv*n) stop("wrong length of radii, needs to be 
             dim(polyeder$vertices)[2]*dim(centers)[2]")
   vertices <- sweep(vertices,2:3,radii,"*")
   vertices <- sweep(vertices,c(1,3),centers,"+")
   dim(vertices) <- c(3,nv*n)
   indices <- sweep(matrix(indices,ni,n),2,((1:n)-1)*nv,"+")
   rgl.triangles(vertices[1,indices],vertices[2,indices],vertices[3,indices],
                 color=colors,alpha=alpha,...)
}
show3d.data <- function(radii,vertices,centers=NULL,minalpha=1,...){
   if(is.null(centers)){
      centers <- matrix(0,3,1)
      n <- 1
   } else {
      dcenters <- dim(centers)
      if(length(dcenters)!=2 || dcenters[1]!=3) stop("centers needs to be NULL or a matrix 
      with dimension (3,n)")
      n <- dcenters[2]
   }
   maxradii <- apply(radii,2,max)
   alpha <- minalpha+(1-minalpha)*sweep(radii,2,maxradii,"/")
   colors <- rgb(abs(vertices[1,]),abs(vertices[2,]),abs(vertices[3,]))
   if(length(alpha)!=n){
      nc <- length(alpha)
      nnc <- n%/%nc+1
      alpha <- rep(alpha,nnc)[1:n]
   }
   nv <- dim(vertices)[2]
   indices <- gettriangles(vertices)$triangles
   ni <- length(indices)
   indices <- sweep(matrix(indices,ni,n),2,((1:n)-1)*nv,"+")
   colors <- matrix(colors,nv,n)
   alpha <- matrix(alpha,nv,n)
   vertices <- array(vertices,c(3,nv,n))
   if(length(radii)!=nv*n) stop("wrong length of radii, needs to be 
             dim(vertices)[2]*dim(centers)[2]")
   vertices <- sweep(vertices,2:3,radii,"*")
   vertices <- sweep(vertices,c(1,3),centers,"+")
   dim(vertices) <- c(3,nv*n)
   rgl.triangles(vertices[1,indices],vertices[2,indices],vertices[3,indices],
                 color=colors[indices],alpha=alpha[indices],...)
#   rgl.points(vertices[1,],vertices[2,],vertices[3,],
#                 color=colors,alpha=alpha,size=2,...)
}
show3d.cdata <- function(radii,polyeder,centers=NULL,minalpha=1,scale=.5,...){
   if(is.null(centers)){
      centers <- matrix(0,3,1)
      n <- 1
   } else {
      dcenters <- dim(centers)
      if(length(dcenters)!=2 || dcenters[1]!=3) stop("centers needs to be NULL or a matrix 
      with dimension (3,n)")
      n <- dcenters[2]
   }
   vertices <- polyeder$vertices
   colors <- rainbow(1000,start=0,end=.7)[pmax(1,pmin(1000,as.integer(1000*radii)))]
   alpha <- minalpha+(1-minalpha)*radii
   nv <- polyeder$nv
   ni <- polyeder$ni*3
   colors <- matrix(colors,nv,n)
   alpha <- matrix(alpha,nv,n)
   vertices <- array(polyeder$vertices,c(3,nv,n))
   indices <- matrix(polyeder$indices,c(ni,n))
   if(length(radii)!=nv*n) stop("wrong length of radii, needs to be 
             dim(polyeder$vertices)[2]*dim(centers)[2]")
   vertices <- sweep(scale*vertices,c(1,3),centers,"+")
   dim(vertices) <- c(3,nv*n)
   indices <- sweep(matrix(indices,ni,n),2,((1:n)-1)*nv,"+")
   rgl.triangles(vertices[1,indices],vertices[2,indices],vertices[3,indices],
                 color=colors[indices],alpha=alpha[indices],...)
}
show3d.odf <- function(radii,polyeder,centers=NULL,minalpha=1,...){
   if(is.null(centers)){
      centers <- matrix(0,3,1)
      n <- 1
   } else {
      dcenters <- dim(centers)
      if(length(dcenters)!=2 || dcenters[1]!=3) stop("centers needs to be NULL or a matrix 
      with dimension (3,n)")
      n <- dcenters[2]
   }
   vertices <- polyeder$vertices
   colors <- rgb(abs(vertices[1,]),abs(vertices[2,]),abs(vertices[3,]))
   alpha <- rep(minalpha,length(radii))
   nv <- polyeder$nv
   ni <- polyeder$ni*3
   colors <- matrix(colors,nv,n)
   alpha <- matrix(alpha,nv,n)
   vertices <- array(polyeder$vertices,c(3,nv,n))
   indices <- matrix(polyeder$indices,c(ni,n))
   if(length(radii)!=nv*n) stop("wrong length of radii, needs to be 
             dim(polyeder$vertices)[2]*dim(centers)[2]")
   vertices <- sweep(vertices,2:3,radii,"*")
   vertices <- sweep(vertices,c(1,3),centers,"+")
   dim(vertices) <- c(3,nv*n)
   indices <- sweep(matrix(indices,ni,n),2,((1:n)-1)*nv,"+")
   rgl.triangles(vertices[1,indices],vertices[2,indices],vertices[3,indices],
                 color=colors[indices],alpha=alpha[indices],...)
}
