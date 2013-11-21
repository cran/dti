dkiIndices2 <- function(object, 
                   mc.cores = getOption( "mc.cores", 2L),
                   verbose = FALSE) {
			 
            if ( verbose) cat( "dkiTensor: entering function", format( Sys.time()), "\n")
            
            ## call history  
            args <- c( object@call, sys.call(-1))

            ## we need this for all the arrays
            ddim <- object@ddim
            nvox <- prod( ddim)
            nvox0 <- sum( object@mask)

            ## perform the DTI indices calculations
            D <- object@D
            z <- dtiind3D( object@D[ 1:6, , , ], object@mask, mc.cores = mc.cores)

            if ( verbose) cat( "dkiTensor: DTI indices calculated", format( Sys.time()), "\n")

            ## perform the DKI indices calculations
            ## DO WE NEED THIS?
            if ( mc.cores > 1) {
              mc.cores.old <- setCores( , reprt = FALSE)
              setCores( mc.cores)
            }
            
            dim( D) <- c( 21, nvox)
            andir <- matrix( 0, 9, nvox)
            lambda <- matrix( 0, 3, nvox)
            
            t1 <- Sys.time()

            zz <- .Fortran("dti3DevAll",
                           as.double( D[ 1:6, object@mask]),
                           as.integer( nvox0),
                           andir = double( 9*nvox0),
                           evalues = double( 3*nvox0),
                           DUP = FALSE,
                           PACKAGE = "dti")[ c( "andir", "evalues")]
            andir[ , object@mask] <- zz$andir
            lambda[ , object@mask] <- zz$evalues 

            t2 <- Sys.time()
            if ( verbose) cat( "dkiTensor: calculation took ", difftime( t2, t1), attr(difftime( t2, t1), "units"), " for", nvox0, "voxel\n")
           
            if ( mc.cores > 1) setCores( mc.cores.old, reprt = FALSE)
            dim( andir) <- c( 3, 3, nvox)
            
            ind <- object@bvalue != 0
            xxx <- dkiDesign( object@gradient[ , ind])
            Tabesh_AD <- xxx[ , 1:6]
            Tabesh_AK <- xxx[ , 7:21]
            Dapp <- Tabesh_AD %*% D[ 1:6, ] ## ngrad*prod(ddim)
            Kapp <- z$md^2 * (Tabesh_AK %*% D[ 7:21, ] / Dapp^2)
            mk <- apply( Kapp, 2, mean)
           
            x1 <- dkiDesign( andir[ , 1, ])
            x2 <- dkiDesign( andir[ , 2, ])
            x3 <- dkiDesign( andir[ , 3, ])

            ## cannot allocate memory for the following:
            ## w1111 <- diag( x1[ , 7:21] %*% D[ 7:21, ])
            w1111 <- numeric( nvox)
            for ( i in 1:nvox) w1111[ i] <- x1[ i, 7:21] %*% D[ 7:21, i]
            w2222 <- numeric( nvox)
            for ( i in 1:nvox) w2222[ i] <- x2[ i, 7:21] %*% D[ 7:21, i]
            w3333 <- numeric( nvox)
            for ( i in 1:nvox) w3333[ i] <- x3[ i, 7:21] %*% D[ 7:21, i]
            
            k1 <- z$md^2 / lambda[ 1, ]^2 * w1111
            k2 <- z$md^2 / lambda[ 2, ]^2 * w2222
            k3 <- z$md^2 / lambda[ 3, ]^2 * w3333


            ## old voxelwise test version
            ##            k1 <- k2 <- k3 <- numeric( nvox)

            ##            if (verbose) pb <- txtProgressBar(0, nvox, style = 3)
            ##            for ( i in 1:nvox) {
            ##              if (verbose) setTxtProgressBar(pb, i)
            ##              xxx <- dkiDesign( andir[ , , i])
            
            ##              w1111 <- xxx[ 1, 7:21] %*% D[ 7:21, i]
            ##              w2222 <- xxx[ 2, 7:21] %*% D[ 7:21, i]
            ##              w3333 <- xxx[ 3, 7:21] %*% D[ 7:21, i]

            ##              k1[ i] <- z$md[ i]^2 / lambda[ 1, i]^2 * w1111
            ##              k2[ i] <- z$md[ i]^2 / lambda[ 2, i]^2 * w2222
            ##              k3[ i] <- z$md[ i]^2 / lambda[ 3, i]^2 * w3333
            
            ##            }
            ##            if (verbose) close(pb)
            

            kbar <- ( k1 + k2 + k3) / 3
            kaxial <- k1
            kradial <- ( k2 + k3) / 2
            fak <- sqrt( 3/2 *( (k1-kbar)^2 + (k2-kbar)^2 + (k3-kbar)^2) / ( k1^2 + k2^2 + k3^2))

            ## we finally got k1, k2, k3, kaxial, kradial, mk, fak
            dim( k1) <- ddim
            dim( k2) <- ddim
            dim( k3) <- ddim
            dim( mk) <- ddim
            dim( kaxial) <- ddim
            dim( kradial) <- ddim
            dim( fak) <- ddim
            
            if ( verbose) cat( "dkiTensor: DKI indices calculated", format( Sys.time()), "\n")
                      
            if ( verbose) cat( "dkiTensor: exiting function", format( Sys.time()), "\n")

            invisible( new("dkiIndices",
                           call = args,
                           fa = array( z$fa, ddim),
                           ga = array( z$ga, ddim),
                           md = array( z$md, ddim),
                           andir = array( z$andir, c( 3, ddim)),
                           bary = array( z$bary, c( 3, ddim)),
                           k1 = k1,
                           k2 = k2,
                           k3 = k3,
                           mk = mk,
                           kaxial = kaxial,
                           kradial = kradial,
                           fak = fak,
                           lambda = lambda,
                           evec = andir,
                           gradient = object@gradient,
                           bvalue = object@bvalue,
                           btb   = object@btb,
                           ngrad = object@ngrad,
                           s0ind = object@s0ind,
                           ddim  = ddim,
                           ddim0 = object@ddim0,
                           voxelext = object@voxelext,
                           orientation = object@orientation,
                           rotation = object@rotation,
                           xind  = object@xind,
                           yind  = object@yind,
                           zind  = object@zind,
                           method = object@method,
                           level = object@level,
                           source = object@source)
                      )
}
