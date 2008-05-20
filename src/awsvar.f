CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cawsvar(y,mask,n1,n2,n3,hakt,lambda,th,
     1                   lth,bi,thn,spmin,lwght,wght)
C   
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   th       estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   thn       \sum  Wi Y / \sum  Wi    (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   
      implicit logical (a-z)
      external kldist,lkern
      real*8 kldist,lkern
      integer n1,n2,n3
      logical aws,mask(1)
      real*8 y(1),th(1),lth(1),bi(1),thn(1),lambda,wght(2),
     1       hakt,lwght(1),spmin,spf
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3
      real*8 thi,bii,sij,swj,swjy,z1,z2,z3,wj,hakt2,hmax2,lthip1
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      ih1=hakt
      aws=lambda.lt.1d40
C
C     first fill lth with log(th)
C
      DO i1=1,n1*n2*n3
         if(.not.mask(i1)) CYCLE
         lth(i1)=log(th(i1))
      END DO
C
C   first calculate location weights
C
      ih3=hakt/wght(2)
      ih2=hakt/wght(1)
      ih1=hakt
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1+1
      clw2=ih2+1
      clw3=ih3+1
      dlw1=ih1+clw1
      dlw2=ih2+clw2
      dlw3=ih3+clw3
      z2=0.d0
      z3=0.d0
      DO j3=1,dlw3
         if(n3.gt.1) THEN
            z3=(clw3-j3)*wght(2)
            z3=z3*z3
            ih2=sqrt(hakt2-z3)/wght(1)
            jind3=(j3-1)*dlw1*dlw2
	 ELSE
	    jind3=0
	 END IF
         DO j2=clw2-ih2,clw2+ih2
            if(n2.gt.1) THEN
               z2=(clw2-j2)*wght(1)
               z2=z3+z2*z2
               ih1=sqrt(hakt2-z2)
               jind2=jind3+(j2-1)*dlw1
	    ELSE
	       jind2=0
	    END IF
            DO j1=clw1-ih1,clw1+ih1
C  first stochastic term
               jind=j1+jind2
               z1=clw1-j1
               lwght(jind)=1-(z1*z1+z2)/hakt2
               if(lwght(jind).gt.0.d0) hmax2=max(hmax2,z2+z1*z1)
            END DO
         END DO
      END DO
      call rchkusr()
      DO i3=1,n3
         DO i2=1,n2
             DO i1=1,n1
	       iind=i1+(i2-1)*n1+(i3-1)*n1*n2
               if(.not.mask(iind)) CYCLE 
C    nothing to do, final estimate is already fixed by control 
               thi=th(iind)
               lthip1=1.d0+lth(iind)
               bii=bi(iind)/lambda
C   scaling of sij outside the loop
               swj=0.d0
               swjy=0.d0
               DO jw3=1,dlw3
	          j3=jw3-clw3+i3
	          if(j3.lt.1.or.j3.gt.n3) CYCLE
		  jwind3=(jw3-1)*dlw1*dlw2
	          jind3=(j3-1)*n1*n2
                  z3=(clw3-jw3)*wght(2)
                  z3=z3*z3
                  if(n2.gt.1) ih2=sqrt(hakt2-z3)/wght(1)
                  DO jw2=clw2-ih2,clw2+ih2
	             j2=jw2-clw2+i2
	             if(j2.lt.1.or.j2.gt.n2) CYCLE
		     jwind2=jwind3+(jw2-1)*dlw1
	             jind2=(j2-1)*n1+jind3
                     z2=(clw2-jw2)*wght(1)
                     z2=z3+z2*z2
                     ih1=sqrt(hakt2-z2)
                     DO jw1=clw1-ih1,clw1+ih1
C  first stochastic term
	                j1=jw1-clw1+i1
	                if(j1.lt.1.or.j1.gt.n1) CYCLE
                        jind=j1+jind2
                        if(.not.mask(jind)) CYCLE 
                        wj=lwght(jw1+jwind2)
                        z1=(clw1-jw1)
                        z1=z2+z1*z1
                        IF (aws) THEN
                           sij=bii*(thi/th(jind)-lthip1+lth(jind))
C           kldist(thi,th(jind))
                           IF (sij.gt.1.d0) CYCLE
			   IF (sij.gt.spmin) THEN
			      wj=wj*(1.d0-spf*(sij-spmin))
                           END IF
                        END IF
                        swj=swj+wj
                        swjy=swjy+wj*y(jind)
                     END DO
                  END DO
               END DO
               thn(iind)=swjy/swj
               bi(iind)=swj
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
