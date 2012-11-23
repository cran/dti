      subroutine afvarest(y,n1,n2,n3,mask,h,vext,sigma)
      implicit logical (a-z)
      integer n1,n2,n3
      real*8 y(n1,n2,n3),sigma(n1,n2,n3),h,vext(2)
      logical mask(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,ih1,ih2,ih3,ni
      real*8 m1,m2,z
      ih1=h
      ih2=h*vext(1)
      ih3=h*vext(2)
      Do i1=1+ih1,n1-ih1
         Do i2=1+ih2,n2-ih2
            Do i3=1+ih3,n3-ih3
               if(mask(i1,i2,i3)) THEN
                  ni=0
                  m1=0.d0
                  m2=0.d0
                  DO j1=i1-ih1,i1+ih1
                     DO j2=i2-ih2,i2+ih2
                        DO j3=i3-ih3,i3+ih3
                           if(mask(j1,j2,j3)) THEN
                              z=y(j1,j2,j3)
                              m1=m1+z
                              m2=m2+z*z
                              ni=ni+1
                           ENDIF
                        END DO
                     END DO
                  END DO
                  m1=m1/ni
                  m2=m2/ni
                  z=m2-m1*m1
                  if(ni.gt.1) THEN
                     sigma(i1,i2,i3)=ni*z/(ni-1)
                  ELSE
                     sigma(i1,i2,i3)=0.d0
                  ENDIF
               ELSE
                  sigma(i1,i2,i3)=0.d0
               ENDIF
            END DO
         END DO
      END DO
      RETURN
      END
