      subroutine replvar(x,ngrad,n1,n2,n3,rind,tind,lind,nind,sigma2,
     1                   mean)
C
C  estimate variances from replicated gradient images
C  returns  df*variance estimate which should be distributed 
C  chisq(df) * variance
C
      implicit logical (a-z)
      integer ngrad,n1,n2,n3,nind
      integer x(ngrad,n1,n2,n3),rind(ngrad),tind(nind),lind(nind)
      real*8 sigma2(n1,n2,n3),mean(n1,n2,n3)
      integer i,j,k,l,m,r,s
      real*8 z
      DO i=1,n1
         DO j=1,n2
            DO k=1,n3
               sigma2(i,j,k)=0.d0
            END DO
         END DO
      END DO
      DO l=1,nind
         if(tind(l).le.1) CYCLE
         DO i=1,n1
            DO j=1,n2
               DO k=1,n3
                  mean(i,j,k)=0.d0
               END DO
            END DO
         END DO
         m=lind(l)
         s=0
         DO r=1,ngrad
            if(rind(r).ne.m) CYCLE
            s=s+1
            DO i=1,n1
               DO j=1,n2
                  DO k=1,n3
                     z=x(r,i,j,k)
                     mean(i,j,k)=mean(i,j,k)+z
                     sigma2(i,j,k)=sigma2(i,j,k)+z*z
                  END DO
               END DO
            END DO
         END DO
         DO i=1,n1
            DO j=1,n2
               DO k=1,n3
                  z=mean(i,j,k)
                  sigma2(i,j,k)=sigma2(i,j,k)-z*z/tind(l)
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
