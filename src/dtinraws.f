      integer function indvar(i,j)
      implicit logical (a-z)
      integer i,j
      indvar=i+(j*(j-1))/2
      RETURN 
      END
      subroutine sihat(th0i,Di,btb,F,nb)
      implicit logical (a-z)
      integer nb
      real*8 th0i,Di(6),btb(6,nb),F(nb)
      integer j,k
      real*8 z
      DO k=1,nb
         z=0.d0
         DO j=1,6
            z=z+btb(j,k)*Di(j)
         END DO
         F(k)=th0i*exp(-z)
      END DO
      RETURN
      END
