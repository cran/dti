      subroutine tensres(th0,D,s,nvox,nb,b,res,rss)
C
C  compute residuals
C
      implicit logical (a-z)
      integer nb,nvox
      real*8 th0(nvox),D(6,nvox),s(nb,nvox),b(6,nb),res(nb,nvox),
     1       rss(nvox)
      integer i,j,k
      real*8 zrss,z
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(nb,nvox,th0,D,s,b,res,rss)
C$OMP& PRIVATE(i,j,k,z,zrss)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,nvox
         zrss = 0.d0
         DO j = 1,nb
            z = 0.d0
            DO k = 1,6
               z = z - b(k,j)*D(k,i)
            END DO
            z = s(j,i) - th0(i)*dexp(z)
            res(j,i) = z
            zrss = zrss + z*z
         END DO
         rss(i) = zrss
      END DO
C$OMP END DO 
C$OMP END PARALLEL
      RETURN
      END
      subroutine ftensor(par,s,nb,b,vinv,gv,fv)
C
C  compute f(par) 
C
      implicit logical (a-z)
      integer nb
      real*8 par(7),s(nb),b(6,nb),vinv(nb),gv(nb),fv
      integer i
      real*8 D(6),th0,rss,res
      th0=par(1)
C      call dblepr("par",3,par,7)
      call rho2D(par(2),D)
      call sihat(th0,D,b,gv,nb)
C      call dblepr("th0",3,th0,1)
C      call dblepr("D",1,D,6)
C      call dblepr("vinv",4,vinv,6)
C      call dblepr("gv",2,gv,nb)
C      call dblepr("si",2,s,nb)
C
C   this gives vector  th0*exp(-b g_i^T D g_i) in gv
C
      rss=0.d0
      DO i=1,nb
         res=s(i)-gv(i)
         rss=rss+res*res*vinv(i)
      END DO
      fv=rss
C      call dblepr("fv",2,fv,1)
C
C   now we have value of the criterion in fv
C
      RETURN
      END
      subroutine gtensor(par,s,nb,b,vinv,gv,fv,grad)
C
C  compute gradient of f(par) 
C
      implicit logical (a-z)
      integer nb
      real*8 par(7),s(nb),b(6,nb),vinv(nb),gv(nb),fv(nb),grad(7)
      integer i
      real*8 D(6),th0,z,z1,res
      th0=par(1)
      call rho2D(par(2),D)
      call sihat(th0,D,b,gv,nb)
C
C   this gives vector  th0*exp(-b g_i^T D g_i) in gv
C
      DO i=1,nb
         res=s(i)-gv(i)
         fv(i)=2.d0*res*vinv(i)
      END DO
C
C     derivative with respect to theta0
C
      z=0.d0
      DO i=1,nb
         z=z+fv(i)*gv(i)
      END DO
      grad(1)=-z/th0
C
C     derivatives with respect to r
C
      DO i=2,7
         grad(i)=0.d0
      END DO
      DO i=1,nb
         z=fv(i)*gv(i)
         z1=2.d0*b(1,i)*par(2)+b(2,i)*par(3)+b(3,i)*par(4)
         grad(2)=grad(2)+z*z1
         z1=2.d0*b(4,i)*par(3)+b(2,i)*par(2)+b(5,i)*par(4)
         grad(3)=grad(3)+z*z1
         z1=2.d0*b(6,i)*par(4)+b(3,i)*par(2)+b(5,i)*par(3)
         grad(4)=grad(4)+z*z1
         z1=2.d0*b(4,i)*par(5)+b(5,i)*par(6)
         grad(5)=grad(5)+z*z1
         z1=2.d0*b(6,i)*par(6)+b(5,i)*par(5)
         grad(6)=grad(6)+z*z1
         z1=2.d0*b(6,i)*par(7)
         grad(7)=grad(7)+z*z1         
      END DO       
C      call dblepr("grad",4,grad,7)
C
C     We now have the gradient in grad
C      
      RETURN
      END
