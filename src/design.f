      subroutine parofor(grad,ngrad,par)
      implicit logical (a-z)
      integer ngrad
      real*8 par(2,ngrad),grad(3,ngrad)
      integer i
      real*8 th,sth,phi,z
      DO i=1,ngrad
         th = acos(grad(3,i))
         sth = sin(th)
         phi = 0.d0
         IF (sth.lt.1d-8) THEN
            th = 0.d0
         ELSE
            z = grad(1,i)/sth
            IF (abs(z).ge.1) THEN
               IF (z.lt.0.d0) THEN
                  phi = 0.d0
               ELSE
                  phi = 3.14159265358979d0
               END IF
            ELSE
               phi = sign(acos(z),grad(2,i))
            END IF
            IF (phi.lt.0) phi = phi+6.2831853071796d0
         END IF
         par(1,i) = th
         par(2,i) = phi
      END DO
      RETURN
      END
      subroutine orofpar(par,ngrad,grad)
      implicit logical (a-z)
      integer ngrad
      real*8 par(2,ngrad),grad(3,ngrad)
      integer i
      real*8 c1,s1,c2,s2
      DO i=1,ngrad
         c1=cos(par(1,i))
         s1=sin(par(1,i))
         c2=cos(par(2,i))
         s2=sin(par(2,i))
         grad(1,i) = s1*c2
         grad(2,i) = s1*s2
         grad(3,i) = c1
      END DO 
      RETURN 
      END
      subroutine optdesi(par,ngrad,grad,value)
      implicit logical (a-z)
      integer ngrad
      real*8 par(2,ngrad),grad(3,ngrad),value
      integer i,j
      real*8 d1,d2,d3,s1,s2,s3,z
      call orofpar(par,ngrad,grad)
      z = 0.d0
      DO i=1,ngrad-1
         DO j=i+1,ngrad
            d1=grad(1,i)-grad(1,j)
            d2=grad(2,i)-grad(2,j)
            d3=grad(3,i)-grad(3,j)
            s1=grad(1,i)+grad(1,j)
            s2=grad(2,i)+grad(2,j)
            s3=grad(3,i)+grad(3,j)
            z=z+1.d0/(d1*d1+d2*d2+d3*d3)+1.d0/(s1*s1+s2*s2+s3*s3)
         END DO
      END DO
      value = z
      RETURN
      END
      subroutine optdesi0(grad,ngrad,value)
      implicit logical (a-z)
      integer ngrad
      real*8 grad(3,ngrad),value
      integer i,j
      real*8 d1,d2,d3,z
      z = 0.d0
      DO i=1,ngrad-1
         DO j=i+1,ngrad
            d1=grad(1,i)-grad(1,j)
            d2=grad(2,i)-grad(2,j)
            d3=grad(3,i)-grad(3,j)
            z=z+1.d0/(d1*d1+d2*d2+d3*d3)
         END DO
      END DO
      value = z
      RETURN
      END
