      subroutine mfun(par,siq,grad,m,lpar,ngrad,z,erg)
      implicit logical (a-z)
      integer m,lpar,ngrad
      real*8 par(lpar),siq(ngrad),grad(3,ngrad),z(ngrad),erg
      integer i,j
      real*8 dotprod3,c1,c2,w,sw,ew,dir(3),sth,z1,z2
      c1 = exp(par(1))
      c2 = exp(par(2))
      sw = 0
      DO i = 1,m
         IF(i.eq.m) THEN
            w = 1.d0 - sw
         ELSE
            ew = exp(par(3*i+2))
            w = (1.d0 - sw)*ew/(1.d0+ew)
            sw = sw + w
         END IF
         sth = sin(par(i*3))
         dir(1) = sth*cos(par(i*3+1))
         dir(2) = sth*sin(par(i*3+1))
         dir(3) = cos(par(i*3))
         DO j = 1,ngrad
            z1 = dotprod3(dir,grad(1,j))
            z(j) = z(j) + w*exp(-c2 - c1*z1*z1)
         END DO
      END DO
      z1 = siq(1) - z(1)
      z2 = z1*z1
      DO j = 2,ngrad
         z1 = siq(j) - z(j)
         z2 = z2 + z1*z1
      END DO
      erg = z2
      RETURN
      END 
      subroutine mfun2(par,siq,grad,m,p,lpar,ngrad,z,erg)
      implicit logical (a-z)
      integer m,lpar,ngrad
      real*8 par(lpar),siq(ngrad),grad(3,ngrad),z(ngrad),p,erg
      integer i,j
      real*8 dotprod3,c1,c2,c12,w,sw,ew,dir(3),sth,z1,z2
      c1 = exp(par(1))
      c2 = exp(par(2))
      c12 = c1/(c1+c2)
      sw = 0
      DO i = 1,m
         IF(i.eq.m) THEN
            w = 1.d0 - sw
         ELSE
            ew = exp(par(3*i+2))
            w = (1.d0 - sw)*ew/(1.d0+ew)
            sw = sw + w
         END IF
         sth = sin(par(i*3))
         dir(1) = sth*cos(par(i*3+1))
         dir(2) = sth*sin(par(i*3+1))
         dir(3) = cos(par(i*3))
         DO j = 1,ngrad
            z1 = dotprod3(dir,grad(1,j))
            z(j) = z(j) + w*exp(-p*log(1 + ( c2 + c1*z1*z1 )/p))
         END DO
      END DO
      z1 = siq(1) - z(1)
      z2 = z1*z1
      DO j = 2,ngrad
         z1 = siq(j) - z(j)
         z2 = z2 + z1*z1
      END DO
      erg = z2
      RETURN
      END 
      subroutine smsi(si,n,snind,nsn,w,sms)
      implicit logical (a-z)
      integer n,nsn,snind(nsn,n)
      real*8 si(n),sms(n),w(nsn)
      integer i,j
      real*8 z,sw
      sw=0.d0
      DO j=1,nsn
         sw=sw+w(j)
      END DO
      DO i=1,n
         z=0.d0
         DO j=1,nsn
            z=z+w(j)*si(snind(j,i))
         END DO
         sms(i)=z/sw
      END DO
      RETURN
      END