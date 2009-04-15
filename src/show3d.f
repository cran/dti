CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  calculate size of vertices for ellipses in show3d.tensor
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine adcradii(vert,nv,tens,ntens,scale,radii)
      implicit logical (a-z)
      integer nv,ntens
      real*8 vert(3,nv),tens(6,ntens),scale,radii(nv,ntens)
      integer i,j
      real*8 x,y,z,s,xx,xy,xz,yy,yz,zz
      DO i=1,nv
         x=vert(1,i)*scale
         y=vert(2,i)*scale
         z=vert(3,i)*scale
         xx=x*x
         xy=2.d0*x*y
         xz=2.d0*x*z
         yy=y*y
         yz=2.d0*y*z
         zz=z*z
         DO j=1,ntens
            s=tens(1,j)*xx+tens(2,j)*xy+tens(3,j)*xz+
     1        tens(4,j)*yy+tens(5,j)*yz+tens(6,j)*zz
            radii(i,j)=s
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  calculate size of vertices for ellipses in show3d.tensor
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ellradii(vert,nv,tens,ntens,scale,radii)
      implicit logical (a-z)
      integer nv,ntens
      real*8 vert(3,nv),tens(6,ntens),scale,radii(nv,ntens)
      integer i,j
      real*8 x,y,z,s,xx,xy,xz,yy,yz,zz,a,b,c,d,e,f,det,c2,b2,e2
      DO j=1,ntens
         a=tens(1,j)
         b=tens(2,j)
         c=tens(3,j)
         d=tens(4,j)
         e=tens(5,j)
         f=tens(6,j)
         b2=b*b
         c2=c*c
         e2=e*e
         det=a*d*f+2.d0*b*c*e-c2*d-b2*f-e2*a
         if(det.le.1d-6) det=1d-6
         tens(1,j)=(d*f-e2)/det
         tens(2,j)=(e*c-b*f)/det
         tens(3,j)=(b*e-d*c)/det
         tens(4,j)=(a*f-c2)/det
         tens(5,j)=(b*c-a*e)/det
         tens(6,j)=(a*d-b2)/det
      END DO
      DO i=1,nv
         x=vert(1,i)/scale
         y=vert(2,i)/scale
         z=vert(3,i)/scale
         xx=x*x
         xy=2.d0*x*y
         xz=2.d0*x*z
         yy=y*y
         yz=2.d0*y*z
         zz=z*z
         DO j=1,ntens
            s=tens(1,j)*xx+tens(2,j)*xy+tens(3,j)*xz+
     1        tens(4,j)*yy+tens(5,j)*yz+tens(6,j)*zz
            radii(i,j)=1.d0/sqrt(s)
         END DO
      END DO
      RETURN
      END
      real*8 function dotprod3(a,b)
      implicit logical (a-z)
      real*8 a(3),b(3)
      dotprod3=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Compute distances of vertices
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine distvert(vert,nvert,ab,distab,ndist)
      implicit logical (a-z)
      integer nvert,ndist,ab(2,ndist)
C  ndist = nvert*(nvert-1)/2
      real*8 vert(3,nvert),distab(ndist)
      integer i,j,k
      real*8 dotprod3,z
      k=1
      DO i=1,nvert-1
         DO j=i+1,nvert
            ab(1,k)=i
            ab(2,k)=j
            z=dotprod3(vert(1,i),vert(1,j))
            distab(k)=acos(max(-1.d0,min(1.d0,z)))
            k=k+1
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Compute triangulation of the sphere form vertex distances
C   ab and distab are assumed to be sorted by distab in R
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine triedges(ab,distab,iab,ndist,abc,ntria)
      implicit logical (a-z)
      integer ndist,ntria,ab(2,ndist),iab(ndist),abc(3,ntria)
      real*8 distab(ndist)
      integer i,it,m1
C
C   Initialization
C         
      DO i=1,ndist
         iab(i)=0
      END DO
      it=0
      m1=1
C  it is index of last triangle
C  m1 is index of first edge that is used once
C  now get first triangle
      ierr=.TRUE.
      DO while(ierr)
         it=it+1
         call findtri(m1,ab,distab,iab,ndist,abc,it,ierr)
         ierr=.FALSE.
         DO i=1,ndist
            if(iab(i).eq.1) THEN
               ierr=.TRUE.
               m1=i
               EXIT
            END IF
         END DO
         if(it.ge.ntria) THEN
            call intpr("Found it triangles",18,it,1)
            RETURN
         END IF
      END DO
      ntria=it
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Add a new triangle
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine findtri(iedge,ab,distab,iab,ndist,abc,it,ierr)
      implicit logical (a-z)
      integer iedge,ndist,it,ab(2,ndist),iab(ndist),abc(3,it)
      real*8 distab(ndist)
      integer a,b,i,j,k,c,ci,cj
      real*8 z,mindist,dab
      logical ierr,checktri
      external checktri
      a=ab(1,iedge)
      b=ab(2,iedge)
      dab=distab(iedge)
      ierr=.FALSE.
      mindist=1d10
      ci=1
      cj=1
      DO i=1,ndist
         if(i.eq.iedge) CYCLE
         if(iab(i).eq.2) CYCLE
         if(ab(1,i).eq.a) THEN
            k=ab(2,i)
            if(it.gt.1.and.checktri(a,b,k,abc,it-1)) CYCLE
            DO j=1,ndist
               if(iab(j).eq.2) CYCLE
               if(ab(1,j).eq.b) THEN
                  if(ab(2,j).eq.k) THEN
                     z=distab(i)*distab(j)*dab
                     if(z.lt.mindist) THEN
                        c=k
                        ci=i
                        cj=j
                        mindist=z
                        ierr=.TRUE.
                     END IF
                  END IF
               END IF    
               if(ab(2,j).eq.b) THEN
                  if(ab(1,j).eq.k) THEN
                     z=distab(i)*distab(j)*dab
                      if(z.lt.mindist) THEN
                        c=k
                        ci=i
                        cj=j
                        mindist=z
                        ierr=.TRUE.
                     END IF
                  END IF
               END IF    
            END DO
         END IF
         if(ab(2,i).eq.a) THEN
            k=ab(1,i)
            if(it.gt.1.and.checktri(a,b,k,abc,it-1)) CYCLE
            DO j=1,ndist
               if(iab(j).eq.2) CYCLE
               if(ab(1,j).eq.b) THEN
                  if(ab(2,j).eq.k) THEN
                     z=distab(i)*distab(j)*dab
                     if(z.lt.mindist) THEN
                        c=k
                        ci=i
                        cj=j
                        mindist=z
                        ierr=.TRUE.
                     END IF
                  END IF
               END IF    
               if(ab(2,j).eq.b) THEN
                  if(ab(1,j).eq.k) THEN
                     z=distab(i)*distab(j)*dab
                     if(z.lt.mindist) THEN
                        c=k
                        ci=i
                        cj=j
                        mindist=z
                        ierr=.TRUE.
                     END IF
                  END IF
               END IF    
            END DO
         END IF
      END DO
      IF(ierr) THEN
         call sortabc(a,b,c,abc(1,it))
         iab(iedge)=iab(iedge)+1
         iab(ci)=iab(ci)+1
         iab(cj)=iab(cj)+1
      END IF
      RETURN
      END
C
C  sort values a,b,c into abc(3)
C
      subroutine sortabc(a,b,c,abc)
C  sort in decreasing order
      integer a,b,c,abc(3)
      IF(a.lt.b) THEN
         IF(a.lt.c) THEN
            abc(1)=a
            IF(b.lt.c) THEN
               abc(2)=b
               abc(3)=c
            ELSE
               abc(2)=c
               abc(3)=b
            END IF
         ELSE 
            abc(1)=c
            abc(2)=a
            abc(3)=b
         END IF
      ELSE
         IF(b.lt.c) THEN
            abc(1)=b
            IF(a.lt.c) THEN
               abc(2)=a
               abc(3)=c
            ELSE
               abc(2)=c
               abc(3)=a
            END IF
         ELSE 
            abc(1)=c
            abc(2)=b
            abc(3)=a
         END IF
      END IF
      RETURN
      END
C
C   check if triangle exists
C
      logical function checktri(a,b,c,abc,it)
      implicit logical (a-z)
      integer a,b,c,it,abc(3,it),s(3)
      integer i
      checktri=.FALSE.
      DO i=1,it
         call sortabc(a,b,c,s)
         IF(s(1).eq.abc(1,i).and.s(2).eq.abc(2,i).and.
     1                           s(3).eq.abc(3,i)) THEN
            checktri=.TRUE.
            CYCLE
         END IF
      END DO
      RETURN 
      END
C
C   compute radii for vertices of a polyeder from radii on gradients
C
      subroutine datinter(gradii,n,grad,ng,vert,nv,nn,dnn,inn,vradii)
      implicit logical (a-z)
      integer n,ng,nv,nn,inn(nn)
      real*8 gradii(ng,n),grad(3,ng),vert(3,nv),dnn(nn),vradii(nv,n)
      integer i,j,k,mininn
      real*8 mindnn,z,mdist,sdnn,sz,dotprod3
      external dotprod3
      DO i=1,nv
         mindnn=1.d0
         mininn=1
         DO j=1,nn
            dnn(j)=dotprod3(vert(1,i),grad(1,j))
            inn(j)=j
            if(dnn(j).lt.mindnn) THEN
               mindnn=dnn(j)
               mininn=j
            END IF
         END DO
         DO j=nn+1,ng
            z=dotprod3(vert(1,i),grad(1,j))
            if(z.gt.mindnn) THEN
               dnn(mininn)=z
               inn(mininn)=j
               mindnn=dnn(1)
               mininn=1
               IF(nn.eq.1) CYCLE
               DO k=2,nn
                  if(dnn(k).lt.mindnn) THEN
                     mindnn=dnn(k)
                     mininn=k
                  END IF
               END DO
            END IF
         END DO
C
C  we now have the nn nearest neighbors in dnn and inn
C  now interpolate
         IF(nn.gt.1) THEN
            mdist=1.d0
            DO j=1,nn-1
               DO k=j+1,nn
                  mdist=min(mdist,dotprod3(grad(1,k),grad(1,j)))
               END DO
            END DO
         ELSE
            mdist=0.d0
         END IF
         sdnn=0.d0
         sz=0.d0
         DO j=1,nn
            z=max(0.d0,dnn(j)-mdist)
            dnn(j)=z
            sdnn=sdnn+z
         END DO
         DO k=1,n
            sz=0.d0
            DO j=1,nn
               sz=sz+gradii(inn(j),k)*dnn(j)              
            END DO
            vradii(i,k)=sz/sdnn
         END DO
      END DO
      RETURN
      END












