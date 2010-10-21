      subroutine initdat0(si,n1,n2,n3,nb,maxvalue)
C
C   set si()==0 for all voxel that have si-values <=0 or > maxvalue
C
      integer n1,n2,n3,nb,maxvalue,si(n1,n2,n3,nb)
      logical mask
      integer i1,i2,i3,k,sii
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               mask=.FALSE.
               DO k=1,nb
                  sii=si(i1,i2,i3,k)
                  if(sii.le.0.or.sii.gt.maxvalue) THEN
                     mask=.TRUE.
                     CYCLE
                  END IF
               END DO
               IF(mask) THEN
                  DO k=1,nb
                     si(i1,i2,i3,k)=0.d0
                  END DO
               END IF
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine initdata(si,n1,n2,n3,nb,maxvalue)
C
C   project all values to (1,maxvalue) to avoid infinite estimates
C
      integer n1,n2,n3,nb,maxvalue,si(n1,n2,n3,nb)
      integer i1,i2,i3,k,sii
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               DO k=1,nb
                  sii=si(i1,i2,i3,k)
                  if(sii.le.0) si(i1,i2,i3,k)=1
                  if(sii.gt.maxvalue) si(i1,i2,i3,k)=maxvalue
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine outlier(si,n,nb,s0ind,ls0,sinew,ind,lind)
C
C   replace physically meaningless Si values by mean S0
C
      implicit logical(a-z)
      integer n,nb,ls0,sinew(n,nb),ind(n),lind
      real*8 si(n,nb)
      logical s0ind(nb)
      integer i,j,ls0m1
      real*8 s0
      logical changed
      ls0m1=ls0-1
      lind=0
      DO i=1,n
         s0=0
         DO j=1,nb
            if(s0ind(j)) THEN
               s0=s0+si(i,j)
               sinew(i,j)=si(i,j)
            END IF
         END DO
         s0=(s0+ls0m1)/ls0
         changed=.FALSE.
         DO j=1,nb
            if(.not.s0ind(j)) THEN 
               if(si(i,j).gt.s0) THEN
                  sinew(i,j)=s0
                  changed=.TRUE.
               ELSE 
                  sinew(i,j)=si(i,j)
               END IF
            END IF
         END DO
         if(changed) THEN
            lind=lind+1
            ind(lind)=i
         END IF
      END DO
      call rchkusr()
      RETURN
      END
      subroutine mcorrlag(res,mask,n1,n2,n3,nv,sigma,mean,scorr,lag)

      implicit logical(a-z)
      integer n1,n2,n3,nv,lag(3)
      real*8 scorr,res(nv,n1,n2,n3),sigma(n1,n2,n3),mean(n1,n2,n3)
      logical mask(n1,n2,n3)
      real*8 vrm,zcorr,z
      integer i1,i2,i3,i4,l1,l2,l3,k,j1,j2,j3
      l1=lag(1)
      l2=lag(2)
      l3=lag(3)
      z=0.d0
      k=0
C  correlation in x
      do i1=1,n1-l1
         j1=i1+l1
         do i2=1,n2-l2
            j2=i2+l2
            do i3=1,n3-l3
               j3=i3+l3
               if (.not.(mask(i1,i2,i3).and.mask(j1,j2,j3))) CYCLE
               zcorr=0.d0
               do i4=1,nv
                  zcorr=zcorr+(res(i4,i1,i2,i3)-mean(i1,i2,i3))*
     1                        (res(i4,j1,j2,j3)-mean(j1,j2,j3))
               enddo
               vrm=sigma(i1,i2,i3)*sigma(j1,j2,j3)
               if(vrm.gt.1e-10) THEN
                  z=z+zcorr/vrm
                  k=k+1
               end if
            enddo
         enddo
      enddo
      if( k.gt.0 ) then
         scorr=z/k/nv
      ELSE
         scorr=0
      END IF
      return
      end

      subroutine msd(res,mask,n1,n2,n3,nv,sigma,mean)
      implicit logical(a-z)
      integer n1,n2,n3,nv
      real*8 sigma(n1,n2,n3),res(nv,n1,n2,n3),mean(n1,n2,n3)
      logical mask(n1,n2,n3)
      integer i1,i2,i3,iv
      real*8 z,resi,zm
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               if(mask(i1,i2,i3)) THEN
                  z=0.d0
                  zm=0.d0
                  DO iv=1,nv
                     resi=res(iv,i1,i2,i3)
                     zm=zm+resi
                     z=z+resi*resi
                  END DO
                  zm=zm/nv
                  z=z/nv-zm*zm
                  sigma(i1,i2,i3)=sqrt(z)
                  mean(i1,i2,i3)=zm
               ELSE
                  sigma(i1,i2,i3)=0.d0
               ENDIF
            END DO
         END DO
      END DO
      RETURN
      END
      
      subroutine mcorr(res,mask,n1,n2,n3,nv,sigma,mean,scorr,l1,l2,l3)

      implicit logical(a-z)
      integer n1,n2,n3,nv,l1,l2,l3,lag(3)
      real*8 scorr(l1,l2,l3),res(nv,n1,n2,n3),sigma(n1,n2,n3),
     1       mean(n1,n2,n3)
      logical mask(n1,n2,n3)
      integer i1,i2,i3
      call msd(res,mask,n1,n2,n3,nv,sigma,mean)
      Do i1=1,l1
         lag(1)=i1-1
         DO i2=1,l2
            lag(2)=i2-1
            DO i3=1,l3
               lag(3)=i3-1
               call mcorrlag(res,mask,n1,n2,n3,nv,sigma,mean,
     1                       scorr(i1,i2,i3),lag)
               call rchkusr()
            END DO
         END DO
      END DO
      return
      end
      subroutine thcorr(w,n1,n2,n3,scorr,l1,l2,l3)

      implicit logical(a-z)
      integer n1,n2,n3,l1,l2,l3,lag(3)
      real*8 scorr(l1,l2,l3),w(n1,n2,n3)
      integer i1,i2,i3
      real*8 z,zcorr
      z=0.d0
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               z=z+w(i1,i2,i3)*w(i1,i2,i3)
            END DO
         END DO
      END DO
      Do i1=1,l1
         lag(1)=i1-1
         DO i2=1,l2
            lag(2)=i2-1
            DO i3=1,l3
               lag(3)=i3-1
               call thcorlag(w,n1,n2,n3,zcorr,lag)
               scorr(i1,i2,i3)=zcorr/z
               call rchkusr()
            END DO
         END DO
      END DO
      return
      end

      subroutine thcorlag(w,n1,n2,n3,scorr,lag)

      implicit logical(a-z)
      integer n1,n2,n3,lag(3)
      real*8 scorr,w(n1,n2,n3)
      integer i1,i2,i3,c1,c2,c3,j1,j2,j3,l1,l2,l3
      real*8 z
      c1=(n1-1)/2
      c2=(n2-1)/2
      c3=(n3-1)/2
      z=0.d0
      Do i1=-c1,c1
         j1=i1+c1+1
         l1=lag(1)-i1+c1+1
         if(l1.lt.1.or.l1.gt.n1) CYCLE
         DO i2=-c2,c2
            j2=i2+c2+1
            l2=lag(2)-i2+c2+1
            if(l2.lt.1.or.l2.gt.n2) CYCLE
            DO i3=-c3,c3
               j3=i3+c3+1
               l3=lag(3)-i3+c3+1
               if(l3.lt.1.or.l3.gt.n3) CYCLE
               z=z+w(j1,j2,j3)*w(l1,l2,l3)
            END DO
         END DO
      END DO
      scorr=z
      return
      end
      subroutine lconnect(segm,n1,n2,n3,i1,i2,i3,ind1,ind2,ind3,
     1                    mask)
C
C   assumes that we search for a connected region in segm==.TRUE.
C   that contains seed voxel (i1,i2,i3)
C   result: mask == .TRUE. if voxel is connected to seed
      implicit logical (a-z)
      integer n1,n2,n3,i1,i2,i3,ind1(1),ind2(1),ind3(1)
      logical final,mask(n1,n2,n3),segm(n1,n2,n3)
      integer j1,j2,j3,k,l1,l2,l3,lind,lind0,ichecked
C     first find pixel close to (i1,i2) with segm(j1,j2)=0
      DO j1=1,n1
         DO j2=1,n2
            DO j3=1,n3
               mask(j1,j2,j3)=.FALSE.
            END DO
         END DO
      END DO
      if(.not.segm(i1,i2,i3)) THEN
         final=.FALSE.
         DO k=1,max(n1,n2,n3)
            DO l1=-k,k
               DO l2=-k,k
                  DO l3=-k,k
                     if(max(abs(l1),abs(l2),abs(l3)).ne.k) CYCLE
                     j1=i1+l1
                     if(j1.lt.1.or.j1.gt.n1) CYCLE
                     j2=i2+l2
                     if(j2.lt.1.or.j2.gt.n2) CYCLE
                     j3=i3+l3
                     if(j3.lt.1.or.j3.gt.n3) CYCLE
                     if(segm(j1,j2,j3)) THEN
                        final=.TRUE.
                        i1=j1
                        i2=j2
                        i3=j3
                     END IF
                     if(final) EXIT
                  END DO
                  if(final) EXIT
               END DO
               if(final) EXIT
            END DO
            if(final) EXIT
         END DO
      END IF
      mask(i1,i2,i3)=.TRUE.
      ind1(1)=i1
      ind2(1)=i2
      ind3(1)=i3
      lind=1
      lind0=1
      ichecked=1
      final=.FALSE.
      DO while(.not.final)
         DO k=ichecked,lind0
            DO l1=-1,1
               DO l2=-1,1
                  DO l3=-1,1
                     if(l1.eq.0.and.l2.eq.0.and.i3.eq.0) CYCLE
                     j1=ind1(k)+l1
                     if(j1.lt.1.or.j1.gt.n1) CYCLE
                     j2=ind2(k)+l2
                     if(j2.lt.1.or.j2.gt.n2) CYCLE
                     j3=ind3(k)+l3
                     if(j3.lt.1.or.j3.gt.n3) CYCLE
                     if(segm(j1,j2,j3).and..not.mask(j1,j2,j3)) THEN
                        mask(j1,j2,j3)=.TRUE.
                        lind=lind+1
                        ind1(lind)=j1
                        ind2(lind)=j2
                        ind3(lind)=j3
                     END IF
                  END DO
               END DO
            END DO 
         END DO
         if(lind.eq.lind0) THEN
            final=.TRUE.
         ELSE
            ichecked=lind0
            lind0=lind
         END IF
      END DO
      RETURN
      END


