      subroutine initdata(si,n1,n2,n3,nb,maxvalue)
C
C   project all values to (1,maxvalue) to avoid infinite estimates
C
      integer n1,n2,n3,nb
      double precision   si(n1,n2,n3,nb),sii,maxvalue
      integer i1,i2,i3,k
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               DO k=1,nb
                  sii=si(i1,i2,i3,k)
                  if(sii.le.0.0) si(i1,i2,i3,k)=1
                  if(sii.gt.maxvalue) si(i1,i2,i3,k)=maxvalue
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine outlier(si,n,nb,s0ind,siind,ls0,sinew,ind)
C
C   replace physically meaningless Si values by mean S0
C
      implicit logical(a-z)
      integer n,nb,ls0,s0ind(ls0),siind(*)
      double precision si(nb,n),sinew(nb,n)
      logical ind(n)
      integer i,j1,j,ls0m1
      double precision s0,sji
      logical changed
      ls0m1=ls0-1
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(s0ind,siind,si,sinew,n,nb,ls0,ind)
C$OMP& FIRSTPRIVATE(ls0m1)
C$OMP& PRIVATE(i,j,changed,s0,sji)
C$OMP DO SCHEDULE(STATIC)
      DO i=1,n
         s0=0
         DO j1=1,ls0
            j=s0ind(j1)
            sji=si(j,i)
            s0=s0+sji
            sinew(j,i)=sji
         END DO
         s0=(s0+ls0m1)/ls0
         changed=.FALSE.
         DO j1=1,nb-ls0
            j=siind(j1)
            sji=si(j,i)
            if(sji.ge.s0) THEN
               sji=s0
               changed=.TRUE.
            END IF
            sinew(j,i)=sji
         END DO
         ind(i)=changed
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
      RETURN
      END
      subroutine outlierp(si,n,nb,s0ind,ls0,siind,lsi,sinew,nb1)
C
C   replace physically meaningless Si values by mean S0
C
      implicit logical(a-z)
      integer n,nb,nb1,ls0,lsi,s0ind(ls0),siind(lsi)
      double precision si(nb,n),sinew(nb1,n)
      integer i,j1,j,ls0m1,changed
      double precision s0,sinn(251),sji
      ls0m1=ls0-1
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(s0ind,siind,si,sinew,n,nb,ls0,nb1,lsi)
C$OMP& FIRSTPRIVATE(ls0m1)
C$OMP& PRIVATE(i,j,changed,s0,sji,sinn)
C$OMP DO SCHEDULE(STATIC)
      DO i=1,n
         s0=0
         DO j1=1,ls0
            j=s0ind(j1)
            sji=si(j,i)
            s0=s0+sji
            sinn(j)=sji
         END DO
         s0=(s0+ls0m1)/ls0
         changed=0
         DO j1=1,lsi
            j=siind(j1)
            sji=si(j,i)
            if(sji.ge.s0) THEN
               sji=s0
               changed=1
            END IF
            sinn(j)=sji
         END DO
         sinn(nb1)=changed
         DO j=1,nb1
            sinew(j,i)=sinn(j)
         END DO         
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
      RETURN
      END
     
      subroutine mcorrlag(res,mask,n1,n2,n3,nv,sigma,mean,scorr,lag)

      implicit logical(a-z)
      integer n1,n2,n3,nv,lag(3)
      double precision scorr,res(nv,n1,n2,n3),sigma(n1,n2,n3),
     1       mean(n1,n2,n3)
      logical mask(n1,n2,n3)
      double precision vrm,zcorr,z,mi,mj
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
               vrm=sigma(i1,i2,i3)*sigma(j1,j2,j3)
               if(vrm.le.1e-10) CYCLE
               mi=mean(i1,i2,i3)
               mj=mean(j1,j2,j3)
               zcorr=(res(1,i1,i2,i3)-mi)*(res(1,j1,j2,j3)-mj)
               do i4=2,nv
                  zcorr=zcorr+(res(i4,i1,i2,i3)-mi)*
     1                         (res(i4,j1,j2,j3)-mj)
               enddo
               z=z+zcorr/vrm
               k=k+1
            enddo
         enddo
      enddo
      if( k.gt.0 ) then
         scorr=z/k/nv
      ELSE
         scorr=0.d0
      END IF
      return
      end

      subroutine msd(res,mask,n,nv,sigma,mean)
      implicit logical(a-z)
      integer n,nv
      double precision sigma(n),res(nv,n),mean(n)
      logical mask(n)
      integer i,iv
      double precision z,resi,zm,sigi
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(res,n,mask,nv,sigma,mean)
C$OMP& PRIVATE(z,iv,i,resi,zm,sigi)
C$OMP DO SCHEDULE(GUIDED)
      DO i=1,n
         sigi=0.d0
         zm=0.d0
         if(mask(i)) THEN
            z=0.d0
            zm=0.d0
            DO iv=1,nv
               resi=res(iv,i)
               zm=zm+resi
               z=z+resi*resi
            END DO
            zm=zm/nv
            sigi=sqrt(z/nv-zm*zm)
         ENDIF
         mean(i)=zm
         sigma(i)=sigi
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(mean,sigma)
      RETURN
      END
      
      subroutine mcorr(res,mask,n1,n2,n3,nv,sigma,mean,scorr,l1,l2,l3)

      implicit logical(a-z)
      integer n1,n2,n3,nv,l1,l2,l3,lag(3),n
      double precision scorr(l1,l2,l3),res(nv,n1,n2,n3),
     1       sigma(n1,n2,n3),mean(n1,n2,n3)
      logical mask(n1,n2,n3)
      integer i1,i2,i3
      double precision sci
      n=n1*n2*n3
      call msd(res,mask,n,nv,sigma,mean)
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(res,mask,n1,n2,n3,nv,sigma,mean,scorr,l1,l2,l3)
C$OMP& PRIVATE(lag,i1,i2,i3,sci)
C$OMP DO SCHEDULE(GUIDED)
      Do i1=1,l1
         lag(1)=i1-1
         DO i2=1,l2
            lag(2)=i2-1
            DO i3=1,l3
               lag(3)=i3-1
               call mcorrlag(res,mask,n1,n2,n3,nv,sigma,mean,sci,lag)
               scorr(i1,i2,i3)=sci
            END DO
         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(scorr)
      return
      end
      subroutine thcorr(w,n1,n2,n3,scorr,l1,l2,l3)

      implicit logical(a-z)
      integer n1,n2,n3,l1,l2,l3,lag(3)
      double precision scorr(l1,l2,l3),w(n1,n2,n3)
      integer i1,i2,i3
      double precision z,zcorr
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
      double precision scorr,w(n1,n2,n3)
      integer i1,i2,i3,c1,c2,c3,j1,j2,j3,l1,l2,l3
      double precision z
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
      integer n1,n2,n3,i1,i2,i3,ind1(*),ind2(*),ind3(*)
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
      subroutine getmask(s0,n1,n2,n3,ns,level,msize,prop,s0m,mask)
      implicit logical (a-z) 
      integer n1,n2,n3,ns,msize
      double precision  s0(n1,n2,n3,ns)
      double precision s0m(n1,n2,n3),prop,level
      logical mask(n1,n2,n3)
      integer i1,i2,i3,j,j1,j2,j3
      double precision z,anz,anz1
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               z=0.d0
               DO j=1,ns
                  z=z+s0(i1,i2,i3,j)
               END DO
               s0m(i1,i2,i3)=z/ns
            END DO
         END DO
      END DO
C
C   thats mean s0
C
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               anz=0
               anz1=0
               DO j1=max(1,i1-msize),min(n1,i1+msize)
                  DO j2=max(1,i2-msize),min(n2,i2+msize)
                     DO j3=max(1,i3-msize),min(n3,i3+msize)
                        if(s0m(j1,j2,j3).gt.level) anz1=anz1+1
                        anz=anz+1
                     END DO
                  END DO
               END DO
               mask(i1,i2,i3)=.FALSE.
               if(anz1/anz.gt.prop) mask(i1,i2,i3)=.TRUE.            
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Calculate exceedence probabilities in awstestprop
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine exceed(x,n,z,nz,exprob)
      implicit logical (a-z)
      integer n,nz
      double precision x(n),z(nz),exprob(nz)
      integer i,j
      double precision sk,zj
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(n,nz,x,z,exprob)
C$OMP& PRIVATE(i,j,sk,zj)
C$OMP DO SCHEDULE(GUIDED)
      DO j=1,nz
         sk=0.d0
         zj=z(j)
         DO i=1,n
            if(x(i).gt.zj) sk=sk+1
         END DO
         exprob(j)=sk/n
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(exprob)
      Return
      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  remove interior points from a cube, and return matrix of coords
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rminteri(x,lx,n1,n2,n3,xcoord,nc,value)
      implicit logical (a-z)
      integer n1,n2,n3,xcoord(3,*),nc,x(n1,n2,n3),value(*)
      logical lx(n1,n2,n3),sel
      integer i1,i2,i3,k,a1,a2,a3,e1,e2,e3
      k=1
      DO i1=1,n1
         a1=max(1,i1-1)
         e1=min(n1,i1+1)
         DO i2=1,n2
            a2=max(1,i2-1)
            e2=min(n2,i2+1)
            DO i3=1,n3
               a3=max(1,i3-1)
               e3=min(n3,i3+1)
               if(.not.lx(i1,i2,i3)) CYCLE
               sel=.FALSE.
               if(.not.lx(a1,i2,i3)) sel=.TRUE.
               if(.not.lx(e1,i2,i3)) sel=.TRUE.
               if(.not.lx(i1,a2,i3)) sel=.TRUE.
               if(.not.lx(i1,e2,i3)) sel=.TRUE.
               if(.not.lx(i1,i2,a3)) sel=.TRUE.
               if(.not.lx(i1,i2,e3)) sel=.TRUE.
C               if(a1.eq.1.or.a2.eq.1.or.a3.eq.1) sel=.TRUE.
C               if(e1.eq.n1.or.e2.eq.n2.or.e3.eq.n3) sel=.TRUE.
               if(sel) THEN
                  xcoord(1,k)=i1
                  xcoord(2,k)=i2
                  xcoord(3,k)=i3
                  value(k)=x(i1,i2,i3)
                  k=k+1
               END IF
            END DO
         END DO
      END DO
      nc=k-1
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Hypergeometric 1F1 NIST HB 13.2.2, 13.2.39, 13.2(iv)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine hg1f1(a,b,z,n,fz)
      implicit logical (a-z)
      integer n
      double precision a,b,z(n),fz(n)
      integer i
      double precision x,y,d,eps,zi,ezi,ai,gofbai
      double precision gammaf
      external gammaf
      eps=1.d-15
      gofbai=gammaf(b)/gammaf(b-a)
      DO i=1,n
         d = 1.d0
         zi = z(i)
         IF(zi.lt.0) THEN
            ezi=exp(zi/2)
            ai=b-a
            if(zi.lt.-1400) THEN
               fz(i) = exp((-a)*log(-zi))*gofbai+5.6e-3+1.9e-3*b
C   add +5.6e-3+1.9e-3*b to keep the function monotone 
               CYCLE
            END IF
         ELSE
            ezi=1.d0
            ai=a
         ENDIF
         x = ezi
         y = ezi
         DO WHILE (abs(y).gt.abs(x)*eps)
            y = -y*(ai+d-1.d0)/(b+d-1.d0)*zi/d
            x = x+y
            d = d+1.d0
         END DO
         fz(i) = ezi*x
      END DO
      RETURN
      END
      
