      subroutine roifiber(fiber,newfiber,sizef,ifiber,mlf,startf,
     1           lengthf,nfibers,roi,nx,ny,nz,vext,sizenf,nnfiber)
C   fiber(sizef,6) - fibers in real coordinates  (in)
C   newfiber(sizef,6) - fibers in real coordinates  (out)
C   ifiber(max(lengthf(nfibers))/2,3) - fiber in voxelcoords
C   startf(nfibers)  -  startpoints of fibers
C   lengthf(nfibers) -  length(fibers)
      implicit logical (a-z)
      integer sizef,nfibers,mlf,ifiber(3,mlf),startf(nfibers),
     1        lengthf(nfibers),nx,ny,nz,sizenf,nnfiber
      real*8 fiber(sizef,6),newfiber(sizef,6),vext(3)
      logical roi(nx,ny,nz)
      integer i,j,k,j1,istart,ilen,nstart
      logical finroi,inroi
      nstart=0
      nnfiber=0
      DO i=1,nfibers
         istart=startf(i)
         ilen=lengthf(i)
         DO j=1,ilen
            DO k=1,3
               j1=istart+j-1
               ifiber(k,j)=fiber(j1,k)/vext(k)+1
            END DO
         END DO
C  we now have the voxelcoords for this fiber
         inroi=finroi(ifiber,ilen,roi,nx,ny,nz)
         if(inroi) THEN
            DO j=1,ilen
               DO k=1,6
                  newfiber(nstart+j,k)=fiber(istart+j-1,k)
               END DO
            END DO
            nnfiber=nnfiber+1
            startf(nnfiber)=nstart+1
            nstart=nstart+ilen
         END IF
      END DO
      sizenf=nstart
      RETURN
      END
      logical function finroi(ifiber,sf,roi,nx,ny,nz)
      implicit logical (a-z)
      integer sf,ifiber(3,sf),nx,ny,nz
      logical roi(nx,ny,nz)
      integer i
      logical in
      in=.FALSE.
      DO i=1,sf
         IF(roi(ifiber(1,i),ifiber(2,i),ifiber(3,i))) THEN
            in=.TRUE.
            CYCLE
         ENDIF
      END DO
      finroi=in
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Compute start-points for fibers
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine fibersta(fibers,nfibers,start,nstart)
C
C     fibers array of fiber segments (2,nfiber,3)
C     start vector of indices where fibers start
C     nstart  number of fibers
      implicit logical (a-z)
      integer nfibers,nstart,start(nfibers)
      real*8 fibers(2,nfibers,3)
      integer i,ns
      real*8 zd1,zd2,zd3
      start(1)=1
      ns=2
      DO i=1,nfibers-1
         zd1=fibers(2,i,1)-fibers(1,i+1,1)
         zd2=fibers(2,i,2)-fibers(1,i+1,2)
         zd3=fibers(2,i,3)-fibers(1,i+1,3)
         if(zd1*zd1+zd2*zd2+zd3*zd3.gt.1d-12) THEN
            start(ns)=2*i+1
            ns=ns+1
         END IF
      END DO
      nstart=ns-1
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Remove fibers that have a distance less than maxd to longer ones
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine reducefi(fibers,nsegm,startf,endf,nfibers,keep,maxd)
      implicit logical (a-z)
      integer nsegm,nfibers,startf(nfibers),endf(nfibers)
      logical keep(nfibers)
      real*8 fibers(3,nsegm),maxd
      integer i,il,is,ilong,ishort,ncounts,sfil,efil,nlong
      real*8 z1,z2,z3,mindist,f1,f2,f3
      DO i=1,nfibers
         keep(i)=.TRUE.
      END DO
      ncounts=0
      nlong=0
      DO ilong=1,nfibers-1
         if(.not.keep(ilong)) CYCLE
         nlong=nlong+1
         sfil=startf(ilong)
         efil=endf(ilong)
         DO ishort=ilong+1,nfibers
            if(.not.keep(ishort)) CYCLE
            keep(ishort)=.FALSE.
            DO is=startf(ishort),endf(ishort)
               f1=fibers(1,is)
               f2=fibers(2,is)
               f3=fibers(3,is)
               DO il=sfil,efil
                  z1=f1-fibers(1,il)
                  z2=f2-fibers(2,il)
                  z3=f3-fibers(3,il)
                  mindist=z1*z1+z2*z2+z3*z3
                  if(mindist.lt.maxd) EXIT
               END DO
               if(mindist.ge.maxd) THEN
                  keep(ishort)=.TRUE.
                  EXIT
               END IF
            END DO
            if(.not.keep(ishort)) THEN
               ncounts=ncounts+1
            END IF
         END DO
         if((nlong/1000)*1000.eq.nlong) THEN
            call intpr("Inspected Fibers",16,nlong,1)
            call intpr("Current Fiber",13,ilong,1)
            call intpr("removed",7,ncounts,1)
         END IF
         call rchkusr()
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Remove fibers that have a distance less than maxd to longer ones
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine reducefe(fibers,nsegm,startf,endf,nfibers,keep,maxd)
      implicit logical (a-z)
      integer nsegm,nfibers,startf(nfibers),endf(nfibers)
      logical keep(nfibers)
      real*8 fibers(3,nsegm),maxd
      integer i,il,is,ilong,ishort,ncounts,sfil,efil,nlong
      real*8 z1,z2,z3,mindist,f1,f2,f3
      DO i=1,nfibers
         keep(i)=.TRUE.
      END DO
      ncounts=0
      nlong=0
      DO ilong=1,nfibers-1
         if(.not.keep(ilong)) CYCLE
         nlong=nlong+1
         sfil=startf(ilong)
         efil=endf(ilong)
         DO ishort=ilong+1,nfibers
            if(.not.keep(ishort)) CYCLE
            keep(ishort)=.FALSE.
            is=startf(ishort)
            mindist=maxd
            f1=fibers(1,is)
            f2=fibers(2,is)
            f3=fibers(3,is)
            DO il=sfil,efil
               z1=f1-fibers(1,il)
               z2=f2-fibers(2,il)
               z3=f3-fibers(3,il)
               mindist=min(mindist,z1*z1+z2*z2+z3*z3)
            END DO
            if(mindist.ge.maxd) THEN
               keep(ishort)=.TRUE.
               CYCLE
            END IF
            is=endf(ishort)
            mindist=maxd
            f1=fibers(1,is)
            f2=fibers(2,is)
            f3=fibers(3,is)
            DO il=sfil,efil
               z1=f1-fibers(1,il)
               z2=f2-fibers(2,il)
               z3=f3-fibers(3,il)
               mindist=min(mindist,z1*z1+z2*z2+z3*z3)
            END DO
            if(mindist.ge.maxd) THEN
               keep(ishort)=.TRUE.
               CYCLE
            END IF
            ncounts=ncounts+1
         END DO
         if((nlong/1000)*1000.eq.nlong) THEN
            call intpr("Inspected Fibers",16,nlong,1)
            call intpr("Current Fiber",13,ilong,1)
            call intpr("removed",7,ncounts,1)
         END IF
         call rchkusr()
      END DO
      RETURN
      END
