      subroutine roifiber(fiber,newfiber,sizef,ifiber,mlf,startf,
     1                    lengthf,nfibers,roi,nx,ny,nz,vext,sizenf)
C   fiber(sizef,6) - fibers in real coordinates  (in)
C   newfiber(sizef,6) - fibers in real coordinates  (out)
C   ifiber(max(lengthf(nfibers))/2,3) - fiber in voxelcoords
C   startf(nfibers)  -  startpoints of fibers
C   lengthf(nfibers) -  length(fibers)
      implicit logical (a-z)
      integer sizef,nfibers,mlf,ifiber(3,mlf),startf(nfibers),
     1        lengthf(nfibers),nx,ny,nz,sizenf
      real*8 fiber(sizef,6),newfiber(sizef,6),vext(3)
      logical roi(nx,ny,nz)
      integer i,j,k,j1,j2,istart,ilen,nstart
      logical finroi,inroi
      nstart=0
      DO i=1,nfibers
         istart=startf(i)
         ilen=lengthf(i)
         DO j=1,ilen
            DO k=1,3
               j1=istart+2*(j-1)
               j2=j1+1
               ifiber(k,j)=(fiber(j1,k)+fiber(j2,k))/vext(k)/2+1
            END DO
         END DO
C  we now have the voxelcoords for this fiber
         inroi=finroi(ifiber,ilen,roi,nx,ny,nz)
         if(inroi) THEN
            DO j=1,2*ilen
               DO k=1,6
                  newfiber(nstart+j,k)=fiber(istart+j-1,k)
               END DO
            END DO
            nstart=nstart+2*ilen
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