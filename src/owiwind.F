C******************************************************************************
C  MODULE OWIWIND
C    Written by s.b. 08/17/2006
C****************************************************************************** 
C 
      MODULE OWIWIND
C
      USE SIZES,ONLY : SZ,MYPROC
      USE GLOBAL,ONLY : NSCREEN, ScreenUnit
#ifdef CMPI
      USE MESSENGER,ONLY : MSG_FINI
#endif         
      IMPLICIT NONE
      character*100 :: title
      character :: part1*12,part2*7,part3*6
      real(SZ), dimension(:,:), allocatable :: uR,vR,pR,uB,vB,pB
      real(SZ), dimension(:), allocatable :: latR,longR,latB,longB
      real(SZ) :: Along, Alat
      real(SZ) :: ramp,rampfrac
      real(SZ) :: uu,vv,PP
      real(SZ) :: Penv
      integer*8 :: date1R,date2R,date1B,date2B
      integer*8 :: date1w,date2w
      integer*8 :: date1p,date2p
      integer :: iLatR,iLongR,iCYMDHR,iMinR
      real(SZ)  :: dxR,dyR,swlatR,swlongR
      integer :: iLatB,iLongB,iCYMDHB,iMinB
      real(SZ)  :: dxB,dyB,swlatB,swlongB
      integer :: iLatw,iLongw,iCYMDHw,iMinw
      real(SZ)  :: dxw,dyw,swlatw,swlongw
      integer :: iLatp,iLongp,iCYMDHp,iMinp
      real(SZ)  :: dxp,dyp,swlatp,swlongp
      integer :: isnapR,updateR
      integer :: isnapB,updateB

      logical :: regionExists

      integer,allocatable :: swpointsR(:,:)
      integer,allocatable :: swpointsB(:,:)

      real(SZ) :: w,w1,w2,w3,w4
      real(SZ),allocatable :: wR(:,:)
      real(SZ),allocatable :: wB(:,:)

      CHARACTER FNAME*60

      integer :: numSets,numBlankSnaps,cntSnaps,numSkipSnaps
      real(SZ) :: windMultiplier


      PUBLIC

C      INTEGER,SAVE,PRIVATE :: REALTYPE, DBLETYPE   
C      REAL(SZ), PRIVATE,ALLOCATABLE :: SENDBUF(:,:), RECVBUF(:,:)
C

C---------------------end of data declarations--------------------------------C


      CONTAINS


C***********************************************************************
C   SOBROUTINE NWS12INIT
C***********************************************************************

      SUBROUTINE NWS12INIT(WVNX,WVNY,PRN,NP,RHOWAT0,G)

      USE SIZES, ONLY : SZ,MYPROC, GBLINPUTDIR

      IMPLICIT NONE
      INTEGER   NP,I
      REAL(SZ)  WVNX(*),WVNY(*),PRN(*)
      REAL(SZ)  RHOWAT0,RHOWATG,G

      allocate(swpointsB(NP,2),wB(NP,4))
      allocate(swpointsR(NP,2),wR(NP,4))
      
! Read meta info -------------------------------------------------

      OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22',STATUS='OLD',
     & ACTION='READ')

      ! Read the number of sets of .pre and .win files.
      ! If numSets = 1 then ADCIRC requires UNIT 221 and 222.
      ! If numSets = 2 then ADCIRC requires UNIT 223 and 224 in addition to 221 and 222.
      !   UNIT 221 and 223 are atmospheric pressure fields
      !   UNIT 222 and 224 are wind velocity fields.
      read(22,*,err=99999) numSets

      if(numSets.NE.1.AND.numSets.NE.2) then
        if (myproc == 0) write(screenunit,1004)
        write(16,1004)
#ifdef CMPI
        call msg_fini()
#endif
        stop
      endif

      ! Read the number of blank snaps to be inserted before OWI winds start
      read(22,*,err=99999) numBlankSnaps

      ! If numBlankSnaps < 0, ABS(numBlankSnaps) snaps in OWI wind files (UNIT221,222,223 and 224) will be skipped.
      if(numBlankSnaps.LT.0) then
        numSkipSnaps = ABS(numBlankSnaps)
        numBlankSnaps = 0
      endif

      ! Read a wind velocity multiplier
      read(22,*,err=99999) windMultiplier

      close(22)

! Read basin pre file header ------------------------------------------------

      OPEN(221,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.221', STATUS='OLD')

      ! Read begining/ending dates of pre file
      read(221,11,err=99999,ADVANCE="NO")part1,part2,part3
 11   FORMAT(a,1x,a,1x,a)
      read(221,*)date1p,date2p
      if (myproc == 0)then
         write(screenunit,*)'date1 in basin pre file = ',date1p
         write(screenunit,*)'date2 in basin pre file = ',date2p
      endif

! Read basin win file header -------------------------------------------------

      OPEN(222,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.222',STATUS='OLD')

      ! Read begining/ending dates of win file
      read(222,'(a,1x,a,1x,a)',err=99999,ADVANCE="NO")part1,part2,part3
      read(222,*)date1w,date2w
      if (myproc == 0)then
         write(screenunit,*)'date1 in basin win file = ',date1w
         write(screenunit,*)'date2 in basin win file = ',date2w
      endif

      if(date1p.ne.date1w.or.date2p.ne.date2w) then
        if (myproc == 0) then
          WRITE(screenunit,1001)
          WRITE(screenunit,1003)
        endif
        WRITE(16,1001)
        WRITE(16,1003)
#ifdef CMPI
        call msg_fini()
#endif
        STOP
      endif

      date1B = date1p
      date2B = date2p

! Check if region scale data exist

      IF(numSets.eq.1) GOTO 100

! Read region pre file header -----------------------------------------------

      OPEN(223,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.223',STATUS='OLD')

      ! Read begining/ending dates of pre file
      read(223,'(a,1x,a,1x,a)',err=99999,ADVANCE="NO")part1,part2,part3
      read(223,*)date1p,date2p
      if (myproc == 0)then
         write(screenunit,*)'date1 in region pre file = ',date1p
         write(screenunit,*)'date2 in region pre file = ',date2p
      endif

! Read region win file header -----------------------------------------------

      OPEN(224,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.224',STATUS='OLD')

      ! Read begining/ending dates of win file
      read(224,'(a,1x,a,1x,a)',err=99999,ADVANCE="NO")part1,part2,part3
      read(224,*)date1w,date2w
      if (myproc == 0)then
         write(screenunit,*)'date1 in region win file = ',date1w
         write(screenunit,*)'date2 in region win file = ',date2w
      endif

      if(date1p.ne.date1w.or.date2p.ne.date2w) then
        if (myproc == 0) then
          WRITE(screenunit,1002)
          WRITE(screenunit,1003)
        endif
        WRITE(16,1002)
        WRITE(16,1003)
#ifdef CMPI
        call msg_fini()
#endif
        STOP
      endif

      date1R = date1p
      date2R = date2p

 100  CONTINUE


! Initialize flags ----------------------------------------------------------
      isnapB = 0
      isnapR = 0
      updateB = 1
      updateR = 1
      cntSnaps = 0

! Skip snaps if necessary ---------------------------------------------------
      do i = 1,numSkipSnaps
        CALL NWS12GET(WVNX,WVNY,PRN,NP,RHOWAT0,G)
      enddo

      RETURN


 10   format(a51,i12,a5,i12)
 1001 FORMAT(1X,'HEADER INFO IN UNIT222 AND UNIT221 DO NOT MATCH')
 1002 FORMAT(1X,'HEADER INFO IN UNIT224 AND UNIT223 DO NOT MATCH')
 1003 FORMAT(1X,'EXECUTION WILL BE TERMINATED')
 1004 FORMAT(//,1X,' NUMBER OF SETS WAS SPECIFIED'//
     &       'INCORRECTLY IN UNIT22.'/
     &       ' IT MUST BE ETEHR OF 1 or 2'/
     &       ' EXECUTION WILL BE TERMINATED.'//)

99999 CONTINUE

#ifdef CMPI
      call msg_fini()
#endif
      STOP 'OWIWIND READ ERROR (1)'

      END SUBROUTINE

C***********************************************************************
C   SOBROUTINE NWS12GET
C***********************************************************************

      SUBROUTINE NWS12GET(WVNX,WVNY,PRN,NP,RHOWAT0,G)
      USE SIZES,ONLY : MYPROC,MNPROC

      IMPLICIT NONE

      INTEGER   NP,I,J,XI,YI
      REAL(SZ)  WVNX(*),WVNY(*),PRN(*)
      REAL(SZ)  RHOWAT0,RHOWATG,G
      CHARACTER*80 PBLJAGF

      RHOWATG=RHOWAT0*G

! Read basin data  ---------------------------------------------------------

      ! Increment counter
      cntSnaps = cntSnaps+1
      
      ! Put a blank snap for the first 'numBlankSnaps' snaps
      if(cntSnaps.LE.numBlankSnaps) then
        do i=1,NP
          WVNX(I)=0.d0
          WVNY(I)=0.d0
          PRN(I)=101300.d0/RHOWATG
        enddo
        IF(NSCREEN.GE.1) THEN
          if(MYPROC.EQ.0) then
            WRITE(screenunit,16) cntSnaps
          endif
          WRITE(16,15) cntSnaps
        ENDIF
        return
      endif
        

      ! Increment counter
      isnapB = isnapB+1

      ! Read grid specifications/date in pressure file
      read (221,11,end=10000,err=99999) 
     &      iLatp,iLongp,dxp,dyp,swlatp,swlongp,iCYMDHp,iMinp

      ! Read grid specifications/date in wind file
      read (222,11,end=10000,err=99999) 
     &      iLatw,iLongw,dxw,dyw,swlatw,swlongw,iCYMDHw,iMinw

      ! Check consistency
      if(iLatp.ne.iLatw.or.iLongp.ne.iLongw.or.dxp.ne.dxw.or. 
     &     dyp.ne.dyw.or.swlatp.ne.swlatw.or.swlongp.ne.swlongw.or. 
     &     iCYMDHp.ne.iCYMDHw.or.iMinp.ne.iMinw) then
        if (myproc == 0) then
          write(screenunit,12)
          write(screenunit,13)
        endif
        write(16,12)
        write(16,13)
#ifdef CMPI
        call msg_fini()
#endif
        stop
      endif

      ! Check if header info has changed from the previous snapshot
      if(isnapB.gt.1) then
        if(iLatp.ne.iLatB.or.iLongp.ne.iLongB.or.dxp.ne.dxB.or.
     &       dyp.ne.dyB.or.swlatp.ne.swlatB.or.
     &       swlongp.ne.swlongB) then
          updateB = 1
        else
          updateB = 0
        endif
      endif
      
      iCYMDHB = iCYMDHp
      iMinB = iMinp
      
      ! Update coordinate mapping coefficients if necessary
      if(updateB.eq.1) then
        call nws12interp_basin(np)
      endif
      
      ! Read pressure
      read(221,22,end=9999,err=99999) ((pB(i,j),i=1,iLongB),j=1,iLatB)

      ! Read u/v components of the wind
      read(222,22,end=9999,err=99999) ((uB(i,j),i=1,iLongB),j=1,iLatB)
      read(222,22,end=9999,err=99999) ((vB(i,j),i=1,iLongB),j=1,iLatB)

! Read region data  --------------------------------------------------------

      regionExists = .FALSE.

      IF(numSets.EQ.1) GOTO 100
      if(iCYMDHB.lt.date1R) goto 100
      if(iCYMDHB.eq.date2R.and.iMinR.ne.0) goto 100
      if(iCYMDHB.gt.date2R) goto 100

      regionExists = .TRUE.

      ! Increment counter         
      isnapR = isnapR+1

      ! Read grid specifications/date in pressure file
      read (223,11,end=10000,err=99999) 
     &     iLatp,iLongp,dxp,dyp,swlatp,swlongp,iCYMDHp,iMinp

      ! Read grid specifications/date in wind file
      read (224,11,end=10000,err=99999) 
     &        iLatw,iLongw,dxw,dyw,swlatw,swlongw,iCYMDHw,iMinw

      if(iLatp.ne.iLatw.or.iLongp.ne.iLongw.or.dxp.ne.dxw.or.
     &     dyp.ne.dyw.or.swlatp.ne.swlatw.or.swlongp.ne.swlongw.or.
     &     iCYMDHp.ne.iCYMDHw.or.iMinp.ne.iMinw) then
        if (myproc == 0) then
          write(screenunit,12)
          write(screenunit,13)
        endif
        write(16,12)
        write(16,13)
#ifdef CMPI
        call msg_fini()
#endif
        stop
      endif

      ! Check if header info has changed from the previous snapshot
      if(isnapR.gt.1) then
        if(iLatp.ne.iLatR.or.iLongp.ne.iLongR.or.dxp.ne.dxR.or. 
     &       dyp.ne.dyR.or.swlatp.ne.swlatR.or.
     &       swlongp.ne.swlongR) then
          updateR = 1
        else
          updateR = 0
        endif

      endif


      iCYMDHR = iCYMDHp
      iMinR = iMinp

      if(iCYMDHB.ne.iCYMDHR.or.iMinB.ne.iMinR) then
        if (myproc == 0) then
          WRITE(screenunit,*) 'SNAPSHOTS NOT SYNCRONIZED'
          WRITE(screenunit,*) '  iCYMDHB=',iCYMDHB, ' iMinB=',iMinB
          WRITE(screenunit,*) '  iCYMDHR=',iCYMDHR, ' iMinR=',iMinR
          WRITE(screenunit,*) 'EXECUTION WILL BE TERMINATED'
        endif
        WRITE(16,*) 'SNAPSHOTS NOT SYNCRONIZED'
        WRITE(16,*) '  iCYMDHB=',iCYMDHB, ' iMinB=',iMinB
        WRITE(16,*) '  iCYMDHR=',iCYMDHR, ' iMinR=',iMinR
        WRITE(16,*) 'EXECUTION WILL BE TERMINATED'
#ifdef CMPI
        call msg_fini()
#endif
        STOP
      endif

      ! Update coordinate mapping coefficients if necessary
      if(updateR.eq.1) then
        call nws12interp_region(np)
      endif

      
      ! Read pressure
      read(223,22,end=9999,err=99999) ((pR(i,j),i=1,iLongR),j=1,iLatR)
      
      ! Read u/v components of the wind
      read(224,22,end=9999,err=99999) ((uR(i,j),i=1,iLongR),j=1,iLatR)
      read(224,22,end=9999,err=99999) ((vR(i,j),i=1,iLongR),j=1,iLatR)

 100  CONTINUE

! Interpolate onto ADCIRC grid and write to file -------------------------
      
      rampfrac = isnapB-1
c      if (rampfrac<36) then
c        ramp = tanh(18d0*rampfrac/36d0)
c      end if
      ramp = 1.0

      IF(NSCREEN.GE.1) THEN
        if(regionExists.EQV..TRUE.) then
          if(MYPROC.EQ.0) then
            WRITE(screenunit,15) iCYMDHB,iMinB
          endif
          WRITE(16,15) iCYMDHB,iMinB
        else
          if(MYPROC.EQ.0) then
            WRITE(screenunit,14) iCYMDHB,iMinB
          endif
          WRITE(16,14) iCYMDHB,iMinB
        endif
      ENDIF
         
      do i=1,NP

        uu=-9999.9D0

        ! BASIN ---------------------------------------------------------

        if (swpointsB(i,1).gt.0) then
          xi = swpointsB(i,1)
          yi = swpointsB(i,2)

          w1=wB(i,1)
          w2=wB(i,2)
          w3=wB(i,3)
          w4=wB(i,4)
          
          uu=w1*uB(xi,yi)+w2*uB(xi+1,yi)+w3*
     &         uB(xi+1,yi+1)+w4*uB(xi,yi+1)
          vv=w1*vB(xi,yi)+w2*vB(xi+1,yi)+w3*
     &         vB(xi+1,yi+1)+w4*vB(xi,yi+1)
          PP=w1*pB(xi,yi)+w2*pB(xi+1,yi)+w3*
     &         pB(xi+1,yi+1)+w4*pB(xi,yi+1)
          
        endif
        
        ! REGION ---------------------------------------------------------
        ! uu, vv and PP will be overwritten if region data exist.

        if ((regionExists).and.(swpointsR(i,1).gt.0)) then
          xi = swpointsR(i,1)
          yi = swpointsR(i,2)
          
          w1=wR(i,1)
          w2=wR(i,2)
          w3=wR(i,3)
          w4=wR(i,4)
          
          uu=w1*uR(xi,yi)+w2*uR(xi+1,yi)+
     &         w3*uR(xi+1,yi+1)+w4*uR(xi,yi+1)
          vv=w1*vR(xi,yi)+w2*vR(xi+1,yi)+
     &         w3*vR(xi+1,yi+1)+w4*vR(xi,yi+1)
          PP=w1*pR(xi,yi)+w2*pR(xi+1,yi)+
     &         w3*pR(xi+1,yi+1)+w4*pR(xi,yi+1)
          
        endif

        ! COPY TO ARRAYS ----------------------------------------------------------

        if(uu.eq.-9999.9D0) then
          WVNX(I)=0.d0
          WVNY(I)=0.d0
          PRN(I)=101300.d0/RHOWATG
        else
          if (rampfrac<36) then
            uu=uu*ramp
            vv=vv*ramp
            PP=Penv-(Penv-PP)*ramp
          endif

          !CONVERT MILLIBARS TO M OF WATER
          PRN(i)  = 100.d0*PP/RHOWATG
 
          ! Apply wind velocity multiplier
          uu = uu * windMultiplier
          vv = vv * windMultiplier

          WVNX(i) = uu
          WVNY(i) = vv

        end if

      enddo

      RETURN
      
9999  continue

      if (myproc == 0) then
      WRITE(screenunit,*)''
      WRITE(screenunit,*)' !!!!!FATAL ERROR!!!!!'
      WRITE(screenunit,*)'EITHER OF UNIT221,222,223,224
     & COULD NOT BE READ'
      WRITE(screenunit,*)' EXECUTION WILL BE TERMINATED'
      WRITE(screenunit,*)''
      endif

      WRITE(16,*) ''
      WRITE(16,*) ' !!!!!FATAL ERROR!!!!!'
      WRITE(16,*) ' EITHER OF UNIT221,222,223,224 COULD NOT BE READ'
      WRITE(16,*) ' EXECUTION WILL BE TERMINATED'
      WRITE(16,*) ''
#ifdef CMPI
      call msg_fini()
#endif
      STOP

10000 continue

      IF(MYPROC.EQ.0) THEN
        WRITE(screenunit,*) ''
        WRITE(screenunit,*) ' !!! WARNING !!!'
        WRITE(screenunit,*) ' EITHER OF UNIT221,222,223,224 RAN OUT'
        WRITE(screenunit,*) ' EXECUTION WILL CONTINUE'
        WRITE(screenunit,*) ''
      ENDIF

      WRITE(16,*) ''
      WRITE(16,*) ' !!! WARNING !!!'
      WRITE(16,*) ' EITHER OF UNIT221,222,223,224 RAN OUT'
      WRITE(16,*) ' EXECUTION WILL CONTINUE'
      WRITE(16,*) ''

      do i=1,NP
        WVNX(I)=0.d0
        WVNY(I)=0.d0
        PRN(I)=101300.d0/RHOWATG
      enddo

      RETURN


 11   format(t6,i4,t16,i4,t23,f6.0,t32,f6.0,
     &     t44,f8.0,t58,f8.0,t69,i10,i2)
 12   format(1X,'SNAPSHOT HEADER IN WIN AND PRE FILES DO NOT MATCH')
 13   format(1X,'EXECUTION WILL BE TERMINATED')
 14   format(/,1X,'PROCESSING BASIN-SCALE WIND DATA',i12,' ',i2)
 15   format(/,1X,'PROCESSING BASIN&REGION-SCALE DATA',i12,' ',i2)
 16   format(/,1X,'INSERTING A BLANK WIND SNAP, COUNT=',i4)
 
 22   format(8f10.0)
      
99999 CONTINUE

#ifdef CMPI
      call msg_fini()
#endif
      STOP 'OWIWIND READ ERROR (2)'

      END SUBROUTINE

C***********************************************************************
C   SOBROUTINE NWS12INTERP_BASIN
C
C   This generates and saves interpolation coefficients for mapping
C   from a basin-scale OWI to a ADCIRC grid.
C
C***********************************************************************

      SUBROUTINE NWS12INTERP_BASIN(NP)
      USE GLOBAL,ONLY : SLAM,SFEA,RAD2DEG

      IMPLICIT NONE

      INTEGER   NP,I,J,K,XI,YI
      REAL(SZ)  adcLat,adcLong

      WRITE(16,*) ''
      WRITE(16,*) 'BASIN-SCALE WIND MAPPING UPDATED'
      WRITE(16,*) ''
        
      iLatB = iLatw
      iLongB = iLongw
      dxB = dxw
      dyB = dyw
      swlatB = swlatw
      swlongB = swlongw
        
      ! Allocate and create matrices
      if(allocated(uB))    deallocate(uB)
      if(allocated(vB))    deallocate(vB)
      if(allocated(pB))    deallocate(pB)
      if(allocated(longB)) deallocate(longB)
      if(allocated(latB))  deallocate(latB)

      allocate(uB(iLongB,iLatB),vB(iLongB,iLatB),pB(iLongB,iLatB))
      allocate(longB(iLongB),latB(iLatB))
        
      ! Generate long&lat on each grid point
      do i=1,iLatB
        latB(i) = swlatB+(i-1)*dyB
      enddo
      do i=1,iLongB
        longB(i) = swlongB+(i-1)*dxB
      enddo
        
      ! Generate interpolation coefficients (south west point and weights)
      do i=1,NP
        adcLat  = RAD2DEG*SFEA(i)
        adcLong = RAD2DEG*SLAM(i)

        if (adcLong>=longB(1).and.adcLong<longB(iLongB) .and. 
     &       adcLat>=latB(1).and.adcLat<latB(iLatB)) then
          do j=1,iLongB-1
            if (adcLong>=longB(j) .and.  
     &           adcLong<longB(j+1)) then
              xi=j
              goto 200
            endif
          enddo
 200      continue
            
          do k=1,iLatB-1
            if (adcLat>=latB(k) .and. 
     &           adcLat<latB(k+1)) then
              yi=k
              goto 300
            endif
          enddo
 300      continue
            
          swpointsB(i,1) = xi
          swpointsB(i,2) = yi
            
          w=(longB(xi+1)-longB(xi))*(latB(yi+1)-latB(yi))
          w1=(longB(xi+1)-adcLong)*(latB(yi+1)-adcLat)
          w2=(adcLong-longB(xi))*(latB(yi+1)-adcLat)
          w3=(adcLong-longB(xi))*(adcLat-latB(yi))
          w4=(longB(xi+1)-adcLong)*(adcLat-latB(yi))
            
          wB(i,1)=w1/w
          wB(i,2)=w2/w
          wB(i,3)=w3/w
          wB(i,4)=w4/w

        else
          swpointsB(i,1) = 0
          swpointsB(i,2) = 0
            
        endif
      enddo
      END SUBROUTINE


C***********************************************************************
C   SOBROUTINE NWS12INTERP_REGION
C
C   This generates and saves interpolation coefficients for mapping
C   from a region-scale OWI to a ADCIRC grid.
C
C***********************************************************************

      SUBROUTINE NWS12INTERP_REGION(NP)
      USE GLOBAL,ONLY : SLAM,SFEA,RAD2DEG
      IMPLICIT NONE

      INTEGER   NP,I,J,K,XI,YI
      REAL(SZ)  adcLat,adcLong

      WRITE(16,*) ''
      WRITE(16,*) 'REGION-SCALE WIND MAPPING UPDATED'
      WRITE(16,*) ''
        
      iLatR = iLatw
      iLongR = iLongw
      dxR = dxw
      dyR = dyw
      swlatR = swlatw
      swlongR = swlongw
 
      ! Allocate and create matrices
      if(allocated(uR))    deallocate(uR)
      if(allocated(vR))    deallocate(vR)
      if(allocated(pR))    deallocate(pR)
      if(allocated(longR)) deallocate(longR)
      if(allocated(latR))  deallocate(latR)

      allocate(uR(iLongR,iLatR),vR(iLongR,iLatR),pR(iLongR,iLatR))
      allocate(longR(iLongR),latR(iLatR))

      ! Generate long&lat on each grid point
      do i=1,iLatR
        latR(i) = swlatR+(i-1)*dyR
      enddo
      do i=1,iLongR
        longR(i) = swlongR+(i-1)*dxR
      enddo

      ! Generate interpolation info (south west point and weights)
      do i=1,NP
        adcLat  = RAD2DEG*SFEA(i)
        adcLong = RAD2DEG*SLAM(i)

        if (adcLong>=longR(1).and.adcLong<longR(iLongR) .and. 
     &       adcLat>=latR(1).and.adcLat<latR(iLatR)) then
          do j=1,iLongR-1
            if (adcLong>=longR(j).and.adcLong<longR(j+1)) then
              xi=j
              goto 201
            endif
          enddo
 201      continue
                  
          do k=1,iLatR-1
            if (adcLat>=latR(k).and.adcLat<latR(k+1)) then
              yi=k
              goto 301
            endif
          enddo
 301      continue

          swpointsR(i,1) = xi
          swpointsR(i,2) = yi

          w=(longR(xi+1)-longR(xi))*(latR(yi+1)-latR(yi))
          w1=(longR(xi+1)-adcLong)*(latR(yi+1)-adcLat)
          w2=(adcLong-longR(xi))*(latR(yi+1)-adcLat)
          w3=(adcLong-longR(xi))*(adcLat-latR(yi))
          w4=(longR(xi+1)-adcLong)*(adcLat-latR(yi))

          wR(i,1)=w1/w
          wR(i,2)=w2/w
          wR(i,3)=w3/w
          wR(i,4)=w4/w
        else
          swpointsR(i,1) = 0
          swpointsR(i,2) = 0
        endif
      enddo
      END SUBROUTINE


      END MODULE OWIWIND
      
      
