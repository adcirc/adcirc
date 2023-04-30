C******************************************************************************
C  MODULE OWI_RAIN
C    Written by Kendra Dresback, University of Oklahoma 12/27/2013
C    Modified from owiwind.F Written by s.b. 08/17/2006
C
C    The values read in from the fort.425/fort.427 contain rainfall
C    intensity values that has to be given in m/s.
C
C******************************************************************************
C
      MODULE OWI_RAIN
C
      USE SIZES,ONLY : SZ,MYPROC
      USE GLOBAL,ONLY : NSCREEN, ScreenUnit
#ifdef CMPI
      USE MESSENGER,ONLY : MSG_FINI
#endif
      IMPLICIT NONE
      character*100 :: title

!      tcm v48.4645 Changed the way the header line is read
!      character :: part1*12,part2*7,part3*6
!      character :: owiheader*80
      character(80) :: owiheader

      real(SZ), dimension(:,:), allocatable :: CRAINR,CRAINB
      real(SZ), dimension(:), allocatable :: latR,longR,latB,longB
      real(SZ) :: Along, Alat
      real(SZ) :: ramp,rampfrac
      real(SZ) :: RAIN1
      real(SZ) :: CRAIN_env
      integer(8):: date1R,date2R,date1B,date2B
      integer :: date1rain,date2rain
      integer :: iLatR,iLongR,iCYMDHR,iMinR
      real(SZ)  :: dxR,dyR,swlatR,swlongR
      integer :: iLatB,iLongB,iCYMDHB,iMinB
      real(SZ)  :: dxB,dyB,swlatB,swlongB
      integer,save :: iLatc,iLongc,iCYMDHc,iMinc
      real(SZ)  :: dxc,dyc,swlatc,swlongc
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

      PUBLIC

C---------------------end of data declarations--------------------------------C


      CONTAINS


C***********************************************************************
C   SUBROUTINE NCRAIN1_INIT
C***********************************************************************

      SUBROUTINE NCRAIN1_INIT(Crainval,NP)

      USE SIZES, ONLY : SZ,MYPROC, GBLINPUTDIR

      IMPLICIT NONE
      INTEGER   NP,I
      REAL(SZ)  Crainval(*)

      allocate(swpointsB(NP,2),wB(NP,4))
      allocate(swpointsR(NP,2),wR(NP,4))

! Read meta info -------------------------------------------------

      OPEN(27,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.27',STATUS='OLD',
     & ACTION='READ')

      ! Read the number of sets of .rain files.
      ! If numSets = 1 then ADCIRC requires UNIT 425
      ! If numSets = 2 then ADCIRC requires UNIT 425 and 427.
      !   UNIT 425 and 427 are rainfall intensity fields

      read(27,*,err=99999) numSets

      if(numSets.NE.1.AND.numSets.NE.2) then
        if (myproc == 0) write(screenunit,1004)
        write(16,1004)
#ifdef CMPI
        call msg_fini()
#endif
        stop
      endif

      ! Read the number of blank snaps to be inserted before OWI rain fields start
      read(27,*,err=99999) numBlankSnaps

      ! If numBlankSnaps < 0, ABS(numBlankSnaps) snaps in OWI rain files (UNIT225 and 227) will be skipped.
      if(numBlankSnaps.LT.0) then
        numSkipSnaps = ABS(numBlankSnaps)
        numBlankSnaps = 0

      ! v48.4628 TCM 10/28/2009 -- Added else to initialize numSkipSnaps to be 0
      else
        numSkipSnaps = 0
      endif

      close(27)

! Read basin rain file header ------------------------------------------------

      OPEN(425,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.425', STATUS='OLD')

      ! Read begining/ending dates of rainfall intensity file
      owiheader(:) = ' '  !set owiheader to blanks before read
      read(425,fmt='(a80)',err=99999) owiheader
      read(owiheader(56:65),'(I10)') date1rain
      read(owiheader(71:80),'(I10)') date2rain

      if (myproc == 0)then
         write(screenunit,*)'date1 in basin rainfall intensity file = ',date1rain
         write(screenunit,*)'date2 in basin rainfall intensity file = ',date2rain
      endif

      date1B = date1rain
      date2B = date2rain

! Check if region scale data exist

      IF(numSets.eq.1) GOTO 100

! Read region rain file header -----------------------------------------------

      OPEN(427,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.427',STATUS='OLD')

      ! Read begining/ending dates of rain file
      owiheader(:) = ' '  !set owiheader to blanks before read
      read(427,fmt='(a80)',err=99999) owiheader
      read(owiheader(56:65),'(I10)') date1rain
      read(owiheader(71:80),'(I10)') date2rain

      if (myproc == 0)then
         write(screenunit,*)'date1 in region rainfall intensity file = ',date1rain
         write(screenunit,*)'date2 in region rainfall intensity file = ',date2rain
      endif

      date1R = date1rain
      date2R = date2rain

 100  CONTINUE


! Initialize flags ----------------------------------------------------------
      isnapB = 0
      isnapR = 0
      updateB = 1
      updateR = 1
      cntSnaps = 0

! Skip snaps if necessary ---------------------------------------------------
      do i = 1,numSkipSnaps
        CALL NCRAIN_GET(Crainval,NP)
      enddo

      RETURN


 10   format(a51,i12,a5,i12)
 1004 FORMAT(//,1X,' NUMBER OF SETS WAS SPECIFIED'//
     &       'INCORRECTLY IN UNIT23.'/
     &       ' IT MUST BE ETEHR OF 1 or 2'/
     &       ' EXECUTION WILL BE TERMINATED.'//)

99999 CONTINUE

#ifdef CMPI
      call msg_fini()
#endif
      STOP 'OWI_RAIN READ ERROR (1)'

      END SUBROUTINE

C***********************************************************************
C   SUBROUTINE NCRAIN_GET
C***********************************************************************

      SUBROUTINE NCRAIN_GET(Crainval,NP)

      USE SIZES,ONLY : MYPROC,MNPROC
      USE GLOBAL,ONLY : Screenunit

      IMPLICIT NONE

      INTEGER   NP,I,J,XI,YI,K,ibar
      REAL(SZ)  Crainval(*)
      REAL(SZ)  CV(4),CVbar
      CHARACTER*80 PBLJAGF


! Read basin data  ---------------------------------------------------------

      ! Increment counter (cntSnaps initialized to zero in ncrain1_init)
      cntSnaps = cntSnaps+1

      ! Put a blank snap for the first 'numBlankSnaps' snaps
      if(cntSnaps.LE.numBlankSnaps) then
        do i=1,NP
          Crainval(I)=0.d0
        enddo
        IF(NSCREEN.GE.1) THEN
          if(MYPROC.EQ.0) then
            WRITE(screenunit,16) cntSnaps
          endif
          WRITE(16,16) cntSnaps    !TCMv48.4629 (Changed format number from 15 to 16)
        ENDIF
        return
      endif


      ! Increment counter
      isnapB = isnapB+1

      ! Read grid specifications/date in rain file
      read (425,11,end=10000,err=99999)
     &      iLatc,iLongc,dxc,dyc,swlatc,swlongc,iCYMDHc,iMinc

      ! Check if header info has changed from the previous snapshot
      if(isnapB.gt.1) then
        if(iLatc.ne.iLatB.or.iLongc.ne.iLongB.or.dxc.ne.dxB.or.
     &       dyc.ne.dyB.or.swlatc.ne.swlatB.or.
     &       swlongc.ne.swlongB) then
          updateB = 1
        else
          updateB = 0
        endif
      endif

      iCYMDHB = iCYMDHc
      iMinB = iMinc

      IF (myproc == 0 ) THEN
        write(screenunit,'("Processing rain data for time ",I10," ",I2)')
     &              iCYMDHc,iMinc
      END IF

      ! Update coordinate mapping coefficients if necessary
      if(updateB.eq.1) then
        call ncrain1_interp_basin(np)
      endif

      ! Read rain concentrations
      read(425,22,end=9999,err=99999) ((CrainB(i,j),i=1,iLongB),j=1,iLatB)

! Read region data  --------------------------------------------------------

      regionExists = .FALSE.

      IF(numSets.EQ.1) GOTO 100
      if(iCYMDHB.lt.date1R) goto 100
      if(iCYMDHB.eq.date2R.and.iMinR.ne.0) goto 100
      if(iCYMDHB.gt.date2R) goto 100

      regionExists = .TRUE.

      ! Increment counter
      isnapR = isnapR+1

      ! Read grid specifications/date in rain file
      read (427,11,end=10000,err=99999)
     &     iLatc,iLongc,dxc,dyc,swlatc,swlongc,iCYMDHc,iMinc

      ! Check if header info has changed from the previous snapshot
      if(isnapR.gt.1) then
        if(iLatc.ne.iLatR.or.iLongc.ne.iLongR.or.dxc.ne.dxR.or.
     &       dyc.ne.dyR.or.swlatc.ne.swlatR.or.
     &       swlongc.ne.swlongR) then
          updateR = 1
        else
          updateR = 0
        endif

      endif


      iCYMDHR = iCYMDHc
      iMinR = iMinc

      IF (myproc == 0 ) THEN
        write(screenunit,'("Processing 2nd rain set at",I10," ",I2)')
     &              iCYMDHc,iMinc
      END IF

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
        call ncrain1_interp_region(np)
      endif

      ! Read rainfall intensity concentrations
      read(427,22,end=9999,err=99999) ((CrainR(i,j),i=1,iLongR),j=1,iLatR)


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

        rain1=-9999.9D0

        ! BASIN ---------------------------------------------------------

        if (swpointsB(i,1).gt.0) then
          xi = swpointsB(i,1)
          yi = swpointsB(i,2)

          w1=wB(i,1)
          w2=wB(i,2)
          w3=wB(i,3)
          w4=wB(i,4)
!         If an rainfall concentration value is less than zero (land flag)
!         then set the value to the average non-zero value before
!         interpolating
          CV(1) = CrainB(xi,yi)
          CV(2) = CrainB(xi+1,yi)
          CV(3) = CrainB(xi+1,yi+1)
          CV(4) = CrainB(xi,yi+1)
          CVbar = 0.d0
          ibar = 0
          do k=1,4
             if(CV(k) .ge. 0.d0) then
                CVbar = CVbar + CV(k)
                ibar = ibar + 1
             endif
          enddo
          if (ibar.gt.0) CVbar = CVbar/real(ibar,sz)
          do k=1,4
             if (CV(k).lt.0.d0) CV(k) = CVbar
          enddo

!         perform the interpolation
          rain1=w1*CV(1)+w2*CV(2)+w3*CV(3)+w4*CV(4)

        endif

        ! REGION ---------------------------------------------------------
        ! rain1 will be overwritten if region data exist.

        if ((regionExists).and.(swpointsR(i,1).gt.0)) then
          xi = swpointsR(i,1)
          yi = swpointsR(i,2)

          w1=wR(i,1)
          w2=wR(i,2)
          w3=wR(i,3)
          w4=wR(i,4)
!         If an rainfall value is less than zero (land flag)
!         then set the value to the average non-zero value before
!         interpolating
          CV(1) = CrainR(xi,yi)
          CV(2) = CrainR(xi+1,yi)
          CV(3) = CrainR(xi+1,yi+1)
          CV(4) = CrainR(xi,yi+1)
          CVbar = 0.d0
          ibar = 0

          do k=1,4
             if(CV(k) .ge. 0.d0) then
                CVbar = CVbar + CV(k)
                ibar = ibar + 1
             endif
          enddo
          if (ibar.gt.0) CVbar = CVbar/real(ibar,sz)
          do k=1,4
             if (CV(k).lt.0.d0) CV(k) = CVbar
          enddo

!         perform the interpolation
          rain1=w1*CV(1)+w2*CV(2)+w3*CV(3)+w4*CV(4)

        endif

        ! COPY TO ARRAYS ----------------------------------------------------------

        if(rain1.eq.-9999.9D0) then
          Crainval(I)=0.d0
        else
          if (rampfrac<36) then
            rain1=Crain_env-(Crain_env-rain1)*ramp
          endif
          Crainval(i) = rain1
        end if

! Convert the rainfall from mm/hr to m/s
        Crainval(i)=Crainval(i)*0.00000027778d0

      enddo

      RETURN

9999  continue

      if (myproc == 0) then
      WRITE(screenunit,*)''
      WRITE(screenunit,*)' !!!!!FATAL ERROR!!!!!'
      WRITE(screenunit,*)'EITHER UNIT 425 or 427',
     & ' COULD NOT BE READ'
      WRITE(screenunit,*)' EXECUTION WILL BE TERMINATED'
      WRITE(screenunit,*)''
      endif

      WRITE(16,*) ''
      WRITE(16,*) ' !!!!!FATAL ERROR!!!!!'
      WRITE(16,*) ' EITHER UNIT 425 or 427 COULD NOT BE READ'
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
        WRITE(screenunit,*) ' EITHER UNIT 425 or 427 RAN OUT'
        WRITE(screenunit,*) ' EXECUTION WILL CONTINUE'
        WRITE(screenunit,*) ''
      ENDIF

      WRITE(16,*) ''
      WRITE(16,*) ' !!! WARNING !!!'
      WRITE(16,*) ' EITHER UNIT 425 or 427 RAN OUT'
      WRITE(16,*) ' EXECUTION WILL CONTINUE'
      WRITE(16,*) ''

      do i=1,NP
        Crainval(I)=0.d0
      enddo

      RETURN


 11   format(t6,i4,t16,i4,t23,f6.0,t32,f6.0,
     &     t44,f8.0,t58,f8.0,t69,i10,i2)
 12   format(1X,'SNAPSHOT HEADER IN RAINFALL FILES DO NOT MATCH')
 13   format(1X,'EXECUTION WILL BE TERMINATED')
 14   format(/,1X,'PROCESSING BASIN-SCALE RAINFALL DATA',i12,' ',i2)
 15   format(/,1X,'PROCESSING BASIN&REGION-SCALE DATA',i12,' ',i2)
 16   format(/,1X,'INSERTING A BLANK RAINFALL SNAP, COUNT=',i4)

 22   format(8f10.0)

99999 CONTINUE

#ifdef CMPI
      call msg_fini()
#endif
      STOP 'OWI_RAIN READ ERROR (2)'

      END SUBROUTINE

C***********************************************************************
C   SUBROUTINE NCRAIN1_INTERP_BASIN
C
C   This generates and saves interpolation coefficients for mapping
C   from a basin-scale OWI to a ADCIRC grid.
C
C***********************************************************************

      SUBROUTINE NCRAIN1_INTERP_BASIN(NP)
      USE CONSTANTS, ONLY: RAD2DEG
      USE MESH, ONLY : SLAM, SFEA

      IMPLICIT NONE

      INTEGER   NP,I,J,K,XI,YI
      REAL(SZ)  adcLat,adcLong

      WRITE(16,*) ''
      WRITE(16,*) 'BASIN-SCALE RAINFALL MAPPING UPDATED'
      WRITE(16,*) ''

      iLatB = iLatc
      iLongB = iLongc
      dxB = dxc
      dyB = dyc
      swlatB = swlatc
      swlongB = swlongc

      ! Allocate and create matrices
      if(allocated(CrainB))    deallocate(CrainB)
      if(allocated(longB)) deallocate(longB)
      if(allocated(latB))  deallocate(latB)

      allocate(CrainB(iLongB,iLatB))
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
C   SUBROUTINE NCRAIN1_INTERP_REGION
C
C   This generates and saves interpolation coefficients for mapping
C   from a region-scale OWI to a ADCIRC grid.
C
C***********************************************************************

      SUBROUTINE NCRAIN1_INTERP_REGION(NP)
      USE CONSTANTS, ONLY: RAD2DEG
      USE MESH, ONLY : SLAM, SFEA
      IMPLICIT NONE

      INTEGER   NP,I,J,K,XI,YI
      REAL(SZ)  adcLat,adcLong

      WRITE(16,*) ''
      WRITE(16,*) 'REGION-SCALE RAINFALL MAPPING UPDATED'
      WRITE(16,*) ''

      iLatR = iLatc
      iLongR = iLongc
      dxR = dxc
      dyR = dyc
      swlatR = swlatc
      swlongR = swlongc

      ! Allocate and create matrices
      if(allocated(CrainR))    deallocate(CrainR)
      if(allocated(longR)) deallocate(longR)
      if(allocated(latR))  deallocate(latR)

      allocate(CrainR(iLongR,iLatR))
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


      END MODULE OWI_RAIN
