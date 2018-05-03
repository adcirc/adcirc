C******************************************************************************
C  MODULE OWIWIND
C    Written by s.b. 08/17/2006
C******************************************************************************
C
      MODULE OWIWIND
C
      USE SIZES,ONLY : SZ,MYPROC
      USE GLOBAL,ONLY : NSCREEN, ScreenUnit, DEBUG, ECHO, INFO,
     &   WARNING, ERROR, screenMessage, logMessage, allMessage,
     &   setMessageSource, unsetMessageSource, scratchMessage,
     &   openFileForRead
#ifdef CMPI
      USE MESSENGER,ONLY : MSG_FINI
#endif
      IMPLICIT NONE
      character*100 :: title

!      tcm v49.44 Changed the way the header line is read
!      character :: part1*12,part2*7,part3*6
      character :: owiheader*80

      real(SZ), dimension(:,:), allocatable :: uR,vR,pR,uB,vB,pB
      real(SZ), dimension(:), allocatable :: latR,longR,latB,longB
      real(SZ) :: Along, Alat
      real(SZ) :: ramp,rampfrac
      real(SZ) :: uu,vv,PP
      real(SZ) :: Penv
      ! start and end dates for region and basin scale data
      integer(8):: date1R,date2R,date1B,date2B
      integer(8):: date1,date2  ! generic start and end dates from OWI header
      integer(8):: date1w,date2w ! start and end dates for wind data
      integer(8):: date1p,date2p ! start and end dates for atm. press. data
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

      CHARACTER FNAME*1024 ! file to read from; used in error messages

      integer :: numSets,numBlankSnaps,cntSnaps,numSkipSnaps
      real(SZ) :: windMultiplier

      integer :: lun       ! fortran logical number to read from; used in
                           ! error messages
      integer :: errorIO   ! zero if the file opened successfully
C      logical :: fileFound ! true if the file is present
      character(len=100) :: errorVar ! name of variable that was being read
                                    ! when an error occurred

      logical :: moving_grid = .false. ! true will spatially interpolate each dataset

      PUBLIC

C---------------------end of data declarations--------------------------------C

      CONTAINS

C***********************************************************************
C   SUBROUTINE NWS12INIT
C***********************************************************************

      SUBROUTINE NWS12INIT(WVNX,WVNY,PRN,NP,RHOWAT0,G)

      USE SIZES, ONLY : SZ,MYPROC, GBLINPUTDIR

      IMPLICIT NONE
      INTEGER   NP,I
      REAL(SZ), intent(out), dimension(:) ::  WVNX,WVNY,PRN
      REAL(SZ)  RHOWAT0,RHOWATG,G
C
      call setMessageSource("nws12init")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      ! Allocate and create matrices  Jie 09/2013
      if(allocated(swpointsB))    deallocate(swpointsB)
      if(allocated(swpointsR))    deallocate(swpointsR)
      if(allocated(wB))           deallocate(wB)
      if(allocated(wR))           deallocate(wR)

      allocate(swpointsB(NP,2),wB(NP,4))
      allocate(swpointsR(NP,2),wR(NP,4))
      errorIO = 0

! Read meta info -------------------------------------------------
C
C     R E A D   F O R T 2 2
C     ---------------------
      lun = 22
      FNAME = TRIM(GBLINPUTDIR)//'/'//'fort.22'
      errorVar = ""
      call openFileForRead(lun, FNAME, errorIO)
      call check_err(errorIO)
      ! Read the number of sets of .pre and .win files from the fort.22.
      ! If numSets = 1 then ADCIRC requires UNIT 221 and 222.
      ! If numSets = 2 then ADCIRC requires UNIT 223 and 224
      !           in addition to 221 and 222.
      !   UNIT 221 and 223 are atmospheric pressure fields
      !   UNIT 222 and 224 are wind velocity fields.
      errorVar = "NWSET" ! used in error msgs
      read(lun,*,err=99999,end=99998,iostat=errorIO) numSets
      if((numSets.NE.1).AND.(numSets.NE.2)) then
         write(scratchMessage,1004) numSets
 1004    format("NWSET was '",I2,"' in unit 22.",
     &   " It must be set to 1 or 2. ADCIRC will stop.")
         call check_err(1)
      endif

      ! Read the number of blank snaps to be inserted before OWI winds start
      errorVar = "NWBS" ! used in error messages
      read(lun,*,err=99999,end=99998,iostat=errorIO) numBlankSnaps
      call check_err(errorIO)
      ! If numBlankSnaps < 0, ABS(numBlankSnaps) snaps in OWI wind files (UNIT221,222,223 and 224) will be skipped.
      if(numBlankSnaps.LT.0) then
         numSkipSnaps = ABS(numBlankSnaps)
         numBlankSnaps = 0

      ! v49.01 TCM 10/28/2009 -- Added else to initialize numSkipSnaps to be 0
      else
         numSkipSnaps = 0
      endif

      ! Read a wind velocity multiplier
      errorVar = "DWM" ! used in error messages
      read(lun,*,err=99999, end=99998, iostat=errorIO) windMultiplier
      call check_err(errorIO)
      close(lun)

! Read basin pre file header ------------------------------------------------

      lun = 221
      FNAME = TRIM(GBLINPUTDIR)//'/'//'fort.221'
      errorVar = ""
      call readHeader()
      date1p = date1
      date2p = date2

! Read basin win file header -------------------------------------------------

      lun = 222
      FNAME = TRIM(GBLINPUTDIR)//'/'//'fort.222'
      errorVar = ""
      call readHeader()
      date1w = date1
      date2w = date2

C
C     Error checking for basin scale data
      if ((date1p.ne.date1w).or.(date2p.ne.date2w)) then
         call allMessage(ERROR,
     &      "Start and end dates of basin met. data do not match.")
         errorVar = ""
         call check_err(1)
      endif
      date1B = date1p
      date2B = date2p

! Check if region scale data exist

      if (numSets.eq.2) then

! Read region pre file header -----------------------------------------------

         lun = 223
         FNAME = TRIM(GBLINPUTDIR)//'/'//'fort.223'
         errorVar = ""
         call readHeader()
         date1p = date1
         date2p = date2

! Read region win file header -----------------------------------------------

         lun = 224
         FNAME = TRIM(GBLINPUTDIR)//'/'//'fort.224'
         errorVar = ""
         call readHeader()
         date1w = date1
         date2w = date2
C
C     Error checking for region scale data
         if(date1p.ne.date1w.or.date2p.ne.date2w) then
            call allMessage(ERROR,
     &         "Start and end dates of region met. data do not match.")
            errorVar = ""
            call check_err(1)
         endif

         date1R = date1p
         date2R = date2p

      endif

! Initialize flags ----------------------------------------------------------
      isnapB = 0
      isnapR = 0
      updateB = 1
      updateR = 1
      cntSnaps = 0

! Skip snaps if necessary ---------------------------------------------------
      do i = 1,numSkipSnaps
         write(scratchMessage,41) i
  41     format("Skipping snap '",I6,"' in OWI wind data.")
         call logMessage(DEBUG,trim(scratchMessage))
         CALL NWS12GET(WVNX,WVNY,PRN,NP,RHOWAT0,G)
      enddo

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN

99998 call allMessage(ERROR,"Unexpectedly reached end-of-file.") ! END jumps here

99999 call check_err(1) !  ERR jumps here
C-----------------------------------------------------------------------
      END SUBROUTINE NWS12INIT
C-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!      S U B R O U T I N E   F I N D   S T O R M   C E N T E R
!-----------------------------------------------------------------------
!     Called to find the lowest pressure in the met field below
!     the threshold pressure.
!-----------------------------------------------------------------------
      subroutine findStormCenter(eyeLatR, eyeLonR, foundEye)
      implicit none
      real(sz), intent(out) :: eyeLatR(3)
      real(sz), intent(out) :: eyeLonR(3)
      logical, intent(out)  :: foundEye
      !Casey 110518: Added for Mark Powell's sector-based wind drag.
      INTEGER  :: EyeLatI
      INTEGER  :: EyeLonI
      REAL(SZ) :: EyeLatTemp
      REAL(SZ) :: EyeLonTemp
      REAL(SZ) :: EyePressure      
      integer :: i, j

      !Casey 110518: Find lon,lat location of eye.
      EyeLatI = 0
      EyeLonI = 0
      EyePressure = 1013.D0
      DO i=1,iLongB
         DO j=1,iLatB
            IF((pB(i,j).LT.EyePressure).AND.(pB(i,j).LT.1000.D0))THEN
               EyeLatI = j
               EyeLonI = i
               EyePressure = pB(i,j)
            ENDIF
         ENDDO
      ENDDO
      IF((EyeLatI.EQ.0).AND.(EyeLonI.EQ.0))THEN
         FoundEye = .FALSE.
      ELSEIF((EyeLatI.EQ.1).OR.(EyeLatI.EQ.iLatB).OR.(EyeLonI.EQ.1).OR.(EyeLonI.EQ.iLongB))THEN
         FoundEye = .FALSE.
      ELSE
         FoundEye = .TRUE.
         EyeLatTemp = swlatB  + (EyeLatI-1)*dyB
         EyeLonTemp = swlongB + (EyeLonI-1)*dxB
         IF((EyeLatTemp.EQ.EyeLatR(3)).AND.(EyeLonTemp.EQ.EyeLonR(3)))THEN
            CONTINUE
         ELSE
            EyeLatR(1) = EyeLatR(2)
            EyeLonR(1) = EyeLonR(2)
            EyeLatR(2) = EyeLatR(3)
            EyeLonR(2) = EyeLonR(3)
            EyeLatR(3) = EyeLatTemp
            EyeLonR(3) = EyeLonTemp
         ENDIF
      ENDIF

      !Casey 110518: Find lon,lat location of eye.
      EyeLatI = 0
      EyeLonI = 0
      EyePressure = 1013.D0
      DO i=1,iLongR
         DO j=1,iLatR
            IF((pR(i,j).LT.EyePressure).AND.(pR(i,j).LT.1000.D0))THEN
               EyeLatI = j
               EyeLonI = i
               EyePressure = pR(i,j)
            ENDIF
         ENDDO
      ENDDO
      IF((EyeLatI.EQ.0).AND.(EyeLonI.EQ.0))THEN
         IF(.NOT.FoundEye)THEN
            FoundEye = .FALSE.
         ELSE
            FoundEye = .TRUE.
         ENDIF
      ELSEIF((EyeLatI.EQ.1).OR.(EyeLatI.EQ.iLatR).OR.(EyeLonI.EQ.1).OR.(EyeLonI.EQ.iLongR))THEN
         IF(.NOT.FoundEye)THEN
            FoundEye = .FALSE.
         ELSE
            FoundEye = .TRUE.
         ENDIF
      ELSE
         FoundEye = .TRUE.
         EyeLatTemp = swlatR  + (EyeLatI-1)*dyR
         EyeLonTemp = swlongR + (EyeLonI-1)*dxR
         IF((EyeLatTemp.EQ.EyeLatR(3)).AND.(EyeLonTemp.EQ.EyeLonR(3)))THEN
            CONTINUE
         ELSE
            EyeLatR(1) = EyeLatR(2)
            EyeLonR(1) = EyeLonR(2)
            EyeLatR(2) = EyeLatR(3)
            EyeLonR(2) = EyeLonR(3)
            EyeLatR(3) = EyeLatTemp
            EyeLonR(3) = EyeLonTemp
         ENDIF
      ENDIF
!-----------------------------------------------------------------------
      end subroutine findStormCenter
!-----------------------------------------------------------------------


C***********************************************************************
C   SUBROUTINE NWS12GET
C***********************************************************************

      SUBROUTINE NWS12GET(WVNX,WVNY,PRN,NP,RHOWAT0,G)
      USE SIZES,ONLY : MYPROC,MNPROC
      IMPLICIT NONE
      INTEGER   NP,I,J,XI,YI
      REAL(SZ), intent(out), dimension(:) ::  WVNX,WVNY,PRN
      REAL(SZ)  RHOWAT0,RHOWATG,G
      CHARACTER*80 PBLJAGF

      !jgf51.52.05: Used to detect mesh nodes outside both met grids
      ! and set those nodes to zero wind velocity and to background pressure.
      LOGICAL :: uvpHasBeenSet 
C
      call setMessageSource("nws12get")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      RHOWATG=RHOWAT0*G

! Set a value for Penv

      Penv=1013d0         !mb

! Read basin data  ---------------------------------------------------------

      ! Increment counter (cntSnaps initialized to zero in nws12init)
      cntSnaps = cntSnaps+1

      ! Put a blank snap for the first 'numBlankSnaps' snaps and then return
      if (cntSnaps.LE.numBlankSnaps) then
         do i=1,NP
            WVNX(I)=0.d0
            WVNY(I)=0.d0
            PRN(I)=101300.d0/RHOWATG
         enddo
         !TCM v49.02 (Changed format number from 15 to 16)
         write(scratchMessage,16) cntSnaps
 16      format('INSERTING A BLANK WIND SNAP, COUNT=',i4)
         call allMessage(INFO,trim(scratchMessage))
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         RETURN
      endif

      ! Increment counter
      isnapB = isnapB+1

      ! Read grid specifications/date in basin pressure file
      errorVar = "grid specifications/date in basin pressure file"
      read (221,11,end=10000,err=9999,iostat=errorIO)
     &      iLatp,iLongp,dxp,dyp,swlatp,swlongp,iCYMDHp,iMinp
      call check_err(errorIO)

      ! Read grid specifications/date in basin wind file
      errorVar = "grid specifications/date in basin wind file"
      read (222,11,end=10000,err=9999,iostat=errorIO)
     &      iLatw,iLongw,dxw,dyw,swlatw,swlongw,iCYMDHw,iMinw
      call check_err(errorIO)

      ! Check consistency
      if(iLatp.ne.iLatw.or.iLongp.ne.iLongw.or.dxp.ne.dxw.or.
     &     dyp.ne.dyw.or.swlatp.ne.swlatw.or.swlongp.ne.swlongw.or.
     &     iCYMDHp.ne.iCYMDHw.or.iMinp.ne.iMinw) then
         call allMessage(ERROR,
     &  "Grid specifications/date in OWI win and pre files must match.")
         errorVar = ""
         call check_err(1)
      endif

      ! Check if header info has changed from the previous snapshot
      if(isnapB.gt.1) then
        if(iLatp.ne.iLatB.or.iLongp.ne.iLongB.or.dxp.ne.dxB.or.
     &       dyp.ne.dyB.or.swlatp.ne.swlatB.or.
     &       swlongp.ne.swlongB) then
          call logMessage(INFO,"Basin scale grid has changed.")
          updateB = 1
        else
          updateB = 0
        endif
      endif

      ! for a moving grid,update it every time  Jie 09/2013
      if (moving_grid) updateB = 1

      iCYMDHB = iCYMDHp
      iMinB = iMinp

      ! Update coordinate mapping coefficients if necessary
      if(updateB.eq.1) then
         call logMessage(INFO,
     &       "Updating basin grid coordinate mapping coefficients.")
         call nws12interp_basin(np)
      endif

      ! Read basin scale atmospheric pressure snapshot
      errorVar = "basin scale atmospheric pressure snapshot"
      read(221,22,end=10000,err=9999,iostat=errorIO)
     &   ((pB(i,j),i=1,iLongB),j=1,iLatB)
      call check_err(errorIO)

      ! Read basin scale snapshot of u/v components of the wind
      errorVar = "basin scale wind u-velocity snapshot"
      read(222,22,end=10000,err=9999,iostat=errorIO)
     &   ((uB(i,j),i=1,iLongB),j=1,iLatB)
      call check_err(errorIO)
      errorVar = "basin scale wind v-velocity snapshot"
      read(222,22,end=10000,err=9999,iostat=errorIO)
     &   ((vB(i,j),i=1,iLongB),j=1,iLatB)
      call check_err(errorIO)


! Read region data  --------------------------------------------------------

      regionExists = .FALSE.

      IF(numSets.EQ.1) GOTO 100
      if(iCYMDHB.lt.date1R) goto 100
      if(iCYMDHB.eq.date2R.and.iMinR.ne.0) goto 100
      if(iCYMDHB.gt.date2R) goto 100

      regionExists = .TRUE.

      ! Increment counter
      isnapR = isnapR+1

      ! Read grid specifications/date in region pressure file
      errorVar = "grid specifications/date in region pressure file"
      read (223,11,end=10000,err=9999,iostat=errorIO)
     &     iLatp,iLongp,dxp,dyp,swlatp,swlongp,iCYMDHp,iMinp
      call check_err(errorIO)

      ! Read grid specifications/date in region wind file
      errorVar = "grid specifications/date in region wind file"
      read (224,11,end=10000,err=9999,iostat=errorIO)
     &        iLatw,iLongw,dxw,dyw,swlatw,swlongw,iCYMDHw,iMinw
      call check_err(errorIO)

      if(iLatp.ne.iLatw.or.iLongp.ne.iLongw.or.dxp.ne.dxw.or.
     &     dyp.ne.dyw.or.swlatp.ne.swlatw.or.swlongp.ne.swlongw.or.
     &     iCYMDHp.ne.iCYMDHw.or.iMinp.ne.iMinw) then
         call allMessage(ERROR,
     &   "Grid specfications/date in OWI win and pre files must match.")
         errorVar = ""
         call check_err(1)
      endif

      ! Check if header info has changed from the previous snapshot
      if(isnapR.gt.1) then
         if(iLatp.ne.iLatR.or.iLongp.ne.iLongR.or.dxp.ne.dxR.or.
     &       dyp.ne.dyR.or.swlatp.ne.swlatR.or.
     &       swlongp.ne.swlongR) then
            call logMessage(INFO,"Region scale grid has changed.")
            updateR = 1
         else
            updateR = 0
         endif

      endif

      ! for a moving grid,update it every time  Jie 09/2013
      if (moving_grid) updateR = 1
      
      iCYMDHR = iCYMDHp
      iMinR = iMinp

      if (iCYMDHB.ne.iCYMDHR.or.iMinB.ne.iMinR) then
         call allMessage(ERROR,"Snapshots not synchronized.")
         write(scratchMessage,51) iCYMDHB, iMinB, iCYMDHR, iMinR
  51     format("Basin snapshot date is '",I10,"' and '",
     &      I2,"' minutes. Region snapshot date is '",I10,
     &      "' and '",I2,"' minutes.")
         errorVar = ""
         call check_err(1)
      endif

      ! Update coordinate mapping coefficients if necessary
      if (updateR.eq.1) then
         call logMessage(INFO,
     &       "Updating region grid coordinate mapping coefficients.")
         call nws12interp_region(np)
      endif

      ! Read pressure
      errorVar = "region scale atmospheric pressure snapshot"
      read(223,22,end=10000,err=9999,iostat=errorIO)
     &   ((pR(i,j),i=1,iLongR),j=1,iLatR)
      call check_err(errorIO)

      ! Read u/v components of the wind
      errorVar = "region scale wind u-velocity snapshot"
      read(224,22,end=10000,err=9999,iostat=errorIO)
     &   ((uR(i,j),i=1,iLongR),j=1,iLatR)
      call check_err(errorIO)
      errorVar = "region scale wind v-velocity snapshot"
      read(224,22,end=10000,err=9999,iostat=errorIO)
     &   ((vR(i,j),i=1,iLongR),j=1,iLatR)
      call check_err(errorIO)


 100  CONTINUE

! Interpolate onto ADCIRC grid and write to file -------------------------

      rampfrac = isnapB-1
c      if (rampfrac<36) then
c        ramp = tanh(18d0*rampfrac/36d0)
c      end if
      ramp = 1.0

      if (regionExists.EQV..TRUE.) then
         write(scratchMessage,15) iCYMDHB,iMinB
 15      format("Processing basin scale wind data ",I12," ",I2,".")
         call allMessage(INFO,trim(scratchMessage))
      else
         write(scratchMessage,14) iCYMDHB,iMinB
 14      format("Processing region scale wind data ",I12," ",I2,".")
         call allMessage(INFO,trim(scratchMessage))
      endif

      do i=1,NP

         ! jgf51.52.05: Fixed an issue with floating point .eq. that
         ! was used to determine whether the value of uu had been set. 
         ! Replaced floating point comparison with a logical variable.
         uvpHasBeenSet = .false.

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

            uvpHasBeenSet = .true.
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

            uvpHasBeenSet = .true.
        endif

        ! COPY TO ARRAYS ---------------------------------------------------

         if (uvpHasBeenSet.eqv..false.) then
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

         endif

      enddo

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN


9999  call check_err(1) ! ERR during read jumps to here

      ! END during read of required data jumps to here
9997  call allMessage(ERROR,"Unexpected end-of-file reached.")
      call check_err(1)

10000 continue  ! END during read of u, v, p data jumps to here
      write(scratchMessage,61) trim(errorVar)
  61  format("Unexpected end-of-file while reading '",A,
     &    "'. Wind speeds set to zero and pressure to 1013mb.")
      call allMessage(WARNING,trim(scratchMessage))

      WVNX=0.d0
      WVNY=0.d0
      PRN=101300.d0/RHOWATG

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN


 11   format(t6,i4,t16,i4,t23,f6.0,t32,f6.0,
     &     t44,f8.0,t58,f8.0,t69,i10,i2)

 22   format(8f10.0)

      END SUBROUTINE NWS12GET

C***********************************************************************
C   SOBROUTINE NWS12INTERP_BASIN
C
C   This generates and saves interpolation coefficients for mapping
C   from a basin-scale OWI to a ADCIRC grid.
C
C***********************************************************************

      SUBROUTINE NWS12INTERP_BASIN(NP)
      USE GLOBAL, ONLY : RAD2DEG
      USE MESH, ONLY : SLAM, SFEA

      IMPLICIT NONE

      INTEGER   NP,I,J,K,XI,YI
      REAL(SZ)  adcLat,adcLong
C
      call setMessageSource("nws12interp_basin")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
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

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
      END SUBROUTINE NWS12INTERP_BASIN


C***********************************************************************
C   SOBROUTINE NWS12INTERP_REGION
C
C   This generates and saves interpolation coefficients for mapping
C   from a region-scale OWI to a ADCIRC grid.
C
C***********************************************************************

      SUBROUTINE NWS12INTERP_REGION(NP)
      USE GLOBAL, ONLY : RAD2DEG
      USE MESH, ONLY : SLAM, SFEA
      IMPLICIT NONE

      INTEGER   NP,I,J,K,XI,YI
      REAL(SZ)  adcLat,adcLong

      call setMessageSource("nws12interp_region")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
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

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
      END SUBROUTINE NWS12INTERP_REGION

C-----------------------------------------------------------------------
C     S U B R O U T I N E   R E A D   H E A D E R
C-----------------------------------------------------------------------
C     jgf50.16 Read the header from the OWI file and sets the value
C     of the date1p, date2p, date1w, or date2w variables as appropriate.
C-----------------------------------------------------------------------
      SUBROUTINE readHeader()
      IMPLICIT NONE
C
      call setMessageSource("readHeader")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      call openFileForRead(lun, FNAME, errorIO)
      call check_err(errorIO)
C
      ! Read begining/ending dates of pre file
!      tcm v49.44 changed the way the header is read
!      read(221,11,err=99999,ADVANCE="NO")part1,part2,part3
! 11   FORMAT(a,1x,a,1x,a)
!      read(221,*)date1p,date2p
      owiheader(:) = ' '  !set owiheader to blanks before read
      errorVar = "owiheader"
      read(lun, fmt='(a80)',
     &   end=99998,err=99999,iostat=errorIO) owiheader
      call check_err(errorIO)
C
      errorVar = "start date"
      read(owiheader(56:65),'(I10)',
     &   end=99998,err=99999,iostat=errorIO) date1
      call check_err(errorIO)
      write(scratchMessage,31) trim(errorVar), trim(FNAME), date1
  31  format("'",A,"' in  '",A,"' is '",I10,"'.")
      call allMessage(INFO,scratchMessage)
C
      errorVar = "end date"
      read(owiheader(71:80),'(I10)',
     &    end=99998,err=99999,iostat=errorIO) date2
      call check_err(errorIO)
      write(scratchMessage,31) trim(errorVar), trim(FNAME), date2
      call allMessage(INFO,scratchMessage)
C
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN

99998 call allMessage(ERROR,"Unexpectedly reached end-of-file.") ! END jumps here
99999 call check_err(1) !  ERR jumps here
C-----------------------------------------------------------------------
      end subroutine readHeader
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   C H E C K  _  E R R
C-----------------------------------------------------------------------
C     jgf50.16 Checks the return value from subroutine calls; if there
C     was an error, it writes a termination message to the screen and to the
C     fort.16 file and terminates ADCIRC.
C-----------------------------------------------------------------------
      subroutine check_err(iret)
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      IMPLICIT NONE
      INTEGER, intent(in) :: iret
      call setMessageSource("check_err")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      if (iret.ne.0) then
         if (errorVar.ne."") then
            write(scratchMessage,888)
     &         trim(errorVar), trim(FNAME), lun
 888        format("Failed to read '",A,"' from '",A,
     &         "' (unit number ",I3,").")
            call allMessage(ERROR,trim(scratchMessage))
         endif
         call allMessage(ERROR,"ADCIRC execution terminated.")
#ifdef CMPI
         call msg_fini()
#endif
         stop
      endif
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C-----------------------------------------------------------------------
      end subroutine check_err
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      END MODULE OWIWIND
C-----------------------------------------------------------------------

