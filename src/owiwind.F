C******************************************************************************
C  MODULE OWIWIND
C    Written by s.b. 08/17/2006
C    Modified by TCM 08/08/2019 -- Added LandFall Scale win/pre (finest level)
C                                   Basin ->   Region ->  LandFall
C                      File Units (221,222)   (223,224)   (217,218)
C                      Added PRdeflt to nws12init and nws12get calls
C                      in order to use background pressure value instead of
C                      hard wired 1013.0 mb.
C                      Added the OWI_FNAMES namelist to the fort.22 file
C                      in order to specify non-default file names for the
C                      fort.221,fort.222,fort.223,fort.224,fort.217, fort.218
C                      BASIN_PRE_FNAME, (added by ast)
C******************************************************************************
C
      MODULE OWIWIND
C
      USE SIZES,ONLY : SZ,MYPROC
      USE GLOBAL,ONLY : NSCREEN, ScreenUnit, DEBUG, ECHO, INFO,
     &   WARNING, ERROR, screenMessage, logMessage, allMessage,
     &   setMessageSource, unsetMessageSource, scratchMessage,
     &   openFileForRead, Flag_ElevError
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
      real(SZ), dimension(:,:), allocatable :: uLF,vLF,pLF
      real(SZ), dimension(:), allocatable :: latLF,longLF
      real(SZ) :: Along, Alat
      real(SZ) :: ramp,rampfrac
      real(SZ) :: uu,vv,PP
      real(SZ) :: Penv
      ! start and end dates for region and basin scale data
      integer(8):: date1R,date2R,date1B,date2B
      integer(8):: date1LF,date2LF
      integer(8):: date1,date2  ! generic start and end dates from OWI header
      integer(8):: date1w,date2w ! start and end dates for wind data
      integer(8):: date1p,date2p ! start and end dates for atm. press. data
      integer :: iLatR,iLongR,iCYMDHR,iMinR
      real(SZ)  :: dxR,dyR,swlatR,swlongR
      integer :: iLatLF,iLongLF,ICYMDHLF,iMinLF
      real(SZ)  :: dxLF,dyLF,swlatLF,swlongLF      
      integer :: iLatB,iLongB,iCYMDHB,iMinB
      real(SZ)  :: dxB,dyB,swlatB,swlongB
      integer :: iLatw,iLongw,iCYMDHw,iMinw
      real(SZ)  :: dxw,dyw,swlatw,swlongw
      integer :: iLatp,iLongp,iCYMDHp,iMinp
      real(SZ)  :: dxp,dyp,swlatp,swlongp
      integer :: isnapR,updateR
      integer :: isnapB,updateB
      integer :: isnapLF,updateLF

      logical :: regionExists
      logical :: LandFallExists

      integer,allocatable :: swpointsR(:,:)
      integer,allocatable :: swpointsB(:,:)
      integer,allocatable :: swpointsLF(:,:)

      real(SZ) :: w,w1,w2,w3,w4
      real(SZ),allocatable :: wR(:,:)
      real(SZ),allocatable :: wB(:,:)
      real(SZ),allocatable :: wLF(:,:)

      CHARACTER FNAME*1024 ! file to read from; used in error messages

      integer :: numSets,numBlankSnaps,cntSnaps,numSkipSnaps
      real(SZ) :: windMultiplier

      integer :: lun       ! fortran logical number to read from; used in
                           ! error messages
      integer :: errorIO   ! zero if the file opened successfully
C      logical :: fileFound ! true if the file is present
      character(len=100) :: errorVar ! name of variable that was being read
                                    ! when an error occurred
C     FileName Specification (Optional)
      CHARACTER(len=1024) :: BASIN_PRE_FNAME,BASIN_WND_FNAME
      CHARACTER(len=1024) :: REGION_PRE_FNAME,REGION_WND_FNAME
      CHARACTER(len=1024) :: LANDFALL_PRE_FNAME,LANDFALL_WND_FNAME

      logical :: moving_grid = .false. ! true will spatially interpolate each dataset
C      logical :: DragLawString = 
      PUBLIC

C---------------------end of data declarations--------------------------------C

      CONTAINS

C***********************************************************************
C   SUBROUTINE NWS12INIT
C***********************************************************************

      SUBROUTINE NWS12INIT(WVNX,WVNY,PRN,NP,RHOWAT0,G,PRdeflt)

      USE SIZES, ONLY : SZ,MYPROC, GBLINPUTDIR

      IMPLICIT NONE
      INTEGER   NP,I
      REAL(SZ), intent(out), dimension(:) ::  WVNX,WVNY,PRN
      REAL(SZ)  RHOWAT0,RHOWATG,G,PRdeflt
      LOGICAL :: found_owifname_nml = .false.   !flag to determine if the OWI_FNAME namelist was present
      character(len=1024) :: tmpstr

C     TCM 20170129 ADDED Optional Namelist for File Names
      NAMELIST /OWI_FNAMES/ BASIN_PRE_FNAME,BASIN_WND_FNAME,
     &             REGION_PRE_FNAME,REGION_WND_FNAME,
     &             LANDFALL_PRE_FNAME,LANDFALL_WND_FNAME

      
C
      call setMessageSource("nws12init")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

C      ! Allocate and create matrices  Jie 09/2013
C      if(allocated(swpointsB))    deallocate(swpointsB)
C      if(allocated(swpointsR))    deallocate(swpointsR)
C      if(allocated(wB))           deallocate(wB)
C      if(allocated(wR))           deallocate(wR)
C
C      allocate(swpointsB(NP,2),wB(NP,4))
C      allocate(swpointsR(NP,2),wR(NP,4))

      errorIO = 0

! Read meta info -------------------------------------------------
C
C     R E A D   F O R T 2 2
C     ---------------------
      lun = 22
      FNAME = TRIM(GBLINPUTDIR)//'/'//'fort.22'
      errorVar = ""
      call openFileForRead(lun, FNAME, errorIO)
      call check_owi_err(errorIO)
      ! Read the number of sets of .pre and .win files from the fort.22.
      ! If numSets = 1 then ADCIRC requires UNIT 221 and 222.
      ! If numSets = 2 then ADCIRC requires UNIT 223 and 224
      !           in addition to 221 and 222.
      ! If numSets = 3 then ADCIRC requires UNIT 217 and 218
      !           in addition to 221, 222, 223, and 224.
      !   UNIT 221 and 223 (217) are atmospheric pressure field
      !   UNIT 222 and 224 (218) are wind velocity fields.

      errorVar = "NWSET" ! used in error msgs
      read(lun,*,err=99999,end=99998,iostat=errorIO) numSets
      if((numSets.NE.1).AND.(numSets.NE.2).AND.(numSets.NE.3)) then
         write(scratchMessage,1004) numSets
 1004    format("NWSET was '",I2,"' in unit 22.",
     &   " It must be set to 1, 2 or 3. ADCIRC will stop.")
         call check_owi_err(errorIO)
      endif

      ! Read the number of blank snaps to be inserted before OWI winds start
      errorVar = "NWBS" ! used in error messages
      read(lun,*,err=99999,end=99998,iostat=errorIO) numBlankSnaps
      call check_owi_err(errorIO)
      ! If numBlankSnaps < 0, ABS(numBlankSnaps) snaps
      ! in OWI wind files (UNIT 221,222,223, 224 and 217 & 218) will be skipped.
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
      call check_owi_err(errorIO)

      ! Check for Optional Namelist for FileNames
      BASIN_PRE_FNAME(:)= "fort.221"
      BASIN_WND_FNAME(:)= "fort.222"
      REGION_PRE_FNAME(:)= "fort.223"
      REGION_WND_FNAME(:)= "fort.224"
      LANDFALL_PRE_FNAME(:)= "fort.217"
      LANDFALL_WND_FNAME(:)= "fort.218"

      errorIO=0
      errorVar = "OWI_FNAMES" !used in error messages
      read(unit=lun,NML=OWI_FNAMES,iostat=errorIO)

      if (errorIO.gt.0) then  !error condition
         backspace(lun)
         read(lun,fmt='(A)') ErrorVar
         tmpstr="Invalid line in namelist = OWI_FNAMES (fort.22) ::" // trim(ErrorVar)
         call allMessage(ERROR,trim(tmpstr) )
         call check_owi_err(1)
         found_owifname_nml = .true.
      elseif(errorIO.eq.0) then
         found_owifname_nml = .true.
         call logMessage(INFO,"Found namelist = OWI_FNAMES in the fort.22 file")
      else
         call logMessage(INFO,"Namelist = OWI_FNAMES not found in the fort.22 file")
         found_owifname_nml = .false.
      endif
      
      close(lun)

! Read basin pre file header ------------------------------------------------

      lun = 221
      FNAME = TRIM(GBLINPUTDIR)//'/'//trim(BASIN_PRE_FNAME)
      errorVar = ""
      call readHeader()
      date1p = date1
      date2p = date2

! Read basin win file header -------------------------------------------------

      lun = 222
      FNAME = TRIM(GBLINPUTDIR)//'/'//trim(BASIN_WND_FNAME)
      errorVar = ""
      call readHeader()
      date1w = date1
      date2w = date2

C
C     Error checking for basin scale data
      if ((date1p.ne.date1w).or.(date2p.ne.date2w)) then
         call allMessage(ERROR,"Start and end dates of basin " //
     &  " (fort.221 and fort.222) meteorological data do not match.")
         errorVar = ""
         call check_owi_err(errorIO)
      endif
      date1B = date1p
      date2B = date2p

! Check if region scale data exist

      if ((numSets.eq.2).or.(numSets.eq.3)) then

! Read region pre file header -----------------------------------------------

         lun = 223
         FNAME = TRIM(GBLINPUTDIR)//'/'//trim(REGION_PRE_FNAME)
         errorVar = ""
         call readHeader()
         date1p = date1
         date2p = date2

! Read region win file header -----------------------------------------------

         lun = 224
         FNAME = TRIM(GBLINPUTDIR)//'/'//trim(REGION_WND_FNAME)
         errorVar = ""
         call readHeader()
         date1w = date1
         date2w = date2
C
C     Error checking for region scale data
         if(date1p.ne.date1w.or.date2p.ne.date2w) then
            call allMessage(ERROR,"Start and end dates of region " //
     &       "(fort.223 and fort.224) meteorological data do not match.")
            errorVar = ""
            call check_owi_err(errorIO)
         endif

         date1R = date1p
         date2R = date2p

      endif

! Check if LandFall scale data exist (this is the finest scale)

      if (numSets.eq.3) then

! Read LandFall pre file header -----------------------------------------------

         lun = 217
         FNAME = TRIM(GBLINPUTDIR)//'/'//trim(LANDFALL_PRE_FNAME)
         errorVar = ""
         call readHeader()
         date1p = date1
         date2p = date2

! Read LandFall win file header -----------------------------------------------

         lun = 218
         FNAME = TRIM(GBLINPUTDIR)//'/'//trim(LANDFALL_WND_FNAME)
         errorVar = ""
         call readHeader()
         date1w = date1
         date2w = date2
C
C     Error checking for LandFall scale data
         if(date1p.ne.date1w.or.date2p.ne.date2w) then
            call allMessage(ERROR,
     &         "Start and end dates of LandFall met. data do not match.")
            errorVar = ""
            call check_owi_err(1)
         endif

         date1LF = date1p
         date2LF = date2p

      endif  !LandFall Data
! ---------------------------------
!
!   Now that we know how many sets of field files to read
!   Allocate the needed arrays
      ! Allocate and create matrices  Jie 09/2013
      if(allocated(swpointsB))    deallocate(swpointsB)
      if(allocated(swpointsR))    deallocate(swpointsR)
      if(allocated(swpointsLF))   deallocate(swpointsLF)
      if(allocated(wB))           deallocate(wB)
      if(allocated(wR))           deallocate(wR)
      if(allocated(wLF))          deallocate(wLF)

      allocate(swpointsB(NP,2),wB(NP,4))
      if (numSets.gt.1) allocate(swpointsR(NP,2),wR(NP,4))
      if (numSets.gt.2) allocate(swpointsLF(NP,2),wLF(NP,4))

!--------------------------------------------
      
! Initialize flags ----------------------------------------------------------
      isnapB = 0
      isnapR = 0
      isnapLF = 0
      updateB = 1
      updateR = 1
      updateLF = 1
      cntSnaps = 0

! Skip snaps if necessary ---------------------------------------------------
      do i = 1,numSkipSnaps
         write(scratchMessage,41) i
  41     format("Skipping snap '",I6,"' in OWI wind data.")
         call logMessage(DEBUG,trim(scratchMessage))
         CALL NWS12GET(WVNX,WVNY,PRN,NP,RHOWAT0,G,PRdeflt)
      enddo

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN

99998 call allMessage(ERROR,"Unexpectedly reached end-of-file.") ! END jumps here

99999 call check_owi_err(errorIO) !  ERR jumps here
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
      REAL(SZ) :: EyePressure, PRdeflt
      integer :: i, j

      !Casey 110518: Find lon,lat location of eye.
      EyeLatI = 0
      EyeLonI = 0
      EyePressure = PRdeflt  !1013.D0  changed to default background pressure
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

      SUBROUTINE NWS12GET(WVNX,WVNY,PRN,NP,RHOWAT0,G,PRdeflt)
      USE SIZES,ONLY : MYPROC,MNPROC,GBLINPUTDIR
      IMPLICIT NONE
      real(SZ), dimension(:) :: eyeLatR(3)
      real(SZ), dimension(:) :: eyeLonR(3)
      INTEGER   NP,I,J,XI,YI
      REAL(SZ), intent(out), dimension(:) ::  WVNX,WVNY,PRN
      REAL(SZ)  RHOWAT0,RHOWATG,G,PRdeflt,EyePressure
      CHARACTER*80 PBLJAGF
      CHARACTER*80 DragLawString
      INTEGER :: EyeLatI,EyeLonI
      INTEGER :: EyeLatTemp,EyeLonTemp

      !jgf51.52.05: Used to detect mesh nodes outside both met grids
      ! and set those nodes to zero wind velocity and to background pressure.
      LOGICAL :: uvpHasBeenSet,FoundEye

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
            PRN(I)=PRdeflt*100.d0/RHOWATG
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
      lun=221
      read (lun,11,end=10000,err=9999,iostat=errorIO)
     &      iLatp,iLongp,dxp,dyp,swlatp,swlongp,iCYMDHp,iMinp
      call check_owi_err(errorIO)

      ! Read grid specifications/date in basin wind file
      lun=222
      errorVar = "grid specifications/date in basin wind file"
      read (lun,11,end=10000,err=9999,iostat=errorIO)
     &      iLatw,iLongw,dxw,dyw,swlatw,swlongw,iCYMDHw,iMinw
      call check_owi_err(errorIO)

      ! Check consistency
      if(iLatp.ne.iLatw.or.iLongp.ne.iLongw.or.dxp.ne.dxw.or.
     &     dyp.ne.dyw.or.swlatp.ne.swlatw.or.swlongp.ne.swlongw.or.
     &     iCYMDHp.ne.iCYMDHw.or.iMinp.ne.iMinw) then
         call allMessage(ERROR,
     &  "Grid specifications/date in OWI win and pre files must match.")
         errorVar = ""
         call check_owi_err(errorIO)
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
      lun=221
      FNAME = TRIM(GBLINPUTDIR)//'/'//'fort.221'
      read(lun,22,end=10000,err=9999,iostat=errorIO)
     &   ((pB(i,j),i=1,iLongB),j=1,iLatB)
      call check_owi_err(errorIO)

      ! Read basin scale snapshot of u/v components of the wind
      errorVar = "basin scale wind u-velocity snapshot"
      lun=222
      FNAME = TRIM(GBLINPUTDIR)//'/'//'fort.222'      
      read(lun,22,end=10000,err=9999,iostat=errorIO)
     &   ((uB(i,j),i=1,iLongB),j=1,iLatB)
      call check_owi_err(errorIO)
      errorVar = "basin scale wind v-velocity snapshot"
      read(lun,22,end=10000,err=9999,iostat=errorIO)
     &   ((vB(i,j),i=1,iLongB),j=1,iLatB)
      call check_owi_err(errorIO)


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
      lun=223
      read (lun,11,end=10000,err=9999,iostat=errorIO)
     &     iLatp,iLongp,dxp,dyp,swlatp,swlongp,iCYMDHp,iMinp
      call check_owi_err(errorIO)

      ! Read grid specifications/date in region wind file
      errorVar = "grid specifications/date in region wind file"
      lun=224    
      read (lun,11,end=10000,err=9999,iostat=errorIO)
     &        iLatw,iLongw,dxw,dyw,swlatw,swlongw,iCYMDHw,iMinw
      call check_owi_err(errorIO)

      if(iLatp.ne.iLatw.or.iLongp.ne.iLongw.or.dxp.ne.dxw.or.
     &     dyp.ne.dyw.or.swlatp.ne.swlatw.or.swlongp.ne.swlongw.or.
     &     iCYMDHp.ne.iCYMDHw.or.iMinp.ne.iMinw) then
         call allMessage(ERROR,
     &   "Grid specfications/date in OWI win and pre files must match.")
         errorVar = ""
         call check_owi_err(errorIO)
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
         call check_owi_err(errorIO)
      endif

      ! Update coordinate mapping coefficients if necessary
      if (updateR.eq.1) then
         call logMessage(INFO,
     &       "Updating region grid coordinate mapping coefficients.")
         call nws12interp_region(np)
      endif

      ! Read pressure
      errorVar = "region scale atmospheric pressure snapshot"
      lun=223   
      read(lun,22,end=10000,err=9999,iostat=errorIO)
     &   ((pR(i,j),i=1,iLongR),j=1,iLatR)
      call check_owi_err(errorIO)

      ! Read u/v components of the wind
      errorVar = "region scale wind u-velocity snapshot"
      lun=224     
      read(lun,22,end=10000,err=9999,iostat=errorIO)
     &   ((uR(i,j),i=1,iLongR),j=1,iLatR)
      call check_owi_err(errorIO)
      errorVar = "region scale wind v-velocity snapshot"
      read(224,22,end=10000,err=9999,iostat=errorIO)
     &   ((vR(i,j),i=1,iLongR),j=1,iLatR)
      call check_owi_err(errorIO)

! --- End region wind/pre reads

! Read LandFall data  --------------------------------------------------------

      LandFallExists = .FALSE.

      IF(numSets.LT.3) GOTO 100
      if(iCYMDHB.lt.date1LF) goto 100
      if(iCYMDHB.eq.date2LF.and.iMinLF.ne.0) goto 100
      if(iCYMDHB.gt.date2LF) goto 100

      LandFallExists = .TRUE.

      ! Increment counter
      isnapLF = isnapLF+1

      ! Read grid specifications/date in LandFall pressure file
      errorVar = "grid specifications/date in LandFall pressure file"
      lun=217
      read (lun,11,end=10000,err=9999,iostat=errorIO)
     &     iLatp,iLongp,dxp,dyp,swlatp,swlongp,iCYMDHp,iMinp
      call check_owi_err(errorIO)

      ! Read grid specifications/date in LandFall wind file
      errorVar = "grid specifications/date in LandFall wind file"
      lun=218
      read (lun,11,end=10000,err=9999,iostat=errorIO)
     &        iLatw,iLongw,dxw,dyw,swlatw,swlongw,iCYMDHw,iMinw
      call check_owi_err(errorIO)

      if(iLatp.ne.iLatw.or.iLongp.ne.iLongw.or.dxp.ne.dxw.or.
     &     dyp.ne.dyw.or.swlatp.ne.swlatw.or.swlongp.ne.swlongw.or.
     &     iCYMDHp.ne.iCYMDHw.or.iMinp.ne.iMinw) then
         call allMessage(ERROR,
     &   "Grid specfications/date in OWI win and pre files must match.")
         errorVar = ""
         call check_owi_err(1)
      endif
      ! Check if header info has changed from the previous snapshot
      if(isnapLF.gt.1) then
         if(iLatp.ne.iLatR.or.iLongp.ne.iLongR.or.dxp.ne.dxR.or.
     &       dyp.ne.dyR.or.swlatp.ne.swlatR.or.
     &       swlongp.ne.swlongR) then
            call logMessage(INFO,"LandFall scale grid has changed.")
            updateLF = 1
         else
            updateLF = 0
         endif

      endif

      ! for a moving grid,update it every time  Jie 09/2013
      if (moving_grid) updateLF = 1

      iCYMDHLF = iCYMDHp
      iMinLF = iMinp

      if (iCYMDHB.ne.iCYMDHLF.or.iMinB.ne.iMinLF) then
         call allMessage(ERROR,"Snapshots not synchronized.")
         write(scratchMessage,81) iCYMDHB, iMinB, iCYMDHLF, iMinLF
  81     format("Basin snapshot date is '",I10,"' and '",
     &      I2,"' minutes. LandFall snapshot date is '",I10,
     &      "' and '",I2,"' minutes.")
         errorVar = ""
         call check_owi_err(1)
      endif

      ! Update coordinate mapping coefficients if necessary
      if (updateLF.eq.1) then
         call logMessage(INFO,
     &       "Updating LandFall grid coordinate mapping coefficients.")
         call nws12interp_landfall(np)
      endif

      ! Read pressure
      errorVar = "LandFall scale atmospheric pressure snapshot"
      lun=217
      read(lun,22,end=10000,err=9999,iostat=errorIO)
     &   ((pLF(i,j),i=1,iLongLF),j=1,iLatLF)
      call check_owi_err(errorIO)

      ! Read u/v components of the wind
      errorVar = "LandFall scale wind u-velocity snapshot"
      lun=218
      read(lun,22,end=10000,err=9999,iostat=errorIO)
     &   ((uLF(i,j),i=1,iLongLF),j=1,iLatLF)
      call check_owi_err(errorIO)
      errorVar = "LandFall scale wind v-velocity snapshot"
      read(lun,22,end=10000,err=9999,iostat=errorIO)
     &   ((vLF(i,j),i=1,iLongLF),j=1,iLatLF)
      call check_owi_err(errorIO)

      ! jgf50.32: Made this if block dependent on the drag law formulation.
      IF (trim(DragLawString) == "Powell" .or.
     &    trim(DragLawString) == "POWELL" .or.
     &    trim(DragLawString) == "powell" ) THEN
         !Casey 110518: Find lon,lat location of eye.
         EyeLatI = 0
         EyeLonI = 0
         EyePressure = PRdeflt  !1013.D0  changed to default background pressure
         DO i=1,iLongLF
            DO j=1,iLatLF
               IF((pLF(i,j).LT.EyePressure).AND.(pLF(i,j).LT.1000.D0))THEN
                  EyeLatI = j
                  EyeLonI = i
                  EyePressure = pLF(i,j)
               ENDIF
            ENDDO
         ENDDO
         IF((EyeLatI.EQ.0).AND.(EyeLonI.EQ.0))THEN
            IF(.NOT.FoundEye)THEN
               FoundEye = .FALSE.
            ELSE
               FoundEye = .TRUE.
            ENDIF
         ELSEIF((EyeLatI.EQ.1).OR.(EyeLatI.EQ.iLatLF).OR.(EyeLonI.EQ.1).OR.(EyeLonI.EQ.iLongLF))THEN
            IF(.NOT.FoundEye)THEN
               FoundEye = .FALSE.
            ELSE
               FoundEye = .TRUE.
            ENDIF
         ELSE
            FoundEye = .TRUE.
            EyeLatTemp = swlatLF  + (EyeLatI-1)*dyLF
            EyeLonTemp = swlongLF + (EyeLonI-1)*dxLF
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
      ENDIF

! --- End LandFall wind/pre reads
!-----------------------
      
      

 100  CONTINUE

! Interpolate onto ADCIRC grid and write to file -------------------------

      rampfrac = isnapB-1
c      if (rampfrac<36) then
c        ramp = tanh(18d0*rampfrac/36d0)
c      end if
      ramp = 1.d0

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

         ! LANDFALL--------------------------------------------------------
         ! uu, vv and PP will be overwritten if region data exist.

         if ((LandFallExists).and.(swpointsLF(i,1).gt.0)) then
            xi = swpointsLF(i,1)
            yi = swpointsLF(i,2)

            w1=wLF(i,1)
            w2=wLF(i,2)
            w3=wLF(i,3)
            w4=wLF(i,4)

            uu=w1*uLF(xi,yi)+w2*uLF(xi+1,yi)+
     &         w3*uLF(xi+1,yi+1)+w4*uLF(xi,yi+1)
            vv=w1*vLF(xi,yi)+w2*vLF(xi+1,yi)+
     &         w3*vLF(xi+1,yi+1)+w4*vLF(xi,yi+1)
            PP=w1*pLF(xi,yi)+w2*pLF(xi+1,yi)+
     &         w3*pLF(xi+1,yi+1)+w4*pLF(xi,yi+1)

            uvpHasBeenSet = .true.
        endif


        ! COPY TO ARRAYS ---------------------------------------------------

         if (uvpHasBeenSet.eqv..false.) then
            WVNX(I)=0.d0
            WVNY(I)=0.d0
            PRN(I)=PRdeflt*100.d0/RHOWATG
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


9999  call check_owi_err(errorIO) ! ERR during read jumps to here

      ! END during read of required data jumps to here
9997  call allMessage(ERROR,"Unexpected end-of-file reached.")
      call check_owi_err(errorIO)

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
C   SUBROUTINE NWS12INTERP_LANDFALL
C
C   This generates and saves interpolation coefficients for mapping
C   from a LandFall-scale OWI to a ADCIRC grid.
C
C***********************************************************************

      SUBROUTINE NWS12INTERP_LANDFALL(NP)
      USE GLOBAL, ONLY : RAD2DEG
      USE MESH, ONLY : SLAM, SFEA
      IMPLICIT NONE

      INTEGER   NP,I,J,K,XI,YI
      REAL(SZ)  adcLat,adcLong

      call setMessageSource("nws12interp_landfall")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      WRITE(16,*) ''
      WRITE(16,*) 'LandFall-SCALE WIND MAPPING UPDATED'
      WRITE(16,*) ''

      iLatLF = iLatw
      iLongLF = iLongw
      dxLF = dxw
      dyLF = dyw
      swlatLF = swlatw
      swlongLF = swlongw

      ! Allocate and create matrices
      if(allocated(uLF))    deallocate(uLF)
      if(allocated(vLF))    deallocate(vLF)
      if(allocated(pLF))    deallocate(pLF)
      if(allocated(longLF)) deallocate(longLF)
      if(allocated(latLF))  deallocate(latLF)

      allocate(uLF(iLongLF,iLatLF),vLF(iLongLF,iLatLF))
      allocate(pLF(iLongLF,iLatLF))
      allocate(longLF(iLongLF),latLF(iLatLF))

      ! Generate long&lat on each grid point
      do i=1,iLatLF
        latLF(i) = swlatLF+(i-1)*dyLF
      enddo
      do i=1,iLongLF
        longLF(i) = swlongLF+(i-1)*dxLF
      enddo

      ! Generate interpolation info (south west point and weights)
      do i=1,NP
        adcLat  = RAD2DEG*SFEA(i)
        adcLong = RAD2DEG*SLAM(i)


        if (adcLong>=longLF(1).and.adcLong<longLF(iLongLF) .and.
     &       adcLat>=latLF(1).and.adcLat<latLF(iLatLF)) then
          do j=1,iLongLF-1
            if (adcLong>=longLF(j).and.adcLong<longLF(j+1)) then
              xi=j
              goto 201
            endif
          enddo
 201      continue

          do k=1,iLatLF-1
            if (adcLat>=latLF(k).and.adcLat<latLF(k+1)) then
              yi=k
              goto 301
            endif
          enddo
 301      continue

          swpointsLF(i,1) = xi
          swpointsLF(i,2) = yi

          w=(longLF(xi+1)-longLF(xi))*(latLF(yi+1)-latLF(yi))
          w1=(longLF(xi+1)-adcLong)*(latLF(yi+1)-adcLat)
          w2=(adcLong-longLF(xi))*(latLF(yi+1)-adcLat)
          w3=(adcLong-longLF(xi))*(adcLat-latLF(yi))
          w4=(longLF(xi+1)-adcLong)*(adcLat-latLF(yi))

          wLF(i,1)=w1/w
          wLF(i,2)=w2/w
          wLF(i,3)=w3/w
          wLF(i,4)=w4/w
        else
          swpointsLF(i,1) = 0
          swpointsLF(i,2) = 0
        endif
      enddo

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
      END SUBROUTINE NWS12INTERP_LandFall

        

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
      call check_owi_err(errorIO)
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
      call check_owi_err(errorIO)
C
      errorVar = "start date"
      read(owiheader(56:65),'(I10)',
     &   end=99998,err=99999,iostat=errorIO) date1
      call check_owi_err(errorIO)
      write(scratchMessage,31) trim(errorVar), trim(FNAME), date1
  31  format("'",A,"' in  '",A,"' is '",I10,"'.")
      call allMessage(INFO,scratchMessage)
C
      errorVar = "end date"
      read(owiheader(71:80),'(I10)',
     &    end=99998,err=99999,iostat=errorIO) date2
      call check_owi_err(errorIO)
      write(scratchMessage,31) trim(errorVar), trim(FNAME), date2
      call allMessage(INFO,scratchMessage)
C
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN

99998 call allMessage(ERROR,"Unexpectedly reached end-of-file.") ! END jumps here
99999 call check_owi_err(errorIO) !  ERR jumps here
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
      subroutine check_owi_err(iret)
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      IMPLICIT NONE
      INTEGER, intent(in) :: iret
      call setMessageSource("check_owi_err")
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
         Flag_ElevError = .true.  !overloading this flag to help kill all mpi proc
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
      end subroutine check_owi_err
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      END MODULE OWIWIND
C-----------------------------------------------------------------------

