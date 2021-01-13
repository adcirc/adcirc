!-----------------------------------------------------------------------
!  MODULE OWIWIND
!-----------------------------------------------------------------------
!> @author Zachary Cobell, The Water Institute, zcobell@thewaterinstitute.org
!> @author Shintaro Bunya
!> 
!> @copyright Dr. R.A. Luettich and Dr. J.J. Westerink
!>
!> @brief This module handles the I/O and interpolation to the mesh for Oceanweather
!> format wind and pressure fields. 
!> 
!> The module is initialized by calling nws12init and subsequently calling 
!> nws12get. The code will place data into the WVNX, WVNY, and PRN arrays
!> 
!> The fort.22 file is used to control behavior. fort.22 must contain at least
!> three lines as follows:
!>     numSets        <-- Number of sets of wind fields
!>     numBlankSnaps  <-- Number of snaps to add/skip at beginning of field
!>     windMultiplier <-- Wind velocity multiplier to apply
!>
!> If numSets is negative, the code will attempt to read the filenames of 
!> the wind files from the fort.22. numSets may be positive 1, 2, or 3 which
!> indicates use of the standard names:
!>
!> @verbatim
!>    Set 1: fort.221, fort.222 [ pressure field, wind field ]
!>    Set 2: fort.223, fort.224 [ pressure field, wind field ]
!>    Set 3: fort.217, fort.218 [ pressure field, wind field ]
!> @endverbatim
!>
!> When negative, there is no enforced limit on the number of wind fields that 
!> may be used. The wind fields should be listed in comma separated format
!> in the fort.22, 1:abs(numSets), one set of wind files per line as follows:
!>
!> @verbatim
!> pressureFile_001.pre, windFile_001.wnd
!> pressureFile_002.pre, windFile_002.wnd
!> ...
!> ...
!> @endverbatim
!>
!> Wind fields should always be listed from coarse to fine. The priority with
!> which the wind fields are interpolated to the mesh is determined by the 
!> reverse order for which they are listed, meaning that the wind field that
!> is listed last is considered highest priority.
!>
!> Wind fields are interpolated onto each node using a bilinear interpolation
!> scheme. If a wind field lies outside of all wind fields, it is set to 
!> PRDEFLT (usually 1013mb) and 0.0 m/s velocity.
!>
!> When using Powell wind drag with Oceanweather wind fields, the findStormCenter
!> subroutine is used to search for the lowest central pressure under 1000mb. In
!> some cases, the storm center may vary randomly enough within an area when using a 
!> simple search like this on grids with high spatial and temporal resolution. If 
!> this happens, it is possible that the storm direction will not be correctly 
!> interpreted. In this case, it usually is best to search only the more coarse
!> wind fields to find the storm center which will prevent these sorts of oscillations
!> Using the &metControl namelist, you can specify the number of Powell drag
!> search domains via:
!> @verbatim
!> &metControl nPowellSearchDomains=[X] /
!> @endverbatim
!> where [X] is the number of domains to search through. If set to 1, only the 
!> most coarse domain will participiate in the search.
!
!-----------------------------------------------------------------------
      MODULE OWIWIND
!-----------------------------------------------------------------------

        USE SIZES,ONLY : SZ,MYPROC
        USE GLOBAL,ONLY : NSCREEN, ScreenUnit, DEBUG, ECHO, INFO,
     &     WARNING, ERROR, screenMessage, logMessage, allMessage,
     &     setMessageSource, unsetMessageSource, scratchMessage,
     &     openFileForRead, Flag_ElevError
#ifdef CMPI
        USE MESSENGER,ONLY : MSG_FINI
#endif  
        IMPLICIT NONE

        integer  :: nPowellSearchDomains = -1 !< Number of search domains for Powell wind drag
        integer  :: numSets !< Number of sets of wind and pressure fields
        integer  :: numBlankSnaps !< Number of blank snaps to prepend to Oceanweather files
        integer  :: numSkipSnaps !< Number of snaps to skip in Oceanweather files
        integer  :: cntSnaps !< Counter for number of snaps that have been processed so far
        real(SZ) :: windMultiplier !< Wind multiplier applied to Oceanweather wind velicities
        logical  :: moving_grid = .false. !< indicates the owi domains use a moving grid (for aswip.F)
        
        TYPE OCEANWEATHER !< Container for Oceanweather domains. Represents a wind and pressure field
            CHARACTER(1024)      :: pressure_file !< filename used for the pressure field
            CHARACTER(1024)      :: wind_file !< filename used for the wind velocity field
            INTEGER(8)           :: startDate !< start date in the domain
            INTEGER(8)           :: endDate !< end date in the domain
            INTEGER(8)           :: date !< current date in domain for data read
            INTEGER              :: minute !< current minute in domain for data read
            INTEGER              :: iLat !< number of points in latitude direction
            INTEGER              :: iLon !< number of points in longitude direction
            INTEGER              :: iounit_pressure !< fortran unit number for pressure file
            INTEGER              :: iounit_wind !< fortran unit number for wind file
            INTEGER              :: iSnap !< current snap number 
            INTEGER              :: iUpdate !< flag to update interpolation weights and indicies if wind field moves
            REAL(SZ)             :: dLat !< spacing in latitude direction
            REAL(SZ)             :: dLon !< spacing in longitude direction
            REAL(SZ)             :: swLon !< southwest longitude
            REAL(SZ)             :: swLat !< southwest latitude
            REAL(SZ)             :: adcircStartTime !< start time of adcirc simulation
            REAL(SZ)             :: windTimeElapsed !< time that has elapsed in the wind fields
            LOGICAL              :: hasData !< set to false when fields are no longer have data
            LOGICAL              :: atEnd !< set to false when the fields are at the end of file
            REAL(SZ),ALLOCATABLE :: p(:,:) !< pressure field
            REAL(SZ),ALLOCATABLE :: u(:,:) !< u-velocity field
            REAL(SZ),ALLOCATABLE :: v(:,:) !< v-velocity field
            REAL(SZ),ALLOCATABLE :: weight(:,:) !< bilinear interpolation weights
            INTEGER,ALLOCATABLE  :: interpPoints(:,:) !< indicies for interpolation
        END TYPE OCEANWEATHER

        TYPE(OCEANWEATHER),ALLOCATABLE :: owi(:) !< Array containing each of the wind and pressure field used

        PUBLIC nws12get,nws12init,findStormCenter,windMultiplier
        PUBLIC nPowellSearchDomains,moving_grid

        PRIVATE nws12interp,cntSnaps,numBlankSnaps,numSkipSnaps
        PRIVATE owi_allocate,owi_deallocate,owi_generateCoordinates
        PRIVATE owi_setFilenames,owi_setDefaultFilenames,owi_setDynamicFilenames
        PRIVATE owi_initializeFileHeaders,owi_allocateInterpolationArrays
        PRIVATE check_owi_err,readHeader,owi_skipSnaps,owi,numSets
        PRIVATE owi_readNextSnap,owi_interpolateToMesh,owi_setEndOfFileBlankSnap

!---------------------end of data declarations--------------------------------

        CONTAINS

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! P U B L I C  F U N C T I O N S
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> Initializes reading data from the oceanweather wind fields
!> @param[inout] wvnx    wind speed, x-direction
!> @param[inout] wvny    wind speed, y-direction
!> @param[inout] prn     atmospheric pressure
!> @param[in]    np      number of nodes in the adcirc mesh for this processor
!> @param[in]    prdeflt background atmospheric pressure
!-----------------------------------------------------------------------
      SUBROUTINE NWS12INIT(WVNX,WVNY,PRN,NP,PRdeflt)
!-----------------------------------------------------------------------
        USE SIZES, ONLY : SZ,MYPROC, GBLINPUTDIR
        IMPLICIT NONE
        REAL(SZ), intent(out) :: WVNX(:),WVNY(:),PRN(:)
        INTEGER,intent(in)    :: NP
        INTEGER               :: I
        REAL(SZ)              :: PRdeflt
        character(len=1024)   :: tmpstr
        character(len=1024)   :: errorVar = ""
        character(len=1024)   :: fname,f1,f2
        integer               :: iounit
        integer               :: errorIO = 0

        call setMessageSource("nws12init")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif  

        !...Begin reading fort.22
        fname = "fort.22"
        iounit = 22
        errorVar = "open"
        call openFileForRead(iounit, TRIM(GBLINPUTDIR)//"/fort.22", errorIO)
        call check_owi_err(errorIO,errorVar,fname,iounit)
        
        errorVar = "reading number of wind sets" ! used in error msgs
        read(iounit,*,err=99999,end=99998,iostat=errorIO) numSets
        call check_owi_err(errorIO,errorVar,fname,iounit)

        if(numSets.EQ.0)then
            call check_owi_err(1,"number of wind sets must not be zero",fname,iounit)
        endif

        ! Read the number of blank snaps to be inserted before OWI winds start
        errorVar = "reading number of blank snaps" ! used in error messages
        read(iounit,*,err=99999,end=99998,iostat=errorIO) numBlankSnaps
        call check_owi_err(errorIO, "read number of blank snaps",fname,iounit)
        ! If numBlankSnaps < 0, ABS(numBlankSnaps) snaps
        ! in OWI wind files (UNIT 221,222,223, 224 and 217 & 218) will be skipped.
        if(numBlankSnaps.LT.0) then
           numSkipSnaps = ABS(numBlankSnaps)
           numBlankSnaps = 0
        else
           numSkipSnaps = 0 
        endif

        ! Read a wind velocity multiplier
        errorVar = "reading wind multiplier" ! used in error messages
        read(iounit,*,err=99999, end=99998, iostat=errorIO) windMultiplier
        call check_owi_err(errorIO,errorVar,fname,iounit)

        errorIO=0

        ! Allocate the owi arrays
        allocate(owi(1:abs(numSets)))
        call owi_setFilenames(owi,iounit,numSets)
        close(iounit)

        ! set the files to contain data
        do i = 1,abs(numSets)
            owi(i)%atEnd = .false.
            owi(i)%hasData = .false.
        enddo

        ! Check the number of powell search domains
        IF(nPowellSearchDomains.GT.abs(numSets))THEN
            nPowellSearchDomains = abs(numSets)
        ELSEIF(nPowellSearchDomains.LE.0)THEN
            nPowellSearchDomains = abs(numSets)
        ENDIF

        !...Initialize the wind objects based on file headers
        call owi_initializeFileHeaders(owi)

        !...Allocate interpolation arrays
        call owi_allocateInterpolationArrays(owi,np)

        !...Set snap counter to zero
        cntSnaps = 0
       
        !...Skip specified snaps 
        call owi_skipSnaps(numSkipSnaps,np,WVNX,WVNY,PRN,PRdeflt)

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.")
#endif  
        call unsetMessageSource()
        RETURN

99998   call allMessage(ERROR,"Unexpectedly reached end-of-file.") ! END jumps here

99999   call check_owi_err(errorIO,errorVar,fname,iounit) !  ERR jumps here
!-----------------------------------------------------------------------
      END SUBROUTINE NWS12INIT
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Reads a wind and pressure snap from the oceanweather file and 
!>        interpolates it to the adcirc mesh
!> @param[inout] wvnx    wind speed, x-direction
!> @param[inout] wvny    wind speed, y-direction
!> @param[inout] prn     atmospheric pressure
!> @param[in]    np      number of nodes in the adcirc mesh for this processor
!> @param[in]    prdeflt background atmospheric pressure
!> @param[out]   hasdata (otional), index of where OWI data is assigned
!-----------------------------------------------------------------------
      SUBROUTINE NWS12GET(WVNX,WVNY,PRN,NP,PRdeflt,HasData)
!-----------------------------------------------------------------------
        USE SIZES,ONLY: MYPROC,MNPROC,GBLINPUTDIR
        USE GLOBAL,ONLY: RHOWAT0, G, NWS
        IMPLICIT NONE
        INTEGER,intent(in)     :: NP
        REAL(SZ),intent(out)   :: WVNX(:),WVNY(:),PRN(:)
        REAL(SZ),intent(in)    :: PRdeflt
        REAL(SZ)               :: RHOWATG,EyePressure
        INTEGER                :: I,J,S,XI,YI
        CHARACTER(1024)        :: filename,errorVar
        INTEGER                :: iounit
        REAL(SZ)               :: uu,vv,pp
        REAL(SZ)               :: w1,w2,w3,w4
        INTEGER                :: errorIO
        LOGICAL                :: uvpHasBeenSet
        LOGICAL,intent(out),optional :: HasData(:)

        call setMessageSource("nws12get")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif  
        errorVar = ""
        RHOWATG=RHOWAT0*G

        ! Increment counter (cntSnaps initialized to zero in nws12init)
        cntSnaps = cntSnaps+1

        ! Put a blank snap for the first 'numBlankSnaps' snaps and then return
        if (cntSnaps.LE.numBlankSnaps) then
           if (abs(NWS).ne.14) then ! WJP: Don't overwrite if NWS = 14
              WVNX(:) = 0D0
              WVNY(:) = 0D0
              PRN(:) = Prdeflt*100D0/RHOWATG
           endif
           if (present(HasData)) then 
              HasData(:) = .false.
           endif
           write(scratchMessage,16) cntSnaps
 16        format('INSERTING A BLANK WIND SNAP, COUNT=',i4)
           call allMessage(INFO,trim(scratchMessage))
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
           call allMessage(DEBUG,"Return.")
#endif  
           call unsetMessageSource()
           RETURN
        endif

        DO S = 1,numSets

          if(s.gt.1)then
              if((owi(1)%date.lt.owi(s)%startDate) .or.
     &          (owi(1)%date.eq.owi(s)%endDate.and.owi(s)%minute.ne.0) .or.
     &          (owi(1)%date.gt.owi(s)%endDate)) then
                    owi(s)%hasData = .false.
                    cycle
              endif
          endif

          if(owi(s)%atEnd)then
                !  Set the blank snap and report to screen if past end
                !  of wind data
                call owi_setEndOfFileBlankSnap(owi(s)%pressure_file,
     &                                         prdeflt,owi(s))
          else
                !  Read the next snap from the wind and pressure files
                call owi_readNextSnap(s,np,prdeflt)

                ! Increment counter
                owi(s)%isnap = owi(s)%isnap + 1
          endif

        ENDDO

        do i=1,NP
           uvpHasBeenSet = .false.

           call owi_interpolateToMesh(i,uu,vv,pp,uvpHasBeenSet)

           if (.not.uvpHasBeenSet) then
              if (abs(NWS).ne.14) then ! WJP: Don't overwrite if NWS = 14
                 WVNX(I)=0.d0
                 WVNY(I)=0.d0
                 PRN(I)=PRdeflt*100.d0/RHOWATG
              endif
           else
              !Convert millibars to m of water
              PRN(i)  = 100.d0*PP/RHOWATG

              ! Apply wind velocity multiplier
              ! and save to arrays
              WVNX(i) = uu * windMultiplier
              WVNY(i) = vv * windMultiplier

           endif
              
           if (present(HasData)) then 
              HasData(i) = uvpHasBeenSet
           endif

        enddo

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.")
#endif  
        call unsetMessageSource()
        RETURN

!-----------------------------------------------------------------------
      END SUBROUTINE NWS12GET
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Called to find the lowest pressure in the met field below
!>        the threshold pressure. Loops over all domains and finds 
!>        the lowest pressure across all domains
!> @param[inout] eyeLat prior 3 locations (y) where the eye was found
!> @param[inout] eyeLon prior 3 locations (x) where the eye was found
!> @param[out]   foundEye set to true if the eye has been found
!-----------------------------------------------------------------------
      subroutine findStormCenter(eyeLat, eyeLon, foundEye, outfoundRMax, outRMaxLon, outRMaxLat)
!-----------------------------------------------------------------------
        use GLOBAL,only:PRBCKGRND
        implicit none
        intrinsic               :: minloc
        real(sz), intent(inout) :: eyeLat(3)
        real(sz), intent(inout) :: eyeLon(3)
        logical, intent(out)    :: foundEye
        logical, intent(out),optional    :: outFoundRMax
        real(sz), intent(out),optional   :: outRmaxLon, outRMaxLat
        integer                 :: EyeLatI
        integer                 :: EyeLonI
        integer                 :: EyeDomainI
        integer                 :: maxWindLonI
        integer                 :: maxWindLatI
        real(SZ)                :: EyeLatTemp
        real(SZ)                :: EyeLonTemp
        real(SZ)                :: EyePressure
        real(SZ)                :: windMag,windMax
        real(SZ)                :: rMaxLon,rMaxLat
        real(SZ)                :: minp
        integer                 :: s,i,j
        integer                 :: minidx(2)
        logical                 :: FoundRMax = .FALSE.

        call setMessageSource("findStormCenter")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif

        EyeLatI = 0
        EyeLonI = 0
        EyeDomainI = 0
        EyeLatTemp = 0D0
        EyeLonTemp = 0D0
        windMax = 0D0
        EyePressure = PRBCKGRND
        
        do S = 1,nPowellSearchDomains
            if(owi(s)%hasData)then
                minidx = minloc(owi(s)%p)
                minp = owi(s)%p(minidx(1),minidx(2))
                if(minp.lt.1000D0.and.minp.lt.eyePressure)then
                    EyeLonI = minidx(1)
                    EyeLatI = minidx(2)
                    EyeDomainI = s
                    EyePressure = owi(s)%p(minidx(1),minidx(2))
                endif

                do i=1,owi(s)%iLon
                    do j=1,owi(s)%iLat
                        windMag = (owi(s)%u(i,j)**2D0 + owi(s)%v(i,j)**2D0)**(1D0/2D0)
                        if(windMag.GT.windMax)then
                            windMax         = windMag
                            maxWindLonI     = i
                            maxWindLatI     = j
                        endif
                    enddo
                enddo

            endif
        enddo
        
        IF(EyeDomainI.NE.0)THEN
            IF((EyeLatI.EQ.0).OR.(EyeLonI.EQ.0).OR.(EyeDomainI.EQ.0).OR.
     &         (EyeLatI.EQ.1).OR.(EyeLatI.EQ.owi(EyeDomainI)%iLat).OR.
     &             (EyeLonI.EQ.1).OR.(EyeLonI.EQ.owi(EyeDomainI)%iLon))THEN
               FoundEye = .FALSE.
            ELSE
               FoundEye = .TRUE.
               EyeLatTemp = owi(EyeDomainI)%swlat+(EyeLatI-1)*owi(EyeDomainI)%dLat
               EyeLonTemp = owi(EyeDomainI)%swlon+(EyeLonI-1)*owi(EyeDomainI)%dLon
               rMaxLat = owi(EyeDomainI)%swlat+(maxWindLatI-1)*owi(EyeDomainI)%dLat
               rMaxLon = owi(EyeDomainI)%swlon+(maxWindLonI-1)*owi(EyeDomainI)%dLon 
               FoundRMax = .TRUE.
               IF((EyeLatTemp.EQ.EyeLat(3)).AND.(EyeLonTemp.EQ.EyeLon(3)))THEN
                  CONTINUE
               ELSE
                  EyeLat(1) = EyeLat(2)
                  EyeLon(1) = EyeLon(2)
                  EyeLat(2) = EyeLat(3)
                  EyeLon(2) = EyeLon(3)
                  EyeLat(3) = EyeLatTemp
                  EyeLon(3) = EyeLonTemp
               ENDIF
            ENDIF
        ELSE
             FoundEye = .FALSE.
             FoundRMax = .FALSE.
        ENDIF

        IF(PRESENT(outFoundRMax)) outFoundRMax = foundRMax
        IF(PRESENT(outRMaxLon)) outRMaxLon = rMaxLon
        IF(PRESENT(outRMaxLat)) outRMaxLat = rMaxLat

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.")
#endif  
        call unsetMessageSource()

!-----------------------------------------------------------------------
      end subroutine findStormCenter
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! P R I V A T E  F U N C T I O N S
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
          
!-----------------------------------------------------------------------
!> @brief reads the next data snapshot from the owi wind and pressure
!> files 
!> @param[in] idx domain index for the owi data
!> @param[in] np number of points in the adcirc mesh
!> @param[in] prbk background pressure
!-----------------------------------------------------------------------
      subroutine owi_readNextSnap(idx,np,prbk)
!-----------------------------------------------------------------------
        implicit none
        integer,intent(in)  :: idx
        integer,intent(in)  :: np
        real(sz),intent(in) :: prbk 
        character(1024)     :: filename
        character(1024)     :: errorVar
        integer             :: i,j
        integer             :: iounit,errorIO
        integer             :: iLatw,iLongw,iLatp,iLongp
        integer(8)          :: iCYMDHw,iCYMDHp
        integer             :: iMinw,iMinp
        real(SZ)            :: dxw,dyw,swlatw,swlongw
        real(SZ)            :: dxp,dyp,swlatp,swlongp
        
        call setMessageSource("owi_readNextSnap")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif
        
        owi(idx)%hasData = .false.

        ! Read grid specifications/date in pressure file
        filename = owi(idx)%pressure_file
        iounit = owi(idx)%iounit_pressure
        write(errorVar,'(A,I0)') "read grid specifications/date in Oceanweather pressure file, domain ",idx
        read (owi(idx)%iounit_pressure,11,end=10000,err=9999,iostat=errorIO)
     &        iLatp,iLongp,dxp,dyp,swlatp,swlongp,iCYMDHp,iMinp
        call check_owi_err(errorIO,errorVar,filename,iounit)

        if(idx.eq.1)then
            if(numSets.GT.1)then
                write(scratchMessage,'("Processing Oceanweather wind data for domain ",I0," through "'//
     &                               ',I0," for time stamp ",I12," ",I2)') idx,numSets,iCYMDHp,iMinp
            else
                write(scratchMessage,'("Processing Oceanweather wind data for time stamp ",I12," ",I2)') 
     &              iCYMDHp,iMinp
            endif
            call allmessage(INFO,scratchMessage)
        endif

        ! Read grid specifications/date in wind file
        filename = owi(idx)%wind_file
        iounit = owi(idx)%iounit_wind
        WRITE(errorVar,'(A,I0)') "read grid specifications/date in Oceanweather wind file, domain ",idx
        read (owi(idx)%iounit_wind,11,end=10000,err=9999,iostat=errorIO)
     &        iLatw,iLongw,dxw,dyw,swlatw,swlongw,iCYMDHw,iMinw
        call check_owi_err(errorIO,errorVar,filename,iounit)

        ! Check consistency
        if(iLatp.ne.iLatw.or.iLongp.ne.iLongw.or.dxp.ne.dxw.or.
     &       dyp.ne.dyw.or.swlatp.ne.swlatw.or.swlongp.ne.swlongw.or.
     &       iCYMDHp.ne.iCYMDHw.or.iMinp.ne.iMinw) then
           call allMessage(ERROR,
     &    "Grid specifications/date in OWI win and pre files must match.")
           errorVar = ""
           call check_owi_err(errorIO,errorVar,filename,iounit)
        endif

        if(owi(idx)%isnap.gt.1) then
          if(iLatp.ne.owi(idx)%iLat.or.iLongp.ne.owi(idx)%iLon.or.
     &         dxp.ne.owi(idx)%dlon.or.dyp.ne.owi(idx)%dlat.or.
     &         swlatp.ne.owi(idx)%swlat.or.swlongp.ne.owi(idx)%swlon) then
            write(scratchMessage,'(A,I0,A)') "Oceanweather domain ",idx," has changed."
            call logMessage(INFO,TRIM(scratchMessage))
            owi(idx)%iupdate = 1
          else
            owi(idx)%iupdate = 0
          endif
        endif

        owi(idx)%date = iCYMDHp
        owi(idx)%minute = iMinp

        ! Update coordinate mapping coefficients if necessary
        if(owi(idx)%iupdate.eq.1.or.moving_grid) then
           owi(idx)%ilon = ilongp
           owi(idx)%ilat = ilatp
           owi(idx)%dlon = dxp
           owi(idx)%dlat = dyp
           owi(idx)%swlat = swlatp
           owi(idx)%swlon = swlongp
           WRITE(errorVar,'(A,I0)') "Updating grid coordinate mapping coefficients, domain ",idx
           call logMessage(INFO,errorVar)
           call nws12interp(idx,np)
        endif

        ! Read basin scale atmospheric pressure snapshot
        filename = owi(idx)%pressure_file
        iounit =  owi(idx)%iounit_pressure
        WRITE(errorVar,'(A,I0)') "reading atmospheric pressure snapshot, domain ",idx
        read(owi(idx)%iounit_pressure,22,end=10000,err=9999,iostat=errorIO)
     &     ((owi(idx)%p(i,j),i=1,owi(idx)%iLon),j=1,owi(idx)%iLat)
        call check_owi_err(errorIO,errorVar,filename,iounit)

        ! Read basin scale snapshot of u/v components of the wind
        filename = owi(idx)%wind_file
        iounit = owi(idx)%iounit_wind
        WRITE(errorVar,'(A,I0)') "reading wind u-velocity snapshot, domain ",idx
        read(owi(idx)%iounit_wind,22,end=10000,err=9999,iostat=errorIO)
     &     ((owi(idx)%u(i,j),i=1,owi(idx)%iLon),j=1,owi(idx)%iLat)
        call check_owi_err(errorIO,errorVar,filename,iounit)
        filename = owi(idx)%wind_file
        iounit = owi(idx)%iounit_wind
        WRITE(errorVar,'(A,I0)') "reading wind v-velocity snapshot, domain ",idx
        read(owi(idx)%iounit_wind,22,end=10000,err=9999,iostat=errorIO)
     &     ((owi(idx)%v(i,j),i=1,owi(idx)%iLon),j=1,owi(idx)%iLat)
        call check_owi_err(errorIO,errorVar,filename,iounit)

        ! Set this flag so interpolation routine knows
        ! to use this domain or not
        owi(idx)%hasData = .true.

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.")
#endif  
        call unsetMessageSource()

        return

! Oceanweather read formats
 11     format(t6,i4,t16,i4,t23,f6.0,t32,f6.0,
     &       t44,f8.0,t58,f8.0,t69,i10,i2)
 22     format(8f10.0)

! Error IO
9999    call check_owi_err(errorIO,errorVar,filename,iounit) ! ERR during read jumps to here

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.")
#endif  
        call unsetMessageSource()
        return

! END during read of required data jumps to here
9997    call allMessage(ERROR,"Unexpected end-of-file reached.")
        call check_owi_err(errorIO,errorVar,filename,iounit)

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.")
#endif  
        call unsetMessageSource()
        return

! END during read of grid, u, v, p data jumps to here
10000   continue  
        call owi_setEndOfFileBlankSnap(filename,prbk,owi(idx))

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.")
#endif  
        call unsetMessageSource()
        return

!-----------------------------------------------------------------------
      end subroutine owi_readNextSnap
!-----------------------------------------------------------------------
      
!-----------------------------------------------------------------------
!> @brief Sets the blank snaps at the end of the file when data has run
!>        out
!> @param[in] filename name of file that is being read
!> @param[in] prbk background pressure 
!> @param[inout] owi oceanweather object array
!-----------------------------------------------------------------------
      subroutine owi_setEndOfFileBlankSnap(filename,prbk,obj)
!-----------------------------------------------------------------------
        implicit none
        CHARACTER(*),INTENT(IN)          :: filename
        REAL(SZ),INTENT(IN)              :: prbk
        TYPE(OCEANWEATHER),INTENT(INOUT) :: obj

        call setMessageSource("owi_setEndOfFileBlankSnap")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif        

        write(scratchMessage,61) trim(filename)
  61    format("End-of-file reached while reading '",A,
     &      "'. Wind speeds set to zero and the pressure field",
     &      " to its background value.")
        call allMessage(WARNING,trim(scratchMessage))

        ! Set no wind forcing for this domain
        obj%u(:,:)  = 0D0
        obj%v(:,:)  = 0D0
        obj%p(:,:)  = prbk
        obj%hasData = .false.
        obj%atEnd   = .true.

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.")
#endif  
        call unsetMessageSource()

        return
!-----------------------------------------------------------------------
      end subroutine owi_setEndOfFileBlankSnap
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief performs the bilinear interpolation with precomputed
!> coefficients for an adcirc node
!> @param[in] node adcirc node index
!> @param[out] u interpolated u-velocity
!> @param[out] v interpolated v-velocity
!> @param[out] p interpolated pressure
!> @param[out] isSet true if this node was interpolated from the owi data
!-----------------------------------------------------------------------
      subroutine owi_interpolateToMesh(node,u,v,p,isSet)
!-----------------------------------------------------------------------
        implicit none
        integer,intent(in)   :: node
        real(sz),intent(out) :: u,v,p
        logical,intent(out)  :: isSet
        integer              :: j,xi,yi
        real(sz)             :: w1,w2,w3,w4

        call setMessageSource("owi_interpolateToMesh")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif

        isSet = .false.

        do j = numSets,1,-1
           if(.not.owi(j)%hasData)cycle
           if (owi(j)%interppoints(node,1).gt.0) then
              xi = owi(j)%interpPoints(node,1)
              yi = owi(j)%interpPoints(node,2)
              
              w1=owi(j)%weight(node,1)
              w2=owi(j)%weight(node,2)
              w3=owi(j)%weight(node,3)
              w4=owi(j)%weight(node,4)

              u=w1*owi(j)%u(xi,yi)+w2*owi(j)%u(xi+1,yi)+w3*
     &           owi(j)%u(xi+1,yi+1)+w4*owi(j)%u(xi,yi+1)
              v=w1*owi(j)%v(xi,yi)+w2*owi(j)%v(xi+1,yi)+w3*
     &           owi(j)%v(xi+1,yi+1)+w4*owi(j)%v(xi,yi+1)
              p=w1*owi(j)%p(xi,yi)+w2*owi(j)%p(xi+1,yi)+w3*
     &           owi(j)%p(xi+1,yi+1)+w4*owi(j)%p(xi,yi+1)

              isSet = .true.
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
              call allMessage(DEBUG,"Return.")
#endif  
              call unsetMessageSource()
              return
           endif
        enddo

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.")
#endif  
        call unsetMessageSource()
!-----------------------------------------------------------------------
       end subroutine owi_interpolateToMesh
!-----------------------------------------------------------------------
          

!-----------------------------------------------------------------------
!> @brief  Allocates the u,v,p fields in an oceanweather object and
!>         deallocates if necessary
!> @param[inout] this array of oceanweather objects
!-----------------------------------------------------------------------
      subroutine owi_allocate(this)
!-----------------------------------------------------------------------
            IMPLICIT NONE
            TYPE(OCEANWEATHER),INTENT(INOUT) :: this
            call setMessageSource("owi_allocate")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG,"Enter.")
#endif
            IF((SIZE(this%p,1).NE.this%iLon).OR.(SIZE(this%p,2).NE.this%iLat))THEN
                CALL owi_deallocate(this)
                ALLOCATE(this%p(this%ilon,this%ilat))
                ALLOCATE(this%u(this%ilon,this%ilat))
                ALLOCATE(this%v(this%ilon,this%ilat))
            ENDIF
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG,"Return.")
#endif
            call unsetMessageSource()
!-----------------------------------------------------------------------
      end subroutine owi_allocate
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Deallocates the u,v,p fields in an oceanweather object
!> @param[inout] this array of oceanweather objects
!-----------------------------------------------------------------------
      subroutine owi_deallocate(this)
!-----------------------------------------------------------------------
            IMPLICIT NONE
            TYPE(OCEANWEATHER),INTENT(INOUT) :: this
            call setMessageSource("owi_deallocate")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG,"Enter.")
#endif
            IF(ALLOCATED(this%p))DEALLOCATE(this%p)
            IF(ALLOCATED(this%u))DEALLOCATE(this%u)
            IF(ALLOCATED(this%v))DEALLOCATE(this%v)
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG,"Return.")
#endif
            call unsetMessageSource()
            RETURN
!-----------------------------------------------------------------------
      end subroutine owi_deallocate
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Generates coordinates for a specified oceanweather domain
!> @param[in] this array of oceanweather objects
!> @param[out] lat latitude coordinates of oceanweather grid
!> @param[out] lon longiude coordinates of oceanweather grid
!-----------------------------------------------------------------------
      subroutine owi_generateCoordinates(this,lat,lon)
!-----------------------------------------------------------------------
            IMPLICIT NONE
            TYPE(OCEANWEATHER),INTENT(IN)      :: this
            REAL(SZ),ALLOCATABLE,INTENT(OUT)   :: lat(:)
            REAL(SZ),ALLOCATABLE,INTENT(OUT)   :: lon(:)
            INTEGER                            :: I
      
            call setMessageSource("owi_generateCoordinates")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG,"Enter.")
#endif
            ALLOCATE(lat(1:this%ilat))
            ALLOCATE(lon(1:this%ilon))

            DO I = 1,this%ilat
                lat(I) = this%swlat+(I-1)*this%dLat
            ENDDO
            DO I = 1,this%ilon
                lon(i) = this%swlon+(I-1)*this%dLon
            ENDDO

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG,"Return.")
#endif
            call unsetMessageSource()

            RETURN
!-----------------------------------------------------------------------
      end subroutine owi_generateCoordinates
!-----------------------------------------------------------------------
        

!-----------------------------------------------------------------------
!> @brief Sets the filenames that will be used for the oceanweather wind 
!>        fields
!> @param[inout] this array of oceanweather objects
!> @param[in] iounit fortran unit number for fort.22
!> @param[inout] nWindFields number of wind fields specified in fort.22
!-----------------------------------------------------------------------
      subroutine owi_setFilenames(this,iounit,nWindFields)
!-----------------------------------------------------------------------
            implicit none
            integer,intent(inout)            :: nWindFields
            integer,intent(in)               :: iounit
            type(oceanweather),intent(inout) :: this(*)

            call setMessageSource("owi_setFilenames")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG,"Enter.")
#endif

            if(nWindFields.LT.0)then
                nWindFields = abs(nWindFields)
                call owi_setDynamicFilenames(this,iounit,nWindFields)
            else
                if(nWindFields.GT.3)then
                  CALL check_owi_err(1,"owi file list must be used for numSets>3. "//
     &                                 "Use negative numSets to read a file list.","fort.22",22)
                 endif
                 call owi_setDefaultFilenames(this)
            endif
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG,"Return.")
#endif
            call unsetMessageSource()
!-----------------------------------------------------------------------
      end subroutine owi_setFilenames
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Sets the dynamic filenames read from the fort.22 file
!> @param[inout] this array of oceanweather objects
!> @param[in] iounit fortran unit number for fort.22
!> @param[in] nWindFields number of wind fields specified in fort.22
!-----------------------------------------------------------------------
      subroutine owi_setDynamicFilenames(this,iounit,nWindFields)
!-----------------------------------------------------------------------
            use sizes,only:gblinputdir
            implicit none
            integer,intent(in)               :: iounit
            integer,intent(in)               :: nWindFields
            type(oceanweather),intent(inout) :: this(*)
            character(1024)                  :: f1,f2
            integer                          :: i

            call setMessageSource("owi_setDynamicFilenames")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG,"Enter.")
#endif
            do i = 1,nWindFields
              read(iounit,*) f1, f2
              owi(i)%pressure_file = TRIM(GBLINPUTDIR)//"/"//TRIM(f1)
              owi(i)%wind_file = TRIM(GBLINPUTDIR)//"/"//TRIM(f2)
              owi(i)%iounit_pressure = 500+i
              owi(i)%iounit_wind = 600+i
            enddo

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG,"Return.")
#endif
            call unsetMessageSource()
!-----------------------------------------------------------------------
      end subroutine owi_setDynamicFilenames
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Skips the specified number of snaps in the oceanweather wind files
!> @param[in] numSnaps number of snaps to skip
!> @param[in] np number of nodes in adcirc mesh
!> @param[inout] wx wind velocity (x)
!> @param[inout] wy wind velocity (y)
!> @param[inout] pr atmospheric pressure
!> @param[in] prbk background atmospheric pressure
!-----------------------------------------------------------------------
      subroutine owi_skipSnaps(numSnaps,np,wx,wy,pr,prbk)
!-----------------------------------------------------------------------
             implicit none
             integer,intent(in)     :: numSnaps,np
             real(sz),intent(inout) :: wx(:),wy(:),pr(:)
             real(sz),intent(in)    :: prbk
             integer                :: i
            call setMessageSource("owi_skipSnaps")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG,"Enter.")
#endif
             do i = 1,numSkipSnaps
                write(scratchMessage,41) i
  41            format("Skipping snap '",I6,"' in OWI wind data.")
                call logMessage(DEBUG,trim(scratchMessage))
                CALL NWS12GET(wx,wy,pr,np,prbk)
             enddo
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG,"Return.")
#endif
            call unsetMessageSource()
!-----------------------------------------------------------------------
      end subroutine owi_skipSnaps
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Sets the default names for the oceanweather files if numSets is
!>        positive. fort.221, fort.222, fort.223, fort.224, fort.217, 
!>        fort.218
!> @param[inout] this array of oceanweather objects
!-----------------------------------------------------------------------
      subroutine owi_setDefaultFilenames(this)
!-----------------------------------------------------------------------
            use sizes,only:gblinputdir
            implicit none
            type(oceanweather),intent(inout) :: this(*)
            call setMessageSource("owi_setDefaultFilenames")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG,"Enter.")
#endif
            this(1)%pressure_file = TRIM(GBLINPUTDIR)//"/fort.221"
            this(1)%wind_file = TRIM(GBLINPUTDIR)//"/fort.222"
            this(1)%iounit_pressure = 221
            this(1)%iounit_wind = 222
            IF(numSets.GT.1)THEN
              this(2)%pressure_file = TRIM(GBLINPUTDIR)//"/fort.223"
              this(2)%wind_file = TRIM(GBLINPUTDIR)//"/fort.224"
              this(2)%iounit_pressure = 223
              this(2)%iounit_wind = 224
            ENDIF
            IF(numSets.GT.2)THEN
              this(3)%pressure_file = TRIM(GBLINPUTDIR)//"/fort.217"
              this(3)%wind_file =  TRIM(GBLINPUTDIR)//"/fort.218"
              this(3)%iounit_pressure = 217
              this(3)%iounit_wind = 218
            ENDIF
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG,"Return.")
#endif
            call unsetMessageSource()
!-----------------------------------------------------------------------
      end subroutine owi_setDefaultFilenames
!-----------------------------------------------------------------------
        

!-----------------------------------------------------------------------
!> @brief Initializes the data in the oceanweather data objects based on
!>        the file headers
!> @param[inout] this array of oceanweather objects
!-----------------------------------------------------------------------
      subroutine owi_initializeFileHeaders(this)
!-----------------------------------------------------------------------
            implicit none
            type(oceanweather),intent(inout) :: this(*)
            integer                          :: i
            integer(8)                       :: date1pressure,date2pressure
            integer(8)                       :: date1wind,date2wind
            character(1024)                  :: errorVar

            call setMessageSource("owi_initializeFileHeaders")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG,"Enter.")
#endif
            DO I = 1,numSets
              CALL readHeader(this(i)%pressure_file,this(i)%iounit_pressure,
     &                        date1pressure,date2pressure)
              CALL readHeader(this(i)%wind_file,this(i)%iounit_wind,
     &                        date1wind,date2wind)

              IF(date1pressure.NE.date1wind.OR.date2pressure.NE.date2wind)THEN
                  CALL allMessage(ERROR,"Wind and pressure dates do not match.")
                  errorVar = ""
                  CALL check_owi_err(1,errorVar,this(i)%wind_file,this(i)%iounit_wind)
              ENDIF

              this(i)%startDate = date1wind
              this(i)%endDate   = date2wind

              IF(this(i)%startDate.NE.this(1)%startdate)THEN
                  CALL allMessage(ERROR,"Wind domains have different time frames")
                  errorVar = ""
                  CALL check_owi_err(1,errorVar,this(i)%wind_file,this(i)%iounit_wind)
              ENDIF

              this(i)%isnap = 0
              this(i)%iupdate = 1
            ENDDO

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG,"Return.")
#endif
            call unsetMessageSource()
!-----------------------------------------------------------------------
      end subroutine owi_initializefileHeaders
!-----------------------------------------------------------------------
     
!-----------------------------------------------------------------------
!> @brief Allocates the arrays used to interpolate oceanweather grids to
!>        the mesh
!> @param[inout] this array of oceanweather objects
!> @param[in] numNodes number of nodes in the adcirc mesh
!-----------------------------------------------------------------------
      subroutine owi_allocateInterpolationArrays(this,numNodes)
!-----------------------------------------------------------------------
           implicit none
           type(oceanweather),intent(inout) :: this(*)
           integer,intent(in)               :: numNodes
           integer                          :: i
           
           call setMessageSource("owi_initializeFileHeaders")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
           call allMessage(DEBUG,"Enter.")
#endif
           do i=1,numSets
              ALLOCATE(this(i)%interpPoints(numNodes,2))
              ALLOCATE(this(i)%weight(numNodes,4))
           enddo

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG,"Return.")
#endif
            call unsetMessageSource()
!-----------------------------------------------------------------------
      end subroutine owi_allocateInterpolationArrays
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief This generates and saves interpolation coefficients for mapping
!>        from a OWI grid to an ADCIRC mesh.
!> @param[in] idx index in the owi object array
!> @param[in] np number of points in the adcirc mesh
!-----------------------------------------------------------------------
      SUBROUTINE NWS12INTERP(idx, NP)
!-----------------------------------------------------------------------
      USE GLOBAL, ONLY : RAD2DEG
      USE MESH, ONLY : SLAM, SFEA

      IMPLICIT NONE

      INTEGER,INTENT(IN)   :: idx
      INTEGER,INTENT(IN)   :: NP
      INTEGER              :: I,J,K,XI,YI
      REAL(SZ)             :: adcLat,adcLong
      REAL(SZ)             :: w,w1,w2,w3,w4
      REAL(SZ),ALLOCATABLE :: lat(:), lon(:)

      call setMessageSource("nws12interp")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      WRITE(16,'(A)') ''
      WRITE(16,'(A,I0)') 'WIND MAPPING UPDATED FOR REGION ',IDX
      WRITE(16,'(A)') ''

      CALL owi_allocate(owi(idx))
      CALL owi_generateCoordinates(owi(idx),lat,lon)

      ! Generate interpolation coefficients (south west point and weights)
      do i=1,NP
        adcLat  = RAD2DEG*SFEA(i)
        adcLong = RAD2DEG*SLAM(i)

        if (adcLong>=lon(1).and.adcLong<lon(owi(idx)%iLon) .and.
     &       adcLat>=lat(1).and.adcLat<lat(owi(idx)%iLat)) then
          do j=1,owi(idx)%iLon-1
            if (adcLong>=lon(j) .and.
     &           adcLong<lon(j+1)) then
              xi=j
              goto 200
            endif
          enddo
 200      continue

          do k=1,owi(idx)%iLat-1
            if (adcLat>=lat(k) .and.
     &           adcLat<lat(k+1)) then
              yi=k
              goto 300
            endif
          enddo
 300      continue

          owi(idx)%interpPoints(i,1) = xi
          owi(idx)%interpPoints(i,2) = yi

          w=(lon(xi+1)-lon(xi))*(lat(yi+1)-lat(yi))
          w1=(lon(xi+1)-adcLong)*(lat(yi+1)-adcLat)
          w2=(adcLong-lon(xi))*(lat(yi+1)-adcLat)
          w3=(adcLong-lon(xi))*(adcLat-lat(yi))
          w4=(lon(xi+1)-adcLong)*(adcLat-lat(yi))

          owi(idx)%weight(i,1)=w1/w
          owi(idx)%weight(i,2)=w2/w
          owi(idx)%weight(i,3)=w3/w
          owi(idx)%weight(i,4)=w4/w

        else
          owi(idx)%interpPoints(i,1) = 0
          owi(idx)%interpPoints(i,2) = 0
        endif
      enddo

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
!-----------------------------------------------------------------------
      END SUBROUTINE NWS12INTERP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Read the header from the OWI file and returns the values of
!>        start and end dates
!> @param[in] filename name of the file to open and read
!> @param[in] lun unit number to open the file as
!> @param[out] date1 start date in the header
!> @param[out] date2 end date in the header
!-----------------------------------------------------------------------
      SUBROUTINE readHeader(filename, lun, date1, date2)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER,INTENT(IN)       :: lun
      CHARACTER(*),INTENT(IN)  :: filename
      INTEGER(8),INTENT(OUT)   :: date1
      INTEGER(8),INTENT(OUT)   :: date2
      CHARACTER(80)            :: owiheader
      CHARACTER(1024)          :: errorVar
      INTEGER                  :: errorIO

      call setMessageSource("readHeader")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      call openFileForRead(lun, filename, errorIO)
      errorVar = "opening Oceanweather file"
      call check_owi_err(errorIO,errorVar,filename,lun)

      ! Read begining/ending dates of pre file
      owiheader(:) = ' '  !set owiheader to blanks before read
      errorVar = "reading Oceanweather header line"
      read(lun, fmt='(a80)',
     &   end=99998,err=99999,iostat=errorIO) owiheader
      call check_owi_err(errorIO,errorVar,filename,lun)

      errorVar = "reading Oceanweather start date"
      read(owiheader(56:65),'(I10)',
     &   end=99998,err=99999,iostat=errorIO) date1
      call check_owi_err(errorIO,errorVar,filename,lun)
      write(scratchMessage,31) TRIM(errorVar), trim(filename), date1
  31  format("'",A,"' in  '",A,"' is '",I10,"'.")
      call allMessage(INFO,scratchMessage)

      errorVar = "reading Oceanweather end date"
      read(owiheader(71:80),'(I10)',
     &    end=99998,err=99999,iostat=errorIO) date2
      call check_owi_err(errorIO,errorVar,filename,lun)
      write(scratchMessage,31) TRIM(errorVar), trim(filename), date2
      call allMessage(INFO,scratchMessage)

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN

99998 call allMessage(ERROR,"Unexpectedly reached end-of-file.")    ! END jumps here
99999 call check_owi_err(errorIO,"Error reading Oceanweather file",filename,lun) !  ERR jumps here
!-----------------------------------------------------------------------
      end subroutine readHeader
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! @brief Checks the return value from subroutine calls; if there
!        was an error, it writes a termination message to the screen 
!        and to the fort.16 file and terminates ADCIRC.
! @param[in] iret return value to check. nonzero indicates error
! @param[in] errorVar process being executed when the error occured
! @param[in] filename file being processed when error occured
! @param[in] unitnumber unit number of the file being processed during error
!-----------------------------------------------------------------------
      subroutine check_owi_err(iret, errorVar, filename, unitnumber)
!-----------------------------------------------------------------------
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      IMPLICIT NONE
      INTEGER,INTENT(IN)               :: iret
      CHARACTER(*),INTENT(IN)          :: errorVar
      CHARACTER(*),INTENT(IN)          :: filename
      INTEGER,INTENT(IN),OPTIONAL      :: unitnumber

      call setMessageSource("check_owi_err")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      if (iret.ne.0) then
         if (errorVar.ne."") then
            write(scratchMessage,888)
     &         trim(errorVar), trim(filename), unitnumber
 888        format("Failed to '",A,"' from '",A,
     &         "' (unit number ",I3,").")
            call allMessage(ERROR,trim(scratchMessage))
            write(6,*) 'iostat =',iret
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
!-----------------------------------------------------------------------
      end subroutine check_owi_err
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      END MODULE OWIWIND
!-----------------------------------------------------------------------

