C******************************************************************************
C  MODULE OWIWIND
C    Written by s.b. 08/17/2006
C    Modified by TCM 08/08/2019 -- Added LandFall Scale win/pre (finest level)
C                                   Basin ->   Region ->  LandFall
C                      File Units (221,222)   (223,224)   (217,218)
C                      Added PRdeflt to nws12init and nws12get calls
C                      in order to use background pressure value instead of
C                      hard wired 1013.0 mb.
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

      integer  :: nPowellSearchDomains = -1
      integer  :: numSets,numBlankSnaps,cntSnaps,numSkipSnaps
      real(SZ) :: windMultiplier
      logical :: moving_grid = .false. ! true will spatially interpolate each dataset
      
      TYPE OCEANWEATHER
          CHARACTER(1024)      :: pressure_file
          CHARACTER(1024)      :: wind_file
          INTEGER(8)           :: startDate
          INTEGER(8)           :: endDate
          INTEGER(8)           :: date
          INTEGER              :: minute
          INTEGER              :: iLat
          INTEGER              :: iLon
          INTEGER              :: iounit_pressure
          INTEGER              :: iounit_wind
          INTEGER              :: iSnap
          INTEGER              :: iUpdate
          REAL(SZ)             :: dLat
          REAL(SZ)             :: dLon
          REAL(SZ)             :: swLon
          REAL(SZ)             :: swLat
          REAL(SZ),ALLOCATABLE :: p(:,:)
          REAL(SZ),ALLOCATABLE :: u(:,:)
          REAL(SZ),ALLOCATABLE :: v(:,:)
          REAL(SZ),ALLOCATABLE :: weight(:,:)
          INTEGER,ALLOCATABLE  :: interpPoints(:,:)
      END TYPE OCEANWEATHER

      TYPE(OCEANWEATHER),ALLOCATABLE :: owi(:) 

      PRIVATE OCEANWEATHER,owi,owi_allocate,owi_deallocate
      PRIVATE owi_generateCoordinates

C---------------------end of data declarations--------------------------------C

      CONTAINS
          
        SUBROUTINE owi_allocate(this)
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
        END SUBROUTINE owi_allocate

        SUBROUTINE owi_deallocate(this)
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
        END SUBROUTINE owi_deallocate

        SUBROUTINE owi_generateCoordinates(this,lat,lon)
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
        END SUBROUTINE owi_generateCoordinates
     
C***********************************************************************
C   SUBROUTINE NWS12INIT
C***********************************************************************

      SUBROUTINE NWS12INIT(WVNX,WVNY,PRN,NP,RHOWAT0,G,PRdeflt)

      USE SIZES, ONLY : SZ,MYPROC, GBLINPUTDIR
      USE GLOBAL,ONLY: PRBCKGRND

      IMPLICIT NONE
      INTEGER               :: NP,I
      INTEGER(8)            :: date1pressure,date2pressure
      INTEGER(8)            :: date1wind,date2wind
      REAL(SZ), intent(out) :: WVNX(:),WVNY(:),PRN(:)
      REAL(SZ)              :: RHOWAT0,RHOWATG,G,PRdeflt
      character(len=1024)   :: tmpstr
      character(len=1024)   :: errorVar = ""
      character(len=1024)   :: fname,f1,f2
      integer               :: iounit
      integer               :: errorIO = 0

      call setMessageSource("nws12init")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

! Read meta info -------------------------------------------------
C
C     R E A D   F O R T 2 2
C     ---------------------
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

      IF(nPowellSearchDomains.GT.numSets)THEN
          nPowellSearchDomains = numSets
      ELSEIF(nPowellSearchDomains.LE.0)THEN
          nPowellSearchDomains = numSets
      ENDIF

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

      allocate(owi(1:abs(numSets)))
      if(numSets.LT.0)then
          numSets = abs(numSets)
          do i = 1,numSets
            read(iounit,*) f1, f2
            owi(i)%pressure_file = TRIM(GBLINPUTDIR)//"/"//TRIM(f1)
            owi(i)%wind_file = TRIM(GBLINPUTDIR)//"/"//TRIM(f2)
            owi(i)%iounit_pressure = 500+i
            owi(i)%iounit_wind = 600+i
          enddo
      else
        IF(numSets.GT.3)THEN
            CALL check_owi_err(1,"owi file list must be used for numSets>3. "//
     &                           "Use negative numSets to read a file list.",fname,iounit)
        ENDIF
        owi(1)%pressure_file = TRIM(GBLINPUTDIR)//"/fort.221"
        owi(1)%wind_file = TRIM(GBLINPUTDIR)//"/fort.222"
        owi(1)%iounit_pressure = 221
        owi(1)%iounit_wind = 222
        IF(numSets.GT.1)THEN
          owi(2)%pressure_file = TRIM(GBLINPUTDIR)//"/fort.223"
          owi(2)%wind_file = TRIM(GBLINPUTDIR)//"/fort.224"
          owi(2)%iounit_pressure = 223
          owi(2)%iounit_wind = 224
        ENDIF
        IF(numSets.GT.2)THEN
          owi(3)%pressure_file = TRIM(GBLINPUTDIR)//"/fort.217"
          owi(3)%wind_file =  TRIM(GBLINPUTDIR)//"/fort.218"
          owi(3)%iounit_pressure = 217
          owi(3)%iounit_wind = 218
        ENDIF
      endif
            
      close(iounit)

      DO I = 1,numSets
        CALL readHeader(owi(i)%pressure_file,owi(i)%iounit_pressure,
     &                  date1pressure,date2pressure)
        CALL readHeader(owi(i)%wind_file,owi(i)%iounit_wind,
     &                  date1wind,date2wind)
        IF(date1pressure.NE.date1wind.OR.date2pressure.NE.date2wind)THEN
            CALL allMessage(ERROR,"Wind and pressure dates do not match.")
            errorVar = ""
            CALL check_owi_err(1,errorVar,owi(i)%wind_file,owi(i)%iounit_wind)
        ENDIF
        owi(i)%startDate = date1wind
        owi(i)%endDate   = date2wind
        IF(owi(i)%startDate.NE.owi(1)%startdate)THEN
            CALL allMessage(ERROR,"Wind domains have different time frames")
            errorVar = ""
            CALL check_owi_err(1,errorVar,owi(i)%wind_file,owi(i)%iounit_wind)
        ENDIF
        ALLOCATE(owi(i)%interpPoints(np,2))
        ALLOCATE(owi(i)%weight(np,4))
        owi(i)%isnap = 0
        owi(i)%iupdate = 1
      ENDDO

      cntSnaps = 0

! Skip snaps if necessary ---------------------------------------------------
      do i = 1,numSkipSnaps
         write(scratchMessage,41) i
  41     format("Skipping snap '",I6,"' in OWI wind data.")
         call logMessage(DEBUG,trim(scratchMessage))
         CALL NWS12GET(WVNX,WVNY,PRN,NP,RHOWAT0,G,PRBCKGRND)
      enddo

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN

99998 call allMessage(ERROR,"Unexpectedly reached end-of-file.") ! END jumps here

99999 call check_owi_err(errorIO,errorVar,fname,iounit) !  ERR jumps here
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
      USE GLOBAL,ONLY:PRBCKGRND
      implicit none
      real(sz), intent(out) :: eyeLatR(3)
      real(sz), intent(out) :: eyeLonR(3)
      logical, intent(out)  :: foundEye
      !Casey 110518: Added for Mark Powell's sector-based wind drag.
      INTEGER  :: EyeLatI
      INTEGER  :: EyeLonI
      REAL(SZ) :: EyeLatTemp
      REAL(SZ) :: EyeLonTemp
      REAL(SZ) :: EyePressure = PRBCKGRND
      integer :: i, j, s

      DO S = 1,nPowellSearchDomains
          DO i=1,owi(i)%iLon
             DO j=1,owi(i)%iLat
                IF((owi(i)%p(i,j).LT.EyePressure).AND.
     &             (owi(i)%p(i,j).LT.1000.D0))THEN
                   EyeLatI = j
                   EyeLonI = i
                   EyePressure = owi(i)%p(i,j)
                ENDIF
             ENDDO
          ENDDO
          IF((EyeLatI.EQ.0).AND.(EyeLonI.EQ.0))THEN
             FoundEye = .FALSE.
          ELSEIF((EyeLatI.EQ.1).OR.(EyeLatI.EQ.owi(i)%iLat).OR.
     &           (EyeLonI.EQ.1).OR.(EyeLonI.EQ.owi(i)%iLon))THEN
             FoundEye = .FALSE.
          ELSE
             FoundEye = .TRUE.
             EyeLatTemp = owi(i)%swlat  + (EyeLatI-1)*owi(i)%dLat
             EyeLonTemp = owi(i)%swlon + (EyeLonI-1)*owi(i)%dLon
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
      ENDDO

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
      INTEGER                ::NP,I,J,S,XI,YI
      REAL(SZ), intent(out), dimension(:) ::  WVNX,WVNY,PRN
      REAL(SZ), intent(in)   :: G, PRdeflt
      REAL(SZ)               :: RHOWAT0,RHOWATG,EyePressure
      INTEGER                :: EyeLatI,EyeLonI
      INTEGER                :: EyeLatTemp,EyeLonTemp
      CHARACTER(1024)        :: filename,errorVar = ""
      INTEGER                :: iounit
      INTEGER                :: iLatw,iLongw,iLatp,iLongp
      REAL(SZ)               :: dxw,dyw,swlatw,swlongw
      REAL(SZ)               :: dxp,dyp,swlatp,swlongp
      REAL(SZ)               :: uu,vv,pp
      REAL(SZ)               :: w1,w2,w3,w4
      INTEGER(8)             :: iCYMDHw,iCYMDHp
      INTEGER                :: iMinw,iMinp
      INTEGER                :: errorIO

      !jgf51.52.05: Used to detect mesh nodes outside both met grids
      ! and set those nodes to zero wind velocity and to background pressure.
      LOGICAL :: uvpHasBeenSet,FoundEye

      call setMessageSource("nws12get")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      RHOWATG=RHOWAT0*G

      ! Increment counter (cntSnaps initialized to zero in nws12init)
      cntSnaps = cntSnaps+1

      ! Put a blank snap for the first 'numBlankSnaps' snaps and then return
      if (cntSnaps.LE.numBlankSnaps) then
         WVNX(:) = 0D0
         WVNY(:) = 0D0
         PRN(:) = Prdeflt*100D0/RHOWATG

         write(scratchMessage,16) cntSnaps
 16      format('INSERTING A BLANK WIND SNAP, COUNT=',i4)
         call allMessage(INFO,trim(scratchMessage))
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         RETURN
      endif

      DO S = 1,numSets

        if(s.gt.1)then
            if(owi(1)%date.lt.owi(s)%startDate) cycle
            if(owi(1)%date.eq.owi(s)%endDate.and.owi(s)%minute.ne.0) cycle
            if(owi(1)%date.gt.owi(s)%endDate) cycle
        endif

        ! Increment counter
        owi(s)%isnap = owi(s)%isnap + 1

        ! Read grid specifications/date in pressure file
        filename = owi(s)%pressure_file
        iounit = owi(s)%iounit_pressure
        WRITE(errorVar,'(A,I0)') "read grid specifications/date in Oceanweather pressure file, domain ",s
        read (owi(s)%iounit_pressure,11,end=10000,err=9999,iostat=errorIO)
     &        iLatp,iLongp,dxp,dyp,swlatp,swlongp,iCYMDHp,iMinp
        call check_owi_err(errorIO,errorVar,filename,iounit)

        ! Read grid specifications/date in wind file
        filename = owi(s)%wind_file
        iounit = owi(s)%iounit_wind
        WRITE(errorVar,'(A,I0)') "read grid specifications/date in Oceanweather wind file, domain ",s
        read (owi(s)%iounit_wind,11,end=10000,err=9999,iostat=errorIO)
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

        ! Check if header info has changed from the previous snapshot
        if(moving_grid)then
            ! for a moving grid,update it every time  Jie 09/2013
            owi(s)%iupdate = 1
        else
            if(owi(s)%isnap.gt.1) then
              if(iLatp.ne.owi(s)%iLat.or.iLongp.ne.owi(s)%iLon.or.
     &             dxp.ne.owi(s)%dlon.or.dyp.ne.owi(s)%dlat.or.
     &             swlatp.ne.owi(s)%swlat.or.swlongp.ne.owi(s)%swlon) then
                call logMessage(INFO,"Basin scale grid has changed.")
                owi(s)%iupdate = 1
              else
                owi(s)%iupdate = 0
              endif
            endif
        endif

        owi(s)%date = iCYMDHp
        owi(s)%minute = iMinp

        ! Update coordinate mapping coefficients if necessary
        if(owi(s)%iupdate.eq.1) then
           owi(s)%ilon = ilongp
           owi(s)%ilat = ilatp
           owi(s)%dlon = dxp
           owi(s)%dlat = dyp
           owi(s)%swlat = swlatp
           owi(s)%swlon = swlongp
           WRITE(errorVar,'(A,I0)') "Updating grid coordinate mapping coefficients, domain ",s
           call logMessage(INFO,errorVar)
           call nws12interp(s,np)
        endif

        ! Read basin scale atmospheric pressure snapshot
        filename = owi(s)%pressure_file
        iounit =  owi(s)%iounit_pressure
        WRITE(errorVar,'(A,I0)') "reading atmospheric pressure snapshot, domain ",s
        read(owi(s)%iounit_pressure,22,end=10000,err=9999,iostat=errorIO)
     &     ((owi(s)%p(i,j),i=1,owi(s)%iLon),j=1,owi(s)%iLat)
        call check_owi_err(errorIO,errorVar,filename,iounit)

        ! Read basin scale snapshot of u/v components of the wind
        filename = owi(s)%wind_file
        iounit = owi(s)%iounit_wind
        WRITE(errorVar,'(A,I0)') "reading wind u-velocity snapshot, domain ",s
        read(owi(s)%iounit_wind,22,end=10000,err=9999,iostat=errorIO)
     &     ((owi(s)%u(i,j),i=1,owi(s)%iLon),j=1,owi(s)%iLat)
        call check_owi_err(errorIO,errorVar,filename,iounit)
        filename = owi(s)%wind_file
        iounit = owi(s)%iounit_wind
        WRITE(errorVar,'(A,I0)') "reading wind v-velocity snapshot, domain ",s
        read(owi(s)%iounit_wind,22,end=10000,err=9999,iostat=errorIO)
     &     ((owi(s)%v(i,j),i=1,owi(s)%iLon),j=1,owi(s)%iLat)
        call check_owi_err(errorIO,errorVar,filename,iounit)
      ENDDO

      do i=1,NP

         ! jgf51.52.05: Fixed an issue with floating point .eq. that
         ! was used to determine whether the value of uu had been set. 
         ! Replaced floating point comparison with a logical variable.
         uvpHasBeenSet = .false.

         do j = numSets,1,-1
            if (owi(j)%interppoints(i,1).gt.0) then
               xi = owi(j)%interpPoints(i,1)
               yi = owi(j)%interpPoints(i,2)
               
               w1=owi(j)%weight(i,1)
               w2=owi(j)%weight(i,2)
               w3=owi(j)%weight(i,3)
               w4=owi(j)%weight(i,4)

               uu=w1*owi(j)%u(xi,yi)+w2*owi(j)%u(xi+1,yi)+w3*
     &            owi(j)%u(xi+1,yi+1)+w4*owi(j)%u(xi,yi+1)
               vv=w1*owi(j)%v(xi,yi)+w2*owi(j)%v(xi+1,yi)+w3*
     &            owi(j)%v(xi+1,yi+1)+w4*owi(j)%v(xi,yi+1)
               pp=w1*owi(j)%p(xi,yi)+w2*owi(j)%p(xi+1,yi)+w3*
     &            owi(j)%p(xi+1,yi+1)+w4*owi(j)%p(xi,yi+1)

               uvpHasBeenSet = .true.
               exit
            endif
         enddo

         if (.not.uvpHasBeenSet) then
            WVNX(I)=0.d0
            WVNY(I)=0.d0
            PRN(I)=PRdeflt*100.d0/RHOWATG
         else
            !Convert millibars to m of water
            PRN(i)  = 100.d0*PP/RHOWATG

            ! Apply wind velocity multiplier
            ! and save to arrays
            WVNX(i) = uu * windMultiplier
            WVNY(i) = vv * windMultiplier

         endif

      enddo

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN


9999  call check_owi_err(errorIO,errorVar,filename,iounit) ! ERR during read jumps to here

      ! END during read of required data jumps to here
9997  call allMessage(ERROR,"Unexpected end-of-file reached.")
      call check_owi_err(errorIO,errorVar,filename,iounit)

10000 continue  ! END during read of u, v, p data jumps to here
      write(scratchMessage,61) trim(errorVar)
  61  format("Unexpected end-of-file while reading '",A,
     &    "'. Wind speeds set to zero and pressure to 1013mb.")
      call allMessage(WARNING,trim(scratchMessage))

      WVNX(:)=0.d0
      WVNY(:)=0.d0
      PRN(:)=101300.d0/RHOWATG

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
C   SOBROUTINE NWS12INTERP
C
C   This generates and saves interpolation coefficients for mapping
C   from a OWI grid to an ADCIRC mesh.
C
C***********************************************************************

      SUBROUTINE NWS12INTERP(idx, NP)
      USE GLOBAL, ONLY : RAD2DEG
      USE MESH, ONLY : SLAM, SFEA

      IMPLICIT NONE

      INTEGER,INTENT(IN)   :: idx
      INTEGER,INTENT(IN)   :: NP
      INTEGER              :: I,J,K,XI,YI
      REAL(SZ)             :: adcLat,adcLong
      REAL(SZ)             :: w,w1,w2,w3,w4
      REAL(SZ),ALLOCATABLE :: lat(:), lon(:)
C
      call setMessageSource("nws12interp")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
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
      END SUBROUTINE NWS12INTERP

C-----------------------------------------------------------------------
C     S U B R O U T I N E   R E A D   H E A D E R
C-----------------------------------------------------------------------
C     jgf50.16 Read the header from the OWI file and sets the value
C     of the date1p, date2p, date1w, or date2w variables as appropriate.
C-----------------------------------------------------------------------
      SUBROUTINE readHeader(filename, lun, date1, date2)
      IMPLICIT NONE
      INTEGER,INTENT(IN)       :: lun
      CHARACTER(*),INTENT(IN)  :: filename
      INTEGER(8),INTENT(OUT)   :: date1
      INTEGER(8),INTENT(OUT)   :: date2
      CHARACTER(80)            :: owiheader
      CHARACTER(1024)          :: errorVar
      INTEGER                  :: errorIO
C
      call setMessageSource("readHeader")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      call openFileForRead(lun, filename, errorIO)
      errorVar = "opening Oceanweather file"
      call check_owi_err(errorIO,errorVar,filename,lun)
C
      ! Read begining/ending dates of pre file
      owiheader(:) = ' '  !set owiheader to blanks before read
      errorVar = "reading Oceanweather header line"
      read(lun, fmt='(a80)',
     &   end=99998,err=99999,iostat=errorIO) owiheader
      call check_owi_err(errorIO,errorVar,filename,lun)
C
      errorVar = "reading Oceanweather start date"
      read(owiheader(56:65),'(I10)',
     &   end=99998,err=99999,iostat=errorIO) date1
      call check_owi_err(errorIO,errorVar,filename,lun)
      write(scratchMessage,31) TRIM(errorVar), trim(filename), date1
  31  format("'",A,"' in  '",A,"' is '",I10,"'.")
      call allMessage(INFO,scratchMessage)
C
      errorVar = "reading Oceanweather end date"
      read(owiheader(71:80),'(I10)',
     &    end=99998,err=99999,iostat=errorIO) date2
      call check_owi_err(errorIO,errorVar,filename,lun)
      write(scratchMessage,31) TRIM(errorVar), trim(filename), date2
      call allMessage(INFO,scratchMessage)
C
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN

99998 call allMessage(ERROR,"Unexpectedly reached end-of-file.")    ! END jumps here
99999 call check_owi_err(errorIO,"Error reading Oceanweather file",filename,lun) !  ERR jumps here
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
      subroutine check_owi_err(iret, errorVar, filename, unitnumber)
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

