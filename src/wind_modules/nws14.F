      module mod_nws14
        use mod_datetime, only: t_datetime, t_timedelta
        use global, only: ERROR, ECHO, WARNING, logMessage,
     &     allMessage, scratchMessage, DEBUG, Found_InterpInset_Namelist, UseInterpInset,
     &     basedatetime, setMessageSource, unsetMessageSource
        use adc_constants, only: rad2deg, g, PRBCKGRND, rhoWat0
        implicit none

        integer :: ub, ubc ! size of weightsp and indp interpolant weights
        real(8),allocatable :: weightsp(:,:,:)
        real(8),allocatable :: inset_weight(:,:,:)
        integer,allocatable  :: indp(:,:,:)
        integer,allocatable :: inset_indp(:,:,:)
        character(len=200) :: Pfile, Wfile, Pinv, Winv, Cinv, Cfile
        character(len=200) :: Wfile1, Pinv1, Winv1
        INTEGER :: PfileNCID, WfileNCID, CfileNCID, Wfile1NCID
        character(len=200) :: Pvar, Uvar, Vvar, Cvar, Tvar, Lonvar, Latvar
        character(len=200) :: Tdim, Londim, Latdim, Tformat

        CHARACTER(LEN=80)  :: InsetLon, InsetLat, InsetLonDim, InsetLatDim
        CHARACTER(LEN=256) :: InsetPvar, InsetUVar, InsetVVar
        type(t_datetime) :: CurDT, refdate, stepdate
        logical :: grb2flag
        integer :: NT, NTC
        real(8) :: rhoWat0g, PRBCKGRND_MH2O

        logical  :: read_NWS14_NetCdf_using_core_0 = .true. ! read NWS14 (grib2,netcdf) from a computed core 0
        logical  :: change_14 ! true when the NWS = 14 data has been updated

        private

        public :: NWS14INIT, NWS14GET, CLOSE_NWS14_FILES, set_change_14, get_change_14, read_NWS14_NetCdf_using_core_0

        contains

        subroutine set_change_14(val)
          implicit none
          logical, intent(in) :: val
          change_14 = val
        end subroutine set_change_14

        logical function get_change_14()
          implicit none
          get_change_14 = change_14
        end function get_change_14

C----------------------------------------------------------------------
      SUBROUTINE NWS14INIT(NWS, WTIME1)
C----------------------------------------------------------------------
C.... WJP Adding in subroutines for GRB2 reading below (NWS = 14)
C.... or netcdf (depending on what is present)
C....
C.... CPB 10/2023: reorganized this subroutine to make it more readable
C.... as well as to improve input for TACC
C----------------------------------------------------------------------
      use mod_datetime, only: operator(+)
      implicit none
      INTEGER,INTENT(IN) :: NWS
      REAL(8),INTENT(IN) :: WTIME1

      ! Initialise the datetime
      CurDT = basedatetime + t_timedelta(minutes=floor(WTIME1/60d0))
      NT = 0
      NTC = 0
      rhoWat0g = rhoWat0 * g
      PRBCKGRND_MH2O = 100.0D0 * PRBCKGRND / rhoWat0g

      ! check if we are using grib2 or netcdf input files:
      CALL NWS14CHECK_FILETYPE()

      ! set filenames appropriately
      CALL NWS14SET_FILENAMES()

      ! make grib2 .inv files (returns if using netcdf)
      CALL NWS14GRB2_MAKE_INV()

      ! read fort.22 to get netCDF variable info (returns if using
      ! grib2)
      CALL NWS14NC_READ_F22(NWS)

      ! Calculate the interpolant weights
      CALL NWS14_CALC_INTERP_WTS()

      call allMessage(ECHO,'Finished init of NWS14')

      IF(Found_InterpInset_Namelist) THEN
        CALL allMessage(ECHO,'Requesting to use insets. If files found,
     &                      insets will be interpolated!')
        CALL INIT_INSET_NC()
        CALL CALC_INSET_WTS()
      ENDIF
C----------------------------------------------------------------------
      END SUBROUTINE NWS14INIT
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE NWS14CHECK_FILETYPE()
C----------------------------------------------------------------------
C     CPB 10/2023: This subroutine simply checks whether we are using
C     grib2 or netCDF winds in NWS14 and sets the grb2flag logical
C     appropriately. It does this by:
C        1. Checking if the pressure grib2 file (fort.221.grb2) exists
C        2. If it does not exist or if we are compiled without grib2
C           then it checks if the pressure netCDF file (fort.221.nc)
C           exists.
C        3. If neither exist we abort with an error message.
C----------------------------------------------------------------------
      implicit none
      logical :: Fexists
      Pfile = 'fort.221.grb2'
      INQUIRE(file=Pfile,exist=Fexists)
#ifdef GRIB2API
      IF (Fexists) THEN

         grb2flag = .TRUE.
         RETURN
      ENDIF
#endif
      Pfile = 'fort.221.nc'
      INQUIRE(file=Pfile,exist=Fexists)
      IF (Fexists) THEN
         grb2flag = .FALSE.
         RETURN
      ELSE
         call allMessage(ERROR,
     &       'Neither .grb2 nor .nc wind files exist. Or, if .grb2 '//
     &       'files exist, did you compile with GRIB2 compiler flags?')
         call nws14Terminate()
      ENDIF

C----------------------------------------------------------------------
      END SUBROUTINE NWS14CHECK_FILETYPE
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE NWS14SET_FILENAMES()
C----------------------------------------------------------------------
C     CPB 10/2023: This subroutine sets the met forcing filenames based
C     on grb2flag as well as the existence or lack thereof of various
C     files. It also opens the netcdf files in read only mode to
C     eliminate opening/closing the files repeatedly
C----------------------------------------------------------------------
#ifdef ADCNETCDF
      use netcdf
      use netcdf_error, only: check_err
#endif
#ifdef CMPI
      USE MESSENGER, ONLY: MSG_BARRIER
#endif
      USE SIZES, ONLY : MYPROC
      USE GLOBAL, ONLY : NCICE
      implicit none
      logical :: Fexists
      ! grib2
      IF (grb2flag) THEN
         Pfile = 'fort.221.grb2'
         Pinv  = 'fort.221.inv'
         Wfile = 'fort.222.grb2'
         Wfile1 = 'fort.222.grb2'
         Winv  = 'fort.222.inv'
         Cfile = 'fort.225.grb2'
         Cinv  = 'fort.225.inv'

         Pvar = 'S'  ! searching for PRMSL or MSLET or PRES:SURFACE
                     ! prev: :PRMSL:mean sea level:'
         Uvar = ':UGRD:10 m above ground:'
         Vvar = ':VGRD:10 m above ground:'
         Cvar = ':ICEC:surface:'
         ! Assign dummy values to ncid variables
         PfileNCID = -9999
         WfileNCID = -9999
         Wfile1NCID = -9999
         IF (NCICE.EQ.14) THEN
            CfileNCID = -9999
         ENDIF
         RETURN
      ELSEIF (.NOT.grb2flag) THEN
         Pfile = 'fort.221.nc'
         Wfile = 'fort.222.nc'
         Wfile1 = 'fort.222.nc'
         Cfile = 'fort.225.nc'
         ! check if we are using two different wind files
         INQUIRE(FILE=Wfile,exist=Fexists)
         IF (.NOT.Fexists) THEN
            ! using 2 wind files
            Wfile = 'fort.222u.nc'
            Wfile1 = 'fort.222v.nc'
         ENDIF
         ! create dummy names for the inventory files
         Pinv = 'dmy'
         Winv = 'dmy'
         Cinv = 'dmy'
         ! create dummy names for the grib2 variables
         Pvar = 'dmy'
         Uvar = 'dmy'
         Vvar = 'dmy'
         Cvar = 'dmy'

#ifdef ADCNETCDF
         IF ( read_NWS14_NetCdf_using_core_0 ) THEN
           IF ( MYPROC .EQ. 0 ) THEN
             call Check_err(NF90_OPEN(Pfile,nf90_nowrite,PfileNCID))
             call Check_err(NF90_OPEN(Wfile,nf90_nowrite,WfileNCID))
             call Check_err(NF90_OPEN(Wfile1,nf90_nowrite,Wfile1NCID))
             IF (NCICE.EQ.14) THEN
                call Check_err(NF90_OPEN(Cfile,nf90_nowrite,CfileNCID))
             ENDIF
           END IF
#ifdef CMPI
         ! wait til proc 0 has opened all the files
           CALL MSG_BARRIER()
#endif
         END IF

#endif

         RETURN
      ELSE
         ! should be unreachable
         CALL nws14Terminate()
      ENDIF
C----------------------------------------------------------------------
      END SUBROUTINE NWS14SET_FILENAMES
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE NWS14GRB2_MAKE_INV()
C----------------------------------------------------------------------
C     CPB 10/2023: This subroutine makes the .inv files for grib2
C     meteorological forcing files. If we are not using grib2 then it
C     just returns
C----------------------------------------------------------------------
      USE SIZES, ONLY : MYPROC
#ifdef GRIB2API
      USE wgrib2api
#endif
      USE GLOBAL, ONLY : NCICE
#ifdef CMPI
      USE MESSENGER, ONLY: MSG_BARRIER
#endif
      USE GL2LOC_MAPPING ,ONLY : BcastToLocal_Int
      implicit none
      logical :: Fexists
      integer :: iret
      ! return if we are using netcdf
      IF (.NOT.grb2flag) THEN
         RETURN
      ENDIF
#ifdef GRIB2API
#ifdef CMPI
      ! only do this on proc 0
      IF (MYPROC.EQ.0) THEN
#endif
         ! pressure file (fort.221.grb2)
         inquire(file=Pinv,exist=Fexists)
         if (.not.Fexists) then
            iret = grb2_mk_inv(Pfile, Pinv)
            if (iret.ne.0) then
               call allMessage(ERROR,
     &                        'Fatal error in reading fort.221.grb2.')
            else
               call allMessage(ECHO,'successfully read fort.221.grb2 '
     &                            //'and wrote out fort.221.inv file')
            endif
         endif
         ! wind file (fort.222.grb2)
         inquire(file=Winv,exist=Fexists)
         if (.not.Fexists.and.iret.eq.0) then
            iret = grb2_mk_inv(Wfile, Winv)
            if (iret.ne.0) then
               call allMessage(ERROR,
     &                         'Fatal error in reading fort.222.grb2.')
            endif
            call allMessage(ECHO,'successfully read fort.222.grb2 '
     &                         //'and wrote out fort.222.inv file')
         endif
         ! ice file (if we are using it) (fort.225.grb2)
         if (NCICE.eq.14) then
            inquire(file=Cinv,exist=Fexists)
            if (.not.Fexists.and.iret.eq.0) then
               iret = grb2_mk_inv(Cfile, Cinv)
               if (iret.ne.0) then
                  call allMessage(ERROR,
     &                         'Fatal error in reading fort.225.grb2.')
               endif
               call allMessage(ECHO,'successfully read fort.225.grb2 '
     &                         //'and wrote out fort.225.inv file')
            endif
         endif
#ifdef CMPI
      ENDIF
#endif
      CALL BcastToLocal_Int(iret)
      if (iret.ne.0) then
         call nws14Terminate()
      endif
      RETURN
#endif
C----------------------------------------------------------------------
      END SUBROUTINE NWS14GRB2_MAKE_INV
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE NWS14NC_READ_F22(NWS)
C----------------------------------------------------------------------
C     CPB 10/2023: This subroutine reads the fort.22 file for NWS14 to
C     get netcdf variable names as well as starting indices for reading
C     in the met data. If we are using grib2 it just returns
C----------------------------------------------------------------------
      USE SIZES, ONLY : MYPROC
      USE GLOBAL, ONLY : NCICE
      implicit none
      INTEGER,INTENT(IN) :: NWS
      logical :: Fexists
      integer :: ierr, PvarInd
      ! return if we are using grib2
      IF (grb2flag) THEN
         RETURN
      ENDIF
      CALL setMessageSource("NWS14NC_READ_F22")
      IF (NWS.EQ.-14) THEN
         ! For inset OWI winds with netcdf files something needs to be
         ! done but I haven't programmed it yet.
         CALL allMessage(ERROR,"NWS = -14 has not yet been set up to "
     &                   //"use netcdf files. Try again with grib2 "
     &                   //"background wind files.")
         CALL nws14Terminate()
      ENDIF
      ! open the fort.22
      OPEN(22,file='fort.22',ACTION='READ',IOSTAT=ierr)
      ! if we had trouble opening then abort
      IF (ierr.NE.0) THEN
         call allMessage(ERROR,'Unable to open fort.22.')
         call nws14Terminate()
      ENDIF
      read(22,*) Tdim    ! Time dimension name
      read(22,*) Tvar    ! Time variable name
      read(22,*) Tformat ! Format of time datestr
      read(22,*) Londim  ! Lon dimension name
      read(22,*) Lonvar  ! Lon var name
      read(22,*) Latdim  ! Lat dimension name
      read(22,*) Latvar  ! Lat var name
      read(22,*) Pvar    ! Pressure var name
      read(22,*) Uvar    ! U10 var name
      read(22,*) Vvar    ! V10 Var name
      if (NCICE.eq.14) then
         read(22,*) Cvar ! Ice var name
      endif
      close(22)
      PvarInd = index(Pvar,'HPa')
      if (PvarInd > 0) then
         rhoWat0g = rhoWat0g/100d0
         Pvar = Pvar(1:PvarInd-1)
         write(16,*) 'Pressure is in Hpa, New Pvar = ',trim(Pvar)
      endif
      ! get starting time index for netcdf file
      CALL READNWS14_NC_StaTime(Pfile, PfileNCID)
      CALL unsetMessageSource()
      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE NWS14NC_READ_F22
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE NWS14_CALC_INTERP_WTS()
C----------------------------------------------------------------------
C     CPB 10/2023: This subroutine calculates the interpolation indices
C     (indp) and weights (weightsp) to interpolate from the grib2 or
C     netcdf input meterological forcing to the ADCIRC grid. NWS = 14
C     was originally set up to use NOAA CFSv2/GFS met forcing products.
C     These meteorological products provide two different resolutions
C     for wind forcing but not for pressure (or ice). To limit
C     pre-processing necessary to use these products, this subroutine
C     checks to see if the wind and pressure have the same dimensions.
C     If so, it only allocates one set of interpolant indices (indp) and
C     weights (weightsp). If they are different then it stores two
C     different sets. If necessary, it also checks the dimensions of the
C     ice forcing and sets the interpolants for ice to be equal to
C     either the wind or pressure forcing depending on the grid.
C----------------------------------------------------------------------
      use mesh, only: NP, SLAM, SFEA, bl_interp, bl_interp2
      use kdtree2_module, only: kdtree2, kdtree2_create
      USE GLOBAL, ONLY : NCICE
      implicit none
      integer :: i, ii, nxp, nyp, indt(4), NTkeep, xi, yi, cc
      real,allocatable,dimension(:,:) :: Pmsl, U10, V10, Cice, lat, lon
      real(8), allocatable ::  lonv(:), latv(:), XY(:,:)
      real(8) :: xx, yy, wt(4)
      logical :: regular_grid
      type(kdtree2), pointer :: tree
      ! determine if wind/pressure/ice grids match or not.
      ! assume grids are the same to start
      ub = 1
      ubc = 1 ! index where we store the ice data (will match either
              ! wind or pressure). Assume match pressure first.
      ! get pressure dimensions
      CALL READNWS14LatLon(Pfile,lat,lon,Pinv,Pvar,PfileNCID)
      nxp = ubound(lon,1); nyp = ubound(lat,2)
      ! get wind dimensions and see if they match pressure
      CALL READNWS14LatLon(Wfile,lat,lon,Winv,Uvar,WfileNCID)
      if (nxp.NE.ubound(lon,1).OR.nyp.NE.ubound(lat,2)) then
         ! if we don't match both lat and lon then we need to store two
         ! sets of interpolation indices/weights
         ub = 2
         nxp = ubound(lon,1)
         nyp = ubound(lat,2)
      endif
      ! If winds and pressure do not match, then check ice file to see
      ! whether we should use wind or pressure interpolants. If wind and
      ! pressure DO match then we assume ice also matches.
      IF (NCICE.EQ.14.AND.ub.EQ.2) THEN
         CALL READNWS14LatLon(Cfile,lat,lon,Cinv,Cvar,CfileNCID)
         IF (nxp.EQ.ubound(lon,1).AND.nyp.EQ.ubound(lat,2)) THEN
            ! match winds
            ubc = 2
         ENDIF
      ENDIF
      ! Allocate the indices and weights
      allocate(indp(ub,4,np),weightsp(ub,4,np))
      ! Make the interpolant weights. If w
      do ii = 1,ub
         if (ii.eq.1) then
            CALL READNWS14LatLon(Pfile,lat,lon,Pinv,Pvar,PfileNCID)
         else
            CALL READNWS14LatLon(Wfile,lat,lon,Winv,Uvar,WfileNCID)
         endif
         nxp = ubound(lon,1); nyp = ubound(lat,2)
         ! Check if structured regular..
         if (abs(lon(1,1) - lon(1,nyp)) < 1d-6) then
            regular_grid = .true.
            allocate(lonv(nxp),latv(nyp))
            lonv = lon(:,1)
            latv = lat(1,:)
         else
            regular_grid = .false.
            ! Convert point matrices to 2 X (NXP*NYP) point vector
            allocate(XY(2,nxp*nyp))
            cc = 0
            do xi = 1,nxp
               do yi = 1,nyp
                  cc = cc + 1
                  XY(1,cc) = lon(xi,yi)
                  XY(2,cc) = lat(xi,yi)
               enddo
            enddo
            ! Create the search tree
            tree => kdtree2_create(XY)
         endif
         ! Loop over each node and interpolate
         do i = 1,np
            xx = rad2deg*slam(i)
            ! Convert our numbers if grids are 0 to 360
            if (maxval(lon).gt.180d0.and.xx < 0d0) xx = xx + 360d0
            yy = rad2deg*sfea(i)
            if (regular_grid) then
               call bl_interp(nxp,lonv,nyp,latv,xx,yy,indt,wt)
            else
               call bl_interp2(nxp,nyp,lon,lat,xx,yy,indt,wt,tree)
            endif
            indp(ii,:,i) = indt
            weightsp(ii,:,i) = wt
         enddo
         if (regular_grid) then
            deallocate(lon,lat,latv,lonv)
         else
            deallocate(lon,lat,XY)
         endif
      enddo
C----------------------------------------------------------------------
      END SUBROUTINE NWS14_CALC_INTERP_WTS
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE NWS14GET(WTIMINC,WVNX2,WVNY2,PRN2,CICE2)
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      use mesh, only: NP
      use mod_datetime, only: operator(+)
      implicit none

      real(8),intent(in) :: WTIMINC
      real(8),intent(out) :: PRN2(NP), WVNX2(NP), WVNY2(NP)
      real(8),intent(out),optional :: CICE2(NP)
      integer :: i, NTkeep
      real,allocatable,dimension(:,:) :: Pmsl, U10, V10, Cice

      !
      CALL setMessageSource('NWS14GET')
      ! Show the date in str_date format
      call logMessage(ECHO,'CurDT strng: '//
     &                CurDT%strftime("%Y-%m-%d %H:%M"))


      ! Get the Pressure
      !':PRMSL:mean sea level:'
      call READNWS14(Pfile,Pmsl,Pinv,Pvar,PfileNCID)

      ! Get the U-winds
      call READNWS14(Wfile,U10,Winv,Uvar,WfileNCID)

      ! Get the V-winds
      call READNWS14(Wfile1,V10,Winv,Vvar,Wfile1NCID)

      ! Get the ice concentration if required
      if (present(CICE2)) then
         NTkeep = NT; NT = NTC
         call READNWS14(Cfile,Cice,Cinv,Cvar,CfileNCID)
         NT = NTkeep
      endif


      ! Doing the interpolation from their grid to ours
      do i = 1,np
         if (indp(1,1,i) < 0) then
            PRN2(i) = PRBCKGRND_MH2O
         else
            PRN2(i) = ( Pmsl(indp(1,1,i),indp(1,2,i))*weightsp(1,1,i) +
     &             Pmsl(indp(1,3,i),indp(1,2,i))*weightsp(1,2,i)      +
     &             Pmsl(indp(1,1,i),indp(1,4,i))*weightsp(1,3,i)      +
     &             Pmsl(indp(1,3,i),indp(1,4,i))*weightsp(1,4,i) )
     &              / rhoWat0g
         endif
         if (indp(ub,1,i) < 0) then
            WVNX2(i) = 0d0; WVNY2(i) = 0d0
         else
            WVNX2(i) = U10(indp(ub,1,i),indp(ub,2,i))*weightsp(ub,1,i) +
     &              U10(indp(ub,3,i),indp(ub,2,i))*weightsp(ub,2,i)    +
     &              U10(indp(ub,1,i),indp(ub,4,i))*weightsp(ub,3,i)    +
     &              U10(indp(ub,3,i),indp(ub,4,i))*weightsp(ub,4,i)
            WVNY2(i) = V10(indp(ub,1,i),indp(ub,2,i))*weightsp(ub,1,i) +
     &              V10(indp(ub,3,i),indp(ub,2,i))*weightsp(ub,2,i)    +
     &              V10(indp(ub,1,i),indp(ub,4,i))*weightsp(ub,3,i)    +
     &              V10(indp(ub,3,i),indp(ub,4,i))*weightsp(ub,4,i)
         endif
      enddo
      deallocate(Pmsl,U10,V10)

       IF(found_InterpInset_Namelist.AND.UseInterpInset) THEN

      ! Insert insets in netcdf format  (Aman Tejaswi)
            CALL INSET_NC_INPUT(PRN2,WVNX2,WVNY2)

       ELSE IF(found_InterpInset_Namelist.AND.(.NOT.UseInterpInset)) THEN

            CALL logMessage(ECHO,'UseInterpInset set to F, Please Check
     &                          Namelist if you want to use inset files')
       ENDIF
!      ! Add WTIMINC on CurDT for next WTIME
      CurDT = CurDT + t_timedelta(minutes=nint(WTIMINC/60d0))
      ! Next time index for the netcdf reading
      NT = NT + 1

      ! If ice concentration is present
      if (present(CICE2)) then
         do i = 1,np
            if (indp(ubc,1,i) < 0) then
               CICE2(i) = 0d0;
            else
               CICE2(i) =
     &             Cice(indp(ubc,1,i),indp(ubc,2,i))*weightsp(ubc,1,i) +
     &             Cice(indp(ubc,3,i),indp(ubc,2,i))*weightsp(ubc,2,i) +
     &             Cice(indp(ubc,1,i),indp(ubc,4,i))*weightsp(ubc,3,i) +
     &             Cice(indp(ubc,3,i),indp(ubc,4,i))*weightsp(ubc,4,i)
            endif
         enddo
         deallocate(Cice)
         ! Next ice time index for the netcdf reading
         NTC = NTC + 1
      endif
      CALL unsetMessageSource()
      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE NWS14GET
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE READNWS14(fileN,data2,invN,var,ncid)
C----------------------------------------------------------------------
C     CPB 10/2023: Added to make NWS = 14 more more readable. Simply
C     calls the appropriate subroutine based on the filetype we are
C     using.
C----------------------------------------------------------------------
      implicit none
      CHARACTER(LEN=200),INTENT(IN) :: fileN
      REAL,ALLOCATABLE,INTENT(OUT) :: data2(:,:)
      CHARACTER(LEN=200),INTENT(IN) :: invN, var
      INTEGER,INTENT(INOUT) :: ncid
      IF (.NOT.grb2flag) THEN
         CALL READNWS14_netCDF(fileN, ncid, var, data2)
      ELSEIF (grb2flag) THEN
         CALL READNWS14_grib2(fileN, invN, var, data2)
      ELSE
         ! should be unreachable
         CALL nws14Terminate()
      ENDIF
C----------------------------------------------------------------------
      END SUBROUTINE READNWS14
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE READNWS14_grib2(fileN,invN,var,data2)
C----------------------------------------------------------------------
C       CPB 10/2023: Reads in grib2 format met forcing for NWS = 14.
C       NOTE: reads in on Proc 0 and broadcasts.
C----------------------------------------------------------------------
#ifdef GRIB2API
      use wgrib2api
#endif
      USE GL2LOC_MAPPING ,ONLY : BcastToLocal_2DRealArray,
     &                           BcastToLocal_Int
      USE SIZES, ONLY : MYPROC
      implicit none
      integer :: iret, iter
      character(len=200),intent(in) :: var, fileN, invN
      real,allocatable,intent(out) :: data2(:,:)
      character(len=200) :: str_date
      integer :: NX, NY
      CALL setMessageSource("READNWS14_grib2")
      IF (MYPROC.EQ.0) THEN
         ! Getting the date in str_date format
         str_date = ':start_FT='//trim(CurDT%strftime("%Y%m%d%H"))//'0000:'
         iret = -1; iter = 0
         do while (iret.le.0.and.iter < 10)
            if (iter > 0) then
               call logMessage(WARNING,'Trying to read again '//var)
               call sleep(5)
            endif
#ifdef GRIB2API
            iret = grb2_inq(fileN,invN,var,str_date,data2=data2)
            if (iret.gt.1.and.iter.eq.0) then
               ! May have two entries bcause of forecast/DA overlap
               call logMessage(DEBUG,'> 1 msg, trying to read '//var)
               iret = grb2_inq(fileN,invN,var,':anl:',
     &                         str_date,data2=data2)
            endif
#endif
            iter = iter + 1
         enddo
         if (iret.gt.0) then
            call logMessage(DEBUG,'Successfully read '//var)
         else
            call ArnoldSchwarzenegger(iret,var,fileN)
         endif
         NX = UBOUND(data2,1)
         NY = UBOUND(data2,2)
      ENDIF
      CALL BcastToLocal_Int(NX)
      CALL BcastToLocal_INT(NY)
      IF (MYPROC.NE.0) THEN
         ALLOCATE( data2(NX,NY) )
      ENDIF
      CALL BcastToLocal_2DRealArray(data2,NX,NY)
      CALL unsetMessageSource()
      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE READNWS14_grib2
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE READNWS14_netCDF(fileN, NC_ID,var,data2)
C----------------------------------------------------------------------
C       CPB 10/2023: Reads in netCDF format met forcing for NWS = 14.
C       NOTE: reads in on Proc 0 and broadcasts.
C----------------------------------------------------------------------
#ifdef ADCNETCDF
      use netcdf
      use netcdf_error, only: check_err
#endif
      USE GL2LOC_MAPPING ,ONLY : BcastToLocal_2DRealArray,
     &                           BcastToLocal_Int
      USE SIZES, ONLY : MYPROC
      implicit none

      character(len=*),intent(in):: fileN
      character(len=200),intent(in) :: var
      real,allocatable,intent(out) :: data2(:,:)
      INTEGER,INTENT(INOUT) :: NC_ID
      integer :: i, Temp_ID, Lat_ID, Lon_ID, NX, NY

      CALL setMessageSource('READNWS14_netCDF')

#ifdef ADCNETCDF
      IF ( read_NWS14_NetCdf_using_core_0 ) THEN
        ! Use core zero to read data
        ! and then broadcast them
        IF (MYPROC.EQ.0) THEN
           call Check_err(NF90_INQ_DIMID(NC_ID,Latdim,Temp_ID))
           call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,len=NY))
           call Check_err(NF90_INQ_DIMID(NC_ID,Londim,Temp_ID))
           call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,len=NX))
        ENDIF
        CALL BcastToLocal_Int(NX)
        CALL BcastToLocal_Int(NY)
        allocate(data2(NX,NY))
        IF (MYPROC.EQ.0) THEN
           call Check_err(NF90_INQ_VARID(NC_ID,var,Temp_ID))
           call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,data2,
     &                  start=[1, 1, NT],count=[NX, NY, 1]))
        ENDIF
        CALL BCastToLocal_2DRealArray(Data2,NX,NY)
      ELSE
        ! Each core read data
        call Check_err(NF90_OPEN(fileN,nf90_nowrite,NC_ID))
        call Check_err(NF90_INQ_DIMID(NC_ID,Latdim,Temp_ID))
        call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,len=NY))
        call Check_err(NF90_INQ_DIMID(NC_ID,Londim,Temp_ID))
        call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,len=NX))

        allocate(data2(NX,NY))
        call Check_err(NF90_INQ_VARID(NC_ID,var,Temp_ID))
        call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,data2,
     &                  start=[1, 1, NT],count=[NX, NY, 1]))

        call Check_err(NF90_CLOSE(NC_ID))
      ENDIF

#endif
      CALL unsetMessageSource()
      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE READNWS14_netCDF
C----------------------------------------------------------------------


C---------------------------------------------------------------------
      SUBROUTINE INSET_NC_INPUT(PRN2,WVNX2,WVNY2)
C----------------------------------------------------------------------
C     Aman Tejaswi (11th Feb 2025)
C     This subroutine uses high res inset netcdf wind and pressure data.
C     The way it is implemented that first the global background wind (ex. GFS,
C     ERA5) and pressure are interpolated on the ADCIRC mesh and then
C     wherever highres wind is available, it interpolates the inset data on the
C     Finite Element grid. The file 'inset_netcdf.in' is used to get the
C     lat lon variable names as well as variable names.
C     Expected Structure of The File Is Like This: (example)
C
C-------------------------------------
C              lonDimName
C              latDimName
C              lonVarName
C              latVarName
C              InsetPresureVarName
C              InsetUWindVarName
C              InsetVWindVarName
C-------------------------------------
C     For now the code assumes that your background wind and your high res
C     inset has same temporal coverage.
C     The code only works with NWS=14 (netcdf). The implementation is controlled
C     by Namelist -- &InterpInset, UseInterpInset=T. The
C     subroutine can be used in the sameway as OWI_NETCDF with and
C     advantage that one doesn't need to write global and regional
C     wind/pressure data in a particular (OWI) format as sometimes it
C     becomes tedious if running decadal simulations.
C----------------------------------------------------------------------
       USE GLOBAL, ONLY : PFile_inset, Wfile_inset, InsetControlFile

       IMPLICIT NONE

       REAL(8), INTENT(INOUT) :: PRN2(:), WVNX2(:), WVNY2(:)

       INTEGER :: TempNCID_P, TempNCID_W


       REAL(8), ALLOCATABLE, DIMENSION(:,:) :: InsetPData
       REAL(8), ALLOCATABLE, DIMENSION(:,:) :: InsetUData, InsetVData

       CALL setMessageSource("INSET_NC_INPUT")
       Pfile_inset= TRIM(ADJUSTL(Pfile_inset))
       Wfile_inset= TRIM(ADJUSTL(Wfile_inset))

         CALL READINSET_NC(Pfile_inset
     &       ,TempNCID_P,InsetLatDim,InsetLat,InsetLonDim,InsetLon,InsetPvar,InsetPData)


         CALL READINSET_NC(Wfile_inset
     &       ,TempNCID_W,InsetLatDim,InsetLat,InsetLonDim,InsetLon,InsetUVar,InsetUData)

         CALL READINSET_NC(Wfile_inset
     &       ,TempNCID_W,InsetLatDim,InsetLat,InsetLonDim,InsetLon,InsetVvar,InsetVData)


         CALL INTERP_INSET_GRID(InsetPData,InsetUData,InsetVData,PRN2,WVNX2,WVNY2)

         CALL unsetMessageSource()

C----------------------------------------------------------------------
      END SUBROUTINE INSET_NC_INPUT

C----------------------------------------------------------------------
      SUBROUTINE INIT_INSET_NC()

       USE SIZES, ONLY : MYPROC, GBLINPUTDIR
       USE GLOBAL, ONLY : InsetControlFile

       IMPLICIT NONE

       LOGICAL :: InsetFileExist
       INTEGER :: IERR


       CALL setMessageSource("INIT_INSET_NC")

       INQUIRE(FILE=TRIM(GBLINPUTDIR)//'/'//TRIM(ADJUSTL(InsetControlFile)), EXIST=InsetFileExist)
       IF(.NOT.InsetFileExist) THEN
           write(scratchMessage,'(A)') "Unable to find InsetControlFile: "//TRIM(InsetControlFile)
           call allMessage(ERROR, scratchMessage)
           call nws14Terminate()
       ENDIF

       ! open the inset_netcdf.in file
       OPEN(2200,FILE=TRIM(GBLINPUTDIR)//'/'//TRIM(ADJUSTL(InsetControlFile)),ACTION='READ',IOSTAT=IERR)
       ! if we had trouble opening then abort
       IF (IERR.NE.0) THEN
          write(scratchMessage,'(A)') 'Unable to open insetControlFile: '//TRIM(InsetControlFile)
          call allMessage(ERROR,scratchMessage)
          call nws14Terminate()
       ENDIF

       READ(2200, '(A)', IOSTAT=IERR) InsetLonDim
       READ(2200, '(A)', IOSTAT=IERR) InsetLatDim
       READ(2200, '(A)', IOSTAT=IERR) InsetLon
       READ(2200, '(A)', IOSTAT=IERR) InsetLat
       READ(2200, '(A)', IOSTAT=IERR) InsetPvar
       READ(2200, '(A)', IOSTAT=IERR) InsetUVar
       READ(2200, '(A)', IOSTAT=IERR) InsetVVar

       CLOSE(2200) !closing the file here..

      END SUBROUTINE INIT_INSET_NC
C----------------------------------------------------------------------




C----------------------------------------------------------------------
C     Aman Tejaswi (11 Feb 2025)
C     Read High Res-Inset Winds in NETCDF format and Interpolate onto
C     ADCIRC grid. These insets can be used particularly with a
C     background global wind (e.g GFS) which is read also in NETCDF. The
C     advantage of using this is specific to GLOBAL MODEL (GSTOFS).

C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE READINSET_NC(InsetFileN,InsetID,InsetLatDim,InsetLat,InsetLonDim,InsetLon,InsetVar,InsetOUT)

#ifdef ADCNETCDF
      use netcdf
      use netcdf_error, only: check_err
#endif

      USE GL2LOC_MAPPING ,ONLY : BcastToLocal_2DRealArray,
     &                           BcastToLocal_Int
      USE SIZES, ONLY : MYPROC, GBLINPUTDIR
      USE GLOBAL, ONLY: scratchMessage, allMessage, ERROR

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)      :: InsetFileN
      CHARACTER(LEN=*), INTENT(IN)      :: InsetVar, InsetLat, InsetLon,InsetLatDim, InsetLonDim
      INTEGER                           :: InsetID
      REAL(8), ALLOCATABLE, INTENT(OUT) :: InsetOUT(:,:)
      INTEGER                           :: InsetLoc, INSET_NX, INSET_NY, VarLoc,alloc_err

         CALL setMessageSource('READINSET_NC')
#ifdef ADCNETCDF

         CALL Check_err(NF90_OPEN(InsetFileN,nf90_nowrite,InsetID))

         CALL Check_err(NF90_INQ_DIMID(InsetID,InsetLonDim,InsetLoc))
         CALL Check_err(NF90_INQUIRE_DIMENSION(InsetID,InsetLoc,len=INSET_NX))

         CALL Check_err(NF90_INQ_DIMID(InsetID,InsetLatDim,InsetLoc))
         CALL Check_err(NF90_INQUIRE_DIMENSION(InsetID,InsetLoc,len=INSET_NY))

         ALLOCATE(InsetOUT(INSET_NX,INSET_NY),STAT=alloc_err)

         IF (alloc_err /= 0) THEN
            WRITE(scratchMessage,'(A)') "Unable to allocate memory for InsetOUT"
            call allMessage(ERROR, scratchMessage)
            call nws14Terminate()
         ENDIF


         CALL Check_err(NF90_INQ_VARID(InsetID,InsetVar,VarLoc))
         CALL Check_err(NF90_GET_VAR(InsetID,VarLoc,InsetOUT,
     &                  start=[1, 1, NT],count=[Inset_NX, Inset_NY, 1]))


         CALL Check_err(NF90_CLOSE(InsetID))

#endif
         CALL unsetMessageSource()

         RETURN

       END SUBROUTINE READINSET_NC

!--------------------------------------------------------------------------
       SUBROUTINE CALC_INSET_WTS()


        use mesh, only: NP, SLAM, SFEA, bl_interp, bl_interp2
        use kdtree2_module, only: kdtree2, kdtree2_create,kdtree2_destroy
        USE GLOBAL, ONLY : Pfile_inset
        IMPLICIT NONE

        integer :: i, ii, nxp, nyp, indt(4), NTkeep, xi, yi,cc,ncid_p
        real,allocatable,dimension(:,:) :: lat, lon
        real(8), allocatable ::  lonv(:), latv(:), XY(:,:)
        real(8) :: xx, yy, wt(4)
        logical :: regular_grid
        type(kdtree2), pointer :: tree
        !
        ! assume grids are the same to start
        ub = 1

        CALL READ_INSET_LatLon(Pfile_inset,NCID_p,lat,InsetLatDim,InsetLat,lon,InsetLonDim,InsetLon,nxp,nyp)
        nxp = ubound(lon,1); nyp = ubound(lat,2)



        if (nxp.NE.ubound(lon,1).OR.nyp.NE.ubound(lat,2)) then
           ! if we don't match both lat and lon then we need to store two
           ! sets of interpolation indices/weights
           ub = 2
           nxp = ubound(lon,1)
           nyp = ubound(lat,2)
        endif

        IF ( .NOT.allocated(inset_indp) .and. .not.allocated(inset_weight)) then
        ! Allocate the indices and weights
        allocate(inset_indp(ub,4,np),inset_weight(ub,4,np))
        endif
        ! Make the interpolant weights. If w
        do ii = 1,ub

           ! Check if structured regular..
           if (abs(lon(1,1) - lon(1,nyp)) < 1d-6) then
              regular_grid = .true.
              if(.not.allocated(lonv).and..not.allocated(latv)) then
              allocate(lonv(nxp),latv(nyp))
               endif
              lonv = lon(:,1)
              latv = lat(1,:)

           else
              regular_grid = .false.
              ! Convert point matrices to 2 X (NXP*NYP) point vector
              if(.not.allocated(xy)) then
              allocate(XY(2,nxp*nyp))
              endif
              cc = 0
              do xi = 1,nxp
                 do yi = 1,nyp
                    cc = cc + 1
                    XY(1,cc) = lon(xi,yi)
                    XY(2,cc) = lat(xi,yi)
                 enddo
              enddo
              ! Create the search tree
              tree => kdtree2_create(XY)
           endif
           ! Loop over each node and interpolate
           do i = 1,np
              xx = rad2deg*slam(i)
              ! Convert our numbers if grids are 0 to 360
              if (maxval(lon).gt.180d0.and.xx < 0d0) xx = xx + 360d0
              yy = rad2deg*sfea(i)
              if (regular_grid) then
                 call bl_interp(nxp,lonv,nyp,latv,xx,yy,indt,wt)
              else
                 call bl_interp2(nxp,nyp,lon,lat,xx,yy,indt,wt,tree)
              endif
              inset_indp(ii,:,i) = indt
              inset_weight(ii,:,i) = wt
           enddo
           if (regular_grid) then
              deallocate(lon,lat,latv,lonv)
           else
              deallocate(lon,lat,XY)
           endif

            if (associated(tree)) then
              call kdtree2_destroy(tree)
            nullify(tree)
            endif

        enddo

       END SUBROUTINE CALC_inset_wts

C------------------------------------------------------------------------
       SUBROUTINE INTERP_INSET_GRID(InsetPData,InsetUData,InsetVData,PRN2,WVNX2,WVNY2)

C---------------------------------------------------------------------------------------
        USE MESH, ONLY : NP
        IMPLICIT NONE

        INTEGER :: i
        REAL(8), INTENT(IN) :: InsetPData(:,:), InsetUData(:,:), InsetVData(:,:)
        REAL(8), INTENT(INOUT) :: PRN2(:), WVNX2(:), WVNY2(:)

        ! Doing the interpolation from inset grid to finite ele grid
        do i = 1,np
         if (inset_indp(1,1,i) < 0) then
            PRN2(i) = PRN2(i)     !Keeping the background global pressure
         else

            PRN2(i) = (InsetPData(inset_indp(1,1,i),inset_indp(1,2,i))*inset_weight(1,1,i) +
     &            InsetPData(inset_indp(1,3,i),inset_indp(1,2,i))*inset_weight(1,2,i)      +
     &            InsetPData(inset_indp(1,1,i),inset_indp(1,4,i))*inset_weight(1,3,i)      +
     &            InsetPData(inset_indp(1,3,i),inset_indp(1,4,i))*inset_weight(1,4,i) )
     &              / rhoWat0g
         endif
         if (inset_indp(ub,1,i) < 0) then
            WVNX2(i) = WVNX2(i); WVNY2(i) = WVNY2(i)  !Keeping the background global wind
         else

            WVNX2(i) =InsetUData(inset_indp(1,1,i),inset_indp(1,2,i))*inset_weight(1,1,i) +
     &             InsetUData(inset_indp(1,3,i),inset_indp(1,2,i))*inset_weight(1,2,i)    +
     &             InsetUData(inset_indp(1,1,i),inset_indp(1,4,i))*inset_weight(1,3,i) +
     &             InsetUData(inset_indp(1,3,i),inset_indp(1,4,i))*inset_weight(1,4,i)


            WVNY2(i) =InsetVData(inset_indp(1,1,i),inset_indp(1,2,i))*inset_weight(1,1,i) +
     &             InsetVData(inset_indp(1,3,i),inset_indp(1,2,i))*inset_weight(1,2,i)    +
     &             InsetVData(inset_indp(1,1,i),inset_indp(1,4,i))*inset_weight(1,3,i)    +
     &             InsetVData(inset_indp(1,3,i),inset_indp(1,4,i))*inset_weight(1,4,i)
         endif
        enddo


        RETURN

C-----------------------------------------------------------------------------------------

       END SUBROUTINE INTERP_INSET_GRID
C-----------------------------------------------------------------------------------


      SUBROUTINE READ_INSET_LATLON(ncfile,NCIDf,lat_ins,
     &                  InsetLatDim,InsetLat,lon_ins,InsetLonDim,InsetLon,nxp,nyp)

#ifdef ADCNETCDF
        use netcdf
        use netcdf_error, only: check_err
#endif
        USE GL2LOC_MAPPING ,ONLY : BcastToLocal_2DRealArray,BcastToLocal_Int
        USE SIZES, ONLY : MYPROC, GBLINPUTDIR

        IMPLICIT NONE

        CHARACTER(len=*), INTENT(IN) :: ncfile  ! NetCDF file name
        REAL, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: lat_ins,lon_ins

        CHARACTER(LEN=80), INTENT(IN)  :: InsetLat, InsetLon, InsetLatDim,InsetLonDim
        INTEGER, INTENT(OUT) :: NXP, NYP
        INTEGER              :: NCIDf

        INTEGER :: i, Temp_ID, Lat_ID, Lon_ID, ndims

        CALL setMessageSource("READ_INSET_LATLON")
#ifdef ADCNETCDF

         CALL Check_err(NF90_OPEN(ncfile,nf90_nowrite,NCIDf))

         call Check_err(NF90_INQ_DIMID(ncidf,InsetLatDim,Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(ncidf,Temp_ID,
     &                                         len=NYP))
         call Check_err(NF90_INQ_DIMID(ncidf,InsetLonDim,Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(ncidf,Temp_ID,
     &                                         len=NXP))

         allocate(lon_ins(NXP,NYP),lat_ins(NXP,NYP))


         call Check_err(NF90_INQ_VARID(ncidf,InsetLat,Lat_ID))
         call Check_err(nf90_inquire_variable(ncidf,Lat_ID,ndims=ndims))
         call Check_err(NF90_INQ_VARID(ncidf,InsetLon,Lon_ID))

         if (ndims.eq.2) then
            call Check_err(NF90_GET_VAR(ncidf,Lat_ID,lat_ins))
            call Check_err(NF90_GET_VAR(ncidf,Lon_ID,lon_ins))
         else
            call Check_err(NF90_GET_VAR(ncidf,Lat_ID,lat_ins(1,:)))
            call Check_err(NF90_GET_VAR(ncidf,Lon_ID,lon_ins(:,1)))
            do i = 2,NXP
               lat_ins(i,:) = lat_ins(1,:)
            enddo
            do i = 2,NYP
               lon_ins(:,i) = lon_ins(:,1)
            enddo
         endif

        CALL Check_err(NF90_CLOSE(NCIDf))
#else
        NXP = 0
        NYP = 0
#endif
        CALL unsetMessageSource()
      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE READ_INSET_LATLON


C----------------------------------------------------------------------
      SUBROUTINE READNWS14LatLon(fileN,lat,lon,invN,var,ncid)
C----------------------------------------------------------------------
C     CPB 10/2023: Added to make NWS = 14 more more readable. Simply
C     calls the appropriate subroutine based on the filetype we are
C     using.
C----------------------------------------------------------------------
      implicit none
      CHARACTER(len=200),INTENT(IN) :: fileN
      REAL,ALLOCATABLE,INTENT(OUT) :: lat(:,:),lon(:,:)
      CHARACTER(LEN=200),INTENT(IN),OPTIONAL :: var, invN
      INTEGER,INTENT(INOUT),OPTIONAL :: ncid
      IF (grb2flag) THEN
         CALL READNWS14LatLon_grib2(fileN,invN,var,lat,lon)
      ELSEIF (.NOT.grb2flag) THEN
         CALL READNWS14LatLon_netCDF(fileN,ncid,lat,lon)
      ELSE
         ! should be unreachable
         CALL nws14Terminate()
      ENDIF
C----------------------------------------------------------------------
      END SUBROUTINE READNWS14LatLon
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE READNWS14LatLon_grib2(fileN,invN,var,lat,lon)
C----------------------------------------------------------------------
C     CPB 10/2023: Reads in the lat and lon from a grib2 format
C     meteorological file. NOTE: reads in on Proc 0 and broadcasts
C----------------------------------------------------------------------
#ifdef GRIB2API
      use wgrib2api
#endif
      USE GL2LOC_MAPPING ,ONLY : BcastToLocal_2DRealArray,
     &                           BcastToLocal_Int
      USE SIZES, ONLY : MYPROC
      implicit none
      integer :: iret, iter
      character(len=200),intent(in) :: var, fileN, invN
      real,allocatable,intent(out) :: lat(:,:), lon(:,: )
      integer :: NX, NY
      character(len=200) :: str_date
      CALL setMessageSource("READNWS14LatLon_grib2")
      ! Getting the date in str_date format
      str_date = ':start_FT='//trim(CurDT%strftime("%Y%m%d%H"))//'0000:'
      ! read in on core 0
      IF (MYPROC.EQ.0) THEN
         iret = -1; iter = 0
         do while (iret.le.0.and.iter < 10)
            if (iter > 0) then
               call logMessage(WARNING,'Trying to read again '//var)
               call sleep(5)
            endif
#ifdef GRIB2API
            iret = grb2_inq(fileN,invN,var,str_date,
     &                      lat=lat,lon=lon)
#endif
            iter = iter + 1
         enddo
         if (iret.gt.0) then
            call logMessage(ECHO,'Successfully read LatLon '//var)
         else
            call ArnoldSchwarzenegger(iret,var,fileN)
         endif
      ENDIF
      ! broadcast data
      IF (MYPROC.EQ.0) THEN
         NX = UBOUND(lon,1)
         NY = UBOUND(lat,2)
      ENDIF
      CALL BcastToLocaL_Int(NX)
      Call BcastToLocal_Int(NY)
      IF (MYPROC.NE.0) THEN
         ALLOCATE( lat(NX,NY), lon(NX,NY) )
      ENDIF
      ! latitude
      CALL BcastToLocal_2DRealArray(lat,NX,NY)
      ! longitude
      CALL BcastToLocal_2DRealArray(lon,NX,NY)
      CALL unsetMessageSource()
      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE READNWS14LatLon_grib2
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE READNWS14LatLon_netCDF(fileN,ncid,lat,lon)
C----------------------------------------------------------------------
C     CPB 10/2023: Reads in the lat and lon from a netCDF format
C     meteorological forcing file. NOTE: reads in on Proc 0 and
C     broadcasts.
C----------------------------------------------------------------------
#ifdef ADCNETCDF
      use netcdf
      use netcdf_error, only: check_err
#endif
      USE GL2LOC_MAPPING ,ONLY : BcastToLocal_2DRealArray,
     &                           BcastToLocal_Int
      USE SIZES, ONLY : MYPROC
      implicit none
      INTEGER,INTENT(INOUT) :: ncid ! already opened ncid tag
      character(len=200),intent(in) :: fileN
      real,allocatable,intent(out) :: lat(:,:), lon(:,:)
      integer :: i, Temp_ID, Lat_ID, Lon_ID, NX, NY, ndims
      CALL setMessageSource("READNWS14LatLon_netCDF")
#ifdef ADCNETCDF
      IF (MYPROC.EQ.0) THEN

         IF ( .NOT. read_NWS14_NetCdf_using_core_0 ) THEN
           call Check_err(NF90_OPEN(FileN,nf90_nowrite,ncid))
         END IF

         call Check_err(NF90_INQ_DIMID(ncid,Latdim,Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(ncid,Temp_ID,
     &                                         len=NY))
         call Check_err(NF90_INQ_DIMID(ncid,Londim,Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(ncid,Temp_ID,
     &                                         len=NX))
         allocate(lon(NX,NY),lat(NX,NY))
         call Check_err(NF90_INQ_VARID(ncid,Latvar,Lat_ID))
         call Check_err(nf90_inquire_variable(ncid,Lat_ID,
     &               ndims=ndims))
         call Check_err(NF90_INQ_VARID(ncid,Lonvar,Lon_ID))
         if (ndims.eq.3) then
            call Check_err(NF90_GET_VAR(ncid,Lat_ID,lat,
     &                     start=[1, 1, NT],count=[NX, NY, 1]))
            call Check_err(NF90_GET_VAR(ncid,Lon_ID,lon,
     &                     start=[1, 1, NT],count=[NX, NY, 1]))
         elseif (ndims.eq.2) then
            call Check_err(NF90_GET_VAR(ncid,Lat_ID,lat))
            call Check_err(NF90_GET_VAR(ncid,Lon_ID,lon))
         else
            call Check_err(NF90_GET_VAR(ncid,Lat_ID,lat(1,:)))
            call Check_err(NF90_GET_VAR(ncid,Lon_ID,lon(:,1)))
            do i = 2,NX
               lat(i,:) = lat(1,:)
            enddo
            do i = 2,NY
               lon(:,i) = lon(:,1)
            enddo
         endif

         IF ( .NOT. read_NWS14_NetCdf_using_core_0 ) THEN
           call Check_err(NF90_CLOSE(ncid))
         END IF

      ENDIF
      CALL BcastToLocal_Int(NX)
      CALL BcastToLocal_Int(NY)
      IF (MYPROC.NE.0) THEN
         ALLOCATE( lon(NX,NY), lat(NX,NY) )
      ENDIF
      CALL BcastToLocal_2DRealArray(lon,NX,NY)
      CALL BcastToLocal_2DRealArray(lat,NX,NY)
#endif
      CALL unsetMessageSource()
      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE READNWS14LatLon_netCDF
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE READNWS14_NC_StaTime(FileN, NC_ID)
C----------------------------------------------------------------------
C       CPB 10/2023: Sets the time index from which we start reading in
C       netCDF meteorological forcing. NOTE: reads in on proc 0 and
C       broadcasts.
C----------------------------------------------------------------------
#ifdef ADCNETCDF
      use netcdf
      use netcdf_error, only: check_err
#endif
      USE GL2LOC_MAPPING ,ONLY : BcastToLocal_Int
      USE SIZES, ONLY : MYPROC
      use mod_datetime, only: operator(==), operator(+)
      implicit none
      integer :: iret, iter
      INTEGER,INTENT(INOUT) :: NC_ID
      character(len=200),intent(in) :: fileN
      character(len=19),allocatable :: TimeStr(:)
      character(len=200) :: str_date
      integer :: i, Temp_ID, TL
      ! CPB 4/20/2023
      REAL(8),ALLOCATABLE :: NCTIMES(:)
      CALL setMessageSource('READNWS14_NC_StaTime')
#ifdef ADCNETCDF
#ifdef CMPI
      IF (MYPROC.EQ.0) THEN
#endif

         IF ( .NOT. read_NWS14_NetCdf_using_core_0 ) THEN
            call Check_err(NF90_OPEN(FileN,nf90_nowrite,NC_ID))
         END IF

         call Check_err(NF90_INQ_DIMID(NC_ID,Tdim,Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,
     &                                         len=TL))
         IF (Tformat(1:1).eq.'%') THEN
            ! in this case the netCDF files have a datetime variable
            ! so we can determine where to start based off of that
            CALL logMessage(ECHO,"fort.22 indicates that a "
     &          //"datetime variable is provided. Starting index "
     &          //"will be calculated.")
            ALLOCATE( TIMESTR(TL) )
            call Check_err(NF90_INQ_VARID(NC_ID,Tvar,Temp_ID))
            call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,TimeStr))
            ! Determine the first timestep index
            str_date = CurDT%strftime(trim(Tformat))
            do NT = 1,TL
               if (trim(str_date).eq.TimeStr(NT)) exit
            enddo
            call logMessage(ECHO,'Starting from '//TimeStr(NT))
            DEALLOCATE( TIMESTR )
         ELSEIF (Tvar.eq.'refdate') THEN
            ! in this case the reference date for the netCDF file is
            ! provided in the fort.22 and the units of the time
            ! variable are in the form of "seconds since refdate".
            ! This allows us to find the start index
            CALL logMessage(ECHO,"fort.22 provides a reference "
     &          //"date. Starting time index will be calculated.")
            ALLOCATE( NCTIMES(TL) )
            refdate = t_datetime(trim(Tformat))
            CALL CHECK_ERR(NF90_INQ_VARID(NC_ID,TDIM,TEMP_ID))
            call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,NCTIMES))
            DO NT = 1,TL
               stepdate = refdate + t_timedelta(seconds=INT(NCTIMES(NT)))
               IF (stepdate==curDT) THEN
                  EXIT
               ENDIF
            ENDDO
            DEALLOCATE( NCTIMES )
         else
            ! If neither of those two things are true we just assume
            ! we start at the beginning of the file
            call logMessage(ECHO,
     &                      'Neither a datetime variable nor a '
     &          //'reference date were provided. Assume met forcing '
     &          //'starts at the beginning of the netCDF file.')
            NT = 1
         endif
         NTC = NT

         IF ( .NOT. read_NWS14_NetCdf_using_core_0 ) THEN
            call Check_err(NF90_CLOSE(NC_ID))
         END IF

#ifdef CMPI
      endif
      CALL BcastToLocal_Int(NT)
      CALL BcastToLocal_Int(NTC)
#endif
#endif
      CALL unsetMessageSource()
      return
C----------------------------------------------------------------------
      END SUBROUTINE READNWS14_NC_StaTime
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE ArnoldSchwarzenegger(iret,var,fileN)
C----------------------------------------------------------------------
#ifdef CMPI
      USE MESSENGER, ONLY : subdomainFatalError
#endif
      IMPLICIT NONE
      integer,intent(in) :: iret
      character(len=200),intent(in) :: var, fileN
C
      if (iret.eq.0) then
         call allMessage(ERROR,'Cound not find message when reading '
     &                         //trim(var)//' in '//trim(fileN)//'.')
      elseif (iret.lt.0) then
         call allMessage(ERROR,'Fatal error when reading '
     &                         //trim(var)//' in '//trim(fileN)//'.')
      elseif (iret.gt.1) then
         call allMessage(ERROR,'Found multiple messages when reading '
     &                         //trim(var)//' in '//trim(fileN)//'.')
      endif
#ifdef CMPI
      subdomainFatalError = .true.
#endif
      call nws14Terminate()
C----------------------------------------------------------------------
      END SUBROUTINE ArnoldSchwarzenegger
C-----------------------------------------------------------------------

C---------------------------------------------------------------------
      SUBROUTINE CLOSE_NWS14_FILES()
        use global, only: NWS, NCICE
        USE SIZES, ONLY : MYPROC

#ifdef ADCNETCDF
        use netcdf
        use netcdf_error, only: check_err
#endif

        IMPLICIT NONE

#ifdef ADCNETCDF
        IF ( NWS == 14 ) THEN
          IF ( (.not. grb2flag) .AND. read_NWS14_NetCdf_using_core_0 )  THEN
             IF ( MYPROC == 0 ) THEN
               call Check_err(NF90_CLOSE(PfileNCID))
               call Check_err(NF90_CLOSE(WfileNCID))
               call Check_err(NF90_CLOSE(Wfile1NCID))

               IF (NCICE.EQ.14) THEN
                  call Check_err(NF90_CLOSE(CfileNCID))
               ENDIF
             ENDIF
          END IF
        END IF
#endif
        RETURN
      END SUBROUTINE CLOSE_NWS14_FILES
C---------------------------------------------------------------------

C----------------------------------------------------------------------
C...  jgf50.38.05: Subroutine to terminate the run cleanly.
C...
C----------------------------------------------------------------------
      SUBROUTINE nws14Terminate()
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      IMPLICIT NONE
      INTEGER,ALLOCATABLE :: dmy(:)
C
      call setMessageSource("nws14Terminate")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      call allMessage(ERROR,"ADCIRC terminating.")
#ifdef CMPI
      call msg_fini()
#endif
      CALL EXIT(1)
C

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE nws14Terminate
C----------------------------------------------------------------------

        end module mod_nws14