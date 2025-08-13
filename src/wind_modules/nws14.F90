!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2023 R.A. Luettich, Jr., J.J. Westerink
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------!
module mod_nws14
   use mod_datetime, only: t_datetime, t_timedelta
   use adc_constants, only: rad2deg, g, PRBCKGRND, rhoWat0
   use global, only: setMessageSource, unsetMessageSource, &
                     allMessage, logMessage, ERROR, WARNING, ECHO, DEBUG, &
                     scratchMessage, found_InterpInset_namelist, useInterpInset, &
                     basedatetime

   implicit none

   integer :: ub, ubc ! size of weightsp and indp interpolant weights
   real(8), allocatable :: weightsp(:, :, :)

   integer, allocatable  :: indp(:, :, :)
   character(len=200) :: Pvar, Uvar, Vvar, Cvar
   character(len=200) :: Pfile, Wfile, Pinv, Winv, Cinv, Cfile
   character(len=200) :: Wfile1
   integer :: PfileNCID, WfileNCID, CfileNCID, Wfile1NCID

   type(t_datetime) :: CurDT
   logical :: grb2flag
   integer :: NT, NTC
   real(8) :: rhoWat0g, PRBCKGRND_MH2O

   logical  :: read_NWS14_NetCdf_using_core_0 = .true. ! read NWS14 (grib2,netcdf) from a computed core 0
   logical  :: change_14 ! true when the NWS = 14 data has been updated

#ifdef ADCNETCDF
   type(t_datetime) :: refdate, stepdate
   character(len=200) :: Lonvar, Latvar, Tvar
   character(len=200) :: Tdim, Londim, Latdim, Tformat
   character(LEN=80)  :: InsetLon, InsetLat, InsetLonDim, InsetLatDim
   character(LEN=256) :: InsetPvar, InsetUVar, InsetVVar
   integer, allocatable :: inset_indp(:, :, :)
   real(8), allocatable :: inset_weight(:, :, :)
#endif

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

!----------------------------------------------------------------------
   subroutine NWS14INIT(NWS, WTIME1)
!----------------------------------------------------------------------
!.... WJP Adding in subroutines for GRB2 reading below (NWS = 14)
!.... or netcdf (depending on what is present)
!....
!.... CPB 10/2023: reorganized this subroutine to make it more readable
!.... as well as to improve input for TACC
!----------------------------------------------------------------------
      use mod_datetime, only: operator(+)
      implicit none
      integer, intent(in) :: NWS
      real(8), intent(in) :: WTIME1

      call setMessageSource('NWS14INIT')

#if !defined(ADCNETCDF) && !defined(GRIB2API)
      write (scratchMessage, '(A,I0,A)') 'Neither GRIB2API nor ADCNETCDF compiler '// &
         'flags set. Cannot use NWS = ', NWS, '. Please recompile '// &
         'with one of these options.'
      call allMessage(ERROR, scratchMessage)
      call nws14Terminate()
#endif

      ! Initialise the datetime
      CurDT = basedatetime + t_timedelta(minutes=floor(WTIME1/60d0))
      NT = 0
      NTC = 0

      rhoWat0g = rhoWat0*g
      PRBCKGRND_MH2O = 100.0d0*PRBCKGRND/rhoWat0g

      ! check if we are using grib2 or netcdf input files:
      call NWS14CHECK_FILETYPE()

      ! set filenames appropriately
      call NWS14SET_FILENAMES()

#ifdef  GRIB2API
      ! make grib2 .inv files
      call NWS14GRB2_MAKE_INV()
#endif

      ! read fort.22 to get netCDF variable info (returns if using
      ! grib2)
#ifdef ADCNETCDF
      call NWS14NC_READ_F22(NWS)
#endif

      ! Calculate the interpolant weights
      call NWS14_CALC_INTERP_WTS()

      call allMessage(ECHO, 'Finished init of NWS14')

      if (Found_InterpInset_Namelist) then
#ifdef ADCNETCDF
         call allMessage(ECHO, 'Requesting to use insets. If files found, insets will be interpolated!')
         call INIT_INSET_NC()
         call CALC_INSET_WTS()
#else
         call allMessage(ERROR, "Inset files in netCDF format requested, "// &
                         "but ADCIRC not compiled with netCDF support."// &
                         "Please recompile with netCDF support.")
         call nws14Terminate()
#endif
      end if

      call unsetMessageSource()
!----------------------------------------------------------------------
   end subroutine NWS14INIT
!----------------------------------------------------------------------

!----------------------------------------------------------------------
   subroutine NWS14GET(WTIMINC, WVNX2, WVNY2, PRN2, CICE2)
!----------------------------------------------------------------------
      use mesh, only: NP
      use mod_datetime, only: operator(+)
      implicit none

      real(8), intent(in) :: WTIMINC
      real(8), intent(inout) :: PRN2(NP), WVNX2(NP), WVNY2(NP)
      real(8), intent(inout), optional :: CICE2(NP)

#if defined(ADCNETCDF) || defined(GRIB2API)

      integer :: i, NTkeep
      real, allocatable, dimension(:, :) :: Pmsl, U10, V10, Cice

      call setMessageSource('NWS14GET')
      ! Show the date in str_date format
      call logMessage(ECHO, 'CurDT strng: '// &
                      CurDT%strftime("%Y-%m-%d %H:%M"))

      ! Get the Pressure
      !':PRMSL:mean sea level:'
      call READNWS14(Pfile, Pmsl, Pinv, Pvar, PfileNCID)

      ! Get the U-winds
      call READNWS14(Wfile, U10, Winv, Uvar, WfileNCID)

      ! Get the V-winds
      call READNWS14(Wfile1, V10, Winv, Vvar, Wfile1NCID)

      ! Get the ice concentration if required
      if (present(CICE2)) then
         NTkeep = NT
         NT = NTC
         call READNWS14(Cfile, Cice, Cinv, Cvar, CfileNCID)
         NT = NTkeep
      end if

      ! Doing the interpolation from their grid to ours
      do i = 1, np
         if (indp(1, 1, i) < 0) then
            PRN2(i) = PRBCKGRND_MH2O
         else
            PRN2(i) = (dble(Pmsl(indp(1, 1, i), indp(1, 2, i)))*weightsp(1, 1, i) + &
                       dble(Pmsl(indp(1, 3, i), indp(1, 2, i)))*weightsp(1, 2, i) + &
                       dble(Pmsl(indp(1, 1, i), indp(1, 4, i)))*weightsp(1, 3, i) + &
                       dble(Pmsl(indp(1, 3, i), indp(1, 4, i)))*weightsp(1, 4, i)) &
                      /rhoWat0g
         end if
         if (indp(ub, 1, i) < 0) then
            WVNX2(i) = 0d0
            WVNY2(i) = 0d0
         else
            WVNX2(i) = dble(U10(indp(ub, 1, i), indp(ub, 2, i)))*weightsp(ub, 1, i) + &
                       dble(U10(indp(ub, 3, i), indp(ub, 2, i)))*weightsp(ub, 2, i) + &
                       dble(U10(indp(ub, 1, i), indp(ub, 4, i)))*weightsp(ub, 3, i) + &
                       dble(U10(indp(ub, 3, i), indp(ub, 4, i)))*weightsp(ub, 4, i)
            WVNY2(i) = dble(V10(indp(ub, 1, i), indp(ub, 2, i)))*weightsp(ub, 1, i) + &
                       dble(V10(indp(ub, 3, i), indp(ub, 2, i)))*weightsp(ub, 2, i) + &
                       dble(V10(indp(ub, 1, i), indp(ub, 4, i)))*weightsp(ub, 3, i) + &
                       dble(V10(indp(ub, 3, i), indp(ub, 4, i)))*weightsp(ub, 4, i)
         end if
      end do
      deallocate (Pmsl, U10, V10)

      if (found_InterpInset_Namelist .and. UseInterpInset) then
#ifdef ADCNETCDF
         ! Insert insets in netcdf format  (Aman Tejaswi)
         call INSET_NC_INPUT(PRN2, WVNX2, WVNY2)
#else
         call allMessage(ERROR, "Inset files in netCDF format requested, "// &
                         "but ADCIRC not compiled with netCDF support."// &
                         "Please recompile with netCDF support.")
         call nws14Terminate()
#endif
      else if (found_InterpInset_Namelist .and. (.not. UseInterpInset)) then
         call logMessage(ECHO, 'UseInterpInset set to F, Please Check Namelist if you want to use inset files')
      end if

      ! Add WTIMINC on CurDT for next WTIME
      CurDT = CurDT + t_timedelta(minutes=nint(WTIMINC/60d0))
      ! Next time index for the netcdf reading
      NT = NT + 1

      ! If ice concentration is present
      if (present(CICE2)) then
         do i = 1, np
            if (indp(ubc, 1, i) < 0) then
               CICE2(i) = 0d0; 
            else
               CICE2(i) = &
                  dble(Cice(indp(ubc, 1, i), indp(ubc, 2, i)))*weightsp(ubc, 1, i) + &
                  dble(Cice(indp(ubc, 3, i), indp(ubc, 2, i)))*weightsp(ubc, 2, i) + &
                  dble(Cice(indp(ubc, 1, i), indp(ubc, 4, i)))*weightsp(ubc, 3, i) + &
                  dble(Cice(indp(ubc, 3, i), indp(ubc, 4, i)))*weightsp(ubc, 4, i)
            end if
         end do
         deallocate (Cice)
         ! Next ice time index for the netcdf reading
         NTC = NTC + 1
      end if

      call unsetMessageSource()

#else
      ! Force compiler not to warn about unused variables
      write (scratchMessage, '(I0)') WTIMINC
      WVNX2 = WVNX2
      WVNY2 = WVNY2
      PRN2 = PRN2
      if (present(CICE2)) then
         CICE2 = CICE2
      end if

      call allMessage(ERROR, 'Neither GRIB2API nor ADCNETCDF compiler '// &
                      'flags set. Cannot use NWS = 14. Please recompile '// &
                      'with one of these options.')
      call nws14Terminate()
#endif

   end subroutine NWS14GET
!----------------------------------------------------------------------

   !----------------------------------------------------------------------
   subroutine NWS14CHECK_FILETYPE()
      !----------------------------------------------------------------------
      !     CPB 10/2023: This subroutine simply checks whether we are using
      !     grib2 or netCDF winds in NWS14 and sets the grb2flag logical
      !     appropriately. It does this by:
      !        1. Checking if the pressure grib2 file (fort.221.grb2) exists
      !        2. If it does not exist or if we are compiled without grib2
      !           then it checks if the pressure netCDF file (fort.221.nc)
      !           exists.
      !        3. If neither exist we abort with an error message.
      !----------------------------------------------------------------------
      implicit none
      logical :: Fexists

      call setMessageSource("NWS14CHECK_FILETYPE")

      Pfile = 'fort.221.grb2'
      inquire (file=Pfile, exist=Fexists)
#ifdef GRIB2API
      if (Fexists) then
         grb2flag = .true.
         call unsetMessageSource()
         return
      end if
#endif
      Pfile = 'fort.221.nc'
      inquire (file=Pfile, exist=Fexists)
      if (Fexists) then
         grb2flag = .false.
         call unsetMessageSource()
         return
      else
         call allMessage(ERROR, &
                         'Neither .grb2 nor .nc wind files exist. Or, if .grb2 '// &
                         'files exist, did you compile with GRIB2 compiler flags?')
         call nws14Terminate()
      end if

      call unsetMessageSource()

   end subroutine NWS14CHECK_FILETYPE
   !----------------------------------------------------------------------

   !----------------------------------------------------------------------
   subroutine NWS14SET_FILENAMES()
      !----------------------------------------------------------------------
      !     CPB 10/2023: This subroutine sets the met forcing filenames based
      !     on grb2flag as well as the existence or lack thereof of various
      !     files. It also opens the netcdf files in read only mode to
      !     eliminate opening/closing the files repeatedly
      !----------------------------------------------------------------------
#ifdef ADCNETCDF
      use netcdf, only: NF90_OPEN, NF90_NOWRITE, NF90_CLOSE
      use netcdf_error, only: check_err
      use SIZES, only: MYPROC
#endif
#ifdef CMPI
      use MESSENGER, only: MSG_BARRIER
#endif

      use GLOBAL, only: NCICE
      implicit none
      logical :: Fexists

      call setMessageSource("NWS14SET_FILENAMES")

      if (grb2flag) then
         Pfile = 'fort.221.grb2'
         Pinv = 'fort.221.inv'
         Wfile = 'fort.222.grb2'
         Wfile1 = 'fort.222.grb2'
         Winv = 'fort.222.inv'
         Cfile = 'fort.225.grb2'
         Cinv = 'fort.225.inv'

         Pvar = 'S' ! searching for PRMSL or MSLET or PRES:SURFACE
         ! prev: :PRMSL:mean sea level:'
         Uvar = ':UGRD:10 m above ground:'
         Vvar = ':VGRD:10 m above ground:'
         Cvar = ':ICEC:surface:'
         ! Assign dummy values to ncid variables
         PfileNCID = -9999
         WfileNCID = -9999
         Wfile1NCID = -9999
         if (NCICE == 14) then
            CfileNCID = -9999
         end if
      elseif (.not. grb2flag) then
         Pfile = 'fort.221.nc'
         Wfile = 'fort.222.nc'
         Wfile1 = 'fort.222.nc'
         Cfile = 'fort.225.nc'
         ! check if we are using two different wind files
         inquire (FILE=Wfile, exist=Fexists)
         if (.not. Fexists) then
            ! using 2 wind files
            Wfile = 'fort.222u.nc'
            Wfile1 = 'fort.222v.nc'
         end if
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
         if (read_NWS14_NetCdf_using_core_0) then
            if (MYPROC == 0) then
               call Check_err(NF90_OPEN(Pfile, nf90_nowrite, PfileNCID))
               call Check_err(NF90_OPEN(Wfile, nf90_nowrite, WfileNCID))
               call Check_err(NF90_OPEN(Wfile1, nf90_nowrite, Wfile1NCID))
               if (NCICE == 14) then
                  call Check_err(NF90_OPEN(Cfile, nf90_nowrite, CfileNCID))
               end if
            end if

#ifdef CMPI
            ! wait til proc 0 has opened all the files
            call MSG_BARRIER()
#endif
         end if

#endif

      else
         ! should be unreachable
         call nws14Terminate()
      end if

      call unsetMessageSource()
   end subroutine NWS14SET_FILENAMES
   !----------------------------------------------------------------------

#ifdef GRIB2API
!----------------------------------------------------------------------
   subroutine NWS14GRB2_MAKE_INV()
!----------------------------------------------------------------------
!     CPB 10/2023: This subroutine makes the .inv files for grib2
!     meteorological forcing files. If we are not using grib2 then it
!     just returns
!----------------------------------------------------------------------
      use wgrib2api, only: grb2_mk_inv
#ifdef CMPI
      use MESSENGER, only: MSG_BARRIER
#endif
      use SIZES, only: MYPROC
      use GL2LOC_MAPPING, only: BcastToLocal_Int
      use GLOBAL, only: NCICE
      implicit none
      integer :: iret
      logical :: Fexists

      call setMessageSource("NWS14GRB2_MAKE_INV")

      if (.not. grb2flag) then
         call unsetMessageSource()
         return
      end if

      iret = 0

      ! only do this on proc 0
      if (MYPROC == 0) then
         ! pressure file (fort.221.grb2)
         inquire (file=Pinv, exist=Fexists)
         if (.not. Fexists) then
            iret = grb2_mk_inv(Pfile, Pinv)
            if (iret /= 0) then
               call allMessage(ERROR, 'Fatal error in reading fort.221.grb2.')
            else
               call allMessage(ECHO, 'successfully read fort.221.grb2 and wrote out fort.221.inv file')
            end if
         end if
      end if
      call BcastToLocal_Int(iret)
      if (iret /= 0) then
         call allMessage(ERROR, 'Error generating grib2 pressure inventory')
         call nws14Terminate()
      end if

      if (MYPROC == 0) then
         ! wind file (fort.222.grb2)
         inquire (file=Winv, exist=Fexists)
         if (.not. Fexists .and. iret == 0) then
            iret = grb2_mk_inv(Wfile, Winv)
            if (iret /= 0) then
               call allMessage(ERROR, 'Fatal error in reading fort.222.grb2.')
            else
               call allMessage(ECHO, 'successfully read fort.222.grb2 and wrote out fort.222.inv file')
            end if
         end if
      end if
      call BcastToLocal_Int(iret)
      if (iret /= 0) then
         call allMessage(ERROR, 'Error generating grib2 wind inventory')
         call nws14Terminate()
      end if

      if (MYPROC == 0) then
         ! ice file (if we are using it) (fort.225.grb2)
         if (NCICE == 14) then
            inquire (file=Cinv, exist=Fexists)
            if (.not. Fexists .and. iret == 0) then
               iret = grb2_mk_inv(Cfile, Cinv)
               if (iret /= 0) then
                  call allMessage(ERROR, 'Fatal error in reading fort.225.grb2.')
               else
                  call allMessage(ECHO, 'successfully read fort.225.grb2 and wrote out fort.225.inv file')
               end if
            end if
         end if
      end if
      call BcastToLocal_Int(iret)
      if (iret /= 0) then
         call allMessage(ERROR, 'Error generating grib2 ice inventory')
         call nws14Terminate()
      end if

      call unsetMessageSource()
   end subroutine NWS14GRB2_MAKE_INV
!----------------------------------------------------------------------
#endif

   !----------------------------------------------------------------------
#ifdef ADCNETCDF
   subroutine NWS14NC_READ_F22(NWS)
!----------------------------------------------------------------------
!     CPB 10/2023: This subroutine reads the fort.22 file for NWS14 to
!     get netcdf variable names as well as starting indices for reading
!     in the met data. If we are using grib2 it just returns
!----------------------------------------------------------------------
      use GLOBAL, only: NCICE
      implicit none
      integer, intent(in) :: NWS
      integer :: ierr, PvarInd

      call setMessageSource("NWS14NC_READ_F22")

      ! return if we are using grib2
      if (grb2flag) then
         call unsetMessageSource()
         return
      end if

      if (NWS == -14) then
         ! For inset OWI winds with netcdf files something needs to be
         ! done but I haven't programmed it yet.
         call allMessage(ERROR, "NWS = -14 has not yet been set up to " &
                         //"use netcdf files. Try again with grib2 " &
                         //"background wind files.")
         call nws14Terminate()
      end if

      ! open the fort.22
      open (22, file='fort.22', ACTION='READ', IOSTAT=ierr)

      ! if we had trouble opening then abort
      if (ierr /= 0) then
         call allMessage(ERROR, 'Unable to open fort.22.')
         call nws14Terminate()
      end if

      read (22, *) Tdim ! Time dimension name
      read (22, *) Tvar ! Time variable name
      read (22, *) Tformat ! Format of time datestr
      read (22, *) Londim ! Lon dimension name
      read (22, *) Lonvar ! Lon var name
      read (22, *) Latdim ! Lat dimension name
      read (22, *) Latvar ! Lat var name
      read (22, *) Pvar ! Pressure var name
      read (22, *) Uvar ! U10 var name
      read (22, *) Vvar ! V10 Var name

      if (NCICE == 14) then
         read (22, *) Cvar ! Ice var name
      end if
      close (22)

      PvarInd = index(Pvar, 'HPa')
      if (PvarInd > 0) then
         rhoWat0g = rhoWat0g/100d0
         Pvar = Pvar(1:PvarInd - 1)
         write (16, *) 'Pressure is in Hpa, New Pvar = ', trim(Pvar)
      end if

      ! get starting time index for netcdf file
      call READNWS14_NC_StaTime(Pfile, PfileNCID)
      call unsetMessageSource()
   end subroutine NWS14NC_READ_F22
#endif
   !----------------------------------------------------------------------

   !----------------------------------------------------------------------
   subroutine NWS14_CALC_INTERP_WTS()
      !----------------------------------------------------------------------
      !     CPB 10/2023: This subroutine calculates the interpolation indices
      !     (indp) and weights (weightsp) to interpolate from the grib2 or
      !     netcdf input meterological forcing to the ADCIRC grid. NWS = 14
      !     was originally set up to use NOAA CFSv2/GFS met forcing products.
      !     These meteorological products provide two different resolutions
      !     for wind forcing but not for pressure (or ice). To limit
      !     pre-processing necessary to use these products, this subroutine
      !     checks to see if the wind and pressure have the same dimensions.
      !     If so, it only allocates one set of interpolant indices (indp) and
      !     weights (weightsp). If they are different then it stores two
      !     different sets. If necessary, it also checks the dimensions of the
      !     ice forcing and sets the interpolants for ice to be equal to
      !     either the wind or pressure forcing depending on the grid.
      !----------------------------------------------------------------------
      use mesh, only: NP, SLAM, SFEA, bl_interp, bl_interp2
      use kdtree2_module, only: kdtree2, kdtree2_create
      use GLOBAL, only: NCICE
      implicit none
      integer :: i, ii, nxp, nyp, indt(4), xi, yi, cc
      real, allocatable, dimension(:, :) :: lat, lon
      real(8), allocatable ::  lonv(:), latv(:), XY(:, :)
      real(8) :: xx, yy, wt(4)
      logical :: regular_grid
      type(kdtree2), pointer :: tree

      ! determine if wind/pressure/ice grids match or not.
      ! assume grids are the same to start
      ub = 1
      ubc = 1 ! index where we store the ice data (will match either
      ! wind or pressure). Assume match pressure first.
      ! get pressure dimensions
      call READNWS14LatLon(Pfile, lat, lon, Pinv, Pvar, PfileNCID)
      nxp = ubound(lon, 1); nyp = ubound(lat, 2)
      ! get wind dimensions and see if they match pressure
      call READNWS14LatLon(Wfile, lat, lon, Winv, Uvar, WfileNCID)
      if (nxp /= ubound(lon, 1) .or. nyp /= ubound(lat, 2)) then
         ! if we don't match both lat and lon then we need to store two
         ! sets of interpolation indices/weights
         ub = 2
         nxp = ubound(lon, 1)
         nyp = ubound(lat, 2)
      end if
      ! If winds and pressure do not match, then check ice file to see
      ! whether we should use wind or pressure interpolants. If wind and
      ! pressure DO match then we assume ice also matches.
      if (NCICE == 14 .and. ub == 2) then
         call READNWS14LatLon(Cfile, lat, lon, Cinv, Cvar, CfileNCID)
         if (nxp == ubound(lon, 1) .and. nyp == ubound(lat, 2)) then
            ! match winds
            ubc = 2
         end if
      end if
      ! Allocate the indices and weights
      allocate (indp(ub, 4, np), weightsp(ub, 4, np))
      ! Make the interpolant weights. If w
      do ii = 1, ub
         if (ii == 1) then
            call READNWS14LatLon(Pfile, lat, lon, Pinv, Pvar, PfileNCID)
         else
            call READNWS14LatLon(Wfile, lat, lon, Winv, Uvar, WfileNCID)
         end if
         nxp = ubound(lon, 1); nyp = ubound(lat, 2)
         ! Check if structured regular..
         if (dble(abs(lon(1, 1) - lon(1, nyp))) < 1d-6) then
            regular_grid = .true.
            allocate (lonv(nxp), latv(nyp))
            lonv = dble(lon(:, 1))
            latv = dble(lat(1, :))
         else
            regular_grid = .false.
            ! Convert point matrices to 2 X (NXP*NYP) point vector
            allocate (XY(2, nxp*nyp))
            cc = 0
            do xi = 1, nxp
               do yi = 1, nyp
                  cc = cc + 1
                  XY(1, cc) = dble(lon(xi, yi))
                  XY(2, cc) = dble(lat(xi, yi))
               end do
            end do
            ! Create the search tree
            tree => kdtree2_create(XY)
         end if
         ! Loop over each node and interpolate
         do i = 1, np
            xx = rad2deg*slam(i)
            ! Convert our numbers if grids are 0 to 360
            if (dble(maxval(lon)) > 180d0 .and. xx < 0d0) xx = xx + 360d0
            yy = rad2deg*sfea(i)
            if (regular_grid) then
               call bl_interp(nxp, lonv, nyp, latv, xx, yy, indt, wt)
            else
               call bl_interp2(nxp, nyp, lon, lat, xx, yy, indt, wt, tree)
            end if
            indp(ii, :, i) = indt
            weightsp(ii, :, i) = wt
         end do
         if (regular_grid) then
            deallocate (lon, lat, latv, lonv)
         else
            deallocate (lon, lat, XY)
         end if
      end do
   end subroutine NWS14_CALC_INTERP_WTS
!----------------------------------------------------------------------

!----------------------------------------------------------------------
#if defined(GRIB2API) || defined(ADCNETCDF)
   subroutine READNWS14(fileN, data2, invN, var, ncid)
!----------------------------------------------------------------------
!     CPB 10/2023: Added to make NWS = 14 more more readable. Simply
!     calls the appropriate subroutine based on the filetype we are
!     using.
!----------------------------------------------------------------------
      implicit none
      character(LEN=200), intent(in) :: fileN
      real, allocatable, intent(inout) :: data2(:, :)
      character(LEN=200), intent(in) :: invN, var
      integer, intent(inout) :: ncid

      if (.not. grb2flag) then
#ifdef ADCNETCDF
         call READNWS14_netCDF(fileN, ncid, var, data2)
#else
         write (scratchMessage, '(5a)') 'netCDF support not compiled inADCIRC.', &
            ' Cannot read NWS=14 lat/lon from netCDF file: ', trim(fileN), &
            ' Please recompile with ADCNETCDF option. var: ', trim(var)
#endif
      else
#ifdef GRIB2API
         call READNWS14_grib2(fileN, invN, var, data2)
#else
         write (scratchMessage, '(7a)') 'GRIB2 support not compiled in ADCIRC.', &
            ' Cannot read NWS=14 lat/lon from GRIB2 file: ', trim(fileN), &
            ' Please recompile with GRIB2API option. invN: ', trim(invN), ' var: ', trim(var)
         call allMessage(ERROR, scratchMessage)
         call nws14Terminate()
#endif
      end if
!----------------------------------------------------------------------
   end subroutine READNWS14
!----------------------------------------------------------------------
#endif

#ifdef GRIB2API
!----------------------------------------------------------------------
   subroutine READNWS14_grib2(fileN, invN, var, data2)
!----------------------------------------------------------------------
!       CPB 10/2023: Reads in grib2 format met forcing for NWS = 14.
!       NOTE: reads in on Proc 0 and broadcasts.
!----------------------------------------------------------------------
      use wgrib2api, only: grb2_inq
      use GL2LOC_MAPPING, only: BcastToLocal_2DRealArray, &
                                BcastToLocal_Int
      use SIZES, only: MYPROC
      implicit none
      integer :: iret, iter
      character(len=200), intent(in) :: var, fileN, invN
      real, allocatable, intent(out) :: data2(:, :)
      character(len=200) :: str_date
      integer :: NX, NY
      call setMessageSource("READNWS14_grib2")
      if (MYPROC == 0) then
         ! Getting the date in str_date format
         str_date = ':start_FT='//trim(CurDT%strftime("%Y%m%d%H"))//'0000:'
         iret = -1; iter = 0
         do while (iret <= 0 .and. iter < 10)
            if (iter > 0) then
               call logMessage(WARNING, 'Trying to read again '//var)
               call sleep(5)
            end if
            iret = grb2_inq(fileN, invN, var, str_date, data2=data2)
            if (iret > 1 .and. iter == 0) then
               ! May have two entries bcause of forecast/DA overlap
               call logMessage(DEBUG, '> 1 msg, trying to read '//var)
               iret = grb2_inq(fileN, invN, var, ':anl:', &
                               str_date, data2=data2)
            end if
            iter = iter + 1
         end do
         if (iret > 0) then
            call logMessage(DEBUG, 'Successfully read '//var)
         else
            call ArnoldSchwarzenegger(iret, var, fileN)
         end if
         NX = ubound(data2, 1)
         NY = ubound(data2, 2)
      end if
      call BcastToLocal_Int(NX)
      call BcastToLocal_INT(NY)
      if (MYPROC /= 0) then
         allocate (data2(NX, NY))
      end if
      call BcastToLocal_2DRealArray(data2, NX, NY)
      call unsetMessageSource()
      return
!----------------------------------------------------------------------
   end subroutine READNWS14_grib2
!----------------------------------------------------------------------
#endif

!----------------------------------------------------------------------
#ifdef ADCNETCDF
   subroutine READNWS14_netCDF(fileN, NC_ID, var, data2)
!----------------------------------------------------------------------
!       CPB 10/2023: Reads in netCDF format met forcing for NWS = 14.
!       NOTE: reads in on Proc 0 and broadcasts.
!----------------------------------------------------------------------

      use netcdf, only: NF90_OPEN, NF90_NOWRITE, NF90_CLOSE, &
                        NF90_INQ_DIMID, NF90_INQUIRE_DIMENSION, &
                        NF90_INQ_VARID, NF90_GET_VAR
      use netcdf_error, only: check_err
      use GL2LOC_MAPPING, only: BcastToLocal_2DRealArray, &
                                BcastToLocal_Int
      use SIZES, only: MYPROC
      implicit none

      character(len=*), intent(in):: fileN
      character(len=200), intent(in) :: var
      real(4), allocatable, intent(inout) :: data2(:, :)
      integer, intent(inout) :: NC_ID
      integer :: Temp_ID, NX, NY

      call setMessageSource('READNWS14_netCDF')

      if (read_NWS14_NetCdf_using_core_0) then
         ! Use core zero to read data
         ! and then broadcast them
         if (MYPROC == 0) then
            call Check_err(NF90_INQ_DIMID(NC_ID, Latdim, Temp_ID))
            call Check_err(NF90_INQUIRE_DIMENSION(NC_ID, Temp_ID, len=NY))
            call Check_err(NF90_INQ_DIMID(NC_ID, Londim, Temp_ID))
            call Check_err(NF90_INQUIRE_DIMENSION(NC_ID, Temp_ID, len=NX))
         end if
         call BcastToLocal_Int(NX)
         call BcastToLocal_Int(NY)
         allocate (data2(NX, NY))
         if (MYPROC == 0) then
            call Check_err(NF90_INQ_VARID(NC_ID, var, Temp_ID))
            call Check_err(NF90_GET_VAR(NC_ID, Temp_ID, data2, &
                                        start=[1, 1, NT], count=[NX, NY, 1]))
         end if
         call BCastToLocal_2DRealArray(Data2, NX, NY)
      else
         ! Each core read data
         call Check_err(NF90_OPEN(fileN, nf90_nowrite, NC_ID))
         call Check_err(NF90_INQ_DIMID(NC_ID, Latdim, Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(NC_ID, Temp_ID, len=NY))
         call Check_err(NF90_INQ_DIMID(NC_ID, Londim, Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(NC_ID, Temp_ID, len=NX))

         allocate (data2(NX, NY))
         call Check_err(NF90_INQ_VARID(NC_ID, var, Temp_ID))
         call Check_err(NF90_GET_VAR(NC_ID, Temp_ID, data2, &
                                     start=[1, 1, NT], count=[NX, NY, 1]))

         call Check_err(NF90_CLOSE(NC_ID))
      end if

      call unsetMessageSource()

!----------------------------------------------------------------------
   end subroutine READNWS14_netCDF
!----------------------------------------------------------------------
#endif

#ifdef ADCNETCDF
!---------------------------------------------------------------------
   subroutine INSET_NC_INPUT(PRN2, WVNX2, WVNY2)
!----------------------------------------------------------------------
!     Aman Tejaswi (11th Feb 2025)
!     This subroutine uses high res inset netcdf wind and pressure data.
!     The way it is implemented that first the global background wind (ex. GFS,
!     ERA5) and pressure are interpolated on the ADCIRC mesh and then
!     wherever highres wind is available, it interpolates the inset data on the
!     Finite Element grid. The file 'inset_netcdf.in' is used to get the
!     lat lon variable names as well as variable names.
!     Expected Structure of The File Is Like This: (example)
!
!-------------------------------------
!              lonDimName
!              latDimName
!              lonVarName
!              latVarName
!              InsetPresureVarName
!              InsetUWindVarName
!              InsetVWindVarName
!-------------------------------------
!     For now the code assumes that your background wind and your high res
!     inset has same temporal coverage.
!     The code only works with NWS=14 (netcdf). The implementation is controlled
!     by Namelist -- &InterpInset, UseInterpInset=T. The
!     subroutine can be used in the sameway as OWI_NETCDF with and
!     advantage that one doesn't need to write global and regional
!     wind/pressure data in a particular (OWI) format as sometimes it
!     becomes tedious if running decadal simulations.
!----------------------------------------------------------------------
      use GLOBAL, only: PFile_inset, Wfile_inset

      implicit none

      real(8), intent(inout) :: PRN2(:), WVNX2(:), WVNY2(:)

      integer :: TempNCID_P, TempNCID_W

      real(8), allocatable, dimension(:, :) :: InsetPData
      real(8), allocatable, dimension(:, :) :: InsetUData, InsetVData

      call setMessageSource("INSET_NC_INPUT")
      Pfile_inset = trim(adjustl(Pfile_inset))
      Wfile_inset = trim(adjustl(Wfile_inset))

      call READINSET_NC(Pfile_inset, TempNCID_P, InsetLatDim, InsetLonDim, InsetPvar, InsetPData)
      call READINSET_NC(Wfile_inset, TempNCID_W, InsetLatDim, InsetLonDim, InsetUVar, InsetUData)
      call READINSET_NC(Wfile_inset, TempNCID_W, InsetLatDim, InsetLonDim, InsetVvar, InsetVData)

      call INTERP_INSET_GRID(InsetPData, InsetUData, InsetVData, PRN2, WVNX2, WVNY2)

      call unsetMessageSource()

!----------------------------------------------------------------------
   end subroutine INSET_NC_INPUT
#endif

!----------------------------------------------------------------------
#ifdef ADCNETCDF
   subroutine INIT_INSET_NC()

      use SIZES, only: GBLINPUTDIR
      use GLOBAL, only: InsetControlFile

      implicit none

      logical :: InsetFileExist
      integer :: IERR

      call setMessageSource("INIT_INSET_NC")

      inquire (FILE=trim(GBLINPUTDIR)//'/'//trim(adjustl(InsetControlFile)), EXIST=InsetFileExist)
      if (.not. InsetFileExist) then
         write (scratchMessage, '(A)') "Unable to find InsetControlFile: "//trim(InsetControlFile)
         call allMessage(ERROR, scratchMessage)
         call nws14Terminate()
      end if

      ! open the inset_netcdf.in file
      open (2200, FILE=trim(GBLINPUTDIR)//'/'//trim(adjustl(InsetControlFile)), ACTION='READ', IOSTAT=IERR)

      ! if we had trouble opening then abort
      if (IERR /= 0) then
         write (scratchMessage, '(A)') 'Unable to open insetControlFile: '//trim(InsetControlFile)
         call allMessage(ERROR, scratchMessage)
         call nws14Terminate()
      end if

      read (2200, '(A)', IOSTAT=IERR) InsetLonDim
      read (2200, '(A)', IOSTAT=IERR) InsetLatDim
      read (2200, '(A)', IOSTAT=IERR) InsetLon
      read (2200, '(A)', IOSTAT=IERR) InsetLat
      read (2200, '(A)', IOSTAT=IERR) InsetPvar
      read (2200, '(A)', IOSTAT=IERR) InsetUVar
      read (2200, '(A)', IOSTAT=IERR) InsetVVar

      close (2200) !closing the file here..

   end subroutine INIT_INSET_NC
#endif
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!     Aman Tejaswi (11 Feb 2025)
!     Read High Res-Inset Winds in NETCDF format and Interpolate onto
!     ADCIRC grid. These insets can be used particularly with a
!     background global wind (e.g GFS) which is read also in NETCDF. The
!     advantage of using this is specific to GLOBAL MODEL (GSTOFS).
!----------------------------------------------------------------------

#ifdef ADCNETCDF
!----------------------------------------------------------------------
   subroutine READINSET_NC(InsetFileN, InsetID, InsetLatDim, InsetLonDim, InsetVar, InsetOUT)

      use netcdf, only: NF90_OPEN, NF90_NOWRITE, NF90_CLOSE, &
                        NF90_INQ_DIMID, NF90_INQUIRE_DIMENSION, &
                        NF90_INQ_VARID, NF90_GET_VAR
      use netcdf_error, only: check_err

      use GL2LOC_MAPPING, only: BcastToLocal_2DRealArray, &
                                BcastToLocal_Int
      use GLOBAL, only: scratchMessage, allMessage, ERROR

      implicit none

      character(LEN=*), intent(in)      :: InsetFileN
      character(LEN=*), intent(in)      :: InsetVar, InsetLatDim, InsetLonDim
      integer, intent(out)              :: InsetID
      real(8), allocatable, intent(out) :: InsetOUT(:, :)
      integer                           :: InsetLoc, INSET_NX, INSET_NY, VarLoc, alloc_err

      call setMessageSource('READINSET_NC')

      call Check_err(NF90_OPEN(InsetFileN, nf90_nowrite, InsetID))

      call Check_err(NF90_INQ_DIMID(InsetID, InsetLonDim, InsetLoc))
      call Check_err(NF90_INQUIRE_DIMENSION(InsetID, InsetLoc, len=INSET_NX))

      call Check_err(NF90_INQ_DIMID(InsetID, InsetLatDim, InsetLoc))
      call Check_err(NF90_INQUIRE_DIMENSION(InsetID, InsetLoc, len=INSET_NY))

      allocate (InsetOUT(INSET_NX, INSET_NY), STAT=alloc_err)

      if (alloc_err /= 0) then
         write (scratchMessage, '(A)') "Unable to allocate memory for InsetOUT"
         call allMessage(ERROR, scratchMessage)
         call nws14Terminate()
      end if

      call Check_err(NF90_INQ_VARID(InsetID, InsetVar, VarLoc))
      call Check_err(NF90_GET_VAR(InsetID, VarLoc, InsetOUT, &
                                  start=[1, 1, NT], count=[Inset_NX, Inset_NY, 1]))

      call Check_err(NF90_CLOSE(InsetID))

      call unsetMessageSource()

   end subroutine READINSET_NC
#endif

!--------------------------------------------------------------------------
#ifdef ADCNETCDF
   subroutine CALC_INSET_WTS()

      use mesh, only: NP, SLAM, SFEA, bl_interp, bl_interp2
      use kdtree2_module, only: kdtree2, kdtree2_create, kdtree2_destroy
      use GLOBAL, only: Pfile_inset
      implicit none

      integer :: i, ii, nxp, nyp, indt(4), xi, yi, cc, ncid_p
      real, allocatable, dimension(:, :) :: lat, lon
      real(8), allocatable ::  lonv(:), latv(:), XY(:, :)
      real(8) :: xx, yy, wt(4)
      logical :: regular_grid
      type(kdtree2), pointer :: tree

      ! assume grids are the same to start
      ub = 1

      call READ_INSET_LatLon(Pfile_inset, NCID_p, lat, InsetLatDim, InsetLat, lon, InsetLonDim, InsetLon, nxp, nyp)
      nxp = ubound(lon, 1); nyp = ubound(lat, 2)

      if (nxp /= ubound(lon, 1) .or. nyp /= ubound(lat, 2)) then
         ! if we don't match both lat and lon then we need to store two
         ! sets of interpolation indices/weights
         ub = 2
         nxp = ubound(lon, 1)
         nyp = ubound(lat, 2)
      end if

      if (.not. allocated(inset_indp) .and. .not. allocated(inset_weight)) then
         ! Allocate the indices and weights
         allocate (inset_indp(ub, 4, np), inset_weight(ub, 4, np))
      end if

      ! Make the interpolant weights. If w
      do ii = 1, ub

         ! Check if structured regular..
         if (dble(abs(lon(1, 1) - lon(1, nyp))) < 1d-6) then
            regular_grid = .true.
            if (.not. allocated(lonv) .and. .not. allocated(latv)) then
               allocate (lonv(nxp), latv(nyp))
            end if
            lonv = dble(lon(:, 1))
            latv = dble(lat(1, :))

         else
            regular_grid = .false.
            ! Convert point matrices to 2 X (NXP*NYP) point vector
            if (.not. allocated(xy)) then
               allocate (XY(2, nxp*nyp))
            end if
            cc = 0
            do xi = 1, nxp
               do yi = 1, nyp
                  cc = cc + 1
                  XY(1, cc) = dble(lon(xi, yi))
                  XY(2, cc) = dble(lat(xi, yi))
               end do
            end do
            ! Create the search tree
            tree => kdtree2_create(XY)
         end if
         ! Loop over each node and interpolate
         do i = 1, np
            xx = rad2deg*slam(i)
            ! Convert our numbers if grids are 0 to 360
            if (dble(maxval(lon)) > 180d0 .and. xx < 0d0) xx = xx + 360d0
            yy = rad2deg*sfea(i)
            if (regular_grid) then
               call bl_interp(nxp, lonv, nyp, latv, xx, yy, indt, wt)
            else
               call bl_interp2(nxp, nyp, lon, lat, xx, yy, indt, wt, tree)
            end if
            inset_indp(ii, :, i) = indt
            inset_weight(ii, :, i) = wt
         end do
         if (regular_grid) then
            deallocate (lon, lat, latv, lonv)
         else
            deallocate (lon, lat, XY)
         end if

         if (associated(tree)) then
            call kdtree2_destroy(tree)
            nullify (tree)
         end if

      end do

   end subroutine CALC_inset_wts
#endif

!------------------------------------------------------------------------
#ifdef ADCNETCDF
   subroutine INTERP_INSET_GRID(InsetPData, InsetUData, InsetVData, PRN2, WVNX2, WVNY2)
      use MESH, only: NP
      implicit none

      integer :: i
      real(8), intent(in) :: InsetPData(:, :), InsetUData(:, :), InsetVData(:, :)
      real(8), intent(inout) :: PRN2(:), WVNX2(:), WVNY2(:)

      ! Doing the interpolation from inset grid to finite ele grid
      do i = 1, np
         if (inset_indp(1, 1, i) < 0) then
            PRN2(i) = PRN2(i) !Keeping the background global pressure
         else

            PRN2(i) = (InsetPData(inset_indp(1, 1, i), inset_indp(1, 2, i))*inset_weight(1, 1, i) + &
                       InsetPData(inset_indp(1, 3, i), inset_indp(1, 2, i))*inset_weight(1, 2, i) + &
                       InsetPData(inset_indp(1, 1, i), inset_indp(1, 4, i))*inset_weight(1, 3, i) + &
                       InsetPData(inset_indp(1, 3, i), inset_indp(1, 4, i))*inset_weight(1, 4, i)) &
                      /rhoWat0g
         end if
         if (inset_indp(ub, 1, i) < 0) then
            WVNX2(i) = WVNX2(i); WVNY2(i) = WVNY2(i) !Keeping the background global wind
         else

            WVNX2(i) = InsetUData(inset_indp(1, 1, i), inset_indp(1, 2, i))*inset_weight(1, 1, i) + &
                       InsetUData(inset_indp(1, 3, i), inset_indp(1, 2, i))*inset_weight(1, 2, i) + &
                       InsetUData(inset_indp(1, 1, i), inset_indp(1, 4, i))*inset_weight(1, 3, i) + &
                       InsetUData(inset_indp(1, 3, i), inset_indp(1, 4, i))*inset_weight(1, 4, i)

            WVNY2(i) = InsetVData(inset_indp(1, 1, i), inset_indp(1, 2, i))*inset_weight(1, 1, i) + &
                       InsetVData(inset_indp(1, 3, i), inset_indp(1, 2, i))*inset_weight(1, 2, i) + &
                       InsetVData(inset_indp(1, 1, i), inset_indp(1, 4, i))*inset_weight(1, 3, i) + &
                       InsetVData(inset_indp(1, 3, i), inset_indp(1, 4, i))*inset_weight(1, 4, i)
         end if
      end do

   end subroutine INTERP_INSET_GRID
#endif
!-----------------------------------------------------------------------------------

#ifdef ADCNETCDF
   subroutine READ_INSET_LATLON(ncfile, NCIDf, lat_ins, &
                                InsetLatDim, InsetLat, lon_ins, InsetLonDim, InsetLon, nxp, nyp)

      use netcdf, only: NF90_OPEN, NF90_NOWRITE, NF90_CLOSE, &
                        NF90_INQ_DIMID, NF90_INQUIRE_DIMENSION, &
                        NF90_INQ_VARID, NF90_GET_VAR, &
                        NF90_INQUIRE_VARIABLE
      use netcdf_error, only: check_err
      use GL2LOC_MAPPING, only: BcastToLocal_2DRealArray, BcastToLocal_Int

      implicit none

      character(len=*), intent(in) :: ncfile ! NetCDF file name
      real, allocatable, dimension(:, :), intent(out) :: lat_ins, lon_ins

      character(LEN=80), intent(in)  :: InsetLat, InsetLon, InsetLatDim, InsetLonDim
      integer, intent(out) :: NXP, NYP
      integer, intent(out) :: NCIDf

      integer :: i, Temp_ID, Lat_ID, Lon_ID, ndims

      call setMessageSource("READ_INSET_LATLON")

      call Check_err(NF90_OPEN(ncfile, nf90_nowrite, NCIDf))

      call Check_err(NF90_INQ_DIMID(ncidf, InsetLatDim, Temp_ID))
      call Check_err(NF90_INQUIRE_DIMENSION(ncidf, Temp_ID, &
                                            len=NYP))
      call Check_err(NF90_INQ_DIMID(ncidf, InsetLonDim, Temp_ID))
      call Check_err(NF90_INQUIRE_DIMENSION(ncidf, Temp_ID, &
                                            len=NXP))

      allocate (lon_ins(NXP, NYP), lat_ins(NXP, NYP))

      call Check_err(NF90_INQ_VARID(ncidf, InsetLat, Lat_ID))
      call Check_err(nf90_inquire_variable(ncidf, Lat_ID, ndims=ndims))
      call Check_err(NF90_INQ_VARID(ncidf, InsetLon, Lon_ID))

      if (ndims == 2) then
         call Check_err(NF90_GET_VAR(ncidf, Lat_ID, lat_ins))
         call Check_err(NF90_GET_VAR(ncidf, Lon_ID, lon_ins))
      else
         call Check_err(NF90_GET_VAR(ncidf, Lat_ID, lat_ins(1, :)))
         call Check_err(NF90_GET_VAR(ncidf, Lon_ID, lon_ins(:, 1)))
         do i = 2, NXP
            lat_ins(i, :) = lat_ins(1, :)
         end do
         do i = 2, NYP
            lon_ins(:, i) = lon_ins(:, 1)
         end do
      end if

      call Check_err(NF90_CLOSE(NCIDf))

      call unsetMessageSource()
      return
   end subroutine READ_INSET_LATLON
!----------------------------------------------------------------------
#endif

!----------------------------------------------------------------------
   subroutine READNWS14LatLon(fileN, lat, lon, invN, var, ncid)
!----------------------------------------------------------------------
!     CPB 10/2023: Added to make NWS = 14 more more readable. Simply
!     calls the appropriate subroutine based on the filetype we are
!     using.
!----------------------------------------------------------------------
      implicit none
      character(200), intent(in) :: fileN
      real, allocatable, intent(out) :: lat(:, :), lon(:, :)
      character(200), intent(in), optional :: var, invN
      integer, intent(inout), optional :: ncid

      if (grb2flag) then
#ifdef GRIB2API
         call READNWS14LatLon_grib2(fileN, invN, var, lat, lon)
#else
         write (scratchMessage, '(7a)') 'GRIB2 support not compiled in ADCIRC.', &
            ' Cannot read NWS=14 lat/lon from GRIB2 file: ', trim(fileN), &
            ' Please recompile with GRIB2API option. invN: ', trim(invN), ' var: ', trim(var)
         call allMessage(ERROR, scratchMessage)
         call nws14Terminate()
#endif
      elseif (.not. grb2flag) then
#ifdef ADCNETCDF
         call READNWS14LatLon_netCDF(fileN, ncid, lat, lon)
#else
         ! Force compiler not to warn about unused variables
         ncid = ncid
         lat = lat
         lon = lon
         write (scratchMessage, '(5a)') 'netCDF support not compiled in ADCIRC.', &
            ' Cannot read NWS=14 lat/lon from netCDF file: ', trim(fileN), &
            ' Please recompile with ADCNETCDF option.'
         call allMessage(ERROR, scratchMessage)
         call nws14Terminate()
#endif
      else
         call nws14Terminate()
      end if
!----------------------------------------------------------------------
   end subroutine READNWS14LatLon
!----------------------------------------------------------------------

#ifdef GRIB2API
!----------------------------------------------------------------------
   subroutine READNWS14LatLon_grib2(fileN, invN, var, lat, lon)
!----------------------------------------------------------------------
!     CPB 10/2023: Reads in the lat and lon from a grib2 format
!     meteorological file. NOTE: reads in on Proc 0 and broadcasts
!----------------------------------------------------------------------
      use wgrib2api, only: grb2_inq
      use GL2LOC_MAPPING, only: BcastToLocal_2DRealArray, &
                                BcastToLocal_Int
      use SIZES, only: MYPROC
      implicit none

      integer :: iret, iter
      character(len=200), intent(in) :: var, fileN, invN
      real, allocatable, intent(out) :: lat(:, :), lon(:, :)
      integer :: NX, NY
      character(len=200) :: str_date
      call setMessageSource("READNWS14LatLon_grib2")
      ! Getting the date in str_date format
      str_date = ':start_FT='//trim(CurDT%strftime("%Y%m%d%H"))//'0000:'
      ! read in on core 0
      if (MYPROC == 0) then
         iret = -1; iter = 0
         do while (iret <= 0 .and. iter < 10)
            if (iter > 0) then
               call logMessage(WARNING, 'Trying to read again '//var)
               call sleep(5)
            end if
            iret = grb2_inq(fileN, invN, var, str_date, &
                            lat=lat, lon=lon)
            iter = iter + 1
         end do
         if (iret > 0) then
            call logMessage(ECHO, 'Successfully read LatLon '//var)
         else
            call ArnoldSchwarzenegger(iret, var, fileN)
         end if
      end if

      ! broadcast data
      if (MYPROC == 0) then
         NX = ubound(lon, 1)
         NY = ubound(lat, 2)
      end if
      call BcastToLocaL_Int(NX)
      call BcastToLocal_Int(NY)
      if (MYPROC /= 0) then
         allocate (lat(NX, NY), lon(NX, NY))
      end if
      ! latitude
      call BcastToLocal_2DRealArray(lat, NX, NY)
      ! longitude
      call BcastToLocal_2DRealArray(lon, NX, NY)
      call unsetMessageSource()
      return
!----------------------------------------------------------------------
   end subroutine READNWS14LatLon_grib2
!----------------------------------------------------------------------
#endif

!----------------------------------------------------------------------
#ifdef ADCNETCDF
   subroutine READNWS14LatLon_netCDF(fileN, ncid, lat, lon)
!----------------------------------------------------------------------
!     CPB 10/2023: Reads in the lat and lon from a netCDF format
!     meteorological forcing file. NOTE: reads in on Proc 0 and
!     broadcasts.
!----------------------------------------------------------------------

      use netcdf, only: NF90_OPEN, NF90_NOWRITE, NF90_CLOSE, &
                        NF90_INQ_DIMID, NF90_INQUIRE_DIMENSION, &
                        NF90_INQ_VARID, NF90_GET_VAR, nf90_inquire_variable
      use netcdf_error, only: check_err

      use GL2LOC_MAPPING, only: BcastToLocal_2DRealArray, &
                                BcastToLocal_Int
      use SIZES, only: MYPROC
      implicit none
      integer, intent(inout) :: ncid ! already opened ncid tag
      character(len=200), intent(in) :: fileN
      real, allocatable, intent(out) :: lat(:, :), lon(:, :)
      integer :: i, Temp_ID, Lat_ID, Lon_ID, NX, NY, ndims
      call setMessageSource("READNWS14LatLon_netCDF")

      if (MYPROC == 0) then

         if (.not. read_NWS14_NetCdf_using_core_0) then
            call Check_err(NF90_OPEN(FileN, nf90_nowrite, ncid))
         end if

         call Check_err(NF90_INQ_DIMID(ncid, Latdim, Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(ncid, Temp_ID, &
                                               len=NY))
         call Check_err(NF90_INQ_DIMID(ncid, Londim, Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(ncid, Temp_ID, &
                                               len=NX))
         allocate (lon(NX, NY), lat(NX, NY))
         call Check_err(NF90_INQ_VARID(ncid, Latvar, Lat_ID))
         call Check_err(nf90_inquire_variable(ncid, Lat_ID, &
                                              ndims=ndims))
         call Check_err(NF90_INQ_VARID(ncid, Lonvar, Lon_ID))
         if (ndims == 3) then
            call Check_err(NF90_GET_VAR(ncid, Lat_ID, lat, &
                                        start=[1, 1, NT], count=[NX, NY, 1]))
            call Check_err(NF90_GET_VAR(ncid, Lon_ID, lon, &
                                        start=[1, 1, NT], count=[NX, NY, 1]))
         elseif (ndims == 2) then
            call Check_err(NF90_GET_VAR(ncid, Lat_ID, lat))
            call Check_err(NF90_GET_VAR(ncid, Lon_ID, lon))
         else
            call Check_err(NF90_GET_VAR(ncid, Lat_ID, lat(1, :)))
            call Check_err(NF90_GET_VAR(ncid, Lon_ID, lon(:, 1)))
            do i = 2, NX
               lat(i, :) = lat(1, :)
            end do
            do i = 2, NY
               lon(:, i) = lon(:, 1)
            end do
         end if

         if (.not. read_NWS14_NetCdf_using_core_0) then
            call Check_err(NF90_CLOSE(ncid))
         end if

      end if
      call BcastToLocal_Int(NX)
      call BcastToLocal_Int(NY)
      if (MYPROC /= 0) then
         allocate (lon(NX, NY), lat(NX, NY))
      end if
      call BcastToLocal_2DRealArray(lon, NX, NY)
      call BcastToLocal_2DRealArray(lat, NX, NY)

      call unsetMessageSource()
!----------------------------------------------------------------------
   end subroutine READNWS14LatLon_netCDF
#endif
!----------------------------------------------------------------------

#ifdef ADCNETCDF
!----------------------------------------------------------------------
   subroutine READNWS14_NC_StaTime(FileN, NC_ID)
!----------------------------------------------------------------------
!       CPB 10/2023: Sets the time index from which we start reading in
!       netCDF meteorological forcing. NOTE: reads in on proc 0 and
!       broadcasts.
!----------------------------------------------------------------------
      use netcdf, only: NF90_OPEN, NF90_NOWRITE, NF90_CLOSE, &
                        NF90_INQ_DIMID, NF90_INQUIRE_DIMENSION, &
                        NF90_INQ_VARID, NF90_GET_VAR
      use netcdf_error, only: check_err
      use GL2LOC_MAPPING, only: BcastToLocal_Int
      use SIZES, only: MYPROC
#ifdef CMPI
      use SIZES, only: MYPROC
#endif

      use mod_datetime, only: operator(==), operator(+)
      implicit none
      integer, intent(inout) :: NC_ID
      character(len=200), intent(in) :: fileN
      character(len=19), allocatable :: TimeStr(:)
      character(len=200) :: str_date
      integer :: Temp_ID, TL
      ! CPB 4/20/2023
      real(8), allocatable :: NCTIMES(:)
      call setMessageSource('READNWS14_NC_StaTime')

      if (MYPROC == 0) then
         if (.not. read_NWS14_NetCdf_using_core_0) then
            call Check_err(NF90_OPEN(FileN, nf90_nowrite, NC_ID))
         end if

         call Check_err(NF90_INQ_DIMID(NC_ID, Tdim, Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(NC_ID, Temp_ID, &
                                               len=TL))
         if (Tformat(1:1) == '%') then
            ! in this case the netCDF files have a datetime variable
            ! so we can determine where to start based off of that
            call logMessage(ECHO, "fort.22 indicates that a " &
                            //"datetime variable is provided. Starting index " &
                            //"will be calculated.")
            allocate (TIMESTR(TL))
            call Check_err(NF90_INQ_VARID(NC_ID, Tvar, Temp_ID))
            call Check_err(NF90_GET_VAR(NC_ID, Temp_ID, TimeStr))
            ! Determine the first timestep index
            str_date = CurDT%strftime(trim(Tformat))
            do NT = 1, TL
               if (trim(str_date) == TimeStr(NT)) exit
            end do
            call logMessage(ECHO, 'Starting from '//TimeStr(NT))
            deallocate (TIMESTR)
         elseif (Tvar == 'refdate') then
            ! in this case the reference date for the netCDF file is
            ! provided in the fort.22 and the units of the time
            ! variable are in the form of "seconds since refdate".
            ! This allows us to find the start index
            call logMessage(ECHO, "fort.22 provides a reference " &
                            //"date. Starting time index will be calculated.")
            allocate (NCTIMES(TL))
            refdate = t_datetime(trim(Tformat))
            call CHECK_ERR(NF90_INQ_VARID(NC_ID, TDIM, TEMP_ID))
            call Check_err(NF90_GET_VAR(NC_ID, Temp_ID, NCTIMES))
            do NT = 1, TL
               stepdate = refdate + t_timedelta(seconds=int(NCTIMES(NT)))
               if (stepdate == curDT) then
                  exit
               end if
            end do
            deallocate (NCTIMES)
         else
            ! If neither of those two things are true we just assume
            ! we start at the beginning of the file
            call logMessage(ECHO, &
                            'Neither a datetime variable nor a ' &
                            //'reference date were provided. Assume met forcing ' &
                            //'starts at the beginning of the netCDF file.')
            NT = 1
         end if
         NTC = NT

         if (.not. read_NWS14_NetCdf_using_core_0) then
            call Check_err(NF90_CLOSE(NC_ID))
         end if

      end if

#ifdef CMPI
      call BcastToLocal_Int(NT)
      call BcastToLocal_Int(NTC)
#endif

      call unsetMessageSource()
!----------------------------------------------------------------------
   end subroutine READNWS14_NC_StaTime
!----------------------------------------------------------------------
#endif

#ifdef GRIB2API
!----------------------------------------------------------------------
   subroutine ArnoldSchwarzenegger(iret, var, fileN)
!----------------------------------------------------------------------
      implicit none
      integer, intent(in) :: iret
      character(len=200), intent(in) :: var, fileN

      if (iret == 0) then
         call allMessage(ERROR, 'Cound not find message when reading ' &
                         //trim(var)//' in '//trim(fileN)//'.')
      elseif (iret < 0) then
         call allMessage(ERROR, 'Fatal error when reading ' &
                         //trim(var)//' in '//trim(fileN)//'.')
      elseif (iret > 1) then
         call allMessage(ERROR, 'Found multiple messages when reading ' &
                         //trim(var)//' in '//trim(fileN)//'.')
      end if
      call nws14Terminate()
!----------------------------------------------------------------------
   end subroutine ArnoldSchwarzenegger
!-----------------------------------------------------------------------
#endif

!---------------------------------------------------------------------
   subroutine CLOSE_NWS14_FILES()

#ifdef ADCNETCDF
      use global, only: NWS, NCICE
      use netcdf, only: NF90_CLOSE
      use netcdf_error, only: check_err
      use SIZES, only: MYPROC
#endif

      implicit none

#ifdef ADCNETCDF
      if (NWS == 14) then
         if ((.not. grb2flag) .and. read_NWS14_NetCdf_using_core_0) then
            if (MYPROC == 0) then
               call Check_err(NF90_CLOSE(PfileNCID))
               call Check_err(NF90_CLOSE(WfileNCID))
               call Check_err(NF90_CLOSE(Wfile1NCID))

               if (NCICE == 14) then
                  call Check_err(NF90_CLOSE(CfileNCID))
               end if
            end if
         end if
      end if
#endif
      return
   end subroutine CLOSE_NWS14_FILES
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!...  jgf50.38.05: Subroutine to terminate the run cleanly.
!...
!----------------------------------------------------------------------
   subroutine nws14Terminate()
#ifdef CMPI
      use MESSENGER, only: MSG_FINI
#endif
      implicit none

      call setMessageSource("nws14Terminate")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      call allMessage(ERROR, "ADCIRC terminating.")
#ifdef CMPI
      call msg_fini()
#endif
      call exit(1)

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
!----------------------------------------------------------------------
   end subroutine nws14Terminate
!----------------------------------------------------------------------

end module mod_nws14
