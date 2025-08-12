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
!-----------------------------------------------------------------------
!  MODULE mod_nws12
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
module mod_nws12

   use SIZES, only: MYPROC
   use GLOBAL, only: NSCREEN, ScreenUnit, DEBUG, ECHO, INFO, &
                     WARNING, ERROR, screenMessage, logMessage, allMessage, &
                     setMessageSource, unsetMessageSource, scratchMessage, &
                     openFileForRead, Flag_ElevError
   use mod_datetime, only: t_datetime
#ifdef CMPI
   use MESSENGER, only: MSG_FINI
#endif
   implicit none

   integer :: nPowellSearchDomains = -1 !< Number of search domains for Powell wind drag
   integer :: numSets !< Number of sets of wind and pressure fields
   integer :: numBlankSnaps !< Number of blank snaps to prepend to Oceanweather files
   integer :: numSkipSnaps !< Number of snaps to skip in Oceanweather files
   integer :: cntSnaps !< Counter for number of snaps that have been processed so far
   real(8) :: windMultiplier !< Wind multiplier applied to Oceanweather wind velicities
   logical :: moving_grid = .false. !< indicates the owi domains use a moving grid (for aswip.F)

   type OCEANWEATHER !< Container for Oceanweather domains. Represents a wind and pressure field
      character(1024) :: pressure_file !< filename used for the pressure field
      character(1024) :: wind_file !< filename used for the wind velocity field
      type(t_datetime) :: startDate !< start date in the domain
      type(t_datetime) :: endDate !< end date in the domain
      type(t_datetime) :: date !< current date in domain for data read
      integer :: iLat !< number of points in latitude direction
      integer :: iLon !< number of points in longitude direction
      integer :: iounit_pressure !< fortran unit number for pressure file
      integer :: iounit_wind !< fortran unit number for wind file
      integer :: iSnap !< current snap number
      integer :: iUpdate !< flag to update interpolation weights and indicies if wind field moves
      real(8) :: dLat !< spacing in latitude direction
      real(8) :: dLon !< spacing in longitude direction
      real(8) :: swLon !< southwest longitude
      real(8) :: swLat !< southwest latitude
      logical :: hasData !< set to false when fields are no longer have data
      logical :: atEnd !< set to false when the fields are at the end of file
      real(8), allocatable :: p(:, :) !< pressure field
      real(8), allocatable :: u(:, :) !< u-velocity field
      real(8), allocatable :: v(:, :) !< v-velocity field
      real(8), allocatable :: weight(:, :) !< bilinear interpolation weights
      integer, allocatable :: interpPoints(:, :) !< indicies for interpolation
   end type OCEANWEATHER

   type(OCEANWEATHER), allocatable :: owi(:) !< Array containing each of the wind and pressure field used

   private

   public nws12get, nws12init, findStormCenter, windMultiplier
   public nPowellSearchDomains, moving_grid

   !---------------------end of data declarations--------------------------------

contains

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
   subroutine NWS12INIT(WVNX, WVNY, PRN, NP, PRdeflt)
      use SIZES, only: GBLINPUTDIR
      implicit none
      real(8), intent(out) :: WVNX(:), WVNY(:), PRN(:)
      integer, intent(in) :: NP
      integer :: I
      real(8), intent(in) :: PRdeflt
      character(len=1024) :: errorVar
      character(len=1024) :: fname
      integer :: iounit
      integer :: errorIO

      call setMessageSource("nws12init")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      !...Begin reading fort.22
      fname = "fort.22"
      iounit = 22
      errorVar = "open"
      errorIO = 0
      call openFileForRead(iounit, trim(GBLINPUTDIR)//"/fort.22", errorIO)
      call check_owi_err(errorIO, errorVar, fname, iounit)

      errorVar = "reading number of wind sets" ! used in error msgs
      read (iounit, *, err=99999, end=99998, iostat=errorIO) numSets
      call check_owi_err(errorIO, errorVar, fname, iounit)

      if (numSets == 0) then
         call check_owi_err(1, "number of wind sets must not be zero", fname, iounit)
      end if

      ! Read the number of blank snaps to be inserted before OWI winds start
      errorVar = "reading number of blank snaps" ! used in error messages
      read (iounit, *, err=99999, end=99998, iostat=errorIO) numBlankSnaps
      call check_owi_err(errorIO, "read number of blank snaps", fname, iounit)
      ! If numBlankSnaps < 0, ABS(numBlankSnaps) snaps
      ! in OWI wind files (UNIT 221,222,223, 224 and 217 & 218) will be skipped.
      if (numBlankSnaps < 0) then
         numSkipSnaps = abs(numBlankSnaps)
         numBlankSnaps = 0
      else
         numSkipSnaps = 0
      end if

      ! Read a wind velocity multiplier
      errorVar = "reading wind multiplier" ! used in error messages
      read (iounit, *, err=99999, end=99998, iostat=errorIO) windMultiplier
      call check_owi_err(errorIO, errorVar, fname, iounit)

      errorIO = 0

      ! Allocate the owi arrays
      allocate (owi(1:abs(numSets)))
      call owi_setFilenames(owi, iounit, numSets)
      close (iounit)

      ! set the files to contain data
      do i = 1, abs(numSets)
         owi(i)%atEnd = .false.
         owi(i)%hasData = .false.
      end do

      ! Check the number of powell search domains
      if (nPowellSearchDomains > abs(numSets)) then
         nPowellSearchDomains = abs(numSets)
      elseif (nPowellSearchDomains <= 0) then
         nPowellSearchDomains = abs(numSets)
      end if

      !...Initialize the wind objects based on file headers
      call owi_initializeFileHeaders(owi)

      !...Allocate interpolation arrays
      call owi_allocateInterpolationArrays(owi, np)

      !...Set snap counter to zero
      cntSnaps = 0

      !...Skip specified snaps
      call owi_skipSnaps(numSkipSnaps, np, WVNX, WVNY, PRN, PRdeflt)

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return

99998 call allMessage(ERROR, "Unexpectedly reached end-of-file.") ! END jumps here

99999 call check_owi_err(errorIO, errorVar, fname, iounit) !  ERR jumps here

   end subroutine NWS12INIT
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
   subroutine NWS12GET(WVNX, WVNY, PRN, NP, PRdeflt, HasData)
      use GLOBAL, only: NWS
      use ADC_CONSTANTS, only: RHOWAT0, G
      use mod_datetime, only: operator(<), operator(==), operator(>)
      implicit none
      integer, intent(in) :: NP
      real(8), intent(out) :: WVNX(:), WVNY(:), PRN(:)
      real(8), intent(in) :: PRdeflt
      real(8) :: RHOWATG
      integer :: I, S
      character(1024) :: errorVar
      real(8) :: uu, vv, pp
      logical :: uvpHasBeenSet
      logical, intent(out), optional :: HasData(:)

      call setMessageSource("nws12get")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif
      errorVar = ""
      RHOWATG = RHOWAT0*G

      ! Increment counter (cntSnaps initialized to zero in nws12init)
      cntSnaps = cntSnaps + 1

      ! Put a blank snap for the first 'numBlankSnaps' snaps and then return
      if (cntSnaps <= numBlankSnaps) then
         if (abs(NWS) /= 14) then ! WJP: Don't overwrite if NWS = 14
            WVNX(:) = 0d0
            WVNY(:) = 0d0
            PRN(:) = Prdeflt*100d0/RHOWATG
         end if
         if (present(HasData)) then
            HasData(:) = .false.
         end if
         write (scratchMessage, 16) cntSnaps
16       format('INSERTING A BLANK WIND SNAP, COUNT=', i4)
         call allMessage(INFO, trim(scratchMessage))
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG, "Return.")
#endif
         call unsetMessageSource()
         return
      end if

      do S = 1, numSets

         if (s > 1) then
            if ((owi(1)%date < owi(s)%startDate) .or. &
                (owi(1)%date == owi(s)%endDate .and. owi(s)%date%minute() /= 0) .or. &
                (owi(1)%date > owi(s)%endDate)) then
               owi(s)%hasData = .false.
               cycle
            end if
         end if

         if (owi(s)%atEnd) then
            !  Set the blank snap and report to screen if past end
            !  of wind data
            call owi_setEndOfFileBlankSnap(owi(s)%pressure_file, &
                                           prdeflt, owi(s))
         else
            !  Read the next snap from the wind and pressure files
            call owi_readNextSnap(s, np, prdeflt)

            ! Increment counter
            owi(s)%isnap = owi(s)%isnap + 1
         end if

      end do

      do i = 1, NP
         uvpHasBeenSet = .false.

         call owi_interpolateToMesh(i, uu, vv, pp, uvpHasBeenSet)

         if (.not. uvpHasBeenSet) then
            if (abs(NWS) /= 14) then ! WJP: Don't overwrite if NWS = 14
               WVNX(I) = 0.d0
               WVNY(I) = 0.d0
               PRN(I) = PRdeflt*100.d0/RHOWATG
            end if
         else
            !             !Convert millibars to m of water
            PRN(i) = 100.d0*PP/RHOWATG

            ! Apply wind velocity multiplier
            ! and save to arrays
            WVNX(i) = uu*windMultiplier
            WVNY(i) = vv*windMultiplier

         end if

         if (present(HasData)) then
            HasData(i) = uvpHasBeenSet
         end if

      end do

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return

   end subroutine NWS12GET
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Called to find the lowest pressure in the met field below
   !>        the threshold pressure. Loops over all domains and finds
   !>        the lowest pressure across all domains
   !> @param[inout] eyeLat prior 3 locations (y) where the eye was found
   !> @param[inout] eyeLon prior 3 locations (x) where the eye was found
   !> @param[out]   foundEye set to true if the eye has been found
   !-----------------------------------------------------------------------
   subroutine findStormCenter(eyeLat, eyeLon, foundEye)
      !-----------------------------------------------------------------------
      use ADC_CONSTANTS, only: PRBCKGRND
      implicit none
      intrinsic :: minloc
      real(8), intent(inout) :: eyeLat(3)
      real(8), intent(inout) :: eyeLon(3)
      logical, intent(out) :: foundEye
      integer :: EyeLatI
      integer :: EyeLonI
      integer :: EyeDomainI
      real(8) :: EyeLatTemp
      real(8) :: EyeLonTemp
      real(8) :: EyePressure
      real(8) :: minp
      integer :: s
      integer :: minidx(2)
      real(8), parameter :: eps = epsilon(1d0)

      call setMessageSource("findStormCenter")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      EyeLatI = 0
      EyeLonI = 0
      EyeDomainI = 0
      EyeLatTemp = 0d0
      EyeLonTemp = 0d0
      EyePressure = PRBCKGRND

      do S = 1, nPowellSearchDomains
         if (owi(s)%hasData) then
            minidx = minloc(owi(s)%p)
            minp = owi(s)%p(minidx(1), minidx(2))
            if (minp < 1000d0 .and. minp < eyePressure) then
               EyeLonI = minidx(1)
               EyeLatI = minidx(2)
               EyeDomainI = s
               EyePressure = owi(s)%p(minidx(1), minidx(2))
            end if
         end if
      end do

      if (EyeDomainI /= 0) then
         if ((EyeLatI == 0) .or. (EyeLonI == 0) .or. (EyeDomainI == 0) .or. &
             (EyeLatI == 1) .or. (EyeLatI == owi(EyeDomainI)%iLat) .or. &
             (EyeLonI == 1) .or. (EyeLonI == owi(EyeDomainI)%iLon)) then
            FoundEye = .false.
         else
            FoundEye = .true.
            EyeLatTemp = owi(EyeDomainI)%swlat + dble(EyeLatI - 1)*owi(EyeDomainI)%dLat
            EyeLonTemp = owi(EyeDomainI)%swlon + dble(EyeLonI - 1)*owi(EyeDomainI)%dLon
            if ((abs(EyeLatTemp - EyeLat(3)) >= eps) .or. (abs(EyeLonTemp - EyeLon(3)) >= eps)) then
               EyeLat(1) = EyeLat(2)
               EyeLon(1) = EyeLon(2)
               EyeLat(2) = EyeLat(3)
               EyeLon(2) = EyeLon(3)
               EyeLat(3) = EyeLatTemp
               EyeLon(3) = EyeLonTemp
            end if
         end if
      else
         FoundEye = .false.
      end if

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()

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
   subroutine owi_readNextSnap(idx, np, prbk)
      use mod_datetime, only: operator(/=)
      implicit none
      integer, intent(in) :: idx
      integer, intent(in) :: np
      real(8), intent(in) :: prbk
      character(1024) :: filename
      character(1024) :: errorVar
      integer :: i, j
      integer :: iounit, errorIO
      integer :: iLatw, iLongw, iLatp, iLongp
      character(20) :: iCYMDHMw, iCYMDHMp
      real(8) :: dxw, dyw, swlatw, swlongw
      real(8) :: dxp, dyp, swlatp, swlongp
      real(8), parameter :: eps = epsilon(1d0)
      type(t_datetime) :: current_date_wind, current_date_press

      call setMessageSource("owi_readNextSnap")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      owi(idx)%hasData = .false.

      ! Read grid specifications/date in pressure file
      filename = owi(idx)%pressure_file
      iounit = owi(idx)%iounit_pressure
      write (errorVar, '(A,I0)') "read grid specifications/date in Oceanweather pressure file, domain ", idx
      read (owi(idx)%iounit_pressure, 11, end=10000, err=9999, iostat=errorIO) &
         iLatp, iLongp, dxp, dyp, swlatp, swlongp, iCYMDHMp
      call check_owi_err(errorIO, errorVar, filename, iounit)
      current_date_press = t_datetime(iCYMDHMp, "%Y%m%d%H%M")

      if (idx == 1) then
         if (numSets > 1) then
            write (scratchMessage, '("Processing Oceanweather wind data for domain ",I0," through "'// &
                   ',I0," for time stamp ",A)') idx, numSets, trim(current_date_press%to_iso_string())
         else
            write (scratchMessage, '("Processing Oceanweather wind data for time stamp ",A)') &
               trim(current_date_press%to_iso_string())
         end if
         call allmessage(INFO, scratchMessage)
      end if

      ! Read grid specifications/date in wind file
      filename = owi(idx)%wind_file
      iounit = owi(idx)%iounit_wind
      write (errorVar, '(A,I0)') "read grid specifications/date in Oceanweather wind file, domain ", idx
      read (owi(idx)%iounit_wind, 11, end=10000, err=9999, iostat=errorIO) &
         iLatw, iLongw, dxw, dyw, swlatw, swlongw, iCYMDHMw
      call check_owi_err(errorIO, errorVar, filename, iounit)
      current_date_wind = t_datetime(iCYMDHMw, "%Y%m%d%H%M")

      ! Check consistency
      if (iLatp /= iLatw .or. &
          iLongp /= iLongw .or. &
          abs(dxp - dxw) > eps .or. &
          abs(dyp - dyw) > eps .or. &
          abs(swlatp - swlatw) > eps .or. &
          abs(swlongp - swlongw) > eps .or. &
          current_date_press /= current_date_wind) then
         call allMessage(ERROR, &
                         "Grid specifications/date in OWI win and pre files must match.")
         errorVar = ""
         call check_owi_err(errorIO, errorVar, filename, iounit)
      end if

      if (owi(idx)%isnap > 1) then
         if (iLatp /= owi(idx)%iLat .or. iLongp /= owi(idx)%iLon .or. &
             abs(dxp - owi(idx)%dlon) > eps .or. abs(dyp - owi(idx)%dlat) > eps .or. &
             abs(swlatp - owi(idx)%swlat) > eps .or. abs(swlongp - owi(idx)%swlon) > eps) then
            write (scratchMessage, '(A,I0,A)') "Oceanweather domain ", idx, " has changed."
            call logMessage(INFO, trim(scratchMessage))
            owi(idx)%iupdate = 1
         else
            owi(idx)%iupdate = 0
         end if
      end if

      owi(idx)%date = current_date_press

      ! Update coordinate mapping coefficients if necessary
      if (owi(idx)%iupdate == 1 .or. moving_grid) then
         owi(idx)%ilon = ilongp
         owi(idx)%ilat = ilatp
         owi(idx)%dlon = dxp
         owi(idx)%dlat = dyp
         owi(idx)%swlat = swlatp
         owi(idx)%swlon = swlongp
         write (errorVar, '(A,I0)') "Updating grid coordinate mapping coefficients, domain ", idx
         call logMessage(INFO, errorVar)
         call nws12interp(idx, np)
      end if

      ! Read basin scale atmospheric pressure snapshot
      filename = owi(idx)%pressure_file
      iounit = owi(idx)%iounit_pressure
      write (errorVar, '(A,I0)') "reading atmospheric pressure snapshot, domain ", idx
      read (owi(idx)%iounit_pressure, 22, end=10000, err=9999, iostat=errorIO) &
         ((owi(idx)%p(i, j), i=1, owi(idx)%iLon), j=1, owi(idx)%iLat)
      call check_owi_err(errorIO, errorVar, filename, iounit)

      ! Read basin scale snapshot of u/v components of the wind
      filename = owi(idx)%wind_file
      iounit = owi(idx)%iounit_wind
      write (errorVar, '(A,I0)') "reading wind u-velocity snapshot, domain ", idx
      read (owi(idx)%iounit_wind, 22, end=10000, err=9999, iostat=errorIO) &
         ((owi(idx)%u(i, j), i=1, owi(idx)%iLon), j=1, owi(idx)%iLat)
      call check_owi_err(errorIO, errorVar, filename, iounit)
      filename = owi(idx)%wind_file
      iounit = owi(idx)%iounit_wind
      write (errorVar, '(A,I0)') "reading wind v-velocity snapshot, domain ", idx
      read (owi(idx)%iounit_wind, 22, end=10000, err=9999, iostat=errorIO) &
         ((owi(idx)%v(i, j), i=1, owi(idx)%iLon), j=1, owi(idx)%iLat)
      call check_owi_err(errorIO, errorVar, filename, iounit)

      ! Set this flag so interpolation routine knows
      ! to use this domain or not
      owi(idx)%hasData = .true.

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()

      return

      ! Oceanweather read formats
11    format(t6, i4, t16, i4, t23, f6.0, t32, f6.0, t44, f8.0, t58, f8.0, t69, a12)
22    format(8f10.0)

      ! Error IO
9999  call check_owi_err(errorIO, errorVar, filename, iounit) ! ERR during read jumps to here

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return

      ! END during read of grid, u, v, p data jumps to here
10000 continue
      call owi_setEndOfFileBlankSnap(filename, prbk, owi(idx))

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return

   end subroutine owi_readNextSnap
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Sets the blank snaps at the end of the file when data has run
   !>        out
   !> @param[in] filename name of file that is being read
   !> @param[in] prbk background pressure
   !> @param[inout] owi oceanweather object array
   !-----------------------------------------------------------------------
   subroutine owi_setEndOfFileBlankSnap(filename, prbk, obj)
      !-----------------------------------------------------------------------
      implicit none
      character(*), intent(in) :: filename
      real(8), intent(in) :: prbk
      type(OCEANWEATHER), intent(inout) :: obj

      call setMessageSource("owi_setEndOfFileBlankSnap")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      write (scratchMessage, 61) trim(filename)
61    format("End-of-file reached while reading '", A, &
             "'. Wind speeds set to zero and the pressure field", &
             " to its background value.")
      call allMessage(WARNING, trim(scratchMessage))

      ! Set no wind forcing for this domain
      obj%u(:, :) = 0d0
      obj%v(:, :) = 0d0
      obj%p(:, :) = prbk
      obj%hasData = .false.
      obj%atEnd = .true.

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()

      return
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
   subroutine owi_interpolateToMesh(node, u, v, p, isSet)
      implicit none
      integer, intent(in) :: node
      real(8), intent(out) :: u, v, p
      logical, intent(out) :: isSet
      integer :: i, j, xi, yi
      real(8) :: ui(4)
      real(8) :: vi(4)
      real(8) :: pi(4)
      real(8), parameter :: flag = -999.0d0
      real(8), parameter :: eps = epsilon(1d0)

      call setMessageSource("owi_interpolateToMesh")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      isSet = .false.

      do j = numSets, 1, -1
         if (.not. owi(j)%hasData) cycle
         if (owi(j)%interppoints(node, 1) > 0) then
            xi = owi(j)%interpPoints(node, 1)
            yi = owi(j)%interpPoints(node, 2)

            ui = [owi(j)%u(xi, yi), owi(j)%u(xi + 1, yi), &
                  owi(j)%u(xi + 1, yi + 1), owi(j)%u(xi, yi + 1)]
            vi = [owi(j)%v(xi, yi), owi(j)%v(xi + 1, yi), &
                  owi(j)%v(xi + 1, yi + 1), owi(j)%v(xi, yi + 1)]
            pi = [owi(j)%p(xi, yi), owi(j)%p(xi + 1, yi), &
                  owi(j)%p(xi + 1, yi + 1), owi(j)%p(xi, yi + 1)]

            if (any(abs(ui - flag) < eps) .or. &
                any(abs(vi - flag) < eps) .or. &
                any(abs(pi - flag) < eps)) cycle

            u = 0.0d0
            v = 0.0d0
            p = 0.0d0

            do i = 1, 4
               u = u + owi(j)%weight(node, i)*ui(i)
               v = v + owi(j)%weight(node, i)*vi(i)
               p = p + owi(j)%weight(node, i)*pi(i)
            end do

            isSet = .true.
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG, "Return.")
#endif
            call unsetMessageSource()
            return
         end if
      end do

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
   end subroutine owi_interpolateToMesh
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief  Allocates the u,v,p fields in an oceanweather object and
   !>         deallocates if necessary
   !> @param[inout] this array of oceanweather objects
   !-----------------------------------------------------------------------
   subroutine owi_allocate(this)
      implicit none
      type(OCEANWEATHER), intent(inout) :: this
      logical :: do_allocate
      call setMessageSource("owi_allocate")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif
      do_allocate = .false.

      if (.not. allocated(this%p)) then
         do_allocate = .true.
      elseif ((size(this%p, 1) /= this%iLon) .or. (size(this%p, 2) /= this%iLat)) then
         do_allocate = .true.
      end if

      if (do_allocate) then
         call owi_deallocate(this)
         allocate (this%p(this%ilon, this%ilat))
         allocate (this%u(this%ilon, this%ilat))
         allocate (this%v(this%ilon, this%ilat))
      end if

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
   end subroutine owi_allocate
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Deallocates the u,v,p fields in an oceanweather object
   !> @param[inout] this array of oceanweather objects
   !-----------------------------------------------------------------------
   subroutine owi_deallocate(this)
      implicit none
      type(OCEANWEATHER), intent(inout) :: this
      call setMessageSource("owi_deallocate")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif
      if (allocated(this%p)) deallocate (this%p)
      if (allocated(this%u)) deallocate (this%u)
      if (allocated(this%v)) deallocate (this%v)
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
   end subroutine owi_deallocate
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Generates coordinates for a specified oceanweather domain
   !> @param[in] this array of oceanweather objects
   !> @param[out] lat latitude coordinates of oceanweather grid
   !> @param[out] lon longiude coordinates of oceanweather grid
   !-----------------------------------------------------------------------
   subroutine owi_generateCoordinates(this, lat, lon)
      implicit none
      type(OCEANWEATHER), intent(in) :: this
      real(8), allocatable, intent(out) :: lat(:)
      real(8), allocatable, intent(out) :: lon(:)
      integer :: I

      call setMessageSource("owi_generateCoordinates")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif
      allocate (lat(1:this%ilat))
      allocate (lon(1:this%ilon))

      do I = 1, this%ilat
         lat(I) = this%swlat + dble(I - 1)*this%dLat
      end do
      do I = 1, this%ilon
         lon(i) = this%swlon + dble(I - 1)*this%dLon
      end do

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()

      return
   end subroutine owi_generateCoordinates
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Sets the filenames that will be used for the oceanweather wind
   !>        fields
   !> @param[inout] this array of oceanweather objects
   !> @param[in] iounit fortran unit number for fort.22
   !> @param[inout] nWindFields number of wind fields specified in fort.22
   !-----------------------------------------------------------------------
   subroutine owi_setFilenames(this, iounit, nWindFields)
      implicit none
      integer, intent(inout) :: nWindFields
      integer, intent(in) :: iounit
      type(oceanweather), intent(inout) :: this(nWindFields)

      call setMessageSource("owi_setFilenames")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      if (nWindFields < 0) then
         nWindFields = abs(nWindFields)
         call owi_setDynamicFilenames(this, iounit, nWindFields)
      else
         if (nWindFields > 3) then
            call check_owi_err(1, "owi file list must be used for numSets>3. "// &
                               "Use negative numSets to read a file list.", "fort.22", 22)
         end if
         call owi_setDefaultFilenames(this)
      end if
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
   end subroutine owi_setFilenames
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Sets the dynamic filenames read from the fort.22 file
   !> @param[inout] this array of oceanweather objects
   !> @param[in] iounit fortran unit number for fort.22
   !> @param[in] nWindFields number of wind fields specified in fort.22
   !-----------------------------------------------------------------------
   subroutine owi_setDynamicFilenames(this, iounit, nWindFields)
      use sizes, only: gblinputdir
      implicit none
      integer, intent(in) :: iounit
      integer, intent(in) :: nWindFields
      type(oceanweather), intent(inout) :: this(nwindFields)
      character(1024) :: f1, f2
      integer :: i

      call setMessageSource("owi_setDynamicFilenames")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif
      do i = 1, nWindFields
         read (iounit, *) f1, f2
         this(i)%pressure_file = trim(GBLINPUTDIR)//"/"//trim(f1)
         this(i)%wind_file = trim(GBLINPUTDIR)//"/"//trim(f2)
         this(i)%iounit_pressure = 500 + i
         this(i)%iounit_wind = 600 + i
      end do

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
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
   subroutine owi_skipSnaps(nSnaps, np, wx, wy, pr, prbk)
      implicit none
      integer, intent(in) :: nSnaps
      integer, intent(in) :: np
      real(8), intent(inout) :: wx(:), wy(:), pr(:)
      real(8), intent(in) :: prbk
      integer :: i
      call setMessageSource("owi_skipSnaps")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif
      do i = 1, nSnaps
         write (scratchMessage, 41) i
41       format("Skipping snap '", I6, "' in OWI wind data.")
         call logMessage(DEBUG, trim(scratchMessage))
         call NWS12GET(wx, wy, pr, np, prbk)
      end do
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
   end subroutine owi_skipSnaps
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Sets the default names for the oceanweather files if numSets is
   !>        positive. fort.221, fort.222, fort.223, fort.224, fort.217,
   !>        fort.218
   !> @param[inout] this array of oceanweather objects
   !-----------------------------------------------------------------------
   subroutine owi_setDefaultFilenames(this)
      use sizes, only: gblinputdir
      implicit none
      type(oceanweather), intent(inout) :: this(numSets)
      call setMessageSource("owi_setDefaultFilenames")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif
      this(1)%pressure_file = trim(GBLINPUTDIR)//"/fort.221"
      this(1)%wind_file = trim(GBLINPUTDIR)//"/fort.222"
      this(1)%iounit_pressure = 221
      this(1)%iounit_wind = 222
      if (numSets > 1) then
         this(2)%pressure_file = trim(GBLINPUTDIR)//"/fort.223"
         this(2)%wind_file = trim(GBLINPUTDIR)//"/fort.224"
         this(2)%iounit_pressure = 223
         this(2)%iounit_wind = 224
      end if
      if (numSets > 2) then
         this(3)%pressure_file = trim(GBLINPUTDIR)//"/fort.217"
         this(3)%wind_file = trim(GBLINPUTDIR)//"/fort.218"
         this(3)%iounit_pressure = 217
         this(3)%iounit_wind = 218
      end if
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
   end subroutine owi_setDefaultFilenames
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Initializes the data in the oceanweather data objects based on
   !>        the file headers
   !> @param[inout] this array of oceanweather objects
   !-----------------------------------------------------------------------
   subroutine owi_initializeFileHeaders(this)
      use mod_datetime, only: operator(/=)
      implicit none
      type(oceanweather), intent(inout) :: this(numSets)
      integer :: i
      type(t_datetime) :: date1pressure, date2pressure
      type(t_datetime) :: date1wind, date2wind
      character(1024) :: errorVar

      call setMessageSource("owi_initializeFileHeaders")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif
      do I = 1, numSets
         call readHeader(this(i)%pressure_file, this(i)%iounit_pressure, &
                         date1pressure, date2pressure)
         call readHeader(this(i)%wind_file, this(i)%iounit_wind, &
                         date1wind, date2wind)

         if (date1pressure /= date1wind .or. date2pressure /= date2wind) then
            call allMessage(ERROR, "Wind and pressure dates do not match.")
            errorVar = ""
            call check_owi_err(1, errorVar, this(i)%wind_file, this(i)%iounit_wind)
         end if

         this(i)%startDate = date1wind
         this(i)%endDate = date2wind

         if (this(i)%startDate /= this(1)%startdate) then
            call allMessage(ERROR, "Wind domains have different time frames")
            errorVar = ""
            call check_owi_err(1, errorVar, this(i)%wind_file, this(i)%iounit_wind)
         end if

         this(i)%isnap = 0
         this(i)%iupdate = 1
      end do

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
   end subroutine owi_initializefileHeaders
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Allocates the arrays used to interpolate oceanweather grids to
   !>        the mesh
   !> @param[inout] this array of oceanweather objects
   !> @param[in] numNodes number of nodes in the adcirc mesh
   !-----------------------------------------------------------------------
   subroutine owi_allocateInterpolationArrays(this, numNodes)
      implicit none
      type(oceanweather), intent(inout) :: this(numSets)
      integer, intent(in) :: numNodes
      integer :: i

      call setMessageSource("owi_initializeFileHeaders")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif
      do i = 1, numSets
         allocate (this(i)%interpPoints(numNodes, 2))
         allocate (this(i)%weight(numNodes, 4))
      end do

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
   end subroutine owi_allocateInterpolationArrays
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief This generates and saves interpolation coefficients for mapping
   !>        from a OWI grid to an ADCIRC mesh.
   !> @param[in] idx index in the owi object array
   !> @param[in] np number of points in the adcirc mesh
   !-----------------------------------------------------------------------
   subroutine NWS12INTERP(idx, NP)
      use ADC_CONSTANTS, only: RAD2DEG
      use MESH, only: SLAM, SFEA

      implicit none

      integer, intent(in) :: idx
      integer, intent(in) :: NP
      integer :: I, J, K, XI, YI
      real(8) :: adcLat, adcLong
      real(8) :: w, w1, w2, w3, w4
      real(8), allocatable :: lat(:), lon(:)

      call setMessageSource("nws12interp")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      write (16, '(A)') ''
      write (16, '(A,I0)') 'WIND MAPPING UPDATED FOR REGION ', IDX
      write (16, '(A)') ''

      call owi_allocate(owi(idx))
      call owi_generateCoordinates(owi(idx), lat, lon)

      xi = -1
      yi = -1

      ! Generate interpolation coefficients (south west point and weights)
      do i = 1, NP
         adcLat = RAD2DEG*SFEA(i)
         adcLong = RAD2DEG*SLAM(i)

         if (adcLong >= lon(1) .and. adcLong < lon(owi(idx)%iLon) .and. &
             adcLat >= lat(1) .and. adcLat < lat(owi(idx)%iLat)) then
            do j = 1, owi(idx)%iLon - 1
               if (adcLong >= lon(j) .and. &
                   adcLong < lon(j + 1)) then
                  xi = j
                  exit
               end if
            end do

            do k = 1, owi(idx)%iLat - 1
               if (adcLat >= lat(k) .and. &
                   adcLat < lat(k + 1)) then
                  yi = k
                  exit
               end if
            end do

            owi(idx)%interpPoints(i, 1) = xi
            owi(idx)%interpPoints(i, 2) = yi

            w = (lon(xi + 1) - lon(xi))*(lat(yi + 1) - lat(yi))
            w1 = (lon(xi + 1) - adcLong)*(lat(yi + 1) - adcLat)
            w2 = (adcLong - lon(xi))*(lat(yi + 1) - adcLat)
            w3 = (adcLong - lon(xi))*(adcLat - lat(yi))
            w4 = (lon(xi + 1) - adcLong)*(adcLat - lat(yi))

            owi(idx)%weight(i, 1) = w1/w
            owi(idx)%weight(i, 2) = w2/w
            owi(idx)%weight(i, 3) = w3/w
            owi(idx)%weight(i, 4) = w4/w

         else
            owi(idx)%interpPoints(i, 1) = 0
            owi(idx)%interpPoints(i, 2) = 0
         end if
      end do

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
   end subroutine NWS12INTERP
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> @brief Read the header from the OWI file and returns the values of
   !>        start and end dates
   !> @param[in] filename name of the file to open and read
   !> @param[in] lun unit number to open the file as
   !> @param[out] date1 start date in the header
   !> @param[out] date2 end date in the header
   !-----------------------------------------------------------------------
   subroutine readHeader(filename, lun, date1, date2)
      implicit none
      integer, intent(in) :: lun
      character(*), intent(in) :: filename
      character(10) :: date1str, date2str
      type(t_datetime), intent(out) :: date1
      type(t_datetime), intent(out) :: date2
      character(80) :: owiheader
      character(1024) :: errorVar
      integer :: errorIO

      call setMessageSource("readHeader")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif
      call openFileForRead(lun, filename, errorIO)
      errorVar = "opening Oceanweather file"
      call check_owi_err(errorIO, errorVar, filename, lun)

      ! Read begining/ending dates of pre file
      owiheader(:) = ' ' !set owiheader to blanks before read
      errorVar = "reading Oceanweather header line"
      read (lun, fmt='(a80)', end=99998, err=99999, iostat=errorIO) owiheader
      call check_owi_err(errorIO, errorVar, filename, lun)

      errorVar = "reading Oceanweather start date"
      read (owiheader(56:65), '(A10)', end=99998, err=99999, iostat=errorIO) date1str
      call check_owi_err(errorIO, errorVar, filename, lun)
      date1 = t_datetime(date1str, "%Y%m%d%H")
      write (scratchMessage, 31) trim(errorVar), trim(filename), trim(date1%to_iso_string())
31    format("'", A, "' in  '", A, "' is '", A, "'.")
      call allMessage(INFO, scratchMessage)

      errorVar = "reading Oceanweather end date"
      read (owiheader(71:80), '(A10)', end=99998, err=99999, iostat=errorIO) date2str
      call check_owi_err(errorIO, errorVar, filename, lun)
      date2 = t_datetime(date2str, "%Y%m%d%H")
      write (scratchMessage, 31) trim(errorVar), trim(filename), trim(date2%to_iso_string())
      call allMessage(INFO, scratchMessage)

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return

99998 call allMessage(ERROR, "Unexpectedly reached end-of-file.") ! END jumps here
99999 call check_owi_err(errorIO, "Error reading Oceanweather file", filename, lun) !  ERR jumps here
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
#ifdef CMPI
      use MESSENGER, only: MSG_FINI
#endif
      implicit none
      integer, intent(in) :: iret
      character(*), intent(in) :: errorVar
      character(*), intent(in) :: filename
      integer, intent(in), optional :: unitnumber

      call setMessageSource("check_owi_err")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      if (iret /= 0) then
         if (errorVar /= "") then
            write (scratchMessage, 888) &
               trim(errorVar), trim(filename), unitnumber
888         format("Failed to '", A, "' from '", A, &
                   "' (unit number ", I3, ").")
            call allMessage(ERROR, trim(scratchMessage))
            write (6, *) 'iostat =', iret
         end if
         call allMessage(ERROR, "ADCIRC execution terminated.")
#ifdef CMPI
         Flag_ElevError = .true. !overloading this flag to help kill all mpi proc
         call msg_fini()
#endif
         stop
      end if
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
   end subroutine check_owi_err
   !-----------------------------------------------------------------------

end module mod_nws12
!-----------------------------------------------------------------------
