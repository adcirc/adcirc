!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2025 R.A. Luettich, Jr., J.J. Westerink
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
!> @brief This module handles the NWS8 forcing, which consists of a traditional
!> Holland vortex model and the CLE15 vortex model (Chavas et. al, 2015).
!>
!> @author Zachary Cobell
!> @author Coleman Blakely
!>
!> @date 2025-04-21
!>
!> @copyright Dr. R.A. Luettich and Dr. J.J. Westerink
!>
!> This module is initialized in three steps. First, the namelist data must be
!> read from an input file, typically the fort.15. Second, the additional parameters
!> that are set in the fort.15 file must be supplied. Lastly, the NWS08INIT function
!> is called.
!>
!> The namelist used for initialization is formatted as follows:
!> @code
!> &nws08Control
!>     vortexModel = "Holland" ! Holland or CLE15
!>     backgroundWindModel = "radialVelocityWeighted" ! radialVelocityWeighted or LC12
!>     BCalc = "limited" ! limited or exact
!>     thetaLatDep = .false. ! theta used in calculation of radius dependent upon latitude
!>     useInflow = .false. ! whether to use inflow angle in translation speed calc
!>     windspeed_averaging_minute = 1 ! 1 or 10
!>     w_cool = 2.0 ! magnitude of the radiative-subsidence rate in the free troposphere (Chavas and Lin, 2016)
!>     CkCd_calc = .false. ! use constant (false) or best fit (true) for CkCd
!>     CkCd = 1.0 ! ratio of exchange coefficients of enthalpy (Ck) and momentum
!>     WindMultiplier = 1.0 ! wind multiplier applied to wind speed (default = 1.0, Holland model only)
!> /
!> @endcode
module mod_nws08

   use mod_logging, only: DEBUG, ERROR, ECHO, allMessage, setMessageSource, unsetMessageSource

   implicit none

   ! Parameters used at runtime to control the vortex model
   integer, parameter :: VORTEX_MODEL_HOLLAND = 11110
   integer, parameter :: VORTEX_MODEL_CLE15 = 11111
   integer, parameter :: BACKGROUND_MODEL_RADIALVELOCITYWEIGHTED = 11120
   integer, parameter :: BACKGROUND_MODEL_LC12 = 11121
   integer, parameter :: BCALC_LIMITED = 11130
   integer, parameter :: BCALC_EXACT = 11131

   !> Weight Ratio for the time interpolation of storm track data
   real(8) :: WTRATIO = 1.0d0

   !> Boundary layer adjustment factor (from user input)
   real(8) :: BLAdj = 1.0d0

   !> Vortex model identifier (Holland or CLE15)
   integer :: vortexModelId

   !> Background model identifier (radialVelocityWeighted or LC12)
   integer :: backgroundWindModelId

   !> Allow "limited" or using "exact" calculation for Holland B parameter
   integer :: BCalcId ! use "exact" or "limited" shape

   !> Theta used in calculation of radius dependent upon latitude
   logical :: thetaLatDep ! theta used in  calculation of

   !> Boolean controlling use of inflow angle in translation speed calc
   logical :: useInflow ! whether to use inflow angle in

   !> Wind speed averaging period for sustained winds in fort.22 file. 1 or 10
   integer :: windspeed_averaging_minute = 1

   ! These next parameters are used in either the CLE15 vortex model
   !     (Chavas et. al, 2015) or the LC12 background wind model (Lin and
   !      Chavas, 2012). See those publications for more information.

   !> Magnitude of the radiative-subsidence rate in the free troposphere (Chavas and Lin, 2016)
   real(8) :: w_cool ! magnitude of the

   !> Use constant (false) or best fit (true) for CkCd
   !> If computed, then:
   !> ckcd = 0.00055*Vmax -
   !>        0.02590*Vmax +
   !>        0.76300
   !> See Chavas et al., 2015
   logical :: CkCd_calc

   !> Ratio of exchange coefficients of enthalpy (Ck) and momentum.
   !> If CkCd_calc is true, this value is computed instead
   real(8) :: CkCd ! ratio of exchange coefficients

   !> Factor of translation speed to remove from V_max
   real(8) :: backgroundWindFac ! factor of translation speed to remove from Vmax

   !> Used for storing best track time w.r.t. start of the year
   real(8), allocatable :: CastTime(:)
   character(len=4), allocatable :: CastType(:)

   !> Variables for storing best track information
   real(8), allocatable :: Lat(:), Lon(:), Spd(:), CPress(:), &
                           RRP(:), RMW(:), TVX(:), TVY(:)

   !> Variables for tracking which time stamp of the best track data to use
   integer :: bestTrackCounter = 0

   !> Variables for calculating nws = 8 winds @ nodes
   real(8), allocatable, private :: RAD(:), DX(:), DY(:), V_r(:), &
                                    THETA(:), XCOOR(:), YCOOR(:)

   !> Seconds since beginning of the year corresponding to time=0 of the simulation
   real(8) :: WindRefTime

   !> Density of water multiplied by g
   real(8) :: RhoWat0g

   !> Wind multiplier applied to wind speed
   real(8) :: WindMultiplier = 1.0d0

   private

   public :: nws08init, nws08get, getVortexModelString, readNws08Namelist, setNws08f15Parameters

contains

   !----------------------------------------------------------------
   !> Convert the vortex model string to an integer ID
   !> @param vortexModelStr The vortex model string
   !> @return The vortex model ID
   !----------------------------------------------------------------
   integer function getVortexModelId(vortexModelStr) result(id)
      use global, only: toLowercase
      use mod_logging, only: allMessage
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none
      character(*), intent(in) :: vortexModelStr

      character(len=len(vortexModelStr)) :: vortexModelStrLower
      character(1024) :: scratchMessage

      vortexModelStrLower = toLowercase(trim(vortexModelStr))

      id = -1

      select case (trim(adjustl(vortexModelStrLower)))
      case ("holland")
         id = VORTEX_MODEL_HOLLAND
      case ("cle15")
         id = VORTEX_MODEL_CLE15
      case default
         write (scratchMessage, '(3a)') "Unrecognized vortex model: '", &
            trim(vortexModelStr), "'"
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message=trim(scratchMessage))
      end select

   end function getVortexModelId

   !----------------------------------------------------------------
   !> Convert the vortex model ID to a string
   !> @return The vortex model string
   !----------------------------------------------------------------
   character(256) function getVortexModelString()
      use mod_logging, only: allMessage
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none

      character(1024) :: scratchMessage

      select case (vortexModelId)
      case (VORTEX_MODEL_HOLLAND)
         getVortexModelString = "Holland"
      case (VORTEX_MODEL_CLE15)
         getVortexModelString = "CLE15"
      case default
         write (scratchMessage, '(a,i0)') "Unrecognized vortex model ID: ", vortexModelId
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message=trim(scratchMessage))
      end select

   end function getVortexModelString

   !----------------------------------------------------------------
   !> Convert the background wind model string to an integer ID
   !> @param backgroundWindModelStr The background wind model string
   !> @return The background wind model ID
   !----------------------------------------------------------------
   integer function getBackgroundWindModelId(backgroundWindModelStr) result(id)
      use global, only: toLowercase
      use mod_logging, only: allMessage
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none
      character(*), intent(in) :: backgroundWindModelStr

      character(len=len(backgroundWindModelStr)) :: backgroundWindModelStrLower
      character(1024) :: scratchMessage

      backgroundWindModelStrLower = toLowercase(trim(backgroundWindModelStr))

      id = -1

      select case (trim(adjustl(backgroundWindModelStrLower)))
      case ("radialvelocityweighted")
         id = BACKGROUND_MODEL_RADIALVELOCITYWEIGHTED
      case ("lc12")
         id = BACKGROUND_MODEL_LC12
      case default
         write (scratchMessage, '(3a)') "Unrecognized background wind model: '", &
            trim(backgroundWindModelStr), "'"
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message=trim(scratchMessage))
      end select

   end function getBackgroundWindModelId

   !----------------------------------------------------------------
   !> Convert the background wind model ID to a string
   !> @return The background wind model string
   !----------------------------------------------------------------
   character(256) function getBackgroundWindModelString()
      use mod_logging, only: allMessage
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none

      character(1024) :: scratchMessage

      select case (backgroundWindModelId)
      case (BACKGROUND_MODEL_RADIALVELOCITYWEIGHTED)
         getBackgroundWindModelString = "radialVelocityWeighted"
      case (BACKGROUND_MODEL_LC12)
         getBackgroundWindModelString = "LC12"
      case default
         write (scratchMessage, '(a,i0)') "Unrecognized background wind model ID: ", backgroundWindModelId
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message=trim(scratchMessage))
      end select

   end function getBackgroundWindModelString

   !----------------------------------------------------------------
   !> Convert the BCalc string to an integer ID
   !> @param BCalcStr The BCalc string
   !> @return The BCalc ID
   !----------------------------------------------------------------
   integer function getBCalcId(BCalcStr) result(id)
      use global, only: toLowercase
      use mod_logging, only: allMessage
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none
      character(*), intent(in) :: BCalcStr

      character(len=len(BCalcStr)) :: BCalcStrLower
      character(1024) :: scratchMessage

      BCalcStrLower = toLowercase(trim(BCalcStr))

      id = -1

      select case (trim(adjustl(BCalcStrLower)))
      case ("exact")
         id = BCALC_EXACT
      case ("limited")
         id = BCALC_LIMITED
      case default
         write (scratchMessage, '(3a)') "Unrecognized BCalc: '", &
            trim(BCalcStr), "'"
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message=trim(scratchMessage))
      end select

   end function getBCalcId

   !----------------------------------------------------------------
   !> Convert the BCalc ID to a string
   !> @return The BCalc string
   !----------------------------------------------------------------
   character(256) function getBCalcString()
      use mod_logging, only: logMessage
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none

      character(1024) :: scratchMessage

      select case (BCalcId)
      case (BCALC_EXACT)
         getBCalcString = "exact"
      case (BCALC_LIMITED)
         getBCalcString = "limited"
      case default
         write (scratchMessage, '(a,i0)') "Unrecognized BCalc ID: ", BCalcId
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message=trim(scratchMessage))
      end select

   end function getBCalcString

   !----------------------------------------------------------------
   !> Read the NWS8 namelist from the specified unit
   !> @param iounit The unit number to read from
   !----------------------------------------------------------------
   subroutine readNws08Namelist(iounit)
      use global, only: toLowercase, logNamelistReadStatus
      use mod_logging, only: screenMessage, logMessage, WARNING
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none
      real(8), parameter :: eps = epsilon(1.0d0)
      integer, intent(in) :: iounit
      integer :: ios, io_stat
      character(len=256) :: vortexModel
      character(len=256) :: backgroundWindModel
      character(len=256) :: BCalc
      character(len=256) :: ios_nml_error_msg
      character(len=256) :: namelistSpecifier
      character(1024) :: scratchMessage

      namelist /nws08Control/ &
         vortexModel, &
         backgroundWindModel, &
         BCalc, &
         thetaLatDep, &
         useInflow, &
         windspeed_averaging_minute, &
         w_cool, &
         CkCd_calc, &
         CkCd, &
         WindMultiplier

      ! Check if the specified unit is open for read
      inquire (UNIT=iounit, IOSTAT=io_stat)
      if (io_stat /= 0) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message="fort.15 file is not open for read.")
      end if

      ! CPB: add nws8Control namelist defaults
      vortexModel = 'Holland'
      backgroundWindModel = 'radialVelocityWeighted'
      BCalc = 'limited'
      thetaLatDep = .false.
      useInflow = .false.
      windspeed_averaging_minute = 1
      w_cool = 2d0
      CkCd_calc = .false.
      CkCd = 1d0
      WindMultiplier = 1.0d0

      ! read the namelist
      namelistSpecifier = 'nws08Control'
      read (unit=iounit, nml=nws08Control, iostat=IOS, iomsg=ios_nml_error_msg)
      call logNameListReadStatus(namelistSpecifier, ios, ios_nml_error_msg)

      write (scratchMessage, '(a,a)') "vortexModel = ", trim(vortexModel)
      call logMessage(ECHO, trim(scratchMessage))
      write (scratchMessage, '(a,a)') "backgroundWindModel = ", trim(backgroundWindModel)
      call logMessage(ECHO, trim(scratchMessage))
      write (scratchMessage, '(a,a)') "BCalc = ", trim(BCalc)
      call logMessage(ECHO, trim(scratchMessage))
      write (scratchMessage, '(a,l1)') "thetaLatDep = ", thetaLatDep
      call logMessage(ECHO, trim(scratchMessage))
      write (scratchMessage, '(a,l1)') "useInflow = ", useInflow
      call logMessage(ECHO, trim(scratchMessage))
      write (scratchMessage, '(a,i0)') "windspeed_averaging_minute = ", windspeed_averaging_minute
      call logMessage(ECHO, trim(scratchMessage))
      write (scratchMessage, '(a,f6.2)') "WindMultiplier = ", WindMultiplier
      call logMessage(ECHO, trim(scratchMessage))

      vortexModelId = getVortexModelId(vortexModel)
      backgroundWindModelId = getBackgroundWindModelId(backgroundWindModel)
      bCalcId = getBCalcId(BCalc)

      if (WindMultiplier - 1.0d0 > eps .and. vortexModelId == VORTEX_MODEL_CLE15) then
         write (scratchMessage, '(a)') "WindMultiplier is only used with the Holland vortex model."
         call screenMessage(WARNING, trim(scratchMessage))
      end if

      ! only print these next lines if we are using the CLE15 vortex
      ! model since they are not used in regular Holland
      if (vortexModelId == VORTEX_MODEL_CLE15) then
         write (scratchMessage, '(a,f6.2,a)') "w_cool  =", w_cool, " mm s^-1"
         call logMessage(ECHO, trim(scratchMessage))
         if (CkCd_calc .eqv. .true.) then
            write (scratchMessage, '(a,a)') "CkCd will be calculated", &
               " using the relationship: CkCd = 0.00055*V^2 - 0.02590*V + 0.763"
            call logMessage(ECHO, trim(scratchMessage))
         else
            write (scratchMessage, '(a,f6.2)') "CkCd  =", CkCd
            call logMessage(ECHO, trim(scratchMessage))
         end if
      end if

      rewind (iounit)

      ! convert w_cool from mm per second to m per second
      w_cool = w_cool/1000d0

   end subroutine readNws08Namelist

   !----------------------------------------------------------------
   !> @brief Set the NWS8 parameters
   !> @param WindRefTime_in The reference time for the wind
   !> @param BLAdj_in The boundary layer adjustment factor
   !> @param g The acceleration due to gravity
   !----------------------------------------------------------------
   subroutine setNws08f15Parameters(WindRefTime_in, BLAdj_in, g)
      use ADC_CONSTANTS, only: rhowat0
      implicit none
      real(8), intent(in) :: WindRefTime_in
      real(8), intent(in) :: BLAdj_in
      real(8), intent(in) :: g

      ! Set up input variables
      WindRefTime = WindRefTime_in
      BLAdj = BLAdj_in
      RhoWat0g = rhowat0*g

   end subroutine setNws08f15Parameters

! ----------------------------------------------------------------
!> @brief Initializes the NWS08 module
!> @param timeloc The current time of the simulation
! ----------------------------------------------------------------
   subroutine NWS08INIT(timeloc)
      use MESH, only: SLAM, SFEA, NP
      use ADC_CONSTANTS, only: RAD2DEG
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE

      implicit none

      real(8), intent(in) :: TIMELOC
      integer :: i

      call setMessageSource("NWS08INIT")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      ! Allocate necessary arrays
      allocate (RAD(NP), DX(NP), DY(NP), XCOOR(NP), YCOOR(NP))
      allocate (V_r(NP), THETA(NP))
      XCOOR = SLAM*RAD2DEG
      YCOOR = SFEA*RAD2DEG

      ! Read in ATCF Best Track/Objective Aid/Wind Radii Format winds
      call readBestTrackData()

      ! Determine the correspondence between the current simulation time and
      ! the fort.22 file.
      i = 2
      do
         if (dble(ceiling(TIMELOC)) >= CastTime(i - 1) .and. dble(floor(TIMELOC)) < CastTime(i)) then
            bestTrackCounter = i
            exit
         else
            i = i + 1
            if (i > size(CastTime)) then
               call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                              message="The Storm Hindcast/Forecast Input File (unit 22) " &
                              //"does not contain times/dates that correspond " &
                              //"to the ADCIRC current model time. " &
                              //"ADCIRC terminating.")
            end if
         end if
      end do

      ! Get factor of forward speed to remove from Vmax
      select case (backgroundWindModelId)
      case (BACKGROUND_MODEL_RADIALVELOCITYWEIGHTED)
         backgroundWindFac = 1.0d0
      case (BACKGROUND_MODEL_LC12)
         backgroundWindFac = 0.62d0
      end select

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()

   end subroutine NWS08INIT

! ----------------------------------------------------------------
!> @brief Computes the wind velocity and pressure at each node using
!> the Holland or CLE15 vortex models
!> @author Coleman Blakely
!> @date 2024-02-01
!>
!> @param[inout] WVNX The x-component of thewind velocity at each node
!> @param[inout] WVNY The y-component of the wind velocity ateach node
!> @param[inout] PRESS The pressure ateach node
!> @param[in] TIMELOC The current time of the simulation
!> @param[inout] FoundEye A logical variable indicating if the eye was found
!> @param[inout] EyeLon The longitude of the storm eye
!> @param[inout] EyeLat The latitude of the storm eye
!>
!> Calls the appropriate subroutines to calculate winds using
!> vortex models that take as input track data in ATCF Best
!> Track/Objective Aid/Wind Radii format. Also applies the
!> appropriate background wind field. Options are controlled by the
!> nws8Control namelist. Current options are:
!>
!> Vortex models:
!>    - Holland Model
!>    - CLE15 Model (Chavas et. al., 2015)
!>
!> Background wind models:
!>    - None (standard Holland)
!>    - LC12 (Lin and Chavas, 2012)
!>
! ----------------------------------------------------------------
   subroutine NWS08GET(WVNX, WVNY, PRESS, TIMELOC, FoundEye, EyeLon, EyeLat)
      use MESH, only: NP
      implicit none
      real(8), intent(IN) :: TIMELOC
      real(8), intent(INOUT), dimension(NP) :: WVNX, WVNY
      real(8), intent(INOUT), dimension(NP) :: PRESS
      logical, intent(INOUT) :: FoundEye
      real(8), intent(INOUT), dimension(3) :: EyeLon, EyeLat

      call setMessageSource("NWS08GET")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif
      select case (vortexModelId)
      case (VORTEX_MODEL_HOLLAND)
         call HollandGet(WVNX, WVNY, PRESS, TIMELOC, FoundEye, EyeLon, EyeLat)
      case (VORTEX_MODEL_CLE15)
         call CLE15GET(WVNX, WVNY, PRESS, TIMELOC, FoundEye, EyeLon, EyeLat)
      end select

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
   end subroutine NWS08GET

!-----------------------------------------------------------------
!> @brief Calculate wind velocity at nodes from the Holland model
!>
!> @param[inout] WVNX The x-component of thewind velocity at each node
!> @param[inout] WVNY The y-component of the windvelocity ateach node
!> @param[in] PRESS The pressureateach node
!> @param[inout] TIMELOC The current time of the simulation
!> @param[inout] FoundEye A logical variable indicatingif the eye was found
!> @param[inout] EyeLon Thelongitude of the storm eye
!> @param[inout] EyeLat The latitude of thestorm eye
!>
!> Subroutine to calculate wind velocity at nodes from
!> the Holland Wind model.
!>
!> The format statement takes into account whether the track data is
!> hindcast/nowcast (BEST) or forecast (OFCL).
!>
!> The first line in the file MUST be a hindcast, since the central
!> pressure and the RMW are carried forward from hindcasts into
!> forecasts. So there needs to be at least one hindcast to carry the data
!> forward.
!>
!> Assumes spherical coordinates (ICS=2 in fort.15 file).
! ----------------------------------------------------------------
   subroutine HollandGet(WVNX, WVNY, PRESS, TIMELOC, FoundEye, EyeLon, EyeLat)
      use MESH, only: NP
      use ADC_CONSTANTS, only: RHOWAT0, G, mb2pa, DEG2RAD, prbckgrnd, one2ten, e, rhoAir
      use GLOBAL, only: CORIF, sphericalDistance
      implicit none
      real(8), intent(in) :: TIMELOC
      real(8), intent(out), dimension(NP) :: WVNX, WVNY
      real(8), intent(out), dimension(NP) :: PRESS
      real(8), intent(inout), dimension(3) :: EyeLon, EyeLat
      logical, intent(inout) :: FoundEye
      integer :: I
      real(8) :: TVX, TVY, RRP, RMW, B
      real(8) :: TransSpdX, TransSpdY
      real(8) :: cpress, lon, lat, spd, rrad, alpha
      real(8) :: ts ! storm translation speed, m/s
      real(8) :: centralPressureDeficit ! difference btw ambient and cpress
      real(8) :: inflowAngle
      real(8) :: non_cyclostrophic_part

      call setMessageSource("HollandGet")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      b = huge(1.0d0)

      ! Get data for this time step.
      call GetHollandStormData(lat, lon, cpress, spd, rrp, rmw, tvx, tvy, TIMELOC)
      !
      ! @jasonfleming: If this is a "CALM" period, set winds to zero
      ! velocity and pressure equal to the background pressure and return.
      if (cpress < 0.d0) then
         press(:) = PRBCKGRND*mb2pa/rhoWat0g ! convert mb to mH2O
         wvnx(:) = 0.d0
         wvny(:) = 0.d0
         EyeLat(1) = EyeLat(2)
         EyeLon(1) = EyeLon(2)
         EyeLat(2) = EyeLat(3)
         EyeLon(2) = EyeLon(3)
         EyeLat(3) = lat
         EyeLon(3) = lon
         FoundEye = .true.
#if defined(WIND_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG, "Return.")
#endif
         call unsetMessageSource()
         return
      end if

      ! jgf50.32: If we are using sector-based wind drag, record the location
      ! of the center of the storm.
      call update_storm_eye_position(lon, lat, EyeLon, EyeLat, FoundEye)

      ! Calculate and limit central pressure deficit; some track files
      ! (e.g., Charley 2004) may have a central pressure greater than the
      ! ambient pressure that this subroutine assumes
      ! jgf51.14: Apply a factor of 100.0 to convert PRBCKGRND from
      ! mb to pa for use in this subroutine.
      centralPressureDeficit = PRBCKGRND*100.d0 - cpress
      if (centralPressureDeficit < 100.d0) then
         centralPressureDeficit = 100.d0
      end if

      ! jgf46.29 Subtract the translational speed of the storm from the
      ! observed max wind speed to avoid distortion in the Holland curve
      ! fit. The translational speed will be added back later.
      ts = sqrt(tvx*tvx + tvy*tvy)
      spd = spd - backgroundWindFac*ts

      ! Convert wind speed from 10 meter altitude (which is what the
      ! NHC forecast contains) to wind speed at the top of the atmospheric
      ! boundary layer (which is what the Holland curve fit requires).
      spd = spd/BLAdj

      ! Calculate Holland parameters and limit the result to its appropriate range.
      if (BCalcId == BCALC_LIMITED) then
         B = rhoAir*e*(spd**2.d0)/(centralPressureDeficit)
         if (B < 1.0d0) B = 1.0d0
         if (B > 2.5d0) B = 2.5d0
         non_cyclostrophic_part = 0d0
      end if

      ! WJP convert RMW to meters
      RMW = RMW*1000d0

      ! Calculate wind velocity and pressure at each node.
      do I = 1, NP
         DX(I) = (XCOOR(I) - lon)*DEG2RAD
         DY(I) = (YCOOR(I) - lat)*DEG2RAD
         if (thetaLatDep .eqv. .true.) then
            THETA(I) = atan2(DY(I), (cos(lat*DEG2RAD)*DX(I)))
         else
            THETA(I) = atan2(DY(I), DX(I))
         end if

         ! Compute the distances based on haversine formula for distance along a sphere
         rad(i) = sphericalDistance(DX(I), DY(I), lat, ycoor(i))

         ! Add option to calculate exact B parameter
         if (BCalcId == BCALC_EXACT) then
            non_cyclostrophic_part = RMW*CORIF(I)
            B = rhoAir*e*spd*(spd + non_cyclostrophic_part)/(centralPressureDeficit)
         end if

         ! Compute pressure field
         alpha = (RMW/RAD(I))**B

         !P.V 11/04/2022 to account for Southern Hemisphere
         PRESS(I) = cpress + centralPressureDeficit*exp(-alpha)

         ! Compute velocity field with absolute around
         ! CORIOLIS for the Southern Hempisphere
         V_r(I) = sqrt(alpha*exp(1.d0 - alpha)*spd*(spd + non_cyclostrophic_part) &
                       + 0.25d0*(RAD(I)**2.d0)*(CORIF(I)**2.d0) &
                       ) - &
                  0.5d0*RAD(I)*abs(CORIF(I))

         ! Apply wind speed mutliplier
         V_r(I) = V_r(I)*WindMultiplier

         ! Convert wind velocity from top of atmospheric boundary layer to
         ! 10m wind velocity.
         V_r(I) = V_r(I)*BLADj

         ! Find the velocity components.
         rrad = RAD(I)/RMW
         inflowAngle = calcInflowAngle(rrad)

         !P.V 11/04/2022 to account for Southern Hemisphere
         WVNX(I) = sign(1.0d0, -lat)*V_r(I)*sin(THETA(I) + sign(1.0d0, lat)*inflowangle)
         WVNY(I) = sign(1.0d0, lat)*V_r(I)*cos(THETA(I) + sign(1.0d0, lat)*inflowangle)

         ! jgf46.31 Determine translation speed that should be added to final
         ! storm wind speed. This is tapered to zero as the storm wind tapers
         ! to zero toward the eye of the storm and at long distances from the
         ! storm.
         call calcTranslationSpeed(V_r(I), spd, TVX, TVY, lat, rrad, TransSpdX, TransSpdY)

         ! jgf46.31 Add the storm translation speed.
         WVNX(I) = WVNX(I) + TransSpdX
         WVNY(I) = WVNY(I) + TransSpdY

         ! jgf46.21 Also convert from 1 minute averaged winds to 10
         ! minute averaged winds (0.88 factor).
         if (windspeed_averaging_minute == 1) then
            WVNX(I) = WVNX(I)*one2ten
            WVNY(I) = WVNY(I)*one2ten
         end if

         !P.V 11/04/2022
         PRESS(I) = max(0.85d5, min(1.1d5, PRESS(I))) !Typhoon Tip 870 hPa ... 12-oct-1979
         WVNX(I) = max(-200.d0, min(200.d0, WVNX(I)))
         WVNY(I) = max(-200.d0, min(200.d0, WVNY(I)))

         !-------------------------------------------
         ! Convert atmospheric pressure (Pascals) to
         ! atmospheric pressure-induced water surface
         ! elevation (meters).
         !-------------------------------------------
         PRESS(I) = PRESS(I)/(RHOWAT0*G)

      end do
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
   end subroutine HollandGet

!-----------------------------------------------------------------------
!> @brief Calculate the translation speed of the storm
!>
!> @param[in] Wind velocity at the specified radius
!> @param[in] Storm speed
!> @param[in] TVX X-component of the storm translation velocity
!> @param[in] TVY Y-component of the storm translation velocity
!> @param[in] Latitude of the storm center
!> @param[in] rrad Current radius divided by the radius to max wind speed
!> @param[out] TransSpdX X-component of the translation speed at the specified radius
!> @param[out] TransSpdY Y-component of the translation speed at the specified radius
!-----------------------------------------------------------------------
   subroutine calcTranslationSpeed(V_r, spd, TVX, TVY, lat, rrad, TransSpdX, TransSpdY)
      use ADC_CONSTANTS, only: deg2rad

      implicit none

      real(8), intent(IN) :: V_r, spd, TVX, TVY, lat, rrad
      real(8), intent(OUT) :: TransSpdX, TransSpdY
      real(8) :: beta, alpha

      call setMessageSource("calcTranslationSpeed")
#ifdef ALL_TRACE
      call allMessage(DEBUG, "Enter.")
#endif
      select case (backgroundWindModelId)
      case (BACKGROUND_MODEL_RADIALVELOCITYWEIGHTED)
         ! translation speed weighted by the relative radial velocity
         TransSpdX = (abs(V_r)/spd)*TVX
         TransSpdY = (abs(V_r)/spd)*TVY
      case (BACKGROUND_MODEL_LC12)
         ! Add (0.3-0.67)vt to the left of the direction of forward motion
         ! and use direction change of (0-20)deg clockwise; note the hemisphere
         ! Get the ts factor, alpha in LC12
         if (rrad <= 1d0) then
            alpha = 0.3d0 + (0.62d0 - 0.3d0)*rrad
         elseif (rrad <= 1.35d0) then
            alpha = 0.62d0 + (0.67d0 - 0.62d0)*(rrad - 1d0)/0.35d0
         elseif (rrad <= 2.2d0) then
            alpha = 0.67d0 + (0.6d0 - 0.67d0)*(rrad - 1.35d0)/0.85d0
         elseif (rrad <= 8d0) then
            alpha = 0.6d0 + (0.5d0 - 0.6d0)*(rrad - 2.2d0)/5.8d0
         elseif (rrad <= 12d0) then
            alpha = 0.5d0 + (0.55d0 - 0.5d0)*(rrad - 8d0)/4d0
         else
            alpha = 0.55d0
         end if
         ! Get the rotation factor, beta in LC12
         if (rrad <= 0.8d0) then
            beta = 0.0d0 + (20d0 - 0.0d0)*rrad/0.8d0
         elseif (rrad <= 1.25d0) then
            beta = 20d0 + (10d0 - 20d0)*(rrad - 0.8d0)/0.45d0
         elseif (rrad <= 1.5d0) then
            beta = 10d0
         elseif (rrad <= 3d0) then
            beta = 10d0 + (15d0 - 10d0)*(rrad - 1.5d0)/1.5d0
         elseif (rrad <= 8d0) then
            beta = 15d0 + (20d0 - 15d0)*(rrad - 3d0)/5d0
         else
            beta = 20d0
         end if
         beta = deg2rad*beta
         TransSpdX = alpha*(TVX*cos(beta) &
                            - sign(1.0d0, lat)*TVY*sin(beta))
         TransSpdY = alpha*(TVY*cos(beta) &
                            + sign(1.0d0, lat)*TVX*sin(beta))
      end select
#ifdef ALL_TRACE
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()

   end subroutine calcTranslationSpeed

!-----------------------------------------------------------------------
!> @brief Calculate the inflow angle based on the radial distance
!> @param[in] RRAD Relative radial distance (fraction of RMW)
!> @return inflowAngle The inflow angle
!-----------------------------------------------------------------------
   real(8) pure function calcInflowAngle(RRAD) result(inflowAngle)
      implicit none
      real(8), intent(IN) :: RRAD ! relative radial distance (fraction of RMW)

      if (useInflow .eqv. .true.) then
         ! Adjustment to account for the inflow toward the center
         ! Bretschneider, 1972. New by Ning
         if (RRAD < 1d0) then
            inflowangle = 0.1745d0*(1.d0 + RRAD)
         elseif (RRAD < 1.2d0) then
            inflowangle = 0.3491d0 + 0.4363d0*(RRAD - 1.0d0)
         else
            inflowangle = 0.4363d0
         end if
         !     original ADCIRC uses inflowangle=0
      else
         inflowAngle = 0d0
      end if

   end function calcInflowAngle

! ----------------------------------------------------------------
!> @brief Returns the storm data by interpolating in time from the data
!> read in from the best track file (fort.22) previously in readBestTrackData
!>
!> @param[out] LatOut The latitude of the storm center
!> @param[out] LonOut The longitude of the storm center
!> @param[out] CPressOut The central pressure of the storm
!> @param[out] SpdOut The maximum wind speed of the storm
!> @param[out] RRPOut The radius of the last closed isobar
!> @param[out] RMWOut The radius of maximum wind speed
!> @param[out] TVXOut The x-component of thestorm translation velocity
!> @param[out] TVYOut The y-component of thestorm translation velocity
!> @param[in] TIMELOC Thecurrent time of the simulation
! ----------------------------------------------------------------
   subroutine GetHollandStormData(LatOut, LonOut, CPressOut, SpdOut, &
                                  RRPOut, RMWOut, TVXOut, TVYOut, TIMELOC)
      use VORTEX, only: uvtrans
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none
      real(8), intent(in) :: TIMELOC
      real(8), intent(out) :: LatOut, LonOut, CPressOut
      real(8), intent(out) :: SpdOut, RRPOut, RMWOut, TVXOut, TVYOut

      integer :: i ! Current array counter for fort.22 file
      call setMessageSource("getHollandStormData")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif
      i = bestTrackCounter

      ! If time exceeds the next hindcast/nowcast/forecast time, increment the
      ! array counter.
      if (TIMELOC > CastTime(i)) then
         i = i + 1
         ! jgf51.14: Check to see that we haven't gone off the end of
         ! meteorological data.
         if (i > size(CastTime)) then
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message='The simulation time has extended beyond the end '// &
                           'of the meteorological dataset.')
         end if
      end if

      ! Interpolate w.r.t. time
      ! @jasonfleming: Add CALM handling.
      LonOut = Lon(i - 1) + WTRATIO*(Lon(i) - lon(i - 1))
      LatOut = Lat(i - 1) + WTRATIO*(Lat(i) - lat(i - 1))
      if ((trim(castType(i)) /= 'CALM') .and. trim(castType(i - 1)) /= 'CALM') then
         WTRATIO = (TIMELOC - CastTime(i - 1))/(CastTime(i) - CastTime(i - 1))
         CPressOut = CPress(i - 1) + WTRATIO*(CPress(i) - CPress(i - 1))
         SpdOut = Spd(i - 1) + WTRATIO*(Spd(i) - Spd(i - 1))
         RRPOut = RRP(i - 1) + WTRATIO*(RRP(i) - RRP(i - 1))
         RMWOut = RMW(i - 1) + WTRATIO*(RMW(i) - RMW(i - 1))
         TVXOut = TVX(i - 1) + WTRATIO*(TVX(i) - TVX(i - 1))
         TVYOut = TVY(i - 1) + WTRATIO*(TVY(i) - TVY(i - 1))
      else
         CPressOut = -99999.d0
         SpdOut = -99999.d0
         RRPOut = -99999.d0
         RMWOut = -99999.d0
         TVXOut = -99999.d0
         TVYOut = -99999.d0
      end if

      bestTrackCounter = i

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()

      return

   end subroutine GetHollandStormData

!----------------------------------------------------------------
!> @brief Reads the best trackdata from the fort.22 file as
!> ATCF Best Track/Objective Aid/Wind Radii Format
!----------------------------------------------------------------
   subroutine readBestTrackData()
      use SIZES, only: GBLINPUTDIR
      use GLOBAL, only: RNDAY, timeconv
      use mod_io, only: openFileForRead
      use VORTEX, only: uvtrans
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none
      integer, allocatable :: iYear(:), iMth(:), iDay(:), iHr(:)
      integer, allocatable :: iLat(:), iLon(:)
      character(1), allocatable :: ns(:), ew(:)
      integer, allocatable :: iSpd(:), iCPress(:), iRRP(:), iRMW(:)
      integer, allocatable :: iFcstInc(:) ! hours
      real(8), allocatable :: FcstInc(:) ! seconds
      integer :: iNowcastCPress, iNowcastRRP, iNowcastRMW
      integer :: i ! Current array counter for fort.22 file
      integer :: nl ! Number of lines in the fort.22 file
      integer :: pl ! populated length of Holland Data array
      integer :: ios ! return code for an i/o operation
      call setMessageSource("Read_Best_Track_Data")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      pl = huge(1)

      ! Determine the number of lines in the file.
      nl = 0
      call openFileForRead(22, trim(GBLINPUTDIR)//'/'//'fort.22', ios)
      if (ios > 0) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message="The symmetric vortex parameter file was not found. " &
                        //"ADCIRC terminating.")
      end if

      do
         read (UNIT=22, FMT='(A170)', end=8888)
         nl = nl + 1
      end do
8888  rewind (22)

      ! Dimension the arrays according to the number of lines in the file,
      ! this will be greater than or equal to the size of the array we need
      ! (probably greater because of the repeated lines that we throw away)
      allocate (iYear(nl), iMth(nl), iDay(nl), iHr(nl), iLat(nl), iLon(nl), &
                iSpd(nl), iCpress(nl), iRRP(nl), iRMW(nl), iFcstInc(nl))
      allocate (ns(nl), ew(nl))
      allocate (Lat(nl), Lon(nl), Spd(nl), CPress(nl), RRP(nl), RMW(nl), FcstInc(nl), TVX(nl), TVY(nl))
      allocate (CastType(nl))
      allocate (CastTime(nl))

      ! Now read the data into the arrays. The first
      ! line must be a hindcast/nowcast.
      i = 1
      do
         ! Get another line of data from the file and check to see if the
         ! line represents a new point in time, or is a repeated time
         ! point. Repeated time points occur in hindcasts for the purpose of
         ! describing winds in the quadrants of the storm. We don't use the
         ! quadrant-by-quadrant wind data. Repeated time data occur in the
         ! forecast because the time data is just the time that the forecast
         ! was made. The important parameter in the forecast file is the
         ! forecast increment.
         read (UNIT=22, FMT=228, end=9999) &
            iYear(i), iMth(i), iDay(i), iHr(i), &
            CastType(i), iFcstInc(i), iLat(i), ns(i), iLon(i), ew(i), &
            iSpd(i), iCPress(i), iRRP(i), iRMW(i)

         ! yr,mo,dy,hr, ,type, inc,  lat,NS,  lon,EW,  spd,   pc,
228      format(8x, i4, i2, i2, i2, 6x, a4, 2x, i3, 1x, i4, a1, 2x, i4, a1, 2x, i3, 2x, i4, 47x, i3, 2x, i3)

         select case (trim(CastType(i)))
         case ("BEST") ! nowcast/hindcast
            !     Check to see if this is a repeated line. If so, go directly to the
            !     next line without any processing.
            if (i > 1) then
               if (iYear(i) == iYear(i - 1) .and. &
                   iMth(i) == iMth(i - 1) .and. &
                   iDay(i) == iDay(i - 1) .and. &
                   iHr(i) == iHr(i - 1)) then
                  cycle
               end if
            end if
            ! Save the central pressure, radius of last closed isobar, and
            ! radius to max wind for use in forecasts
            iNowcastCPress = iCPress(i)
            iNowcastRMW = iRMW(i)
            iNowcastRRP = iRRP(i)

            ! Determine the time of this hindcast in seconds since the beginning
            ! of the year.
            CastTime(i) = TimeConv(iYear(i), iMth(i), iDay(i), iHr(i), 0, 0.d0)

            ! Determine the CastTime in seconds since the beginning of the simulation.
            CastTime(i) = CastTime(i) - WindRefTime
            FcstInc(i) = dble(iFcstInc(i))

         case ("OFCL") ! forecast
            !  Check to see if this is a repeated line (i.e., a forecast that
            !  coincides with the nowcast, or a repeated forecast). If so, go
            !  directly to the next line without any processing.
            if (i > 1) then
               if ((iFcstInc(i) == 0 .and. &
                    (iYear(i) == iYear(i - 1) .and. &
                     iMth(i) == iMth(i - 1) .and. &
                     iDay(i) == iDay(i - 1) &
                     .and. iHr(i) == iHr(i - 1))) &
                   .or. (iFcstInc(i) /= 0 .and. iFcstInc(i) == iFcstInc(i - 1))) cycle
            end if
            FcstInc(i) = dble(iFcstInc(i))

            ! Determine the time of this forecast in seconds since the beginning
            ! of the year.
            if (iFcstInc(i) == 0) then
               CastTime(i) = TimeConv(iYear(i), iMth(i), iDay(i), iHr(i), 0, 0.d0)
               CastTime(i) = CastTime(i) - WindRefTime
            else
               FcstInc(i) = FcstInc(i)*3600.d0 ! convert hours to seconds
               CastTime(i) = CastTime(i - 1) + (FcstInc(i) - FcstInc(i - 1))
            end if

            ! jgf48.4637 If the forecast values of central pressure or RMW are
            ! zero (and they will be if unless the user has filled them in, because
            ! the NHC does not forecast these parameters), exit with a fatal error.
            if ((iCPress(i) == 0) .or. (iRMW(i) == 0)) then
               call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                              message='The storm hindcast/forecast input file (unit 22) '// &
                              'contains invalid data for central pressure or Rmax.')
            end if

            ! @jasonfleming: Adding a new type to allow the analyst to add lines
            ! that do nothing but produce zero winds and background barometric
            ! pressure. These lines can have a date/time like a BEST line or
            ! a date/time and forecast period like an OFCL line.
         case ("CALM")
            call allMessage(ECHO, 'The fort.22 file contains at least one "CALM" line.')
            if (i > 1) then
               if ((iFcstInc(i) == 0 .and. (iYear(i) == iYear(i - 1) .and. iMth(i) == iMth(i - 1) &
                                            .and. iDay(i) == iDay(i - 1) .and. iHr(i) == iHr(i - 1))) .or. &
                   (iFcstInc(i) /= 0 .and. iFcstInc(i) == iFcstInc(i - 1))) cycle
            end if
            FcstInc(i) = dble(iFcstInc(i))

            if (iFcstInc(i) == 0) then
               CastTime(i) = TimeConv(iYear(i), iMth(i), iDay(i), iHr(i), 0, 0.d0)
               CastTime(i) = CastTime(i) - WindRefTime
            else
               FcstInc(i) = FcstInc(i)*3600.d0 ! convert hours to seconds
               CastTime(i) = CastTime(i - 1) + (FcstInc(i) - FcstInc(i - 1))
            end if

         case DEFAULT ! unrecognized
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message='Only "BEST", "OFCL", or "CALM" are allowed '// &
                           'in the 5th column of fort.22.')
         end select

         ! Convert integers to reals.
         Lat(i) = dble(iLat(i))
         Lon(i) = dble(iLon(i))
         Spd(i) = dble(iSpd(i))
         CPress(i) = dble(iCPress(i))
         RRP(i) = dble(iRRP(i))
         RMW(i) = dble(iRMW(i))

         ! Convert units.
         Lat(i) = Lat(i)/10.d0 ! convert 10ths of degs to degs
         Lon(i) = Lon(i)/10.d0 ! convert 10ths of degs to degs
         if (ew(i) == 'W') then
            Lon(i) = -1.d0*Lon(i) ! negative b/c WEST longitude
         end if
         if (ns(i) == 'S') then
            Lat(i) = -1.d0*Lat(i) ! negative b/c SOUTH latitude
         end if
         CPress(i) = CPress(i)*100.d0 ! convert mbar to Pa
         RRP(i) = RRP(i)*1.852000003180799d0*1000.0d0 ! convert nm to m
         RMW(i) = RMW(i)*1.852000003180799d0 ! convert nm to km
         Spd(i) = Spd(i)*0.51444444d0 ! convert kts to m/s

         ! Save the number of non-repeated lines from the fort.22 file, this
         ! is the populated length of the array.
         pl = i

         !Increment array counter
         i = i + 1

      end do
9999  close (22)

      ! Calculate storm translation velocities based on change in position, then
      ! convert degrees/time to m/s
      ! initialize storm translation speeds
      TVX = 0.0d0
      TVY = 0.0d0
      do i = 2, pl
         ! Calculate storm translation velocities based on change in position,
         ! approximate u and v translation velocities
         call uvTrans(lat(i - 1), lon(i - 1), &
                      lat(i), lon(i), &
                      CastTime(i - 1), CastTime(i), &
                      tvx(i), tvy(i))

      end do
      ! extrapolate later storm translation speeds back to the
      ! first storm position
      if (pl >= 2) then
         TVX(1) = TVX(2)
         TVY(1) = TVY(2)
      end if
      !
      ! @jasonfleming: Check to see if there is enough data to cover
      ! the whole run and bomb out immediately if there isn't.
      if (castTime(pl) < RNDAY*86400.d0) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message='The fort.22 file ends before RNDAY.')
      end if

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()

      return
   end subroutine readBestTrackData

! ----------------------------------------------------------------
!> @brief Calculate the wind velocity ateach node based on the CLE15 wind model
!>
!> @author Coleman Blakely
!> @date Feb. 2024
!>
!> @param[in] WVNX The x-component of the wind velocity
!> @param[in] WVNY The y-component of thewind velocity
!> @param[in] PRESS Theatmospheric pressure
!> @param[in] TIMELOC The current time ofthe simulation
!> @param[in] FoundEye Flag indicating if the eye of the stormwas found
!> @param[in] EyeLon The longitude of thestorm center
!> @param[in] EyeLat The latitudeof thestorm center
!>
!> Calculates wind velocity at nodes based on the wind model
!> described in Chavas et. al., 2015. This model consists of an
!> inner vortex model described by:
!>
!>   (M_i/M_m)^(2-(C_k/C_d)) =
!>      (2(r/r_m)^2)/(2-(C_k/C_d)+(C_k/C_d)(r/r_m)^2)        (1)
!>
!>   M_m = (r_m)*(V_m) + (1/2)*(f)*(r_m)^2                   (2)
!> Where,
!>   M_i = angular momentum at radius r
!>   M_m = angular momentum at radius of max winds
!>   C_k/C_d = ratio of exchange coefficients of enthalpy and
!>             momentum
!>   r = radius
!>   r_m = radius of max winds
!>   V_m = max wind speed
!>
!> As well as an outer vortex model described by:
!>
!>   (dM_o/dr) = gamma*((r*V)^2)(r_0^2-r^2)                  (3)
!>
!>   gamma = (2*C_d)/w_cool                                  (4)
!> Where,
!>   M_o = angular momentum at radius r
!>   C_d = surface drag coefficient
!>   w_cool = free-tropospheric air subsidence rate
!>   r_0 = radius of diminishing winds (@r=r_0, V=0)
!>
!> Note that angular momentum is calculated as M = r*V+(1/2)*f*r^2
!>
!> To find the complete wind profile a radius r_a is found where the
!> angular momentum of the inner model and outer model as well as
!> their derivatives match. Rather than solve for the angular
!> momentum, this subroutine instead iteratively solves for V in both
!> the inner and outer models with given r_m, V_m and estimates of
!> r_0 until M_o(r_m,V_m,r_0,r_a) = M_i(r_m,V_m,r_0,r_a) = M_a.
!>
!> Chavas et. al., 2015: https://doi.org/10.1175/JAS-D-15-0014.1
!>
! ----------------------------------------------------------------
   subroutine CLE15GET(WVNX, WVNY, PRESS, TIMELOC, FoundEye, EyeLon, EyeLat)
      use MESH, only: NP
      use ADC_CONSTANTS, only: RHOWAT0, G, mb2pa, deg2rad, e, one2ten, prbckgrnd, rhoAir
      use GLOBAL, only: CORIF, sphericalDistance

      implicit none

      real(8), intent(in) :: TIMELOC
      real(8), intent(out), dimension(NP) :: WVNX, WVNY
      real(8), intent(out), dimension(NP) :: PRESS
      real(8), intent(inout), dimension(3) :: EyeLon, EyeLat
      logical, intent(inout) :: FoundEye
      integer :: I
      real(8) :: TVX, TVY, RRP, RMW, B
      real(8) :: TransSpdX, TransSpdY
      real(8) :: cpress, lon, lat, spd, rrad
      real(8) :: ts ! storm translation speed, m/s
      real(8) :: centralPressureDeficit ! difference btw ambient and cpress
      real(8) :: inflowAngle
      real(8), allocatable, dimension(:) :: V_cle15
      real(8) :: r0

      call setMessageSource("CLE15Get")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      b = huge(1d0)

      ! Get data for this time step.
      call GetHollandStormData(lat, lon, cpress, spd, rrp, rmw, tvx, tvy, TIMELOC)

      ! If this is a "CALM" period, set winds to zero
      ! velocity and pressure equal to the background pressure and return.
      if (cpress < 0.d0) then
         press(:) = PRBCKGRND*mb2pa/rhoWat0g ! convert mb to mH2O
         wvnx(:) = 0.d0
         wvny(:) = 0.d0
         EyeLat(1) = EyeLat(2)
         EyeLon(1) = EyeLon(2)
         EyeLat(2) = EyeLat(3)
         EyeLon(2) = EyeLon(3)
         EyeLat(3) = lat
         EyeLon(3) = lon
         FoundEye = .true.
#if defined(WIND_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG, "Return.")
#endif
         call unsetMessageSource()
         return
      end if

      ! If we are using sector-based wind drag, record the location
      ! of the center of the storm.
      call update_storm_eye_position(lat, lon, EyeLat, EyeLon, FoundEye)

      ! Calculate and limit central pressure deficit; some track files
      ! (e.g., Charley 2004) may have a central pressure greater than the
      ! ambient pressure that this subroutine assumes
      ! Apply a factor of 100.0 to convert PRBCKGRND from mb to pa for use in this subroutine.
      centralPressureDeficit = PRBCKGRND*100.d0 - cpress
      if (centralPressureDeficit < 100.d0) then
         centralPressureDeficit = 100.d0
      end if

      ! jgf46.29 Subtract the translational speed of the storm from the
      ! observed max wind speed to avoid distortion in the Holland curve
      ! fit. The translational speed will be added back later.
      ts = sqrt(tvx*tvx + tvy*tvy)
      spd = spd - backgroundWindFac*ts

      ! Convert wind speed from 10 meter altitude (which is what the
      ! NHC forecast contains) to wind speed at the top of the atmospheric
      ! boundary layer (which is what the Holland curve fit requires).
      spd = spd/BLAdj

      ! Calculate Holland parameters and limit the result to its appropriate
      ! range.
      if (BCalcId == BCALC_LIMITED) then
         B = rhoAir*e*(spd**2.d0)/(centralPressureDeficit)
         if (B < 1.0d0) B = 1.0d0
         if (B > 2.5d0) B = 2.5d0
      end if

      ! WJP convert RMW to meters
      RMW = RMW*1000d0

      ! Get CLE15 profile for this time steop
      call get_CLE15_profile(spd, RMW, lat, V_cle15, r0)

      ! Calculate wind velocity and pressure at each node.
      do I = 1, NP
         DX(I) = (XCOOR(I) - lon)*DEG2RAD
         DY(I) = (YCOOR(I) - lat)*DEG2RAD
         if (thetaLatDep .eqv. .true.) then
            THETA(I) = atan2(DY(I), (cos(lat*DEG2RAD)*DX(I)))
         else
            THETA(I) = atan2(DY(I), DX(I))
         end if

         ! RJW v48.45
         ! Compute the distances based on haversine formula for distance along a sphere
         rad(i) = sphericalDistance(DX(I), DY(I), lat, ycoor(i))

         ! Add option to calculate exact B parameter
         if (BCalcId == BCALC_EXACT) then
            B = rhoAir*e*(spd**2.d0 + spd*RMW*CORIF(I))/(centralPressureDeficit)
         end if

         ! Compute pressure field
         !P.V 11/04/2022 to account for Southern Hemisphere
         PRESS(I) = cpress + (centralPressureDeficit)*exp(-(RMW/RAD(I))**B)

         !P.V 11/04/2022
         PRESS(I) = PRESS(I)/(RHOWAT0*G)

         ! CLE15 model: complete TC wind profile
         if (RAD(I) > r0) then
            WVNX(I) = 0.0d0
            WVNY(I) = 0.0d0
            cycle
         else
            V_r(I) = V_cle15(int(RAD(I)/1000.d0) + 1)
         end if

         ! Convert wind velocity from top of atmospheric boundary layer to
         ! 10m wind velocity.
         V_r(I) = V_r(I)*BLADj

         ! Find the velocity components.
         rrad = RAD(I)/RMW
         inflowAngle = calcInflowAngle(rrad)

         !P.V 11/04/2022 to account for Southern Hemisphere
         WVNX(I) = sign(1.0d0, -lat)*V_r(I)*sin(THETA(I) + sign(1.0d0, lat)*inflowangle)
         WVNY(I) = sign(1.0d0, lat)*V_r(I)*cos(THETA(I) + sign(1.0d0, lat)*inflowangle)

         ! jgf46.31 Determine translation speed that should be added to final
         ! storm wind speed. This is tapered to zero as the storm wind tapers
         ! to zero toward the eye of the storm and at long distances from the
         ! storm.
         call calcTranslationSpeed(V_r(I), spd, TVX, TVY, lat, rrad, TransSpdX, TransSpdY)

         ! jgf46.31 Add the storm translation speed.
         WVNX(I) = WVNX(I) + TransSpdX
         WVNY(I) = WVNY(I) + TransSpdY

         ! jgf46.21 Also convert from 1 minute averaged winds to 10
         ! minute averaged winds (0.88 factor).
         if (windspeed_averaging_minute == 1) then
            WVNX(I) = WVNX(I)*one2ten
            WVNY(I) = WVNY(I)*one2ten
         end if

      end do
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()

   end subroutine CLE15GET

   !----------------------------------------------------------------
   !> @brief Get the CLE15 wind profile
   !>
   !> @param[in] Vm The maximum wind speed
   !> @param[in] rm The radius of maximum wind
   !> @param[in] lat0 The latitude of the storm center
   !> @param[out] Vout The wind profile
   !> @param[out] r0out The radius of diminishing winds
   !-----------------------------------------------------------------
   subroutine get_CLE15_profile(Vm, rm, lat0, Vout, r0out)
      use adc_constants, only: omega, DEG2RAD
      implicit none
      real(8), intent(IN) :: Vm ! max wind (m/s)
      real(8), intent(IN) :: rm ! radius of max wind (m)
      real(8), intent(IN) :: lat0 ! center of storm (deg)
      real(8), allocatable, intent(OUT) :: Vout(:) ! cle15 wind profile (m/s)
      real(8), intent(OUT) :: r0out ! radius of diminishing winds (m)
      ! local variables
      integer :: I, r_N, num_its, rmerge1, rmerge2, rmerge, rclose
      real(8) :: r0, r0_out_bound, r0_in_bound, r0_mid_bound, tmp, Vdiff
      logical :: NEED_INTERATION
      real(8), dimension(:), allocatable :: r, Vinner, Vouter
      real(8) :: f

      r0 = 1000d3
      r0_out_bound = max(2000d3, 30d0*rm)
      r0_in_bound = 2d0*rm
      r0_mid_bound = 0.5d0*(r0_out_bound + r0_in_bound)
      rclose = huge(1)

      NEED_INTERATION = .true.
      ! calculate latitude dependent coriolis w.r.t center of storm
      f = 2.0d0*omega*sin(abs(lat0)*DEG2RAD)

      num_its = 0
      do while (NEED_INTERATION)

         num_its = num_its + 1

         ! Note that we switched these to subroutine calls.
         ! gfortran < 15 struggles with properly handling the
         ! allocate-on-assign without triggering some warnings
         ! The result is that the compiler ends up inlining
         ! the subroutines instead and the assembly ends up
         ! a bit nicer anyway .. ZC
         call get_inner_wind(Vm, rm, f, r0_mid_bound, Vinner)
         call get_outer_wind(f, r0_mid_bound, Vouter)

         r_N = int(r0_mid_bound/1000d0) + 1
         rmerge1 = -1
         rmerge2 = -1

         do I = 1, r_N, 1

            if (Vinner(I) >= Vouter(I)) then
               if (rmerge1 == -1) then
                  rmerge1 = I
               else
                  rmerge2 = I
               end if
            end if

         end do

         if (rmerge1 == -1 .and. rmerge2 == -1) then
            ! No cross point
            ! ALLOCATE (Vdiff (r_N))
            ! Vdiff = abs(Vinner-Vouter)
            tmp = 9999999999d0

            do I = 1, r_N, 1
               Vdiff = abs(Vinner(I) - Vouter(I))
               if (Vdiff <= tmp) then
                  tmp = Vdiff
                  rclose = I
               end if
            end do

            if (abs(Vinner(rclose) - Vouter(rclose)) < 0.1d0) then
               ! condition to quit while:
               ! 1. wind speed at the closest point < .5
               NEED_INTERATION = .false.
               rmerge = rclose
               r0 = r0_mid_bound
            else
               r0_out_bound = r0_mid_bound
               r0_mid_bound = 0.5d0*(r0_out_bound + r0_in_bound)
               if (r0_out_bound - r0_in_bound <= 1 .and. num_its >= 50) then
                  ! condition to quit while:
                  ! 1. r0_out_bound - r0_in_bound <= 3 r grid point
                  ! 2. loop more than 50 times
                  NEED_INTERATION = .false.
                  rmerge = rclose
                  r0 = r0_mid_bound
               end if
            end if

            ! DEALLOCATE(Vdiff)

         else
            ! At least one cross point

            if (rmerge1 /= -1 .and. rmerge2 == -1) then
               ! condition to quit while:
               ! 1. perfect merge
               NEED_INTERATION = .false.
               rmerge = rmerge1
               r0 = r0_mid_bound

            else

               if (rmerge2 - rmerge1 <= 0) then
                  ! condition to quit while:
                  ! 1. left (rmerge1) and right (rmerge2) cross points
                  ! happen within 3 r grids
                  NEED_INTERATION = .false.
                  rmerge = rmerge1
                  r0 = r0_mid_bound
               else
                  r0_in_bound = r0_mid_bound
                  r0_mid_bound = 0.5d0*(r0_out_bound + r0_in_bound)
                  !condition to quit while:
                  !1. r0_out_bound - r0_in_bound <= 3 r grid point
                  if (r0_out_bound - r0_in_bound <= 1 .and. num_its >= 50) then
                     NEED_INTERATION = .false.
                     rmerge = rmerge1
                     r0 = r0_mid_bound
                  end if
               end if

            end if

         end if

      end do

      r_N = int(r0/1000d0) + 1
      r_N = size(Vinner)
      allocate (r(r_N))
      allocate (Vout(r_N))

      do I = 1, rmerge
         r(I) = dble((I - 1)*1000)
         Vout(I) = Vinner(I)
      end do
      do I = rmerge + 1, r_N
         r(I) = dble((I - 1)*1000)
         Vout(I) = Vouter(I)
      end do

      r0out = r0

      do I = 1, r_N
         if (r(I) < rm) then
            Vout(I) = Vout(I)*(r(I)/rm)**0.25d0
         end if
      end do
   end subroutine get_CLE15_profile

   !----------------------------------------------------------------------
   !> @brief Compute the CkCd value based on the maximum wind speed
   !>
   !> @param Vm The maximum wind speed
   !> @return CkCd_computed The computed CkCd value
   !-----------------------------------------------------------------------
   real(8) pure function compute_ckcd(Vm) result(CkCd_computed)
      implicit none
      real(8), intent(in) :: Vm
      CkCd_computed = min(max(0.00055d0*Vm**2d0 - 0.0259d0*Vm + 0.763d0, 0.4d0), 1.4d0)
   end function compute_ckcd

   !----------------------------------------------------------------------
   !> @brief Get the inner wind profile
   !>
   !> @param Vm The maximum wind speed
   !> @param rm The radius of maximum wind
   !> @param f The Coriolis parameter
   !> @param r0 The radius of diminishing winds
   !> @return res The inner wind profile
   !-----------------------------------------------------------------------
   subroutine get_inner_wind(Vm, rm, f, r0, res)
      implicit none
      real(8), intent(in)   :: Vm, rm, f, r0
      real(8), dimension(:), allocatable, intent(out) :: res
      real(8), dimension(:), allocatable :: r
      ! local
      real(8) :: CkCd_this
      integer :: I, r_N

      r_N = int(r0/1000d0) + 1

      allocate (r(r_N))
      allocate (res(r_N))

      ! if desired calculated CkCd from the best fit relationship
      ! found in
      if (CkCd_calc) then
         CkCd_this = compute_ckcd(Vm)
      else
         CkCd_this = CkCd
      end if

      r = dble([(I, I=0, int(r0), 1000)])
      res = (2d0*(r/rm)**2/(2d0 - CkCd_this + CkCd_this*(r/rm)**2))**(1d0/(2d0 - CkCd_this))* &
            (rm*Vm + 0.5d0*f*rm**2d0)
      res = (res - 0.5d0*f*r**2d0)/r
      res(1) = 0d0

   end subroutine get_inner_wind

   !----------------------------------------------------------------------
   !> @brief Get the outer wind profile
   !>
   !> @param f The Coriolis parameter
   !> @param r0 The radius of diminishing winds
   !> @return res The outer wind profile
   !-----------------------------------------------------------------------
   subroutine get_outer_wind(f, r0, res)
      implicit none
      real(8), intent(in)   :: f, r0
      real(8), dimension(:), allocatable, intent(out) :: res
      real(8), dimension(:), allocatable :: r
      ! local
      real(8) :: Cd
      integer :: I, J, r_N
      real(8) :: dr

      r_N = int(r0/1d3) + 1

      allocate (r(r_N))
      allocate (res(r_N))

      r = dble([(I, I=0, int(r0), 1000)])
      dr = r(2) - r(1)
      res(:) = 0d0

      J = 1
      do I = r_N - 1, 2, -1

         if (res(I) <= 6d0) then
            Cd = 6.16d-4
         else if (res(I) >= 35.4d0) then
            Cd = 2.4d-3
         else
            Cd = 5.19d-5*res(I) + 2.614d-4
         end if

         res(I - 1) = res(I) - (2*Cd/w_cool*r(I)*res(I)**2/(r0**2 - r(I)**2) - res(I)/r(I) - f)*dr

         if (res(I - 1) > 300d0) then
            res(I - 1) = 300d0
            J = I - 1
            exit
         end if

      end do

      if (J > 1) then
         res(:J) = 300d0
      end if

   end subroutine get_outer_wind

   !----------------------------------------------------------------------
   !> @brief Update the storm eye position
   !>
   !> @param lat Latitude of the storm center
   !> @param lon Longitude of the storm center
   !> @param EyeLon Array of longitudes of the storm eye
   !> @param EyeLat Array of latitudes of the storm eye
   !> @param FoundEye Logical flag indicating if the eye has been found
   !>
   !> This subroutine updates the position of the storm eye if the current
   !> latitude and longitude differ from the last recorded position. It
   !> shifts the previous positions in the arrays and updates the current
   !> position.
   !-----------------------------------------------------------------------
   subroutine update_storm_eye_position(lat, lon, EyeLon, EyeLat, FoundEye)
      implicit none
      real(8), parameter :: eps = epsilon(1.d0)
      real(8), intent(in) :: lat, lon
      real(8), intent(inout), dimension(3) :: EyeLon, EyeLat
      logical, intent(inout) :: FoundEye

      if ((abs(lat - EyeLat(3)) > eps) .or. (abs(lon - EyeLon(3)) > eps)) then
         EyeLat(1) = EyeLat(2)
         EyeLon(1) = EyeLon(2)
         EyeLat(2) = EyeLat(3)
         EyeLon(2) = EyeLon(3)
         EyeLat(3) = lat
         EyeLon(3) = lon
         FoundEye = .true.
      end if
   end subroutine update_storm_eye_position

end module mod_nws08
