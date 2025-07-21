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
!
! Fortran module for computing the full luna-solar equilibrium tides of
! different orders of approximation,
!
! Ref:
!   J,-M. Hervouet, Free surface flows: modelling with the finite
!      element method. John Wiley & Sons, Ltd, 2007
!
!   Z. Kowalik and J. L. Luick, Modern theory and practice of tide
!      analysis and tidal power, Austides consulting, 2019
!      https://austides.com/downloads/
!
! DW/AC/ZC, 2024
!
module mod_tidepotential
   use mod_astronomic, only: MassRatioSunEarth, MassRatioMoonEarth, EarthRadiusAU, EarthRadiuskm
   use mod_moon_sun_coors, only: t_moon_sun
   use mod_ephemerides, only: t_ephemerides
   use mod_datetime, only: t_datetime
   implicit none
   !
   !  UseFullTIPFormula = T/F (default F)
   !  TIPOrder = 2 (P2), 3(p3), >= 4(exact)  (default 2)
   !  TIPStartDate = 'BaseDate'/'YYYY-MM-DD HH:mm:SS'
   !  MoonSunPostionMethod = 'JM'/'External' (defaul 'JM')
   !           'JM' = Jean Meeus's approach
   !           'External' = User supplied epherides (Moon's, Sun's RA, DEC and
   !                        distances)
   !  MoonSunCoordFile = filename of an external file storting the
   !                        coordinates of the Moon and the Sun
   !  UniformResMoonSunTimeData = T/F (deafult F). Logical flag
   !                indicating the time intervals of the MoonSunCoordFile.
   !                T = Uniform interval (locating the intervals
   !                    from which the Moon/Sun coordinates are
   !                    interpolate is trivial).
   !                F = non uniform interval (using a binary search
   !                    to locate the intervals from which the
   !                    Moon/Sun coordiantes are interpolated).
   !  IncludeNutation = T/F, Include/exclude the nutation in JM formula
   !  k2value = Love number
   !  h2value = Love number
   !
   real(8), parameter :: massratio(2) = [MassRatioMoonEarth, MassRatioSunEarth]
   real(8), parameter :: significant_radiusearth(2) = [EarthRadiuskm(1), EarthRadiusAu(1)]
   real(8), parameter :: exponent_radiusearth(2) = [EarthRadiuskm(2), EarthRadiusAu(2)]

   integer, parameter :: ComputeMethod_FT = 0
   integer, parameter :: ComputeMethod_JM = 1
   integer, parameter :: ComputeMethod_External = 2

   type t_tidePotential
      private
      logical :: m_UseFullTIPFormula = .false.
      integer :: m_TIPOrder = 2
      type(t_datetime) :: m_TIPStartDate
      integer :: m_MoonSunPositionComputeMethod = ComputeMethod_JM
      logical :: m_UniformResMoonSunTimeData = .false.
      logical :: m_IncludeNutation = .true.
      real(8) :: m_k2value = 0.302d0
      real(8) :: m_h2value = 0.609d0
      real(8) :: m_JDE_BEG
      real(8) :: m_JDE_CURRENT
      real(8), dimension(:), allocatable :: m_AH
      real(8), dimension(:), allocatable :: m_cosZA
      real(8), dimension(:), allocatable :: m_workarr
      real(8), dimension(:), allocatable :: m_workarr1
      real(8), dimension(:), allocatable :: m_CSFEA
      real(8), dimension(:), allocatable :: m_SSFEA
      real(8), dimension(:), allocatable :: m_S2SFEA
      type(t_moon_sun) :: m_moon_sun_position
      type(t_ephemerides) :: m_ephemerides

   contains
      procedure, pass(self), public :: compute => compute_full_tip
      procedure, pass(self), public :: active => useFullTIPFormula
      procedure, pass(self), private :: COMP_FULL_TIP_SUB0
      procedure, pass(self), private :: init => tidalPotentialConstructor
      procedure, pass(self), private :: init_full_tip
   end type t_tidePotential

   interface t_tidePotential
      module procedure createTidalPotential
   end interface t_tidePotential

   type(t_tidePotential) :: tidePotential

   private

   public :: tidePotential, t_tidePotential

contains

   !> @brief Initializes a tidal potential object with the given parameters.
   !>
   !> This subroutine constructs a tidal potential object and initializes its member variables
   !> with the provided input parameters.
   !>
   !> @param self The tidal potential object to be initialized.
   !> @param np The number of points in the computational grid.
   !> @param rnday The number of days in the simulation
   !> @param in_UseFullTIPFormula A logical value indicating whether to use the full TIP formula.
   !> @param in_TIPOrder The order of the TIP formula.
   !> @param in_TIPStartDate The start date for the TIP calculations.
   !> @param in_MoonSunPositionComputeMethod The method used to compute the positions of the Moon and Sun.
   !> @param in_MoonSunCoordFile The file containing the coordinates of the Moon and Sun.
   !> @param in_IncludeNutation A logical value indicating whether to include nutation in the calculations.
   !> @param in_k2value The k2 love number used in the TIP calculations.
   !> @param in_h2value The h2 love number used in the TIP calculations.
   !>
   !> @note This subroutine initializes the tidal potential object by setting its member variables
   !>       with the provided input parameters. It also calls the INIT_FULL_TIP subroutine to initialize
   !>       the full TIP calculations and creates an ephemri object for computing the positions of the
   !>       Moon and Sun from an external file.
   !>
   subroutine tidalPotentialConstructor(self, np, rnday, in_UseFullTIPFormula, in_TIPOrder, in_TIPStartDate, &
                                        in_MoonSunPositionComputeMethod, in_MoonSunCoordFile, &
                                        in_IncludeNutation, in_k2value, in_h2value)
      use global, only: allMessage, WARNING
      implicit none

      class(t_tidePotential), intent(INOUT) :: self
      integer, intent(IN) :: np
      real(8), intent(IN) :: rnday
      logical, intent(IN) :: in_UseFullTIPFormula
      integer, intent(IN) :: in_TIPOrder
      character(LEN=*), intent(IN) :: in_TIPStartDate
      character(LEN=*), intent(IN) :: in_MoonSunPositionComputeMethod
      character(LEN=*), intent(IN) :: in_MoonSunCoordFile
      logical, intent(IN) :: in_IncludeNutation
      real(8), intent(IN) :: in_k2value, in_h2value
      character(len=len_trim(in_MoonSunPositionComputeMethod)) :: moonSunPositionString

      self%m_UseFullTIPFormula = in_UseFullTIPFormula

      if (.not. self%m_UseFullTIPFormula) return

      self%m_TIPOrder = in_TIPOrder
      self%m_TIPStartDate = t_datetime(in_TIPStartDate, &
                                      ["%Y-%m-%d %H:%M:%S    ", &
                                       "%Y-%m-%dT%H:%M:%S    ", &
                                       "%Y-%m-%d %H:%M:%S UTC", &
                                       "%Y-%m-%dT%H:%M:%SZ   ", &
                                       "%Y-%m-%dT%H:%M       ", &
                                       "%Y/%m/%d %H:%M:%S    ", &
                                       "%Y/%m/%d %H:%M       ", &
                                       "%Y%m%d%H%M%S         ", &
                                       "%Y/%m/%d             ", &
                                       "%Y-%m-%d             "] &
                                  )
      self%m_IncludeNutation = in_IncludeNutation
      self%m_k2value = in_k2value
      self%m_h2value = in_h2value

      moonSunPositionString = StringUpper(in_MoonSunPositionComputeMethod)
      if (trim(adjustl(moonSunPositionString)) == 'JM') then
         self%m_MoonSunPositionComputeMethod = ComputeMethod_JM
      elseif (trim(adjustl(MoonSunPositionString)) == 'EXTERNAL') then
         self%m_MoonSunPositionComputeMethod = ComputeMethod_External
      else
         call allMessage(WARNING, "Invalid MoonSunPositionComputeMethod: "//trim(adjustl(in_MoonSunPositionComputeMethod)))
         call allMessage(WARNING, "Defaulting to JM")
         self%m_MoonSunPositionComputeMethod = ComputeMethod_JM
      end if

      call self%INIT_FULL_TIP(np)

      if (self%m_MoonSunPositionComputeMethod == ComputeMethod_External) then
         self%m_ephemerides = t_ephemerides(rnday, in_moonsuncoordfile)
      end if

   end subroutine tidalPotentialConstructor

   !> @brief Creates a tidal potential object.
   !>
   !> @param np Number of processors.
   !> @param in_rnday Number of days in the tidal potential record.
   !> @param in_UseFullTIPFormula Flag indicating whether to use the full Tidal Inverse Problem (TIP) formula.
   !> @param in_TIPOrder Order of the TIP formula.
   !> @param in_TIPStartDate Start date of the TIP formula.
   !> @param in_MoonSunPositionComputeMethod Method for computing the positions of the Moon and Sun.
   !> @param in_MoonSunCoordFile File containing the coordinates of the Moon and Sun.
   !> @param in_IncludeNutation Flag indicating whether to include nutation in the tidal potential calculations.
   !> @param in_k2value Value of the k2 Love number.
   !> @param in_h2value Value of the h2 Love number.
   !>
   !> @return tp Tidal potential object.
   type(t_tidePotential) function createTidalPotential(np, in_rnday, in_UseFullTIPFormula, in_TIPOrder, in_TIPStartDate, &
                                                       in_MoonSunPositionComputeMethod, in_MoonSunCoordFile, &
                                                       in_IncludeNutation, in_k2value, in_h2value) result(tp)
      implicit none

      integer, intent(IN) :: np
      real(8), intent(IN) :: in_rnday
      logical, intent(IN) :: in_UseFullTIPFormula
      integer, intent(IN) :: in_TIPOrder
      character(LEN=*), intent(IN) :: in_TIPStartDate
      character(LEN=*), intent(IN) :: in_MoonSunPositionComputeMethod
      character(LEN=*), intent(IN) :: in_MoonSunCoordFile
      logical, intent(IN) :: in_IncludeNutation
      real(8), intent(IN) :: in_k2value, in_h2value

      call tp%init(np, in_rnday, in_UseFullTIPFormula, in_TIPOrder, in_TIPStartDate, &
                   in_MoonSunPositionComputeMethod, in_MoonSunCoordFile, &
                   in_IncludeNutation, in_k2value, in_h2value)

   end function createTidalPotential

   !> @brief Returns the value of the UseFullTIPFormula flag.
   !> @param self The tidal potential object.
   !> @return The value of the UseFullTIPFormula flag. (Used in self%active)
   logical function useFullTIPFormula(self)
      class(t_tidePotential), intent(IN) :: self
      useFullTIPFormula = self%m_UseFullTIPFormula
   end function useFullTIPFormula

   ! subroutine computing the tidal potential. It is a private
   ! subroutine called by the higher level subroutines
   function COMP_FULL_TIP_SUB0(self, tocgmst, np, slam, MoonSunCoor) result(tipval)
      use ADC_CONSTANTS, only: DEG2RAD
      use mod_astronomic, only: EarthRadiusm
      implicit none

      class(t_tidePotential), intent(INOUT) :: self
      integer, intent(IN) :: np
      real(8), intent(IN) :: tocgmst ! Time in seconds since TIPStartDate !
      real(8), intent(IN) :: MoonSunCoor(3, 2)
      real(8), intent(IN) :: SLAM(:)

      integer :: IOBJ
      real(8) :: RA, DEC, DELTA
      real(8) :: radius_div_Delta, KP, C0

      real(8) :: TIPVAL(np)

      TIPVAL = 0.d0
      do IOBJ = 1, 2
         ! Hour angle !
         self%m_AH = DEG2RAD*tocgmst + SLAM - MoonSunCoor(1, IOBJ)
         ! AH = modulo( tocgmst + SLAM*RAD2DEG, 360.D0 )*DEG2RAD -  MoonSunCoor(1,IOBJ)
         RA = MoonSunCoor(1, IOBJ) ! right ascendsion
         DEC = MoonSunCoor(2, IOBJ) ! declination
         DELTA = MoonSunCoor(3, IOBJ) ! distance (km for the Moon,
         ! AU for the sun

         ! Ratio of the radius of the earth and the Moon/Sun
         radius_div_delta = significant_radiusearth(IOBJ)/DELTA
         radius_div_delta = radius_div_delta*exponent_radiusearth(IOBJ)

         !c cosZA = cos( Z )
         self%m_cosZA = self%m_SSFEA*sin(DEC)
         self%m_cosZA = self%m_cosZA + self%m_CSFEA*cos(DEC)*cos(self%m_AH)

         select case (self%m_TIPOrder)
         case (2, 3)
            ! P2  term. Degree 2 or 3 tidal potential !
            self%m_workarr = 3.d0*self%m_cosZA*self%m_cosZA
            self%m_workarr = self%m_workarr - 1.d0
            self%m_workarr = 0.5d0*self%m_workarr

            KP = AINTPOWER(radius_div_delta, 3)
            KP = KP*EarthRadiusm(1)*EarthRadiusm(2)
            KP = KP*massratio(IOBJ)

            TIPVAL = TIPVAL + KP*self%m_workarr

            ! if include P3 term
            if (self%m_TIPOrder == 3) then
               self%m_workarr = 5.d0*self%m_cosZA*self%m_cosZA*self%m_cosZA
               self%m_workarr = self%m_workarr - 3.d0*self%m_cosZA
               self%m_workarr = 0.5d0*self%m_workarr

               KP = AINTPOWER(radius_div_delta, 4)
               KP = KP*EarthRadiusm(1)*EarthRadiusm(2)
               KP = KP*massratio(IOBJ)
               ! Offset values is zero for this term c!

               TIPVAL = TIPVAL + KP*self%m_workarr
            end if
         case (22)
            !   P2 - Only dirunal + semi-dirunal. A special case of
            !   the degree 2 tidal potential,
            !
            ! Dirunal
            self%m_workarr = (3.d0/4.d0)*self%m_S2SFEA*sin(2.d0*DEC)*cos(self%m_AH)

            ! Semi-Diurnal
            self%m_workarr1 = self%m_CSFEA
            self%m_workarr = self%m_workarr + (3.d0/4.d0)*self%m_workarr1*self%m_workarr1*cos(DEC)*cos(DEC)*cos(2.d0*self%m_AH)

            KP = AINTPOWER(radius_div_delta, 3)
            KP = KP*EarthRadiusm(1)*EarthRadiusm(2)
            KP = KP*massratio(IOBJ)

            TIPVAL = TIPVAL + KP*self%m_workarr
         case DEFAULT
            ! Exact form without any truncation !
            self%m_workarr = 1.d0 + AINTPOWER(radius_div_delta, 2)
            self%m_workarr = self%m_workarr - 2.d0*radius_div_delta*self%m_cosZA

            self%m_workarr = sqrt(self%m_workarr)
            self%m_workarr = 1.d0/self%m_workarr

            self%m_workarr = self%m_workarr - radius_div_delta*self%m_cosZA

            ! Offset value, use Proudman's approach, see Kowalik,
            ! page 14.
            C0 = 1.d0 + AINTPOWER(radius_div_delta, 2)
            C0 = sqrt(C0 + 2.d0*radius_div_delta) - &
                 sqrt(C0 - 2.d0*radius_div_delta)
            C0 = -0.5d0*C0/radius_div_delta

            ! NOTE:
            ! The constant of integration C0 equals to
            !
            ! C0 = -1/2*(1/e)*(sqrt(1 + e**2  + 2*e) - sqrt(1 + e**2 - 2*e))
            ! e =  a/r = radius_div_delta
            !
            ! For radius_div_delta (a/r) << 1, it can be shown
            ! through the use of the Taylor series that
            !
            !    sqrt(1 + (a/r)^2 + 2*(a/r)) \sim  1 + (a/r) + High order terms
            !    sqrt(1 + (a/r)^2 + 2*(a/r)) \sim  1 - (a/r) + High order terms
            !
            ! As a result, a good approximate value of the offset is
            !
            !    C0 \sim -1
            !
            self%m_workarr = self%m_workarr + C0

            KP = AINTPOWER(radius_div_delta, 1)
            KP = KP*EarthRadiusm(1)*EarthRadiusm(2)
            KP = KP*massratio(IOBJ)

            TIPVAL = TIPVAL + KP*self%m_workarr
         end select
      end do

      TIPVAL = (1.d0 + self%m_k2value - self%m_h2value)*TIPVAL

   end function COMP_FULL_TIP_SUB0

   function compute_full_tip(self, TimeLoc, NP, SLAM) result(tip)
      use ADC_CONSTANTS, only: sec2day, DEG2RAD
      use global, only: setMessageSource, unsetMessageSource, allMessage
#ifdef ADCNETCDF
      use mod_ephemerides, only: HEAVENLY_OBJS_COORDS_FROM_TABLE
#else
      use global, only: ERROR
#ifdef CMPI
      use messenger, only: msg_fini
#endif
#endif

#ifdef ALL_TRACE
      use global, only: DEBUG
#endif

      implicit none

      class(t_tidePotential), intent(INOUT) :: self
      real(8), intent(IN) :: TimeLoc ! Time in seconds since TIPStartDate !
      integer, intent(IN) :: NP
      real(8), intent(IN), dimension(:) :: SLAM
      real(8) :: tip(np)

      ! local !
      real(8) :: JDELoc, tocgmst
      real(8) :: MoonSunCoor(3, 2)
      integer :: IERR

      call setMessageSource("comp_full_tip")
#if defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      ! Julian day
      JDELoc = self%m_JDE_BEG + TimeLoc*sec2day

      ! Compute the postions of the Moon and Sun !
      if (self%m_MoonSunPositionComputeMethod == ComputeMethod_JM) then
         ! use algorithms in Jean Meeus's book !
         call self%m_moon_sun_position%HEAVENLY_OBJS_COORDS_JM(MoonSunCoor(:, 1), MoonSunCoor(:, 2), JDELoc, IERR)
      elseif (self%m_MoonSunPositionComputeMethod == ComputeMethod_External) then
         ! interpolate from an external look up table !
#ifdef ADCNETCDF
         call self%m_ephemerides%HEAVENLY_OBJS_COORDS_FROM_TABLE(MoonSunCoor, JDELoc, IERR, self%m_UniformResMoonSunTimeData)
#else
         call allMessage(ERROR, "Must compile with netCDF to use TIP from table")
#ifdef CMPI
         call msg_fini()
#endif
         call exit(1)
#endif
      else
         IERR = 1
      end if
      call check_tip_err(IERR)

      ! Convert RA, DEC from deg --> rad
      MoonSunCoor(1:2, :) = MoonSunCoor(1:2, :)*DEG2RAD

      ! Sidereal time at Greenwich !
      tocgmst = self%m_moon_sun_position%GMST_DEG_FN(JDELoc)

      tip = self%COMP_FULL_TIP_SUB0(tocgmst, np, slam, MoonSunCoor)

#if defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()

   end function compute_full_tip

   subroutine check_tip_err(IERR)
      use global, only: screenMessage, ERROR
      implicit none
      integer, intent(IN) :: IERR

      if (IERR > 1) then
         select case (IERR)
         case (1)
            call screenMessage(ERROR, "Error in the calculation"// &
                               " the postion of the Moon and Sun. ADCIRC is not "// &
                               " built with the -DADCNETCDF flag.")
         case (2)
            call screenMessage(ERROR, "Error in the calculation"// &
                               " the postion of the Moon and Sun. Date not within"// &
                               "  database.")
         end select
      end if

   end subroutine check_tip_err

   function StringUpper(string) result(outString)
      implicit none

      character(LEN=*), intent(IN) :: string
      character(LEN=len_trim(string)) :: outString

      integer :: I, asciinum

      outString = string; 
      do I = 1, len(string)
         asciinum = iachar(string(I:I))
         select case (asciinum)
         case (97:122)
            outString(I:I) = char(asciinum - 32)
         end select
      end do
   end function StringUpper

   ! A^{N}
   elemental function AINTPOWER(A, N) result(AP)
      implicit none

      real(8) :: AP
      real(8), intent(IN) :: A
      integer, intent(IN) :: N

      integer :: I

      AP = 1.d0
      do I = 1, N
         AP = AP*A
      end do

   end function AINTPOWER

   subroutine INIT_FULL_TIP(self, NP)
      use mesh, only: SFEA
      implicit none

      class(t_tidePotential), intent(INOUT) :: self
      integer, intent(IN) :: NP

      ! Get  correspond Julian days !
      self%m_JDE_BEG = self%m_TIPStartDate%julian_day()
      self%m_JDE_CURRENT = self%m_JDE_BEG

      if (self%m_MoonSunPositionComputeMethod == ComputeMethod_JM) then
         call self%m_moon_sun_position%SET_NUTATION(self%m_IncludeNutation)
      end if

      call ALLOCATEWORKARR(self%m_AH, NP)
      call ALLOCATEWORKARR(self%m_cosZA, NP)
      call ALLOCATEWORKARR(self%m_workarr, NP)
      call ALLOCATEWORKARR(self%m_CSFEA, NP)
      call ALLOCATEWORKARR(self%m_SSFEA, NP)

      if (self%m_TIPOrder == 22) then
         call ALLOCATEWORKARR(self%m_S2SFEA, NP)
      end if

      if (self%m_TIPOrder > 3) then
         call ALLOCATEWORKARR(self%m_workarr1, NP)
      end if

      self%m_CSFEA = cos(SFEA)
      self%m_SSFEA = sin(SFEA)
      if (self%m_TIPOrder == 22) then
         self%m_S2SFEA = sin(2.d0*SFEA)
      end if

      return
   contains

      subroutine ALLOCATEWORKARR(arr, N)
         implicit none

         integer, intent(in) :: N
         real(8), allocatable, intent(inout) :: arr(:)

         if (allocated(arr)) then
            deallocate (arr)
         end if
         allocate (arr(1:N))
         arr = 0.d0

         return
      end subroutine ALLOCATEWORKARR

   end subroutine INIT_FULL_TIP

end module mod_tidepotential
