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
!     M O D U L E  G L O B A L _ W E I R _ V A R I A B L E S
!-----------------------------------------------------------------------
!     THIS MODULE CONTAINS THE VARIABLES USED THROUGHOUT WEIR
!     CALCULATIONS AND THE DIFFERENT TIME LEVELS ASSOCIATED WITH
!     A TIME VARYING WEIR.
!-----------------------------------------------------------------------
module mod_weir_data
   implicit none
   real(8), allocatable :: BARINHT1(:) !...Internal Barrier elevation at previous timestep
   real(8), allocatable :: BARINHT2(:) !...Internal Barrier elevation at current timestep
   real(8), allocatable :: BARLANHT1(:) !...External Barrier elevation at previous timestep
   real(8), allocatable :: BARLANHT2(:) !...External Barrier elevation at current timestep
   real(8), parameter :: BARMIN = 0.04d0
   real(8), parameter :: BARMIN64 = 0.04d0
   real(8), parameter :: SUBMIN64 = 0.04d0

   logical :: EXT_TVW = .false. !...If there are external barrier TVW present, avoids checking an unallocated array
   logical :: INT_TVW = .false. !...If there are internal barrier TVW present, avoids checking an unallocated array
   logical :: FOUND_TVW_NML = .false. !...If the namelist in the fort.15 was found

   !-----------------------------------------------------------------------
contains
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E   A L L O C A T E _ W E I R S
   !-----------------------------------------------------------------------
   !  THIS ROUTINE WILL ALLOCATE THE NECESSARY VARIABLES AND INITIALIZE
   !  THE TIME VARYING ELEVATIONS TO THEIR FORT.14 ELEVATIONS AT TWO TIME
   !  LEVELS
   !-----------------------------------------------------------------------
   subroutine ALLOCATE_WEIRS()
      use mod_logging, only: allMessage, setMessageSource, &
                             unsetMessageSource
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
      use mod_logging, only: DEBUG
#endif
      use BOUNDARIES, only: NVEL, BARINHT, BARLANHT
      implicit none

      call setMessageSource("ALLOCATE_WEIRS")
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif

      allocate (BARINHT1(NVEL))
      allocate (BARINHT2(NVEL))
      allocate (BARLANHT1(NVEL))
      allocate (BARLANHT2(NVEL))

      BARINHT1(1:NVEL) = BARINHT(1:NVEL)
      BARINHT2(1:NVEL) = BARINHT(1:NVEL)
      BARLANHT1(1:NVEL) = BARLANHT(1:NVEL)
      BARLANHT2(1:NVEL) = BARLANHT(1:NVEL)

#if defined(WEIR_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()

   end subroutine ALLOCATE_WEIRS
   !-----------------------------------------------------------------------

   subroutine step_weir_elevation(NFLUXIB, NFLUXB, BARINHT2, BARLANHT2, BARINHT1, BARLANHT1)

      integer, intent(in) :: NFLUXIB
      integer, intent(in) :: NFLUXB
      real(8), intent(in) :: BARINHT2(:)
      real(8), intent(in) :: BARLANHT2(:)
      real(8), intent(inout) :: BARINHT1(:)
      real(8), intent(inout) :: BARLANHT1(:)

      if (NFLUXIB == 1) BARINHT1 = BARINHT2
      if (NFLUXB == 1) BARLANHT1 = BARLANHT2
   end subroutine step_weir_elevation

   !-----------------------------------------------------------------------
end module mod_weir_data
!-----------------------------------------------------------------------