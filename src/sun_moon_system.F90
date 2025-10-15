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
! Driver subroutine
module mod_moon_sun_coors
   use mod_astronomic, only: t_astronomic_values

   implicit none

   type t_moon_sun
      logical, private :: include_nutation = .true.
      type(t_astronomic_values) :: astronomic_values
   contains
      procedure, public, pass(self) :: set_nutation
      procedure, public, pass(self) :: heavenly_objs_coords_jm
      procedure, public, pass(self) :: gmst_deg_fn
   end type t_moon_sun

   private

   public :: t_moon_sun

contains

   subroutine set_nutation(self, nutation)
      implicit none
      class(t_moon_sun), intent(inout) :: self
      logical, intent(in) :: NUTATION
      self%INCLUDE_NUTATION = NUTATION
   end subroutine set_nutation

   !
   ! Compute the geocentric coodinates of the Moon, Sun
   !
   ! Output:
   !   MOON_POS = Moon coordinates and distance
   !            = (/ RA, DEC, Delta /)
   !   SUN_POS  = Sun coordinates and distance
   !            = (/ RA, DEC, Delta /)
   ! Input:
   !   JD = Julian days
   !
   ! NOTE:
   !   - By default, the nutation is accounted for and tthe output
   !     coodinates are the apparent coordinates.
   !
   !   - To ignore the nutation and get the geometrical
   !     coordinates, call the subroutine
   !
   !         CALL SET_NUTATION( .FALSE. )
   !
   !     prior to the use of this subroutine.
   !
   subroutine HEAVENLY_OBJS_COORDS_JM(self, MOON_POS, SUN_POS, JD, IERR)
      use mod_moon_position, only: MOON_COORDINATES
      use mod_sun_coordinates, only: SUN_COORDINATES
      implicit none
      class(t_moon_sun), intent(inout) :: self
      real(8), intent(IN) :: JD
      real(8), intent(OUT) :: MOON_POS(3), SUN_POS(3)
      integer, intent(OUT) :: IERR

      call self%astronomic_values%compute_astronomic_values(JD)
      call MOON_COORDINATES(MOON_POS(1), MOON_POS(2), MOON_POS(3), JD, self%astronomic_values, self%INCLUDE_NUTATION)
      call SUN_COORDINATES(SUN_POS(1), SUN_POS(2), SUN_POS(3), JD, self%astronomic_values, self%INCLUDE_NUTATION)

      IERR = 0
   end subroutine HEAVENLY_OBJS_COORDS_JM

   real(8) function GMST_DEG_FN(self, JDE) result(GMST)
      use mod_astronomic, only: GMST_DEG
      implicit none
      class(t_moon_sun), intent(IN) :: self
      real(8), intent(IN) :: JDE
      gmst = GMST_DEG(JDE, self%INCLUDE_NUTATION, self%astronomic_values)
   end function GMST_DEG_FN

end module mod_moon_sun_coors
