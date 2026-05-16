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
!-----------------------------------------------------------------------
!  MODULE RAMP
!-----------------------------------------------------------------------
!> @brief Hyperbolic-tangent ramp functions for forcing and boundary conditions
!>
!> This module centralizes the temporal ramp applied to ADCIRC forcings
!> (boundary elevation, flux, tidal potential, meteorology, wave radiation
!> stress). The ramp gradually scales a forcing from zero to full strength
!> over a user-specified duration to avoid shocking the model at startup.
!>
!> Two formulations are provided:
!>  - A normalized ramp that divides by tanh(2) so the ramp reaches
!>    exactly 1.0 at the end of the ramp duration. This is the default.
!>  - The original ramp, which is not normalized and therefore
!>    only reaches tanh(2) (~0.964) at the end of the ramp duration,
!>    and will reach nearly 1.0 at ~130% of the ramp duration.
!>
!> The sign of the ramp duration selects the formulation: a positive
!> duration uses the normalized ramp, while a negative duration is a
!> sentinel that reproduces the original, un-normalized behavior.
!-----------------------------------------------------------------------
module mod_ramp
   implicit none

   private

   public :: ramp_function

contains

!-----------------------------------------------------------------------
!> @brief Evaluate the temporal ramp factor for a forcing
!>
!> Returns the ramp multiplier (nominally between 0.0 and 1.0) to apply to
!> a forcing at a given time. The sign of @p duration selects the
!> formulation: a positive duration uses the normalized ramp, while a
!> negative duration is a sentinel that selects the legacy, un-normalized
!> ramp. The sign carries no other meaning, so its magnitude is passed on
!> and both helper functions only ever see a positive duration.
!>
!> @param[in] duration     ramp duration in days; negative selects legacy behavior
!> @param[in] current_time elapsed time in seconds
!> @return    ramp multiplier to apply to the forcing
!-----------------------------------------------------------------------
   pure real(8) function ramp_function(duration, current_time)
      real(8), intent(in) :: duration
      real(8), intent(in) :: current_time

      if (duration < 0d0) then
         ramp_function = ramp_function_original(abs(duration), current_time)
      elseif (duration > 0d0) then
         ramp_function = ramp_function_normalized(duration, current_time)
      else
         ramp_function = 1.0d0
      end if
   end function ramp_function

!-----------------------------------------------------------------------
!> @brief Normalized hyperbolic-tangent ramp
!>
!> Computes a tanh ramp normalized by tanh(2) so that it reaches exactly
!> 1.0 at the end of the ramp duration, capped at 1.0 thereafter.
!>
!> @param[in] duration     ramp duration in days (assumed positive)
!> @param[in] current_time elapsed time in seconds
!> @return    ramp multiplier in the range [.., 1.0]
!-----------------------------------------------------------------------
   pure real(8) function ramp_function_normalized(duration, current_time)
      implicit none
      real(8), intent(in) :: duration
      real(8), intent(in) :: current_time
      real(8), parameter :: sec2day = 1.0d0/86400.0d0
      real(8), parameter :: tan2_inv = 1.0d0/tanh(2.0d0)
      ramp_function_normalized = min(tanh((2.d0*current_time*sec2day)/duration)*tan2_inv, 1.0d0)
   end function ramp_function_normalized

!-----------------------------------------------------------------------
!> @brief Original (legacy) hyperbolic-tangent ramp
!>
!> Computes the historical, un-normalized tanh ramp. Because it is not
!> normalized, it only reaches tanh(2) (~0.964) at the end of the ramp
!> duration rather than 1.0. Retained to preserve legacy behavior.
!>
!> @param[in] duration     ramp duration in days (assumed positive)
!> @param[in] current_time elapsed time in seconds
!> @return    ramp multiplier
!-----------------------------------------------------------------------
   pure real(8) function ramp_function_original(duration, current_time)
      implicit none
      real(8), intent(in) :: duration
      real(8), intent(in) :: current_time
      real(8), parameter :: sec2day = 1.0d0/86400.0d0
      ramp_function_original = tanh((2.d0*current_time*sec2day)/duration)
   end function ramp_function_original

end module mod_ramp
