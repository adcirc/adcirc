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
!  MODULE MOMENTUM3D
!-----------------------------------------------------------------------
!> @brief 3D Vertical Sigma (3DVS) momentum equation solver
!>
!> This module contains the 3DVS momentum equation solver which calls
!> VSSOL from vsmy.F. It is kept separate from momentum.F because
!> vsmy.F has dependencies (write_output, globalio, etc.) that are
!> not needed by adcprep.
!-----------------------------------------------------------------------
module momentum3d

   implicit none

   private
   public :: solve3DVSMomentumEq

contains

   !--------------------------------------------------------------------
   ! SUBROUTINE SOLVE3DVSMOMENTUMEQ
   !--------------------------------------------------------------------
   !> @brief Solves the 3DVS momentum equation
   !>
   !> Loads barotropic pressure terms (water level, atmospheric pressure,
   !> tidal potential) averaged between time levels s and s+1, then calls
   !> VSSOL to solve for velocity at the new time level.
   !>
   !> @param[in] IT      Current timestep number
   !> @param[in] TimeLoc Current simulation time
   !--------------------------------------------------------------------
   subroutine solve3DVSMomentumEq(IT, TimeLoc)
      use global, only: MOM_LV_X, ETA1, ETA2, PR1, PR2, TIP1, TIP2, &
                        NWS, CTIP
      use mesh, only: NP
      use ADC_CONSTANTS, only: G
      implicit none

      integer, intent(in) :: IT
      real(8), intent(in) :: TimeLoc
      integer :: I

      interface
         subroutine VSSOL(IT, TimeLoc)
            implicit none
            integer, intent(in) :: IT
            real(8), intent(in) :: TimeLoc
         end subroutine VSSOL
      end interface

      ! Load barotropic pressure terms (water level, atm pressure, tidal potential)
      ! averaged between time levels s and s+1
      do I = 1, NP
         MOM_LV_X(I) = ETA1(I) + ETA2(I)
         if (NWS /= 0) MOM_LV_X(I) = MOM_LV_X(I) + PR1(I) + PR2(I)
         if (CTIP) MOM_LV_X(I) = MOM_LV_X(I) - TIP1(I) - TIP2(I)
         MOM_LV_X(I) = G*MOM_LV_X(I)/2.d0
      end do

      ! Solve for velocity at new time level
      call VSSOL(IT, TimeLoc)

   end subroutine solve3DVSMomentumEq

end module momentum3d
