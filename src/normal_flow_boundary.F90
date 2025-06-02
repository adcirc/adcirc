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
!===============================================================================
!> \file normal_flow_boundary.F90
!> \brief Module for converting and formulating normal flow (i.e. river) boundary conditions
!>        in a rotated coordinate system. Only used when ICS=-20,-22 etc..
!>
!> This module provides routines to
!>   - Convert normal flow (QNIN1, QNIN2) into rotated coordinate flux
!>     components and formulate the rotated flux forcing terms (QN_R, QN2_R)
!>     for the ADCIRC boundary.
!>
!> It uses global variables for boundary indices, rotated normals, and mesh scaling
!> factors to compute the correct projected fluxes along normal flow boundaries.
!===============================================================================

module mod_normal_flow_boundary

   implicit none

   private

   public :: generate_rotated_normal_flow_boundary_data

contains

   !=====================================================================
   !> \brief Convert normal fluxes into rotated coordinate fluxes.
   !>
   !>  For each boundary node J where LBCODEI(J) == 2, 12, 22 or 32, this routine:
   !>    - Reads the two normal flux arrays QN_IN
   !>    - Computes the rotated normal vectors (NX_R, NY_R).
   !>    - Solves for rotated flux components QX_R, QY_R
   !>      such that the flux is projected correctly onto the rotated axes.
   !>
   !> \param[inout]  QN  Array of first normal flow rates at boundary nodes (size NVEL)
   !>                    which will be modified to contain the rotated fluxes.
   !>
   !> \note Requires that Q and Q2 be previously allocated to size MNVEL.
   !>       If they are not allocated, the subroutine returns immediately.
   !=====================================================================
   subroutine generate_rotated_normal_flow_boundary_data(QN)
      use boundaries, only: LBCODEI, NVEL, NBV, CSII, SIII, CSII_OLD, SIII_OLD
      use mesh, only: SFMX, SFCX, SFCY, YCSFAC, RVELF
      use global, only: IFSPROTS

      implicit none

      real(8), intent(inout) :: QN(NVEL)
      real(8) :: QX_R, QY_R, TAUX_R, TAUY_R
      real(8) :: QX_ROT, QY_ROT
      real(8) :: VO(2), VR(2)
      real(8) :: NX_R, NY_R
      integer :: J

      do J = 1, NVEL
         select case (LBCODEI(J))
         case (2, 12, 22, 32, 52)
            ! Compute rotated‐coordinate normals
            NX_R = CSII_OLD(J)/SFMX(NBV(J))
            NY_R = SIII_OLD(J)

            TAUX_R = -1.0d0*NY_R
            TAUY_R = NX_R

            ! Solve for QX1_R and QY1_R
            QY_R = QN(J)/(((-1.0d0*TAUY_R)/TAUX_R)*NX_R + NY_R)
            QX_R = -TAUY_R*QY_R/TAUX_R

            if (IFSPROTS == 1) then
               VO = [QX_R, QY_R]
               VR = matmul(RVELF(1:2, 1:2, J), VO)
               QX_ROT = VR(1); 
               QY_ROT = VR(2); 
            else
               QX_ROT = QX_R
               QY_ROT = QY_R
            end if

            ! Generate the new flux forcing terms
            QN(J) = SFCX(NBV(J))*QX_ROT*CSII(J) + SFCY(NBV(J))*QY_ROT*SIII(J)*YCSFAC(NBV(J))
         end select
      end do

   end subroutine generate_rotated_normal_flow_boundary_data

end module mod_normal_flow_boundary
