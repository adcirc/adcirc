!===============================================================================
!> \file river.F90
!> \brief Module for converting and formulating river flux boundary conditions
!>        in a rotated coordinate system. Only used when ICS=-20,-22 etc..
!>
!> This module provides routines to
!>   - Convert normal flow (QNIN1, QNIN2) into rotated coordinate flux
!>     components and formulate the rotated flux forcing terms (QN_R, QN2_R)
!>     for the ADCIRC boundary.
!>
!> It uses global variables for boundary indices, rotated normals, and mesh scaling
!> factors to compute the correct projected fluxes along river boundaries.
!===============================================================================

module mod_river_flux_boundary

   implicit none

   private

   public :: generate_rotated_river_boundary_data

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
   subroutine generate_rotated_river_boundary_data(QN)
      use boundaries, only: LBCODEI, NVEL, NBV, CSII, SIII
      use mesh, only: SFMX, SFCX, SFCY, YCSFAC
      implicit none
      real(8), intent(inout) :: QN(NVEL)
      real(8) :: QX_R, QY_R, TAUX_R, TAUY_R
      real(8) :: NX_R, NY_R
      integer :: J

      do J = 1, NVEL
         select case (LBCODEI(J))
         case (2, 12, 22, 32)
            ! Compute rotated‐coordinate normals
            NX_R = CSII(J)/SFMX(NBV(J))
            NY_R = SIII(J)
            TAUX_R = -1.0d0*NY_R
            TAUY_R = NX_R

            ! Solve for QX1_R and QY1_R
            QY_R = QN(J)/(((-1.0d0*TAUY_R)/TAUX_R)*NX_R + NY_R)
            QX_R = -TAUY_R*QY_R/TAUX_R

            ! Generate the new flux forcing terms
            QN(J) = SFCX(NBV(J))*QX_R*CSII(J) + SFCY(NBV(J))*QY_R*SIII(J)*YCSFAC(NBV(J))
         case default
            QN(J) = 0.d0
         end select
      end do

   end subroutine generate_rotated_river_boundary_data

end module mod_river_flux_boundary
