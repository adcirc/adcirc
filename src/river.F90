!===============================================================================
!> \file river.F90
!> \brief Module for converting and formulating river flux boundary conditions
!>        in a rotated coordinate system. Only used when ICS=-20,-22 etc..
!>
!> This module provides routines to:
!>   - Initialize river forcing nodes.
!>   - Convert normal flow (QNIN1, QNIN2) into rotated coordinate flux components.
!>   - Formulate the rotated flux forcing terms (QN_R, QN2_R) for the ADCIRC boundary.
!>
!> It uses global variables for boundary indices, rotated normals, and mesh scaling
!> factors to compute the correct projected fluxes along river boundaries.
!===============================================================================

module read_river
   use global, only: QN1, QN2, QN0, QX1_R, QY1_R, QX2_R, QY2_R, QN_R, QN2_R, &
                     ScreenUnit, setMessageSource, allMessage, NFFR, NX_R, &
                     NY_R, TAUX_R, TAUY_R, QNIN1, QNIN2
   use boundaries, only: NETA, NFLUXF, NOPE, NVEL, LBCODEI, NPEBC, CSII, &
                         SIII, NVELL, NBV
   use sizes, only: MNPROC, GLOBALDIR, MNVEL
   use mesh, only: SFac, SFMX, SFMY, SFCT, SFCX, SFCY, YCSFAC, TANPHI

   implicit none

   integer :: I, J, K
   integer, allocatable, save :: FORCENODES(:)
   real(8), allocatable        :: Q(:), Q2(:)
   public

contains

!=====================================================================
!> \brief Convert normal fluxes into rotated coordinate fluxes.
!
!  For each boundary node J where LBCODEI(J) == 22 or 32, this routine:
!    - Reads the two normal flux arrays QNIN1(J) and QNIN2(J).
!    - Computes the rotated normal vectors (NX_R(J), NY_R(J)).
!    - Solves for rotated flux components QX1_R, QY1_R, QX2_R, QY2_R
!      such that the flux is projected correctly onto the rotated axes.
!
!> \param[in]  QNIN1  Array of first normal flow rates at boundary nodes (size NVEL).
!> \param[in]  QNIN2  Array of second normal flow rates at boundary nodes (size NVEL).
!>
!> \note Requires that Q and Q2 be previously allocated to size MNVEL.
!>       If they are not allocated, the subroutine returns immediately.
!=====================================================================

   subroutine convert_qn(QNIN1, QNIN2)
      implicit none
      real(8), intent(in) :: QNIN1(NVEL), QNIN2(NVEL)
      integer :: J

      ! Initialize all relevant arrays to zero
      if (.not. allocated(Q)) return
      if (.not. allocated(Q2)) return
      Q = 0.0d0
      Q2 = 0.0d0
      QX1_R = 0.0d0
      QY1_R = 0.0d0
      QX2_R = 0.0d0
      QY2_R = 0.0d0

      do J = 1, NVEL
         if (LBCODEI(J) == 22 .or. LBCODEI(J) == 32) then
            Q(J) = QNIN1(J)
            Q2(J) = QNIN2(J)

            ! Compute rotated‐coordinate normals
            NX_R(J) = CSII(J)/SFMX(NBV(J))
            NY_R(J) = SIII(J)
            TAUX_R(J) = -1.0d0*NY_R(J)
            TAUY_R(J) = NX_R(J)

            ! Solve for QX1_R and QY1_R
            QY1_R(J) = Q(J)/(((-1.0d0*TAUY_R(J))/TAUX_R(J))*NX_R(J) + NY_R(J))
            QX1_R(J) = -TAUY_R(J)*QY1_R(J)/TAUX_R(J)

            ! Solve for QX2_R and QY2_R
            QY2_R(J) = Q2(J)/(((-1.0d0*TAUY_R(J))/TAUX_R(J))*NX_R(J) + NY_R(J))
            QX2_R(J) = -TAUY_R(J)*QY2_R(J)/TAUX_R(J)
         end if
      end do
   end subroutine convert_qn

!=====================================================================
!> \brief Formulate rotated flux forcing terms for boundary conditions.
!
!  Using the rotated flux components QX1_R, QY1_R, QX2_R, QY2_R and mesh
!  scaling factors, this routine computes:
!    - QN_R(P):   the rotated normal forcing term for the first flux.
!    - QN2_R(P):  the rotated normal forcing term for the second flux.
!
!  These arrays are used directly by ADCIRC when applying boundary
!  conditions on river inflows
!
!> \param[out] QN_R   Output array of rotated normal fluxes (size NVEL).
!> \param[out] QN2_R  Output array of second rotated normal fluxes (size NVEL).
!
!> \note Requires that FORCENODES(P) contains the mesh node index for each P.
!=====================================================================

   subroutine formulate_qforce(QN_R, QN2_R)
      use mesh, only: SFCX, SFCY, YCSFAC
      use boundaries, only: NVEL, LBCODEI

      implicit none
      real(8), intent(out) :: QN_R(NVEL), QN2_R(NVEL)
      integer :: P

      do P = 1, NVEL
         if (LBCODEI(P) == 22 .or. LBCODEI(P) == 32) then
            QN_R(P) = SFCX(FORCENODES(P))*QX1_R(P)*CSII(P) + &
                      SFCY(FORCENODES(P))*QY1_R(P)*SIII(P)* &
                      YCSFAC(FORCENODES(P))

            QN2_R(P) = SFCX(FORCENODES(P))*QX2_R(P)*CSII(P) + &
                       SFCY(FORCENODES(P))*QY2_R(P)*SIII(P)* &
                       YCSFAC(FORCENODES(P))
         end if
      end do
   end subroutine formulate_qforce

!=====================================================================
!> \brief Initialize river forcing nodes and allocate arrays.
!
!  This routine allocates the temporary arrays Q, Q2, and FORCENODES
!  (size MNVEL). It then iterates over all boundary indices,
!  storing the node index NBV(J) into FORCENODES(J) for J where
!  LBCODEI(J) == 22 or 32 (i.e., normal flow/radiation boundary).
!
!> \note Must be called once before convert_qn or formulate_qforce.
!=====================================================================

   subroutine init_river()
      use boundaries, only: NVEL, NBV

      implicit none

      allocate (Q(MNVEL))
      allocate (Q2(MNVEL))
      allocate (FORCENODES(MNVEL))

      do J = 1, NVEL
         if (LBCODEI(J) == 22 .or. LBCODEI(J) == 32) then
            FORCENODES(J) = NBV(J)
         end if
      end do
   end subroutine init_river

end module read_river
