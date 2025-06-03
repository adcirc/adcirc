module read_river
  use global, only: QN1, QN2, QN0, QX1_R, QY1_R, QX2_R, QY2_R, QN_R, QN2_R, &
                   ScreenUnit, setMessageSource, allMessage, NFFR, NX_R, &
                   NY_R, TAUX_R, TAUY_R, QNIN1, QNIN2
  use boundaries, only: NETA, NFLUXF, NOPE, NVEL, LBCODEI, NPEBC, CSII, &
                       SIII, NVELL, NBV
  use sizes,      only: MNPROC, GLOBALDIR, MNVEL
  use mesh,       only: SFac, SFMX, SFMY, SFCT, SFCX, SFCY, YCSFAC, TANPHI

  implicit none

  integer :: I, J, K
  real(8), allocatable, save :: FORCENODES(:)
  real(8), allocatable        :: Q(:), Q2(:)
  public

contains

  subroutine convert_qn(QNIN1, QNIN2)
    implicit none
    real(8), intent(in) :: QNIN1(NVEL), QNIN2(NVEL)
    integer :: J

    ! Initialize all relevant arrays to zero
    if (.not. allocated(Q))    return
    if (.not. allocated(Q2))   return
    Q   = 0.0d0
    Q2  = 0.0d0
    QX1_R = 0.0d0
    QY1_R = 0.0d0
    QX2_R = 0.0d0
    QY2_R = 0.0d0

    do J = 1, NVEL
      if (LBCODEI(J) == 22 .or. LBCODEI(J) == 32) then
        Q(J)  = QNIN1(J)
        Q2(J) = QNIN2(J)

        ! Compute rotated‐coordinate normals
        NX_R(J)  = CSII(J) / SFMX(NBV(J))
        NY_R(J)  = SIII(J)
        TAUX_R(J) = -1.0d0 * NY_R(J)
        TAUY_R(J) =  NX_R(J)

        ! Solve for QX1_R and QY1_R
        QY1_R(J) = Q(J) / ( ((-1.0d0 * TAUY_R(J)) / TAUX_R(J)) * NX_R(J) + NY_R(J) )
        QX1_R(J) = -TAUY_R(J) * QY1_R(J) / TAUX_R(J)

        ! Solve for QX2_R and QY2_R
        QY2_R(J) = Q2(J) / ( ((-1.0d0 * TAUY_R(J)) / TAUX_R(J)) * NX_R(J) + NY_R(J) )
        QX2_R(J) = -TAUY_R(J) * QY2_R(J) / TAUX_R(J)
      end if
    end do
  end subroutine convert_qn


  subroutine formulate_qforce(QN_R, QN2_R)
    use mesh,       only: SFac, SFMX, SFMY, SFCT, SFCX, SFCY, YCSFAC, TANPHI
    use boundaries, only: NBOU, NVELL, NVEL, NBV, LBCODEI

    implicit none
    real(8), intent(out) :: QN_R(NVEL), QN2_R(NVEL)
    integer :: P

    do P = 1, NVEL
      if (LBCODEI(P) == 22 .or. LBCODEI(P) == 32) then
        QN_R(P) = SFCX(FORCENODES(P)) * QX1_R(P) * CSII(P) + &
                  SFCY(FORCENODES(P)) * QY1_R(P) * SIII(P) * &
                  YCSFAC(FORCENODES(P))

        QN2_R(P) = SFCX(FORCENODES(P)) * QX2_R(P) * CSII(P) + &
                   SFCY(FORCENODES(P)) * QY2_R(P) * SIII(P) * &
                   YCSFAC(FORCENODES(P))
      end if
    end do
  end subroutine formulate_qforce


  subroutine init_river()
    use mesh,       only: SFac, SFMX, SFMY, SFCT, SFCX, SFCY, YCSFAC, TANPHI
    use boundaries, only: NBOU, NVELL, NVEL, NBV

    implicit none

    allocate(Q   (MNVEL))
    allocate(Q2  (MNVEL))
    allocate(FORCENODES(MNVEL))

    do J = 1, NVEL
      if (LBCODEI(J) == 22 .or. LBCODEI(J) == 32) then
        FORCENODES(J) = NBV(J)
      end if
    end do
  end subroutine init_river

end module read_river

