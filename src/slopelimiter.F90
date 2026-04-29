module SLOPELIMITERS

   use DG, only: ze, slopeflag, el_count, eletab, phi_corner
   use global, only: NOFF, nodecode
   use SIZES, only: MNE, mnei
   use MESH, only: NE, NP, NM

   implicit none

   integer, parameter, private :: sz = 8
   private
   public :: slopelimiter

contains

   subroutine SLOPELIMITER(IRK)

      implicit none
      integer, intent(in) :: IRK

      if (slopeflag == 6) then
         call SlopeLimiter6(IRK)
      end if

   end subroutine SLOPELIMITER

   ! namo - copy of SL5 but got rid of QX,QY
   subroutine SLOPELIMITER6(IRK)
!.....Use appropriate modules

#ifdef CMPI
      use MESSENGER, only: updateR
#endif

      implicit none

      integer, intent(in) :: IRK
!.....Declare local variables

      integer, parameter :: blocksize = 8
      integer ::  block_start, block_end
      integer ::  LL, INC1, INC2, INC3,  I, J, jj, kk, k, bb

      real(SZ), dimension(blocksize, 3) :: ZEC, ZEVERTEX, DIF
      real(SZ), dimension(blocksize) :: SUMLOC, SUMDIF, SIGNDIF
      integer, dimension(blocksize) :: KDP
      real(SZ), dimension(blocksize, 3) :: zemin1, zemax1

      real(sz) :: REDFAC(blocksize), REDMAX(blocksize), bound, div(blocksize)
      real(sz) :: ze_dg(MNEI)
      real(SZ), allocatable :: ZE_MIN1(:), ZE_MAX1(:), QX_MIN1(:), QX_MAX1(:)
      real(SZ), allocatable :: QY_MIN1(:), QY_MAX1(:)
      integer :: N1, N2, N3, NO_NBORS, NBOR_EL

      allocate (ZE_MIN1(NP), ZE_MAX1(NP), QX_MIN1(NP))
      allocate (QY_MIN1(NP), QY_MAX1(NP), QX_MAX1(NP))
!     FIND THE MAXIMUM AND MINIMUM OF EACH VARIABLE OVER ALL ELEMENTS
!     SHARING A NODE

      bound = 1.0d-5

      do I = 1, NP
         ZE_MIN1(I) = 99999.d0
         ZE_MAX1(I) = -99999.d0
         QX_MIN1(I) = 99999.d0
         QX_MAX1(I) = -99999.d0
         QY_MIN1(I) = 99999.d0
         QY_MAX1(I) = -99999.d0

         NO_NBORS = EL_COUNT(I)

         do J = 1, NO_NBORS
            NBOR_EL = ELETAB(I, 1 + J)

            !IF(ncele(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

            ZE_DG(J) = ZE(1, NBOR_EL, IRK + 1)

            if (ZE_DG(J) < ZE_MIN1(I)) then
               ZE_MIN1(I) = ZE_DG(J)
            end if
            if (ZE_DG(J) > ZE_MAX1(I)) then
               ZE_MAX1(I) = ZE_DG(J)
            end if
         end do
      end do
!     LOOP OVER ELEMENTS TO CALCULATE NEW VERTEX VALUES
!

#ifdef CMPI

      call UPDATER(ZE_MIN1, ZE_MAX1, QX_MIN1, 3)
      call UPDATER(QX_MAX1, QY_MIN1, QY_MAX1, 3)
#endif

      bb = 1

      do j = 1, NE, blocksize
         block_start = j
         block_end = min(NE, j + blocksize - 1)

!$omp simd
         simd1: do i = block_start, block_end
            jj = i - block_start + 1
            !IF(ncele(I).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
            N1 = NM(I, 1)
            N2 = NM(I, 2)
            N3 = NM(I, 3)

            ZEC(jj, 1) = ZE(1, I, IRK + 1)
            ZEC(jj, 2) = ZE(2, I, IRK + 1)
            ZEC(jj, 3) = ZE(3, I, IRK + 1)

            ZEMAX1(jj, 1) = ZE_MAX1(N1)
            ZEMIN1(jj, 1) = ZE_MIN1(N1)
            ZEMAX1(jj, 2) = ZE_MAX1(N2)
            ZEMIN1(jj, 2) = ZE_MIN1(N2)
            ZEMAX1(jj, 3) = ZE_MAX1(N3)
            ZEMIN1(jj, 3) = ZE_MIN1(N3)

!     COMPUTE THE VERTEX VALUES

            ZEVERTEX(jj, 1) = ZEC(jj, 1)
            ZEVERTEX(jj, 2) = ZEC(jj, 1)
            ZEVERTEX(jj, 3) = ZEC(jj, 1)
            do KK = 2, 3
               ZEVERTEX(jj, 1) = ZEVERTEX(jj, 1) + PHI_CORNER(KK, 1, 1)*ZEC(jj, KK)
               ZEVERTEX(jj, 2) = ZEVERTEX(jj, 2) + PHI_CORNER(KK, 2, 1)*ZEC(jj, KK)
               ZEVERTEX(jj, 3) = ZEVERTEX(jj, 3) + PHI_CORNER(KK, 3, 1)*ZEC(jj, KK)
            end do

!     RESET THE VERTEX VALUE TO BE LESS THAN OR EQUAL TO THE MAX AND
!     GREATER THAN OR EQUAL TO THE MIN AT THAT VERTEX
!
            ZEVERTEX(jj, 1) = DMAX1(DMIN1(ZEVERTEX(jj, 1), ZEMAX1(jj, 1)), ZEMIN1(jj, 1))
            ZEVERTEX(jj, 2) = DMAX1(DMIN1(ZEVERTEX(jj, 2), ZEMAX1(jj, 2)), ZEMIN1(jj, 2))
            ZEVERTEX(jj, 3) = DMAX1(DMIN1(ZEVERTEX(jj, 3), ZEMAX1(jj, 3)), ZEMIN1(jj, 3))
         end do simd1

!     LOOP OVER THE VERTICES 3 TIMES
!     IF THE VALUE AT THE VERTEX IS ABOVE (BELOW) THE MAX (MIN) AT THAT
!     VERTEX THEN SUBTRACT OFF THE DIFFERENCE AND ADD IT TO THE OTHER
!     VERTICES
!
         do LL = 1, 3
!$omp simd
            simd2: do i = block_start, block_end
               jj = i - block_start + 1
               SUMLOC(jj) = (ZEVERTEX(jj, 1) + ZEVERTEX(jj, 2) + ZEVERTEX(jj, 3))/3.0d0
               SUMDIF(jj) = (SUMLOC(jj) - ZEC(jj, 1))*3.0d0
               SIGNDIF(jj) = DSIGN(1.d0, SUMDIF(jj))
               DIF(jj, 1) = (ZEVERTEX(jj, 1) - ZEC(jj, 1))*SIGNDIF(jj)
               DIF(jj, 2) = (ZEVERTEX(jj, 2) - ZEC(jj, 1))*SIGNDIF(jj)
               DIF(jj, 3) = (ZEVERTEX(jj, 3) - ZEC(jj, 1))*SIGNDIF(jj)
               INC1 = 0
               if (DIF(jj, 1) > 0) INC1 = 1
               INC2 = 0
               if (DIF(jj, 2) > 0) INC2 = 1
               INC3 = 0
               if (DIF(jj, 3) > 0) INC3 = 1
               KDP(jj) = INC1 + INC2 + INC3
            end do simd2
!
!$omp simd
            simd3: do i = block_start, block_end
               jj = i - block_start + 1
               do K = 1, 3
                  DIV(jj) = DMAX1(1.d0, DFLOAT(KDP(jj)))
                  if (DIF(jj, K) > 0) then
                     REDFAC(jj) = SUMDIF(jj)*SIGNDIF(jj)/DIV(jj)
                     KDP(jj) = KDP(jj) - 1
                  else
                     REDFAC(jj) = 0
                  end if
                  if (SIGNDIF(jj) > 0) then
                     REDMAX(jj) = ZEVERTEX(jj, K) - ZEMIN1(jj, K)
                  else
                     REDMAX(jj) = ZEMAX1(jj, K) - ZEVERTEX(jj, K)
                  end if
                  REDFAC(jj) = DMIN1(REDFAC(jj), REDMAX(jj))
                  SUMDIF(jj) = SUMDIF(jj) - REDFAC(jj)*SIGNDIF(jj)
                  ZEVERTEX(jj, K) = ZEVERTEX(jj, K) - REDFAC(jj)*SIGNDIF(jj)
               end do
            end do simd3
         end do
!$omp simd
         simd4: do i = block_start, block_end
            jj = i - block_start + 1
            !IF(ncele(I).EQ.1) then ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

            ZE(2, I, IRK + 1) = -1.d0/6.d0*(ZEVERTEX(jj, 1) + ZEVERTEX(jj, 2)) &
                                + 1.d0/3.d0*ZEVERTEX(jj, 3)
            ZE(3, I, IRK + 1) = -.5d0*ZEVERTEX(jj, 1) + .5d0*ZEVERTEX(jj, 2)
            !endif
         end do simd4

      end do

      return
   end subroutine SLOPELIMITER6
end module slopelimiters
