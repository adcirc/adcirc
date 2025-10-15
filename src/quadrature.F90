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
! but WITHOUT ANY WARRANTY without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------!
module QUADRATURETRI

   implicit none

   private

   integer :: pcub, ncubs
   real(8), pointer :: cub2d(:, :)
   real(8), allocatable :: brscub2d(:, :)
   real(8), allocatable :: workvec(:)

   type MATARR
      integer :: mdim
      real(8), pointer :: ARRVAL(:, :)
   end type MATARR

   type VECARR
      integer :: mdim
      real(8), pointer :: VECVAL(:)
   end type VECARR

   public :: allocmatarr, allocvecarr, GetDefaultCub2D, FetchDefaultNcub2D, CompElmMsfh, &
             MATARR, VECARR

contains

   !C=================== High-level subroutines ====================
   subroutine AllocMatArr(msarr, iidim)
      implicit none

      integer, intent(in) :: iidim
      type(MATARR), intent(inout) :: msarr

      msarr%mdim = iidim
      if (associated(msarr%arrval)) then
         deallocate (msarr%arrval)
         nullify (msarr%arrval)
      end if

      allocate (msarr%arrval(iidim, iidim))

      return
   end subroutine AllocMatArr

   subroutine AllocVecArr(msvec, iidim)
      implicit none

      integer, intent(in) :: iidim
      type(VECARR), intent(inout) :: msvec

      msvec%mdim = iidim
      if (associated(msvec%vecval)) then
         deallocate (msvec%vecval)
         nullify (msvec%vecval)
      end if

      allocate (msvec%vecval(iidim))

      return
   end subroutine AllocVecArr

   subroutine GetDefaultCub2D(p)
      implicit none

      integer, intent(in) :: p

      if (associated(cub2d)) then
         deallocate (cub2d)
      end if
      if (allocated(brscub2d)) then
         deallocate (brscub2d)
      end if
      if (allocated(workvec)) then
         deallocate (workvec)
      end if

      pcub = p
      call GetQuadrature2DTri(cub2d, pcub)
      ncubs = ubound(cub2d, 1)
      allocate (brscub2d(ncubs, 3))
      call rsbarycentriccoords(brscub2d, &
                               cub2d(:, 1), cub2d(:, 2), ncubs)
      allocate (workvec(ncubs))
      return
   end subroutine GetDefaultCub2D

   subroutine FetchDefaultNcub2d(nn)
      implicit none

      integer, intent(out) :: nn

      nn = ncubs
      return
   end subroutine FetchDefaultNcub2d

   subroutine InterpDefaultCub2D(fxy, fv)
      implicit none

      real(8), dimension(:), intent(in) :: fv
      real(8), dimension(:), intent(out) :: fxy

      call barycentricinterp(fxy, fv, brscub2d, ncubs)
   end subroutine InterpDefaultCub2D

   !c
   !c Compute elemental:
   !c
   !c S = (1/(A/2)) \int f_{h} phi_{i} phi_{j}
   !c
   subroutine CompElmMsfh(S, fv)
      implicit none

      real(8), intent(out) :: S(3, 3)
      real(8), intent(in) :: fv(:)
      integer :: II, JJ

      workvec = 0.d0
      call InterpDefaultCub2D(workvec, fv)
      do II = 1, 3
         do JJ = 1, 3
            S(II, JJ) = sum(workvec*cub2d(:, 3)* &
                            brscub2d(:, II)*brscub2d(:, JJ))
         end do
      end do

      return
   end subroutine CompElmMsfh
   !============================================================

   !
   ! Linear interpolation on triangle
   !
   ! given \lambda | (r,s)  return  f(xy)
   subroutine barycentricinterp(fxy, fv, brs, nrs)
      implicit none

      integer, intent(in) :: nrs
      real(8), dimension(:), intent(out) :: fxy
      real(8), intent(in) :: brs(:, :), fv(:)

      fxy(1:nrs) = fv(1)*brs(1:nrs, 1) + &
                   fv(2)*brs(1:nrs, 2) + fv(3)*brs(1:nrs, 3)
      return
   end subroutine barycentricinterp

   !c RS to lambda
   subroutine rsbarycentriccoords(brs, r, s, nrs)
      implicit none

      integer, intent(in) :: nrs
      real(8), intent(in) :: r(1:nrs), s(1:nrs)
      real(8), intent(inout) :: brs(:, :)

      brs(1:nrs, 1) = -0.5d0*(r(1:nrs) + s(1:nrs))
      brs(1:nrs, 2) = 0.5d0*(r(1:nrs) + 1.d0)
      brs(1:nrs, 3) = 0.5d0*(s(1:nrs) + 1.d0)
      return
   end subroutine rsbarycentriccoords

   !c
   !c ref element:   -1 < x,y < 1, x + y = 0
   subroutine GetQuadrature2DTri(cub2d, pcub)

      implicit none

      integer, intent(in) :: pcub
      real(8), pointer, intent(out) :: cub2d(:, :)

      select case (pcub)
      case (1)
         allocate (cub2d(1, 3))
         cub2d(1, :) = [-3.333333333333330d-01, -3.333333333333330d-01, 2.000000000000000d+00]
      case (2)
         allocate (cub2d(3, 3))
         cub2d(1, :) = [-6.666666666666670d-01, -6.666666666666670d-01, 6.666666666666670d-01]
         cub2d(2, :) = [3.333333333333330d-01, -6.666666666666670d-01, 6.666666666666670d-01]
         cub2d(3, :) = [-6.666666666666670d-01, 3.333333333333330d-01, 6.666666666666670d-01]
      case (3)
         allocate (cub2d(6, 3))
         cub2d(1, :) = [-8.168475729804580d-01, -8.168475729804580d-01, 2.199034873106440d-01]
         cub2d(2, :) = [6.336951459609170d-01, -8.168475729804590d-01, 2.199034873106440d-01]
         cub2d(3, :) = [-8.168475729804590d-01, 6.336951459609170d-01, 2.199034873106440d-01]
         cub2d(4, :) = [-1.081030181680700d-01, -1.081030181680700d-01, 4.467631793560230d-01]
         cub2d(5, :) = [-7.837939636638600d-01, -1.081030181680700d-01, 4.467631793560230d-01]
         cub2d(6, :) = [-1.081030181680700d-01, -7.837939636638600d-01, 4.467631793560230d-01]
      case (4)
         allocate (cub2d(6, 3))
         cub2d(1, :) = [-8.168475729804580d-01, -8.168475729804580d-01, 2.199034873106440d-01]
         cub2d(2, :) = [6.336951459609170d-01, -8.168475729804590d-01, 2.199034873106440d-01]
         cub2d(3, :) = [-8.168475729804590d-01, 6.336951459609170d-01, 2.199034873106440d-01]
         cub2d(4, :) = [-1.081030181680700d-01, -1.081030181680700d-01, 4.467631793560230d-01]
         cub2d(5, :) = [-7.837939636638600d-01, -1.081030181680700d-01, 4.467631793560230d-01]
         cub2d(6, :) = [-1.081030181680700d-01, -7.837939636638600d-01, 4.467631793560230d-01]
      case (5)
         allocate (cub2d(7, 3))
         cub2d(1, :) = [-3.333333333333330d-01, -3.333333333333330d-01, 4.500000000000000d-01]
         cub2d(2, :) = [-5.971587178977000d-02, -5.971587178977000d-02, 2.647883055770120d-01]
         cub2d(3, :) = [-8.805682564204600d-01, -5.971587178977000d-02, 2.647883055770120d-01]
         cub2d(4, :) = [-5.971587178977000d-02, -8.805682564204600d-01, 2.647883055770120d-01]
         cub2d(5, :) = [-7.974269853530870d-01, -7.974269853530870d-01, 2.518783610896540d-01]
         cub2d(6, :) = [5.948539707061750d-01, -7.974269853530870d-01, 2.518783610896540d-01]
         cub2d(7, :) = [-7.974269853530870d-01, 5.948539707061750d-01, 2.518783610896540d-01]
      case (6)
         allocate (cub2d(12, 3))
         cub2d(1, :) = [-5.014265096581790d-01, -5.014265096581790d-01, 2.335725514527590d-01]
         cub2d(2, :) = [2.853019316358000d-03, -5.014265096581790d-01, 2.335725514527590d-01]
         cub2d(3, :) = [-5.014265096581790d-01, 2.853019316358000d-03, 2.335725514527590d-01]
         cub2d(4, :) = [-8.738219710169960d-01, -8.738219710169960d-01, 1.016898127404140d-01]
         cub2d(5, :) = [7.476439420339910d-01, -8.738219710169960d-01, 1.016898127404140d-01]
         cub2d(6, :) = [-8.738219710169960d-01, 7.476439420339910d-01, 1.016898127404140d-01]
         cub2d(7, :) = [-3.792950979324310d-01, -8.937099003103660d-01, 1.657021512367470d-01]
         cub2d(8, :) = [-8.937099003103660d-01, -3.792950979324310d-01, 1.657021512367470d-01]
         cub2d(9, :) = [2.730049982427970d-01, -8.937099003103660d-01, 1.657021512367470d-01]
         cub2d(10, :) = [-8.937099003103660d-01, 2.730049982427970d-01, 1.657021512367470d-01]
         cub2d(11, :) = [2.730049982427970d-01, -3.792950979324310d-01, 1.657021512367470d-01]
         cub2d(12, :) = [-3.792950979324310d-01, 2.730049982427970d-01, 1.657021512367470d-01]
      case DEFAULT
         print *, "Error: Quadrature for pcub = ", pcub, " has not yet been implemented"
      end select

   end subroutine GetQuadrature2DTri

end module QUADRATURETRI
