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
!---------------------------------------------------------------------
module INTERNALTIDE
   !---------------------------------------------------------------------
   !     CPB 03/2023: This module was created to make the
   !     apply2dinternalwavedrag subroutine in nodalattr.F more readable.
   !     For now I have simply moved the calculation of the de-tided
   !     velocities here. I plan to move more of the subroutine into this
   !     module to hopefully make it a little clearer what is happening in
   !     that subroutine.
   !---------------------------------------------------------------------
   use GLOBAL, only: DEBUG, ECHO, INFO, WARNING, ERROR, &
                     setMessageSource, unsetMessageSource, allMessage, &
                     scratchMessage, logMessage
   implicit none

   private

   ! stored velocity samples for moving average calculation
   real(8), allocatable :: UAV(:, :), VAV(:, :)
   ! averaged velocities so we don't calculate them every single time
   ! step
   real(8), allocatable :: UBar(:), VBar(:)
   ! weights for filters
   real(8), allocatable :: wts(:)

   public :: UNTIDE, MunkHPFilter, UBar, VBar
   !---------------------------------------------------------------------
contains
   !---------------------------------------------------------------------
   !-----------------------------------------------------------------------
   subroutine UNTIDE(U_i, V_i, TimeStep)
      !---------------------------------------------------------------------
      !     This subroutine takes as inputs the current time step velocity,
      !     updates the average velocity vectors (if necessary) and outputs
      !     the 25 hour lagged average velocity based on a lagged, 25 hour
      !     filter with sampling frequency of 12 minutes. It replaces a good
      !     chunk of code in the apply2dinternalwavedrag to make that
      !     subroutine more readable as well as to eliminate some global
      !     variables
      !---------------------------------------------------------------------
      use MESH, only: NP
      use GLOBAL, only: DTDP
      implicit none
      real(8), intent(IN), dimension(:) :: U_i, V_i ! current velocity
      integer, intent(IN) :: TimeStep ! current time step
      real(8), parameter :: filtL = 25d0*3600d0 ! filter length (s)
      real(8), parameter :: Fs = 12d0*60d0 ! sampling interval (s)
      integer, save :: NS = 1 ! number of samples
      logical, save :: first_call = .true.
      ! need indices to check if the current and previous time step are
      ! in the same 12 min windo2
      integer :: L, Lm
      integer :: ii, kk ! for loops
      ! for populating the UAV, VAV matrices at the start of the run
      integer, save :: ISTA

      if (first_call) then
         first_call = .false.
         NS = floor(filtL/Fs)
         ista = 1
         allocate (UAV(NP, NS), VAV(NP, NS))
         UAV = 0d0
         VAV = 0d0
         allocate (UBar(NP), VBar(NP))
         UBar = 0d0
         VBar = 0d0
      end if

      L = floor(dble(TimeStep)*DTDP/Fs)
      Lm = floor(dble(TimeStep - 1)*DTDP/Fs)
      if (L > Lm) then
         if (ISTA > NS) then
            do ii = 1, NP
               do kk = 1, NS - 1
                  UAV(ii, kk) = UAV(ii, kk + 1)
                  VAV(ii, kk) = VAV(ii, kk + 1)
               end do
               UAV(ii, NS) = U_i(ii)
               VAV(ii, NS) = V_i(ii)
            end do
            do ii = 1, NP
               UBar(ii) = sum(UAV(ii, 1:NS))/dble(NS)
               VBar(ii) = sum(VAV(ii, 1:NS))/dble(NS)
            end do
         else
            do ii = 1, NP
               UAV(ii, ISTA) = U_i(ii)
               VAV(ii, ISTA) = V_i(ii)
               UBar(ii) = sum(UAV(ii, 1:ISTA))/dble(ISTA)
               VBar(ii) = sum(VAV(ii, 1:ISTA))/dble(ISTA)
            end do
            ISTA = ISTA + 1
         end if
      end if
      !-----------------------------------------------------------------------
   end subroutine UNTIDE
   !-----------------------------------------------------------------------
   subroutine MunkHPFilter(U_i, V_i, TimeStep)
      !---------------------------------------------------------------------
      !     This subroutine applies a high-pass filter derived from the
      !     so-called Munk "Tide Killer" filter to the velocity field. The
      !     coefficients for the original, low pass filter can be found at:
      !
      !        https://www.sonel.org/Filters-for-the-daily-mean-sea.html
      !
      !     The high-pass filter is derived from the normalized, nonrecursive
      !     low-pass filter as:
      !
      !        W^{HP}_0 = 1-W^{LP}_0
      !        W^{HP}_k = -W^{LP}_k   (k not equal to 0)
      !
      !     The filter is applied as:
      !
      !        y_n = \sum_{k=-m}^m W_k * x_{n+k}
      !
      !     Note that the output of the filter is lagged by 25 hours from the
      !     current time step (which is fine since we are going to apply this
      !     to the internal wave drag and 25 hours is approximately twice the
      !     semi-diurnal period).
      !
      !     Written by: Coleman Blakely 3/2023
      !---------------------------------------------------------------------
      use MESH, only: NP
      use GLOBAL, only: DTDP
      implicit none
      real(8), intent(IN), dimension(:) :: U_i, V_i ! current velocity
      integer, intent(IN) :: TimeStep ! current time step
      real(8), parameter :: filtL = 49d0*3600d0 ! filter length (s)
      real(8), parameter :: T = 60d0*60d0 ! sampling interval (s)
      integer, save :: NS = 1 ! number of samples
      logical, save :: first_call = .true.
      ! need indices to check if the current and previous time step are
      ! in the same 1 hour window
      integer :: L, Lm
      integer :: ii, kk ! for loops
      ! for populating the UAV, VAV matrices at the start of the run
      integer, save :: ISTA

      if (first_call) then
         first_call = .false.
         NS = floor(filtL/T)
         ista = 1
         allocate (UAV(NP, NS), VAV(NP, NS))
         UAV = 0d0
         VAV = 0d0
         allocate (UBar(NP), VBar(NP))
         UBar = 0d0
         VBar = 0d0
         call CalcMunkWeights()
      end if
      ! check if this and the previous timestep are in the same hour
      L = floor(dble(TimeStep)*DTDP/T)
      Lm = floor(dble(TimeStep - 1)*DTDP/T)
      if (L > Lm) then
         ! check if we are >= 49 hours into the run
         if (ISTA > NS) then
            do ii = 1, NP
               UAV(ii, 1:NS - 1) = UAV(ii, 2:NS)
               VAV(ii, 1:NS - 1) = VAV(ii, 2:NS)
               UAV(ii, NS) = U_i(ii)
               VAV(ii, NS) = V_i(ii)
            end do
            UBar = 0d0
            VBar = 0d0
            do ii = 1, NP
               do kk = 1, NS
                  UBar(ii) = UBar(ii) + wts(kk)*UAV(ii, kk)
                  VBar(ii) = VBar(ii) + wts(kk)*VAV(ii, kk)
               end do
            end do
         else
            ! if we do not have a 49 hr record yet just use a lagged
            ! average
            do ii = 1, NP
               UAV(ii, ISTA) = U_i(ii)
               VAV(ii, ISTA) = V_i(ii)
               UBar(ii) = U_i(ii) - sum(UAV(ii, 1:ISTA))/dble(ISTA)
               VBar(ii) = V_i(ii) - sum(VAV(ii, 1:ISTA))/dble(ISTA)
            end do
            ISTA = ISTA + 1
         end if
      end if
      !---------------------------------------------------------------------
   end subroutine MunkHPFilter
   !---------------------------------------------------------------------
   !---------------------------------------------------------------------
   subroutine CalcMunkWeights()
      !---------------------------------------------------------------------
      !     This subroutine sets up the filter weights for use in the
      !     high-pass filter derived from the so-called Munk "Tide Killer"
      !     low-pass filter.
      !---------------------------------------------------------------------
      implicit none
      real(8), dimension(25) :: LPwts ! original weights (one sided)
      integer, parameter :: NS = 49 ! length of HP filter needed
      real(8) :: K ! for normalizing the LP filter
      integer :: ii ! for loops

      ! allocate weights
      allocate (wts(NS))
      ! define low-pass filter weights (not normalized and one-sided)
      LPwts = [395287d0, 386839d0, 370094d0, 354118d0, 338603d0, 325633d0, 314959d0, &
               300054d0, 278167d0, 251492d0, 234033d0, 219260d0, 208050d0, 195518d0, &
               180727d0, 165525d0, 146225d0, 122665d0, 101603d0, 85349d0, 72261d0, &
               60772d0, 47028d0, 30073d0, 13307d0]
      ! find sum to normalize low pass filter
      K = LPwts(1)
      do ii = 2, 25
         K = K + 2d0*LPwts(ii)
      end do
      ! normalize low pass filter weights
      do ii = 1, 25
         LPwts(ii) = LPwts(ii)/K
      end do
      ! turn one-sided low-pass filter into two-sided high-pass filter
      do ii = 1, 24
         ! -m to -1
         wts(ii) = -LPwts(26 - ii)
         ! 1 to m
         wts(ii + 25) = -LPwts(ii + 1)
      end do
      ! 0 (center point)
      wts(25) = 1 - LPwts(1)
      !---------------------------------------------------------------------
   end subroutine CalcMunkWeights
   !---------------------------------------------------------------------
end module INTERNALTIDE
!---------------------------------------------------------------------
