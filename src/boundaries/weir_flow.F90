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
!-----------------------------------------------------------------------
! WEIR_FLOW.F90
!   Written by Zachary Cobell, 2013/01/04
!              The Water Institute
!              zcobell@thewaterinstitute.org
!              (225)421-2761
!
!   CONTAINS THE NEW ROUTINES FOR THE SPECIFICATION OF TIME VARYING WEIRS AND
!   LAND BOUNDARY CONDITIONS. ALL TYPE 3,13,23,4,24,64,5,25 BOUNDARY CONDITIONS
!   HAVE THEIR ELEVATION VALUES SET USING THIS ROUNTINE. ALSO, WEIR OVERTOPPING
!   IS COMPUTED HERE USING EITHER THE ORIGINAL IMPLEMENTATION (JJW) OR THE
!   STANDARD FORMULATION (SHINTARO BUNYA, CHECKS FOR WET EDGE)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     M O D U L E  W E I R _ F L O W
!-----------------------------------------------------------------------
!   THIS MODULE CONTAINS THE ROUTINES FOR THE WEIR OVERTOPPING FLUX
!   CALCULATIONS FOR INTERNAL AND EXTERNAL WEIR TYPE BOUNDARIES.
!-----------------------------------------------------------------------
module mod_weir_flow
   use ADC_CONSTANTS, only: G
   use GLOBAL, only: ETA2, RAMPINTFLUX
   use mod_logging, only: setMessageSource, unsetMessageSource, DEBUG, allMessage
   use BOUNDARIES, only: LBCODEI, NBV
   use mod_time_varying_weir_boundary, only: COMPUTE_BARRIER_HEIGHT, &
                                             BAR_DEG
   use mod_weir_data, only: BARMIN, BARLANHT1, BARLANHT2, &
                            BARINHT1, BARINHT2, EXT_TVW, INT_TVW, &
                            BARMIN64

   private

   public :: BARINHT1, BARINHT2, BARLANHT1, BARLANHT2
   public :: COMPUTE_EXTERNAL_BOUNDARY_FLUX
   public :: COMPUTE_INTERNAL_BOUNDARY_FLUX
   public :: COMPUTE_INTERNAL_BOUNDARY64_FLUX
   public :: COMPUTE_CROSS_BARRIER_PIPE_FLUX
   public :: SET_SUBMERGED64_AT

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     F U N C T I O N
   !       C O M P U T E _ E X T E R N A L _ B O U N A R Y _ F L U X
   !-----------------------------------------------------------------------
   !> THIS SUBROUTINE COMPUTES THE OVERTOPPING FLUX OUT OF THE DOMAIN FOR
   !> TYPE 3,13,23 BOUNDARY CONDITIONS
   !>
   !> @param BARRIER_INDEX - The index of the barrier
   !> @param TIMELOC - The time location of the simulation
   !> @param FLUX - The flux out of the domain
   !-----------------------------------------------------------------------
   real(8) function COMPUTE_EXTERNAL_BOUNDARY_FLUX(BARRIER_INDEX, TIMELOC, TVW) result(FLUX)
      use BOUNDARIES, only: BARLANCFSP, BARLANHT

      implicit none

      integer, intent(in) :: BARRIER_INDEX
      real(8), intent(in) :: TIMELOC
      real(8), intent(inout) :: TVW(:)

      integer :: NNBB
      real(8) :: RBARWL

      call setMessageSource("COMPUTE_EXTERNAL_BOUNDARY_FLUX")
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif

      FLUX = 0.0d0

      NNBB = NBV(BARRIER_INDEX)
      if (EXT_TVW) then
         if (BAR_DEG(BARRIER_INDEX)) then
            call COMPUTE_BARRIER_HEIGHT(BARRIER_INDEX, &
                                        TIMELOC, BARLANHT1(BARRIER_INDEX), &
                                        BARLANHT2(BARRIER_INDEX))
         end if
         TVW(NNBB) = BARLANHT2(BARRIER_INDEX) - &
                     BARLANHT(BARRIER_INDEX)
      end if

      RBARWL = 2.d0*(ETA2(NNBB) - &
                     BARLANHT2(BARRIER_INDEX))/3.d0
      if (RBARWL > BARMIN) then
         FLUX = -RampIntFlux &
                *BARLANCFSP(BARRIER_INDEX)*RBARWL* &
                (RBARWL*G)**0.5d0
      else
         FLUX = 0.0d0
      end if

#if defined(WEIR_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()

   end function COMPUTE_EXTERNAL_BOUNDARY_FLUX
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     F U N C T I O N
   !       C O M P U T E _ I N T E R N A L _ B O U N D A R Y _ F L U X
   !-----------------------------------------------------------------------
   !> THIS SUBROUTINE CONTAINS THE CALCULATIONS OF OVERTOPPING FLUX
   !> PASSED BETWEEN THE THE WEIR PAIR NODES. THIS IS THE VERSION OF
   !> THE CODE THAT CHECKS IF THE WEIR HAS A WET EDGE BEFORE PASSING
   !> FLUX OVER THE WEIR
   !>
   !> @param BARRIER_INDEX - The index of the barrier
   !> @param BOUNDARYNODE - The boundary node
   !> @param BOUNDARY - The boundary
   !> @param TIMELOC - The time location of the simulation
   !> @return FLUX - The flux in/out of the domain
   !-----------------------------------------------------------------------
   real(8) function COMPUTE_INTERNAL_BOUNDARY_FLUX(BARRIER_INDEX, BOUNDARYNODE, &
                                                   BOUNDARY, TIMELOC, NIBNODECODE, &
                                                   TVW) result(FLUX)

      use BOUNDARIES, only: BARINCFSB, BARINCFSP, NVELL, IBCONN, BARINHT
      use GLOBAL, only: NODECODE

      implicit none

      integer, intent(in), target :: BARRIER_INDEX
      integer, intent(in), target :: BOUNDARY
      integer, intent(in), target :: BOUNDARYNODE
      integer, intent(inout), target :: NIBNODECODE(:)
      real(8), intent(in) :: TIMELOC
      real(8), intent(inout) :: TVW(:)

      integer :: NNBB1, NNBB2
      integer :: NNBB1WN, NNBB2WN
      integer :: FLOWDIR
      integer, pointer :: I, J, K
      real(8) :: RBARWL1, RBARWL2
      real(8) :: RBARWL1F, RBARWL2F

      call setMessageSource("COMPUTE_INTERNAL_BOUNDARY_FLUX")
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif

      !..SIMPLIFY VARIABLES
      I => BARRIER_INDEX
      J => BOUNDARYNODE
      K => BOUNDARY

      NNBB1 = NBV(I) ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
      NNBB2 = IBCONN(I) ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
      NNBB1WN = 0 ! COUNT NUMBER OF WET NEIGHBORS
      NNBB2WN = 0 ! COUNT NUMBER OF WET NEIGHBORS
      FLOWDIR = 0 ! DIRECTION OF FLOW
      FLUX = 0.d0

      !...............CHECK TO SEE IF THERE IS A WET EDGE
      if ((J == 1) .or. (J == NVELL(K) + 1)) then
         NNBB1WN = NNBB1WN + NODECODE(NBV(I + 1))
         NNBB2WN = NNBB2WN + NODECODE(IBCONN(I + 1))
      else
         if (J == NVELL(K) .or. (J == NVELL(K)*2)) then
            NNBB1WN = NNBB1WN + NODECODE(NBV(I - 1))
            NNBB2WN = NNBB2WN + NODECODE(IBCONN(I - 1))
         else
            NNBB1WN = NNBB1WN + NODECODE(NBV(I - 1))
            NNBB1WN = NNBB1WN + NODECODE(NBV(I + 1))
            NNBB2WN = NNBB2WN + NODECODE(IBCONN(I + 1))
            NNBB2WN = NNBB2WN + NODECODE(IBCONN(I - 1))
         end if
      end if

      !..CHECK TO SEE IF THE BARRIER ELEVATION NEEDS UPDATING
      if (INT_TVW) then
         if (BAR_DEG(I)) then
            call COMPUTE_BARRIER_HEIGHT(I, TIMELOC, BARINHT1(I), BARINHT2(I))
         end if
         TVW(NNBB1) = BARINHT2(I) - BARINHT(I)
      end if

      !..GET WATER LEVEL ABOVE WEIR ON EACH SIDE TO COMPUTE HEAD
      RBARWL1 = ETA2(NNBB1) - BARINHT2(I)
      RBARWL2 = ETA2(NNBB2) - BARINHT2(I)

      RBARWL1F = 2.d0*RBARWL1/3.d0
      RBARWL2F = 2.d0*RBARWL2/3.d0
      FLUX = 0.d0

      if ((RBARWL1 < 0.d0) .and. (RBARWL2 < 0.d0)) then
         !..WATER LEVEL ON BOTH SIDES OF BARRIER BELOW BARRIER -> CASE 1
         FLUX = 0.d0
      elseif (abs(RBARWL1 - RBARWL2) < 0.01d0) then
         !..WATER LEVEL EQUAL ON BOTH SIDES OF BARRIER
         !..TO WITHIN TOLERANCE BARMIN -> CASE 2
         FLUX = 0.d0
      elseif ((RBARWL1 > RBARWL2) .and. (RBARWL1 > BARMIN)) then
         !...WATER LEVEL GREATER ON THIS SIDE OF THE BARRIER AND IS SUCH
         !...THAT OVERTOPPING IS OCCURING
         !...THUS THIS SIDE IS THE SOURCE SIDE OF THE FLOW ACROSS THE BARRIER
         !...NOTE THAT WE DO NOT FORCE THE SOURCE SIDE OF THE BARRIER TO
         !...REMAIN WET. ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
         !...BARRIER HAS BEEN DRIED. IF IT HAS WE SHUT DOWN THE FLOW ACROSS
         !...THE BARRIER
         !...ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE BARRIER HAS AT
         !...LEAST ONE WET EDGE. IF NOT, WE SHUT DOWN THE FLOW ACROSS
         !...THE BARRIER. shintaro v46.28.sb05.05 11/01/2006
         if (RBARWL2 > RBARWL1F) then
            !..OUTWARD SUBCRITICAL FLOW -> CASE 3
            if (NODECODE(NNBB1) == 0 .or. NNBB1WN == 0) then
               FLUX = 0.0d0
            else
               FLUX = -RampIntFlux*RBARWL2*BARINCFSB(I)*(2.d0*G*(RBARWL1 - RBARWL2))**0.5d0
            end if
         else
            !..OUTWARD SUPERCRITICAL FLOW -> CASE 4
            if (NODECODE(NNBB1) == 0 .or. NNBB1WN == 0) then
               FLUX = 0.0d0
            else
               FLUX = -RampIntFlux*BARINCFSP(I)*RBARWL1F*(RBARWL1F*G)**0.5d0
            end if
         end if
      elseif ((RBARWL2 > RBARWL1) .and. (RBARWL2 > BARMIN)) then
         !..WATER LEVEL LOWER ON THIS SIDE OF BARRIER AND IS SUCH
         !..THAT OVERTOPPING IS OCCURING
         !..THUS THIS IS THE RECEIVING SIDE OF THE FLOW ACROSS THE BARRIER
         !..NOTE THAT WE DO FORCE THE RECEIVING SIDE OF THE BARRIER TO
         !..REMAIN WET WHEN THERE IS FLOW ACROSS THE BARRIER.
         !..ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
         !..BARRIER HAS BEEN DRIED. IF IT HAS, WE SHUT DOWN THE FLOW ACROSS
         !..THE BARRIER
         !..ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE BARRIER HAS AT
         !..LEAST ONE WET EDGE. IF NOT, WE SHUT DOWN THE FLOW ACROSS
         !..THE BARRIER. shintaro v46.28.sb05.05 11/01/2006
         if (RBARWL1 > RBARWL2F) then
            !..INWARD SUBCRITICAL FLOW -> CASE 5
            if (NODECODE(NNBB2) == 0 .or. NNBB2WN == 0) then
               FLUX = 0.0d0
            else
               FLUX = RampIntFlux*RBARWL1*BARINCFSB(I)*(2.0d0*G*(RBARWL2 - RBARWL1))**0.5d0
               NIBNODECODE(NNBB1) = 1
            end if
         else
            !..INWARD SUPERCRITICAL FLOW -> CASE 6
            if (NODECODE(NNBB2) == 0 .or. NNBB2WN == 0) then
               FLUX = 0.0d0
            else
               FLUX = RampIntFlux*BARINCFSP(I)*RBARWL2F*(RBARWL2F*G)**0.5d0
               NIBNODECODE(NNBB1) = 1
            end if
         end if
      else
         FLUX = 0
      end if

#if defined(WEIR_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()

   end function COMPUTE_INTERNAL_BOUNDARY_FLUX
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E  S E T _ S U B M E R G E D 64 _ A T
   !-----------------------------------------------------------------------
   !  THIS SUBROUTINE SETS THE SUBMERGED FLAG
   !-----------------------------------------------------------------------
   subroutine SET_SUBMERGED64_AT( &
      BARRIER_INDEX, BOUNDARYNODE, BOUNDARY, TIMELOC)

      use BOUNDARIES, only: NVELL, IBCONN, BARINHT, ISSUBMERGED64, NBV
      use GLOBAL, only: TVW, NODECODE
      use MESH, only: LBArray_Pointer, DP
      implicit none
      integer, intent(IN), target   :: BARRIER_INDEX
      integer, intent(IN), target   :: BOUNDARY
      integer, intent(IN), target   :: BOUNDARYNODE
      real(8), intent(IN)         :: TIMELOC

      integer                     :: NNBB1, NNBB2, NNBBL, NNBBH
      integer                     :: NNBB1WN, NNBB2WN
      integer                     :: FLOWDIR
      integer, pointer             :: I, J, K
      real(8)                    :: RBARWL1, RBARWL2, RBARWLL
      real(8)                    :: RBARWL1F, RBARWL2F
      logical                    :: CHECK
      integer                    :: ISSUBMP

      call setMessageSource("SET_SUBMERGED64_AT")
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif

!.....SIMPLIFY VARIABLES
      I => BARRIER_INDEX
      J => BOUNDARYNODE
      K => BOUNDARY

      NNBB1 = NBV(I) ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
      NNBB2 = IBCONN(I) ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
      NNBB1WN = 0 ! COUNT NUMBER OF WET NEIGHBORS
      NNBB2WN = 0 ! COUNT NUMBER OF WET NEIGHBORS
      FLOWDIR = 0 ! DIRECTION OF FLOW

      ! Check if the node ids on both sides of the barrier are valid
      ! Either of them can be invalid ocasionally but not always
      ! when the barrier is splitted for paralell computation
      if ((NNBB1 <= 0) .or. (NNBB2 <= 0)) then
         ISSUBMERGED64(I) = 0
         if (NNBB2 > 0) then
            if (LBArray_Pointer(NNBB2) > 0) then
               ISSUBMERGED64(LBArray_Pointer(NNBB2)) = 0
            end if
         end if
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG, "Return")
#endif
         call unsetMessageSource()
         return
      end if

!...CHECK TO SEE IF THERE IS A WET EDGE
!     This check is needed ocasionally but not always
!     when the barrier is splitted for paralell computation
      CHECK = .false.
      if ((J == 1) .or. (J == NVELL(K) + 1)) then
         if ((NBV(I + 1) > 0) .and. (IBCONN(I + 1) > 0)) then
            NNBB1WN = NNBB1WN + NODECODE(NBV(I + 1))
            NNBB2WN = NNBB2WN + NODECODE(IBCONN(I + 1))
            CHECK = .true.
         end if
      else
         if (J == NVELL(K) .or. (J == NVELL(K)*2)) then
            if ((NBV(I - 1) > 0) .and. (IBCONN(I - 1) > 0)) then
               NNBB1WN = NNBB1WN + NODECODE(NBV(I - 1))
               NNBB2WN = NNBB2WN + NODECODE(IBCONN(I - 1))
               CHECK = .true.
            end if
         else
            if ((NBV(I - 1) > 0) .and. (IBCONN(I - 1) > 0) .and. &
                (NBV(I + 1) > 0) .and. (IBCONN(I + 1) > 0)) then
               NNBB1WN = NNBB1WN + NODECODE(NBV(I - 1))
               NNBB1WN = NNBB1WN + NODECODE(NBV(I + 1))
               NNBB2WN = NNBB2WN + NODECODE(IBCONN(I + 1))
               NNBB2WN = NNBB2WN + NODECODE(IBCONN(I - 1))
               CHECK = .true.
            end if
         end if
      end if

      if (.not. CHECK) then
         ISSUBMERGED64(I) = 0
         if (LBArray_Pointer(NNBB2) > 0) then
            ISSUBMERGED64(LBArray_Pointer(NNBB2)) = 0
         end if
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG, "Return")
#endif
         call unsetMessageSource()
         return
      end if

!.....CHECK TO SEE IF THE BARRIER ELEVATION NEEDS UPDATING
      if (INT_TVW) then
         if (BAR_DEG(I)) then
            call COMPUTE_BARRIER_HEIGHT(I, TIMELOC, BARINHT1(I), &
                                        BARINHT2(I))
         end if
         TVW(NNBB1) = BARINHT2(I) - BARINHT(I)
      end if

!..GET WATER LEVEL ABOVE WEIR ON EACH SIDE TO COMPUTE HEAD
!  DEFINE COMPILER FLAG -DAVERAGEWEIRFLOW TO USE BARAVGWT,
!  WHICH GENERALLY IS SET TO ZERO, AND IS HARD CODED
!  TO ZERO HERE. READ_INPUT.F CONTAINS THE SPECIFICATION OF
!  BARAVGWT. WITH BARAVGWT SET TO ZERO, THERE IS NO NEED FOR
!  IBSTART EITHER.
#if AVERAGEWEIRFLOW
      if (IBSTART == 0) then
         RBARWL1AVG(I) = ETA2(NNBB1) - BARINHT2(I)
         RBARWL2AVG(I) = ETA2(NNBB2) - BARINHT2(I)
         IBSTART = 1
      else
         RBARWL1AVG(I) = (ETA2(NNBB1) - BARINHT2(I) + BARAVGWT &
                          *RBARWL1)/(1.d0 + BARAVGWT)
         RBARWL2AVG(I) = (ETA2(NNBB2) - BARINHT2(I) + BARAVGWT &
                          *RBARWL2)/(1.d0 + BARAVGWT)
      end if
      RBARWL1 = RBARWL1AVG(I)
      RBARWL2 = RBARWL2AVG(I)
#else
!.....STREAMLINE THE PROCESS SINCE BARAVGWT IS ZERO BY DEFAULT
      RBARWL1 = ETA2(NNBB1) - BARINHT2(I)
      RBARWL2 = ETA2(NNBB2) - BARINHT2(I)
#endif
!.....RBARWL on the lower ground side SB
      if (DP(NNBB1) >= DP(NNBB2)) then
         RBARWLL = RBARWL1
         NNBBL = NNBB1
         NNBBH = NNBB2
      else
         RBARWLL = RBARWL2
         NNBBL = NNBB2
         NNBBH = NNBB1
      end if

      RBARWL1F = 2.d0*RBARWL1/3.d0
      RBARWL2F = 2.d0*RBARWL2/3.d0
      ISSUBMP = ISSUBMERGED64(I)
      ISSUBMERGED64(I) = 0
      ISSUBMERGED64(LBArray_Pointer(NNBB2)) = 0

      if (RBARWLL >= BARMIN64) then
         !..THIS IS THE MAJOR DIFFERENCE OF IBTYPE=64 FROM IBTYPE=24.
         !  WATER LEVEL ON BOTH SIDES OF BARRIER ABOVE BARRIER -> CASE 0
         ISSUBMERGED64(I) = 1
         ISSUBMERGED64(LBArray_Pointer(NNBB2)) = 1
      end if

#if defined(WEIR_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()
      return
      !-----------------------------------------------------------------------
   end subroutine SET_SUBMERGED64_AT
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     F U N C T I O N
   !       C O M P U T E _ I N T E R N A L _ B O U N D A R Y 6 4 _ F L U X
   !-----------------------------------------------------------------------
   !> THIS SUBROUTINE CONTAINS THE CALCULATIONS OF OVERTOPPING FLUX
   !> PASSED BETWEEN THE THE WEIR PAIR NODES. THIS IS THE VERSION OF
   !> THE CODE THAT CHECKS IF THE WEIR HAS A WET EDGE BEFORE PASSING
   !> FLUX OVER THE WEIR.
   !> THIS IS A SUBROUTINE THAT WORKS WITH IBTYPE=64. THIS RETURNS ZERO
   !> WHEN A WEIR IS SUBMERGED ON BOTH SIDES. SB
   !>
   !> @param BARRIER_INDEX - The index of the barrier
   !> @param BOUNDARYNODE - The boundary node
   !> @param BOUNDARY - The boundary
   !> @param TIMELOC - The time location of the simulation
   !> @param ISFRONT - A flag to indicate if the boundary is a front
   !> @return FLUX - The flux in/out of the domain
   !-----------------------------------------------------------------------
   real(8) function COMPUTE_INTERNAL_BOUNDARY64_FLUX(BARRIER_INDEX, BOUNDARYNODE, &
                                                     BOUNDARY, TIMELOC, TVW) result(FLUX)

      use BOUNDARIES, only: BARINCFSB, BARINCFSP, &
                            NVELL, IBCONN, BARINHT, ISSUBMERGED64, NBV
      use GLOBAL, only: NODECODE
      use MESH, only: LBArray_Pointer, DP

      implicit none

      integer, intent(in), target :: BARRIER_INDEX
      integer, intent(in), target :: BOUNDARY
      integer, intent(in), target :: BOUNDARYNODE
      real(8), intent(in) :: TIMELOC
      real(8), intent(inout) :: TVW(:)

      integer :: NNBB1, NNBB2
      integer :: NNBB1WN, NNBB2WN
      integer :: FLOWDIR
      integer, pointer :: I, J, K
      real(8) :: RBARWL1, RBARWL2, RBARWLL
      real(8) :: RBARWL1F, RBARWL2F
      logical :: CHECK

      call setMessageSource("COMPUTE_INTERNAL_BOUNDARY64_FLUX")
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif

      !..SIMPLIFY VARIABLES
      I => BARRIER_INDEX
      J => BOUNDARYNODE
      K => BOUNDARY

      NNBB1 = NBV(I) ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
      NNBB2 = IBCONN(I) ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
      NNBB1WN = 0 ! COUNT NUMBER OF WET NEIGHBORS
      NNBB2WN = 0 ! COUNT NUMBER OF WET NEIGHBORS
      FLOWDIR = 0 ! DIRECTION OF FLOW
      FLUX = 0.0d0

      ! Check if the node ids on both sides of the barrier are valid
      ! Either of them can be invalid ocasionally but not always
      ! when the barrier is splitted for paralell computation
      if ((NNBB1 <= 0) .or. (NNBB2 <= 0)) then
         FLUX = 0.d0
         ISSUBMERGED64(I) = 0
         if (NNBB2 > 0) then
            if (LBArray_Pointer(NNBB2) > 0) then
               ISSUBMERGED64(LBArray_Pointer(NNBB2)) = 0
            end if
         end if
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG, "Return")
#endif
         call unsetMessageSource()
         return
      end if

      !..CHECK TO SEE IF THERE IS A WET EDGE
      !    This check is needed ocasionally but not always
      !    when the barrier is splitted for paralell computation
      CHECK = .false.
      if ((J == 1) .or. (J == NVELL(K) + 1)) then
         if ((NBV(I + 1) > 0) .and. (IBCONN(I + 1) > 0)) then
            NNBB1WN = NNBB1WN + NODECODE(NBV(I + 1))
            NNBB2WN = NNBB2WN + NODECODE(IBCONN(I + 1))
            CHECK = .true.
         end if
      else
         if (J == NVELL(K) .or. (J == NVELL(K)*2)) then
            if ((NBV(I - 1) > 0) .and. (IBCONN(I - 1) > 0)) then
               NNBB1WN = NNBB1WN + NODECODE(NBV(I - 1))
               NNBB2WN = NNBB2WN + NODECODE(IBCONN(I - 1))
               CHECK = .true.
            end if
         else
            if ((NBV(I - 1) > 0) .and. (IBCONN(I - 1) > 0) .and. &
                (NBV(I + 1) > 0) .and. (IBCONN(I + 1) > 0)) then
               NNBB1WN = NNBB1WN + NODECODE(NBV(I - 1))
               NNBB1WN = NNBB1WN + NODECODE(NBV(I + 1))
               NNBB2WN = NNBB2WN + NODECODE(IBCONN(I + 1))
               NNBB2WN = NNBB2WN + NODECODE(IBCONN(I - 1))
               CHECK = .true.
            end if
         end if
      end if

      if (.not. CHECK) then
         FLUX = 0.d0
         ISSUBMERGED64(I) = 0
         if (LBArray_Pointer(NNBB2) > 0) then
            ISSUBMERGED64(LBArray_Pointer(NNBB2)) = 0
         end if
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG, "Return")
#endif
         call unsetMessageSource()
         return
      end if

      !..CHECK TO SEE IF THE BARRIER ELEVATION NEEDS UPDATING
      if (INT_TVW) then
         if (BAR_DEG(I)) then
            call COMPUTE_BARRIER_HEIGHT(I, TIMELOC, BARINHT1(I), BARINHT2(I))
         end if
         TVW(NNBB1) = BARINHT2(I) - BARINHT(I)
      end if

      !..GET WATER LEVEL ABOVE WEIR ON EACH SIDE TO COMPUTE HEAD
      RBARWL1 = ETA2(NNBB1) - BARINHT2(I)
      RBARWL2 = ETA2(NNBB2) - BARINHT2(I)

      !..RBARWL on the lower ground side SB
      if (DP(NNBB1) >= DP(NNBB2)) then
         RBARWLL = RBARWL1
      else
         RBARWLL = RBARWL2
      end if

      RBARWL1F = 2.d0*RBARWL1/3.d0
      RBARWL2F = 2.d0*RBARWL2/3.d0
      FLUX = 0.d0
      ISSUBMERGED64(I) = 0
      ISSUBMERGED64(LBArray_Pointer(NNBB2)) = 0

      if (RBARWLL >= BARMIN64) then
         !..THIS IS THE MAJOR DIFFERENCE OF IBTYPE=64 FROM IBTYPE=24.
         !  WATER LEVEL ON BOTH SIDES OF BARRIER ABOVE BARRIER -> CASE 0
         FLUX = 0.d0
         ISSUBMERGED64(I) = 1
         ISSUBMERGED64(LBArray_Pointer(NNBB2)) = 1
      elseif ((RBARWL1 < BARMIN64) .and. (RBARWL2 < BARMIN64)) then
         !..WATER LEVEL ON BOTH SIDES OF BARRIER BELOW BARRIER -> CASE 1
         FLUX = 0.d0
      elseif (abs(RBARWL1 - RBARWL2) < 0.01d0) then
         !..WATER LEVEL EQUAL ON BOTH SIDES OF BARRIER
         !..TO WITHIN TOLERANCE BARMIN -> CASE 2
         FLUX = 0.d0
      elseif ((RBARWL1 > RBARWL2) .and. (RBARWL1 > BARMIN64)) then
         !..WATER LEVEL GREATER ON THIS SIDE OF THE BARRIER AND IS SUCH
         !..THAT OVERTOPPING IS OCCURING
         !..THUS THIS SIDE IS THE SOURCE SIDE OF THE FLOW ACROSS THE BARRIER
         !..NOTE THAT WE DO NOT FORCE THE SOURCE SIDE OF THE BARRIER TO
         !..REMAIN WET. ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
         !..BARRIER HAS BEEN DRIED. IF IT HAS WE SHUT DOWN THE FLOW ACROSS
         !..THE BARRIER
         !..ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE BARRIER HAS AT
         !..LEAST ONE WET EDGE. IF NOT, WE SHUT DOWN THE FLOW ACROSS
         !..THE BARRIER. shintaro v46.28.sb05.05 11/01/2006
         if (RBARWL2 > RBARWL1F) then
            !..OUTWARD SUBCRITICAL FLOW -> CASE 3
            if (NODECODE(NNBB1) == 0 .or. NNBB1WN == 0) then
               FLUX = 0.0d0
            else
               FLUX = -RampIntFlux*RBARWL2*BARINCFSB(I)*(2.d0*G*(RBARWL1 - RBARWL2))**0.5d0
            end if
         else
            !..OUTWARD SUPERCRITICAL FLOW -> CASE 4
            if (NODECODE(NNBB1) == 0 .or. NNBB1WN == 0) then
               FLUX = 0.0d0
            else
               FLUX = -RampIntFlux*BARINCFSP(I)*RBARWL1F*(RBARWL1F*G)**0.5d0
            end if
         end if
      elseif ((RBARWL2 > RBARWL1) .and. (RBARWL2 > BARMIN64)) then
         !...WATER LEVEL LOWER ON THIS SIDE OF BARRIER AND IS SUCH
         !...THAT OVERTOPPING IS OCCURING
         !...THUS THIS IS THE RECEIVING SIDE OF THE FLOW ACROSS THE BARRIER
         !...NOTE THAT WE DO FORCE THE RECEIVING SIDE OF THE BARRIER TO
         !...REMAIN WET WHEN THERE IS FLOW ACROSS THE BARRIER.
         !...ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
         !...BARRIER HAS BEEN DRIED. IF IT HAS, WE SHUT DOWN THE FLOW ACROSS
         !...THE BARRIER
         !...ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE BARRIER HAS AT
         !...LEAST ONE WET EDGE. IF NOT, WE SHUT DOWN THE FLOW ACROSS
         !...THE BARRIER. shintaro v46.28.sb05.05 11/01/2006
         if (RBARWL1 > RBARWL2F) then
            !..INWARD SUBCRITICAL FLOW -> CASE 5
            if (NODECODE(NNBB2) == 0 .or. NNBB2WN == 0) then
               FLUX = 0.0d0
            else
               FLUX = RampIntFlux*RBARWL1*BARINCFSB(I)*(2.0d0*G*(RBARWL2 - RBARWL1))**0.5d0
            end if
         else
            !..INWARD SUPERCRITICAL FLOW -> CASE 6
            if (NODECODE(NNBB2) == 0 .or. NNBB2WN == 0) then
               FLUX = 0.0d0
            else
               FLUX = RampIntFlux*BARINCFSP(I)*RBARWL2F*(RBARWL2F*G)**0.5d0
            end if
         end if
      else
         FLUX = 0.0d0
      end if

#if defined(WEIR_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()

   end function COMPUTE_INTERNAL_BOUNDARY64_FLUX
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     F U N C T I O N
   !       C O M P U T E _ C R O S S _ B A R R I E R _ P I P E _ F L U X
   !-----------------------------------------------------------------------
   !> THIS ROUTINE COMPUTE THE DISCHARGE ACROSS TYPE 5,25 BOUNDARY
   !> CONDITIONS (CROSS BARRIER PIPES)
   !>
   !> @param IDX - The index of the barrier
   !> @param TIMELOC - The time location of the simulation
   !> @param FLUX - The flux out of the domain
   !-----------------------------------------------------------------------
   real(8) function COMPUTE_CROSS_BARRIER_PIPE_FLUX(IDX, NIBNODECODE) result(FLUX)
      use BOUNDARIES, only: NBV, IBCONN, PIPEHT, PIPEDIAM, PIPECOEF
      use ADC_CONSTANTS, only: G, PI
      use GLOBAL, only: NODECODE, ETA2, RampIntFlux

      implicit none

      integer, intent(in) :: IDX
      integer, intent(inout) :: NIBNODECODE(:)

      integer :: NNBB1
      integer :: NNBB2

      real(8) :: RBARWL1
      real(8) :: RBARWL2

      call setMessageSource("COMPUTE_CROSS_BARRIER_PIPE_FLUX")
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif

      NNBB1 = NBV(IDX) ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
      NNBB2 = IBCONN(IDX) ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER

      RBARWL1 = ETA2(NNBB1) - PIPEHT(IDX)
      RBARWL2 = ETA2(NNBB2) - PIPEHT(IDX)

      if ((RBARWL1 < 0.d0) .and. (RBARWL2 < 0.d0)) then
         !..WATER LEVEL ON BOTH SIDES OF BARRIER BELOW PIPE -> CASE 1
         FLUX = 0.d0
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG, "Return")
#endif
         call unsetMessageSource()
         return
      end if
      if (abs(RBARWL1 - RBARWL2) < BARMIN) then
         !..WATER LEVEL EQUAL ON BOTH SIDES OF PIPE -> CASE 2
         FLUX = 0.d0
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG, "Return")
#endif
         call unsetMessageSource()
         return
      end if
      if ((RBARWL1 > RBARWL2) .and. (RBARWL1 > BARMIN)) then
         !..WATER LEVEL GREATER ON THIS SIDE OF THE PIPE AND IS SUCH
         !..THAT OUTWARD DISCHARGE IS OCCURING
         !..THUS THIS SIDE IS THE SOURCE SIDE OF THE FLOW THROUGH THE PIPE
         !..NOTE THAT WE DO NOT FORCE THE SOURCE SIDE OF THE PIPE TO
         !..REMAIN WET. ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
         !..PIPE HAS BEEN DRIED. IF IT HAS, WE SHUT DOWN THE FLOW ACROSS
         !..THE PIPE
         if (RBARWL2 <= 0.d0) then
            !..OUTWARD FREE DISCHARGE -> CASE 3
            if (NODECODE(NNBB1) == 0) then
               FLUX = 0.0d0
            else
               FLUX = -RampIntFlux*0.25d0*PI*(PIPEDIAM(IDX))**2*(2.d0*G*RBARWL1/(1.d0 + PIPECOEF(IDX)))**0.5d0
            end if
         else
            !..OUTWARD SUBMERGED DISCHARGE -> CASE 4
            if (NODECODE(NNBB1) == 0) then
               FLUX = 0.0d0
            else
               FLUX = -RampIntFlux*0.25d0*PI*(PIPEDIAM(IDX))**2*(2.d0*G*(RBARWL1 - RBARWL2)/PIPECOEF(IDX))**0.5d0
            end if
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG, "Return")
#endif
            call unsetMessageSource()
            return
         end if
      end if
      if ((RBARWL2 > RBARWL1) .and. (RBARWL2 > BARMIN)) then
         !..WATER LEVEL LOWER ON THIS SIDE OF PIPE AND IS SUCH
         !..THAT INWARD DISCHARGE IS OCCURING
         !..THUS THIS IS THE RECEIVING SIDE OF THE FLOW THROUGH THE  PIPE
         !..NOTE THAT WE DO FORCE THE RECEIVING SIDE OF THE PIPE TO
         !..REMAIN WET WHEN THERE IS FLOW ACROSS THE PIPE.
         !..ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
         !..PIPE HAS BEEN DRIED. IF IT HAS, WE SHUT DOWN THE FLOW THROUGH
         !..THE PIPE
         if (RBARWL1 <= 0) then
            !..INWARD FREE DISCHARGE -> CASE 5
            if (NODECODE(NNBB2) == 0) then
               FLUX = 0.0d0
            else
               FLUX = RampIntFlux*0.25d0*PI*(PIPEDIAM(IDX))**2*(2.d0*G*RBARWL2/(1.d0 + PIPECOEF(IDX)))**0.5d0
               NIBNODECODE(NNBB1) = 1
            end if
         else
            !..INWARD SUBMERGED DISCHARGE -> CASE 6
            if (NODECODE(NNBB2) == 0) then
               FLUX = 0.0d0
            else
               FLUX = RampIntFlux*0.25d0*PI*(PIPEDIAM(IDX))**2*(2.d0*G*(RBARWL2 - RBARWL1)/PIPECOEF(IDX))**0.5d0
               NIBNODECODE(NNBB1) = 1
            end if
         end if
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG, "Return")
#endif
         call unsetMessageSource()
         return
      end if

#if defined(WEIR_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()

   end function COMPUTE_CROSS_BARRIER_PIPE_FLUX
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module mod_weir_flow
!-----------------------------------------------------------------------

