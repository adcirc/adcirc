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
#include "logging_macros.h"
!-----------------------------------------------------------------------
! WEIR_BOUNDARY.F
!   Written by Zachary Cobell, 2013/01/04
!              ICEYE
!              zach.cobell@iceye.com
!
!   CONTAINS THE NEW ROUTINES FOR THE SPECIFICATION OF TIME VARYING WEIRS AND
!   LAND BOUNDARY CONDITIONS. ALL TYPE 3,13,23,4,24,64,5,25 BOUNDARY CONDITIONS
!   HAVE THEIR ELEVATION VALUES SET USING THIS ROUNTINE. ALSO, WEIR OVERTOPPING
!   IS COMPUTED HERE USING EITHER THE ORIGINAL IMPLEMENTATION (JJW) OR THE
!   STANDARD FORMULATION (SHINTARO BUNYA, CHECKS FOR WET EDGE)
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!     M O D U L E  G L O B A L _ W E I R _ V A R I A B L E S
!-----------------------------------------------------------------------
!     THIS MODULE CONTAINS THE VARIABLES USED THROUGHOUT WEIR
!     CALCULATIONS AND THE DIFFERENT TIME LEVELS ASSOCIATED WITH
!     A TIME VARYING WEIR.
!-----------------------------------------------------------------------
module WEIR
   implicit none
   private
   real(8), allocatable :: BARINHT1(:) !...Internal Barrier elevation at previous timestep
   real(8), allocatable :: BARINHT2(:) !...Internal Barrier elevation at current timestep
   real(8), allocatable :: BARLANHT1(:) !...External Barrier elevation at previous timestep
   real(8), allocatable :: BARLANHT2(:) !...External Barrier elevation at current timestep
   real(8), parameter :: BARMIN = 0.04d0 !...Minimum water depth above weir to appply weir formula
   real(8) :: BARMIN64 = 0.04d0 !...Minimum water depth above weir to appply weir formula along vertical element walls
   real(8) :: BARMIN64_SUBM = -1.d0 !...Minimum water depth above weir to consider submerged at vertical element walls when the weir is already submerged (negative: 1.2*H0, set in WEIR_SETUP)
   real(8) :: BARMIN64_NOSUBM = -1.d0 !...Minimum water depth above weir to consider submerged at vertical element walls when the weir is not submerged (negative: 1.2*H0, set in WEIR_SETUP)
   real(8) :: BARMIN64_CORNER = -1.d0 !...Wall-top depth below which a submerged VEW pair whose wall-top node has fewer than 3 wet neighbours (or a steep edge slope) is un-submerged (negative: 1.2*H0, set in WEIR_SETUP; was hard-coded 4*H0)
   real(8) :: BARMIN64_SLIM = 0.5d0 !...Minimum water depth above weir to consider submerged at vertical element walls no matter the slope
   logical :: WEIR_ENERGY_HEAD = .false. !...Same for the internal barrier weirs (IBTYPE 4/24): total energy head in the weir formula
   !...WEIR_SEMI_IMPLICIT (lumped GWCE only): internal barrier pairs (IBTYPE 4/24, and IBTYPE 64 pairs in
   !   weir mode) with both nodes wet get their weir flux with both water levels at the new time level
   !   (SEMI_WEIR_SOLVE), instead of from the previous step's levels.
   logical :: WEIR_SEMI_IMPLICIT = .false.
   !...VEW_CREST_HYDRAULICS (VEW1DChannelControl, default off): crest hydraulics at IBTYPE=64 VEWs.
   !   Turns on, together:
   !   - the recovered crest velocity in the advective terms of coupled pairs (vew1d.F90);
   !   - the total energy head (static head + vn**2/(2g), raw nodal velocity normal to the boundary)
   !     in the VEW weir formula;
   !   - the switch between the coupled and the weir mode described below.
   !   The switch uncouples a submerged VEW pair (weir formula) when the crest Froude number
   !   Fr = (q_nb/h_t)/sqrt(g*h_t) exceeds VEW_FROUDE_CRIT; q_nb is the wall-normal discharge averaged
   !   over the foot's neighbours that are not on a VEW, h_t the depth over the crest. Recouple when
   !   Fr < VEW_FROUDE_RECOUPLE*VEW_FROUDE_CRIT (default 0.8) and the head across the wall face is below
   !   0.5*(VEW_FROUDE_RECOUPLE*VEW_FROUDE_CRIT)**2*h_t.
   !   Replaces the depth thresholds BARMIN64_SUBM/NOSUBM (setting both is an error, see WEIR_SETUP).
   !   On walls with condensed tops the test acts on the lee pair only, the pair whose
   !   opposite foot's neighbours stand higher; the upstream pair is left to the depth test and stays
   !   coupled, so the wall top keeps the upstream level and the lee pair is the one weir across the
   !   wall. When a pair leaves coupled mode, its flux history (QN0, QN1) starts from the discharge of
   !   the foot's neighbours, capped at the critical discharge over the crest, and its two levels are
   !   split, conserving volume, to the head at which the weir formula gives that discharge: a coupled
   !   pair has zero flux and one shared level, so otherwise the whole discharge would switch on in
   !   one step from zero head.
   logical :: VEW_CREST_HYDRAULICS = .false.
   real(8) :: VEW_FROUDE_CRIT = 0.5d0 !...Critical crest Froude number for the coupled/weir switch
   real(8) :: VEW_FROUDE_RECOUPLE = 0.8d0 !...Recouple below VEW_FROUDE_RECOUPLE*VEW_FROUDE_CRIT (hysteresis)
   real(8) :: BARSLIM64_ELEM = 99999.d0 !...Water level elemental slope threhold for vertical element walls. A weir is considered not submerged if the water level slope is greater than this value
   real(8) :: BARSLIM64_EDGE = 99999.d0 !...Water level edge slope threhold for vertical element walls. A weir is considered not submerged if the water level slope is greater than this value

!...Averaging information in time for boundaries.
!   This was previously zeroed out in the code, but the
!   option is maintained.
#ifdef AVERAGEWEIRFLOW
   real(8), allocatable     :: RBARWL1AVG(:), RBARWL2AVG(:)
   real(8), allocatable     :: RPIPEWL1AVG(:), RPIPEWL2AVG(:)
   real(8), parameter       :: BARAVGWT = 0.0d0
   integer                 :: IBSTART
#endif

   logical :: EXT_TVW = .false. !...If there are external barrier TVW present, avoids checking an unallocated array
   logical :: INT_TVW = .false. !...If there are internal barrier TVW present, avoids checking an unallocated array
   logical :: FOUND_TVW_NML = .false. !...If the namelist in the fort.15 was found

   public :: BARINHT1, BARINHT2, BARLANHT1, BARLANHT2, ALLOCATE_WEIRS, &
             BARMIN64_SUBM, BARMIN64_NOSUBM, BARMIN64_SLIM, BARSLIM64_ELEM, BARSLIM64_EDGE, &
             BARMIN64_CORNER, BARMIN64, INT_TVW, EXT_TVW, found_tvw_nml, BARMIN, WEIR_ENERGY_HEAD, &
             VEW_CREST_HYDRAULICS, VEW_FROUDE_CRIT, VEW_FROUDE_RECOUPLE, WEIR_SEMI_IMPLICIT

   !-----------------------------------------------------------------------
contains
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E   A L L O C A T E _ W E I R S
   !-----------------------------------------------------------------------
   !  THIS ROUTINE WILL ALLOCATE THE NECESSARY VARIABLES AND INITIALIZE
   !  THE TIME VARYING ELEVATIONS TO THEIR FORT.14 ELEVATIONS AT TWO TIME
   !  LEVELS
   !-----------------------------------------------------------------------
   subroutine ALLOCATE_WEIRS()
      use mod_logging, only: allMessage, t_log_scope, init_log_scope
      use BOUNDARIES, only: NVEL, BARINHT, BARLANHT
      implicit none

      LOG_SCOPE_TRACED("ALLOCATE_WEIRS", WEIR_BOUNDARY_TRACING)

      allocate (BARINHT1(NVEL))
      allocate (BARINHT2(NVEL))
      allocate (BARLANHT1(NVEL))
      allocate (BARLANHT2(NVEL))

#ifdef AVERAGEWEIRFLOW
      allocate (RBARWL1AVG(1:NVEL))
      allocate (RBARWL2AVG(1:NVEL))
      allocate (RPIPEWL1AVG(1:NVEL))
      allocate (RPIPEWL2AVG(1:NVEL))
#endif

      BARINHT1(1:NVEL) = BARINHT(1:NVEL)
      BARINHT2(1:NVEL) = BARINHT(1:NVEL)
      BARLANHT1(1:NVEL) = BARLANHT(1:NVEL)
      BARLANHT2(1:NVEL) = BARLANHT(1:NVEL)

      !-----------------------------------------------------------------------
   end subroutine ALLOCATE_WEIRS
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
end module WEIR
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     M O D U L E  T I M E _ V A R Y I N G _ W E I R _ B O U N D A R Y
!-----------------------------------------------------------------------
!  THIS MODULE CONTAINS ALL THE ROUTINES NECESSARY FOR SPECIFICATION OF
!  TIME VARYING WEIR BOUNDARIES. THESE ROUTINES ARE USED TO SETUP AND
!  CALCULATE THE ELEVATIONS THROUGHOUT THE SIMULATION
!-----------------------------------------------------------------------
module TIME_VARYING_WEIR_BOUNDARY
   use SIZES, only: MYPROC, LOCALDIR, GLOBALDIR
   use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
   use mod_logging, only: logMessage, allMessage, ScreenMessage, &
                          t_log_scope, init_log_scope, ERROR, INFO, ECHO, &
                          DEBUG, WARNING, screenUnit
   use KDTREE2_MODULE, only: KDTREE2, KDTREE2_CREATE, KDTREE2_N_NEAREST, KDTREE2_RESULT, KDTREE2_DESTROY
   use WEIR, only: BARINHT1, BARINHT2, BARLANHT1, BARLANHT2, ALLOCATE_WEIRS, found_tvw_nml, &
                   BARMIN64_SUBM, BARMIN64_NOSUBM, BARMIN64_SLIM, BARSLIM64_ELEM, BARSLIM64_EDGE, &
                   BARMIN64_CORNER, BARMIN64, INT_TVW, EXT_TVW, VEW_CREST_HYDRAULICS

   implicit none

   private

   real(8), allocatable :: BARHT_FINAL(:) !...Final elevation for time varying boundary
   real(8), allocatable :: BAR_LOCATIONS(:, :) !...Location array for boundary conditions
   real(8), allocatable :: GHOST_LOCATIONS(:, :) !...Location array for the ghost nodes
   real(8), allocatable :: BAR_DEG_START(:) !...Start time for time varying boundary
   real(8), allocatable :: BAR_DEG_END(:) !...End time for time varying boundary
   real(8), allocatable :: BAR_ETA_MAX(:) !...Used if the boundary changes at a critical water surface elevation
   real(8), allocatable :: BAR_FAILURE_START(:) !...Used to track the first if ETA_MAX has been exceeded
   real(8), allocatable :: BAR_FAILURE_DURATION(:) !...Amount of time it takes an ETA_MAX barrier to fail to ZF
   logical, allocatable :: BAR_DEG(:) !...T/F this node is a time varying boundary
   integer, allocatable :: BAR_VARYTYPE(:) !...Type of variation to apply to this boundary node
   type(KDTREE2), pointer :: BARRIER_SEARCHTREE !...Search tree for boundary nodes
   type(KDTREE2), pointer :: GHOST_SEARCHTREE !...Search tree for ghost nodes

   !...........CHARACTERIZES A SECTION OF THE SCHEDULE
   type BSCHED
      real(8) :: BAR_DEG_START !...Time to begin this change
      real(8) :: BAR_DEG_END !...Time to end this change
      real(8) :: BARHT_FINAL !...Elevation at the end time for this section
      !   Flag Values:
      !     -99990 = Nodal elevation (minimum value)
      !     -99993 = fort.14 weir elevation
      !     -99991 = Add BARHT_DELTA
      !     -99992 = Subtract BARHT_DELTA
      real(8) :: BARHT_DELTA
   end type BSCHED

   !...........CHARACTERIZES THE SCHEDULE. BASED ON NUMBER OF FILES
   type BARRIER_SCHEDULE_T
      integer :: SID !...Unique ID corresponding to this schedule
      integer :: NSECTIONS !...Number of sections in the schedule
      type(BSCHED), allocatable :: SECTION(:) !...List of all the schedule sections
   end type BARRIER_SCHEDULE_T

   !...........SCHEDULE THAT POINTS TO MAIN SET THAT HOUSES THE ACTUAL DATA
   !           AND IS BASED UPON THE INDIVIDUAL BOUNDARY NODE
   type BARRIER_SCHEDULE_P
      character(1024), pointer :: MYFILE !...File used for this schedule
      integer :: SID !...SID that has a pointer back to SCHEDULE(:)
      integer :: MYSEC !...The current section this node is operating in
      integer :: LOOP !...Should this schedule loop?
      integer :: NLOOPS !...How many times should it loop
      integer :: LOOPSCOMPLETE !...How many loops we've completed
      integer, pointer :: NSECTIONS !...How many sections are in the schedule
      real(8) :: PREVEND !...Previous schedule end
      real(8) :: BARHT_START !...Starting barrier height
      real(8) :: OFFSET !...Offset from starting time in schedule
      type(BSCHED), pointer :: SECTION(:) !...A section of the schedule containing a variation
   end type BARRIER_SCHEDULE_P
   type(BARRIER_SCHEDULE_T), allocatable, target :: SCHEDULE(:)
   type(BARRIER_SCHEDULE_P), allocatable :: BAR_SCHEDULE(:)
   character(1024), allocatable, target :: SCHEDULE_LIST(:)

   !...........THESE VARIABLES ARE PART OF THE NAMELIST THAT IS READ IN
   character(200) :: ScheduleFile
   real(8) :: X1, X2, Y1, Y2
   real(8) :: ZF, ETA_MAX
   real(8) :: TimeStartDay, TimeStartHour, TimeStartMin
   real(8) :: TimeStartSec
   real(8) :: TimeEndDay, TimeEndHour, TimeEndMin, TimeEndSec
   real(8) :: FailureDurationDay, FailureDurationHour
   real(8) :: FailureDurationMin, FailureDurationSec
   real(8) :: HOTADD
   real(8) :: SearchRadius !...Read for VaryType=4 lines, which are handled in mod_tvv
   integer :: VaryType, HOT
   integer :: LOOP, NLOOPS

   integer :: NSCHEDULES

   !...........THE TIME VARYING WEIR NAMELIST
   namelist /TimeVaryingWeir/ &
      X1, Y1, X2, Y2, &
      VaryType, ZF, &
      ETA_MAX, &
      TimeStartDay, TimeStartHour, TimeStartMin, TimeStartSec, &
      TimeEndDay, TimeEndHour, TimeEndMin, TimeEndSec, &
      FailureDurationDay, FailureDurationHour, &
      FailureDurationMin, FailureDurationSec, &
      ScheduleFile, HOT, LOOP, NLOOPS, SearchRadius

   public :: WEIR_SETUP, ALLOCATE_TIMEVARYINGWEIRS, &
             FIND_BOUNDARY_NODES, COMPUTE_BARRIER_HEIGHT, &
             COMPUTE_BARRIER_HEIGHT_LINEAR, COMPUTE_BARRIER_HEIGHT_ETAMAX, &
             COMPUTE_BARRIER_HEIGHT_SCHEDULE, &
             BARHT_FINAL, BAR_DEG_START, BAR_DEG_END, BAR_DEG, &
             BAR_ETA_MAX, BAR_FAILURE_START, BAR_FAILURE_DURATION, &
             BAR_VARYTYPE, BAR_SCHEDULE, SCHEDULE_LIST, NSCHEDULES, &
             BARRIER_SEARCHTREE, GHOST_SEARCHTREE, GHOST_LOCATIONS

   !-----------------------------------------------------------------------
contains
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E   W E I R _ S E T U P
   !-----------------------------------------------------------------------
   !  THIS ROUTINE IS CALLED FROM adcirc.F TO ARRANGE THE WEIRS ONCE
   !  RESNODE HAS BEEN ESTABLISHED
   !-----------------------------------------------------------------------
   subroutine WEIR_SETUP()
      use GLOBAL, only: H0
      use SIZES, only: MNPROC
      use BOUNDARIES, only: NFLUXIB, NFLUXIBP, NFLUXB
      implicit none
      character(1024) :: ScratchMessage

      LOG_SCOPE_TRACED("WEIR_SETUP", WEIR_BOUNDARY_TRACING)

      if ((NFLUXIB == 1) .or. (NFLUXIBP == 1) .or. (NFLUXB == 1)) then

         !                   !...Allocate the new weir arrays at both time levels
         call ALLOCATE_WEIRS()

         !                   !...Parse the time varying weir file on the local processor
         if (found_tvw_nml) then
            call PARSE_TIME_VARYING_WEIR_INFO()
         end if

      end if

      !               !...Set the minimum water depth above weir to appply weir formula along vertical element walls
      BARMIN64 = 1.2d0*H0 ! This is equivalent to HOFF in wetdry.F
      !               !...The submergence thresholds default to the same value unless set in the
      !               !   VEW1DChannelControl namelist (negative means not set)
      !               !   With VEW_CREST_HYDRAULICS the Froude number decides between the coupled and the
      !               !   weir mode, so the depth thresholds are not used (mutually exclusive): a pair
      !               !   can couple as soon as the water on its lower ground side is above the crest.
      if (VEW_CREST_HYDRAULICS) then
         !...serial only for now: the per-wall state (one discharge per wall, weir mode of the
         !   lee pair) is formed from the pairs a subdomain holds, and a subdomain can hold one face
         !   of a wall without the other
         if (MNPROC > 1) call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                                        message="VEW_CREST_HYDRAULICS is serial only for now")
         if (BARMIN64_SUBM >= 0.d0 .or. BARMIN64_NOSUBM >= 0.d0) then
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message="BARMIN64_SUBM/BARMIN64_NOSUBM and VEW_CREST_HYDRAULICS "// &
                           "are mutually exclusive; remove the thresholds from VEW1DChannelControl")
         end if
         BARMIN64_NOSUBM = 0.d0
         BARMIN64_SUBM = 0.d0
      end if
      if (BARMIN64_NOSUBM < 0.d0) BARMIN64_NOSUBM = 1.2d0*H0
      if (BARMIN64_SUBM < 0.d0) BARMIN64_SUBM = 1.2d0*H0
      if (BARMIN64_CORNER < 0.d0) BARMIN64_CORNER = 1.2d0*H0
      write (ScratchMessage, '(a,f15.7,a,f15.7,a,f15.7)') "VEW submergence thresholds: BARMIN64_SUBM=", &
         BARMIN64_SUBM, " BARMIN64_NOSUBM=", BARMIN64_NOSUBM, " BARMIN64_CORNER=", BARMIN64_CORNER
      call allMessage(INFO, ScratchMessage)

      return

      !-----------------------------------------------------------------------
   end subroutine WEIR_SETUP
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E
   !       A L L O C A T E _ T I M E V A R Y I N G W E I R S
   !-----------------------------------------------------------------------
   !  THIS SUBROUTINE SETS UP THE BOUNDARY ARRAYS USED THROUGHOUT THE
   !  OTHER ROUTINES FOR PROCESSING. ALLOCATES ARRAYS USED GLOBALLY FOR
   !  WEIR TYPE BOUNDARY CONDITIONS AND WEIRS WITH PIPES
   !-----------------------------------------------------------------------
   subroutine ALLOCATE_TIMEVARYINGWEIRS()

      use MESH, only: X, Y
      use BOUNDARIES, only: NVEL, LBCODEI, NBV
#ifdef CMPI
      use MESH, only: NP
      use MESSENGER, only: RESNODE
#endif

      implicit none

      integer :: I
      integer :: NWEIR
      integer :: IDX
#ifdef CMPI
      integer :: NGHOST
#endif

      LOG_SCOPE_TRACED("ALLOCATE_TIMEVARYINGWEIRS", WEIR_BOUNDARY_TRACING)

      allocate (BARHT_FINAL(NVEL))
      allocate (BAR_DEG_START(NVEL))
      allocate (BAR_DEG_END(NVEL))
      allocate (BAR_DEG(NVEL))
      allocate (BAR_ETA_MAX(NVEL))
      allocate (BAR_FAILURE_START(NVEL))
      allocate (BAR_FAILURE_DURATION(NVEL))
      allocate (BAR_VARYTYPE(NVEL))
      allocate (BAR_SCHEDULE(NVEL))

      BARHT_FINAL(1:NVEL) = 0d0
      BAR_DEG_START(1:NVEL) = 0d0
      BAR_DEG_END(1:NVEL) = 0d0
      BAR_DEG(1:NVEL) = .false.
      BAR_ETA_MAX(1:NVEL) = 0d0
      BAR_FAILURE_START(1:NVEL) = -1d0
      BAR_FAILURE_DURATION(1:NVEL) = -1d0

      !...............COUNT THE NUMBER OF WEIR STYLE BOUNDARIES
      NWEIR = 0
      do I = 1, NVEL
         select case (LBCODEI(I))
         case (3, 13, 23, 4, 24, 64, 5, 25)
            NWEIR = NWEIR + 1
         end select
      end do

      !...............CONSTRUCT ARRAYS OF X,Y FOR BOUNDARY NODE LOCATIONS
      allocate (BAR_LOCATIONS(3, NWEIR))
      IDX = 0
      do I = 1, NVEL
         select case (LBCODEI(I))
         case (3, 13, 23, 4, 24, 64, 5, 25)
            IDX = IDX + 1
            BAR_LOCATIONS(1, IDX) = X(NBV(I))
            BAR_LOCATIONS(2, IDX) = Y(NBV(I))
            BAR_LOCATIONS(3, IDX) = dble(I)
         end select
      end do

      !...............BUILD SEARCH TREE
      BARRIER_SEARCHTREE => KDTREE2_CREATE( &
                            BAR_LOCATIONS(1:2, :), REARRANGE=.true., SORT=.true.)

#ifdef CMPI
!...............DO THE SAME PROCESS FOR THE GHOST NODES
      NGHOST = 0
      do I = 1, NP
         if (.not. RESNODE(I)) NGHOST = NGHOST + 1
      end do

      allocate (GHOST_LOCATIONS(3, NGHOST))
      IDX = 0
      do I = 1, NP
         if (.not. RESNODE(I)) then
            IDX = IDX + 1
            GHOST_LOCATIONS(1, IDX) = X(I)
            GHOST_LOCATIONS(2, IDX) = Y(I)
            GHOST_LOCATIONS(3, IDX) = dble(I)
         end if
      end do

      GHOST_SEARCHTREE => KDTREE2_CREATE( &
                          GHOST_LOCATIONS(1:2, :), REARRANGE=.true., SORT=.true.)
#endif

      return

      !-----------------------------------------------------------------------
   end subroutine ALLOCATE_TIMEVARYINGWEIRS
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E  F I N D _ B O U N D A R Y _ N O D E S
   !-----------------------------------------------------------------------
   !  THIS SUBROUTINE FINDS THE NODE THAT MATCHS THAT WHICH IS SPECIFIED
   !  AS A TIME VARYING BOUNDARY POINT IN THE INPUT FILE. THE NEAREST
   !  BOUNDARY NODE IS LOCATED. TYPE 4/24/5/25 BOUNDARIES ARE SPECIFIED WITH BOTH
   !  NODES. TYPE 3/13/23 BOUNDARIES ARE SPECIFIED WITH A SINGLE NODE. THIS SCHEME
   !  IS USED SO THAT MESH NUMBERING AND BOUNDARY CONDITION ORDERING IN THE FORT.14
   !  CANNOT INVALIDATE THE TIME VARYING WEIR INPUT FILE. CHECKS ARE IN PLACE TO ENSURE
   !  THE SAME BOUNDARY NODE IS NOT LOCATED TWICE.
   !-----------------------------------------------------------------------
   subroutine FIND_BOUNDARY_NODES(LON, LAT, IDX)
      use MESH, only: ICS, SLAM0, SFEA0, &
                      DRVSPCOORSROTS, CYLINDERMAP
      use GLOBAL, only: DEG2RAD, IFSPROTS

      implicit none
      real(8), intent(IN) :: LAT
      real(8), intent(IN) :: LON
      real(8) :: LAT_TEMP, LON_TEMP
      real(8) :: LATR, LONR
      real(8) :: EPS
      real(8) :: X, Y
      integer, intent(OUT) :: IDX
      integer, parameter :: SEARCHDEPTH = 2
      type(KDTREE2_RESULT) :: KDRESULTS(SEARCHDEPTH)
      character(1024) :: ScratchMessage

      LOG_SCOPE_TRACED("FIND_BOUNDARY_NODES", WEIR_BOUNDARY_TRACING)

      EPS = epsilon(1.0d0)

      if (ICS == 1) then
         X = LON
         Y = LAT
      elseif (ICS /= 1) then
         LAT_TEMP = LAT*DEG2RAD
         LON_TEMP = LON*DEG2RAD
         if (IFSPROTS == 1) then
            call DRVSPCOORSROTS(lonr, latr, LON_TEMP, LAT_TEMP)
         else
            latr = LAT_TEMP; lonr = LON_TEMP
         end if
         call CYLINDERMAP(X, Y, lonr, latr, SLAM0, SFEA0, ICS)
      end if

      call KDTREE2_N_NEAREST(TP=BARRIER_SEARCHTREE, &
                             QV=[X, Y], NN=SEARCHDEPTH, RESULTS=KDRESULTS)

      !...............THIS NEEDS SOME WORK, HOWEVER, IM NOT SURE HOW TO HANDLE THE POSSIBILITY
      !               THAT A NODE HAS TWO BOUNDARY CONDITIONS BESIDES MAKING IT A FATAL ERROR
      if (abs(KDRESULTS(1)%DIS - KDRESULTS(2)%DIS) <= EPS) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message="MULTIPLE LOCATIONS FOUND FOR TIME VARYING "// &
                        "BOUNDARY NODES")
      end if

      if (abs(KDRESULTS(1)%DIS) > EPS) then
#ifdef CMPI
!...................CHECK TO MAKE SURE IT IS NOT INVALID BECAUSE IT
!                   IS A GHOST NODE
         call KDTREE2_N_NEAREST(TP=GHOST_SEARCHTREE, QV=[X, Y], &
                                NN=SEARCHDEPTH, RESULTS=KDRESULTS)
         if (abs(KDRESULTS(1)%DIS) > EPS) then
!.......................THIS IS NOT A GHOST NODE EITHER.

            write (ScratchMessage, '(A,F0.9,A,F0.9,A)') &
               "SPECIFIED NODE LOCATION X=", &
               LAT, " Y=", LON, &
               " NOT FOUND IN LIST OF BOUNDARY NODES."
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message=trim(ScratchMessage))
         else
            IDX = -1
            return
         end if

#else
         write (ScratchMessage, '(A,F0.9,A,F0.9,A)') &
            "SPECIFIED NODE LOCATION X=", &
            LAT, " Y=", LON, &
            " NOT FOUND IN LIST OF BOUNDARY NODES."
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message=trim(ScratchMessage))

#endif
      end if

      IDX = int(BAR_LOCATIONS(3, KDRESULTS(1)%IDX))

      return

      !-----------------------------------------------------------------------
   end subroutine FIND_BOUNDARY_NODES
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E  C O M P U T E _ B A R R I E R _ H E I G H T
   !-----------------------------------------------------------------------
   !  THIS SUBROUTINE WILL DETERMINE IF THE BARRIER HEIGHTS ARE SUPPOSED TO CHANGE AND
   !  WHEN. IT IS CALLED FOR ALL VARIABLE WEIR TYPE BOUNDARIES AND RETURNS THE INITIAL
   !  VALUE FOR THE BOUNDARY IF NO ACTION IS REQUIRED. BOUNDARY IS NOT ALLOWED TO
   !  DECREASE BELOW THE ELEVATION OF THE SURROUNDING TOPOGRAPHY. THIS ESSENTIALLY
   !  FUNCTIONS AS A WRAPPER FOR ALL THE BOUNDARY VARIATION ROUTINES.
   !-----------------------------------------------------------------------
   subroutine COMPUTE_BARRIER_HEIGHT(MYIDX, TIMELOC, &
                                     BAR_HEIGHT_CURRENT, BAR_HEIGHT)

      use BOUNDARIES, only: NBV

      implicit none

      integer, intent(IN) :: MYIDX
      character(20) :: VC
      real(8), intent(IN) :: TIMELOC
      real(8), intent(IN) :: BAR_HEIGHT_CURRENT
      real(8), intent(OUT) :: BAR_HEIGHT
      real(8), parameter :: eps = epsilon(1d0)
      character(1024) :: ScratchMessage

      LOG_SCOPE_TRACED("COMPUTE_BARRIER_HEIGHT", WEIR_BOUNDARY_TRACING)

      !...............SET DEFAULT RETURN VALUE SO WE CAN DUCK OUT AT ANY TIME
      BAR_HEIGHT = BAR_HEIGHT_CURRENT

      !...............SELECT THE APPROPRIATE TYPE OF CHANGE AND CALL THE ASSOCIATED
      !               ROUTINE. FUTURE DEVELOPMENT CAN TAKE PLACE HERE BY ADDING NEW
      !               SUBROUTINES THAT ARE CALLED IN A SIMLIAR WAY. FOR EXAMPLE,
      !               A SUBROUTINE MIGHT BE WRITTEN THAT COMPUTES WAVE FORCES AND
      !               ALTERS THE HEIGHT OF THE BARRIER BASED UPON THAT.
      select case (BAR_VARYTYPE(MYIDX))
      case (1)
         call COMPUTE_BARRIER_HEIGHT_LINEAR(MYIDX, &
                                            TIMELOC, BAR_HEIGHT_CURRENT, BAR_HEIGHT)
      case (2)
         call COMPUTE_BARRIER_HEIGHT_ETAMAX(MYIDX, &
                                            TIMELOC, BAR_HEIGHT_CURRENT, BAR_HEIGHT)
      case (3)
         call COMPUTE_BARRIER_HEIGHT_SCHEDULE(MYIDX, &
                                              TIMELOC, BAR_HEIGHT_CURRENT, BAR_HEIGHT)
      case DEFAULT
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message="Invalid Barrier Variation")
      end select

      !...............WRITE SCREEN/UNIT 16 INFORMATION AT START/END
      select case (BAR_VARYTYPE(MYIDX))
      case (1, 2)
         if (abs(BAR_DEG_START(MYIDX) - TIMELOC) <= eps) then
            if (BAR_VARYTYPE(MYIDX) == 1) then
               VC = 'LINEAR'
            elseif (BAR_VARYTYPE(MYIDX) == 2) then
               VC = 'ETA_MAX'
            end if
#ifdef CMPI
            write (ScratchMessage, 1900) trim(VC), 'BEGAN', &
               TIMELOC, NBV(MYIDX), MYPROC
#else
            write (ScratchMessage, 1901) trim(VC), 'BEGAN', &
               TIMELOC, NBV(MYIDX)
#endif
            call allMessage(INFO, scratchMessage)
         end if
         if (abs(BAR_DEG_END(MYIDX) - TIMELOC) <= eps) then
            if (BAR_VARYTYPE(MYIDX) == 1) then
               VC = 'LINEAR'
            elseif (BAR_VARYTYPE(MYIDX) == 2) then
               VC = 'ETA_MAX'
            end if
#ifdef CMPI
            write (ScratchMessage, 1900) trim(VC), 'CONCLUDED', &
               TIMELOC, NBV(MYIDX), MYPROC
#else
            write (ScratchMessage, 1901) trim(VC), 'CONCLUDED', &
               TIMELOC, NBV(MYIDX)
#endif
            call allMessage(INFO, ScratchMessage)
         end if
      end select
#ifdef CMPI
1900  format('INFO: ', A, ' TIME VARYING BOUNDARY ', A, &
             ' CHANGING AT TIME = ', E15.8, ' AT NODE = ', &
             I7, ' ON MYPROC = ', I4)
#else
1901  format('INFO: ', A, ' TIME VARYING BOUNDARY ', A, &
             ' CHANGING AT TIME = ', E15.8, ' AT NODE = ', &
             I7)
#endif

      return

      !-----------------------------------------------------------------------
   end subroutine COMPUTE_BARRIER_HEIGHT
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E
   !       C O M P U T E _ B A R R I E R _ H E I G H T _ L I N E A R
   !-----------------------------------------------------------------------
   !  THIS SUBROUTINE COMPUTES THE NEW BARRIER ELEVATION AT THE
   !  SPECIFIED BOUNDARY NODE. THE INTERPOLATION FROM STARTING
   !  ELEVATION (FORT.14) AND FINAL ELEVATION IS LINEAR
   !  OVER THE SPECIFIED TIME AND DOES NOT STOP ONCE
   !  IT HAS BEGUN.
   !-----------------------------------------------------------------------
   subroutine COMPUTE_BARRIER_HEIGHT_LINEAR(MYIDX, &
                                            TIMELOC, BAR_HEIGHT_CURRENT, BAR_HEIGHT, BARHT_START)

      use MESH, only: DP
      use BOUNDARIES, only: BARINHT, BARLANHT, LBCODEI, &
                            NBV, IBCONN

      implicit none

      integer, intent(IN) :: MYIDX
      real(8), intent(IN) :: TIMELOC
      real(8), intent(IN) :: BAR_HEIGHT_CURRENT
      real(8), intent(OUT) :: BAR_HEIGHT
      real(8), intent(IN), optional :: BARHT_START
      real(8) :: BAR_DZ, BAR_DT
      real(8) :: BAR_START
      real(8) :: DEPTH
      real(8), parameter :: eps = epsilon(1d0)
      character(1024) :: ScratchMessage

      LOG_SCOPE_TRACED("COMPUTE_BARRIER_HEIGHT_LINEAR", WEIR_BOUNDARY_TRACING)

      !...............SET DEFAULT RETURN VALUE
      BAR_HEIGHT = BAR_HEIGHT_CURRENT

      !...............GET OUT OF HERE AT THE FIRST OPPORTUNITY
      if (BAR_VARYTYPE(MYIDX) == 1) then
         if ((BAR_DEG_START(MYIDX) > TIMELOC) .or. &
             (BAR_DEG_END(MYIDX) < TIMELOC)) then
            return
         end if
      end if

      !...............COMPUTE THE NEW BOUNDARY HEIGHT
      if (present(BARHT_START)) then
         !                   !...Schedule style boundaries specify
         !                   !   their own initial elevation, so dont use
         !                   !   the elevation in the fort.14 weirs
         BAR_DZ = BARHT_START - BARHT_FINAL(MYIDX)
         BAR_START = BARHT_START
      else
         select case (LBCODEI(MYIDX))
         case (3, 13, 23)
            BAR_DZ = BARLANHT(MYIDX) - &
                     BARHT_FINAL(MYIDX)
            BAR_START = BARLANHT(MYIDX)
         case (4, 24, 64, 5, 25)
            BAR_DZ = BARINHT(MYIDX) - BARHT_FINAL(MYIDX)
            BAR_START = BARINHT(MYIDX)
         case DEFAULT
            call allMessage(WARNING, &
                            "INVALID BOUNDARY CONDITION SPECIFIED")
            return
         end select
      end if

      if (abs(BAR_DEG_END(MYIDX) - TIMELOC) < eps) then

         BAR_HEIGHT = BARHT_FINAL(MYIDX)
      else
         BAR_DT = BAR_DEG_END(MYIDX) - BAR_DEG_START(MYIDX)
         BAR_HEIGHT = BAR_START - &
                      (BAR_DZ*(1d0 - (BAR_DEG_END(MYIDX) - &
                                      TIMELOC)/BAR_DT))
      end if

      !...............CHECK AGAINST BATHYMETRY, DO NOT LET WEIR ELEVATION DECREASE BELOW
      !               PREVAILING GROUND. CHECK NODE ACROSS BOUNDARY AS WELL AS LOCAL NODE
      DEPTH = max(-DP(NBV(MYIDX)), -DP(IBCONN(MYIDX)))
      if (BAR_HEIGHT < DEPTH) then
         !                   !...WE'RE TOO LOW, SHUT DOWN BARRIER, SET TO MINIMUM
         if (BAR_VARYTYPE(MYIDX) == 3) then
            BAR_HEIGHT = DEPTH
         else
            if (BARHT_FINAL(MYIDX) >= DEPTH) then
               !                           !...SET TO USER SPECIFIED MIN
               BAR_HEIGHT = BARHT_FINAL(MYIDX)
            else
#ifdef CMPI
               write (ScratchMessage, 1904) TIMELOC, NBV(MYIDX), &
                  MYPROC

#else
               write (ScratchMessage, 1905) TIMELOC, NBV(MYIDX)
#endif
               call allMessage(INFO, ScratchMessage)
               BAR_HEIGHT = DEPTH
            end if
            BAR_DEG(MYIDX) = .false.
         end if
      end if

#ifdef CMPI
1904  format('WARNING: BARRIER CANNOT DECREASE TO ', &
             'A VALUE BELOW PREVAILING GROUND.', /, &
             '      BARRIER SET TO MINIMUM VALUE AT ', &
             'TIME = ', E15.8, ' NODE = ', I7, &
             ' ON MYPROC = ', I4)
#else
1905  format('WARNING: BARRIER CANNOT DECREASE TO ', &
             'A VALUE BELOW PREVAILING GROUND.', /, &
             '      BARRIER SET TO MINIMUM VALUE AT ', &
             ' TIME = ', E15.8, ' NODE = ', I7)
#endif

      return

      !-----------------------------------------------------------------------
   end subroutine COMPUTE_BARRIER_HEIGHT_LINEAR
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E
   !       C O M P U T E _ B A R R I E R _ H E I G H T _ E T A M A X
   !-----------------------------------------------------------------------
   !  THIS SUBROUTINE CHECKS TO SEE IF THE BARRIER HAS REACHED
   !  ITS ETAMAX FAILURE CRITERIA. IF IT HAS, THE LINEAR BARRIER
   !  HEIGHT CHANGE ROUTINE IS CALLED WITH THE VARIABLES CREATED
   !  WHEN THIS ROUTINE FIRST DECIDED THAT THE ETAMAX CONDITION
   !  WAS MET
   !-----------------------------------------------------------------------
   subroutine COMPUTE_BARRIER_HEIGHT_ETAMAX(MYIDX, TIMELOC, &
                                            BAR_HEIGHT_CURRENT, BAR_HEIGHT)

      use GLOBAL, only: ETA2, NNODECODE
      use BOUNDARIES, only: NBV, IBCONN

      implicit none

      integer, intent(IN) :: MYIDX
      real(8), intent(IN) :: TIMELOC
      real(8), intent(IN) :: BAR_HEIGHT_CURRENT
      real(8), intent(OUT) :: BAR_HEIGHT
      real(8), parameter :: eps = epsilon(1d0)

      LOG_SCOPE_TRACED("COMPUTE_BARRIER_HEIGHT_ETAMAX", WEIR_BOUNDARY_TRACING)

      !...............SET DEFAULT RETURN VALUE
      BAR_HEIGHT = BAR_HEIGHT_CURRENT

      !...............CHECK IF THE ETA_MAX HAS BEEN EXCEEDED ALREADY
      if (abs(BAR_FAILURE_START(MYIDX) - (-1d0)) < eps) then
         if (ETA2(NBV(MYIDX)) > BAR_ETA_MAX(MYIDX) .and. NNODECODE(NBV(MYIDX)) == 1) then
            !.......................WE NEED TO INITIALIZE OUR START TIME FOR THE FAILURE
            BAR_DEG_START(MYIDX) = TIMELOC
            BAR_DEG_END(MYIDX) = TIMELOC + &
                                 BAR_FAILURE_DURATION(MYIDX)
            BAR_FAILURE_START(MYIDX) = 1d0
         elseif (ETA2(IBCONN(MYIDX)) > BAR_ETA_MAX(MYIDX) .and. NNODECODE(IBCONN(MYIDX)) == 1) then
            !.......................WE NEED TO INITIALIZE OUR START TIME FOR THE FAILURE
            BAR_DEG_START(MYIDX) = TIMELOC
            BAR_DEG_END(MYIDX) = TIMELOC + &
                                 BAR_FAILURE_DURATION(MYIDX)
            BAR_FAILURE_START(MYIDX) = 1d0
         else
            return
         end if
      end if

      !...............IF WE ARE HERE, ETA_MAX WAS PREVIOUSLY EXCEEDED AND THE BARRIER
      !               IS DEGRADING AS PRESCRIBED
      call COMPUTE_BARRIER_HEIGHT_LINEAR(MYIDX, &
                                         TIMELOC, BAR_HEIGHT_CURRENT, BAR_HEIGHT)
      return

      !-----------------------------------------------------------------------
   end subroutine COMPUTE_BARRIER_HEIGHT_ETAMAX
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E
   !       C O M P U T E _ B A R R I E R _ H E I G H T _ S C H E D U L E
   !-----------------------------------------------------------------------
   !  THIS SUBROUTINE COMPUTES THE VARYING WEIR ELEVATION BASED UPON A
   !  INPUT FILE THAT CONTAINS A SCHEDULE. THIS SCHEDULE CAN BE USED TO
   !  HELP REPRESENT THINGS LIKE GATE OPERATIONS AND OTHER DYNAMIC
   !  PROCESSES IN THE MODEL THAT BEHAVE IN A KNOWN REPETITIVE WAY
   !-----------------------------------------------------------------------
   subroutine COMPUTE_BARRIER_HEIGHT_SCHEDULE(MYIDX, TIMELOC, &
                                              BAR_HEIGHT_CURRENT, BAR_HEIGHT)
      use MESH, only: DP
      use BOUNDARIES, only: IBCONN, NBV, LBCODEI, BARLANHT, BARINHT
      implicit none
      integer, intent(IN) :: MYIDX
      real(8), intent(IN) :: TIMELOC
      real(8), intent(IN) :: BAR_HEIGHT_CURRENT
      real(8), intent(OUT) :: BAR_HEIGHT

      real(8) :: PREVEND
      real(8) :: BARHT_START
      real(8) :: OFFSET
      integer :: MYSEC
      integer :: NNBB1, NNBB2
      real(8), parameter :: eps = epsilon(1d0)
      character(1024) :: ScratchMessage

      LOG_SCOPE_TRACED("COMPUTE_BARRIER_HEIGHT_SCHEDULE", WEIR_BOUNDARY_TRACING)

      !...............SET DEFAULT RETURN VALUE
      BAR_HEIGHT = BAR_HEIGHT_CURRENT

      !...............Check if we have reached the offset time yet
      !               (Time added to beginning of a schedule)
      if (TIMELOC < BAR_SCHEDULE(MYIDX)%OFFSET) then
         return
      end if

      !...............SIMPLIFY THESE VARIABLES
      MYSEC = BAR_SCHEDULE(MYIDX)%MYSEC
      PREVEND = BAR_SCHEDULE(MYIDX)%PREVEND
      OFFSET = BAR_SCHEDULE(MYIDX)%OFFSET

      !...............IF THIS IS THE FIRST TIMESTEP IN THIS SECTINO OF THE
      !               SCHEDULE, THEN SET UP THE VARIABLES NEEDED
      if (abs(BAR_SCHEDULE(MYIDX)%SECTION(MYSEC)%BAR_DEG_START - &
              (TIMELOC - PREVEND - OFFSET)) < eps) then

         !...................THIS IS THE BARRIER HEIGHT AT THE START OF THIS
         !                   PORTION OF THE SCHEDULE. SAVE IT FOR LINEAR
         !                   INTERPOLATION
         BARHT_START = BARINHT2(MYIDX)
         BAR_SCHEDULE(MYIDX)%BARHT_START = BARHT_START

         !...................GET THE NODE NUMBER
         NNBB1 = NBV(MYIDX)

         !...................INFORM SCREEN AND UNIT 16 WE ARE HERE
#ifdef CMPI
         write (ScratchMessage, 101) MYSEC, NNBB1, TIMELOC, MYPROC
#else
         write (ScratchMessage, 102) MYSEC, NNBB1, TIMELOC
#endif
         call allMessage(INFO, ScratchMessage)

         !...................ALTER THE VARIABLES USED IN COMPUTE_BARRIER_HEIGHT_LINEAR
         !                   SO IT CAN CORRECTLY INTERPOLATE NEW BARRIER VALUES
         BAR_DEG_START(MYIDX) = &
            BAR_SCHEDULE(MYIDX)%SECTION(MYSEC)%BAR_DEG_START &
            + PREVEND + OFFSET
         BAR_DEG_END(MYIDX) = &
            BAR_SCHEDULE(MYIDX)%SECTION(MYSEC)%BAR_DEG_END &
            + PREVEND + OFFSET

         !...................CHECK THE ZF VALUE FOR FLAGS. THESE FLAGS ARE ADDED
         !                   FOR CONVENIENCE TO SET THE BARRIER TO PREDETERMINED
         !                   VALUES BASED UPON MESH GEOMETRY.
         select case (int(BAR_SCHEDULE(MYIDX)% &
                          SECTION(MYSEC)%BARHT_FINAL))
         case (-99990)
            !                           !...Set to fort.14 nodal elevation
            select case (LBCODEI(MYIDX))
            case (3, 13, 23)
               BARHT_FINAL(MYIDX) = -DP(NNBB1)
            case (4, 14, 24, 64, 5, 25)
               NNBB2 = IBCONN(MYIDX)
               BARHT_FINAL(MYIDX) = &
                  min(-DP(NNBB1), -DP(NNBB2))
            end select
         case (-99991)
            !                           !...Add BARHT_DELTA to present height
            select case (LBCODEI(MYIDX))
            case (3, 13, 23)
               BARHT_FINAL(MYIDX) = &
                  BARLANHT2(MYIDX) + &
                  BAR_SCHEDULE(MYIDX)% &
                  SECTION(MYSEC)%BARHT_DELTA
            case (4, 14, 24, 64, 5, 25)
               BARHT_FINAL(MYIDX) = &
                  BARINHT2(MYIDX) + &
                  BAR_SCHEDULE(MYIDX)% &
                  SECTION(MYSEC)% &
                  BARHT_DELTA
            end select
         case (-99992)
            !                           !...Subtract BARHT_DELTA from present height
            select case (LBCODEI(MYIDX))
            case (3, 13, 23)
               BARHT_FINAL(MYIDX) = &
                  BARLANHT2(MYIDX) - &
                  BAR_SCHEDULE(MYIDX)% &
                  SECTION(MYSEC)%BARHT_DELTA
            case (4, 14, 24, 64, 5, 25)
               BARHT_FINAL(MYIDX) = &
                  BARINHT2(MYIDX) - &
                  BAR_SCHEDULE(MYIDX)% &
                  SECTION(MYSEC)% &
                  BARHT_DELTA
            end select
         case (-99993)
            !                           !...Set to weir elevation in fort.14
            select case (LBCODEI(MYIDX))
            case (3, 13, 23)
               BARHT_FINAL(MYIDX) = BARLANHT(MYIDX)
            case (4, 14, 24, 64, 5, 25)
               BARHT_FINAL(MYIDX) = BARINHT(MYIDX)
            end select
         case DEFAULT
            BARHT_FINAL(MYIDX) = BAR_SCHEDULE(MYIDX)% &
                                 SECTION(MYSEC)%BARHT_FINAL
         end select

         !...................USE THE NORMAL LINEAR INTERPOLATION ROUTINE TO
         !                   GET THE NEW BARRIER ELEVATION
         call COMPUTE_BARRIER_HEIGHT_LINEAR(MYIDX, TIMELOC, &
                                            BAR_HEIGHT_CURRENT, BAR_HEIGHT, BARHT_START)

         !...................THIS IS THE END OF A SCHEDULE SECTION, WE NEED TO
         !                   INCREMENT COUNTERS AND PREPARE FOR THE NEXT SECTION
         !                   OF THE SCHEDULE
      elseif (abs(BAR_SCHEDULE(MYIDX)%SECTION(MYSEC)% &
                  BAR_DEG_END - (TIMELOC - PREVEND - OFFSET)) < eps) then

         !...................GET NODE NUMBER AND INFORM SCREEN AND UNIT 16
         NNBB1 = NBV(MYIDX)
#ifdef CMPI
         write (ScratchMessage, 103) MYSEC, NNBB1, TIMELOC, MYPROC
#else
         write (ScratchMessage, 104) MYSEC, NNBB1, TIMELOC
#endif
         call allMessage(INFO, ScratchMessage)

         !...................CALL THE LAST LINEAR INTERPOLATION
         BARHT_START = BAR_SCHEDULE(MYIDX)%BARHT_START
         call COMPUTE_BARRIER_HEIGHT_LINEAR(MYIDX, TIMELOC, &
                                            BAR_HEIGHT_CURRENT, BAR_HEIGHT, BARHT_START)

         !...................INCREMENT SCHEDULE SECTION
         BAR_SCHEDULE(MYIDX)%MYSEC = MYSEC + 1

         !...................CHECK TO SEE IF THIS IS THE END OF THE SCHEDULE.
         !                   IF THE LOOP FLAG HAS BEEN SET TO -1, WE LOOP THE
         !                   SCHEDULE INFINITELY. IF LOOP SET TO 1, WE WILL
         !                   INCREMENT OUR LOOP COUNTER TO CHECK AND SEE IF WE
         !                   HAVE DONE THE CORRECT NUMBER OF SCHEDULE LOOPS
         if (BAR_SCHEDULE(MYIDX)%MYSEC > &
             BAR_SCHEDULE(MYIDX)%NSECTIONS) then

            if (BAR_SCHEDULE(MYIDX)%LOOP == 1) then
               if (BAR_SCHEDULE(MYIDX)%LOOPSCOMPLETE < &
                   BAR_SCHEDULE(MYIDX)%NLOOPS) then
                  BAR_SCHEDULE(MYIDX)%MYSEC = 1
                  BAR_SCHEDULE(MYIDX)%LOOPSCOMPLETE = &
                     BAR_SCHEDULE(MYIDX)% &
                     LOOPSCOMPLETE + 1
                  BAR_SCHEDULE(MYIDX)%PREVEND = TIMELOC &
                                                - OFFSET
               else
                  BAR_DEG(MYIDX) = .false.
               end if
            elseif (BAR_SCHEDULE(MYIDX)%LOOP == -1) then
               BAR_SCHEDULE(MYIDX)%MYSEC = 1
               BAR_SCHEDULE(MYIDX)%PREVEND = TIMELOC &
                                             - OFFSET
            else
               !...............................THIS NODE SHOULD NO LONGER CHANGE, SO
               !                               SET THE APPROPRIATE FLAG
               BAR_DEG(MYIDX) = .false.
            end if
         end if

         !...............THE CASE FOR THE MIDDLE OF THE SCHEDULE. CONTINUE WITH
         !               THE INTERPOLATION
      elseif ((BAR_SCHEDULE(MYIDX)%SECTION(MYSEC)%BAR_DEG_START &
               < TIMELOC - PREVEND - OFFSET) .and. &
              (BAR_SCHEDULE(MYIDX)%SECTION(MYSEC)%BAR_DEG_END &
               > TIMELOC - PREVEND - OFFSET)) then
         BARHT_START = BAR_SCHEDULE(MYIDX)%BARHT_START
         call COMPUTE_BARRIER_HEIGHT_LINEAR(MYIDX, TIMELOC, &
                                            BAR_HEIGHT_CURRENT, BAR_HEIGHT, BARHT_START)
      end if

      return

#ifdef CMPI
101   format("INFO: SCHEDULE STYLE BOUNDARY BEGAN SEGMENT ", &
             I4, " AT NODE = ", I9, " AT TIME = ", E15.8, &
             "ON MYPROC = ", I4)
#else
102   format("INFO: SCHEDULE STYLE BOUNDARY BEGAN SEGMENT ", &
             I4, " AT NODE = ", I9, " AT TIME = ", E15.8)
#endif

#ifdef CMPI
103   format("INFO: SCHEDULE STYLE BOUNDARY CONCLUDED", &
             " SEGMENT ", I4, " AT NODE = ", I9, " AT TIME = ", &
             E15.8, " ON MYPROC = ", I4)
#else
104   format("INFO: SCHEDULE STYLE BOUNDARY CONDLUDED", &
             " SEGMENT ", I4, " AT NODE = ", I9, " AT TIME = ", &
             E15.8)
#endif
      !-----------------------------------------------------------------------
   end subroutine COMPUTE_BARRIER_HEIGHT_SCHEDULE
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E
   !       P A R S E _ T I M E _ V A R Y I N G _ W E I R _ I N F O
   !-----------------------------------------------------------------------
   !  THIS SUBROUTINE READS THE NAMELIST STYLE INPUT FILE FROM THE USER AND
   !  SETS UP THE NECESSARY INFORMATION FOR USE DURING THE SIMLUATION
   !  AS WELL AS PERFORMING VARIOUS SANITY CHECKS ON THE INPUT TO AVOID
   !  ERRORS DURING THE SIMULATION.
   !-----------------------------------------------------------------------
   subroutine PARSE_TIME_VARYING_WEIR_INFO()

      use GLOBAL, only: ITHS, DTDP, IHOT, tvw_file
      use mod_io, only: openFileForRead
      use SIZES, only: INPUTDIR
      use BOUNDARIES, only: LBCODEI, NBV, IBCONN

      implicit none

      character(2000) :: InputString, modifiedString
      character(1024), allocatable :: SCHEDULE_LIST_RAW(:)
      real(8) :: OFFSET
      integer :: IDX, IDX2
      integer :: I, IOS
      integer :: NTIMEVARYINGWEIRS
      integer :: NSCHEDULES_RAW
      character(1024) :: ScratchMessage

      LOG_SCOPE_TRACED("PARSE_TIME_VARYING_WEIR_INFO", WEIR_BOUNDARY_TRACING)

      call openFileForRead(99, trim(INPUTDIR)//'/'// &
                           trim(tvw_file), IOS, required=.false.)
      if (IOS /= 0) then
         NTIMEVARYINGWEIRS = 0
#ifdef CMPI
         write (ScratchMessage, 101) MYPROC
#else
         write (ScratchMessage, 102)
#endif
         call allMessage(INFO, ScratchMessage)
         return
      end if
      read (99, *) NTIMEVARYINGWEIRS
      if (NTIMEVARYINGWEIRS > 0) then
         call ALLOCATE_TIMEVARYINGWEIRS()
#ifdef CMPI
         write (ScratchMessage, 103) MYPROC
#else
         write (ScratchMessage, 104)
#endif
         call allMessage(INFO, ScratchMessage)
      else
#ifdef CMPI
         write (ScratchMessage, 105) MYPROC
#else
         write (ScratchMessage, 106)
#endif
         call allMessage(INFO, ScratchMessage)
         return
      end if

      open (UNIT=98, FILE=trim(LOCALDIR)//'/namelist.scratch', &
            STATUS='REPLACE', ACTION='WRITE')

      !...............READ IN UNFORMATTED NAMELIST, WRITE TO FORMATTED NAMELIST
      do I = 1, NTIMEVARYINGWEIRS
         !...................MODIFY LINE FOR NAMELIST FORMATTING
         read (99, '(A)') InputString
         modifiedString = "&TimeVaryingWeir "// &
                          trim(adjustl(InputString))//" /"
         write (98, '(A)') trim(adjustl(modifiedString))
      end do
      close (98)

      !...............OPEN FORMATTED NAMELIST FILE, READ IN NUMBER OF SCHEDULES
      !               WE WILL NEED TO PARSE
      open (UNIT=98, FILE=trim(LOCALDIR)//'/namelist.scratch', &
            ACTION="READ")
      NSCHEDULES = 0
      NSCHEDULES_RAW = 0
      do I = 1, NTIMEVARYINGWEIRS
         call NULLIFY_TVW_NML()
         read (98, NML=TimeVaryingWeir, IOSTAT=IOS, ERR=200)
         if (VARYTYPE == 4) cycle ! crest table, read by mod_tvv
         if (.not. ISNULL(S=SCHEDULEFILE)) then
            NSCHEDULES_RAW = NSCHEDULES_RAW + 1
         end if
      end do
      if (NSCHEDULES_RAW > 0) then
         rewind (98)
         IDX = 0
         allocate (SCHEDULE_LIST_RAW(1:NSCHEDULES_RAW))
         !...................READ BACK IN A LIST CONTAINING THE SCHEDULE FILE NAMES
         do I = 1, NTIMEVARYINGWEIRS
            call NULLIFY_TVW_NML()
            read (98, NML=TimeVaryingWeir, IOSTAT=IOS)
            if (VARYTYPE == 4) cycle
            if (.not. ISNULL(S=SCHEDULEFILE)) then
               IDX = IDX + 1
               SCHEDULE_LIST_RAW(IDX) = SCHEDULEFILE
            end if
         end do
         !...................FIND OUT HOW MANY UNIQUE SCHEDULES THERE ARE
         call ALPHABETIZE(SCHEDULE_LIST_RAW)
         call UNIQUE_NAMES(SCHEDULE_LIST_RAW, NSCHEDULES, &
                           SCHEDULE_LIST)
         allocate (SCHEDULE(1:NSCHEDULES))
         do I = 1, NSCHEDULES
            call PARSE_SCHEDULE(SCHEDULE_LIST(I), &
                                SCHEDULE(I))
         end do
      end if

      close (98)
      !...............OPEN FORMATTED NAMELIST FILE, READ IN
      open (UNIT=98, FILE=trim(LOCALDIR)//'/namelist.scratch', &
            ACTION="READ")
      do I = 1, NTIMEVARYINGWEIRS
         !...................READ IN NAMELIST LINE
         call NULLIFY_TVW_NML()
         read (98, NML=TimeVaryingWeir, ERR=200, &
               IOSTAT=IOS)

         !...................Time-varying crest VEW lines are handled by mod_tvv
         if (VARYTYPE == 4) cycle

         if (VARYTYPE == 3) then
            if (ISNULL(TimeStartDay)) TimeStartDay = 0d0
            if (ISNULL(TimeStartHour)) TimeStartHour = 0d0
            if (ISNULL(TimeStartMin)) TimeStartMin = 0d0
            if (ISNULL(TimeStartSec)) TimeStartSec = 0d0
            if (ISNULL(I=LOOP)) LOOP = 0
            if (ISNULL(I=NLOOPS)) NLOOPS = 0
            OFFSET = TimeStartDay*86400d0 + &
                     TimeStartHour*3600d0 + &
                     TimeStartMin*60d0 + &
                     TimeStartSec
         else
            OFFSET = 0d0
         end if

         !...................BEGIN SANITY CHECK ON WHAT HAS BEEN SPECIFIED IN NAMELIST
         select case (IHOT)
         case (17, 67, 68, 367, 368, 567, 568)
            if (ISNULL(I=HOT) .or. (HOT == 0)) then
               HOTADD = 0d0
            elseif (HOT == 1) then
               HOTADD = DTDP*dble(ITHS)
               OFFSET = OFFSET + HOTADD
            else
               call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                              message="INCORRECT HOT START VALUE SPECIFIED. HOT=1 OR "// &
                              "HOT=0 FOR HOT START RELATIVE TIME "// &
                              "VARYING WEIRS.")
            end if
         case DEFAULT
            HOTADD = 0d0
         end select
         if (ISNULL(I=VARYTYPE)) then
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message="YOU MUST SPECIFY VARYTYPE= IN THE FORT.15 INPUT "// &
                           "FILE.")
         end if
         select case (VARYTYPE)
         case (1, 2, 3)
            if (VARYTYPE == 1) then
               if (ISNULL(X1) .or. ISNULL(Y1) .or. &
                   ISNULL(ZF)) then
                  call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                                 message="YOU MUST SPECIFY X1=, Y1=, and ZF= ")
               end if
               if (ISNULL(TimeStartday) .and. &
                   ISNULL(TimeStartHour) .and. &
                   ISNULL(TimeStartSec)) then
                  call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                                 message="TIME VARYING BOUNDARY START TIME"// &
                                 " MUST BE SPECIFIED.")
               end if
               if (ISNULL(TimeEndday) .and. &
                   ISNULL(TimeEndHour) .and. &
                   ISNULL(TimeEndSec)) then
                  call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                                 message="TIME VARYING BOUNDARY END TIME"// &
                                 " MUST BE SPECIFIED.")
               end if
               if (ISNULL(TimeStartDay)) TimeStartDay = 0d0
               if (ISNULL(TimeStartHour)) TimeStartHour = 0d0
               if (ISNULL(TimeStartMin)) TimeStartMin = 0d0
               if (ISNULL(TimeStartSec)) TimeStartSec = 0d0
               if (ISNULL(TimeEndDay)) TimeEndDay = 0d0
               if (ISNULL(TimeEndHour)) TimeEndHour = 0d0
               if (ISNULL(TimeEndMin)) TimeEndMin = 0d0
               if (ISNULL(TimeEndSec)) TimeEndSec = 0d0
            elseif (VARYTYPE == 2) then
               if (ISNULL(X1) .or. ISNULL(Y1) .or. &
                   ISNULL(ZF)) then
                  call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                                 message="YOU MUST SPECIFY X1=, Y1=, and ZF= ")
               end if
               if ((ISNULL(FAILUREDURATIONDAY) .and. &
                    ISNULL(FAILUREDURATIONHOUR) .and. &
                    ISNULL(FAILUREDURATIONMIN) .and. &
                    ISNULL(FAILUREDURATIONSEC)) .or. &
                   ISNULL(ETA_MAX)) then
                  call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                                 message="YOU MUST SPECIFY A FAILURE DURATION AND "// &
                                 "MAXIMUM WATER SURFACE BEFORE ELEVATION CHANGE.")
               end if
               if (ISNULL(FAILUREDURATIONDAY)) &
                  FAILUREDURATIONDAY = 0d0
               if (ISNULL(FAILUREDURATIONHOUR)) &
                  FAILUREDURATIONHOUR = 0d0
               if (ISNULL(FAILUREDURATIONMIN)) &
                  FAILUREDURATIONMIN = 0d0
               if (ISNULL(FAILUREDURATIONSEC)) &
                  FAILUREDURATIONSEC = 0d0
            elseif (VARYTYPE == 3) then
               if (ISNULL(X1) .or. ISNULL(Y1)) then
                  call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                                 message="YOU MUST SPECIFY X1= and Y1=")
               end if
               if (ISNULL(S=SCHEDULEFILE)) then
                  call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                                 message="YOU MUST SPECIFY SCHEDULEFILE= FOR VARYTYPE=3.")
               end if
               if (ISNULL(TimeStartDay)) TimeStartDay = 0d0
               if (ISNULL(TimeStartHour)) TimeStartHour = 0d0
               if (ISNULL(TimeStartMin)) TimeStartMin = 0d0
               if (ISNULL(TimeStartSec)) TimeStartSec = 0d0
               if (ISNULL(TimeEndDay)) TimeEndDay = 0d0
               if (ISNULL(TimeEndHour)) TimeEndHour = 0d0
               if (ISNULL(TimeEndMin)) TimeEndMin = 0d0
               if (ISNULL(TimeEndSec)) TimeEndSec = 0d0
            end if
            !.......................END SANITY CHECK

            !.......................BEGIN ASSIGNING BOUNDARY CONDITIONS THEIR ATTRIBUTES
            call FIND_BOUNDARY_NODES(X1, Y1, IDX)
            if (IDX == -1) cycle
            select case (LBCODEI(IDX))
            case (3, 13, 23)
               !...........................A ONE SIDED STYLE WEIR, WE DONT NEED X2 AND Y2
               call ASSIGN_TVW_TIMING(VARYTYPE, IDX)
               EXT_TVW = .true.
               if (VARYTYPE == 1) then
                  if (BAR_DEG_END(IDX) <= 0d0) then
                     call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                                    message="TIME VARYING BOUNDARY END TIME MUST BE "// &
                                    "GREATER THAN ZERO.")
                  end if
               elseif (VARYTYPE == 2) then
                  if (BAR_FAILURE_DURATION(IDX) <= 0d0) then
                     call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                                    message="FAILURE_DURATION MUST BE GREATER THAN ZERO")
                  end if
               elseif (VARYTYPE == 3) then
                  BAR_SCHEDULE(IDX)%OFFSET = OFFSET
               end if

            case (4, 24, 64, 5, 25)
               !...........................A TWO SIDED STYLE WEIR, WE NEED X2 AND Y2
               call FIND_BOUNDARY_NODES(X2, Y2, IDX2)
               if (IDX2 == -1) then
                  call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                                 message="BOUNDARY CONDITION IMPROPERLY SPLIT INTO "// &
                                 "GHOST NODE SPACE!")
               end if
               INT_TVW = .true.

               !...........................WE NEED TO CHECK CONNECTIVITY ACCROSS WEIR TO
               !                           MAKE SURE THIS PAIR IS VALID
               if (NBV(IDX) /= IBCONN(IDX2)) then
                  !...............................THIS IS NOT THE CONNECTIVITY WE EXPECTED, GOOD NIGHT
                  write (ScratchMessage, '(A,F0.9,A,F0.6,A,I0,A,F0.9,A,F0.6,A,I0)') &
                     "CONNECTIVITY ACROSS INTERNAL TIME VARYING BOUNDARY IS INVALID. "// &
                     "X1=", X1, " Y1=", Y1, " LOCAL NODE=", NBV(IDX), &
                     " X2=", X2, " Y2=", Y2, " LOCAL NODE=", NBV(IDX2)
                  call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                                 message=trim(ScratchMessage))
               end if
               call ASSIGN_TVW_TIMING(VARYTYPE, IDX, IDX2)
               if (VARYTYPE == 1) then
                  if (BAR_DEG_START(IDX) > &
                      BAR_DEG_END(IDX)) then
                     call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                                    message="INVALID BARRIER DEGREDATION TIME.")
                  end if
               elseif (VARYTYPE == 2) then
                  if ((BAR_FAILURE_DURATION(IDX) <= 0d0) &
                      .or. (BAR_FAILURE_DURATION(IDX2) <= &
                            0d0)) then
                     call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                                    message="FAILURE_DURATION MUST BE GREATER THAN ZERO.")
                  end if
               elseif (VARYTYPE == 3) then
                  BAR_SCHEDULE(IDX)%OFFSET = OFFSET
                  BAR_SCHEDULE(IDX2)%OFFSET = OFFSET
               end if
            case DEFAULT
               !...........................WE SHOULD NOT BE HERE?
               call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                              message="INVALID BOUNDARY CONDITION DETECTED.")
            end select

         case DEFAULT
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message="INVALID WEIR VARIATION SPECIFIED.")
         end select

      end do

      !...............FREE THE MEMORY USED FOR THE KDTREE2
      call KDTREE2_DESTROY(BARRIER_SEARCHTREE)
      deallocate (BAR_LOCATIONS)
#ifdef CMPI
      call KDTREE2_DESTROY(GHOST_SEARCHTREE)
      deallocate (GHOST_LOCATIONS)
#endif

      return
200   call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                     message="PROBLEM READING TIME VARYING WEIR FILE.")

      !...............GENERAL MESSAGES FOR TIME VARYING WEIRS
#ifdef CMPI
101   format("Time varying weir file was not found. All weirs" &
             , " will be static on MYPROC = ", I6)
103   format("Time varying weir file was found. Time varying", &
             " weirs have been specified on MYPROC = ", I6)
105   format("Time varying weir file was found. No time ", &
             "varying weirs specified on MYPROC = ", I6)
#else
102   format("Time varying weir file was not found. All weirs" &
             , " will be static.")

104   format("Time varying weir file was found. Time varying", &
             " weirs have been specified.")

106   format("Time varying weir file was found. No time ", &
             "varying weirs specified.")
#endif
      !-----------------------------------------------------------------------
   end subroutine PARSE_TIME_VARYING_WEIR_INFO
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E  P A R S E _ S C H E D U L E
   !-----------------------------------------------------------------------
   !  THIS SUBROUTINE READS A FILE THAT CONTAINS THE NAMELIST STYPE INPUT
   !  FOR A WEIR THAT VARIES BASED UPON A SCHEDULE THAT IS SPECIFIED IN THE
   !  TIME VARYING WEIR INPUT FILE
   !-----------------------------------------------------------------------

   subroutine PARSE_SCHEDULE(MYFILE, MYSCHED)
      use GLOBAL, only: DTDP, IHOT, ITHS
      implicit none
      character(*), intent(IN) :: MYFILE
      type(BARRIER_SCHEDULE_T), intent(OUT) :: MYSCHED

      character(500) :: origLine, modLine
      real(8) :: BAR_DEG_START, BAR_DEG_END
      real(8) :: DELTA, HOTSEC
      integer :: I
      logical :: exists
      real(8), parameter :: eps = epsilon(1d0)

      namelist /SCHEDULE/ &
         TimeStartDay, TimeStartHour, TimeStartMin, &
         TimeStartSec, TimeEndDay, TimeEndHour, &
         TimeEndMin, TimeEndSec, ZF, DELTA, HOT

      LOG_SCOPE_TRACED("PARSE_SCHEDULE", WEIR_BOUNDARY_TRACING)

      inquire (FILE=trim(GLOBALDIR)//"/"//trim(MYFILE), &
               EXIST=EXISTS)
      if (EXISTS) then
         open (FILE=trim(GLOBALDIR)//"/"//trim(MYFILE), &
               UNIT=99, ACTION="READ")
         read (99, *) MYSCHED%NSECTIONS
         allocate (MYSCHED%SECTION(1:MYSCHED%NSECTIONS))
         open (FILE=trim(LOCALDIR)//"/namelist2.scratch", &
               UNIT=97, ACTION="WRITE")
         do I = 1, MYSCHED%NSECTIONS
            read (99, '(A)') OrigLine
            write (modLine, '(3A)') "&Schedule ", &
               trim(OrigLine), " /"
            write (97, '(A)') trim(modLine)
         end do
         close (97)
         close (99)
         open (FILE=trim(LOCALDIR)//"/namelist2.scratch", &
               UNIT=97, ACTION="READ")
         do I = 1, MYSCHED%NSECTIONS
            TimeStartDay = 0d0
            TimeStartHour = 0d0
            TimeStartMin = 0d0
            TimeStartSec = 0d0
            TimeEndDay = 0d0
            TimeEndHour = 0d0
            TimeEndMin = 0d0
            TimeEndSec = 0d0
            HOT = 0
            ZF = -99999d0
            DELTA = -99999d0
            read (97, NML=SCHEDULE)
            if (HOT == 1) then
               select case (IHOT)
               case (17, 67, 68, 367, 368, 567, 568)
                  HOTSEC = DTDP*dble(ITHS)
               case DEFAULT
                  HOTSEC = 0d0
               end select
            else
               HOTSEC = 0d0
            end if
            BAR_DEG_START = TimeStartDay*86400d0 + &
                            TimeStartHour*3600d0 + &
                            TimeStartMin*60d0 + &
                            TimeStartSec + &
                            HOTSEC
            BAR_DEG_END = TimeEndDay*86400d0 + &
                          TimeEndHour*3600d0 + &
                          TimeEndMin*60d0 + &
                          TimeStartSec + &
                          HOTSEC
            MYSCHED%SECTION(I)%BAR_DEG_START = BAR_DEG_START
            MYSCHED%SECTION(I)%BAR_DEG_END = BAR_DEG_END
            MYSCHED%SECTION(I)%BARHT_FINAL = ZF
            MYSCHED%SECTION(I)%BARHT_DELTA = DELTA
            if ((BAR_DEG_START < 0d0) .or. &
                (abs(BAR_DEG_END - 0d0) < eps)) then
               call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                              message="You must specify start and end for each point "// &
                              "in barrier schedule.")
            end if
            if (ISNULL(ZF)) then
               call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                              message="You must specify a flag or final elevation for "// &
                              "BARHT_FINAL")
            elseif ((abs(ZF - (-99991d0)) < eps) .or. &
                    (abs(ZF - (-99992d0)) < eps)) then
               if (ISNULL(DELTA)) then
                  call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                                 message="If specifying a relative value for ZF, "// &
                                 "you must specify DELTA = ")
               end if
            end if
         end do
         close (97, STATUS="DELETE")
      end if
      return
      !-----------------------------------------------------------------------
   end subroutine PARSE_SCHEDULE
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E  A S S I G N _ T V W _ T I M I N G
   !-----------------------------------------------------------------------
   !  THIS SUBROUTINE ASSIGNS THE TIMING PRESENT IN THE JUST READ IN
   !  NAMELIST TO THE CORRECT LOCATIONS IN THE BARRIER VARIATION
   !  ARRAYS
   !-----------------------------------------------------------------------
   subroutine ASSIGN_TVW_TIMING(VAR, INDEX1, INDEX2)
      integer, intent(IN) :: VAR
      integer, intent(IN) :: INDEX1
      integer, intent(IN), optional :: INDEX2
      integer :: I

      LOG_SCOPE_TRACED("ASSIGN_TVW_TIMING", WEIR_BOUNDARY_TRACING)

      if (VAR == 1) then
         BAR_DEG_START(INDEX1) = &
            TimeStartDay*86400d0 + &
            TimeStartHour*3600d0 + &
            TimeStartMin*60d0 + &
            TimeStartSec + &
            HOTADD
         BAR_DEG_END(INDEX1) = &
            TimeEndDay*86400d0 + &
            TimeEndHour*3600d0 + &
            TimeEndMin*60d0 + &
            TimeEndSec + &
            HOTADD
         BAR_DEG(INDEX1) = .true.
         BAR_VARYTYPE(INDEX1) = VARYTYPE
         BARHT_FINAL(INDEX1) = ZF
         if (present(INDEX2)) then
            BAR_DEG_START(INDEX2) = BAR_DEG_START(INDEX1)
            BAR_DEG_END(INDEX2) = BAR_DEG_END(INDEX1)
            BAR_DEG(INDEX2) = .true.
            BAR_VARYTYPE(INDEX2) = VARYTYPE
            BARHT_FINAL(INDEX2) = ZF
         end if
      elseif (VAR == 2) then
         BAR_DEG(INDEX1) = .true.
         BAR_VARYTYPE(INDEX1) = VARYTYPE
         BAR_ETA_MAX(INDEX1) = ETA_MAX
         BARHT_FINAL(INDEX1) = ZF
         BAR_FAILURE_DURATION(INDEX1) = &
            FAILUREDURATIONDAY*86400d0 + &
            FAILUREDURATIONHOUR*3600d0 + &
            FAILUREDURATIONMIN*60d0 + &
            FAILUREDURATIONSEC
         if (present(INDEX2)) then
            BAR_FAILURE_DURATION(INDEX2) = &
               BAR_FAILURE_DURATION(INDEX1)
            BAR_DEG(INDEX2) = .true.
            BAR_VARYTYPE(INDEX2) = VARYTYPE
            BAR_ETA_MAX(INDEX2) = ETA_MAX
            BARHT_FINAL(INDEX2) = ZF
         end if
      elseif (VAR == 3) then
         BAR_DEG(INDEX1) = .true.
         BAR_VARYTYPE(INDEX1) = VARYTYPE
         do I = 1, NSCHEDULES
            if (trim(adjustl(SCHEDULEFILE)) == &
                trim(adjustl(SCHEDULE_LIST(I)))) then
               BAR_SCHEDULE(INDEX1)%MYFILE => &
                  SCHEDULE_LIST(I)
               BAR_SCHEDULE(INDEX1)%SID = I
               BAR_SCHEDULE(INDEX1)%MYSEC = 1
               BAR_SCHEDULE(INDEX1)%LOOP = LOOP
               BAR_SCHEDULE(INDEX1)%NLOOPS = NLOOPS
               BAR_SCHEDULE(INDEX1)%NSECTIONS => &
                  SCHEDULE(I)%NSECTIONS
               BAR_SCHEDULE(INDEX1)%PREVEND = 0d0
               BAR_SCHEDULE(INDEX1)%SECTION => &
                  SCHEDULE(I)%SECTION
               if (present(INDEX2)) then
                  BAR_DEG(INDEX2) = .true.
                  BAR_VARYTYPE(INDEX2) = VARYTYPE
                  BAR_SCHEDULE(INDEX2)%MYFILE => &
                     SCHEDULE_LIST(I)
                  BAR_SCHEDULE(INDEX2)%SID = I
                  BAR_SCHEDULE(INDEX2)%MYSEC = 1
                  BAR_SCHEDULE(INDEX2)%LOOP = LOOP
                  BAR_SCHEDULE(INDEX2)%NLOOPS = NLOOPS
                  BAR_SCHEDULE(INDEX2)%NSECTIONS => &
                     SCHEDULE(I)%NSECTIONS
                  BAR_SCHEDULE(INDEX2)%PREVEND = 0d0
                  BAR_SCHEDULE(INDEX2)%SECTION => &
                     SCHEDULE(I)%SECTION
               end if
               exit
            end if
            if (I == NSCHEDULES) then
               call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                              message="Schedule not found.")
            end if
         end do
      end if

      return

      !-----------------------------------------------------------------------
   end subroutine ASSIGN_TVW_TIMING
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     F U N C T I O N  I S N U L L
   !-----------------------------------------------------------------------
   !  THIS FUNCTION CHECKS TO SEE IF A NAMELIST INPUT VARIABLE WAS MODIFIED
   !  WHEN THE NAMELIST WAS READ, MEANING IT WAS PRESENT IN THE NAMELIST
   !  THAT WAS READ IN. RETURNS TRUE IF THE VARIABLE WAS NOT ALTERED NAD
   !  FALSE IF IT WAS.
   !-----------------------------------------------------------------------
   pure logical function ISNULL(R, I, S)
      implicit none
      real(8), intent(IN), optional :: R
      integer, intent(IN), optional :: I
      character(*), intent(IN), optional :: S
      real(8), parameter :: EPS = epsilon(1d0)

      ISNULL = .false.

      if (present(R)) then
         if (abs(R + 99999d0) <= EPS) then
            ISNULL = .true.
         end if
      elseif (present(I)) then
         if (I == -99999) then
            ISNULL = .true.
         end if
      elseif (present(S)) then
         if (trim(adjustl(S)) == "NOFILE") then
            ISNULL = .true.
         end if
      end if
      !-----------------------------------------------------------------------
   end function ISNULL
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E  N U L L I F Y _ T V W _ N M L
   !-----------------------------------------------------------------------
   !  THIS SUBROUTINE WILL SET THE NAMELIST VARIABLES TO THEIR NULL VALUES
   !  SO THEY CAN BE CHECKED FOR MODIFICATION LATER
   !-----------------------------------------------------------------------
   subroutine NULLIFY_TVW_NML()
      implicit none

      LOG_SCOPE_TRACED("NULLIFY_TVW_NML", WEIR_BOUNDARY_TRACING)

      X1 = -99999d0
      X2 = -99999d0
      Y1 = -99999d0
      Y2 = -99999d0
      ZF = -99999d0
      HOT = -99999
      ETA_MAX = -99999d0
      TimeStartDay = -99999d0
      TimeStartHour = -99999d0
      TimeStartMin = -99999d0
      TimeStartSec = -99999d0
      TimeEndDay = -99999d0
      TimeEndHour = -99999d0
      TimeEndMin = -99999d0
      TimeEndSec = -99999d0
      FailureDurationDay = -99999d0
      FailureDurationHour = -99999d0
      FailureDurationMin = -99999d0
      FailureDurationSec = -99999d0
      VARYTYPE = -99999
      LOOP = -99999
      NLOOPS = -99999
      SearchRadius = -99999d0
      ScheduleFile = "NOFILE"

      return
      !-----------------------------------------------------------------------
   end subroutine NULLIFY_TVW_NML
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E  A L P H A B E T I Z E
   !-----------------------------------------------------------------------
   !  THIS IS AN ALPHABETICAL BUBBLE SORT TO PLACE THE NAMES IN THE INPUT
   !  ARRAY IN ALPHABETICAL ORDER SO THAT IT IS EASIER TO FIND THE TOTAL
   !  NUMBER OF UNIQUE STRINGS THAT HAVE BEEN INPUT TO THE CODE
   !-----------------------------------------------------------------------
   pure subroutine ALPHABETIZE(MYCHAR)
      implicit none
      character(1024), intent(INOUT) :: MYCHAR(:)
      character(1024) :: A, B
      integer :: I, J, K

      do I = 1, size(MYCHAR)
         do J = I + 1, size(MYCHAR)
            A = adjustl(MYCHAR(I))
            B = adjustl(MYCHAR(J))
            SORTING: do K = 1, len_trim(A)
               if (K > len_trim(B)) then
                  !...Flip, B should come before A
                  MYCHAR(J) = A
                  MYCHAR(I) = B
                  exit SORTING
               else
                  if (A(K:K) > B(K:K)) then
                     !...Flip, B should come before A
                     MYCHAR(J) = A
                     MYCHAR(I) = B
                     exit SORTING
                  elseif (A(K:K) < B(K:K)) then
                     !...No Flip. Order is correct
                     exit SORTING
                  end if
               end if
            end do SORTING
         end do
      end do
   end subroutine ALPHABETIZE
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E  U N I Q U E _ N A M E S
   !-----------------------------------------------------------------------
   !  THIS SUBROUTINE WILL TAKE AN ALPHABETICAL LIST OF CHARACTER STRINGS
   !  AS INPUT AND RETURN A LIST WITH ONLY UNIQUE STRINGS AND A COUNT OF
   !  THE NUMBER OF UNIQUE STRINGS
   !-----------------------------------------------------------------------
   pure subroutine UNIQUE_NAMES(RAW_LIST, NUNIQUE, UNIQUE_LIST)
      implicit none
      character(1024), intent(IN) :: RAW_LIST(:)
      character(1024), intent(OUT), allocatable :: UNIQUE_LIST(:)
      integer, intent(OUT) :: NUNIQUE

      character(1024) :: PREV, CURR
      integer :: I, J

      if (size(RAW_LIST) < 2) then
         NUNIQUE = 1
         allocate (UNIQUE_LIST(1:size(RAW_LIST)))
         UNIQUE_LIST = RAW_LIST
         return
      end if

      PREV = RAW_LIST(1)
      NUNIQUE = 1
      do I = 2, size(RAW_LIST)
         CURR = RAW_LIST(I)
         if (trim(PREV) /= trim(CURR)) then
            NUNIQUE = NUNIQUE + 1
            PREV = CURR
         end if
      end do
      allocate (UNIQUE_LIST(1:NUNIQUE))
      J = 1
      PREV = RAW_LIST(1)
      UNIQUE_LIST(1) = PREV
      do I = 2, size(RAW_LIST)
         CURR = RAW_LIST(I)
         if (trim(PREV) /= trim(CURR)) then
            J = J + 1
            PREV = CURR
            UNIQUE_LIST(J) = CURR
         end if
      end do

      !-----------------------------------------------------------------------
   end subroutine UNIQUE_NAMES
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
end module TIME_VARYING_WEIR_BOUNDARY
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     M O D U L E  W E I R _ F L U X
!-----------------------------------------------------------------------
!   THIS MODULE CONTAINS THE ROUTINES FOR THE WEIR OVERTOPPING FLUX
!   CALCULATIONS FOR INTERNAL AND EXTERNAL WEIR TYPE BOUNDARIES.
!-----------------------------------------------------------------------
module WEIR_FLUX
   use ADC_CONSTANTS, only: G
   use GLOBAL, only: ETA2, RAMPINTFLUX
   use mod_logging, only: t_log_scope, init_log_scope, DEBUG, allMessage
   use BOUNDARIES, only: LBCODEI, NBV
   use TIME_VARYING_WEIR_BOUNDARY, only: COMPUTE_BARRIER_HEIGHT, &
                                         BAR_DEG
   use WEIR, only: BARMIN, EXT_TVW, BARLANHT1, BARLANHT2, BARINHT1, BARINHT2, &
                   INT_TVW, EXT_TVW, BARSLIM64_ELEM, BARSLIM64_EDGE, &
                   BARMIN64_SUBM, BARMIN64_SLIM, BARMIN64_NOSUBM, BARMIN64_CORNER, &
                   BARMIN64, BARMIN64_SLIM, WEIR_ENERGY_HEAD, &
                   VEW_CREST_HYDRAULICS, VEW_FROUDE_CRIT, VEW_FROUDE_RECOUPLE, WEIR_SEMI_IMPLICIT

   implicit none

   private

   !...Semi-implicit weir pairs (WEIR_SEMI_IMPLICIT): per boundary entry, active flag, nodal boundary
   !   length (lumped), partner entry
   logical, allocatable :: SW_ACT(:)
   real(8), allocatable :: SW_LEN(:)
   integer, allocatable :: SW_PARTNER(:)
   logical :: SW_ANY = .false.

   public :: COMPUTE_EXTERNAL_BOUNDARY_FLUX, &
             COMPUTE_INTERNAL_BOUNDARY_FLUX, &
             COMPUTE_CROSS_BARRIER_PIPE_FLUX, &
             COMPUTE_INTERNAL_BOUNDARY64_FLUX, &
             SET_SUBMERGED64_AT, RESTART_VEW64_FLUXES, &
             SEMI_WEIR_PREPARE, SEMI_WEIR_SOLVE, SW_ACT, SW_ANY

   !-----------------------------------------------------------------------
contains
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E
   !       C O M P U T E _ E X T E R N A L _ B O U N A R Y _ F L U X
   !-----------------------------------------------------------------------
   !  THIS SUBROUTINE COMPUTES THE OVERTOPPING FLUX OUT OF THE DOMAIN FOR
   !  TYPE 3,13,23 BOUNDARY CONDITIONS
   !-----------------------------------------------------------------------
   subroutine COMPUTE_EXTERNAL_BOUNDARY_FLUX( &
      BARRIER_INDEX, TIMELOC, FLUX)

      use GLOBAL, only: TVW
      use BOUNDARIES, only: BARLANCFSP, BARLANHT
      implicit none
      integer, intent(IN) :: BARRIER_INDEX
      integer :: NNBB
      real(8), intent(IN) :: TIMELOC
      real(8), intent(OUT) :: FLUX
      real(8) :: RBARWL

      LOG_SCOPE_TRACED("COMPUTE_EXTERNAL_BOUNDARY_FLUX", WEIR_BOUNDARY_TRACING)

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
      return

      !-----------------------------------------------------------------------
   end subroutine COMPUTE_EXTERNAL_BOUNDARY_FLUX
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E
   !       C O M P U T E _ I N T E R N A L _ B O U N D A R Y _ F L U X
   !-----------------------------------------------------------------------
   !     THIS SUBROUTINE CONTAINS THE CALCULATIONS OF OVERTOPPING FLUX
   !     PASSED BETWEEN THE THE WEIR PAIR NODES. THIS IS THE VERSION OF
   !     THE CODE THAT CHECKS IF THE WEIR HAS A WET EDGE BEFORE PASSING
   !     FLUX OVER THE WEIR
   !-----------------------------------------------------------------------
   subroutine COMPUTE_INTERNAL_BOUNDARY_FLUX( &
      BARRIER_INDEX, BOUNDARYNODE, BOUNDARY, TIMELOC, FLUX)

      use BOUNDARIES, only: BARINCFSB, BARINCFSP, &
                            NVELL, IBCONN, BARINHT, &
                            CSII, SIII
      use GLOBAL, only: NIBNODECODE, TVW, NODECODE, UU1, VV1
      implicit none
      integer, intent(IN), target :: BARRIER_INDEX
      integer, intent(IN), target :: BOUNDARY
      integer, intent(IN), target :: BOUNDARYNODE
      real(8), intent(IN) :: TIMELOC
      real(8), intent(OUT) :: FLUX

      integer :: NNBB1, NNBB2
      integer :: NNBB1WN, NNBB2WN
      integer :: FLOWDIR
      integer, pointer :: I, J, K
      real(8) :: RBARWL1, RBARWL2
      real(8) :: RBARWL1F, RBARWL2F
      real(8) :: RBARWE1, RBARWE2, VN1, VN2 !...Heads for the flow direction and the formula

      LOG_SCOPE_TRACED("COMPUTE_INTERNAL_BOUNDARY_FLUX", WEIR_BOUNDARY_TRACING)

      !...............SIMPLIFY VARIABLES
      I => BARRIER_INDEX
      J => BOUNDARYNODE
      K => BOUNDARY

      !..............CALL THE ORIGINAL IMPLEMENTATION IF THE USER HAS REQUESTED
      !              IT VIA THE COMPILER FLAG -DORIGWEIR. IF NOT, PROCEED WITH
      !              THE DEFAULT FORMULATION
#ifdef ORIGWEIR
      call COMPUTE_INTERNAL_BOUNDARY_FLUX_ORIG(I, J, K, &
                                               TIMELOC, FLUX)
      return
#endif

      NNBB1 = NBV(I) ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
      NNBB2 = IBCONN(I) ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
      NNBB1WN = 0 ! COUNT NUMBER OF WET NEIGHBORS
      NNBB2WN = 0 ! COUNT NUMBER OF WET NEIGHBORS
      FLOWDIR = 0 ! DIRECTION OF FLOW

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

      !..............CHECK TO SEE IF THE BARRIER ELEVATION NEEDS UPDATING
      if (INT_TVW) then
         if (BAR_DEG(I)) then
            call COMPUTE_BARRIER_HEIGHT(I, TIMELOC, BARINHT1(I) &
                                        , BARINHT2(I))
         end if
         TVW(NNBB1) = BARINHT2(I) - BARINHT(I)
      end if

      !..............GET WATER LEVEL ABOVE WEIR ON EACH SIDE TO COMPUTE HEAD
      !              DEFINE COMPILER FLAG -DAVERAGEWEIRFLOW TO USE BARAVGWT,
      !              WHICH GENERALLY IS SET TO ZERO, AND IS HARD CODED
      !              TO ZERO HERE. READ_INPUT.F CONTAINS THE SPECIFICATION OF
      !              BARAVGWT. WITH BARAVGWT SET TO ZERO, THERE IS NO NEED FOR
      !              IBSTART EITHER.
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
!..............ZC - STREAMLINE THE PROCESS SINCE BARAVGWT IS ZERO BY DEFAULT
      RBARWL1 = ETA2(NNBB1) - BARINHT2(I)
      RBARWL2 = ETA2(NNBB2) - BARINHT2(I)
#endif
      !..............Total energy head (WEIR_ENERGY_HEAD): add the velocity head of the raw nodal
      !              velocity normal to the weir (CSII/SIII of this entry; only vn**2 is used).
      !              The static heads are kept for the checks against 0/BARMIN and for the
      !              downstream depth in the subcritical formula.
      RBARWE1 = RBARWL1
      RBARWE2 = RBARWL2
      if (WEIR_ENERGY_HEAD) then
         VN1 = UU1(NNBB1)*CSII(I) + VV1(NNBB1)*SIII(I)
         VN2 = UU1(NNBB2)*CSII(I) + VV1(NNBB2)*SIII(I)
         RBARWE1 = RBARWL1 + VN1*VN1/(2.d0*G)
         RBARWE2 = RBARWL2 + VN2*VN2/(2.d0*G)
      end if

      RBARWL1F = 2.d0*RBARWE1/3.d0
      RBARWL2F = 2.d0*RBARWE2/3.d0
      FLUX = 0.d0

      if ((RBARWL1 < 0.d0) .and. (RBARWL2 < 0.d0)) then
         !...............WATER LEVEL ON BOTH SIDES OF BARRIER BELOW BARRIER -> CASE 1
         FLUX = 0.d0
      elseif (abs(RBARWE1 - RBARWE2) < 0.01d0) then
         !...............WATER LEVEL EQUAL ON BOTH SIDES OF BARRIER
         !................TO WITHIN TOLERANCE BARMIN -> CASE 2
         FLUX = 0.d0
      elseif ((RBARWE1 > RBARWE2) .and. (RBARWL1 > BARMIN)) then
         !...............WATER LEVEL GREATER ON THIS SIDE OF THE BARRIER AND IS SUCH
         !................THAT OVERTOPPING IS OCCURING
         !................THUS THIS SIDE IS THE SOURCE SIDE OF THE FLOW ACROSS THE BARRIER
         !................NOTE THAT WE DO NOT FORCE THE SOURCE SIDE OF THE BARRIER TO
         !................REMAIN WET. ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
         !................BARRIER HAS BEEN DRIED. IF IT HAS WE SHUT DOWN THE FLOW ACROSS
         !................THE BARRIER
         !................ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE BARRIER HAS AT
         !................LEAST ONE WET EDGE. IF NOT, WE SHUT DOWN THE FLOW ACROSS
         !................THE BARRIER. shintaro v46.28.sb05.05 11/01/2006
         if (RBARWL2 > RBARWL1F) then
            !..................OUTWARD SUBCRITICAL FLOW -> CASE 3
            if (NODECODE(NNBB1) == 0 .or. NNBB1WN == 0) then
               FLUX = 0.0d0
            else
               FLUX = -RampIntFlux*RBARWL2*BARINCFSB(I)* &
                      (2.d0*G*(RBARWE1 - RBARWL2))**0.5d0
            end if
         else
            !..................OUTWARD SUPERCRITICAL FLOW -> CASE 4
            if (NODECODE(NNBB1) == 0 .or. NNBB1WN == 0) then
               FLUX = 0.0d0
            else
               FLUX = -RampIntFlux*BARINCFSP(I)*RBARWL1F* &
                      (RBARWL1F*G)**0.5d0
            end if
         end if
      elseif ((RBARWE2 > RBARWE1) .and. (RBARWL2 > BARMIN)) then
         !...............WATER LEVEL LOWER ON THIS SIDE OF BARRIER AND IS SUCH
         !................THAT OVERTOPPING IS OCCURING
         !................THUS THIS IS THE RECEIVING SIDE OF THE FLOW ACROSS THE BARRIER
         !................NOTE THAT WE DO FORCE THE RECEIVING SIDE OF THE BARRIER TO
         !................REMAIN WET WHEN THERE IS FLOW ACROSS THE BARRIER.
         !................ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
         !................BARRIER HAS BEEN DRIED. IF IT HAS, WE SHUT DOWN THE FLOW ACROSS
         !................THE BARRIER
         !................ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE BARRIER HAS AT
         !................LEAST ONE WET EDGE. IF NOT, WE SHUT DOWN THE FLOW ACROSS
         !................THE BARRIER. shintaro v46.28.sb05.05 11/01/2006
         if (RBARWL1 > RBARWL2F) then
            !..................INWARD SUBCRITICAL FLOW -> CASE 5
            if (NODECODE(NNBB2) == 0 .or. NNBB2WN == 0) then
               FLUX = 0.0d0
            else
               FLUX = RampIntFlux*RBARWL1*BARINCFSB(I)* &
                      (2.0d0*G*(RBARWE2 - RBARWL1))**0.5d0
               NIBNODECODE(NNBB1) = 1
            end if
         else
            !..................INWARD SUPERCRITICAL FLOW -> CASE 6
            if (NODECODE(NNBB2) == 0 .or. NNBB2WN == 0) then
               FLUX = 0.0d0
            else
               FLUX = RampIntFlux*BARINCFSP(I)*RBARWL2F* &
                      (RBARWL2F*G)**0.5d0
               NIBNODECODE(NNBB1) = 1
            end if
         end if
      else
         FLUX = 0
      end if

      return
      !-----------------------------------------------------------------------
   end subroutine COMPUTE_INTERNAL_BOUNDARY_FLUX
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E  S E T _ S U B M E R G E D 64 _ A T
   !-----------------------------------------------------------------------
   !  THIS SUBROUTINE SETS THE SUBMERGED FLAG
   !-----------------------------------------------------------------------
   subroutine SET_SUBMERGED64_AT( &
      BARRIER_INDEX, BOUNDARYNODE, BOUNDARY, TIMELOC)

      use BOUNDARIES, only: NVELL, IBCONN, BARINHT, ISSUBMERGED64, NBV, LBCODEI, CSII, SIII, BARINCFSB
      use GLOBAL, only: TVW, NODECODE, H0, &
                        IFSFM, IFNLFA, UU1, VV1, QN0, QN1
      use ADC_CONSTANTS, only: G
      use MESH, only: LBArray_Pointer, DP, AREAS, NM, &
                      NEITAB, NNEIGH, &
                      NNEIGHELE, NEITABELE, &
                      SFacEle, SFMYEle, SFMXEle, FDXE, FDYE, X, Y
      implicit none
      integer, intent(IN), target :: BARRIER_INDEX
      integer, intent(IN), target :: BOUNDARY
      integer, intent(IN), target :: BOUNDARYNODE
      real(8), intent(IN) :: TIMELOC

      integer :: NNBB1, NNBB2, NNBBL, NNBBH
      integer :: FLOWDIR
      integer, pointer :: I, J, K
      real(8) :: RBARWL1, RBARWL2, RBARWLL, RBARWLH
      real(8) :: RBARWL1F, RBARWL2F
      integer :: ISSUBMP

      real(8) :: sfdxfac, sfdyfac, sfacavg, SFmxAvg, SFmyAvg
      real(8) :: fdx1, fdx2, fdx3, fdy1, fdy2, fdy3
      real(8) :: etaN1, etaN2, etaN3
      real(8) :: detadxa, detadya, etaslp
      integer :: ineiele, jnei, nm1, nm2, nm3, nctot, wetcnt

      integer :: NNBBNEI, NWETNEI, II, JJ
      logical :: UNSUBMERGE
      real(8) :: X1, X2, Y1, Y2, ET1, ET2, LEN, SLP, HTOT
      real(8) :: QNB, HTOP, FRCREST, FROFF
      real(8) :: UNB, VNB, DXNB, DYNB, QCARRY, ETANEAR, ETAFAR
      integer :: NB, NCNT, IBN, IFOOT, ITOP, OPPFOOT

      LOG_SCOPE_TRACED("SET_SUBMERGED64_AT", WEIR_BOUNDARY_TRACING)

      !...............SIMPLIFY VARIABLES
      I => BARRIER_INDEX
      J => BOUNDARYNODE
      K => BOUNDARY

      NNBB1 = NBV(I) ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
      NNBB2 = IBCONN(I) ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
      FLOWDIR = 0 ! DIRECTION OF FLOW

      !               ! Check if the node ids on both sides of the barrier are valid
      !               ! Either of them can be invalid ocasionally but not always
      !               ! when the barrier is splitted for paralell computation
      if ((NNBB1 <= 0) .or. (NNBB2 <= 0)) then
         ISSUBMERGED64(I) = 0
         if (NNBB2 > 0) then
            if (LBArray_Pointer(NNBB2) > 0) then
               ISSUBMERGED64(LBArray_Pointer(NNBB2)) = 0
            end if
         end if
         return
      end if

      !...............CHECK TO SEE IF THE BARRIER ELEVATION NEEDS UPDATING
      if (INT_TVW) then
         if (BAR_DEG(I)) then
            call COMPUTE_BARRIER_HEIGHT(I, TIMELOC, BARINHT1(I) &
                                        , BARINHT2(I))
         end if
         TVW(NNBB1) = BARINHT2(I) - BARINHT(I)
      end if

      !...............GET WATER LEVEL ABOVE WEIR ON EACH SIDE TO COMPUTE HEAD
      !               DEFINE COMPILER FLAG -DAVERAGEWEIRFLOW TO USE BARAVGWT,
      !               WHICH GENERALLY IS SET TO ZERO, AND IS HARD CODED
      !               TO ZERO HERE. READ_INPUT.F CONTAINS THE SPECIFICATION OF
      !               BARAVGWT. WITH BARAVGWT SET TO ZERO, THERE IS NO NEED FOR
      !               IBSTART EITHER.
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
!...............ZC - STREAMLINE THE PROCESS SINCE BARAVGWT IS ZERO BY DEFAULT
      RBARWL1 = ETA2(NNBB1) - BARINHT2(I)
      RBARWL2 = ETA2(NNBB2) - BARINHT2(I)
#endif
      !...............RBARWL on the lower ground side SB
      if (DP(NNBB1) >= DP(NNBB2)) then
         RBARWLL = RBARWL1
         RBARWLH = RBARWL2
         NNBBL = NNBB1
         NNBBH = NNBB2
      else
         RBARWLL = RBARWL2
         RBARWLH = RBARWL1
         NNBBL = NNBB2
         NNBBH = NNBB1
      end if

      RBARWL1F = 2.d0*RBARWL1/3.d0
      RBARWL2F = 2.d0*RBARWL2/3.d0
      ISSUBMP = ISSUBMERGED64(I)
      ISSUBMERGED64(I) = 0
      ISSUBMERGED64(LBArray_Pointer(NNBB2)) = 0

      !...............Compute the water surface gradient of the higher ground side
      ETASLP = 0.d0
      if (NODECODE(NNBBH) == 1 .and. RBARWLH < BARMIN64_SLIM) then
         WETCNT = 0
         do INEIELE = 1, NNEIGHELE(NNBBH)
            JNEI = NEITABELE(NNBBH, INEIELE)

            NM1 = NM(JNEI, 1)
            NM2 = NM(JNEI, 2)
            NM3 = NM(JNEI, 3)

            NCTOT = NODECODE(NM1) + NODECODE(NM2) + NODECODE(NM3)
            if (NCTOT /= 3) cycle

            WETCNT = WETCNT + 1

            SFacAvg = SFacEle(JNEI)
            SFmxAvg = SFMXEle(JNEI)
            SFmyAvg = SFMYEle(JNEI)
            sfdxfac = dble(1 - IFSFM)*SFacAvg + dble(IFSFM)*SFmxAvg; 
            sfdyfac = dble(1 - IFSFM)*1.0d0 + dble(IFSFM)*SFmyAvg; 
            FDX1 = FDXE(1, JNEI)*sfdxfac !b1
            FDX2 = FDXE(2, JNEI)*sfdxfac !b2
            FDX3 = FDXE(3, JNEI)*sfdxfac !b3
            FDY1 = FDYE(1, JNEI)*sfdyfac !a1
            FDY2 = FDYE(2, JNEI)*sfdyfac !a2
            FDY3 = FDYE(3, JNEI)*sfdyfac !a3

            ETAN1 = ETA2(NM1)
            ETAN2 = ETA2(NM2)
            ETAN3 = ETA2(NM3)

            DETADXA = ETAN1*FDX1 + ETAN2*FDX2 + ETAN3*FDX3
            DETADYA = ETAN1*FDY1 + ETAN2*FDY2 + ETAN3*FDY3
            ETASLP = ETASLP + sqrt(DETADXA**2 + DETADYA**2)/AREAS(JNEI)
         end do
         if (WETCNT /= 0) then
            ETASLP = ETASLP/dble(WETCNT)
         else
            ETASLP = 99999.d0
         end if
      end if

      !               ! Consider submerged if the water depth over the weir is greater than a threshold
      !               ! and the averaged elemental water surface slope on the higher ground side is less
      !               ! than a threshold
      if (((ISSUBMP == 0 .and. RBARWLL >= BARMIN64_NOSUBM) .or. &
           (ISSUBMP == 1 .and. RBARWLL >= BARMIN64_SUBM)) .and. &
          (ETASLP <= BARSLIM64_ELEM)) then
         ISSUBMERGED64(I) = 1
         if (NNBB2 > 0) then
            ISSUBMERGED64(LBArray_Pointer(NNBB2)) = 1
         end if
      end if

      !...............VEW_CREST_HYDRAULICS: crest Froude number from the discharge of the foot's (lower
      !               ground side node's) neighbours that are not on a VEW, so that the decision does
      !               not depend on the pair's own (shared) velocity. Hysteresis: uncouple above
      !               VEW_FROUDE_CRIT; recouple only below 0.8*VEW_FROUDE_CRIT and once the head across
      !               the wall face has come down to the velocity head of that Froude number.
      !               Only while the wall-top node is wet: a dry top (opening gate) is left to the
      !               depth test, so that coupling can wet it.
      if (VEW_CREST_HYDRAULICS .and. ISSUBMERGED64(I) /= 0 .and. NODECODE(NNBBH) /= 0) then
         HTOP = ETA2(NNBBH) + dble(IFNLFA)*DP(NNBBH)
         QNB = 0.d0
         NCNT = 0
         do JJ = 2, NNeigh(NNBBL)
            NB = NeiTab(NNBBL, JJ)
            if (NB == NNBBL .or. NODECODE(NB) == 0) cycle
            IBN = LBArray_Pointer(NB)
            if (IBN > 0) then
               if (LBCODEI(IBN) == 64) cycle
            end if
            QNB = QNB + (UU1(NB)*CSII(I) + VV1(NB)*SIII(I))*(ETA2(NB) + dble(IFNLFA)*DP(NB))
            NCNT = NCNT + 1
         end do
         !...walls with condensed tops: the Froude test acts on the lee pair only, the pair whose
         !   opposite foot's neighbours stand higher; the upstream pair stays coupled
         OPPFOOT = OPPOSITE_FOOT(I)
         if (OPPFOOT > 0) then
            ETANEAR = NEIGHBOR_LEVEL(NNBBL)
            ETAFAR = NEIGHBOR_LEVEL(OPPFOOT)
            if (.not. (ETANEAR > -9000.d0 .and. ETAFAR > -9000.d0 .and. ETAFAR > ETANEAR)) NCNT = 0
         end if
         if (NCNT > 0 .and. HTOP > H0) then
            FRCREST = abs(QNB/dble(NCNT))/HTOP/sqrt(G*HTOP)
            FROFF = VEW_FROUDE_RECOUPLE*VEW_FROUDE_CRIT
            if ((ISSUBMP /= 0 .and. FRCREST > VEW_FROUDE_CRIT) .or. &
                (ISSUBMP == 0 .and. (FRCREST >= FROFF .or. &
                                     abs(ETA2(NNBB1) - ETA2(NNBB2)) >= 0.5d0*FROFF**2*HTOP))) then
               ISSUBMERGED64(I) = 0
               if (NNBB2 > 0) then
                  ISSUBMERGED64(LBArray_Pointer(NNBB2)) = 0
               end if
            end if
         end if
      end if

      !...............Check if the water suraface slope along the weir is greater than a threshold
      !...............If the slope is greater than a threshold, then the weir is not submerged
      HTOT = ETA2(NNBB1) + dble(IFNLFA)*DP(NNBB1)
      if (ISSUBMERGED64(I) /= 0 .and. &
          NODECODE(NNBB1) /= 0 .and. &
          HTOT < BARMIN64_CORNER) then
         UNSUBMERGE = .false.
         X1 = X(NNBB1)
         Y1 = Y(NNBB1)
         ET1 = ETA2(NNBB1)
         if (J > 1) then
            NNBBNEI = NBV(I - 1)
            X2 = X(NNBBNEI)
            Y2 = Y(NNBBNEI)
            ET2 = ETA2(NNBBNEI)
            LEN = sqrt((X1 - X2)**2 + (Y1 - Y2)**2)
            SLP = abs(ET1 - ET2)/LEN
            if (SLP > BARSLIM64_EDGE) then
               UNSUBMERGE = .true.
            end if
         end if
         if (J < NVELL(K)) then
            NNBBNEI = NBV(I + 1)
            X2 = X(NNBBNEI)
            Y2 = Y(NNBBNEI)
            ET2 = ETA2(NNBBNEI)
            LEN = sqrt((X1 - X2)**2 + (Y1 - Y2)**2)
            SLP = abs(ET1 - ET2)/LEN
            if (SLP > BARSLIM64_EDGE) then
               UNSUBMERGE = .true.
            end if
         end if
         if (.not. UNSUBMERGE) then
            NWETNEI = 0
            do JJ = 2, NNeigh(NNBB1)
               II = NeiTab(NNBB1, JJ)
               if (NODECODE(II) /= 0) then
                  NWETNEI = NWETNEI + 1
               end if
            end do
            if (NWETNEI < 3) then
               UNSUBMERGE = .true.
            end if
         end if
         if (UNSUBMERGE) then
            ISSUBMERGED64(I) = 0
            if (NNBB2 > 0) then
               ISSUBMERGED64(LBArray_Pointer(NNBB2)) = 0
            end if
         end if
      end if

      !...............VEW_CREST_HYDRAULICS: the pair has just left coupled mode; carry the discharge into
      !               its flux history and split its levels to the head that gives that discharge
      if (VEW_CREST_HYDRAULICS .and. ISSUBMP /= 0 .and. &
          ISSUBMERGED64(I) == 0 .and. NNBB2 > 0) then
         if (NODECODE(NNBB1) /= 0 .and. NODECODE(NNBB2) /= 0) then
            QNB = 0.d0
            UNB = 0.d0
            VNB = 0.d0
            DXNB = 0.d0
            DYNB = 0.d0
            NCNT = 0
            do JJ = 2, NNeigh(NNBBL)
               NB = NeiTab(NNBBL, JJ)
               if (NB == NNBBL) cycle
               DXNB = DXNB + X(NB) - X(NNBBL)
               DYNB = DYNB + Y(NB) - Y(NNBBL)
               if (NODECODE(NB) == 0) cycle
               IBN = LBArray_Pointer(NB)
               if (IBN > 0) then
                  if (LBCODEI(IBN) == 64) cycle
               end if
               QNB = QNB + (UU1(NB)*CSII(I) + VV1(NB)*SIII(I))*(ETA2(NB) + dble(IFNLFA)*DP(NB))
               UNB = UNB + UU1(NB)
               VNB = VNB + VV1(NB)
               NCNT = NCNT + 1
            end do
            if (NCNT > 0) then
               !...into the domain at the foot when the flow leaves the foot into its bed side
               !   (at most the critical discharge over the crest: the neighbours' flow can be far more
               !   than passes a nearly dry top, e.g. a wave reflecting off a closed gate)
               HTOP = max(ETA2(NNBBH) + dble(IFNLFA)*DP(NNBBH), 0.d0)
               QCARRY = min(abs(QNB/dble(NCNT)), HTOP*sqrt(G*HTOP))
               if (UNB*DXNB + VNB*DYNB < 0.d0) QCARRY = -QCARRY
               IFOOT = merge(I, LBArray_Pointer(NNBB2), NNBBL == NNBB1)
               ITOP = merge(LBArray_Pointer(NNBB2), I, NNBBL == NNBB1)
               QN0(IFOOT) = QCARRY
               QN1(IFOOT) = QCARRY
               QN0(ITOP) = -QCARRY
               QN1(ITOP) = -QCARRY
               if (HTOP > H0 .and. BARINCFSB(I) > 0.d0) then
                  call SPLIT_PAIR_LEVELS(NNBBH, NNBBL, sign((QCARRY/(BARINCFSB(I)*HTOP))**2/(2.d0*G), QCARRY))
               end if
            end if
         end if
      end if

      return
      !-----------------------------------------------------------------------
   end subroutine SET_SUBMERGED64_AT
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !  NEIGHBOR_LEVEL: mean water level of the wet neighbours of node NF that are not on a VEW, or
   !  -99999 if there are none.
   !-----------------------------------------------------------------------
   real(8) function NEIGHBOR_LEVEL(NF)
      use BOUNDARIES, only: LBCODEI
      use GLOBAL, only: NODECODE
      use MESH, only: LBArray_Pointer, NEITAB, NNEIGH
      implicit none
      integer, intent(IN) :: NF
      integer :: JJ, NB, IBN, CNT

      NEIGHBOR_LEVEL = 0.d0
      CNT = 0
      do JJ = 2, NNeigh(NF)
         NB = NeiTab(NF, JJ)
         if (NB == NF .or. NODECODE(NB) == 0) cycle
         IBN = LBArray_Pointer(NB)
         if (IBN > 0) then
            if (LBCODEI(IBN) == 64) cycle
         end if
         NEIGHBOR_LEVEL = NEIGHBOR_LEVEL + ETA2(NB)
         CNT = CNT + 1
      end do
      if (CNT > 0) then
         NEIGHBOR_LEVEL = NEIGHBOR_LEVEL/dble(CNT)
      else
         NEIGHBOR_LEVEL = -99999.d0
      end if
   end function NEIGHBOR_LEVEL

   !-----------------------------------------------------------------------
   !  OPPOSITE_FOOT: the foot node on the other face of the wall from boundary entry I, through the
   !  condensed wall tops; 0 if the entry's top is not in a two-node group.
   !-----------------------------------------------------------------------
   integer function OPPOSITE_FOOT(I)
      use BOUNDARIES, only: NVEL, NBV, IBCONN
      use MESH, only: LBArray_Pointer
      use NodalAttributes, only: LoadCondensedNodes, NListCondensedNodes, ListCondensedNodes, &
                                 NNodesListCondensedNodes
      implicit none
      integer, intent(IN) :: I
      integer, allocatable, save :: OPP(:)
      integer :: K, L, M, T1, T2, E1, E2, F1, F2

      if (.not. allocated(OPP)) then
         allocate (OPP(NVEL))
         OPP = 0
         if (LoadCondensedNodes) then
            do K = 1, NListCondensedNodes
               T1 = 0
               T2 = 0
               M = 0
               do L = 1, NNodesListCondensedNodes(K)
                  if (ListCondensedNodes(K, L) <= 0) cycle
                  M = M + 1
                  if (M == 1) T1 = ListCondensedNodes(K, L)
                  if (M == 2) T2 = ListCondensedNodes(K, L)
               end do
               if (M /= 2) cycle
               E1 = LBArray_Pointer(T1)
               E2 = LBArray_Pointer(T2)
               if (E1 <= 0 .or. E2 <= 0) cycle
               F1 = merge(IBCONN(E1), NBV(E1), NBV(E1) == T1)
               F2 = merge(IBCONN(E2), NBV(E2), NBV(E2) == T2)
               OPP(E1) = F2
               OPP(LBArray_Pointer(F1)) = F2
               OPP(E2) = F1
               OPP(LBArray_Pointer(F2)) = F1
            end do
         end if
      end if
      OPPOSITE_FOOT = OPP(I)
   end function OPPOSITE_FOOT

   !-----------------------------------------------------------------------
   !  SPLIT_PAIR_LEVELS (VEW_CREST_HYDRAULICS): raise the wall-top side and lower the foot side of a
   !  pair that has just left coupled mode so that their levels differ by HEAD (top minus foot),
   !  conserving volume. Each side moves as its coupled group: VEW-coupled pairs and condensed
   !  nodes (wet members), at both time levels of the water level.
   !-----------------------------------------------------------------------
   subroutine SPLIT_PAIR_LEVELS(NTOP, NFOOT, HEAD)
      use BOUNDARIES, only: NBOU, NVELL, NBV, IBCONN, LBCODEI, ISSUBMERGED64
      use GLOBAL, only: NODECODE, ETA1, ETA2, H2, IFNLFA
      use MESH, only: NP, DP, TotalArea
      use NodalAttributes, only: LoadCondensedNodes, NListCondensedNodes, ListCondensedNodes, &
                                 NNodesListCondensedNodes
      implicit none
      integer, intent(IN) :: NTOP, NFOOT
      real(8), intent(IN) :: HEAD
      integer, allocatable :: root(:)
      integer :: I, I0, J, K, L, NV, NA, RT, RF
      real(8) :: AT, AF, DT1, DF1

      allocate (root(NP))
      root = [(J, J=1, NP)]
      I0 = 0
      do K = 1, NBOU
         NV = NVELL(K)
         select case (LBCODEI(I0 + 1))
         case (64)
            do J = 1, NV
               I = I0 + J
               if (ISSUBMERGED64(I) == 0) cycle
               if (NBV(I) <= 0 .or. IBCONN(I) <= 0) cycle
               if (NODECODE(NBV(I)) /= 0 .and. NODECODE(IBCONN(I)) /= 0) call join(NBV(I), IBCONN(I))
            end do
            I0 = I0 + 2*NV
         case (4, 24, 5, 25)
            I0 = I0 + 2*NV
         case default
            I0 = I0 + NV
         end select
      end do
      if (LoadCondensedNodes) then
         do K = 1, NListCondensedNodes
            NA = 0
            do L = 1, NNodesListCondensedNodes(K)
               J = ListCondensedNodes(K, L)
               if (J <= 0) cycle
               if (NODECODE(J) == 0) cycle
               if (NA == 0) then
                  NA = J
               else
                  call join(NA, J)
               end if
            end do
         end do
      end if
      RT = find(NTOP)
      RF = find(NFOOT)
      if (RT == RF) return ! still tied through another pair: one level, nothing to split
      AT = 0.d0
      AF = 0.d0
      do J = 1, NP
         if (NODECODE(J) == 0) cycle
         if (find(J) == RT) AT = AT + TotalArea(J)
         if (find(J) == RF) AF = AF + TotalArea(J)
      end do
      if (AT <= 0.d0 .or. AF <= 0.d0) return
      DT1 = HEAD*AF/(AT + AF)
      DF1 = -HEAD*AT/(AT + AF)
      do J = 1, NP
         if (NODECODE(J) == 0) cycle
         if (find(J) == RT) then
            ETA1(J) = ETA1(J) + DT1
            ETA2(J) = ETA2(J) + DT1
            H2(J) = DP(J) + dble(IFNLFA)*ETA2(J)
         elseif (find(J) == RF) then
            ETA1(J) = ETA1(J) + DF1
            ETA2(J) = ETA2(J) + DF1
            H2(J) = DP(J) + dble(IFNLFA)*ETA2(J)
         end if
      end do

   contains

      recursive integer function find(N) result(R)
         integer, intent(IN) :: N
         if (root(N) == N) then
            R = N
         else
            root(N) = find(root(N))
            R = root(N)
         end if
      end function find

      subroutine join(A, B)
         integer, intent(IN) :: A, B
         integer :: RA, RB
         RA = find(A)
         RB = find(B)
         if (RA /= RB) root(max(RA, RB)) = min(RA, RB)
      end subroutine join
   end subroutine SPLIT_PAIR_LEVELS
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E
   !       C O M P U T E _ I N T E R N A L _ B O U N D A R Y 6 4 _ F L U X
   !-----------------------------------------------------------------------
   !     THIS SUBROUTINE CONTAINS THE CALCULATIONS OF OVERTOPPING FLUX
   !     PASSED BETWEEN THE THE WEIR PAIR NODES. THIS IS THE VERSION OF
   !     THE CODE THAT CHECKS IF THE WEIR HAS A WET EDGE BEFORE PASSING
   !     FLUX OVER THE WEIR.
   !     THIS IS A SUBROUTINE THAT WORKS WITH IBTYPE=64. THIS RETURNS ZERO
   !     WHEN A WEIR IS SUBMERGED ON BOTH SIDES. SB
   !-----------------------------------------------------------------------
   subroutine COMPUTE_INTERNAL_BOUNDARY64_FLUX( &
      BARRIER_INDEX, BOUNDARYNODE, BOUNDARY, TIMELOC, FLUX)

      use BOUNDARIES, only: BARINCFSB, BARINCFSP, &
                            NVELL, NBV, IBCONN, BARINHT, ISSUBMERGED64, BNDEDGE3RD, &
                            CSII, SIII
      use GLOBAL, only: TVW, NODECODE, UU1, VV1
      use MESH, only: LBArray_Pointer, DP
      implicit none
      integer, intent(IN), target :: BARRIER_INDEX
      integer, intent(IN), target :: BOUNDARY
      integer, intent(IN), target :: BOUNDARYNODE
      real(8), intent(IN) :: TIMELOC
      real(8), intent(OUT) :: FLUX

      integer :: NNBB1, NNBB2
      integer :: NNBB1WN, NNBB2WN
      integer :: FLOWDIR
      integer, pointer :: I, J, K
      integer :: I1, I2
      integer :: N11, N12, N21, N22
      integer :: FORCEDWET
      real(8) :: RBARWL1, RBARWL2, RBARWLL
      real(8) :: RBARWL1F, RBARWL2F
      real(8) :: RBARWE1, RBARWE2 !...Heads used for the flow direction and in the weir formula
      real(8) :: VN1, VN2
      logical :: CHECK

      LOG_SCOPE_TRACED("COMPUTE_INTERNAL_BOUNDARY64_FLUX", WEIR_BOUNDARY_TRACING)

      !...............SIMPLIFY VARIABLES
      I => BARRIER_INDEX
      J => BOUNDARYNODE
      K => BOUNDARY

      FORCEDWET = 0

      I2 = LBArray_Pointer(IBCONN(I)) ! I2 is the pointer to the node on the opposite side of the barrier

      NNBB1 = NBV(I) ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
      NNBB2 = IBCONN(I) ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
      NNBB1WN = 0 ! COUNT NUMBER OF WET NEIGHBORS
      NNBB2WN = 0 ! COUNT NUMBER OF WET NEIGHBORS
      FLOWDIR = 0 ! DIRECTION OF FLOW

      !               ! Check if the node ids on both sides of the barrier are valid
      !               ! Either of them can be invalid ocasionally but not always
      !               ! when the barrier is splitted for paralell computation
      if ((NNBB1 <= 0) .or. (NNBB2 <= 0) .or. (I2 <= 0)) then
         FLUX = 0.d0
         return
      end if

      !...............CHECK TO SEE IF THERE IS A WET EDGE
      !                 This check is needed ocasionally but not always
      !                 when the barrier is splitted for paralell computation
      CHECK = .false.
      !             IF( (J.EQ.1).OR.(J.EQ.NVELL(K)+1) )THEN
      !                 IF((NBV(I+1)>0).AND.(IBCONN(I+1)>0)) THEN
      !                     NNBB1WN = NNBB1WN + NODECODE(NBV(I+1))
      !                     NNBB2WN = NNBB2WN + NODECODE(IBCONN(I+1))
      !                     CHECK = .TRUE.
      !                 ENDIF
      !             ELSE
      !                 IF( J.EQ.NVELL(K).OR.(J.EQ.NVELL(K)*2) ) THEN
      !                     IF((NBV(I-1)>0).AND.(IBCONN(I-1)>0)) THEN
      !                         NNBB1WN = NNBB1WN + NODECODE(NBV(I-1))
      !                         NNBB2WN = NNBB2WN + NODECODE(IBCONN(I-1))
      !                         CHECK = .TRUE.
      !                     ENDIF
      !                 ELSE
      !                     IF((NBV(I-1)>0).AND.(IBCONN(I-1)>0).AND.
      !  &                     (NBV(I+1)>0).AND.(IBCONN(I+1)>0)) THEN
      !                         NNBB1WN = NNBB1WN + NODECODE(NBV(I-1))
      !                         NNBB1WN = NNBB1WN + NODECODE(NBV(I+1))
      !                         NNBB2WN = NNBB2WN + NODECODE(IBCONN(I+1))
      !                         NNBB2WN = NNBB2WN + NODECODE(IBCONN(I-1))
      !                         CHECK = .TRUE.
      !                     ENDIF
      !                 ENDIF
      !             ENDIF

      ! ------
      if ((J == 1) .or. (J == NVELL(K) + 1)) then
         if ((NBV(I + 1) > 0) .and. (IBCONN(I + 1) > 0)) then
            I1 = I + 1
            N11 = NBV(I1)
            N12 = BNDEDGE3RD(I1 - 1)
            N21 = NBV(I2 - 1)
            N22 = BNDEDGE3RD(I2 - 1)
            if (N12 /= 0 .and. N22 /= 0) then
               if (NODECODE(N11) /= 0 .and. NODECODE(N12) /= 0 .and. &
                   NODECODE(N21) /= 0 .and. NODECODE(N22) /= 0) then
                  NNBB1WN = NNBB1WN + 1
                  NNBB2WN = NNBB2WN + 1
                  CHECK = .true.
               end if
            end if
         end if
      else
         if (J == NVELL(K) .or. (J == NVELL(K)*2)) then
            if ((NBV(I - 1) > 0) .and. (IBCONN(I - 1) > 0)) then
               I1 = I - 1
               N11 = NBV(I1)
               N12 = BNDEDGE3RD(I1)
               N21 = NBV(I2 + 1)
               N22 = BNDEDGE3RD(I2)
               if (N12 /= 0 .and. N22 /= 0) then
                  if (NODECODE(N11) /= 0 .and. NODECODE(N12) /= 0 .and. &
                      NODECODE(N21) /= 0 .and. NODECODE(N22) /= 0) then
                     NNBB1WN = NNBB1WN + 1
                     NNBB2WN = NNBB2WN + 1
                     CHECK = .true.
                  end if
               end if
            end if
         else
            if ((NBV(I - 1) > 0) .and. (IBCONN(I - 1) > 0) .and. &
                (NBV(I + 1) > 0) .and. (IBCONN(I + 1) > 0)) then
               I1 = I - 1
               N11 = NBV(I1)
               N12 = BNDEDGE3RD(I1)
               N21 = NBV(I2 + 1)
               N22 = BNDEDGE3RD(I2)
               if (N12 /= 0 .and. N22 /= 0) then
                  if (NODECODE(N11) /= 0 .and. NODECODE(N12) /= 0 .and. &
                      NODECODE(N21) /= 0 .and. NODECODE(N22) /= 0) then
                     NNBB1WN = NNBB1WN + 1
                     NNBB2WN = NNBB2WN + 1
                     CHECK = .true.
                  end if
               end if
               I1 = I + 1
               N11 = NBV(I1)
               N12 = BNDEDGE3RD(I1 - 1)
               N21 = NBV(I2 - 1)
               N22 = BNDEDGE3RD(I2 - 1)
               if (N12 /= 0 .and. N22 /= 0) then
                  if (NODECODE(N11) /= 0 .and. NODECODE(N12) /= 0 .and. &
                      NODECODE(N21) /= 0 .and. NODECODE(N22) /= 0) then
                     NNBB1WN = NNBB1WN + 1
                     NNBB2WN = NNBB2WN + 1
                     CHECK = .true.
                  end if
               end if
            end if
         end if
      end if

      if (.not. CHECK) then
         FLUX = 0.d0
         return
      end if

      !..............CHECK TO SEE IF THE BARRIER ELEVATION NEEDS UPDATING
      if (INT_TVW) then
         if (BAR_DEG(I)) then
            call COMPUTE_BARRIER_HEIGHT(I, TIMELOC, BARINHT1(I) &
                                        , BARINHT2(I))
         end if
         TVW(NNBB1) = BARINHT2(I) - BARINHT(I)
      end if

      !..............GET WATER LEVEL ABOVE WEIR ON EACH SIDE TO COMPUTE HEAD
      !              DEFINE COMPILER FLAG -DAVERAGEWEIRFLOW TO USE BARAVGWT,
      !              WHICH GENERALLY IS SET TO ZERO, AND IS HARD CODED
      !              TO ZERO HERE. READ_INPUT.F CONTAINS THE SPECIFICATION OF
      !              BARAVGWT. WITH BARAVGWT SET TO ZERO, THERE IS NO NEED FOR
      !              IBSTART EITHER.
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
!..............ZC - STREAMLINE THE PROCESS SINCE BARAVGWT IS ZERO BY DEFAULT
      RBARWL1 = ETA2(NNBB1) - BARINHT2(I)
      RBARWL2 = ETA2(NNBB2) - BARINHT2(I)
#endif
      !..............RBARWL on the lower ground side SB
      if (DP(NNBB1) >= DP(NNBB2)) then
         RBARWLL = RBARWL1
      else
         RBARWLL = RBARWL2
      end if

      !..............Total energy head (VEW_CREST_HYDRAULICS): add the velocity head of the raw nodal
      !              velocity normal to the boundary, vn = u*cos + v*sin of the boundary normal.
      !              The static heads RBARWL1/RBARWL2 are kept for the checks against BARMIN64 and
      !              for the downstream depth in the subcritical formula.
      RBARWE1 = RBARWL1
      RBARWE2 = RBARWL2
      if (VEW_CREST_HYDRAULICS) then
         VN1 = UU1(NNBB1)*CSII(I) + VV1(NNBB1)*SIII(I)
         VN2 = UU1(NNBB2)*CSII(I2) + VV1(NNBB2)*SIII(I2)
         RBARWE1 = RBARWL1 + VN1*VN1/(2.d0*G)
         RBARWE2 = RBARWL2 + VN2*VN2/(2.d0*G)
      end if

      RBARWL1F = 2.d0*RBARWE1/3.d0
      RBARWL2F = 2.d0*RBARWE2/3.d0
      FLUX = 0.d0

      if (ISSUBMERGED64(I) /= 0) then
         !...............THIS IS THE MAJOR DIFFERENCE OF IBTYPE=64 FROM IBTYPE=24.
         !               WATER LEVEL ON BOTH SIDES OF BARRIER ABOVE BARRIER -> CASE 0
         FLUX = 0.d0
      elseif ((RBARWL1 < BARMIN64) .and. (RBARWL2 < BARMIN64)) then
         !...............WATER LEVEL ON BOTH SIDES OF BARRIER BELOW BARRIER -> CASE 1
         FLUX = 0.d0
      elseif (abs(RBARWE1 - RBARWE2) < 0.01d0) then
         !...............WATER LEVEL EQUAL ON BOTH SIDES OF BARRIER
         !................TO WITHIN TOLERANCE BARMIN -> CASE 2
         FLUX = 0.d0
      elseif ((RBARWE1 > RBARWE2) .and. (RBARWL1 > BARMIN64)) then
         !...............WATER LEVEL GREATER ON THIS SIDE OF THE BARRIER AND IS SUCH
         !................THAT OVERTOPPING IS OCCURING
         !................THUS THIS SIDE IS THE SOURCE SIDE OF THE FLOW ACROSS THE BARRIER
         !................NOTE THAT WE DO NOT FORCE THE SOURCE SIDE OF THE BARRIER TO
         !................REMAIN WET. ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
         !................BARRIER HAS BEEN DRIED. IF IT HAS WE SHUT DOWN THE FLOW ACROSS
         !................THE BARRIER
         !................ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE BARRIER HAS AT
         !................LEAST ONE WET EDGE. IF NOT, WE SHUT DOWN THE FLOW ACROSS
         !................THE BARRIER. shintaro v46.28.sb05.05 11/01/2006
         if (RBARWL2 > RBARWL1F) then
            !..................OUTWARD SUBCRITICAL FLOW -> CASE 3
            if (NODECODE(NNBB1) == 0 .or. NNBB1WN == 0 .or. &
                FORCEDWET == 1) then
               FLUX = 0.0d0
            else
               FLUX = -RampIntFlux*RBARWL2*BARINCFSB(I)* &
                      (2.d0*G*(RBARWE1 - RBARWL2))**0.5d0
            end if
         else
            !..................OUTWARD SUPERCRITICAL FLOW -> CASE 4
            if (NODECODE(NNBB1) == 0 .or. NNBB1WN == 0 .or. &
                FORCEDWET == 1) then
               FLUX = 0.0d0
            else
               FLUX = -RampIntFlux*BARINCFSP(I)*RBARWL1F* &
                      (RBARWL1F*G)**0.5d0
            end if
         end if
      elseif ((RBARWE2 > RBARWE1) .and. (RBARWL2 > BARMIN64)) then
         !...............WATER LEVEL LOWER ON THIS SIDE OF BARRIER AND IS SUCH
         !................THAT OVERTOPPING IS OCCURING
         !................THUS THIS IS THE RECEIVING SIDE OF THE FLOW ACROSS THE BARRIER
         !................NOTE THAT WE DO FORCE THE RECEIVING SIDE OF THE BARRIER TO
         !................REMAIN WET WHEN THERE IS FLOW ACROSS THE BARRIER.
         !................ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
         !................BARRIER HAS BEEN DRIED. IF IT HAS, WE SHUT DOWN THE FLOW ACROSS
         !................THE BARRIER
         !................ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE BARRIER HAS AT
         !................LEAST ONE WET EDGE. IF NOT, WE SHUT DOWN THE FLOW ACROSS
         !................THE BARRIER. shintaro v46.28.sb05.05 11/01/2006
         if (RBARWL1 > RBARWL2F) then
            !..................INWARD SUBCRITICAL FLOW -> CASE 5
            if (NODECODE(NNBB2) == 0 .or. NNBB2WN == 0 .or. &
                FORCEDWET == 1) then
               FLUX = 0.0d0
            else
               FLUX = RampIntFlux*RBARWL1*BARINCFSB(I)* &
                      (2.0d0*G*(RBARWE2 - RBARWL1))**0.5d0
               !                       ! WRITE(*,*) "Receiving flux, subcritical"
               !                       ! NIBNODECODE(NNBB1)=1
            end if
         else
            !..................INWARD SUPERCRITICAL FLOW -> CASE 6
            if (NODECODE(NNBB2) == 0 .or. NNBB2WN == 0 .or. &
                FORCEDWET == 1) then
               FLUX = 0.0d0
            else
               FLUX = RampIntFlux*BARINCFSP(I)*RBARWL2F* &
                      (RBARWL2F*G)**0.5d0
               !                       ! WRITE(*,*) "Receiving flux, supercritical"
               !                       ! NIBNODECODE(NNBB1)=1
            end if
         end if
      else
         FLUX = 0.d0
      end if

      return
      !-----------------------------------------------------------------------
   end subroutine COMPUTE_INTERNAL_BOUNDARY64_FLUX

   !-----------------------------------------------------------------------
   !  SEMI_WEIR_PREPARE (WEIR_SEMI_IMPLICIT): mark the boundary entries of the pairs solved at the
   !  new time level this step (IBTYPE 4/24 pairs and IBTYPE 64 pairs in weir mode, both nodes wet),
   !  and on the first call set each entry's lumped boundary length (3/4 of BndLen2O3 of each adjacent
   !  edge on the same side of the barrier) and its partner entry. Called before the GWCE load vector
   !  takes the flux-boundary forcing; the active entries are skipped there.
   !-----------------------------------------------------------------------
   subroutine SEMI_WEIR_PREPARE()
      use BOUNDARIES, only: NBOU, NVELL, NVEL, NBV, IBCONN, LBCODEI, ISSUBMERGED64, BndLen2O3, &
                            NFLUXIB64_GBL
      use GLOBAL, only: NODECODE, ILUMP, CPRECOR
      use SIZES, only: MNPROC
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none
      integer :: I0, K, J, NV, IA, IB

      SW_ANY = .false.
      if (.not. WEIR_SEMI_IMPLICIT) return
      if (.not. allocated(SW_ACT)) then
         if (ILUMP == 0) call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
            message="WEIR_SEMI_IMPLICIT needs the lumped GWCE (IM digit 6 = 2)")
         !...in parallel runs for IBTYPE 4/24 only for now: an IBTYPE=64 pair's coupled group (condensed
         !   wall top and coupled feet) can be split between subdomains
         if (MNPROC > 1 .and. NFLUXIB64_GBL > 0) call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
            message="WEIR_SEMI_IMPLICIT is serial only for now on meshes with IBTYPE=64 boundaries")
         if (CPRECOR) call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
            message="WEIR_SEMI_IMPLICIT does not support the predictor-corrector GWCE yet")
         allocate (SW_ACT(NVEL), SW_LEN(NVEL), SW_PARTNER(NVEL))
         SW_LEN = 0.d0
         SW_PARTNER = 0
         I0 = 0
         do K = 1, NBOU
            NV = NVELL(K)
            select case (LBCODEI(I0 + 1))
            case (4, 24, 64)
               do J = 1, NV
                  SW_PARTNER(I0 + J) = I0 + NV + J
                  SW_PARTNER(I0 + NV + J) = I0 + J
               end do
               do J = 2, NV ! edges within each side
                  SW_LEN(I0 + J - 1) = SW_LEN(I0 + J - 1) + 0.75d0*BndLen2O3(I0 + J - 1)
                  SW_LEN(I0 + J) = SW_LEN(I0 + J) + 0.75d0*BndLen2O3(I0 + J - 1)
                  SW_LEN(I0 + NV + J - 1) = SW_LEN(I0 + NV + J - 1) + 0.75d0*BndLen2O3(I0 + NV + J - 1)
                  SW_LEN(I0 + NV + J) = SW_LEN(I0 + NV + J) + 0.75d0*BndLen2O3(I0 + NV + J - 1)
               end do
               I0 = I0 + 2*NV
            case (5, 25)
               I0 = I0 + 2*NV
            case default
               I0 = I0 + NV
            end select
         end do
      end if
      SW_ACT = .false.
      I0 = 0
      do K = 1, NBOU
         NV = NVELL(K)
         select case (LBCODEI(I0 + 1))
         case (4, 24, 64)
            do J = 1, NV
               IA = I0 + J
               IB = I0 + NV + J
               if (NBV(IA) <= 0 .or. IBCONN(IA) <= 0) cycle
               if (NODECODE(NBV(IA)) == 0 .or. NODECODE(IBCONN(IA)) == 0) cycle
               if (LBCODEI(IA) == 64) then
                  if (ISSUBMERGED64(IA) /= 0 .or. ISSUBMERGED64(IB) /= 0) cycle
               end if
               SW_ACT(IA) = .true.
               SW_ACT(IB) = .true.
               SW_ANY = .true.
            end do
            I0 = I0 + 2*NV
         case (5, 25)
            I0 = I0 + 2*NV
         case default
            I0 = I0 + NV
         end select
      end do
   end subroutine SEMI_WEIR_PREPARE

   !-----------------------------------------------------------------------
   !  SEMI_WEIR_SOLVE (WEIR_SEMI_IMPLICIT): after the lumped GWCE division without the active pairs'
   !  forcing (ETAS = increment, COEFDT = diagonal after the VEW / condensed-node summation), solve
   !  each active pair's flux with both levels at the new time level and correct the increments of
   !  the coupled groups of its two nodes. The pair forcing is f = c+ q(n+1) - c- q(n-1),
   !  c+- = 1/(2 dt) +- tau0/2, through one length per pair; the root of
   !  G(f) = f + c- q(n-1) - c+ q(eta_A* - a_A f, eta_B* + a_B f) is found by false position (Illinois).
   !  The solved flux is written to QN2 (positive into the domain) for the momentum equations.
   !-----------------------------------------------------------------------
   subroutine SEMI_WEIR_SOLVE(ETAS, COEFDT, ETA1, TAU0V, DTS, QN0, QN2)
      use BOUNDARIES, only: NBOU, NVELL, NVEL, NBV, IBCONN, LBCODEI, ISSUBMERGED64, &
                            BARINCFSB, BARINCFSP, CSII, SIII
      use GLOBAL, only: NODECODE, UU1, VV1, RampIntFlux
      use MESH, only: NP
      use NodalAttributes, only: LoadCondensedNodes, NListCondensedNodes, ListCondensedNodes, &
                                 NNodesListCondensedNodes
      implicit none
      real(8), intent(INOUT) :: ETAS(:)
      real(8), intent(IN) :: COEFDT(:), ETA1(:), TAU0V(:), QN0(:)
      real(8), intent(IN) :: DTS
      real(8), intent(INOUT) :: QN2(:)
      integer, allocatable :: root(:), kcnt(:)
      real(8), allocatable :: deta(:)
      integer :: I, I2, J, K, L, NA, NB, RA, RB, I0, NV, IT, SIDE
      real(8) :: ell, aa, ab, ea, eb, qold, cp, cm, tau, zc, csb, csp, hva, hvb, va, vb
      real(8) :: fl, fa, fb, ga, gb, fm, gm, qnew, s

      if (.not. SW_ANY) return
      allocate (root(NP), kcnt(NP), deta(NP))
      root = [(J, J=1, NP)]
      kcnt = 0
      deta = 0.d0
      !...coupled groups: VEW-coupled pairs and condensed nodes (wet members), as in the GWCE summation
      I0 = 0
      do K = 1, NBOU
         NV = NVELL(K)
         select case (LBCODEI(I0 + 1))
         case (64)
            do J = 1, NV
               I = I0 + J
               if (ISSUBMERGED64(I) == 0) cycle
               if (NBV(I) <= 0 .or. IBCONN(I) <= 0) cycle
               if (NODECODE(NBV(I)) /= 0 .and. NODECODE(IBCONN(I)) /= 0) call join(NBV(I), IBCONN(I))
            end do
            I0 = I0 + 2*NV
         case (4, 24, 5, 25)
            I0 = I0 + 2*NV
         case default
            I0 = I0 + NV
         end select
      end do
      if (LoadCondensedNodes) then
         do K = 1, NListCondensedNodes
            NA = 0
            do L = 1, NNodesListCondensedNodes(K)
               J = ListCondensedNodes(K, L)
               if (J <= 0) cycle
               if (NODECODE(J) == 0) cycle
               if (NA == 0) then
                  NA = J
               else
                  call join(NA, J)
               end if
            end do
         end do
      end if
      !...number of active pairs on each group
      do I = 1, NVEL
         if (.not. SW_ACT(I)) cycle
         I2 = SW_PARTNER(I)
         if (I2 < I) cycle
         RA = find(NBV(I))
         RB = find(IBCONN(I))
         kcnt(RA) = kcnt(RA) + 1
         kcnt(RB) = kcnt(RB) + 1
      end do
      !...pair solves
      do I = 1, NVEL
         if (.not. SW_ACT(I)) cycle
         I2 = SW_PARTNER(I)
         if (I2 < I) cycle
         NA = NBV(I)
         NB = IBCONN(I)
         RA = find(NA)
         RB = find(NB)
         if (COEFDT(NA) <= 0.d0 .or. COEFDT(NB) <= 0.d0) cycle
         ell = 0.5d0*(SW_LEN(I) + SW_LEN(I2))
         if (ell <= 0.d0) cycle
         aa = ell*dble(kcnt(RA))/COEFDT(NA)
         ab = ell*dble(kcnt(RB))/COEFDT(NB)
         ea = ETA1(NA) + ETAS(NA)
         eb = ETA1(NB) + ETAS(NB)
         zc = BARINHT2(I)
         csb = BARINCFSB(I)
         csp = BARINCFSP(I)
         hva = 0.d0
         hvb = 0.d0
         if ((LBCODEI(I) == 64 .and. VEW_CREST_HYDRAULICS) .or. &
             (LBCODEI(I) /= 64 .and. WEIR_ENERGY_HEAD)) then
            va = UU1(NA)*CSII(I) + VV1(NA)*SIII(I)
            vb = UU1(NB)*CSII(I2) + VV1(NB)*SIII(I2)
            hva = va*va/(2.d0*G)
            hvb = vb*vb/(2.d0*G)
         end if
         qold = 0.5d0*(QN0(I2) - QN0(I)) ! A -> B, from the into-the-domain fluxes of both sides
         tau = 0.5d0*(TAU0V(NA) + TAU0V(NB))
         cp = 0.5d0/DTS + 0.5d0*tau
         cm = 0.5d0/DTS - 0.5d0*tau
         !...bracket the root of the increasing G around the level-equalising f
         fl = (ea - eb)/(aa + ab)
         gm = gfun(fl)
         s = abs(fl) + abs(cm*qold) + 1.d-12
         if (gm > 0.d0) then
            fb = fl
            gb = gm
            fa = min(fl, -cm*qold) - s
            ga = gfun(fa)
            do IT = 1, 60
               if (ga <= 0.d0) exit
               s = 2.d0*s
               fa = fa - s
               ga = gfun(fa)
            end do
         else
            fa = fl
            ga = gm
            fb = max(fl, -cm*qold) + s
            gb = gfun(fb)
            do IT = 1, 60
               if (gb >= 0.d0) exit
               s = 2.d0*s
               fb = fb + s
               gb = gfun(fb)
            end do
         end if
         !...false position with the Illinois modification
         SIDE = 0
         fm = fa
         do IT = 1, 100
            if (gb == ga) exit
            fm = (fa*gb - fb*ga)/(gb - ga)
            gm = gfun(fm)
            if (abs(gm) <= 1.d-14*(abs(fm) + abs(cm*qold) + 1.d-30) .or. abs(fb - fa) <= 1.d-14*(abs(fm) + 1.d-30)) exit
            if (gm > 0.d0) then
               fb = fm
               gb = gm
               if (SIDE == 1) ga = 0.5d0*ga
               SIDE = 1
            else
               fa = fm
               ga = gm
               if (SIDE == -1) gb = 0.5d0*gb
               SIDE = -1
            end if
         end do
         qnew = (fm + cm*qold)/cp
         deta(RA) = deta(RA) - ell*fm/COEFDT(NA)
         deta(RB) = deta(RB) + ell*fm/COEFDT(NB)
         QN2(I) = -qnew
         QN2(I2) = qnew
      end do
      !...correct the increments of every wet member of each group
      do J = 1, NP
         if (NODECODE(J) == 0) cycle
         RA = find(J)
         if (deta(RA) /= 0.d0) ETAS(J) = ETAS(J) + deta(RA)
      end do
      deallocate (root, kcnt, deta)

   contains

      real(8) function gfun(f)
         real(8), intent(IN) :: f
         gfun = f + cm*qold - cp*qlaw(ea - aa*f - zc, eb + ab*f - zc)
      end function gfun

      !...signed weir discharge A -> B for heads ha, hb over the crest (ADCIRC IBTYPE 24/64 law)
      real(8) function qlaw(ha, hb)
         real(8), intent(IN) :: ha, hb
         real(8) :: eu, hd, hu
         if (ha + hva >= hb + hvb) then
            hu = ha
            eu = ha + hva
            hd = hb
            qlaw = 1.d0
         else
            hu = hb
            eu = hb + hvb
            hd = ha
            qlaw = -1.d0
         end if
         if (hu <= 0.d0 .or. eu <= 0.d0) then
            qlaw = 0.d0
         elseif (hd > 2.d0*eu/3.d0) then
            qlaw = qlaw*RampIntFlux*csb*hd*sqrt(2.d0*G*max(eu - hd, 0.d0))
         else
            qlaw = qlaw*RampIntFlux*csp*(2.d0*eu/3.d0)*sqrt(2.d0*eu*G/3.d0)
         end if
      end function qlaw

      recursive integer function find(i) result(ri)
         integer, intent(IN) :: i
         if (root(i) == i) then
            ri = i
         else
            root(i) = find(root(i))
            ri = root(i)
         end if
      end function find

      subroutine join(i, j)
         integer, intent(IN) :: i, j
         integer :: ri, rj
         ri = find(i)
         rj = find(j)
         if (ri /= rj) root(max(ri, rj)) = min(ri, rj)
      end subroutine join

   end subroutine SEMI_WEIR_SOLVE

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E
   !       C O M P U T E _ C R O S S _ B A R R I E R _ P I P E _ F L U X
   !-----------------------------------------------------------------------
   !  THIS ROUTINE COMPUTE THE DISCHARGE ACROSS TYPE 5,25 BOUNDARY
   !  CONDITIONS (CROSS BARRIER PIPES)
   !-----------------------------------------------------------------------
   subroutine COMPUTE_CROSS_BARRIER_PIPE_FLUX(IDX, FLUX)
      use BOUNDARIES, only: NBV, IBCONN, PIPEHT, PIPEDIAM, PIPECOEF
      use ADC_CONSTANTS, only: G, PI
      use GLOBAL, only: NODECODE, NIBNODECODE, ETA2, &
                        RampIntFlux

      implicit none
      integer, intent(IN) :: IDX
      real(8), intent(OUT) :: FLUX

      integer :: NNBB1
      integer :: NNBB2

      real(8) :: RBARWL1
      real(8) :: RBARWL2

      LOG_SCOPE_TRACED("COMPUTE_CROSS_BARRIER_PIPE_FLUX", WEIR_BOUNDARY_TRACING)

      NNBB1 = NBV(IDX) ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
      NNBB2 = IBCONN(IDX) ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER

#ifdef AVERAGEWEIRFLOW
      if (IBSTART == 0) then
         RPIPEWL1AVG(IDX) = ETA2(NNBB1) - PIPEHT(IDX)
         RPIPEWL2AVG(IDX) = ETA2(NNBB2) - PIPEHT(IDX)
         IBSTART = 1
      else
         RPIPEWL1AVG(IDX) = (ETA2(NNBB1) - PIPEHT(IDX) + BARAVGWT &
                             *RPIPEWL1AVG(IDX))/(1.d0 + BARAVGWT)
         RPIPEWL2AVG(IDX) = (ETA2(NNBB2) - PIPEHT(IDX) + BARAVGWT &
                             *RPIPEWL2AVG(IDX))/(1.d0 + BARAVGWT)
      end if
      RBARWL1 = RPIPEWL1AVG(IDX)
      RBARWL2 = RPIPEWL2AVG(IDX)
#else
      RBARWL1 = ETA2(NNBB1) - PIPEHT(IDX)
      RBARWL2 = ETA2(NNBB2) - PIPEHT(IDX)
#endif
      if ((RBARWL1 < 0.d0) .and. (RBARWL2 < 0.d0)) then
         !...............WATER LEVEL ON BOTH SIDES OF BARRIER BELOW PIPE -> CASE 1
         FLUX = 0.d0
         return
      end if
      if (abs(RBARWL1 - RBARWL2) < BARMIN) then
         !...............WATER LEVEL EQUAL ON BOTH SIDES OF PIPE -> CASE 2
         FLUX = 0.d0
         return
      end if
      if ((RBARWL1 > RBARWL2) .and. (RBARWL1 > BARMIN)) then
         !...............WATER LEVEL GREATER ON THIS SIDE OF THE PIPE AND IS SUCH
         !................THAT OUTWARD DISCHARGE IS OCCURING
         !................THUS THIS SIDE IS THE SOURCE SIDE OF THE FLOW THROUGH THE PIPE
         !................NOTE THAT WE DO NOT FORCE THE SOURCE SIDE OF THE PIPE TO
         !................REMAIN WET. ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
         !................PIPE HAS BEEN DRIED. IF IT HAS, WE SHUT DOWN THE FLOW ACROSS
         !................THE PIPE
         if (RBARWL2 <= 0.d0) then
            !..................OUTWARD FREE DISCHARGE -> CASE 3
            if (NODECODE(NNBB1) == 0) then
               FLUX = 0.0d0
            else
               FLUX = -RampIntFlux &
                      *0.25d0*PI*(PIPEDIAM(IDX))**2 &
                      *(2.d0*G*RBARWL1/(1.d0 + PIPECOEF(IDX)))**0.5d0
            end if
         else
            !..................OUTWARD SUBMERGED DISCHARGE -> CASE 4
            if (NODECODE(NNBB1) == 0) then
               FLUX = 0.0d0
            else
               FLUX = -RampIntFlux &
                      *0.25d0*PI*(PIPEDIAM(IDX))**2 &
                      *(2.d0*G*(RBARWL1 - RBARWL2)/PIPECOEF(IDX))**0.5d0
            end if
            return
         end if
      end if
      if ((RBARWL2 > RBARWL1) .and. (RBARWL2 > BARMIN)) then
         !...............WATER LEVEL LOWER ON THIS SIDE OF PIPE AND IS SUCH
         !................THAT INWARD DISCHARGE IS OCCURING
         !................THUS THIS IS THE RECEIVING SIDE OF THE FLOW THROUGH THE  PIPE
         !................NOTE THAT WE DO FORCE THE RECEIVING SIDE OF THE PIPE TO
         !................REMAIN WET WHEN THERE IS FLOW ACROSS THE PIPE.
         !................ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
         !................PIPE HAS BEEN DRIED. IF IT HAS, WE SHUT DOWN THE FLOW THROUGH
         !................THE PIPE
         if (RBARWL1 <= 0) then
            !..................INWARD FREE DISCHARGE -> CASE 5
            if (NODECODE(NNBB2) == 0) then
               FLUX = 0.0d0
            else
               FLUX = RampIntFlux &
                      *0.25d0*PI*(PIPEDIAM(IDX))**2 &
                      *(2.d0*G*RBARWL2/(1.d0 + PIPECOEF(IDX)))**0.5d0
               NIBNODECODE(NNBB1) = 1
            end if
         else
            !..................INWARD SUBMERGED DISCHARGE -> CASE 6
            if (NODECODE(NNBB2) == 0) then
               FLUX = 0.0d0
            else
               FLUX = RampIntFlux &
                      *0.25d0*PI*(PIPEDIAM(IDX))**2 &
                      *(2.d0*G*(RBARWL2 - RBARWL1)/PIPECOEF(IDX))**0.5d0
               NIBNODECODE(NNBB1) = 1
            end if
         end if
         return
      end if

   end subroutine COMPUTE_CROSS_BARRIER_PIPE_FLUX

#ifdef ORIGWEIR
!-----------------------------------------------------------------------
!     S U B R O U T I N E
!       C O M P U T E _ I N T E R N A L _ B O U N D A R Y _ F L U X _ O R I G
!-----------------------------------------------------------------------
!  THIS IS THE ORIGINAL ADCIRC OVERTOPPING ROUTINE THAT ONLY CHECKS TO
!  SEE IF THE NODE IS WET BEFORE PASSING FLUX ACROSS THE WEIR
!-----------------------------------------------------------------------
   subroutine COMPUTE_INTERNAL_BOUNDARY_FLUX_ORIG( &
      BARRIER_INDEX, BOUNDARYNODE, BOUNDARY, TIMELOC, FLUX)

      use BOUNDARIES, only: BARINCFSB, BARINCFSP, NIBNODECODE, &
                            NODECODE, NVELL, NODECODE, IBCONN
      implicit none
      integer, intent(IN), target   :: BARRIER_INDEX
      integer, intent(IN), target   :: BOUNDARY
      integer, intent(IN), target   :: BOUNDARYNODE
      real(8), intent(IN)         :: TIMELOC
      real(8), intent(OUT)        :: FLUX

      integer                     :: NNBB1, NNBB2
      integer                     :: NNBB1WN, NNBB2WN
      integer                     :: FLOWDIR
      integer, pointer             :: I, J, K
      real(8)                    :: RBARWL1, RBARWL2
      real(8)                    :: RBARWL1F, RBARWL2F

      LOG_SCOPE_TRACED("COMPUTE_INTERNAL_BOUNDARY_FLUX", WEIR_BOUNDARY_TRACING)

!...............SIMPLIFY VARIABLES
      I => BARRIER_INDEX
      J => BOUNDARYNODE
      K => BOUNDARY

      NNBB1 = NBV(I) ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
      NNBB2 = IBCONN(I) ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
      if (INT_TVW) then
         if (BAR_DEG(I)) then
            call COMPUTE_BARRIER_HEIGHT(I, TIMELOC, &
                                        BARINHT1(I), BARINHT2(I))
         end if
      end if
      if (IBSTART == 0) then
         RBARWL1AVG(I) = ETA2(NNBB1) - BARINHT2(I)
         RBARWL2AVG(I) = ETA2(NNBB2) - BARINHT2(I)
         IBSTART = 1
      else
         RBARWL1AVG(I) = (ETA2(NNBB1) - BARINHT2(I) + BARAVGWT &
                          *RBARWL1AVG(I))/(1.d0 + BARAVGWT)
         RBARWL2AVG(I) = (ETA2(NNBB2) - BARINHT2(I) + BARAVGWT &
                          *RBARWL2AVG(I))/(1.d0 + BARAVGWT)
      end if
      RBARWL1 = RBARWL1AVG(I)
      RBARWL2 = RBARWL2AVG(I)
      RBARWL1F = 2.d0*RBARWL1/3.d0
      RBARWL2F = 2.d0*RBARWL2/3.d0
      FLUX = 0.d0
      if ((RBARWL1 < 0.d0) .and. (RBARWL2 < 0.d0)) then
!...............WATER LEVEL ON BOTH SIDES OF BARRIER BELOW BARRIER -> CASE 1
         FLUX = 0.d0
         return
      end if
      if (abs(RBARWL1 - RBARWL2) < 0.01d0) then
!...............WATER LEVEL EQUAL ON BOTH SIDES OF BARRIER
!................TO WITHIN TOLERANCE BARMIN -> CASE 2
         FLUX = 0.d0
         return
      end if
      if ((RBARWL1 > RBARWL2) .and. (RBARWL1 > BARMIN)) then
!...............WATER LEVEL GREATER ON THIS SIDE OF THE BARRIER AND IS SUCH
!................THAT OVERTOPPING IS OCCURING
!................THUS THIS SIDE IS THE SOURCE SIDE OF THE FLOW ACROSS THE BARRIER
!................NOTE THAT WE DO NOT FORCE THE SOURCE SIDE OF THE BARRIER TO
!................REMAIN WET. ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
!................BARRIER HAS BEEN DRIED. IF IT HAS WE SHUT DOWN THE FLOW ACROSS
!................THE BARRIER
         if (RBARWL2 > RBARWL1F) then
!..................OUTWARD SUBCRITICAL FLOW -> CASE 3
            if (NODECODE(NNBB1) == 0) then
               FLUX = 0.0d0
            else
               FLUX = -RampIntFlux*RBARWL2*BARINCFSB(I)* &
                      (2.d0*G*(RBARWL1 - RBARWL2))**0.5d0
            end if
         else
!..................OUTWARD SUPERCRITICAL FLOW -> CASE 4
            if (NODECODE(NNBB1) == 0) then
               FLUX = 0.0d0
            else
               FLUX = -RampIntFlux*BARINCFSP(I)*RBARWL1F* &
                      (RBARWL1F*G)**0.5d0
            end if
         end if
         return
      end if
      if ((RBARWL2 > RBARWL1) .and. (RBARWL2 > BARMIN)) then
!...............WATER LEVEL LOWER ON THIS SIDE OF BARRIER AND IS SUCH
!................THAT OVERTOPPING IS OCCURING
!................THUS THIS IS THE RECEIVING SIDE OF THE FLOW ACROSS THE BARRIER
!................NOTE THAT WE DO FORCE THE RECEIVING SIDE OF THE BARRIER TO
!................REMAIN WET WHEN THERE IS FLOW ACROSS THE BARRIER.
!................ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
!................BARRIER HAS BEEN DRIED. IF IT HAS, WE SHUT DOWN THE FLOW ACROSS
!................THE BARRIER
         if (RBARWL1 > RBARWL2F) then
!..................INWARD SUBCRITICAL FLOW -> CASE 5
            if (NODECODE(NNBB2) == 0) then
               FLUX = 0.0d0
            else
               FLUX = RampIntFlux*RBARWL1*BARINCFSB(I)* &
                      (2.0d0*G*(RBARWL2 - RBARWL1))**0.5d0
               NIBNODECODE(NNBB1) = 1
            end if
         else
!..................INWARD SUPERCRITICAL FLOW -> CASE 6
            if (NODECODE(NNBB2) == 0) then
               FLUX = 0.0d0
            else
               FLUX = RampIntFlux*BARINCFSP(I)*RBARWL2F* &
                      (RBARWL2F*G)**0.5d0
               NIBNODECODE(NNBB1) = 1
            end if
         end if
      end if

      return

!-----------------------------------------------------------------------
   end subroutine COMPUTE_INTERNAL_BOUNDARY_FLUX_ORIG
!-----------------------------------------------------------------------
#endif

   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------
   !     S U B R O U T I N E
   !       R E S T A R T _ V E W 6 4 _ F L U X E S
   !-----------------------------------------------------------------------
   !  Hot start: rebuild the over-wall fluxes of the IBTYPE=64 (VEW)
   !  boundary entries, which are not stored in the hotstart file. Called
   !  from HOTSTART before the first time step; on return QN1 and QN2 hold
   !  the fluxes of the two time levels before the restart, as the time
   !  step expects them before its QN0 <- QN1 <- QN2 shift.
   !
   !  The flux of a time level is computed from the water level of the
   !  previous level. The hotstart file holds the last two levels
   !  (ETA1 = n-1, ETA2 = n), so
   !    QN2 (level n)   = flux(ETA1) * ramp(n), exact;
   !    QN1 (level n-1) = flux of the level n-2 water level, which is not
   !                      stored. It is extrapolated linearly from
   !                      flux(ETA1) and flux(ETA2), or set to flux(ETA1)
   !                      when either is zero or they differ in sign
   !                      (e.g. across a change of submergence).
   !  The fluxes are computed with the same routines as in the time step,
   !  including the submergence flags; the wet/dry state of level n is used
   !  for both levels.
   !-----------------------------------------------------------------------
   subroutine RESTART_VEW64_FLUXES(TIMELOC, RAMP_PREV, RAMP_CURR)
      use GLOBAL, only: ETA1, QN1, QN2, NODECODE, NNODECODE
      use BOUNDARIES, only: NBOU, NVELL, NVEL, ISSUBMERGED64
      implicit none
      real(8), intent(IN) :: TIMELOC !...time of the hot start
      real(8), intent(IN) :: RAMP_PREV !...internal flux ramp at level n-1
      real(8), intent(IN) :: RAMP_CURR !...internal flux ramp at level n
      real(8), allocatable :: G1(:), G2(:), ETA2SAVE(:)
      integer, allocatable :: NODECODESAVE(:)
      real(8) :: RAMPSAVE, GX
      logical :: TVWSAVE
      integer :: I

      LOG_SCOPE_TRACED("RESTART_VEW64_FLUXES", WEIR_BOUNDARY_TRACING)

      allocate (G1(NVEL), G2(NVEL))
      G1 = 0d0
      G2 = 0d0

      !...Unramped fluxes; no time-varying weir updates while rebuilding
      RAMPSAVE = RAMPINTFLUX
      RAMPINTFLUX = 1d0
      TVWSAVE = INT_TVW
      INT_TVW = .false.

      !...Wet/dry state of the hot start: HOTSTART reads it into NNODECODE,
      !   and NODECODE is set from it only in the first time step. The flux
      !   routines test NODECODE, which here still holds its initial value,
      !   so that dry nodes would pass weir flow.
      NODECODESAVE = NODECODE
      NODECODE = NNODECODE

      !...Flux from the level n-1 water level
      ETA2SAVE = ETA2
      ETA2 = ETA1
      ISSUBMERGED64 = 0
      call FLUXES64(G1)
      !...Flux from the level n water level
      ETA2 = ETA2SAVE
      call FLUXES64(G2)

      RAMPINTFLUX = RAMPSAVE
      INT_TVW = TVWSAVE
      NODECODE = NODECODESAVE

      do I = 1, NVEL
         if (LBCODEI(I) /= 64) cycle
         QN2(I) = RAMP_CURR*G1(I)
         GX = G1(I)
         if (G1(I) /= 0d0 .and. G2(I) /= 0d0 .and. &
             sign(1d0, G1(I)) == sign(1d0, G2(I))) then
            GX = 2d0*G1(I) - G2(I)
            if (sign(1d0, GX) /= sign(1d0, G1(I))) GX = G1(I)
         end if
         QN1(I) = RAMP_PREV*GX
      end do

   contains

      !...Submergence flags and fluxes of all IBTYPE=64 entries, as in TIMESTEP
      subroutine FLUXES64(G)
         real(8), intent(INOUT) :: G(:)
         integer :: II, J, K
         II = 0
         do K = 1, NBOU
            select case (LBCODEI(II + 1))
            case (4, 24, 5, 25)
               II = II + NVELL(K)*2
            case (64)
               do J = 1, NVELL(K)
                  II = II + 1
                  call SET_SUBMERGED64_AT(II, J, K, TIMELOC)
               end do
               II = II + NVELL(K)
            case DEFAULT
               II = II + NVELL(K)
            end select
         end do
         II = 0
         do K = 1, NBOU
            select case (LBCODEI(II + 1))
            case (64)
               do J = 1, NVELL(K)*2
                  II = II + 1
                  call COMPUTE_INTERNAL_BOUNDARY64_FLUX(II, J, K, TIMELOC, G(II))
               end do
            case (4, 24, 5, 25)
               II = II + NVELL(K)*2
            case DEFAULT
               II = II + NVELL(K)
            end select
         end do
      end subroutine FLUXES64

      !-----------------------------------------------------------------------
   end subroutine RESTART_VEW64_FLUXES
   !-----------------------------------------------------------------------

end module WEIR_FLUX
!-----------------------------------------------------------------------
