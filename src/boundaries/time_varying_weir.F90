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
!     M O D U L E  T I M E _ V A R Y I N G _ W E I R _ B O U N D A R Y
!-----------------------------------------------------------------------
!  THIS MODULE CONTAINS ALL THE ROUTINES NECESSARY FOR SPECIFICATION OF
!  TIME VARYING WEIR BOUNDARIES. THESE ROUTINES ARE USED TO SETUP AND
!  CALCULATE THE ELEVATIONS THROUGHOUT THE SIMULATION
!-----------------------------------------------------------------------
module mod_time_varying_weir_boundary
   use SIZES, only: MYPROC, LOCALDIR, GLOBALDIR
   use mod_terminate, only: terminate
   use mod_logging, only: SCREENUNIT, logMessage, allMessage, &
                          setMessageSource, unsetMessageSource, ERROR, INFO, ECHO, &
                          DEBUG, WARNING, ScratchMessage, ScreenMessage
   use KDTREE2_MODULE, only: KDTREE2, KDTREE2_CREATE, KDTREE2_RESULT, &
                             KDTREE2_N_NEAREST, KDTREE2_DESTROY
   use mod_weir_data, only: BARINHT1, BARINHT2, BARLANHT1, BARLANHT2, &
                            BARMIN, BARMIN64, SUBMIN64, &
                            EXT_TVW, INT_TVW, FOUND_TVW_NML, &
                            ALLOCATE_WEIRS

   implicit none

   real(8), allocatable :: BARHT_FINAL(:) !...Final elevation for time varying boundary
   real(8), allocatable :: BAR_LOCATIONS(:, :) !...Location array for boundary conditions
   real(8), allocatable :: BAR_DEG_START(:) !...Start time for time varying boundary
   real(8), allocatable :: BAR_DEG_END(:) !...End time for time varying boundary
   real(8), allocatable :: BAR_ETA_MAX(:) !...Used if the boundary changes at a critical water surface elevation
   real(8), allocatable :: BAR_FAILURE_START(:) !...Used to track the first if ETA_MAX has been exceeded
   real(8), allocatable :: BAR_FAILURE_DURATION(:) !...Amount of time it takes an ETA_MAX barrier to fail to ZF
   logical, allocatable :: BAR_DEG(:) !...T/F this node is a time varying boundary
   integer, allocatable :: BAR_VARYTYPE(:) !...Type of variation to apply to this boundary node
   type(KDTREE2), pointer :: BARRIER_SEARCHTREE !...Search tree for boundary nodes

#ifdef CMPI
   real(8), allocatable :: GHOST_LOCATIONS(:, :) !...Location array for the ghost nodes
   type(KDTREE2), pointer :: GHOST_SEARCHTREE !...Search tree for ghost nodes
#endif

   !..CHARACTERIZES A SECTION OF THE SCHEDULE
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
   end type

   !..CHARACTERIZES THE SCHEDULE. BASED ON NUMBER OF FILES
   type BARRIER_SCHEDULE_T
      integer :: SID !...Unique ID corresponding to this schedule
      integer :: NSECTIONS !...Number of sections in the schedule
      type(BSCHED), allocatable :: SECTION(:) !...List of all the schedule sections
   end type

   !..SCHEDULE THAT POINTS TO MAIN SET THAT HOUSES THE ACTUAL DATA
   !  AND IS BASED UPON THE INDIVIDUAL BOUNDARY NODE
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
   end type
   type(BARRIER_SCHEDULE_T), allocatable, target :: SCHEDULE(:)
   type(BARRIER_SCHEDULE_P), allocatable :: BAR_SCHEDULE(:)
   character(1024), allocatable, target :: SCHEDULE_LIST(:)

   !..THESE VARIABLES ARE PART OF THE NAMELIST THAT IS READ IN
   character(200) :: ScheduleFile
   real(8) :: X1, X2, Y1, Y2
   real(8) :: ZF, ETA_MAX
   real(8) :: TimeStartDay, TimeStartHour, TimeStartMin
   real(8) :: TimeStartSec
   real(8) :: TimeEndDay, TimeEndHour, TimeEndMin, TimeEndSec
   real(8) :: FailureDurationDay, FailureDurationHour
   real(8) :: FailureDurationMin, FailureDurationSec
   real(8) :: HOTADD
   integer :: VaryType, HOT
   integer :: LOOP, NLOOPS

   integer :: NSCHEDULES

   !..THE TIME VARYING WEIR NAMELIST
   namelist /TimeVaryingWeir/ &
      X1, Y1, X2, Y2, &
      VaryType, ZF, &
      ETA_MAX, &
      TimeStartDay, TimeStartHour, TimeStartMin, TimeStartSec, &
      TimeEndDay, TimeEndHour, TimeEndMin, TimeEndSec, &
      FailureDurationDay, FailureDurationHour, &
      FailureDurationMin, FailureDurationSec, &
      ScheduleFile, HOT, LOOP, NLOOPS

   public :: WEIR_SETUP, PARSE_TIME_VARYING_WEIR_INFO
   public :: BAR_DEG, COMPUTE_BARRIER_HEIGHT

   private

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
      use BOUNDARIES, only: NFLUXIB, NFLUXIBP, NFLUXB
      implicit none

      call setMessageSource("WEIR_SETUP")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif
      if ((NFLUXIB == 1) .or. (NFLUXIBP == 1) .or. (NFLUXB == 1)) then

         !...Allocate the new weir arrays at both time levels
         call ALLOCATE_WEIRS()

         !...Parse the time varying weir file on the local processor
         if (found_tvw_nml) then
            call PARSE_TIME_VARYING_WEIR_INFO()
         end if

      end if

#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()
      return

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
#ifdef CMPI
      use MESH, only: NP
#endif
      use MESH, only: X, Y
      use BOUNDARIES, only: NVEL, LBCODEI, NBV
#ifdef CMPI
      use MESSENGER, only: RESNODE
#endif

      implicit none

      integer :: I
      integer :: NWEIR
      integer :: IDX
#ifdef CMPI
      integer :: NGHOST
#endif

      call setMessageSource("ALLOCATE_TIMEVARYINGWEIRS")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif

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

#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()
      return

   end subroutine ALLOCATE_TIMEVARYINGWEIRS
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E  F I N D _ B O U N D A R Y _ N O D E S
   !-----------------------------------------------------------------------
   !> THIS SUBROUTINE FINDS THE NODE THAT MATCHS THAT WHICH IS SPECIFIED
   !> AS A TIME VARYING BOUNDARY POINT IN THE INPUT FILE. THE NEAREST
   !> BOUNDARY NODE IS LOCATED. TYPE 4/24/5/25 BOUNDARIES ARE SPECIFIED WITH BOTH
   !> NODES. TYPE 3/13/23 BOUNDARIES ARE SPECIFIED WITH A SINGLE NODE. THIS SCHEME
   !> IS USED SO THAT MESH NUMBERING AND BOUNDARY CONDITION ORDERING IN THE FORT.14
   !> CANNOT INVALIDATE THE TIME VARYING WEIR INPUT FILE. CHECKS ARE IN PLACE TO ENSURE
   !> THE SAME BOUNDARY NODE IS NOT LOCATED TWICE.
   !> @param LAT - Latitude of the boundary node
   !> @param LON - Longitude of the boundary node
   !> @param IDX - The index of the boundary node
   !-----------------------------------------------------------------------
   subroutine FIND_BOUNDARY_NODES(LON, LAT, IDX)
      use MESH, only: ICS, SLAM0, SFEA0, &
                      DRVSPCOORSROTS, CYLINDERMAP
      use GLOBAL, only: DEG2RAD, IFSPROTS

      implicit none
      real(8), intent(in) :: LAT
      real(8), intent(in) :: LON
      real(8) :: LAT_TEMP, LON_TEMP
      real(8) :: LATR, LONR
      real(8) :: EPS
      real(8) :: X, Y
      integer, intent(out) :: IDX
      integer, parameter :: SEARCHDEPTH = 2
      type(KDTREE2_RESULT) :: KDRESULTS(SEARCHDEPTH)

      call setMessageSource("FIND_BOUNDARY_NODES")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif

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
                             QV=(/X, Y/), NN=SEARCHDEPTH, RESULTS=KDRESULTS)

      ! THIS NEEDS SOME WORK, HOWEVER, IM NOT SURE HOW TO HANDLE THE POSSIBILITY
      ! THAT A NODE HAS TWO BOUNDARY CONDITIONS BESIDES MAKING IT A FATAL ERROR
      if (abs(KDRESULTS(1)%DIS - KDRESULTS(2)%DIS) <= EPS) then
         call allMessage(ERROR, &
                         "MULTIPLE LOCATIONS FOUND FOR TIME VARYING "// &
                         "BOUNDARY NODES")
         call terminate()
      end if

      if (abs(KDRESULTS(1)%DIS) > EPS) then
#ifdef CMPI
!........CHECK TO MAKE SURE IT IS NOT INVALID BECAUSE IT
!                   IS A GHOST NODE
         call KDTREE2_N_NEAREST(TP=GHOST_SEARCHTREE, QV=(/X, Y/), &
                                NN=SEARCHDEPTH, RESULTS=KDRESULTS)
         if (abs(KDRESULTS(1)%DIS) > EPS) then
!........THIS IS NOT A GHOST NODE EITHER.

            write (ScratchMessage, '(A,F0.9,A,F0.9,A)') &
               "SPECIFIED NODE LOCATION X=", &
               LAT, " Y=", LON, &
               " NOT FOUND IN LIST OF BOUNDARY NODES."
            call allMessage(ERROR, ScratchMessage)
            call terminate()
         else
            IDX = -1
#if defined(TVW_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG, "Return")
#endif
            call unsetMessageSource()
            return
         end if

#else
         write (ScratchMessage, '(A,F0.9,A,F0.9,A)') &
            "SPECIFIED NODE LOCATION X=", &
            LAT, " Y=", LON, &
            " NOT FOUND IN LIST OF BOUNDARY NODES."
         call allMessage(ERROR, ScratchMessage)
         call terminate()

#endif
      end if

      IDX = int(BAR_LOCATIONS(3, KDRESULTS(1)%IDX))

#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()

      return

   end subroutine FIND_BOUNDARY_NODES
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E  C O M P U T E _ B A R R I E R _ H E I G H T
   !-----------------------------------------------------------------------
   !> THIS SUBROUTINE WILL DETERMINE IF THE BARRIER HEIGHTS ARE SUPPOSED TO CHANGE AND
   !> WHEN. IT IS CALLED FOR ALL VARIABLE WEIR TYPE BOUNDARIES AND RETURNS THE INITIAL
   !> VALUE FOR THE BOUNDARY IF NO ACTION IS REQUIRED. BOUNDARY IS NOT ALLOWED TO
   !> DECREASE BELOW THE ELEVATION OF THE SURROUNDING TOPOGRAPHY. THIS ESSENTIALLY
   !> FUNCTIONS AS A WRAPPER FOR ALL THE BOUNDARY VARIATION ROUTINES.
   !>
   !> @param MYIDX - The index of the boundary node
   !> @param TIMELOC - The current time
   !> @param BAR_HEIGHT_CURRENT - The current height of the boundary
   !> @param BAR_HEIGHT - The new height of the boundary
   !-----------------------------------------------------------------------
   subroutine COMPUTE_BARRIER_HEIGHT(MYIDX, TIMELOC, &
                                     BAR_HEIGHT_CURRENT, BAR_HEIGHT)

      use BOUNDARIES, only: NBV

      implicit none

      real(8), parameter :: EPS = epsilon(1.0d0)

      integer, intent(in) :: MYIDX
      character(20) :: VC
      real(8), intent(in) :: TIMELOC
      real(8), intent(in) :: BAR_HEIGHT_CURRENT
      real(8), intent(out) :: BAR_HEIGHT

      call setMessageSource("COMPUTE_BARRIER_HEIGHT")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif

      !...SET DEFAULT RETURN VALUE SO WE CAN DUCK OUT AT ANY TIME
      BAR_HEIGHT = BAR_HEIGHT_CURRENT

      !..SELECT THE APPROPRIATE TYPE OF CHANGE AND CALL THE ASSOCIATED
      !  ROUTINE. FUTURE DEVELOPMENT CAN TAKE PLACE HERE BY ADDING NEW
      !  SUBROUTINES THAT ARE CALLED IN A SIMLIAR WAY. FOR EXAMPLE,
      !  A SUBROUTINE MIGHT BE WRITTEN THAT COMPUTES WAVE FORCES AND
      !  ALTERS THE HEIGHT OF THE BARRIER BASED UPON THAT.
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
         call allMessage(ERROR, "Invalid Barrier " &
                         //"Variation")
         call terminate()
      end select

      !..WRITE SCREEN/UNIT 16 INFORMATION AT START/END
      select case (BAR_VARYTYPE(MYIDX))
      case (1, 2)
         if (abs(BAR_DEG_START(MYIDX) - TIMELOC) <= EPS) then
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
         if (abs(BAR_DEG_END(MYIDX) - TIMELOC) <= EPS) then
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

#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()
      return

   end subroutine COMPUTE_BARRIER_HEIGHT
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E
   !       C O M P U T E _ B A R R I E R _ H E I G H T _ L I N E A R
   !-----------------------------------------------------------------------
   !> THIS SUBROUTINE COMPUTES THE NEW BARRIER ELEVATION AT THE
   !> SPECIFIED BOUNDARY NODE. THE INTERPOLATION FROM STARTING
   !> ELEVATION (FORT.14) AND FINAL ELEVATION IS LINEAR
   !> OVER THE SPECIFIED TIME AND DOES NOT STOP ONCE
   !> IT HAS BEGUN.
   !>
   !> @param MYIDX - The index of the boundary node
   !> @param TIMELOC - The current time
   !> @param BAR_HEIGHT_CURRENT - The current height of the boundary
   !> @param BAR_HEIGHT - The new height of the boundary
   !> @param[optional] BARHT_START - The starting height of the boundary
   !-----------------------------------------------------------------------
   subroutine COMPUTE_BARRIER_HEIGHT_LINEAR(MYIDX, &
                                            TIMELOC, BAR_HEIGHT_CURRENT, BAR_HEIGHT, BARHT_START)

      use MESH, only: DP
      use BOUNDARIES, only: BARINHT, BARLANHT, LBCODEI, &
                            NBV, IBCONN

      implicit none

      real(8), parameter :: EPS = epsilon(1.0d0)

      integer, intent(in) :: MYIDX
      real(8), intent(in) :: TIMELOC
      real(8), intent(in) :: BAR_HEIGHT_CURRENT
      real(8), intent(out) :: BAR_HEIGHT
      real(8), intent(in), optional :: BARHT_START
      real(8) :: BAR_DZ, BAR_DT
      real(8) :: BAR_START
      real(8) :: DEPTH

      call setMessageSource("COMPUTE_BARRIER_HEIGHT_LINEAR")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif

      !..SET DEFAULT RETURN VALUE
      BAR_HEIGHT = BAR_HEIGHT_CURRENT

      !..GET OUT OF HERE AT THE FIRST OPPORTUNITY
      if (BAR_VARYTYPE(MYIDX) == 1) then
         if ((BAR_DEG_START(MYIDX) > TIMELOC) .or. &
             (BAR_DEG_END(MYIDX) < TIMELOC)) then
#if defined(TVW_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG, "Return")
#endif
            call unsetMessageSource()
            return
         end if
      end if

      !..COMPUTE THE NEW BOUNDARY HEIGHT
      if (present(BARHT_START)) then
         !...Schedule style boundaries specify
         !   their own initial elevation, so dont use
         !   the elevation in the fort.14 weirs
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
#if defined(TVW_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG, "Return")
#endif
            call unsetMessageSource()
            return
         end select
      end if

      if (abs(BAR_DEG_END(MYIDX) - TIMELOC) < EPS) then
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
         !...WE'RE TOO LOW, SHUT DOWN BARRIER, SET TO MINIMUM
         if (BAR_VARYTYPE(MYIDX) == 3) then
            BAR_HEIGHT = DEPTH
         else
            if (BARHT_FINAL(MYIDX) >= DEPTH) then
               !...SET TO USER SPECIFIED MIN
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

#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()
      return

   end subroutine COMPUTE_BARRIER_HEIGHT_LINEAR
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E
   !       C O M P U T E _ B A R R I E R _ H E I G H T _ E T A M A X
   !-----------------------------------------------------------------------
   !> THIS SUBROUTINE CHECKS TO SEE IF THE BARRIER HAS REACHED
   !> ITS ETAMAX FAILURE CRITERIA. IF IT HAS, THE LINEAR BARRIER
   !> HEIGHT CHANGE ROUTINE IS CALLED WITH THE VARIABLES CREATED
   !> WHEN THIS ROUTINE FIRST DECIDED THAT THE ETAMAX CONDITION
   !> WAS MET
   !>
   !> @param MYIDX - The index of the boundary node
   !> @param TIMELOC - The current time
   !> @param BAR_HEIGHT_CURRENT - The current height of the boundary
   !> @param BAR_HEIGHT - The new height of the boundary
   !-----------------------------------------------------------------------
   subroutine COMPUTE_BARRIER_HEIGHT_ETAMAX(MYIDX, TIMELOC, &
                                            BAR_HEIGHT_CURRENT, BAR_HEIGHT)

      use GLOBAL, only: ETA2, NNODECODE
      use BOUNDARIES, only: NBV, IBCONN

      implicit none

      real(8), parameter :: EPS = epsilon(1.0d0)

      integer, intent(in) :: MYIDX
      real(8), intent(in) :: TIMELOC
      real(8), intent(in) :: BAR_HEIGHT_CURRENT
      real(8), intent(out) :: BAR_HEIGHT

      call setMessageSource("COMPUTE_BARRIER_HEIGHT_ETAMAX")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif

      !...............SET DEFAULT RETURN VALUE
      BAR_HEIGHT = BAR_HEIGHT_CURRENT

      !..CHECK IF THE ETA_MAX HAS BEEN EXCEEDED ALREADY
      if (abs(BAR_FAILURE_START(MYIDX) + 1d0) <= EPS) then
         if (ETA2(NBV(MYIDX)) > BAR_ETA_MAX(MYIDX) .and. NNODECODE(NBV(MYIDX)) == 1) then
            !..WE NEED TO INITIALIZE OUR START TIME FOR THE FAILURE
            BAR_DEG_START(MYIDX) = TIMELOC
            BAR_DEG_END(MYIDX) = TIMELOC + &
                                 BAR_FAILURE_DURATION(MYIDX)
            BAR_FAILURE_START(MYIDX) = 1d0
         elseif (ETA2(IBCONN(MYIDX)) > BAR_ETA_MAX(MYIDX) .and. NNODECODE(IBCONN(MYIDX)) == 1) then
            !..WE NEED TO INITIALIZE OUR START TIME FOR THE FAILURE
            BAR_DEG_START(MYIDX) = TIMELOC
            BAR_DEG_END(MYIDX) = TIMELOC + &
                                 BAR_FAILURE_DURATION(MYIDX)
            BAR_FAILURE_START(MYIDX) = 1d0
         else
#if defined(TVW_TRACE) || defined(ALL_TRACE)
            call allMessage(DEBUG, "Return")
#endif
            call unsetMessageSource()
            return
         end if
      end if

      !..IF WE ARE HERE, ETA_MAX WAS PREVIOUSLY EXCEEDED AND THE BARRIER
      !  IS DEGRADING AS PRESCRIBED
      call COMPUTE_BARRIER_HEIGHT_LINEAR(MYIDX, &
                                         TIMELOC, BAR_HEIGHT_CURRENT, BAR_HEIGHT)
#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()
      return

   end subroutine COMPUTE_BARRIER_HEIGHT_ETAMAX
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E
   !       C O M P U T E _ B A R R I E R _ H E I G H T _ S C H E D U L E
   !-----------------------------------------------------------------------
   !> THIS SUBROUTINE COMPUTES THE VARYING WEIR ELEVATION BASED UPON A
   !> INPUT FILE THAT CONTAINS A SCHEDULE. THIS SCHEDULE CAN BE USED TO
   !> HELP REPRESENT THINGS LIKE GATE OPERATIONS AND OTHER DYNAMIC
   !> PROCESSES IN THE MODEL THAT BEHAVE IN A KNOWN REPETITIVE WAY
   !>
   !> @param MYIDX - The index of the boundary node
   !> @param TIMELOC - The current time
   !> @param BAR_HEIGHT_CURRENT - The current height of the boundary
   !> @param BAR_HEIGHT - The new height of the boundary
   !-----------------------------------------------------------------------
   subroutine COMPUTE_BARRIER_HEIGHT_SCHEDULE(MYIDX, TIMELOC, &
                                              BAR_HEIGHT_CURRENT, BAR_HEIGHT)
      use MESH, only: DP
      use BOUNDARIES, only: IBCONN, NBV, LBCODEI, BARLANHT, BARINHT
      implicit none

      real(8), parameter :: EPS = epsilon(1.0d0)

      integer, intent(in) :: MYIDX
      real(8), intent(in) :: TIMELOC
      real(8), intent(in) :: BAR_HEIGHT_CURRENT
      real(8), intent(out) :: BAR_HEIGHT

      real(8) :: PREVEND
      real(8) :: BARHT_START
      real(8) :: OFFSET
      integer :: MYSEC
      integer :: NNBB1, NNBB2

      call setMessageSource("COMPUTE_BARRIER_HEIGHT_SCHEDULE")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif

      !..SET DEFAULT RETURN VALUE
      BAR_HEIGHT = BAR_HEIGHT_CURRENT

      !..Check if we have reached the offset time yet
      !               (Time added to beginning of a schedule)
      if (TIMELOC < BAR_SCHEDULE(MYIDX)%OFFSET) then
#if defined(TVW_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG, "Return")
#endif
         call unsetMessageSource()
         return
      end if

      !..SIMPLIFY THESE VARIABLES
      MYSEC = BAR_SCHEDULE(MYIDX)%MYSEC
      PREVEND = BAR_SCHEDULE(MYIDX)%PREVEND
      OFFSET = BAR_SCHEDULE(MYIDX)%OFFSET

      !..IF THIS IS THE FIRST TIMESTEP IN THIS SECTINO OF THE
      !               SCHEDULE, THEN SET UP THE VARIABLES NEEDED
      if (abs(BAR_SCHEDULE(MYIDX)%SECTION(MYSEC)%BAR_DEG_START - (TIMELOC - PREVEND - OFFSET)) <= EPS) then

         !..THIS IS THE BARRIER HEIGHT AT THE START OF THIS
         !                   PORTION OF THE SCHEDULE. SAVE IT FOR LINEAR
         !                   INTERPOLATION
         BARHT_START = BARINHT2(MYIDX)
         BAR_SCHEDULE(MYIDX)%BARHT_START = BARHT_START

         !..GET THE NODE NUMBER
         NNBB1 = NBV(MYIDX)

         !..INFORM SCREEN AND UNIT 16 WE ARE HERE
#ifdef CMPI
         write (ScratchMessage, 101) MYSEC, NNBB1, TIMELOC, MYPROC
#else
         write (ScratchMessage, 102) MYSEC, NNBB1, TIMELOC
#endif
         call allMessage(INFO, ScratchMessage)

         !..ALTER THE VARIABLES USED IN COMPUTE_BARRIER_HEIGHT_LINEAR
         !  SO IT CAN CORRECTLY INTERPOLATE NEW BARRIER VALUES
         BAR_DEG_START(MYIDX) = &
            BAR_SCHEDULE(MYIDX)%SECTION(MYSEC)%BAR_DEG_START &
            + PREVEND + OFFSET
         BAR_DEG_END(MYIDX) = &
            BAR_SCHEDULE(MYIDX)%SECTION(MYSEC)%BAR_DEG_END &
            + PREVEND + OFFSET

         !..CHECK THE ZF VALUE FOR FLAGS. THESE FLAGS ARE ADDED
         !  FOR CONVENIENCE TO SET THE BARRIER TO PREDETERMINED
         !  VALUES BASED UPON MESH GEOMETRY.
         select case (int(BAR_SCHEDULE(MYIDX)% &
                          SECTION(MYSEC)%BARHT_FINAL))
         case (-99990)
            !...Set to fort.14 nodal elevation
            select case (LBCODEI(MYIDX))
            case (3, 13, 23)
               BARHT_FINAL(MYIDX) = -DP(NNBB1)
            case (4, 14, 24, 64, 5, 25)
               NNBB2 = IBCONN(MYIDX)
               BARHT_FINAL(MYIDX) = &
                  min(-DP(NNBB1), -DP(NNBB2))
            end select
         case (-99991)
            !...Add BARHT_DELTA to present height
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
            !...Subtract BARHT_DELTA from present height
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
            !...Set to weir elevation in fort.14
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

         !..USE THE NORMAL LINEAR INTERPOLATION ROUTINE TO
         !  GET THE NEW BARRIER ELEVATION
         call COMPUTE_BARRIER_HEIGHT_LINEAR(MYIDX, TIMELOC, &
                                            BAR_HEIGHT_CURRENT, BAR_HEIGHT, BARHT_START)

         !..THIS IS THE END OF A SCHEDULE SECTION, WE NEED TO
         !  INCREMENT COUNTERS AND PREPARE FOR THE NEXT SECTION
         !  OF THE SCHEDULE
      elseif (abs(BAR_SCHEDULE(MYIDX)%SECTION(MYSEC)%BAR_DEG_END - (TIMELOC - PREVEND - OFFSET)) <= EPS) then

         !..GET NODE NUMBER AND INFORM SCREEN AND UNIT 16
         NNBB1 = NBV(MYIDX)
#ifdef CMPI
         write (ScratchMessage, 103) MYSEC, NNBB1, TIMELOC, MYPROC
#else
         write (ScratchMessage, 104) MYSEC, NNBB1, TIMELOC
#endif
         call allMessage(INFO, ScratchMessage)

         !..CALL THE LAST LINEAR INTERPOLATION
         BARHT_START = BAR_SCHEDULE(MYIDX)%BARHT_START
         call COMPUTE_BARRIER_HEIGHT_LINEAR(MYIDX, TIMELOC, &
                                            BAR_HEIGHT_CURRENT, BAR_HEIGHT, BARHT_START)

         !..INCREMENT SCHEDULE SECTION
         BAR_SCHEDULE(MYIDX)%MYSEC = MYSEC + 1

         !..CHECK TO SEE IF THIS IS THE END OF THE SCHEDULE.
         !  IF THE LOOP FLAG HAS BEEN SET TO -1, WE LOOP THE
         !  SCHEDULE INFINITELY. IF LOOP SET TO 1, WE WILL
         !  INCREMENT OUR LOOP COUNTER TO CHECK AND SEE IF WE
         !  HAVE DONE THE CORRECT NUMBER OF SCHEDULE LOOPS
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
               !..THIS NODE SHOULD NO LONGER CHANGE, SO
               !  SET THE APPROPRIATE FLAG
               BAR_DEG(MYIDX) = .false.
            end if
         end if

         !..THE CASE FOR THE MIDDLE OF THE SCHEDULE. CONTINUE WITH
         !  THE INTERPOLATION
      elseif ((BAR_SCHEDULE(MYIDX)%SECTION(MYSEC)%BAR_DEG_START &
               < TIMELOC - PREVEND - OFFSET) .and. &
              (BAR_SCHEDULE(MYIDX)%SECTION(MYSEC)%BAR_DEG_END &
               > TIMELOC - PREVEND - OFFSET)) then
         BARHT_START = BAR_SCHEDULE(MYIDX)%BARHT_START
         call COMPUTE_BARRIER_HEIGHT_LINEAR(MYIDX, TIMELOC, &
                                            BAR_HEIGHT_CURRENT, BAR_HEIGHT, BARHT_START)
      end if

#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()
      return
#ifdef CMPI
101   format("INFO: SCHEDULE STYLE BOUNDARY BEGAN SEGMENT ", &
             I4, " AT NODE = ", I9, " AT TIME = ", E15.8, &
             "ON MYPROC = ", I4)
103   format("INFO: SCHEDULE STYLE BOUNDARY CONCLUDED", &
             " SEGMENT ", I4, " AT NODE = ", I9, " AT TIME = ", &
             E15.8, " ON MYPROC = ", I4)
#else
102   format("INFO: SCHEDULE STYLE BOUNDARY BEGAN SEGMENT ", &
             I4, " AT NODE = ", I9, " AT TIME = ", E15.8)
104   format("INFO: SCHEDULE STYLE BOUNDARY CONDLUDED", &
             " SEGMENT ", I4, " AT NODE = ", I9, " AT TIME = ", &
             E15.8)
#endif

   end subroutine COMPUTE_BARRIER_HEIGHT_SCHEDULE
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E
   !       P A R S E _ T I M E _ V A R Y I N G _ W E I R _ I N F O
   !-----------------------------------------------------------------------
   !> THIS SUBROUTINE READS THE NAMELIST STYLE INPUT FILE FROM THE USER AND
   !> SETS UP THE NECESSARY INFORMATION FOR USE DURING THE SIMLUATION
   !> AS WELL AS PERFORMING VARIOUS SANITY CHECKS ON THE INPUT TO AVOID
   !> ERRORS DURING THE SIMULATION.
   !-----------------------------------------------------------------------
   subroutine PARSE_TIME_VARYING_WEIR_INFO()

      use GLOBAL, only: ITHS, DTDP, IHOT, &
                        tvw_file, openFileForRead
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

      call setMessageSource("PARSE_TIME_VARYING_WEIR_INFO")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif

      call openFileForRead(99, trim(INPUTDIR)//'/'// &
                           trim(tvw_file), IOS)
      if (IOS /= 0) then
         NTIMEVARYINGWEIRS = 0
#ifdef CMPI
         write (ScratchMessage, 101) MYPROC
#else
         write (ScratchMessage, 102)
#endif
         call allMessage(INFO, ScratchMessage)
#if defined(TVW_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG, "Return")
#endif
         call unsetMessageSource()
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
#if defined(TVW_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG, "Return")
#endif
         call unsetMessageSource()
         return
      end if

      open (UNIT=98, FILE=trim(LOCALDIR)//'/namelist.scratch', &
            STATUS='REPLACE', ACTION='WRITE')

      !..READ IN UNFORMATTED NAMELIST, WRITE TO FORMATTED NAMELIST
      do I = 1, NTIMEVARYINGWEIRS
         !..MODIFY LINE FOR NAMELIST FORMATTING
         read (99, '(A)') InputString
         modifiedString = "&TimeVaryingWeir "// &
                          trim(adjustl(InputString))//" /"
         write (98, '(A)') trim(adjustl(modifiedString))
      end do
      close (98)

      !..OPEN FORMATTED NAMELIST FILE, READ IN NUMBER OF SCHEDULES
      !  WE WILL NEED TO PARSE
      open (UNIT=98, FILE=trim(LOCALDIR)//'/namelist.scratch', &
            ACTION="READ")
      NSCHEDULES = 0
      NSCHEDULES_RAW = 0
      do I = 1, NTIMEVARYINGWEIRS
         call NULLIFY_TVW_NML()
         read (98, NML=TimeVaryingWeir, IOSTAT=IOS, ERR=200)
         if (.not. ISNULL(S=SCHEDULEFILE)) then
            NSCHEDULES_RAW = NSCHEDULES_RAW + 1
         end if
      end do
      if (NSCHEDULES_RAW > 0) then
         rewind (98)
         IDX = 0
         allocate (SCHEDULE_LIST_RAW(1:NSCHEDULES_RAW))
         !...READ BACK IN A LIST CONTAINING THE SCHEDULE FILE NAMES
         do I = 1, NTIMEVARYINGWEIRS
            call NULLIFY_TVW_NML()
            read (98, NML=TimeVaryingWeir, IOSTAT=IOS)
            if (.not. ISNULL(S=SCHEDULEFILE)) then
               IDX = IDX + 1
               SCHEDULE_LIST_RAW(IDX) = SCHEDULEFILE
            end if
         end do
         !..FIND OUT HOW MANY UNIQUE SCHEDULES THERE ARE
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
      !..OPEN FORMATTED NAMELIST FILE, READ IN
      open (UNIT=98, FILE=trim(LOCALDIR)//'/namelist.scratch', &
            ACTION="READ")
      do I = 1, NTIMEVARYINGWEIRS
         !..READ IN NAMELIST LINE
         call NULLIFY_TVW_NML()
         read (98, NML=TimeVaryingWeir, ERR=200, &
               IOSTAT=IOS)

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

         !..BEGIN SANITY CHECK ON WHAT HAS BEEN SPECIFIED IN NAMELIST
         select case (IHOT)
         case (17, 67, 68, 367, 368, 567, 568)
            if (ISNULL(I=HOT) .or. (HOT == 0)) then
               HOTADD = 0d0
            elseif (HOT == 1) then
               HOTADD = DTDP*ITHS
               OFFSET = OFFSET + HOTADD
               write (*, *) DTDP, ITHS
               write (*, *) HOTADD, OFFSET
            else
               call allMessage(ERROR, &
                               "INCORRECT"// &
                               " HOT START VALUE SPECIFIED. HOT=1 OR "// &
                               "HOT=0 FOR HOT START RELATIVE TIME "// &
                               "VARYING WEIRS.")
               call terminate()
            end if
         case DEFAULT
            HOTADD = 0d0
         end select
         if (ISNULL(I=VARYTYPE)) then
            call allMessage(ERROR, &
                            "YOU MUST "// &
                            "SPECIFY VARYTYPE= IN THE FORT.15 INPUT "// &
                            "FILE.")
            call terminate()
         end if
         select case (VARYTYPE)
         case (1, 2, 3)
            if (VARYTYPE == 1) then
               if (ISNULL(X1) .or. ISNULL(Y1) .or. &
                   ISNULL(ZF)) then
                  call allMessage(ERROR, &
                                  "YOU MUST SPECIFY X1=, Y1=, and "// &
                                  "ZF= ")
                  call terminate()
               end if
               if (ISNULL(TimeStartday) .and. &
                   ISNULL(TimeStartHour) .and. &
                   ISNULL(TimeStartSec)) then
                  call allMessage(ERROR, &
                                  "TIME VARYING BOUNDARY START TIME"// &
                                  " MUST BE SPECIFIED.")
                  call terminate()
               end if
               if (ISNULL(TimeEndday) .and. &
                   ISNULL(TimeEndHour) .and. &
                   ISNULL(TimeEndSec)) then
                  call allMessage(ERROR, &
                                  "TIME VARYING BOUNDARY END TIME"// &
                                  " MUST BE SPECIFIED.")
                  call terminate()
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
                  call allMessage(ERROR, &
                                  "YOU MUST SPECIFY X1=, Y1=, "// &
                                  "and ZF= ")
                  call terminate()
               end if
               if ((ISNULL(FAILUREDURATIONDAY) .and. &
                    ISNULL(FAILUREDURATIONHOUR) .and. &
                    ISNULL(FAILUREDURATIONMIN) .and. &
                    ISNULL(FAILUREDURATIONSEC)) .or. &
                   ISNULL(ETA_MAX)) then
                  call allMessage(ERROR, &
                                  " YOU MUST SPECIFY A FAILURE"// &
                                  " DURATION AND MAXIMUM WATER"// &
                                  " SURFACE BEFORE ELEVATION"// &
                                  " CHANGE.")
                  call terminate()
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
                  call allMessage(ERROR, &
                                  "YOU MUST SPECIFY X1= and Y1=")
                  call terminate()
               end if
               if (ISNULL(S=SCHEDULEFILE)) then
                  call allMessage(ERROR, &
                                  "YOU MUST SPECIFY SCHEDULEFILE= "// &
                                  "FOR VARYTYPE=3.")
                  call terminate()
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
            !..END SANITY CHECK

            !..BEGIN ASSIGNING BOUNDARY CONDITIONS THEIR ATTRIBUTES
            call FIND_BOUNDARY_NODES(X1, Y1, IDX)
            if (IDX == -1) cycle
            select case (LBCODEI(IDX))
            case (3, 13, 23)
               !..A ONE SIDED STYLE WEIR, WE DONT NEED X2 AND Y2
               call ASSIGN_TVW_TIMING(VARYTYPE, IDX)
               EXT_TVW = .true.
               if (VARYTYPE == 1) then
                  if (BAR_DEG_END(IDX) <= 0d0) then
                     call allMessage(ERROR, &
                                     "TIME VARYING BOUNDARY END "// &
                                     "TIME MUST BE GREATER THAN "// &
                                     "ZERO.")
                     call terminate()
                  end if
               elseif (VARYTYPE == 2) then
                  if (BAR_FAILURE_DURATION(IDX) <= 0d0) then
                     call allMessage(ERROR, &
                                     "FAILURE_DURATION MUST BE "// &
                                     "GREATER THAN ZERO")
                     call terminate()
                  end if
               elseif (VARYTYPE == 3) then
                  BAR_SCHEDULE(IDX)%OFFSET = OFFSET
               end if

            case (4, 24, 64, 5, 25)
               !..A TWO SIDED STYLE WEIR, WE NEED X2 AND Y2
               call FIND_BOUNDARY_NODES(X2, Y2, IDX2)
               if (IDX2 == -1) then
                  call allMessage(ERROR, "BOUNDARY CONDITION IMPROPERLY "// &
                                  "SPLIT INTO GHOST NODE SPACE!")
                  call terminate()
               end if
               INT_TVW = .true.

               !..WE NEED TO CHECK CONNECTIVITY ACCROSS WEIR TO
               !                           MAKE SURE THIS PAIR IS VALID
               if (NBV(IDX) /= IBCONN(IDX2)) then
                  !..THIS IS NOT THE CONNECTIVITY WE EXPECTED, GOOD NIGHT
                  write (ScratchMessage, '(A)') &
                     "CONNECTIVITY ACROSS INTERNAL"// &
                     " TIME VARYING BOUNDARY IS INVALID."
                  call allMessage(ERROR, ScratchMessage)
                  write (ScratchMessage, &
                         '(A,F0.9,A,F0.6,A,I0)') &
                     "  X1= ", X1, &
                     " Y1=", Y1, " LOCAL NODE=", NBV(IDX)
                  call allMessage(ERROR, ScratchMessage)
                  write (ScratchMessage, &
                         '(A,F0.9,A,F0.6,A,I0)') &
                     " X2= ", X2, &
                     " Y2=", Y2, " LOCAL NODE=", NBV(IDX2)
                  call allMessage(ERROR, ScratchMessage)
                  call terminate()
               end if
               call ASSIGN_TVW_TIMING(VARYTYPE, IDX, IDX2)
               if (VARYTYPE == 1) then
                  if (BAR_DEG_START(IDX) > &
                      BAR_DEG_END(IDX)) then
                     call allMessage(ERROR, &
                                     "INVALID BARRIER DEGREDATION"// &
                                     "TIME.")
                     call terminate()
                  end if
               elseif (VARYTYPE == 2) then
                  if ((BAR_FAILURE_DURATION(IDX) <= 0d0) &
                      .or. (BAR_FAILURE_DURATION(IDX2) <= &
                            0d0)) then
                     call allMessage(ERROR, &
                                     "FAILURE_DURATION MUST BE "// &
                                     "GREATER THAN ZERO.")
                     call terminate()
                  end if
               elseif (VARYTYPE == 3) then
                  BAR_SCHEDULE(IDX)%OFFSET = OFFSET
                  BAR_SCHEDULE(IDX2)%OFFSET = OFFSET
               end if
            case DEFAULT
               !...........................WE SHOULD NOT BE HERE?
               call allMessage(ERROR, &
                               "INVALID BOUNDARY CONDITION DETECTED.")
               call terminate()
            end select

         case DEFAULT
            call allMessage(ERROR, "INVALID WEIR VARIATION"// &
                            "SPECIFIED.")
            call terminate()
         end select

      end do

      !..FREE THE MEMORY USED FOR THE KDTREE2
      call KDTREE2_DESTROY(BARRIER_SEARCHTREE)
      deallocate (BAR_LOCATIONS)
#ifdef CMPI
      call KDTREE2_DESTROY(GHOST_SEARCHTREE)
      deallocate (GHOST_LOCATIONS)
#endif

#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()
      return
200   call allMessage(ERROR, "PROBLEM READING TIME VARYING WEIR FILE.")
      call terminate()

      !..GENERAL MESSAGES FOR TIME VARYING WEIRS
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

   end subroutine PARSE_TIME_VARYING_WEIR_INFO
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E  P A R S E _ S C H E D U L E
   !-----------------------------------------------------------------------
   !> THIS SUBROUTINE READS A FILE THAT CONTAINS THE NAMELIST STYPE INPUT
   !> FOR A WEIR THAT VARIES BASED UPON A SCHEDULE THAT IS SPECIFIED IN THE
   !> TIME VARYING WEIR INPUT FILE
   !>
   !> @param MYFILE - The file that contains the schedule
   !> @param MYSCHED - The schedule that is read in
   !-----------------------------------------------------------------------
   subroutine PARSE_SCHEDULE(MYFILE, MYSCHED)
      use GLOBAL, only: DTDP, IHOT, ITHS
      implicit none
      character(*), intent(in) :: MYFILE
      type(BARRIER_SCHEDULE_T), intent(out) :: MYSCHED

      real(8), parameter :: eps = epsilon(1d0)

      character(500) :: origLine, modLine
      real(8) :: BAR_DEG_START, BAR_DEG_END
      real(8) :: DELTA, HOTSEC
      integer :: I
      logical :: exists

      namelist /SCHEDULE/ &
         TimeStartDay, TimeStartHour, TimeStartMin, &
         TimeStartSec, TimeEndDay, TimeEndHour, &
         TimeEndMin, TimeEndSec, ZF, DELTA, HOT

      call setMessageSource("PARSE_SCHEDULE")
#if defined(TVW_TRACE)|| defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif

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
                  HOTSEC = DTDP*ITHS
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
                (abs(BAR_DEG_END) <= eps)) then
               call allMessage(ERROR, &
                               "You must specify start and end for "// &
                               "each point in barrier schedule.")
               call terminate()
            end if
            if (ISNULL(ZF)) then
               call allMessage(ERROR, "You must specicy"// &
                               " a flag or final elevation for "// &
                               "BARHT_FINAL")
               call terminate()
            elseif (abs(ZF + 99991d0) <= eps .or. abs(ZF + 99992d0) <= eps) then
               if (ISNULL(DELTA)) then
                  call allMessage(ERROR, "If specifying a " &
                                  //"relative value for ZF, you must " &
                                  //"specify DELTA = ")
                  call terminate()
               end if
            end if
         end do
         close (97, STATUS="DELETE")
      end if
#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()
      return

   end subroutine PARSE_SCHEDULE
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E  A S S I G N _ T V W _ T I M I N G
   !-----------------------------------------------------------------------
   !>  THIS SUBROUTINE ASSIGNS THE TIMING PRESENT IN THE JUST READ IN
   !>  NAMELIST TO THE CORRECT LOCATIONS IN THE BARRIER VARIATION
   !>  ARRAYS
   !>
   !> @param VAR - The type of time varying weir
   !> @param INDEX1 - The index of the first node
   !> @param INDEX2 - The index of the second node
   !-----------------------------------------------------------------------
   subroutine ASSIGN_TVW_TIMING(VAR, INDEX1, INDEX2)
      integer, intent(in) :: VAR
      integer, intent(in) :: INDEX1
      integer, intent(in), optional :: INDEX2
      integer :: I

      call setMessageSource("ASSIGN_TVW_TIMING")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif
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
               call allMessage(ERROR, "Schedule not found.")
               call terminate()
            end if
         end do
      end if

#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()
      return

   end subroutine ASSIGN_TVW_TIMING
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     F U N C T I O N  I S N U L L
   !-----------------------------------------------------------------------
   !> THIS FUNCTION CHECKS TO SEE IF A NAMELIST INPUT VARIABLE WAS MODIFIED
   !> WHEN THE NAMELIST WAS READ, MEANING IT WAS PRESENT IN THE NAMELIST
   !> THAT WAS READ IN. RETURNS TRUE IF THE VARIABLE WAS NOT ALTERED NAD
   !> FALSE IF IT WAS.
   !>
   !> @param R - The real variable to check
   !> @param I - The integer variable to check
   !> @param S - The character variable to check
   !> @return .true. if the variable was not modified, .false. if it was
   !-----------------------------------------------------------------------
   logical function ISNULL(R, I, S)
      implicit none
      real(8), intent(in), optional :: R
      integer, intent(in), optional :: I
      character(*), intent(in), optional :: S
      real(8) :: EPS

      call setMessageSource("ISNULL")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif

      if (present(R)) then
         EPS = epsilon(1.0d0)
         ISNULL = .false.
         if (abs(R + 99999d0) <= EPS) then
            ISNULL = .true.
         end if
      elseif (present(I)) then
         ISNULL = .false.
         if (I == -99999) then
            ISNULL = .true.
         end if
      elseif (present(S)) then
         ISNULL = .false.
         if (trim(adjustl(S)) == "NOFILE") then
            ISNULL = .true.
         end if
      else
         call allMessage(ERROR, "No check specified for null.")
         call terminate()
      end if
#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()

   end function ISNULL
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E  N U L L I F Y _ T V W _ N M L
   !-----------------------------------------------------------------------
   !> THIS SUBROUTINE WILL SET THE NAMELIST VARIABLES TO THEIR NULL VALUES
   !> SO THEY CAN BE CHECKED FOR MODIFICATION LATER
   !-----------------------------------------------------------------------
   subroutine NULLIFY_TVW_NML()
      implicit none

      call setMessageSource("NULLIFY_TVW_NML")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif

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
      ScheduleFile = "NOFILE"

#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()
      return
      !-----------------------------------------------------------------------
   end subroutine NULLIFY_TVW_NML
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E  A L P H A B E T I Z E
   !-----------------------------------------------------------------------
   !> THIS IS AN ALPHABETICAL BUBBLE SORT TO PLACE THE NAMES IN THE INPUT
   !> ARRAY IN ALPHABETICAL ORDER SO THAT IT IS EASIER TO FIND THE TOTAL
   !> NUMBER OF UNIQUE STRINGS THAT HAVE BEEN INPUT TO THE CODE
   !>
   !> @param MYCHAR - The character array to sort
   !-----------------------------------------------------------------------
   subroutine ALPHABETIZE(MYCHAR)
      implicit none
      character(*), intent(inout) :: MYCHAR(:)
      character(1024) :: A, B
      integer :: I, J, K

      call setMessageSource("ALPHABETIZE")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif
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
#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()
      return

   end subroutine ALPHABETIZE
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E  U N I Q U E _ N A M E S
   !-----------------------------------------------------------------------
   !>  THIS SUBROUTINE WILL TAKE AN ALPHABETICAL LIST OF CHARACTER STRINGS
   !>  AS INPUT AND RETURN A LIST WITH ONLY UNIQUE STRINGS AND A COUNT OF
   !>  THE NUMBER OF UNIQUE STRINGS
   !>
   !> @param RAW_LIST - The list of strings to sort
   !> @param NUNIQUE - The number of unique strings
   !> @param UNIQUE_LIST - The list of unique strings
   !-----------------------------------------------------------------------
   subroutine UNIQUE_NAMES(RAW_LIST, NUNIQUE, UNIQUE_LIST)
      implicit none
      character(*), intent(in) :: RAW_LIST(:)
      character(*), intent(out), allocatable :: UNIQUE_LIST(:)
      integer, intent(out) :: NUNIQUE

      character(1024) :: PREV, CURR
      integer :: I, J

      call setMessageSource("UNIQUE_NAMES")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter")
#endif
      if (size(RAW_LIST) < 2) then
         NUNIQUE = 1
         allocate (UNIQUE_LIST(1:size(RAW_LIST)))
         UNIQUE_LIST = RAW_LIST
#if defined(TVW_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG, "Return")
#endif
         call unsetMessageSource()
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

#if defined(TVW_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return")
#endif
      call unsetMessageSource()
      return

   end subroutine UNIQUE_NAMES
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module mod_time_varying_weir_boundary
!-----------------------------------------------------------------------

