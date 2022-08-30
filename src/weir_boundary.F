C-----------------------------------------------------------------------
C WEIR_BOUNDARY.F
C   Written by Zachary Cobell, 2013/01/04
C              The Water Institute
C              zcobell@thewaterinstitute.org
C              (225)421-2761
C
C   CONTAINS THE NEW ROUTINES FOR THE SPECIFICATION OF TIME VARYING WEIRS AND
C   LAND BOUNDARY CONDITIONS. ALL TYPE 3,13,23,4,24,5,25 BOUNDARY CONDITIONS
C   HAVE THEIR ELEVATION VALUES SET USING THIS ROUNTINE. ALSO, WEIR OVERTOPPING
C   IS COMPUTED HERE USING EITHER THE ORIGINAL IMPLEMENTATION (JJW) OR THE
C   STANDARD FORMULATION (SHINTARO BUNYA, CHECKS FOR WET EDGE)
C
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C     M O D U L E  G L O B A L _ W E I R _ V A R I A B L E S 
C-----------------------------------------------------------------------
C     THIS MODULE CONTAINS THE VARIABLES USED THROUGHOUT WEIR 
C     CALCULATIONS AND THE DIFFERENT TIME LEVELS ASSOCIATED WITH
C     A TIME VARYING WEIR.
C-----------------------------------------------------------------------
        MODULE WEIR
      use AMPI_LUN_VIRTUALIZED, only : AMPI_LUN
      use AMPI_LUN_migratable
            USE GLOBAL,ONLY:SZ
            IMPLICIT NONE
            REAL(SZ),ALLOCATABLE     :: BARINHT1(:)    !...Internal Barrier elevation at previous timestep
            REAL(SZ),ALLOCATABLE     :: BARINHT2(:)    !...Internal Barrier elevation at current timestep
            REAL(SZ),ALLOCATABLE     :: BARLANHT1(:)   !...External Barrier elevation at previous timestep
            REAL(SZ),ALLOCATABLE     :: BARLANHT2(:)   !...External Barrier elevation at current timestep
            REAL(SZ), PARAMETER      :: BARMIN=0.04D0

            !...Averaging information in time for boundaries. 
            !   This was previously zeroed out in the code, but the
            !   option is maintained.
#ifdef AVERAGEWEIRFLOW 
            REAL(SZ),ALLOCATABLE     :: RBARWL1AVG(:),RBARWL2AVG(:)
            REAL(SZ),ALLOCATABLE     :: RPIPEWL1AVG(:),RPIPEWL2AVG(:)
            REAL(SZ),PARAMETER       :: BARAVGWT = 0.0D0
            INTEGER                  :: IBSTART 
#endif            
            
            LOGICAL                  :: EXT_TVW = .FALSE. !...If there are external barrier TVW present, avoids checking an unallocated array
            LOGICAL                  :: INT_TVW = .FALSE. !...If there are internal barrier TVW present, avoids checking an unallocated array
            LOGICAL                  :: FOUND_TVW_NML = .FALSE.    !...If the namelist in the fort.15 was found

C-----------------------------------------------------------------------
            CONTAINS
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   A L L O C A T E _ W E I R S
C-----------------------------------------------------------------------
C  THIS ROUTINE WILL ALLOCATE THE NECESSARY VARIABLES AND INITIALIZE
C  THE TIME VARYING ELEVATIONS TO THEIR FORT.14 ELEVATIONS AT TWO TIME
C  LEVELS
C-----------------------------------------------------------------------
            SUBROUTINE ALLOCATE_WEIRS()
                USE GLOBAL,ONLY:allMessage,setMessageSource,
     &              unsetMessageSource,DEBUG
                USE BOUNDARIES,ONLY: NVEL,BARINHT,BARLANHT
                IMPLICIT NONE

                CALL setMessageSource("ALLOCATE_WEIRS")
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Enter")
#endif                

                ALLOCATE(BARINHT1(NVEL))
                ALLOCATE(BARINHT2(NVEL))
                ALLOCATE(BARLANHT1(NVEL))
                ALLOCATE(BARLANHT2(NVEL))

#ifdef AVERAGEWEIRFLOW
                ALLOCATE(RBARWL1AVG(1:NVEL))
                ALLOCATE(RBARWL2AVG(1:NVEL))
                ALLOCATE(RPIPEWL1AVG(1:NVEL))
                ALLOCATE(RPIPEWL2AVG(1:NVEL))
#endif                

                BARINHT1(1:NVEL) = BARINHT(1:NVEL)
                BARINHT2(1:NVEL) = BARINHT(1:NVEL)
                BARLANHT1(1:NVEL) = BARLANHT(1:NVEL)
                BARLANHT2(1:NVEL) = BARLANHT(1:NVEL)

#if defined(WEIR_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Return")
#endif                
                CALL unsetMessageSource()

                RETURN

C-----------------------------------------------------------------------
            END SUBROUTINE ALLOCATE_WEIRS
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
        END MODULE WEIR
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     M O D U L E  T I M E _ V A R Y I N G _ W E I R _ B O U N D A R Y
C-----------------------------------------------------------------------
C  THIS MODULE CONTAINS ALL THE ROUTINES NECESSARY FOR SPECIFICATION OF
C  TIME VARYING WEIR BOUNDARIES. THESE ROUTINES ARE USED TO SETUP AND
C  CALCULATE THE ELEVATIONS THROUGHOUT THE SIMULATION
C-----------------------------------------------------------------------
        MODULE TIME_VARYING_WEIR_BOUNDARY
      use AMPI_LUN_VIRTUALIZED, only : AMPI_LUN
      use AMPI_LUN_migratable
            USE SIZES, ONLY: SZ,MYPROC,LOCALDIR,GLOBALDIR
            USE ADCIRC_MOD, ONLY: ADCIRC_TERMINATE
            USE GLOBAL, ONLY: SCREENUNIT,logMessage,allMessage,
     &          setMessageSource,unsetMessageSource,ERROR,INFO,ECHO,
     &          DEBUG,WARNING,ScratchMessage,ScreenMessage
            USE KDTREE2_MODULE
            USE WEIR

            IMPLICIT NONE

            REAL(SZ),ALLOCATABLE     :: BARHT_FINAL(:)          !...Final elevation for time varying boundary
            REAL(SZ),ALLOCATABLE     :: BAR_LOCATIONS(:,:)      !...Location array for boundary conditions
            REAL(SZ),ALLOCATABLE     :: GHOST_LOCATIONS(:,:)    !...Location array for the ghost nodes
            REAL(SZ),ALLOCATABLE     :: BAR_DEG_START(:)        !...Start time for time varying boundary
            REAL(SZ),ALLOCATABLE     :: BAR_DEG_END(:)          !...End time for time varying boundary
            REAL(SZ),ALLOCATABLE     :: BAR_ETA_MAX(:)          !...Used if the boundary changes at a critical water surface elevation
            REAL(SZ),ALLOCATABLE     :: BAR_FAILURE_START(:)    !...Used to track the first if ETA_MAX has been exceeded
            REAL(SZ),ALLOCATABLE     :: BAR_FAILURE_DURATION(:) !...Amount of time it takes an ETA_MAX barrier to fail to ZF
            LOGICAL,ALLOCATABLE      :: BAR_DEG(:)              !...T/F this node is a time varying boundary
            INTEGER,ALLOCATABLE      :: BAR_VARYTYPE(:)         !...Type of variation to apply to this boundary node
            TYPE(KDTREE2),POINTER    :: BARRIER_SEARCHTREE      !...Search tree for boundary nodes
            TYPE(KDTREE2),POINTER    :: GHOST_SEARCHTREE        !...Search tree for ghost nodes

C...........CHARACTERIZES A SECTION OF THE SCHEDULE
            TYPE BSCHED
                REAL(SZ) :: BAR_DEG_START             !...Time to begin this change
                REAL(SZ) :: BAR_DEG_END               !...Time to end this change
                REAL(SZ) :: BARHT_FINAL               !...Elevation at the end time for this section
                                                      !   Flag Values:
                                                      !     -99990 = Nodal elevation (minimum value)
                                                      !     -99993 = fort.14 weir elevation
                                                      !     -99991 = Add BARHT_DELTA
                                                      !     -99992 = Subtract BARHT_DELTA
                REAL(SZ) :: BARHT_DELTA
            END TYPE

C...........CHARACTERIZES THE SCHEDULE. BASED ON NUMBER OF FILES
            TYPE BARRIER_SCHEDULE_T
                INTEGER                         :: SID        !...Unique ID corresponding to this schedule
                INTEGER                         :: NSECTIONS  !...Number of sections in the schedule
                TYPE(BSCHED),ALLOCATABLE        :: SECTION(:) !...List of all the schedule sections
            END TYPE

C...........SCHEDULE THAT POINTS TO MAIN SET THAT HOUSES THE ACTUAL DATA
C           AND IS BASED UPON THE INDIVIDUAL BOUNDARY NODE
            TYPE BARRIER_SCHEDULE_P
                CHARACTER(1024),POINTER  :: MYFILE        !...File used for this schedule
                INTEGER                  :: SID           !...SID that has a pointer back to SCHEDULE(:)
                INTEGER                  :: MYSEC         !...The current section this node is operating in
                INTEGER                  :: LOOP          !...Should this schedule loop?
                INTEGER                  :: NLOOPS        !...How many times should it loop
                INTEGER                  :: LOOPSCOMPLETE !...How many loops we've completed
                INTEGER,POINTER          :: NSECTIONS     !...How many sections are in the schedule
                REAL(SZ)                 :: PREVEND       !...Previous schedule end
                REAL(SZ)                 :: BARHT_START   !...Starting barrier height
                REAL(SZ)                 :: OFFSET        !...Offset from starting time in schedule
                TYPE(BSCHED),POINTER     :: SECTION(:)    !...A section of the schedule containing a variation
            END TYPE    
            TYPE(BARRIER_SCHEDULE_T),ALLOCATABLE,TARGET   :: SCHEDULE(:)
            TYPE(BARRIER_SCHEDULE_P),ALLOCATABLE     :: BAR_SCHEDULE(:)
            CHARACTER(1024),ALLOCATABLE,TARGET  :: SCHEDULE_LIST(:)
            
C...........THESE VARIABLES ARE PART OF THE NAMELIST THAT IS READ IN
            CHARACTER(200)   :: ScheduleFile
            REAL(SZ) :: X1,X2,Y1,Y2
            REAL(SZ) :: ZF,ETA_MAX
            REAL(SZ) :: TimeStartDay,TimeStartHour,TimeStartMin
            REAL(SZ) :: TimeStartSec
            REAL(SZ) :: TimeEndDay,TimeEndHour,TimeEndMin,TimeEndSec
            REAL(SZ) :: FailureDurationDay,FailureDurationHour
            REAL(SZ) :: FailureDurationMin,FailureDurationSec
            REAL(SZ) :: HOTADD
            INTEGER  :: VaryType,HOT
            INTEGER  :: LOOP,NLOOPS

            INTEGER  :: NSCHEDULES

C...........THE TIME VARYING WEIR NAMELIST
            NAMELIST /TimeVaryingWeir/ 
     &          X1,Y1,X2,Y2,
     &          VaryType,ZF,
     &          ETA_MAX,
     &          TimeStartDay,TimeStartHour,TimeStartMin,TimeStartSec,
     &          TimeEndDay,TimeEndHour,TimeEndMin,TimeEndSec,
     &          FailureDurationDay,FailureDurationHour,
     &          FailureDurationMin,FailureDurationSec,
     &          ScheduleFile,HOT,LOOP,NLOOPS

C...........INTERNAL USE ONLY SUBROUTINES/VARIABLES
            PRIVATE ISNULL
            PRIVATE COMPUTE_BARRIER_HEIGHT_LINEAR
            PRIVATE COMPUTE_BARRIER_HEIGHT_ETAMAX
            PRIVATE COMPUTE_BARRIER_HEIGHT_SCHEDULE
            PRIVATE NULLIFY_TVW_NML
            PRIVATE ALPHABETIZE
            PRIVATE UNIQUE_NAMES
            PRIVATE BARRIER_SCHEDULE_P,BARRIER_SCHEDULE_T
            PRIVATE BSCHED
            PRIVATE BAR_SCHEDULE
            PRIVATE TimeVaryingWeir
            PRIVATE X1,X2,Y1,Y2,ZF,ETA_MAX,VaryType,HOT
            PRIVATE TimeStartDay,TimeStartHour,TimeStartMin,TimeStartSec
            PRIVATE TimeEndDay,TimeEndHour,TimeEndMin,TimeEndSec
            PRIVATE FailureDurationDay,FailureDurationHour
            PRIVATE FailureDurationMin,FailureDurationSec
            PRIVATE LOOP,NLOOPS

C-----------------------------------------------------------------------
            CONTAINS
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   W E I R _ S E T U P
C-----------------------------------------------------------------------
C  THIS ROUTINE IS CALLED FROM adcirc.F TO ARRANGE THE WEIRS ONCE
C  RESNODE HAS BEEN ESTABLISHED
C-----------------------------------------------------------------------
            SUBROUTINE WEIR_SETUP()
                USE BOUNDARIES,ONLY:NFLUXIB,NFLUXIBP,NFLUXB
                IMPLICIT NONE
                
                CALL setMessageSource("WEIR_SETUP")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Enter")
#endif                
                IF((NFLUXIB.EQ.1).OR.(NFLUXIBP.EQ.1).OR.(NFLUXB.EQ.1))THEN
        
                    !...Allocate the new weir arrays at both time levels
                    CALL ALLOCATE_WEIRS()

                    !...Parse the time varying weir file on the local processor
                    IF(found_tvw_nml)THEN
                        CALL PARSE_TIME_VARYING_WEIR_INFO()
                    ENDIF

                ENDIF

#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Return")
#endif                
                CALL unsetMessageSource()
                RETURN

C-----------------------------------------------------------------------
            END SUBROUTINE WEIR_SETUP   
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E  
C       A L L O C A T E _ T I M E V A R Y I N G W E I R S
C-----------------------------------------------------------------------
C  THIS SUBROUTINE SETS UP THE BOUNDARY ARRAYS USED THROUGHOUT THE
C  OTHER ROUTINES FOR PROCESSING. ALLOCATES ARRAYS USED GLOBALLY FOR
C  WEIR TYPE BOUNDARY CONDITIONS AND WEIRS WITH PIPES
C-----------------------------------------------------------------------
            SUBROUTINE ALLOCATE_TIMEVARYINGWEIRS()

                USE MESH, ONLY: X,Y,NP
                USE BOUNDARIES, ONLY: NVEL,LBCODEI,NBV
#ifdef CMPI
                USE MESSENGER,ONLY: RESNODE
#endif                

                IMPLICIT NONE

                INTEGER :: I
                INTEGER :: NWEIR
                INTEGER :: NGHOST
                INTEGER :: IDX

                CALL setMessageSource("ALLOCATE_TIMEVARYINGWEIRS")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Enter")
#endif                

                ALLOCATE(BARHT_FINAL(NVEL))
                ALLOCATE(BAR_DEG_START(NVEL))
                ALLOCATE(BAR_DEG_END(NVEL))
                ALLOCATE(BAR_DEG(NVEL))
                ALLOCATE(BAR_ETA_MAX(NVEL))
                ALLOCATE(BAR_FAILURE_START(NVEL))
                ALLOCATE(BAR_FAILURE_DURATION(NVEL))
                ALLOCATE(BAR_VARYTYPE(NVEL))
                ALLOCATE(BAR_SCHEDULE(NVEL))

                BARHT_FINAL(1:NVEL) = 0D0
                BAR_DEG_START(1:NVEL) = 0D0
                BAR_DEG_END(1:NVEL) = 0D0
                BAR_DEG(1:NVEL) = .FALSE.
                BAR_ETA_MAX(1:NVEL) = 0D0
                BAR_FAILURE_START(1:NVEL) = -1D0
                BAR_FAILURE_DURATION(1:NVEL) = -1D0

C...............COUNT THE NUMBER OF WEIR STYLE BOUNDARIES
                NWEIR = 0
                DO I = 1,NVEL
                    SELECT CASE(LBCODEI(I))
                        CASE(3,13,23,4,24,5,25)
                            NWEIR = NWEIR + 1
                    END SELECT
                ENDDO

C...............CONSTRUCT ARRAYS OF X,Y FOR BOUNDARY NODE LOCATIONS
                ALLOCATE(BAR_LOCATIONS(3,NWEIR))
                IDX = 0
                DO I = 1,NVEL
                    SELECT CASE(LBCODEI(I))
                        CASE(3,13,23,4,24,5,25)
                            IDX = IDX + 1
                            BAR_LOCATIONS(1,IDX) = X(NBV(I))
                            BAR_LOCATIONS(2,IDX) = Y(NBV(I))
                            BAR_LOCATIONS(3,IDX) = DBLE(I)
                    END SELECT
                ENDDO

C...............BUILD SEARCH TREE
                BARRIER_SEARCHTREE => KDTREE2_CREATE(
     &               BAR_LOCATIONS(1:2,:),REARRANGE=.TRUE.,SORT=.TRUE.)

#ifdef CMPI
C...............DO THE SAME PROCESS FOR THE GHOST NODES
                NGHOST=0
                DO I = 1,NP
                    IF(.NOT.RESNODE(I))NGHOST=NGHOST+1
                ENDDO

                ALLOCATE(GHOST_LOCATIONS(3,NGHOST))
                IDX = 0
                DO I = 1,NP
                    IF(.NOT.RESNODE(I))THEN
                        IDX = IDX + 1
                        GHOST_LOCATIONS(1,IDX) = X(I)
                        GHOST_LOCATIONS(2,IDX) = Y(I)
                        GHOST_LOCATIONS(3,IDX) = DBLE(I)
                    ENDIF
                ENDDO

                GHOST_SEARCHTREE => KDTREE2_CREATE(
     &              GHOST_LOCATIONS(1:2,:),REARRANGE=.TRUE.,SORT=.TRUE.)
#endif                


#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Return")
#endif                
                CALL unsetMessageSource()
                RETURN

C-----------------------------------------------------------------------
            END SUBROUTINE ALLOCATE_TIMEVARYINGWEIRS
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E  F I N D _ B O U N D A R Y _ N O D E S
C-----------------------------------------------------------------------
C  THIS SUBROUTINE FINDS THE NODE THAT MATCHS THAT WHICH IS SPECIFIED
C  AS A TIME VARYING BOUNDARY POINT IN THE INPUT FILE. THE NEAREST
C  BOUNDARY NODE IS LOCATED. TYPE 4/24/5/25 BOUNDARIES ARE SPECIFIED WITH BOTH
C  NODES. TYPE 3/13/23 BOUNDARIES ARE SPECIFIED WITH A SINGLE NODE. THIS SCHEME
C  IS USED SO THAT MESH NUMBERING AND BOUNDARY CONDITION ORDERING IN THE FORT.14
C  CANNOT INVALIDATE THE TIME VARYING WEIR INPUT FILE. CHECKS ARE IN PLACE TO ENSURE
C  THE SAME BOUNDARY NODE IS NOT LOCATED TWICE.
C-----------------------------------------------------------------------
            SUBROUTINE FIND_BOUNDARY_NODES(LON,LAT,IDX)
                USE MESH, ONLY: ICS, SLAM0, SFEA0,
     &              DRVSPCOORSROTS, CYLINDERMAP
                USE BOUNDARIES, ONLY: NBV
                USE GLOBAL, ONLY: DEG2RAD, IFSPROTS
#ifdef CMPI                
                USE MESSENGER,ONLY: RESNODE
#endif                
                IMPLICIT NONE
                REAL(SZ),INTENT(IN)  :: LAT
                REAL(SZ),INTENT(IN)  :: LON
                REAL(SZ)             :: LAT_TEMP,LON_TEMP
                REAL(SZ)             :: LATR,LONR
                REAL(SZ)             :: EPS
                REAL(SZ)             :: X,Y
                INTEGER,INTENT(OUT)  :: IDX
                INTEGER,PARAMETER    :: SEARCHDEPTH = 2
                TYPE(KDTREE2_RESULT) :: KDRESULTS(SEARCHDEPTH)
                
                CALL setMessageSource("FIND_BOUNDARY_NODES")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Enter")
#endif                

                EPS = EPSILON(1.0D0)

                IF(ICS.EQ.1)THEN
                    X = LON
                    Y = LAT
                ELSEIF(ICS.NE.1)THEN
                    LAT_TEMP = LAT*DEG2RAD
                    LON_TEMP = LON*DEG2RAD
                    IF ( IFSPROTS.eq. 1) THEN
                       CALL DRVSPCOORSROTS( lonr, latr, LON_TEMP, LAT_TEMP )
                    ELSE
                       latr = LAT_TEMP ;  lonr = LON_TEMP 
                    END IF
                    CALL CYLINDERMAP( X, Y, lonr, latr, SLAM0, SFEA0, ICS)
                ENDIF

                CALL KDTREE2_N_NEAREST(TP=BARRIER_SEARCHTREE,
     &                 QV=(/X,Y/),NN=SEARCHDEPTH,RESULTS=KDRESULTS)

C...............THIS NEEDS SOME WORK, HOWEVER, IM NOT SURE HOW TO HANDLE THE POSSIBILITY
C               THAT A NODE HAS TWO BOUNDARY CONDITIONS BESIDES MAKING IT A FATAL ERROR
                IF(ABS(KDRESULTS(1)%DIS-KDRESULTS(2)%DIS).LE.EPS)THEN
                    CALL allMessage(ERROR,
     &                  "MULTIPLE LOCATIONS FOUND FOR TIME VARYING "//
     &                  "BOUNDARY NODES")
                    CALL ADCIRC_TERMINATE()
                ENDIF

                IF(ABS(KDRESULTS(1)%DIS).GT.EPS)THEN
#ifdef CMPI                
C...................CHECK TO MAKE SURE IT IS NOT INVALID BECAUSE IT
C                   IS A GHOST NODE
                    CALL KDTREE2_N_NEAREST(TP=GHOST_SEARCHTREE,QV=(/X,Y/),
     &                  NN=SEARCHDEPTH,RESULTS=KDRESULTS)
                    IF(ABS(KDRESULTS(1)%DIS).GT.EPS)THEN
C.......................THIS IS NOT A GHOST NODE EITHER. 

                        WRITE(ScratchMessage,'(A,F0.9,A,F0.9,A)') 
     &                      "SPECIFIED NODE LOCATION X=",
     &                      LAT," Y=",LON,
     &                      " NOT FOUND IN LIST OF BOUNDARY NODES."
                        CALL allMessage(ERROR,ScratchMessage)
                        CALL ADCIRC_TERMINATE()
                    ELSE
                        IDX=-1
#if defined(TVW_TRACE) || defined(ALL_TRACE) 
                        CALL allMessage(DEBUG,"Return")
#endif                
                        CALL unsetMessageSource()
                        RETURN
                    ENDIF    

#else
                    WRITE(ScratchMessage,'(A,F0.9,A,F0.9,A)') 
     &                  "SPECIFIED NODE LOCATION X=",
     &                  LAT," Y=",LON,
     &                  " NOT FOUND IN LIST OF BOUNDARY NODES."
                    CALL allMessage(ERROR,ScratchMessage)
                    CALL ADCIRC_TERMINATE()

#endif                    
                ENDIF

                IDX = INT(BAR_LOCATIONS(3,KDRESULTS(1)%IDX))

#if defined(TVW_TRACE) || defined(ALL_TRACE) 
                CALL allMessage(DEBUG,"Return")
#endif                
                CALL unsetMessageSource()

                RETURN

C-----------------------------------------------------------------------
            END SUBROUTINE FIND_BOUNDARY_NODES
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E  C O M P U T E _ B A R R I E R _ H E I G H T
C-----------------------------------------------------------------------
C  THIS SUBROUTINE WILL DETERMINE IF THE BARRIER HEIGHTS ARE SUPPOSED TO CHANGE AND
C  WHEN. IT IS CALLED FOR ALL VARIABLE WEIR TYPE BOUNDARIES AND RETURNS THE INITIAL
C  VALUE FOR THE BOUNDARY IF NO ACTION IS REQUIRED. BOUNDARY IS NOT ALLOWED TO
C  DECREASE BELOW THE ELEVATION OF THE SURROUNDING TOPOGRAPHY. THIS ESSENTIALLY
C  FUNCTIONS AS A WRAPPER FOR ALL THE BOUNDARY VARIATION ROUTINES.
C-----------------------------------------------------------------------
            SUBROUTINE COMPUTE_BARRIER_HEIGHT(MYIDX,TIMELOC,
     &           BAR_HEIGHT_CURRENT,BAR_HEIGHT)

                USE BOUNDARIES,ONLY: NBV
                
                IMPLICIT NONE

                INTEGER,INTENT(IN)   :: MYIDX
                CHARACTER(20)        :: VC
                REAL(SZ),INTENT(IN)  :: TIMELOC
                REAL(SZ),INTENT(IN)  :: BAR_HEIGHT_CURRENT
                REAL(SZ),INTENT(OUT) :: BAR_HEIGHT

                CALL setMessageSource("COMPUTE_BARRIER_HEIGHT")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Enter")
#endif                

C...............SET DEFAULT RETURN VALUE SO WE CAN DUCK OUT AT ANY TIME
                BAR_HEIGHT = BAR_HEIGHT_CURRENT

C...............SELECT THE APPROPRIATE TYPE OF CHANGE AND CALL THE ASSOCIATED
C               ROUTINE. FUTURE DEVELOPMENT CAN TAKE PLACE HERE BY ADDING NEW
C               SUBROUTINES THAT ARE CALLED IN A SIMLIAR WAY. FOR EXAMPLE, 
C               A SUBROUTINE MIGHT BE WRITTEN THAT COMPUTES WAVE FORCES AND 
C               ALTERS THE HEIGHT OF THE BARRIER BASED UPON THAT.
                SELECT CASE(BAR_VARYTYPE(MYIDX))
                    CASE(1)
                        CALL COMPUTE_BARRIER_HEIGHT_LINEAR(MYIDX,
     &                      TIMELOC,BAR_HEIGHT_CURRENT,BAR_HEIGHT)
                    CASE(2)
                        CALL COMPUTE_BARRIER_HEIGHT_ETAMAX(MYIDX,
     &                      TIMELOC,BAR_HEIGHT_CURRENT,BAR_HEIGHT)
                    CASE(3)
                        CALL COMPUTE_BARRIER_HEIGHT_SCHEDULE(MYIDX,
     &                      TIMELOC,BAR_HEIGHT_CURRENT,BAR_HEIGHT)
                    CASE DEFAULT
                        CALL allMessage(ERROR,"Invalid Barrier "
     &                       //"Variation")
                        CALL ADCIRC_TERMINATE()
                END SELECT

C...............WRITE SCREEN/UNIT 16 INFORMATION AT START/END
                SELECT CASE(BAR_VARYTYPE(MYIDX))
                CASE(1,2)
                    IF(BAR_DEG_START(MYIDX).EQ.TIMELOC)THEN
                        IF(BAR_VARYTYPE(MYIDX).EQ.1)THEN
                            VC='LINEAR'
                        ELSEIF(BAR_VARYTYPE(MYIDX).EQ.2)THEN
                            VC='ETA_MAX'
                        ENDIF
#ifdef CMPI
                        WRITE(ScratchMessage,1900) TRIM(VC),'BEGAN',
     &                      TIMELOC,NBV(MYIDX),MYPROC
#else
                        WRITE(ScratchMessage,1901) TRIM(VC),'BEGAN',
     &                      TIMELOC,NBV(MYIDX)
#endif
                        CALL allMessage(INFO,scratchMessage)
                    ENDIF
                    IF(BAR_DEG_END(MYIDX).EQ.TIMELOC)THEN
                        IF(BAR_VARYTYPE(MYIDX).EQ.1)THEN
                            VC='LINEAR'
                        ELSEIF(BAR_VARYTYPE(MYIDX).EQ.2)THEN
                            VC='ETA_MAX'
                        ENDIF
#ifdef CMPI
                        WRITE(ScratchMessage,1900) TRIM(VC),'CONCLUDED',
     &                      TIMELOC,NBV(MYIDX),MYPROC
#else
                        WRITE(ScratchMessage,1901) TRIM(VC),'CONCLUDED',
     &                      TIMELOC,NBV(MYIDX)
#endif
                        CALL allMessage(INFO,ScratchMessage)
                    ENDIF
                END SELECT
 1900           FORMAT('INFO: ',A,' TIME VARYING BOUNDARY ',A,
     &                      ' CHANGING AT TIME = ',E15.8,' AT NODE = ',
     &                      I7,' ON MYPROC = ',I4)
 1901           FORMAT('INFO: ',A,' TIME VARYING BOUNDARY ',A,
     &                      ' CHANGING AT TIME = ',E15.8,' AT NODE = ',
     &                      I7)

#if defined(TVW_TRACE) || defined(ALL_TRACE) 
                CALL allMessage(DEBUG,"Return")
#endif                
                CALL unsetMessageSource()
                RETURN

C-----------------------------------------------------------------------
            END SUBROUTINE COMPUTE_BARRIER_HEIGHT
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E 
C       C O M P U T E _ B A R R I E R _ H E I G H T _ L I N E A R
C-----------------------------------------------------------------------
C  THIS SUBROUTINE COMPUTES THE NEW BARRIER ELEVATION AT THE
C  SPECIFIED BOUNDARY NODE. THE INTERPOLATION FROM STARTING
C  ELEVATION (FORT.14) AND FINAL ELEVATION IS LINEAR
C  OVER THE SPECIFIED TIME AND DOES NOT STOP ONCE
C  IT HAS BEGUN.
C-----------------------------------------------------------------------
            SUBROUTINE COMPUTE_BARRIER_HEIGHT_LINEAR(MYIDX,
     &          TIMELOC,BAR_HEIGHT_CURRENT,BAR_HEIGHT,BARHT_START)

                USE GLOBAL,ONLY: DT,SCREENUNIT
                USE MESH,ONLY: DP
                USE BOUNDARIES, ONLY: BARINHT,BARLANHT,LBCODEI,
     &                                NBV,IBCONN

                IMPLICIT NONE

                INTEGER,INTENT(IN)     :: MYIDX
                REAL(SZ),INTENT(IN)    :: TIMELOC
                REAL(SZ),INTENT(IN)    :: BAR_HEIGHT_CURRENT
                REAL(SZ),INTENT(OUT)   :: BAR_HEIGHT
                REAL(SZ),INTENT(IN),OPTIONAL :: BARHT_START
                REAL(SZ)               :: BAR_DZ,BAR_DT
                REAL(SZ)               :: BAR_START
                REAL(SZ)               :: DEPTH



                CALL setMessageSource("COMPUTE_BARRIER_HEIGHT_LINEAR")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Enter")
#endif                

C...............SET DEFAULT RETURN VALUE
                BAR_HEIGHT = BAR_HEIGHT_CURRENT


C...............GET OUT OF HERE AT THE FIRST OPPORTUNITY
                IF(BAR_VARYTYPE(MYIDX).EQ.1)THEN
                    IF( (BAR_DEG_START(MYIDX).GT. TIMELOC).OR.
     &                  (BAR_DEG_END(MYIDX)  .LT. TIMELOC)    )THEN
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                         CALL allMessage(DEBUG,"Return")
#endif
                         CALL unsetMessageSource()
                         RETURN
                    ENDIF
                ENDIF

C...............COMPUTE THE NEW BOUNDARY HEIGHT
                IF(PRESENT(BARHT_START))THEN
                    !...Schedule style boundaries specify
                    !   their own initial elevation, so dont use
                    !   the elevation in the fort.14 weirs
                    BAR_DZ = BARHT_START - BARHT_FINAL(MYIDX)
                    BAR_START = BARHT_START
                ELSE    
                    SELECT CASE(LBCODEI(MYIDX))
                        CASE(3,13,23)
                            BAR_DZ = BARLANHT(MYIDX) - 
     &                               BARHT_FINAL(MYIDX)
                            BAR_START = BARLANHT(MYIDX)
                        CASE(4,24,5,25)
                            BAR_DZ = BARINHT(MYIDX) - BARHT_FINAL(MYIDX)
                            BAR_START = BARINHT(MYIDX)
                        CASE DEFAULT
                            CALL allMessage(WARNING,
     &                          "INVALID BOUNDARY CONDITION SPECIFIED")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                            CALL allMessage(DEBUG,"Return")
#endif                        
                            CALL unsetMessageSource()
                            RETURN
                    END SELECT
                ENDIF    
                
                IF(BAR_DEG_END(MYIDX).EQ.TIMELOC)THEN

                    BAR_HEIGHT = BARHT_FINAL(MYIDX)
                ELSE
                    BAR_DT = BAR_DEG_END(MYIDX) - BAR_DEG_START(MYIDX)
                    BAR_HEIGHT = BAR_START -
     &                  ( BAR_DZ * ( 1D0 - ( BAR_DEG_END(MYIDX) -
C
     &                  TIMELOC ) / BAR_DT ))
                ENDIF

C...............CHECK AGAINST BATHYMETRY, DO NOT LET WEIR ELEVATION DECREASE BELOW
C               PREVAILING GROUND. CHECK NODE ACROSS BOUNDARY AS WELL AS LOCAL NODE
                DEPTH = MAX(-DP(NBV(MYIDX)),-DP(IBCONN(MYIDX)))
                IF(BAR_HEIGHT.LT.DEPTH)THEN
                    !...WE'RE TOO LOW, SHUT DOWN BARRIER, SET TO MINIMUM
                    IF(BAR_VARYTYPE(MYIDX).EQ.3)THEN
                        BAR_HEIGHT = DEPTH
                    ELSE
                        IF(BARHT_FINAL(MYIDX).GE.DEPTH)THEN
                            !...SET TO USER SPECIFIED MIN
                            BAR_HEIGHT = BARHT_FINAL(MYIDX)
                        ELSE    
#ifdef CMPI
                            WRITE(ScratchMessage,1904) TIMELOC,NBV(MYIDX),
     &                          MYPROC

#else
                            WRITE(ScratchMessage,1905) TIMELOC,NBV(MYIDX)
#endif
                            CALL allMessage(INFO,ScratchMessage)
                            BAR_HEIGHT = DEPTH
                        ENDIF   
                        BAR_DEG(MYIDX) = .FALSE.
                    ENDIF
                ENDIF    

 1904           FORMAT('WARNING: BARRIER CANNOT DECREASE TO ',
     &                 'A VALUE BELOW PREVAILING GROUND.',/,
     &                 '      BARRIER SET TO MINIMUM VALUE AT ',
     &                 'TIME = ',E15.8,' NODE = ',I7,
     &                 ' ON MYPROC = ',I4)

 1905           FORMAT('WARNING: BARRIER CANNOT DECREASE TO ',
     &                 'A VALUE BELOW PREVAILING GROUND.',/,
     &                 '      BARRIER SET TO MINIMUM VALUE AT ',
     &                 ' TIME = ',E15.8,' NODE = ',I7)

#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Return")
#endif                
                CALL unsetMessageSource()
                RETURN

C-----------------------------------------------------------------------
            END SUBROUTINE COMPUTE_BARRIER_HEIGHT_LINEAR
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E
C       C O M P U T E _ B A R R I E R _ H E I G H T _ E T A M A X
C-----------------------------------------------------------------------
C  THIS SUBROUTINE CHECKS TO SEE IF THE BARRIER HAS REACHED
C  ITS ETAMAX FAILURE CRITERIA. IF IT HAS, THE LINEAR BARRIER
C  HEIGHT CHANGE ROUTINE IS CALLED WITH THE VARIABLES CREATED
C  WHEN THIS ROUTINE FIRST DECIDED THAT THE ETAMAX CONDITION
C  WAS MET
C-----------------------------------------------------------------------
            SUBROUTINE COMPUTE_BARRIER_HEIGHT_ETAMAX(MYIDX,TIMELOC,
     &              BAR_HEIGHT_CURRENT,BAR_HEIGHT)

                USE GLOBAL,ONLY: DT,SCREENUNIT,ETA2,NNODECODE
                USE MESH,ONLY: DP
                USE BOUNDARIES, ONLY: BARINHT, BARLANHT, LBCODEI,NBV,IBCONN

                IMPLICIT NONE

                INTEGER,INTENT(IN)     :: MYIDX
                REAL(SZ),INTENT(IN)    :: TIMELOC
                REAL(SZ),INTENT(IN)    :: BAR_HEIGHT_CURRENT
                REAL(SZ),INTENT(OUT)   :: BAR_HEIGHT
                REAL(SZ)               :: BAR_DZ,BAR_DT
                REAL(SZ)               :: BAR_START
                REAL(SZ)               :: DEPTH

                CALL setMessageSource("COMPUTE_BARRIER_HEIGHT_ETAMAX")
#if defined(TVW_TRACE) || defined(ALL_TRACE) 
                CALL allMessage(DEBUG,"Enter")
#endif               

C...............SET DEFAULT RETURN VALUE
                BAR_HEIGHT = BAR_HEIGHT_CURRENT

C...............CHECK IF THE ETA_MAX HAS BEEN EXCEEDED ALREADY
                IF(BAR_FAILURE_START(MYIDX).EQ.-1D0)THEN
                    IF(ETA2(NBV(MYIDX)).GT.BAR_ETA_MAX(MYIDX).AND.NNODECODE(NBV(MYIDX)).EQ.1)THEN
C.......................WE NEED TO INITIALIZE OUR START TIME FOR THE FAILURE
                        BAR_DEG_START(MYIDX) = TIMELOC
                        BAR_DEG_END(MYIDX) = TIMELOC +
     &                      BAR_FAILURE_DURATION(MYIDX)
                        BAR_FAILURE_START(MYIDX) = 1D0
                    ELSEIF(ETA2(IBCONN(MYIDX)).GT.BAR_ETA_MAX(MYIDX).AND.NNODECODE(IBCONN(MYIDX)).EQ.1)THEN
C.......................WE NEED TO INITIALIZE OUR START TIME FOR THE FAILURE
                        BAR_DEG_START(MYIDX) = TIMELOC
                        BAR_DEG_END(MYIDX) = TIMELOC +
     &                      BAR_FAILURE_DURATION(MYIDX)
                        BAR_FAILURE_START(MYIDX) = 1D0
                    ELSE
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                        CALL allMessage(DEBUG,"Return")
#endif                
                        CALL unsetMessageSource()
                        RETURN
                    ENDIF
                ENDIF

C...............IF WE ARE HERE, ETA_MAX WAS PREVIOUSLY EXCEEDED AND THE BARRIER
C               IS DEGRADING AS PRESCRIBED
                CALL COMPUTE_BARRIER_HEIGHT_LINEAR(MYIDX,
     &              TIMELOC,BAR_HEIGHT_CURRENT,BAR_HEIGHT)
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Return")
#endif                
                CALL unsetMessageSource()
                RETURN

C-----------------------------------------------------------------------
            END SUBROUTINE COMPUTE_BARRIER_HEIGHT_ETAMAX
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E  
C       C O M P U T E _ B A R R I E R _ H E I G H T _ S C H E D U L E
C-----------------------------------------------------------------------
C  THIS SUBROUTINE COMPUTES THE VARYING WEIR ELEVATION BASED UPON A 
C  INPUT FILE THAT CONTAINS A SCHEDULE. THIS SCHEDULE CAN BE USED TO 
C  HELP REPRESENT THINGS LIKE GATE OPERATIONS AND OTHER DYNAMIC 
C  PROCESSES IN THE MODEL THAT BEHAVE IN A KNOWN REPETITIVE WAY
C-----------------------------------------------------------------------
            SUBROUTINE COMPUTE_BARRIER_HEIGHT_SCHEDULE(MYIDX,TIMELOC,
     &              BAR_HEIGHT_CURRENT,BAR_HEIGHT)
                USE MESH,ONLY:DP
                USE BOUNDARIES,ONLY:IBCONN,NBV,LBCODEI,BARLANHT,BARINHT
                IMPLICIT NONE
                INTEGER,INTENT(IN)   :: MYIDX
                REAL(SZ),INTENT(IN)  :: TIMELOC
                REAL(SZ),INTENT(IN)  :: BAR_HEIGHT_CURRENT
                REAL(SZ),INTENT(OUT) :: BAR_HEIGHT

                REAL(SZ)             :: PREVEND
                REAL(SZ)             :: BARHT_START
                REAL(SZ)             :: OFFSET
                INTEGER              :: MYSEC
                INTEGER              :: NNBB1,NNBB2

                CALL setMessageSource("COMPUTE_BARRIER_HEIGHT_SCHEDULE")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Enter")
#endif

C...............SET DEFAULT RETURN VALUE
                BAR_HEIGHT = BAR_HEIGHT_CURRENT

C...............Check if we have reached the offset time yet
C               (Time added to beginning of a schedule)
                IF(TIMELOC.LT.BAR_SCHEDULE(MYIDX)%OFFSET)THEN
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                    CALL allMessage(DEBUG,"Return")
#endif
                    CALL unsetMessageSource()
                    RETURN
                ENDIF    

C...............SIMPLIFY THESE VARIABLES
                MYSEC = BAR_SCHEDULE(MYIDX)%MYSEC
                PREVEND = BAR_SCHEDULE(MYIDX)%PREVEND
                OFFSET = BAR_SCHEDULE(MYIDX)%OFFSET

C...............IF THIS IS THE FIRST TIMESTEP IN THIS SECTINO OF THE
C               SCHEDULE, THEN SET UP THE VARIABLES NEEDED
                IF(BAR_SCHEDULE(MYIDX)%SECTION(MYSEC)%BAR_DEG_START.EQ.
     &              TIMELOC-PREVEND-OFFSET)THEN
                    
C...................THIS IS THE BARRIER HEIGHT AT THE START OF THIS
C                   PORTION OF THE SCHEDULE. SAVE IT FOR LINEAR
C                   INTERPOLATION
                    BARHT_START = BARINHT2(MYIDX)
                    BAR_SCHEDULE(MYIDX)%BARHT_START = BARHT_START

C...................GET THE NODE NUMBER
                    NNBB1 = NBV(MYIDX)

C...................INFORM SCREEN AND UNIT 16 WE ARE HERE
#ifdef CMPI     
                    WRITE(ScratchMessage,101)MYSEC,NNBB1,TIMELOC,MYPROC
#else
                    WRITE(ScratchMessage,102)MYSEC,NNBB1,TIMELOC
#endif
                    CALL allMessage(INFO,ScratchMessage)

C...................ALTER THE VARIABLES USED IN COMPUTE_BARRIER_HEIGHT_LINEAR
C                   SO IT CAN CORRECTLY INTERPOLATE NEW BARRIER VALUES
                    BAR_DEG_START(MYIDX) = 
     &                  BAR_SCHEDULE(MYIDX)%SECTION(MYSEC)%BAR_DEG_START
     &                  + PREVEND + OFFSET
                    BAR_DEG_END(MYIDX) = 
     &                  BAR_SCHEDULE(MYIDX)%SECTION(MYSEC)%BAR_DEG_END
     &                  + PREVEND + OFFSET


C...................CHECK THE ZF VALUE FOR FLAGS. THESE FLAGS ARE ADDED
C                   FOR CONVENIENCE TO SET THE BARRIER TO PREDETERMINED
C                   VALUES BASED UPON MESH GEOMETRY.
                    SELECT CASE(INT(BAR_SCHEDULE(MYIDX)%
     &                          SECTION(MYSEC)%BARHT_FINAL))
                        CASE(-99990)
                            !...Set to fort.14 nodal elevation 
                            SELECT CASE(LBCODEI(MYIDX))
                                CASE(3,13,23)
                                    BARHT_FINAL(MYIDX) = -DP(NNBB1)
                                CASE(4,14,24,5,25)
                                    NNBB2 = IBCONN(MYIDX)
                                    BARHT_FINAL(MYIDX) = 
     &                                 MIN(-DP(NNBB1),-DP(NNBB2))
                            END SELECT        
                        CASE(-99991)
                            !...Add BARHT_DELTA to present height
                            SELECT CASE(LBCODEI(MYIDX))
                                CASE(3,13,23)
                                    BARHT_FINAL(MYIDX) = 
     &                                  BARLANHT2(MYIDX) + 
     &                                  BAR_SCHEDULE(MYIDX)%
     &                                  SECTION(MYSEC)%BARHT_DELTA
                                CASE(4,14,24,5,25)
                                    BARHT_FINAL(MYIDX) = 
     &                                  BARINHT2(MYIDX) + 
     &                                  BAR_SCHEDULE(MYIDX)%
     &                                  SECTION(MYSEC)%
     &                                  BARHT_DELTA
                            END SELECT
                        CASE(-99992)
                            !...Subtract BARHT_DELTA from present height
                            SELECT CASE(LBCODEI(MYIDX))
                                CASE(3,13,23)
                                    BARHT_FINAL(MYIDX) = 
     &                                  BARLANHT2(MYIDX) - 
     &                                  BAR_SCHEDULE(MYIDX)%
     &                                  SECTION(MYSEC)%BARHT_DELTA
                                CASE(4,14,24,5,25)
                                    BARHT_FINAL(MYIDX) = 
     &                                  BARINHT2(MYIDX) - 
     &                                  BAR_SCHEDULE(MYIDX)%
     &                                  SECTION(MYSEC)%
     &                                  BARHT_DELTA
                            END SELECT
                        CASE(-99993)
                            !...Set to weir elevation in fort.14
                            SELECT CASE(LBCODEI(MYIDX))
                                CASE(3,13,23)
                                    BARHT_FINAL(MYIDX) = BARLANHT(MYIDX)
                                CASE(4,14,24,5,25)
                                    BARHT_FINAL(MYIDX) = BARINHT(MYIDX)
                            END SELECT
                        CASE DEFAULT
                            BARHT_FINAL(MYIDX) = BAR_SCHEDULE(MYIDX)%
     &                          SECTION(MYSEC)%BARHT_FINAL
                    END SELECT

C...................USE THE NORMAL LINEAR INTERPOLATION ROUTINE TO
C                   GET THE NEW BARRIER ELEVATION
                    CALL COMPUTE_BARRIER_HEIGHT_LINEAR(MYIDX,TIMELOC,
     &                  BAR_HEIGHT_CURRENT,BAR_HEIGHT,BARHT_START)


C...................THIS IS THE END OF A SCHEDULE SECTION, WE NEED TO 
C                   INCREMENT COUNTERS AND PREPARE FOR THE NEXT SECTION
C                   OF THE SCHEDULE
                ELSEIF(BAR_SCHEDULE(MYIDX)%SECTION(MYSEC)%
     &              BAR_DEG_END.EQ.TIMELOC-PREVEND-OFFSET)THEN       

C...................GET NODE NUMBER AND INFORM SCREEN AND UNIT 16
                    NNBB1 = NBV(MYIDX)
#ifdef CMPI     
                    WRITE(ScratchMessage,103)MYSEC,NNBB1,TIMELOC,MYPROC
#else
                    WRITE(ScratchMessage,104)MYSEC,NNBB1,TIMELOC
#endif
                    CALL allMessage(INFO,ScratchMessage)

C...................CALL THE LAST LINEAR INTERPOLATION                 
                    BARHT_START = BAR_SCHEDULE(MYIDX)%BARHT_START
                    CALL COMPUTE_BARRIER_HEIGHT_LINEAR(MYIDX,TIMELOC,
     &                  BAR_HEIGHT_CURRENT,BAR_HEIGHT,BARHT_START)

C...................INCREMENT SCHEDULE SECTION
                    BAR_SCHEDULE(MYIDX)%MYSEC = MYSEC + 1

C...................CHECK TO SEE IF THIS IS THE END OF THE SCHEDULE.
C                   IF THE LOOP FLAG HAS BEEN SET TO -1, WE LOOP THE
C                   SCHEDULE INFINITELY. IF LOOP SET TO 1, WE WILL 
C                   INCREMENT OUR LOOP COUNTER TO CHECK AND SEE IF WE
C                   HAVE DONE THE CORRECT NUMBER OF SCHEDULE LOOPS
                    IF(BAR_SCHEDULE(MYIDX)%MYSEC.GT.
     &                  BAR_SCHEDULE(MYIDX)%NSECTIONS)THEN
                            
                            IF(BAR_SCHEDULE(MYIDX)%LOOP.EQ.1)THEN
                                IF(BAR_SCHEDULE(MYIDX)%LOOPSCOMPLETE.LT.
     &                              BAR_SCHEDULE(MYIDX)%NLOOPS)THEN
                                    BAR_SCHEDULE(MYIDX)%MYSEC = 1
                                    BAR_SCHEDULE(MYIDX)%LOOPSCOMPLETE = 
     &                                  BAR_SCHEDULE(MYIDX)%
     &                                  LOOPSCOMPLETE + 1
                                    BAR_SCHEDULE(MYIDX)%PREVEND = TIMELOC 
     &                                  - OFFSET
                                ELSE
                                    BAR_DEG(MYIDX) = .FALSE.
                                ENDIF
                            ELSEIF(BAR_SCHEDULE(MYIDX)%LOOP.EQ.-1)THEN
                                BAR_SCHEDULE(MYIDX)%MYSEC = 1
                                BAR_SCHEDULE(MYIDX)%PREVEND = TIMELOC 
     &                              - OFFSET
                            ELSE
C...............................THIS NODE SHOULD NO LONGER CHANGE, SO
C                               SET THE APPROPRIATE FLAG
                                BAR_DEG(MYIDX) = .FALSE.
                            ENDIF    
                    ENDIF

C...............THE CASE FOR THE MIDDLE OF THE SCHEDULE. CONTINUE WITH
C               THE INTERPOLATION
                ELSEIF((BAR_SCHEDULE(MYIDX)%SECTION(MYSEC)%BAR_DEG_START
     &                  .LT.TIMELOC-PREVEND-OFFSET).AND.
     &                 (BAR_SCHEDULE(MYIDX)%SECTION(MYSEC)%BAR_DEG_END
     &                  .GT.TIMELOC-PREVEND-OFFSET))THEN
                    BARHT_START = BAR_SCHEDULE(MYIDX)%BARHT_START
                    CALL COMPUTE_BARRIER_HEIGHT_LINEAR(MYIDX,TIMELOC,
     &                  BAR_HEIGHT_CURRENT,BAR_HEIGHT,BARHT_START)
                ENDIF 
                

#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Return")
#endif
                CALL unsetMessageSource()
                RETURN
101             FORMAT("INFO: SCHEDULE STYLE BOUNDARY BEGAN SEGMENT ",
     &                 I4," AT NODE = ",I9," AT TIME = ",E15.8,
     &                 "ON MYPROC = ",I4)
102             FORMAT("INFO: SCHEDULE STYLE BOUNDARY BEGAN SEGMENT ",
     &                 I4," AT NODE = ",I9," AT TIME = ",E15.8)
103             FORMAT("INFO: SCHEDULE STYLE BOUNDARY CONCLUDED",
     &                 " SEGMENT ",I4," AT NODE = ",I9," AT TIME = ",
     &                 E15.8," ON MYPROC = ",I4)
104             FORMAT("INFO: SCHEDULE STYLE BOUNDARY CONDLUDED",
     &                 " SEGMENT ",I4," AT NODE = ",I9," AT TIME = ",
     &                 E15.8)          

C-----------------------------------------------------------------------
            END SUBROUTINE COMPUTE_BARRIER_HEIGHT_SCHEDULE
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E
C       P A R S E _ T I M E _ V A R Y I N G _ W E I R _ I N F O
C-----------------------------------------------------------------------
C  THIS SUBROUTINE READS THE NAMELIST STYLE INPUT FILE FROM THE USER AND
C  SETS UP THE NECESSARY INFORMATION FOR USE DURING THE SIMLUATION
C  AS WELL AS PERFORMING VARIOUS SANITY CHECKS ON THE INPUT TO AVOID 
C  ERRORS DURING THE SIMULATION.
C-----------------------------------------------------------------------
            SUBROUTINE PARSE_TIME_VARYING_WEIR_INFO()

                USE GLOBAL,ONLY: ITHS,DTDP,IHOT,
     &              tvw_file,openFileForRead
                USE SIZES,ONLY: INPUTDIR
                USE BOUNDARIES,ONLY: LBCODEI,NBV,IBCONN

                IMPLICIT NONE
                
                CHARACTER(2000)  :: InputString,modifiedString
                CHARACTER(1024),ALLOCATABLE :: SCHEDULE_LIST_RAW(:)
                REAL(SZ) :: OFFSET
                INTEGER  :: IDX,IDX2
                INTEGER  :: I,IOS
                INTEGER  :: NTIMEVARYINGWEIRS
                INTEGER  :: NSCHEDULES_RAW


                CALL setMessageSource("PARSE_TIME_VARYING_WEIR_INFO")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Enter")
#endif                
                
                CALL openFileForRead(99,TRIM(INPUTDIR)//'/'//
     &              TRIM(tvw_file),IOS)
                IF(IOS.NE.0)THEN
                    NTIMEVARYINGWEIRS = 0
#ifdef CMPI         
                    WRITE(ScratchMessage,101)MYPROC
#else
                    WRITE(ScratchMessage,102)
#endif              
                    CALL allMessage(INFO,ScratchMessage)
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                    CALL allMessage(DEBUG,"Return")
#endif                
                    CALL unsetMessageSource()
                    RETURN
                ENDIF
                READ(AMPI_LUN(99),*) NTIMEVARYINGWEIRS
                IF(NTIMEVARYINGWEIRS.GT.0)THEN
                    CALL ALLOCATE_TIMEVARYINGWEIRS()
#ifdef CMPI
                    WRITE(ScratchMessage,103)MYPROC
#else
                    WRITE(ScratchMessage,104)
#endif
                    CALL allMessage(INFO,ScratchMessage)
                ELSE
#ifdef CMPI     
                    WRITE(ScratchMessage,105)MYPROC
#else               
                    WRITE(ScratchMessage,106)
#endif
                    CALL allMessage(INFO,ScratchMessage)
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                    CALL allMessage(DEBUG,"Return")
#endif                
                    CALL unsetMessageSource()
                    RETURN
                ENDIF

                OPEN(UNIT=AMPI_LUN(98),FILE=TRIM(LOCALDIR)//'/namelist.scratch',
     &            STATUS='REPLACE',ACTION='WRITE')
                call AMPI_LUN_REGISTER(UNIT=AMPI_LUN(98),FILE=TRIM(LOCALDIR)//'/namelist.scratch',
     &            STATUS='REPLACE',ACTION='WRITE')
                    
C...............READ IN UNFORMATTED NAMELIST, WRITE TO FORMATTED NAMELIST
                DO I = 1,NTIMEVARYINGWEIRS
C...................MODIFY LINE FOR NAMELIST FORMATTING
                    READ(AMPI_LUN(99),'(A)') InputString
                    modifiedString = "&TimeVaryingWeir "//
     &                  TRIM(ADJUSTL(InputString))//" /"
                    WRITE(AMPI_LUN(98),'(A)') TRIM(ADJUSTL(modifiedString))
                ENDDO  
                CLOSE(AMPI_LUN(98))

C...............OPEN FORMATTED NAMELIST FILE, READ IN NUMBER OF SCHEDULES
C               WE WILL NEED TO PARSE
                OPEN(UNIT=AMPI_LUN(98),FILE=TRIM(LOCALDIR)//'/namelist.scratch',
     &                 ACTION="READ")
                call AMPI_LUN_REGISTER(UNIT=AMPI_LUN(98),FILE=TRIM(LOCALDIR)//'/namelist.scratch',
     &                 ACTION="READ")
                NSCHEDULES = 0
                NSCHEDULES_RAW = 0
                DO I = 1,NTIMEVARYINGWEIRS
                    CALL NULLIFY_TVW_NML()
                    READ(AMPI_LUN(98),NML=TimeVaryingWeir,IOSTAT=IOS,ERR=200)
                    IF(.NOT.ISNULL(S=SCHEDULEFILE))THEN
                        NSCHEDULES_RAW = NSCHEDULES_RAW + 1
                    ENDIF
                ENDDO
                IF(NSCHEDULES_RAW.GT.0)THEN
                    REWIND(AMPI_LUN(98))
                    IDX = 0
                    ALLOCATE(SCHEDULE_LIST_RAW(1:NSCHEDULES_RAW))
C...................READ BACK IN A LIST CONTAINING THE SCHEDULE FILE NAMES
                    DO I = 1,NTIMEVARYINGWEIRS
                        CALL NULLIFY_TVW_NML()
                        READ(AMPI_LUN(98),NML=TimeVaryingWeir,IOSTAT=IOS)
                        IF(.NOT.ISNULL(S=SCHEDULEFILE))THEN
                            IDX = IDX + 1
                            SCHEDULE_LIST_RAW(IDX) = SCHEDULEFILE
                        ENDIF
                    ENDDO
C...................FIND OUT HOW MANY UNIQUE SCHEDULES THERE ARE
                    CALL ALPHABETIZE(SCHEDULE_LIST_RAW)
                    CALL UNIQUE_NAMES(SCHEDULE_LIST_RAW,NSCHEDULES,
     &                  SCHEDULE_LIST)
                    ALLOCATE(SCHEDULE(1:NSCHEDULES))
                    DO I = 1,NSCHEDULES
                        CALL PARSE_SCHEDULE(SCHEDULE_LIST(I),
     &                      SCHEDULE(I))
                    ENDDO    
                ENDIF    

                CLOSE(AMPI_LUN(98))
C...............OPEN FORMATTED NAMELIST FILE, READ IN 
                OPEN(UNIT=AMPI_LUN(98),FILE=TRIM(LOCALDIR)//'/namelist.scratch',
     &                 ACTION="READ")
                call AMPI_LUN_REGISTER(UNIT=AMPI_LUN(98),FILE=TRIM(LOCALDIR)//'/namelist.scratch',
     &                 ACTION="READ")
                DO I = 1,NTIMEVARYINGWEIRS
C...................READ IN NAMELIST LINE
                    CALL NULLIFY_TVW_NML()
                    READ(AMPI_LUN(98),NML=TimeVaryingWeir,ERR=200,
     &                  IOSTAT=IOS)

                    IF(VARYTYPE.EQ.3)THEN
                        IF(ISNULL(TimeStartDay))TimeStartDay=0D0
                        IF(ISNULL(TimeStartHour))TimeStartHour=0D0
                        IF(ISNULL(TimeStartMin))TimeStartMin=0D0
                        IF(ISNULL(TimeStartSec))TimeStartSec=0D0
                        IF(ISNULL(I=LOOP))LOOP=0
                        IF(ISNULL(I=NLOOPS))NLOOPS=0
                        OFFSET = TimeStartDay*86400D0 + 
     &                           TimeStartHour*3600D0 + 
     &                           TimeStartMin*60D0    + 
     &                           TimeStartSec
                    ELSE
                        OFFSET = 0D0
                    ENDIF


C...................BEGIN SANITY CHECK ON WHAT HAS BEEN SPECIFIED IN NAMELIST
                    SELECT CASE(IHOT)
                    CASE(17,67,68,367,368,567,568)
                        IF(ISNULL(I=HOT).OR.(HOT.EQ.0))THEN
                            HOTADD = 0D0
                        ELSEIF(HOT.EQ.1)THEN
                            HOTADD = DTDP*ITHS
                            OFFSET = OFFSET + HOTADD
                            WRITE(*,*) DTDP,ITHS
                            WRITE(*,*) HOTADD,OFFSET
                        ELSE
                            CALL allMessage(ERROR,
     &                       "INCORRECT"//
     &                       " HOT START VALUE SPECIFIED. HOT=1 OR "//
     &                       "HOT=0 FOR HOT START RELATIVE TIME "//
     &                       "VARYING WEIRS.")
                             CALL ADCIRC_TERMINATE()
                        ENDIF
                    CASE DEFAULT
                        HOTADD = 0D0
                    END SELECT
                    IF(ISNULL(I=VARYTYPE))THEN
                        CALL allMessage(ERROR,
     &                    "YOU MUST "//
     &                    "SPECIFY VARYTYPE= IN THE FORT.15 INPUT "//
     &                    "FILE.")
                        CALL ADCIRC_TERMINATE()
                    ENDIF
                    SELECT CASE(VARYTYPE)
                    CASE(1,2,3)
                        IF(VARYTYPE.EQ.1)THEN
                            IF(ISNULL(X1).OR.ISNULL(Y1).OR.
     &                          ISNULL(ZF))THEN
                                CALL allMessage(ERROR,
     &                              "YOU MUST SPECIFY X1=, Y1=, and "//
     &                              "ZF= ")
                                CALL ADCIRC_TERMINATE()
                            ENDIF    
                            IF(ISNULL(TimeStartday).AND.
     &                          ISNULL(TimeStartHour).AND.
     &                          ISNULL(TimeStartSec))THEN
                                CALL allMessage(ERROR,
     &                            "TIME VARYING BOUNDARY START TIME"//
     &                            " MUST BE SPECIFIED.")
                                CALL ADCIRC_TERMINATE()
                            ENDIF
                            IF(ISNULL(TimeEndday).AND.
     &                          ISNULL(TimeEndHour).AND.
     &                          ISNULL(TimeEndSec))THEN
                                CALL allMessage(ERROR,
     &                              "TIME VARYING BOUNDARY END TIME"//
     &                              " MUST BE SPECIFIED.")
                                CALL ADCIRC_TERMINATE()
                            ENDIF
                            IF(ISNULL(TimeStartDay))TimeStartDay=0D0
                            IF(ISNULL(TimeStartHour))TimeStartHour=0D0
                            IF(ISNULL(TimeStartMin))TimeStartMin=0D0
                            IF(ISNULL(TimeStartSec))TimeStartSec=0D0
                            IF(ISNULL(TimeEndDay))TimeEndDay=0D0
                            IF(ISNULL(TimeEndHour))TimeEndHour=0D0
                            IF(ISNULL(TimeEndMin))TimeEndMin=0D0
                            IF(ISNULL(TimeEndSec))TimeEndSec=0D0
                        ELSEIF(VARYTYPE.EQ.2)THEN
                            IF(ISNULL(X1).OR.ISNULL(Y1).OR.
     &                          ISNULL(ZF))THEN
                                CALL allMessage(ERROR,
     &                              "YOU MUST SPECIFY X1=, Y1=, "//
     &                              "and ZF= ")
                                CALL ADCIRC_TERMINATE()
                            ENDIF    
                            IF((ISNULL(FAILUREDURATIONDAY).AND.
     &                          ISNULL(FAILUREDURATIONHOUR).AND.
     &                          ISNULL(FAILUREDURATIONMIN).AND.
     &                          ISNULL(FAILUREDURATIONSEC)).OR.
     &                          ISNULL(ETA_MAX))THEN
                                CALL allMessage(ERROR,
     &                              " YOU MUST SPECIFY A FAILURE"//
     &                              " DURATION AND MAXIMUM WATER"//
     &                              " SURFACE BEFORE ELEVATION"//
     &                              " CHANGE.")
                                CALL ADCIRC_TERMINATE()
                            ENDIF   
                            IF(ISNULL(FAILUREDURATIONDAY))
     &                          FAILUREDURATIONDAY=0D0
                            IF(ISNULL(FAILUREDURATIONHOUR))
     &                          FAILUREDURATIONHOUR=0D0
                            IF(ISNULL(FAILUREDURATIONMIN))
     &                          FAILUREDURATIONMIN=0D0
                            IF(ISNULL(FAILUREDURATIONSEC))
     &                          FAILUREDURATIONSEC=0D0
                        ELSEIF(VARYTYPE.EQ.3)THEN
                            IF(ISNULL(X1).OR.ISNULL(Y1))THEN
                                CALL allMessage(ERROR,
     &                          "YOU MUST SPECIFY X1= and Y1=")
                                CALL ADCIRC_TERMINATE()
                            ENDIF    
                            IF(ISNULL(S=SCHEDULEFILE))THEN
                                CALL allMessage(ERROR,
     &                              "YOU MUST SPECIFY SCHEDULEFILE= "//
     &                              "FOR VARYTYPE=3.")
                                CALL ADCIRC_TERMINATE()
                            ENDIF
                            IF(ISNULL(TimeStartDay))TimeStartDay=0D0
                            IF(ISNULL(TimeStartHour))TimeStartHour=0D0
                            IF(ISNULL(TimeStartMin))TimeStartMin=0D0
                            IF(ISNULL(TimeStartSec))TimeStartSec=0D0
                            IF(ISNULL(TimeEndDay))TimeEndDay=0D0
                            IF(ISNULL(TimeEndHour))TimeEndHour=0D0
                            IF(ISNULL(TimeEndMin))TimeEndMin=0D0
                            IF(ISNULL(TimeEndSec))TimeEndSec=0D0
                        ENDIF
C.......................END SANITY CHECK

C.......................BEGIN ASSIGNING BOUNDARY CONDITIONS THEIR ATTRIBUTES
                        CALL FIND_BOUNDARY_NODES(X1,Y1,IDX)
                        IF(IDX.EQ.-1)CYCLE
                        SELECT CASE(LBCODEI(IDX))
                        CASE(3,13,23)
C...........................A ONE SIDED STYLE WEIR, WE DONT NEED X2 AND Y2
                            CALL ASSIGN_TVW_TIMING(VARYTYPE,IDX)
                            EXT_TVW = .TRUE.
                            IF(VARYTYPE.EQ.1)THEN
                                IF(BAR_DEG_END(IDX).LE.0D0)THEN
                                   CALL allMessage(ERROR, 
     &                                  "TIME VARYING BOUNDARY END "//
     &                                  "TIME MUST BE GREATER THAN "//
     &                                  "ZERO.")
                                    CALL ADCIRC_TERMINATE()
                                ENDIF
                            ELSEIF(VARYTYPE.EQ.2)THEN
                                IF(BAR_FAILURE_DURATION(IDX).LE.0D0)THEN
                                    CALL allMessage(ERROR, 
     &                                  "FAILURE_DURATION MUST BE "//
     &                                  "GREATER THAN ZERO")
                                    CALL ADCIRC_TERMINATE()
                                ENDIF
                            ELSEIF(VARYTYPE.EQ.3)THEN
                                BAR_SCHEDULE(IDX)%OFFSET = OFFSET   
                            ENDIF

                        CASE(4,24,5,25)
C...........................A TWO SIDED STYLE WEIR, WE NEED X2 AND Y2
                            CALL FIND_BOUNDARY_NODES(X2,Y2,IDX2)
                            IF(IDX2.EQ.-1)THEN
                                CALL allMessage(ERROR,"BOUNDARY CONDITION IMPROPERLY "//
     &                                          "SPLIT INTO GHOST NODE SPACE!")
                                CALL ADCIRC_TERMINATE()
                            ENDIF
                            INT_TVW = .TRUE.

C...........................WE NEED TO CHECK CONNECTIVITY ACCROSS WEIR TO
C                           MAKE SURE THIS PAIR IS VALID
                            IF(NBV(IDX).NE.IBCONN(IDX2))THEN
C...............................THIS IS NOT THE CONNECTIVITY WE EXPECTED, GOOD NIGHT
                                WRITE(ScratchMessage,'(A)') 
     &                              "CONNECTIVITY ACROSS INTERNAL"//
     &                              " TIME VARYING BOUNDARY IS INVALID."
                                CALL allMessage(ERROR,ScratchMessage)
                                WRITE(ScratchMessage,
     &                               '(A,F0.9,A,F0.6,A,I0)')
     &                              "  X1= ",X1,
     &                              " Y1=",Y1," LOCAL NODE=",NBV(IDX)
                                CALL allMessage(ERROR,ScratchMessage)
                                WRITE(ScratchMessage,
     &                              '(A,F0.9,A,F0.6,A,I0)')
     &                              " X2= ",X2,
     &                              " Y2=",Y2," LOCAL NODE=",NBV(IDX2)
                                CALL allMessage(ERROR,ScratchMessage)
                                CALL ADCIRC_TERMINATE()
                            ENDIF
                            CALL ASSIGN_TVW_TIMING(VARYTYPE,IDX,IDX2)
                            IF(VARYTYPE.EQ.1)THEN
                                IF(BAR_DEG_START(IDX).GT.
     &                              BAR_DEG_END(IDX))THEN
                                    CALL allMessage(ERROR, 
     &                                  "INVALID BARRIER DEGREDATION"//
     &                                  "TIME.")
                                    CALL ADCIRC_TERMINATE()
                                ENDIF
                            ELSEIF(VARYTYPE.EQ.2)THEN
                                IF((BAR_FAILURE_DURATION(IDX).LE.0D0)
     &                             .OR.(BAR_FAILURE_DURATION(IDX2).LE.
     &                              0D0))THEN
                                    CALL allMessage(ERROR,
     &                                 "FAILURE_DURATION MUST BE "//
     &                                 "GREATER THAN ZERO.")
                                    CALL ADCIRC_TERMINATE()
                                ENDIF
                            ELSEIF(VARYTYPE.EQ.3)THEN
                                BAR_SCHEDULE(IDX)%OFFSET = OFFSET
                                BAR_SCHEDULE(IDX2)%OFFSET = OFFSET
                            ENDIF
                        CASE DEFAULT
C...........................WE SHOULD NOT BE HERE?
                            CALL allMessage(ERROR,
     &                          "INVALID BOUNDARY CONDITION DETECTED.")
                            CALL ADCIRC_TERMINATE()
                        END SELECT

                    CASE DEFAULT
                        CALL allMessage(ERROR,"INVALID WEIR VARIATION"//
     &                      "SPECIFIED.")
                        CALL ADCIRC_TERMINATE()
                    END SELECT
                
                ENDDO

C...............FREE THE MEMORY USED FOR THE KDTREE2
                CALL KDTREE2_DESTROY(BARRIER_SEARCHTREE)
                DEALLOCATE(BAR_LOCATIONS)
#ifdef CMPI                
                CALL KDTREE2_DESTROY(GHOST_SEARCHTREE)
                DEALLOCATE(GHOST_LOCATIONS)
#endif

#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Return")
#endif                
                CALL unsetMessageSource()
                RETURN
  200           CALL allMessage(ERROR,"PROBLEM READING TIME 
     &              VARYING WEIR FILE.")
                CALL ADCIRC_TERMINATE()

C...............GENERAL MESSAGES FOR TIME VARYING WEIRS
  101           FORMAT("Time varying weir file was not found. All weirs"
     &              ," will be static on MYPROC = ",I6)
  102           FORMAT("Time varying weir file was not found. All weirs"
     &              ," will be static.")
  103           FORMAT("Time varying weir file was found. Time varying",
     &              " weirs have been specified on MYPROC = ",I6)
  104           FORMAT("Time varying weir file was found. Time varying",
     &              " weirs have been specified.")
  105           FORMAT("Time varying weir file was found. No time ",
     &              "varying weirs specified on MYPROC = ",I6)
  106           FORMAT("Time varying weir file was found. No time ",
     &              "varying weirs specified.")

C-----------------------------------------------------------------------
            END SUBROUTINE PARSE_TIME_VARYING_WEIR_INFO
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E  P A R S E _ S C H E D U L E
C-----------------------------------------------------------------------
C  THIS SUBROUTINE READS A FILE THAT CONTAINS THE NAMELIST STYPE INPUT 
C  FOR A WEIR THAT VARIES BASED UPON A SCHEDULE THAT IS SPECIFIED IN THE
C  TIME VARYING WEIR INPUT FILE
C-----------------------------------------------------------------------

            SUBROUTINE PARSE_SCHEDULE(MYFILE,MYSCHED)
                USE GLOBAL,ONLY:DTDP,IHOT,ITHS
                IMPLICIT NONE
                CHARACTER(*),INTENT(IN) :: MYFILE
                TYPE(BARRIER_SCHEDULE_T),INTENT(OUT) :: MYSCHED
                
                CHARACTER(500) :: origLine,modLine
                REAL(SZ) :: BAR_DEG_START,BAR_DEG_END
                REAL(SZ) :: DELTA,HOTSEC
                INTEGER  :: I
                LOGICAL  :: exists

                NAMELIST /SCHEDULE/ 
     &              TimeStartDay,TimeStartHour,TimeStartMin,
     &              TimeStartSec,TimeEndDay,TimeEndHour,
     &              TimeEndMin,TimeEndSec,ZF,DELTA,HOT
                
                CALL setMessageSource("PARSE_SCHEDULE")
#if defined(TVW_TRACE)|| defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Enter")
#endif                

                INQUIRE(FILE=TRIM(GLOBALDIR)//"/"//TRIM(MYFILE),
     &              EXIST=EXISTS)
                IF(EXISTS)THEN
                    OPEN(FILE=TRIM(GLOBALDIR)//"/"//TRIM(MYFILE),
     &                  UNIT=99,ACTION="READ")
                    READ(AMPI_LUN(99),*) MYSCHED%NSECTIONS
                    ALLOCATE(MYSCHED%SECTION(1:MYSCHED%NSECTIONS))
                    OPEN(FILE=TRIM(LOCALDIR)//"/namelist2.scratch",
     &                  UNIT=97,ACTION="WRITE")
                    DO I = 1,MYSCHED%NSECTIONS
                        READ(AMPI_LUN(99),'(A)') OrigLine
                        WRITE(modLine,'(3A)') "&Schedule ",
     &                       TRIM(OrigLine)," /"
                        WRITE(AMPI_LUN(97),'(A)') TRIM(modLine)
                    ENDDO
                    CLOSE(AMPI_LUN(97))
                    CLOSE(AMPI_LUN(99))
                    OPEN(FILE=TRIM(LOCALDIR)//"/namelist2.scratch",
     &                  UNIT=97,ACTION="READ")
                    DO I = 1,MYSCHED%NSECTIONS
                        TimeStartDay  = 0D0
                        TimeStartHour = 0D0
                        TimeStartMin  = 0D0
                        TimeStartSec  = 0D0
                        TimeEndDay    = 0D0
                        TimeEndHour   = 0D0
                        TimeEndMin    = 0D0
                        TimeEndSec    = 0D0
                        HOT           = 0
                        ZF            = -99999D0
                        DELTA         = -99999D0
                        READ(AMPI_LUN(97),NML=SCHEDULE)
                        IF(HOT.EQ.1)THEN
                            SELECT CASE(IHOT)
                                CASE(17,67,68,367,368,567,568)
                                    HOTSEC = DTDP*ITHS
                                CASE DEFAULT
                                    HOTSEC = 0D0
                            END SELECT        
                        ELSE
                            HOTSEC = 0D0
                        ENDIF    
                        BAR_DEG_START = TimeStartDay*86400D0 + 
     &                                  TimeStartHour*3600D0 + 
     &                                  TimeStartMin*60D0    + 
     &                                  TimeStartSec         + 
     &                                  HOTSEC
                        BAR_DEG_END  =  TimeEndDay*86400D0   +
     &                                  TimeEndHour*3600D0   +
     &                                  TimeEndMin*60D0      +
     &                                  TimeStartSec         +
     &                                  HOTSEC
                        MYSCHED%SECTION(I)%BAR_DEG_START = BAR_DEG_START
                        MYSCHED%SECTION(I)%BAR_DEG_END   = BAR_DEG_END
                        MYSCHED%SECTION(I)%BARHT_FINAL   = ZF
                        MYSCHED%SECTION(I)%BARHT_DELTA   = DELTA
                        IF((BAR_DEG_START.LT.0D0).OR.
     &                     (BAR_DEG_END.EQ.0D0))THEN
                            CALL allMessage(ERROR,
     &                         "You must specify start and end for "//
     &                         "each point in barrier schedule.")
                            CALL ADCIRC_Terminate()
                        ENDIF
                        IF(ISNULL(ZF))THEN
                            CALL allMessage(ERROR,"You must specicy"//
     &                      " a flag or final elevation for "//
     &                      "BARHT_FINAL")
                            CALL ADCIRC_Terminate()
                        ELSEIF((ZF.EQ.-99991).OR.
     &                         (ZF.EQ.-99992))THEN
                            IF(ISNULL(DELTA))THEN
                                CALL allMessage(ERROR,"If specifying a "
     &                              //"relative value for ZF, you must "
     &                              //"specify DELTA = ")
                                CALL ADCIRC_Terminate()
                            ENDIF    
                        ENDIF    
                    ENDDO
                    CLOSE(AMPI_LUN(97),STATUS="DELETE")
                 ENDIF   
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Return")
#endif                
                CALL unsetMessageSource()
                RETURN
C-----------------------------------------------------------------------
            END SUBROUTINE PARSE_SCHEDULE    
C-----------------------------------------------------------------------
                
C-----------------------------------------------------------------------
C     S U B R O U T I N E  A S S I G N _ T V W _ T I M I N G
C-----------------------------------------------------------------------
C  THIS SUBROUTINE ASSIGNS THE TIMING PRESENT IN THE JUST READ IN 
C  NAMELIST TO THE CORRECT LOCATIONS IN THE BARRIER VARIATION
C  ARRAYS
C-----------------------------------------------------------------------
            SUBROUTINE ASSIGN_TVW_TIMING(VAR,INDEX1,INDEX2)
                INTEGER,INTENT(IN)          :: VAR
                INTEGER,INTENT(IN)          :: INDEX1
                INTEGER,INTENT(IN),OPTIONAL :: INDEX2
                INTEGER                     :: I
                
                CALL setMessageSource("ASSIGN_TVW_TIMING")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Enter")
#endif                
                IF(VAR.EQ.1)THEN
                    BAR_DEG_START(INDEX1) = 
     &                  TimeStartDay*86400D0   +
     &                  TimeStartHour*3600D0   +
     &                  TimeStartMin*60D0      +
     &                  TimeStartSec           +
     &                  HOTADD
                    BAR_DEG_END(INDEX1) =   
     &                  TimeEndDay*86400D0     +
     &                  TimeEndHour*3600D0     +
     &                  TimeEndMin*60D0        +
     &                  TimeEndSec             +
     &                  HOTADD
                    BAR_DEG(INDEX1) = .TRUE.
                    BAR_VARYTYPE(INDEX1) = VARYTYPE
                    BARHT_FINAL(INDEX1) = ZF
                    IF(PRESENT(INDEX2))THEN
                        BAR_DEG_START(INDEX2) = BAR_DEG_START(INDEX1)
                        BAR_DEG_END(INDEX2) =  BAR_DEG_END(INDEX1) 
                        BAR_DEG(INDEX2) = .TRUE.
                        BAR_VARYTYPE(INDEX2) = VARYTYPE
                        BARHT_FINAL(INDEX2) = ZF
                    ENDIF
                ELSEIF(VAR.EQ.2)THEN
                    BAR_DEG(INDEX1) = .TRUE.
                    BAR_VARYTYPE(INDEX1) = VARYTYPE
                    BAR_ETA_MAX(INDEX1) = ETA_MAX
                    BARHT_FINAL(INDEX1) = ZF
                    BAR_FAILURE_DURATION(INDEX1) =
     &                  FAILUREDURATIONDAY*86400D0 +
     &                  FAILUREDURATIONHOUR*3600D0 +
     &                  FAILUREDURATIONMIN*60D0    +
     &                  FAILUREDURATIONSEC
                    IF(PRESENT(INDEX2))THEN
                        BAR_FAILURE_DURATION(INDEX2) = 
     &                      BAR_FAILURE_DURATION(INDEX1)              
                        BAR_DEG(INDEX2) = .TRUE.
                        BAR_VARYTYPE(INDEX2) = VARYTYPE
                        BAR_ETA_MAX(INDEX2) = ETA_MAX
                        BARHT_FINAL(INDEX2) = ZF
                    ENDIF
                ELSEIF(VAR.EQ.3)THEN
                    BAR_DEG(INDEX1) = .TRUE.
                    BAR_VARYTYPE(INDEX1) = VARYTYPE
                    DO I = 1,NSCHEDULES
                        IF(TRIM(ADJUSTL(SCHEDULEFILE)).EQ.
     &                      TRIM(ADJUSTL(SCHEDULE_LIST(I))))THEN
                            BAR_SCHEDULE(INDEX1)%MYFILE => 
     &                          SCHEDULE_LIST(I)
                            BAR_SCHEDULE(INDEX1)%SID = I
                            BAR_SCHEDULE(INDEX1)%MYSEC = 1
                            BAR_SCHEDULE(INDEX1)%LOOP = LOOP
                            BAR_SCHEDULE(INDEX1)%NLOOPS = NLOOPS
                            BAR_SCHEDULE(INDEX1)%NSECTIONS => 
     &                          SCHEDULE(I)%NSECTIONS
                            BAR_SCHEDULE(INDEX1)%PREVEND = 0D0
                            BAR_SCHEDULE(INDEX1)%SECTION => 
     &                          SCHEDULE(I)%SECTION
                            IF(PRESENT(INDEX2))THEN
                                BAR_DEG(INDEX2) = .TRUE.
                                BAR_VARYTYPE(INDEX2) = VARYTYPE
                                BAR_SCHEDULE(INDEX2)%MYFILE => 
     &                              SCHEDULE_LIST(I)
                                BAR_SCHEDULE(INDEX2)%SID = I
                                BAR_SCHEDULE(INDEX2)%MYSEC = 1
                                BAR_SCHEDULE(INDEX2)%LOOP = LOOP
                                BAR_SCHEDULE(INDEX2)%NLOOPS = NLOOPS
                                BAR_SCHEDULE(INDEX2)%NSECTIONS => 
     &                              SCHEDULE(I)%NSECTIONS
                                BAR_SCHEDULE(INDEX2)%PREVEND = 0D0
                                BAR_SCHEDULE(INDEX2)%SECTION => 
     &                              SCHEDULE(I)%SECTION
                            ENDIF
                            EXIT
                        ENDIF
                        IF(I.EQ.NSCHEDULES)THEN
                            CALL allMessage(ERROR,"Schedule not found.")
                            CALL ADCIRC_Terminate()
                        ENDIF    
                    ENDDO
                ENDIF
                
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Return")
#endif                
                CALL unsetMessageSource()
                RETURN

C-----------------------------------------------------------------------
            END SUBROUTINE ASSIGN_TVW_TIMING
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     F U N C T I O N  I S N U L L
C-----------------------------------------------------------------------
C  THIS FUNCTION CHECKS TO SEE IF A NAMELIST INPUT VARIABLE WAS MODIFIED
C  WHEN THE NAMELIST WAS READ, MEANING IT WAS PRESENT IN THE NAMELIST
C  THAT WAS READ IN. RETURNS TRUE IF THE VARIABLE WAS NOT ALTERED NAD 
C  FALSE IF IT WAS.
C-----------------------------------------------------------------------
            LOGICAL FUNCTION ISNULL(R,I,S)
                IMPLICIT NONE
                REAL(SZ),INTENT(IN),OPTIONAL     :: R
                INTEGER,INTENT(IN),OPTIONAL      :: I
                CHARACTER(*),INTENT(IN),OPTIONAL :: S
                REAL(SZ)            :: EPS
                
                CALL setMessageSource("ISNULL")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Enter")
#endif                

                IF(PRESENT(R))THEN
                    EPS = EPSILON(1.0D0)
                    ISNULL = .FALSE.
                    IF(ABS(R+99999D0).LE.EPS)THEN
                        ISNULL = .TRUE.
                    ENDIF
                ELSEIF(PRESENT(I))THEN
                    ISNULL = .FALSE.
                    IF(I.EQ.-99999)THEN
                        ISNULL = .TRUE.
                    ENDIF
                ELSEIF(PRESENT(S))THEN
                    ISNULL = .FALSE.
                    IF(TRIM(ADJUSTL(S)).EQ."NOFILE")THEN
                        ISNULL = .TRUE.
                    ENDIF    
                ELSE
                    CALL allMessage(ERROR,"No check specified for 
     &                   null.")
                    CALL ADCIRC_TERMINATE()
                ENDIF
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Return")
#endif                
                CALL unsetMessageSource()
C-----------------------------------------------------------------------
            END FUNCTION ISNULL
C-----------------------------------------------------------------------
            

C-----------------------------------------------------------------------
C     S U B R O U T I N E  N U L L I F Y _ T V W _ N M L
C-----------------------------------------------------------------------
C  THIS SUBROUTINE WILL SET THE NAMELIST VARIABLES TO THEIR NULL VALUES
C  SO THEY CAN BE CHECKED FOR MODIFICATION LATER
C-----------------------------------------------------------------------
            SUBROUTINE NULLIFY_TVW_NML()
                IMPLICIT NONE

                CALL setMessageSource("NULLIFY_TVW_NML")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Enter")
#endif                

                X1                         = -99999D0
                X2                         = -99999D0
                Y1                         = -99999D0
                Y2                         = -99999D0
                ZF                         = -99999D0
                HOT                        = -99999
                ETA_MAX                    = -99999D0
                TimeStartDay               = -99999D0
                TimeStartHour              = -99999D0
                TimeStartMin               = -99999D0
                TimeStartSec               = -99999D0
                TimeEndDay                 = -99999D0
                TimeEndHour                = -99999D0
                TimeEndMin                 = -99999D0
                TimeEndSec                 = -99999D0
                FailureDurationDay         = -99999D0
                FailureDurationHour        = -99999D0
                FailureDurationMin         = -99999D0
                FailureDurationSec         = -99999D0
                VARYTYPE                   = -99999
                LOOP                       = -99999
                NLOOPS                     = -99999
                ScheduleFile               = "NOFILE"

#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Return")
#endif                
                CALL unsetMessageSource()
                RETURN
C-----------------------------------------------------------------------
            END SUBROUTINE NULLIFY_TVW_NML
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E  A L P H A B E T I Z E
C-----------------------------------------------------------------------
C  THIS IS AN ALPHABETICAL BUBBLE SORT TO PLACE THE NAMES IN THE INPUT 
C  ARRAY IN ALPHABETICAL ORDER SO THAT IT IS EASIER TO FIND THE TOTAL
C  NUMBER OF UNIQUE STRINGS THAT HAVE BEEN INPUT TO THE CODE
C-----------------------------------------------------------------------
            SUBROUTINE ALPHABETIZE(MYCHAR)
                IMPLICIT NONE
                CHARACTER(*),INTENT(INOUT) :: MYCHAR(:)
                CHARACTER(1024) :: A,B
                INTEGER :: I,J,K                       

                CALL setMessageSource("ALPHABETIZE")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Enter")
#endif                
                DO I = 1,SIZE(MYCHAR)
                    DO J = I+1,SIZE(MYCHAR)
                        A = ADJUSTL(MYCHAR(I))
                        B = ADJUSTL(MYCHAR(J))
                        SORTING: DO K = 1,LEN_TRIM(A)
                            IF(K.GT.LEN_TRIM(B))THEN
                                !...Flip, B should come before A
                                MYCHAR(J) = A
                                MYCHAR(I) = B
                                EXIT SORTING
                            ELSE
                                IF(A(K:K).GT.B(K:K))THEN
                                    !...Flip, B should come before A
                                    MYCHAR(J) = A
                                    MYCHAR(I) = B
                                    EXIT SORTING
                                ELSEIF(A(K:K).LT.B(K:K))THEN
                                    !...No Flip. Order is correct
                                    EXIT SORTING
                                ENDIF
                            ENDIF
                        ENDDO SORTING
                    ENDDO
                ENDDO
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Return")
#endif               
                CALL unsetMessageSource()
                RETURN
C-----------------------------------------------------------------------
            END SUBROUTINE ALPHABETIZE
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E  U N I Q U E _ N A M E S
C-----------------------------------------------------------------------
C  THIS SUBROUTINE WILL TAKE AN ALPHABETICAL LIST OF CHARACTER STRINGS
C  AS INPUT AND RETURN A LIST WITH ONLY UNIQUE STRINGS AND A COUNT OF 
C  THE NUMBER OF UNIQUE STRINGS
C-----------------------------------------------------------------------
            SUBROUTINE UNIQUE_NAMES(RAW_LIST,NUNIQUE,UNIQUE_LIST)
                IMPLICIT NONE
                CHARACTER(*),INTENT(IN)              :: RAW_LIST(:)
                CHARACTER(*),INTENT(OUT),ALLOCATABLE :: UNIQUE_LIST(:)
                INTEGER,INTENT(OUT)                  :: NUNIQUE

                CHARACTER(1024)                      :: PREV,CURR
                INTEGER                              :: I,J

                CALL setMessageSource("UNIQUE_NAMES")
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Enter")
#endif                
                IF(SIZE(RAW_LIST).LT.2)THEN
                    NUNIQUE = 1
                    ALLOCATE(UNIQUE_LIST(1:SIZE(RAW_LIST)))
                    UNIQUE_LIST = RAW_LIST
#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Return")
#endif                
                    CALL unsetMessageSource()
                    RETURN
                ENDIF    

                PREV = RAW_LIST(1)
                NUNIQUE = 1
                DO I = 2,SIZE(RAW_LIST)
                    CURR = RAW_LIST(I)
                    IF(TRIM(PREV).NE.TRIM(CURR))THEN
                        NUNIQUE = NUNIQUE + 1
                        PREV = CURR
                    ENDIF
                ENDDO
                ALLOCATE(UNIQUE_LIST(1:NUNIQUE))
                J = 1
                PREV = RAW_LIST(1)
                UNIQUE_LIST(1) = PREV
                DO I = 2,SIZE(RAW_LIST)
                    CURR = RAW_LIST(I)
                    IF(TRIM(PREV).NE.TRIM(CURR))THEN
                        J = J + 1
                        PREV = CURR
                        UNIQUE_LIST(J) = CURR
                    ENDIF
                ENDDO 

#if defined(TVW_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Return")
#endif                
                CALL unsetMessageSource()
                RETURN    
C-----------------------------------------------------------------------
            END SUBROUTINE UNIQUE_NAMES
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
        END MODULE TIME_VARYING_WEIR_BOUNDARY
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     M O D U L E  W E I R _ F L U X
C-----------------------------------------------------------------------
C   THIS MODULE CONTAINS THE ROUTINES FOR THE WEIR OVERTOPPING FLUX
C   CALCULATIONS FOR INTERNAL AND EXTERNAL WEIR TYPE BOUNDARIES.
C-----------------------------------------------------------------------
        MODULE WEIR_FLUX
      use AMPI_LUN_VIRTUALIZED, only : AMPI_LUN
      use AMPI_LUN_migratable
            USE GLOBAL,ONLY:G,ETA2,RAMPINTFLUX,
     &          setMessageSource,unsetMessageSource,DEBUG,allMessage    
            USE BOUNDARIES,ONLY: LBCODEI,NBV
            USE SIZES,ONLY:SZ
            USE TIME_VARYING_WEIR_BOUNDARY,ONLY:COMPUTE_BARRIER_HEIGHT,
     &              BAR_DEG
            USE WEIR


C-----------------------------------------------------------------------
            CONTAINS
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E  
C       C O M P U T E _ E X T E R N A L _ B O U N A R Y _ F L U X
C-----------------------------------------------------------------------
C  THIS SUBROUTINE COMPUTES THE OVERTOPPING FLUX OUT OF THE DOMAIN FOR 
C  TYPE 3,13,23 BOUNDARY CONDITIONS
C-----------------------------------------------------------------------
            SUBROUTINE COMPUTE_EXTERNAL_BOUNDARY_FLUX(
     &          BARRIER_INDEX,TIMELOC,FLUX)

                USE GLOBAL,ONLY:TVW
                USE BOUNDARIES,ONLY: BARLANCFSP,BARLANHT
                IMPLICIT NONE
                INTEGER,INTENT(IN)     :: BARRIER_INDEX
                INTEGER                :: NNBB
                REAL(SZ),INTENT(IN)    :: TIMELOC
                REAL(SZ),INTENT(OUT)   :: FLUX
                REAL(SZ)               :: RBARWL

                CALL setMessageSource("COMPUTE_EXTERNAL_BOUNDARY_FLUX")
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Enter")
#endif                

                NNBB=NBV(BARRIER_INDEX)
                IF(EXT_TVW)THEN
                    IF(BAR_DEG(BARRIER_INDEX))THEN
                        CALL COMPUTE_BARRIER_HEIGHT(BARRIER_INDEX,
     &                      TIMELOC,BARLANHT1(BARRIER_INDEX),
     &                      BARLANHT2(BARRIER_INDEX))
                    ENDIF
                    TVW(NNBB) = BARLANHT2(BARRIER_INDEX)-
     &                  BARLANHT(BARRIER_INDEX)
                ENDIF    
                RBARWL=2.D0*(ETA2(NNBB) -
     &              BARLANHT2(BARRIER_INDEX))/3.D0
                IF(RBARWL.GT.BARMIN) THEN
                   FLUX = -RampIntFlux
     &                  *BARLANCFSP(BARRIER_INDEX)*RBARWL*
     &                 (RBARWL*G)**0.5D0
                ELSE
                   FLUX = 0.0D0
                ENDIF
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
                CALL allMessage(DEBUG,"Return")
#endif
                CALL unsetMessageSource()
                RETURN

C-----------------------------------------------------------------------
            END SUBROUTINE COMPUTE_EXTERNAL_BOUNDARY_FLUX
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E 
C       C O M P U T E _ I N T E R N A L _ B O U N D A R Y _ F L U X
C-----------------------------------------------------------------------
C     THIS SUBROUTINE CONTAINS THE CALCULATIONS OF OVERTOPPING FLUX
C     PASSED BETWEEN THE THE WEIR PAIR NODES. THIS IS THE VERSION OF
C     THE CODE THAT CHECKS IF THE WEIR HAS A WET EDGE BEFORE PASSING
C     FLUX OVER THE WEIR
C-----------------------------------------------------------------------
            SUBROUTINE COMPUTE_INTERNAL_BOUNDARY_FLUX(
     &          BARRIER_INDEX,BOUNDARYNODE,BOUNDARY,TIMELOC,FLUX)

                USE BOUNDARIES,ONLY:BARINCFSB,BARINCFSP,
     &              NVELL,IBCONN,BARINHT
                USE GLOBAL,ONLY:NIBNODECODE,TVW,USE_TVW,NODECODE
                IMPLICIT NONE
                INTEGER,INTENT(IN),TARGET   :: BARRIER_INDEX
                INTEGER,INTENT(IN),TARGET   :: BOUNDARY
                INTEGER,INTENT(IN),TARGET   :: BOUNDARYNODE
                REAL(SZ),INTENT(IN)         :: TIMELOC
                REAL(SZ),INTENT(OUT)        :: FLUX

                INTEGER                     :: NNBB1,NNBB2
                INTEGER                     :: NNBB1WN,NNBB2WN
                INTEGER                     :: FLOWDIR
                INTEGER,POINTER             :: I,J,K
                REAL(SZ)                    :: RBARWL1,RBARWL2
                REAL(SZ)                    :: RBARWL1F,RBARWL2F

                CALL setMessageSource("COMPUTE_INTERNAL_BOUNDARY_FLUX")
#if defined(WEIR_TRACE) || defined(ALL_TRACE) 
                CALL allMessage(DEBUG,"Enter")
#endif                


C...............SIMPLIFY VARIABLES
                I => BARRIER_INDEX
                J => BOUNDARYNODE
                K => BOUNDARY

C..............CALL THE ORIGINAL IMPLEMENTATION IF THE USER HAS REQUESTED
C              IT VIA THE COMPILER FLAG -DORIGWEIR. IF NOT, PROCEED WITH
C              THE DEFAULT FORMULATION
#ifdef ORIGWEIR   
               CALL COMPUTE_INTERNAL_BOUNDARY_FLUX_ORIG(I,J,K,
     &              TIMELOC,FLUX)
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
               CALL allMessage(DEBUG,"Return")
#endif                
               CALL unsetMessageSource()
               RETURN
#endif

                NNBB1=NBV(I)     ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
                NNBB2=IBCONN(I)  ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
                NNBB1WN = 0      ! COUNT NUMBER OF WET NEIGHBORS
                NNBB2WN = 0      ! COUNT NUMBER OF WET NEIGHBORS
                FLOWDIR = 0      ! DIRECTION OF FLOW

C...............CHECK TO SEE IF THERE IS A WET EDGE
                IF( (J.EQ.1).OR.(J.EQ.NVELL(K)+1) )THEN
                    NNBB1WN = NNBB1WN + NODECODE(NBV(I+1))
                    NNBB2WN = NNBB2WN + NODECODE(IBCONN(I+1))
                ELSE
                    IF( J.EQ.NVELL(K).OR.(J.EQ.NVELL(K)*2) ) THEN
                        NNBB1WN = NNBB1WN + NODECODE(NBV(I-1))
                        NNBB2WN = NNBB2WN + NODECODE(IBCONN(I-1))
                    ELSE
                        NNBB1WN = NNBB1WN + NODECODE(NBV(I-1))
                        NNBB1WN = NNBB1WN + NODECODE(NBV(I+1))
                        NNBB2WN = NNBB2WN + NODECODE(IBCONN(I+1))
                        NNBB2WN = NNBB2WN + NODECODE(IBCONN(I-1))
                    ENDIF
               ENDIF

C..............CHECK TO SEE IF THE BARRIER ELEVATION NEEDS UPDATING
               IF(INT_TVW)THEN
                   IF(BAR_DEG(I))THEN
                       CALL COMPUTE_BARRIER_HEIGHT(I,TIMELOC,BARINHT1(I)
     &                  ,BARINHT2(I))
                   ENDIF
                   TVW(NNBB1) = BARINHT2(I)-BARINHT(I)
               ENDIF

C..............GET WATER LEVEL ABOVE WEIR ON EACH SIDE TO COMPUTE HEAD
C              DEFINE COMPILER FLAG -DAVERAGEWEIRFLOW TO USE BARAVGWT,
C              WHICH GENERALLY IS SET TO ZERO, AND IS HARD CODED
C              TO ZERO HERE. READ_INPUT.F CONTAINS THE SPECIFICATION OF 
C              BARAVGWT. WITH BARAVGWT SET TO ZERO, THERE IS NO NEED FOR
C              IBSTART EITHER.
#if AVERAGEWEIRFLOW
               IF(IBSTART.EQ.0)THEN
                  RBARWL1AVG(I)=ETA2(NNBB1)-BARINHT2(I)
                  RBARWL2AVG(I)=ETA2(NNBB2)-BARINHT2(I)
                  IBSTART=1
               ELSE
                  RBARWL1AVG(I)=(ETA2(NNBB1)-BARINHT2(I)+BARAVGWT
     &                 *RBARWL1)/(1.D0+BARAVGWT)
                  RBARWL2AVG(I)=(ETA2(NNBB2)-BARINHT2(I)+BARAVGWT
     &                 *RBARWL2)/(1.D0+BARAVGWT)
               ENDIF
               RBARWL1=RBARWL1AVG(I)
               RBARWL2=RBARWL2AVG(I)
#else               
C..............ZC - STREAMLINE THE PROCESS SINCE BARAVGWT IS ZERO BY DEFAULT
               RBARWL1=ETA2(NNBB1)-BARINHT2(I)
               RBARWL2=ETA2(NNBB2)-BARINHT2(I)
#endif               
               RBARWL1F=2.D0*RBARWL1/3.D0
               RBARWL2F=2.D0*RBARWL2/3.D0
               FLUX=0.D0
   
               IF((RBARWL1.LT.0.D0).AND.(RBARWL2.LT.0.D0)) THEN
C...............WATER LEVEL ON BOTH SIDES OF BARRIER BELOW BARRIER -> CASE 1
                  FLUX=0.D0
               ELSEIF(ABS(RBARWL1-RBARWL2).LT.0.01D0) THEN
C...............WATER LEVEL EQUAL ON BOTH SIDES OF BARRIER
C................TO WITHIN TOLERANCE BARMIN -> CASE 2
                  FLUX=0.D0
               ELSEIF((RBARWL1.GT.RBARWL2).AND.(RBARWL1.GT.BARMIN)) THEN
C...............WATER LEVEL GREATER ON THIS SIDE OF THE BARRIER AND IS SUCH
C................THAT OVERTOPPING IS OCCURING
C................THUS THIS SIDE IS THE SOURCE SIDE OF THE FLOW ACROSS THE BARRIER
C................NOTE THAT WE DO NOT FORCE THE SOURCE SIDE OF THE BARRIER TO
C................REMAIN WET. ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
C................BARRIER HAS BEEN DRIED. IF IT HAS WE SHUT DOWN THE FLOW ACROSS
C................THE BARRIER
C................ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE BARRIER HAS AT
C................LEAST ONE WET EDGE. IF NOT, WE SHUT DOWN THE FLOW ACROSS
C................THE BARRIER. shintaro v46.28.sb05.05 11/01/2006
                  IF(RBARWL2.GT.RBARWL1F) THEN
C..................OUTWARD SUBCRITICAL FLOW -> CASE 3
                     IF(NODECODE(NNBB1).EQ.0.OR.NNBB1WN.EQ.0) THEN
                        FLUX=0.0D0
                     ELSE
                        FLUX=-RampIntFlux*RBARWL2*BARINCFSB(I)*
     &                     (2.D0*G*(RBARWL1-RBARWL2))**0.5D0
                     ENDIF
                  ELSE
C..................OUTWARD SUPERCRITICAL FLOW -> CASE 4
                     IF(NODECODE(NNBB1).EQ.0.OR.NNBB1WN.EQ.0) THEN
                        FLUX=0.0D0
                     ELSE
                        FLUX=-RampIntFlux*BARINCFSP(I)*RBARWL1F*
     &                    (RBARWL1F*G)**0.5D0
                     ENDIF
                  ENDIF
               ELSEIF((RBARWL2.GT.RBARWL1).AND.(RBARWL2.GT.BARMIN)) THEN
C...............WATER LEVEL LOWER ON THIS SIDE OF BARRIER AND IS SUCH
C................THAT OVERTOPPING IS OCCURING
C................THUS THIS IS THE RECEIVING SIDE OF THE FLOW ACROSS THE BARRIER
C................NOTE THAT WE DO FORCE THE RECEIVING SIDE OF THE BARRIER TO
C................REMAIN WET WHEN THERE IS FLOW ACROSS THE BARRIER.
C................ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
C................BARRIER HAS BEEN DRIED. IF IT HAS, WE SHUT DOWN THE FLOW ACROSS
C................THE BARRIER
C................ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE BARRIER HAS AT
C................LEAST ONE WET EDGE. IF NOT, WE SHUT DOWN THE FLOW ACROSS
C................THE BARRIER. shintaro v46.28.sb05.05 11/01/2006
                  IF(RBARWL1.GT.RBARWL2F) THEN
C..................INWARD SUBCRITICAL FLOW -> CASE 5
                     IF(NODECODE(NNBB2).EQ.0.OR.NNBB2WN.EQ.0) THEN
                        FLUX=0.0D0
                     ELSE
                        FLUX=RampIntFlux*RBARWL1*BARINCFSB(I)*
     &                     (2.0D0*G*(RBARWL2-RBARWL1))**0.5D0
                        NIBNODECODE(NNBB1)=1
                     ENDIF
                  ELSE
C..................INWARD SUPERCRITICAL FLOW -> CASE 6
                     IF(NODECODE(NNBB2).EQ.0.OR.NNBB2WN.EQ.0) THEN
                        FLUX=0.0D0
                     ELSE
                        FLUX=RampIntFlux*BARINCFSP(I)*RBARWL2F*
     &                     (RBARWL2F*G)**0.5D0
                        NIBNODECODE(NNBB1)=1
                     ENDIF
                  ENDIF
               ELSE
                  FLUX=0
               ENDIF
               
#if defined(WEIR_TRACE) || defined(ALL_TRACE) 
               CALL allMessage(DEBUG,"Return")
#endif
               CALL unsetMessageSource()
               RETURN
C-----------------------------------------------------------------------
            END SUBROUTINE COMPUTE_INTERNAL_BOUNDARY_FLUX
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E
C       C O M P U T E _ C R O S S _ B A R R I E R _ P I P E _ F L U X
C-----------------------------------------------------------------------
C  THIS ROUTINE COMPUTE THE DISCHARGE ACROSS TYPE 5,25 BOUNDARY
C  CONDITIONS (CROSS BARRIER PIPES)
C-----------------------------------------------------------------------
            SUBROUTINE COMPUTE_CROSS_BARRIER_PIPE_FLUX(IDX,TIMELOC,FLUX)
                USE BOUNDARIES,ONLY:NBV,IBCONN,PIPEHT,PIPEDIAM,PIPECOEF
                USE GLOBAL,ONLY: G,PI,NODECODE,NIBNODECODE,ETA2,
     &                           RampIntFlux           

                IMPLICIT NONE
                INTEGER,INTENT(IN)     :: IDX
                REAL(SZ),INTENT(IN)    :: TIMELOC
                REAL(SZ),INTENT(OUT)   :: FLUX

                INTEGER                :: NNBB1
                INTEGER                :: NNBB2

                REAL(SZ)               :: RBARWL1
                REAL(SZ)               :: RBARWL2
                
                CALL setMessageSource("COMPUTE_CROSS_BARRIER_PIPE_FLUX")
#if defined(WEIR_TRACE) || defined(ALL_TRACE) 
                CALL allMessage(DEBUG,"Enter")
#endif                

                NNBB1=NBV(IDX)     ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
                NNBB2=IBCONN(IDX)  ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER

#ifdef AVERAGEWEIRFLOW
                IF(IBSTART.EQ.0)THEN
                  RPIPEWL1AVG(IDX)=ETA2(NNBB1)-PIPEHT(IDX)
                  RPIPEWL2AVG(IDX)=ETA2(NNBB2)-PIPEHT(IDX)
                  IBSTART=1
                ELSE
                  RPIPEWL1AVG(IDX)=(ETA2(NNBB1)-PIPEHT(IDX)+BARAVGWT
     &                 *RPIPEWL1AVG(IDX))/(1.D0+BARAVGWT)
                  RPIPEWL2AVG(IDX)=(ETA2(NNBB2)-PIPEHT(IDX)+BARAVGWT
     &                 *RPIPEWL2AVG(IDX))/(1.D0+BARAVGWT)
                ENDIF
                RBARWL1=RPIPEWL1AVG(IDX)
                RBARWL2=RPIPEWL2AVG(IDX)
#else
                RBARWL1=ETA2(NNBB1)-PIPEHT(IDX)
                RBARWL2=ETA2(NNBB2)-PIPEHT(IDX)
#endif                
                IF((RBARWL1.LT.0.D0).AND.(RBARWL2.LT.0.D0)) THEN
C...............WATER LEVEL ON BOTH SIDES OF BARRIER BELOW PIPE -> CASE 1
                  FLUX=0.D0
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
                  CALL allMessage(DEBUG,"Return")
#endif                
                  CALL unsetMessageSource()
                  RETURN
                ENDIF
                IF(ABS(RBARWL1-RBARWL2).LT.BARMIN) THEN
C...............WATER LEVEL EQUAL ON BOTH SIDES OF PIPE -> CASE 2
                  FLUX=0.D0
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
                  CALL allMessage(DEBUG,"Return")
#endif                
                  CALL unsetMessageSource()
                  RETURN 
                ENDIF
                IF((RBARWL1.GT.RBARWL2).AND.(RBARWL1.GT.BARMIN)) THEN
C...............WATER LEVEL GREATER ON THIS SIDE OF THE PIPE AND IS SUCH
C................THAT OUTWARD DISCHARGE IS OCCURING
C................THUS THIS SIDE IS THE SOURCE SIDE OF THE FLOW THROUGH THE PIPE
C................NOTE THAT WE DO NOT FORCE THE SOURCE SIDE OF THE PIPE TO
C................REMAIN WET. ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
C................PIPE HAS BEEN DRIED. IF IT HAS, WE SHUT DOWN THE FLOW ACROSS
C................THE PIPE
                  IF(RBARWL2.LE.0.D0) THEN
C..................OUTWARD FREE DISCHARGE -> CASE 3
                     IF(NODECODE(NNBB1).EQ.0) THEN
                        FLUX=0.0D0
                     ELSE
                        FLUX=-RampIntFlux
     &                       *0.25D0*PI*(PIPEDIAM(IDX))**2
     &                       *(2.D0*G*RBARWL1/(1.D0+PIPECOEF(IDX)))**0.5D0
                     ENDIF
                  ELSE
C..................OUTWARD SUBMERGED DISCHARGE -> CASE 4
                     IF(NODECODE(NNBB1).EQ.0) THEN
                        FLUX=0.0D0
                     ELSE
                        FLUX=-RampIntFlux
     &                       *0.25D0*PI*(PIPEDIAM(IDX))**2
     &                    *(2.D0*G*(RBARWL1-RBARWL2)/PIPECOEF(IDX))**0.5D0
                     ENDIF
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
                     CALL allMessage(DEBUG,"Return")
#endif                
                     CALL unsetMessageSource()
                     RETURN
                  ENDIF
               ENDIF
               IF((RBARWL2.GT.RBARWL1).AND.(RBARWL2.GT.BARMIN)) THEN
C...............WATER LEVEL LOWER ON THIS SIDE OF PIPE AND IS SUCH
C................THAT INWARD DISCHARGE IS OCCURING
C................THUS THIS IS THE RECEIVING SIDE OF THE FLOW THROUGH THE  PIPE
C................NOTE THAT WE DO FORCE THE RECEIVING SIDE OF THE PIPE TO
C................REMAIN WET WHEN THERE IS FLOW ACROSS THE PIPE.
C................ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
C................PIPE HAS BEEN DRIED. IF IT HAS, WE SHUT DOWN THE FLOW THROUGH
C................THE PIPE
                  IF(RBARWL1.LE.0) THEN
C..................INWARD FREE DISCHARGE -> CASE 5
                     IF(NODECODE(NNBB2).EQ.0) THEN
                        FLUX=0.0D0
                     ELSE
                        FLUX=RampIntFlux
     &                       *0.25D0*PI*(PIPEDIAM(IDX))**2
     &                       *(2.D0*G*RBARWL2/(1.D0+PIPECOEF(IDX)))**0.5D0
                        NIBNODECODE(NNBB1)=1
                     ENDIF
                  ELSE
C..................INWARD SUBMERGED DISCHARGE -> CASE 6
                     IF(NODECODE(NNBB2).EQ.0) THEN
                        FLUX=0.0D0
                     ELSE
                        FLUX=RampIntFlux
     &                       *0.25D0*PI*(PIPEDIAM(IDX))**2
     &                   *(2.D0*G*(RBARWL2-RBARWL1)/PIPECOEF(IDX))**0.5D0
                        NIBNODECODE(NNBB1)=1
                     ENDIF
                  ENDIF
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
                  CALL allMessage(DEBUG,"Return")
#endif                
                  CALL unsetMessageSource()
                  RETURN
               ENDIF

#if defined(WEIR_TRACE) || defined(ALL_TRACE) 
               CALL allMessage(DEBUG,"Return")
#endif
               CALL unsetMessageSource()

            END SUBROUTINE


#ifdef ORIGWEIR
C-----------------------------------------------------------------------
C     S U B R O U T I N E 
C       C O M P U T E _ I N T E R N A L _ B O U N D A R Y _ F L U X _ O R I G
C-----------------------------------------------------------------------
C  THIS IS THE ORIGINAL ADCIRC OVERTOPPING ROUTINE THAT ONLY CHECKS TO 
C  SEE IF THE NODE IS WET BEFORE PASSING FLUX ACROSS THE WEIR
C-----------------------------------------------------------------------
            SUBROUTINE COMPUTE_INTERNAL_BOUNDARY_FLUX_ORIG(
     &          BARRIER_INDEX,BOUNDARYNODE,BOUNDARY,TIMELOC,FLUX)


                USE BOUNDARIES,ONLY:BARINCFSB,BARINCFSP,NIBNODECODE,
     &              NODECODE,NVELL,NODECODE,IBCONN
                IMPLICIT NONE
                INTEGER,INTENT(IN),TARGET   :: BARRIER_INDEX
                INTEGER,INTENT(IN),TARGET   :: BOUNDARY
                INTEGER,INTENT(IN),TARGET   :: BOUNDARYNODE
                REAL(SZ),INTENT(IN)         :: TIMELOC
                REAL(SZ),INTENT(OUT)        :: FLUX

                INTEGER                     :: NNBB1,NNBB2
                INTEGER                     :: NNBB1WN,NNBB2WN
                INTEGER                     :: FLOWDIR
                INTEGER,POINTER             :: I,J,K
                REAL(SZ)                    :: RBARWL1,RBARWL2
                REAL(SZ)                    :: RBARWL1F,RBARWL2F
                
                
                CALL setMessageSource("COMPUTE_INTERNAL_BOUNDARY_FLUX")
#if defined(WEIR_TRACE) || defined(ALL_TRACE) 
                CALL allMessage(DEBUG,"Enter")
#endif                


C...............SIMPLIFY VARIABLES
                I => BARRIER_INDEX
                J => BOUNDARYNODE
                K => BOUNDARY

                NNBB1=NBV(I)     ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
                NNBB2=IBCONN(I)  ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
                IF(INT_TVW)THEN
                    IF(BAR_DEG(I))THEN
                        CALL COMPUTE_BARRIER_HEIGHT(I,TIMELOC,
     &                      BARINHT1(I),BARINHT2(I))
                    ENDIF
                ENDIF
                IF(IBSTART.EQ.0)THEN
                  RBARWL1AVG(I)=ETA2(NNBB1)-BARINHT2(I)
                  RBARWL2AVG(I)=ETA2(NNBB2)-BARINHT2(I)
                  IBSTART=1
                ELSE
                  RBARWL1AVG(I)=(ETA2(NNBB1)-BARINHT2(I)+BARAVGWT
     &                 *RBARWL1AVG(I))/(1.D0+BARAVGWT)
                  RBARWL2AVG(I)=(ETA2(NNBB2)-BARINHT2(I)+BARAVGWT
     &                 *RBARWL2AVG(I))/(1.D0+BARAVGWT)
                ENDIF
                RBARWL1=RBARWL1AVG(I)
                RBARWL2=RBARWL2AVG(I)
                RBARWL1F=2.D0*RBARWL1/3.D0
                RBARWL2F=2.D0*RBARWL2/3.D0
                FLUX=0.D0
                IF((RBARWL1.LT.0.D0).AND.(RBARWL2.LT.0.D0)) THEN
C...............WATER LEVEL ON BOTH SIDES OF BARRIER BELOW BARRIER -> CASE 1
                  FLUX=0.D0
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
                  CALL allMessage(DEBUG,"Return")
#endif                
                  CALL unsetMessageSource()
                  RETURN
                ENDIF
                IF(ABS(RBARWL1-RBARWL2).LT.0.01D0) THEN
C...............WATER LEVEL EQUAL ON BOTH SIDES OF BARRIER
C................TO WITHIN TOLERANCE BARMIN -> CASE 2
                  FLUX=0.D0
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
                  CALL allMessage(DEBUG,"Return")
#endif                
                  CALL unsetMessageSource()
                  RETURN
                ENDIF
                IF((RBARWL1.GT.RBARWL2).AND.(RBARWL1.GT.BARMIN)) THEN
C...............WATER LEVEL GREATER ON THIS SIDE OF THE BARRIER AND IS SUCH
C................THAT OVERTOPPING IS OCCURING
C................THUS THIS SIDE IS THE SOURCE SIDE OF THE FLOW ACROSS THE BARRIER
C................NOTE THAT WE DO NOT FORCE THE SOURCE SIDE OF THE BARRIER TO
C................REMAIN WET. ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
C................BARRIER HAS BEEN DRIED. IF IT HAS WE SHUT DOWN THE FLOW ACROSS
C................THE BARRIER
                  IF(RBARWL2.GT.RBARWL1F) THEN
C..................OUTWARD SUBCRITICAL FLOW -> CASE 3
                     IF(NODECODE(NNBB1).EQ.0) THEN
                        FLUX=0.0D0
                     ELSE
                        FLUX=-RampIntFlux*RBARWL2*BARINCFSB(I)*
     &                     (2.D0*G*(RBARWL1-RBARWL2))**0.5D0
                     ENDIF
                  ELSE
C..................OUTWARD SUPERCRITICAL FLOW -> CASE 4
                     IF(NODECODE(NNBB1).EQ.0) THEN
                        FLUX=0.0D0
                     ELSE
                        FLUX=-RampIntFlux*BARINCFSP(I)*RBARWL1F*
     &                    (RBARWL1F*G)**0.5D0
                     ENDIF
                  ENDIF
#if defined(WEIR_TRACE) || defined(ALL_TRACE)
                  CALL allMessage(DEBUG,"Return")
#endif                
                  CALL unsetMessageSource()
                  RETURN
               ENDIF
               IF((RBARWL2.GT.RBARWL1).AND.(RBARWL2.GT.BARMIN)) THEN
C...............WATER LEVEL LOWER ON THIS SIDE OF BARRIER AND IS SUCH
C................THAT OVERTOPPING IS OCCURING
C................THUS THIS IS THE RECEIVING SIDE OF THE FLOW ACROSS THE BARRIER
C................NOTE THAT WE DO FORCE THE RECEIVING SIDE OF THE BARRIER TO
C................REMAIN WET WHEN THERE IS FLOW ACROSS THE BARRIER.
C................ALSO WE CHECK TO SEE IF THE SOURCE SIDE OF THE
C................BARRIER HAS BEEN DRIED. IF IT HAS, WE SHUT DOWN THE FLOW ACROSS
C................THE BARRIER
                  IF(RBARWL1.GT.RBARWL2F) THEN
C..................INWARD SUBCRITICAL FLOW -> CASE 5
                     IF(NODECODE(NNBB2).EQ.0) THEN
                        FLUX=0.0D0
                     ELSE
                        FLUX=RampIntFlux*RBARWL1*BARINCFSB(I)*
     &                     (2.0D0*G*(RBARWL2-RBARWL1))**0.5D0
                        NIBNODECODE(NNBB1)=1
                     ENDIF
                  ELSE
C..................INWARD SUPERCRITICAL FLOW -> CASE 6
                     IF(NODECODE(NNBB2).EQ.0) THEN
                        FLUX=0.0D0
                     ELSE
                        FLUX=RampIntFlux*BARINCFSP(I)*RBARWL2F*
     &                     (RBARWL2F*G)**0.5D0
                        NIBNODECODE(NNBB1)=1
                     ENDIF
                  ENDIF
               ENDIF

#if defined(WEIR_TRACE) || defined(ALL_TRACE) 
               CALL allMessage(DEBUG,"Return")
#endif
               CALL unsetMessageSource()

               RETURN


C-----------------------------------------------------------------------
            END SUBROUTINE COMPUTE_INTERNAL_BOUNDARY_FLUX_ORIG
C-----------------------------------------------------------------------
#endif

C-----------------------------------------------------------------------
        END MODULE
C-----------------------------------------------------------------------

