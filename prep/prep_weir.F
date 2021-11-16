C******************************************************************************
C MODULE WEIR_BOUNDARY
C   Written by Zachary Cobell, 2012/03/12
C
C   CONTAINS THE NEW ROUTINES FOR THE SPECIFICATION OF TIME VARYING WEIRS AND
C   LAND BOUNDARY CONDITIONS. ALL TYPE 3,13,23,4,24,5,25 BOUNDARY CONDITIONS
C   HAVE THEIR VALUES SET USING THIS ROUNTINE.
C
C******************************************************************************
C
C
        MODULE PREP_WEIR    
            USE PRE_GLOBAL,ONLY: NVELL,IBCONNR,LBCODE,NBOU,MNVEL,NBVV,X,
     &                           Y,NVEL,NPROC,NNODG 
            USE KDTREE2_MODULE
            USE PRESIZES,ONLY : SZ
            USE GLOBAL,ONLY: SCREENUNIT

            IMPLICIT NONE
            
            REAL(SZ),ALLOCATABLE       :: BAR_LOCATIONS(:,:)
            INTEGER,ALLOCATABLE        :: LNBV(:)
            INTEGER,ALLOCATABLE        :: LIBCONN(:)
            INTEGER,ALLOCATABLE        :: LBCODEI(:)
            INTEGER,ALLOCATABLE        :: NODES_TVW(:)
            INTEGER,ALLOCATABLE        :: NWEIRBNDRY(:)
            INTEGER                    :: NTIMEVARYINGWEIR
            TYPE(KDTREE2),POINTER      :: BARRIER_SEARCHTREE
            TYPE(KDTREE2),POINTER      :: GLOBAL_SEARCHTREE
            CHARACTER(400),ALLOCATABLE :: TIMEVARYINGWEIRMSSG(:)

            CONTAINS
            
            SUBROUTINE ALLOCATE_WEIRBOUNDARY()
                
                IMPLICIT NONE
                
                INTEGER :: I
                INTEGER :: JG
                INTEGER :: K
                INTEGER :: NWEIR
                INTEGER :: IDX,IDX2,IDX3
                REAL(SZ),ALLOCATABLE :: XY(:,:)

                ALLOCATE(XY(1:2,1:NNODG))
                ALLOCATE(NWEIRBNDRY(NPROC))
                NWEIRBNDRY(:) = 0

                XY(1,:) = X(:)
                XY(2,:) = Y(:)
                
                GLOBAL_SEARCHTREE => KDTREE2_CREATE(XY,REARRANGE=.TRUE.,
     &               SORT=.TRUE.)

            END SUBROUTINE ALLOCATE_WEIRBOUNDARY
            
            SUBROUTINE FIND_BOUNDARY_NODES(LAT,LON,IDX)
                IMPLICIT NONE
                REAL(SZ),INTENT(IN)  :: LAT
                REAL(SZ),INTENT(IN)  :: LON
                REAL(SZ)             :: EPS
                REAL(SZ)             :: X,Y
                INTEGER,INTENT(OUT)  :: IDX
                INTEGER,PARAMETER    :: SEARCHDEPTH = 1
                TYPE(KDTREE2_RESULT) :: KDRESULTS(SEARCHDEPTH)

                EPS = EPSILON(1.0D0)
                
                X = LAT
                Y = LON

                CALL KDTREE2_N_NEAREST(TP=GLOBAL_SEARCHTREE,
     &                 QV=(/X,Y/),NN=SEARCHDEPTH,RESULTS=KDRESULTS)

                IF(KDRESULTS(1)%DIS.GT.EPS)THEN
                    WRITE(*,'(A,F0.9,A,F0.9)') "FATAL ERROR: "//
     &                  "GLOBAL NODE NOT FOUND FOR "//
     &                  "X=",X," Y=",Y
                    CALL EXIT(1)
                ENDIF    

                IDX = KDRESULTS(1)%IDX 

                RETURN

            END SUBROUTINE FIND_BOUNDARY_NODES

            SUBROUTINE PARSE_TIME_VARYING_WEIR_INFO()
                USE GLOBAL,ONLY: USE_TVW,TVW_FILE
                IMPLICIT NONE
                CHARACTER(2000)         :: InputString
                CHARACTER(2000)         :: modifiedString
                CHARACTER(200)          :: ScheduleFile
                REAL(SZ) :: X1,X2,Y1,Y2
                REAL(SZ) :: ETA_MAX,TimeStartDay,TimeStartSec
                REAL(SZ) :: TimeStartHour,TimeEndDay,TimeEndHour
                REAL(SZ) :: TimeEndSec,ZF
                INTEGER  :: VARYTYPE
                INTEGER  :: IDX
                INTEGER  :: IOS
                INTEGER  :: I
                LOGICAL  :: exists
                
                NAMELIST /TimeVaryingWeir/ X1,Y1,X2,Y2,VaryType,ZF,
     &                                     ETA_MAX,TimeStartDay,
     &                                     TimeStartHour,TimeStartSec,
     &                                     TimeEndDay,TimeEndHour,
     &                                     TimeEndSec

                CALL ALLOCATE_WEIRBOUNDARY()
C...............OPEN THE TIME VARYING WEIR FILE IF FOUND
                INQUIRE(FILE=TRIM(TVW_FILE),EXIST=exists)
                IF(.NOT.exists)THEN
                    WRITE(*,'(A)') "INFO: No time varying weir file "//
     &                  "was found. Not processed."
                    USE_TVW = .FALSE.
                    TVW_FILE= 'none'
                    RETURN
                ENDIF    
                OPEN(FILE=TRIM(TVW_FILE),UNIT=98,ACTION="READ")        
                READ(98,*) NTIMEVARYINGWEIR
                ALLOCATE(TIMEVARYINGWEIRMSSG(NTIMEVARYINGWEIR))
                ALLOCATE(NODES_TVW(NTIMEVARYINGWEIR))
                NODES_TVW(:) = -1
C...............WRITE NAMELIST INPUT TO SCRATCH FILE
                OPEN(UNIT=99,FILE="namelist.scratch",ACTION="WRITE")
                DO I = 1,NTIMEVARYINGWEIR
                    READ(98,'(A)') InputString
                    TIMEVARYINGWEIRMSSG(I) = ADJUSTL(InputString)
                    WRITE(modifiedString,'(A)') "&TimeVaryingWeir "//
     &                  TRIM(ADJUSTL(InputString))//" /"
                    WRITE(99,'(A)') TRIM(modifiedString)
                ENDDO
                CLOSE(99)

                OPEN(UNIT=99,FILE="namelist.scratch",ACTION="READ",
     &               IOSTAT=IOS,ERR=200)
                DO I = 1,NTIMEVARYINGWEIR
                    X1 = -99999D0
                    Y1 = -99999D0
                    READ(99,NML=TimeVaryingWeir,IOSTAT=IOS)
                    CALL FIND_BOUNDARY_NODES(X1,Y1,IDX)
                    NODES_TVW(I) = IDX
                ENDDO    
                CLOSE(99,STATUS="DELETE")
                RETURN
200             WRITE(SCREENUNIT,'(A)') "FATAL ERROR: ERROR WHILE "//
     &              "READING TIME VARYING BOUNDARY INFORMATION!"
                CALL EXIT(1)

            
            END SUBROUTINE PARSE_TIME_VARYING_WEIR_INFO

            LOGICAL FUNCTION ISNULL(CHECKNUMBER)
                IMPLICIT NONE
                REAL(SZ),INTENT(IN) :: CHECKNUMBER
                REAL(SZ)            :: EPS
                EPS = EPSILON(1.0D0)
                ISNULL = .FALSE.
                IF(ABS(CHECKNUMBER+99999D0).LE.EPS)THEN
                    ISNULL = .TRUE.
                ENDIF
            END FUNCTION ISNULL    

        END MODULE PREP_WEIR

            
            
            
            
