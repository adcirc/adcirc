
module wind_drag_module

  ! jgf50.32: The wind drag limit was originally set to 0.002d0
  ! in v46 ... the higher value of 0.0035d0 was used in ~v48 and
  ! is now the default.
  REAL(8) :: WindDragLimit = 0.0035d0   ! maximum wind drag value

! Casey 110518: Added for Mark Powell's sectorial wind drag.
  REAL(8)  :: EyeLat(3) = 0.D0
  REAL(8)  :: EyeLon(3) = 0.D0
  LOGICAL  :: FoundEye
  REAL(8),allocatable :: PWeight(:,:)
  LOGICAL,allocatable :: foundSector(:)

  ! jgf50.32: If DragLawString is not set in the fort.15 input file,
  ! it will remain 'default'; this will cause 'Garratt' to be used
  ! in FUNCTION WindDrag, and 'IceCube' to be used in
  ! FUNCTION WindIceDrag.
  INTEGER,PARAMETER :: DRAGLAW_GARRATT = 99990
  INTEGER,PARAMETER :: DRAGLAW_POWELL  = 99991
  INTEGER,PARAMETER :: DRAGLAW_SWELL   = 99992
#ifndef POWELL
  CHARACTER(len=10) :: DragLawString = 'default' ! controls wind drag formulation
  INTEGER :: dragLawType = DRAGLAW_GARRATT
#else
      ! jgf50.32 Provide backward compatibility with hardcoded POWELL
      ! support.
#warning The flag -DPOWELL is deprecated and will be removed in a future release. Please use the &metcontrol name list functionality to avoid issues

      CHARACTER(len=10) :: DragLawString = 'Powell'
      INTEGER :: dragLawType = DRAGLAW_POWELL
#endif
  logical :: garatt, swell, powell

  !...Function interface for wind drag
  !   Used so that the selection of wind drag
  !   is only performed once. Each wind drag
  !   function has the same function signiture
  ABSTRACT INTERFACE
    FUNCTION p_computeWindDrag(WindSpeed,NodeNumber) RESULT(Drag)
      REAL(8),INTENT(IN)  :: WindSpeed
      INTEGER,INTENT(IN)  :: NodeNumber
      REAL(8)             :: Drag
    END FUNCTION p_computeWindDrag
  END INTERFACE
  PROCEDURE(p_computeWindDrag),POINTER :: f_computeWindDrag => null()

  contains

    subroutine initialize_wind_drag_module()

      use global, only: nws
      use mesh, only: np
#ifdef CMPI
      use messenger, only: msg_fini
#endif

      implicit none

      CALL checkWindDragType()
      CALL mapWindDragFunctionPointer()

      if (NWS==0) return

      ! Make logicals here if NWS is not 0
      garatt = .false.; swell = .false.; powell = .false.
      SELECT CASE(TRIM(DragLawString))
      CASE("Garratt","GARRATT","garratt","default")
        garatt = .true.
        write(16,*) 'Garratt drag law used'
#ifdef CSWAN
         CASE("Swell","SWELL","swell")
             swell   = .true.
             write(16,*) 'Swell drag law used'
#endif
      CASE("Powell","POWELL","powell")
        if (NWS.eq.12.or.NWS.eq.29.or.NWS.eq.30) then
          powell  = .true.
          allocate(Pweight(3,np)); Pweight = 0d0
          allocate(foundSector(np)); foundSector = .false.
          write(16,*) 'Powell drag law used'
        else
          garatt = .true.
          write(16,*) 'Garratt drag law used because NWS ',&
                      'incompatible with Powell'
        endif
      CASE DEFAULT
        WRITE(16,*) 'ERROR: Wind drag law not recognized:'
        WRITE(16,'(A10)') DragLawString
        WRITE(16,*) 'Execution will now be terminated.'
#ifdef CMPI
            call msg_fini()
#endif
        CALL EXIT(1)

      END SELECT

    end subroutine initialize_wind_drag_module

!----------------------------------------------------------------------
!     Sets an integer for the drag law to be used later in checks
!     This is simpler than checking a string
!----------------------------------------------------------------------
    SUBROUTINE checkWindDragType()
!----------------------------------------------------------------------
      USE GLOBAL,ONLY: toLowercase
      IMPLICIT NONE

      CHARACTER(LEN(dragLawString)) :: dragLawLowercase

      draglawlowercase = toLowercase(dragLawString)
      SELECT CASE(TRIM(draglawlowercase))
      CASE("default","garratt")
        draglawtype = DRAGLAW_GARRATT
      CASE("powell")
        draglawtype = DRAGLAW_POWELL
      CASE("swell")
        draglawtype = DRAGLAW_SWELL
      CASE DEFAULT
        draglawtype = DRAGLAW_GARRATT
      END SELECT
!----------------------------------------------------------------------
    END SUBROUTINE checkWindDragType
!----------------------------------------------------------------------

!   ----------------------------------------------------------------
!>  @brief Maps the function pointer used in the calculation of
!>  wind drag. This is used so that the logical statements required
!>  to select the wind drag are only executed once at startup
!   ----------------------------------------------------------------
    SUBROUTINE mapWindDragFunctionPointer()
!   ----------------------------------------------------------------
      use global, only: allMessage, error, warning
#ifdef CMPI
      use messenger, only: msg_fini
#endif
      IMPLICIT NONE
#ifdef POWELL
          f_computeWindDrag => PowellWindDrag
#else
      SELECT CASE(DragLawType)
      CASE(DRAGLAW_GARRATT)
        f_computeWindDrag => CappedGarrattWindDrag
      CASE(DRAGLAW_POWELL)
        f_computeWindDrag => PowellWindDrag
      CASE(DRAGLAW_SWELL)
#ifndef CSWAN
        call allMessage(WARNING,"SWAN Model not enabled."//&
                       "Drag formulation will be Garratt.")
        f_computeWindDrag => CappedGarrattWindDrag
#else
              f_computeWindDrag => SwanWindDrag
#endif
      CASE DEFAULT
        call allMessage(ERROR,"Wind drag law not recognized"//TRIM(DragLawString))
#ifdef CMPI
              call msg_fini()
#endif
        CALL EXIT(1)
      END SELECT
#endif
!     ----------------------------------------------------------------
    END SUBROUTINE mapWindDragFunctionPointer
!   ----------------------------------------------------------------

!   ---------------------------------------------------------------
!>  @brief Computes wind drag using a function with the signiture:
!>  result = function(windSpeed,nodeNumber)
!>  @param[in] WindSpeed wind magnitude
!>  @param[in] NodeNumber adcirc mesh node index
!>  @result wind drag coefficient
!   ----------------------------------------------------------------
    REAL(8) FUNCTION WindDrag(WindSpeed, NodeNumber)
      USE GLOBAL,ONLY: windDragOut
      IMPLICIT NONE
      REAL(8) :: WindSpeed
      INTEGER  :: NodeNumber
      WindDrag = f_computeWindDrag(windSpeed,NodeNumber)
      windDragOut(nodeNumber) = windDrag
      RETURN
!   ----------------------------------------------------------------
    END FUNCTION WindDrag
!   ----------------------------------------------------------------

!   ----------------------------------------------------------------
!>  @brief Function applying the Garratt wind drag formulation.
!>
!>  Wind drag is calculated using Garratt, 1977 by:
!>  \f[
!>     d = \frac{\left( 0.75 + 0.067 \left| \vec{w} \right| \right)}{1000}
!>  \f]
!>
!>  @param[in] WindSpeed wind magnitude
!>  @param[in] NodeNumber adcirc mesh node index
!>  @result wind drag coefficient
!   ----------------------------------------------------------------
    FUNCTION GarrattWindDrag(WindSpeed,NodeNumber) RESULT(Drag)
!   ----------------------------------------------------------------
      IMPLICIT NONE
      REAL(8),INTENT(IN)  :: WindSpeed
      INTEGER,INTENT(IN)  :: NodeNumber
      REAL(8)             :: Drag
      Drag = 0.001d0*(0.75d0+0.067d0*WindSpeed)
    END FUNCTION GarrattWindDrag

!   ----------------------------------------------------------------
!>  @brief Applies the Garratt Wind Drag capped at WindDragLimit.
!>  By default the cap applied is 0.0035, however this can be adjusted
!>  using the metControl namelist parameter WindDragLimit
!>  @param[in] WindSpeed wind magnitude
!>  @param[in] NodeNumber adcirc mesh node index
!>  @result wind drag coefficient
!   ----------------------------------------------------------------
    FUNCTION CappedGarrattWindDrag(WindSpeed,NodeNumber) RESULT(WindDrag)
!   ----------------------------------------------------------------
      IMPLICIT NONE
      REAL(8),INTENT(IN)  :: WindSpeed
      INTEGER,INTENT(IN)  :: NodeNumber
      REAL(8)             :: WindDrag
      WindDrag = MIN(GarrattWindDrag(WindSpeed,NodeNumber),WindDragLimit)
!   ----------------------------------------------------------------
    END FUNCTION CappedGarrattWindDrag
!   ----------------------------------------------------------------

!   ----------------------------------------------------------------
!>  @brief Applies Mark Powell's sector-based wind drag. This
!>  applies different wind drag values based upon the sector of the
!>  storm a node lies within. The is done to account for the direction
!>  of wave travel and the coincident direction of the wind vector
!>  @param[in] WindSpeed wind magnitude
!>  @param[in] NodeNumber adcirc mesh node index
!>  @result wind drag coefficient
!   ----------------------------------------------------------------
    FUNCTION PowellWindDrag(WindSpeed,NodeNumber) RESULT(PowellDrag)
!   ----------------------------------------------------------------
      USE GLOBAL, ONLY: DEG2RAD
      USE MESH, ONLY : SFEA, SLAM
      IMPLICIT NONE
      INTRINSIC           :: ATAN2
      INTRINSIC           :: MOD
      REAL(8),INTENT(IN)  :: WindSpeed
      INTEGER,INTENT(IN)  :: NodeNumber
      INTEGER             :: IS
      INTEGER             :: IW1,IW2,IW3
      LOGICAL             :: FoundSector
      REAL(8)             :: Dir1
      REAL(8)             :: Dir2
      REAL(8)             :: Drag(3)
      REAL(8)             :: NodeDirection
      REAL(8)             :: NodeLat
      REAL(8)             :: NodeLon
      REAL(8)             :: StormDirection
      REAL(8)             :: Weight(3)
      REAL(8)             :: PowellDrag
      REAL(8)             :: rawGarrattWindDrag

      !.. Check whether we have previously found the eye.
      IF((EyeLat(1)==0.D0).OR.(EyeLon(1)==0.D0).OR.&
         (EyeLat(2)==0.D0).OR.(EyeLon(2)==0.D0).OR.&
         (EyeLat(3)==0.D0).OR.(EyeLon(3)==0.D0))THEN
        FoundEye = .FALSE.
      ENDIF
      !.. If the eye has not been found then apply Garratt.
      IF(.NOT.FoundEye)THEN
        PowellDrag = CappedGarrattWindDrag(WindSpeed,NodeNumber)
        FoundSector = .FALSE.
        Weight = 0.D0
      ELSE
        !.. We have found the eye!
        !.. Compute the storm's direction.  We use the convention of northward
        !.. being zero degrees, and then cycling clockwise.
        Dir1 = ATAN2(EyeLat(2)-EyeLat(1),EyeLon(2)-EyeLon(1))
        Dir2 = ATAN2(EyeLat(3)-EyeLat(2),EyeLon(3)-EyeLon(2))
        StormDirection = ATAN2(SIN(Dir1)+SIN(Dir2),COS(Dir1)+COS(Dir2))
        StormDirection = StormDirection/DEG2RAD
        StormDirection = 90.D0 - StormDirection
        IF(StormDirection<0.D0)THEN
          StormDirection = StormDirection + 360.D0
        ENDIF
        !.. Compute the direction of node relative to the eye.
        NodeLon = SLAM(NodeNumber)/DEG2RAD
        NodeLat = SFEA(NodeNumber)/DEG2RAD
        NodeDirection = ATAN2(NodeLat-EyeLat(3),NodeLon-EyeLon(3))
        NodeDirection = NodeDirection/DEG2RAD
        NodeDirection = 90.D0 - NodeDirection
        IF(NodeDirection<0.D0)THEN
          NodeDirection = NodeDirection + 360.D0
        ENDIF
        !.. Compute the weights for the sector in which the node is located.
        !..     0- 40 : Transition from left (3) to right (1).
        !..    40-130 : Right (1).
        !..   130-170 : Transition from right (1) to rear (2).
        !..   170-220 : Rear (2).
        !..   220-260 : Transition from rear (2) to left (3).
        !..   260-  0 : Left (3).
        FoundSector = .FALSE.
        DO IS=1,6
          IF(IS==1)THEN
            Dir1 = MOD(StormDirection+  0.D0,360.D0)
            Dir2 = MOD(StormDirection+ 40.D0,360.D0)
          ELSEIF(IS==2)THEN
            Dir1 = MOD(StormDirection+ 40.D0,360.D0)
            Dir2 = MOD(StormDirection+130.D0,360.D0)
          ELSEIF(IS==3)THEN
            Dir1 = MOD(StormDirection+130.D0,360.D0)
            Dir2 = MOD(StormDirection+170.D0,360.D0)
          ELSEIF(IS==4)THEN
            Dir1 = MOD(StormDirection+170.D0,360.D0)
            Dir2 = MOD(StormDirection+220.D0,360.D0)
          ELSEIF(IS==5)THEN
            Dir1 = MOD(StormDirection+220.D0,360.D0)
            Dir2 = MOD(StormDirection+260.D0,360.D0)
          ELSEIF(IS==6)THEN
            Dir1 = MOD(StormDirection+260.D0,360.D0)
            Dir2 = MOD(StormDirection+  0.D0,360.D0)
          ENDIF
          IF(Dir1>Dir2)THEN
            IF((Dir1<=NodeDirection).AND.(NodeDirection<=360.D0))THEN
              Dir2 = Dir2 + 360.D0
            ELSEIF((0.D0<=NodeDirection).AND.(NodeDirection<=Dir2))THEN
              Dir1 = Dir1 - 360.D0
            ENDIF
          ENDIF
          IF((Dir1<=NodeDirection).AND.(NodeDirection<=Dir2))THEN
            IF(MOD(IS,2)==0)THEN
              FoundSector = .TRUE.
              Weight = 0.D0
              Weight(IS/2) = 1.D0
            ELSE
              FoundSector = .TRUE.
              Weight = 0.D0
              IF(IS==1)THEN
                IW1 = 3
                IW2 = 1
              ELSEIF(IS==3)THEN
                IW1 = 1
                IW2 = 2
              ELSEIF(IS==5)THEN
                IW1 = 2
                IW2 = 3
              ENDIF
              Weight(IW1) = 1.D0 + (0.D0-1.D0)/(Dir2-Dir1)*(NodeDirection-Dir1)
              Weight(IW2) = 0.D0 + (1.D0-0.D0)/(Dir2-Dir1)*(NodeDirection-Dir1)
            ENDIF
          ENDIF
        ENDDO


        !.. Apply garratt wind drag in case of emergency
        IF(.NOT.FoundSector)THEN
          PowellDrag = CappedGarrattWindDrag(WindSpeed,NodeNumber)
          Weight = 0.D0
          RETURN
        ELSE
          rawGarrattWindDrag = GarrattWindDrag(WindSpeed,NodeNumber)
          !.. Determine the wind drag for sector 1 (right).
          Drag(1) = rawGarrattWindDrag
          IF(WindSpeed<=35.0D0)THEN
            IF(Drag(1)>0.0020D0)THEN
              Drag(1) = 0.0020D0
            ENDIF
          ELSEIF(WindSpeed<=45.0D0)THEN
            Drag(1) = 0.0020D0 + (0.0030D0-0.0020D0)/(45.0D0-35.0D0)*(WindSpeed-35.0D0)
          ELSE
            Drag(1) = 0.0030D0
          ENDIF
          !.. Determine the wind drag for sector 2 (rear).
          Drag(2) = rawGarrattWindDrag
          IF(WindSpeed<=35.0D0)THEN
            IF(Drag(2)>0.0020D0)THEN
              Drag(2) = 0.0020D0
            ENDIF
          ELSEIF(WindSpeed<=45.0D0)THEN
            Drag(2) = 0.0020D0 + (0.0010D0-0.0020D0)/(45.0D0-35.0D0)*(WindSpeed-35.0D0)
          ELSE
            Drag(2) = 0.0010D0
          ENDIF
          !.. Determine the wind drag for sector 3 (left).
          Drag(3) = rawGarrattWindDrag
          IF(Drag(3)>0.0018D0)THEN
            IF(WindSpeed<=25.0D0)THEN
              Drag(3) = 0.0018D0
            ELSEIF(WindSpeed<=30.0D0)THEN
              Drag(3) = 0.0018D0 + (0.0045D0-0.0018D0)/(30.0D0-25.0D0)*(WindSpeed-25.0D0)
            ELSEIF(WindSpeed<=45.0D0)THEN
              Drag(3) = 0.0045D0 + (0.0010D0-0.0045D0)/(45.0D0-30.0D0)*(WindSpeed-30.0D0)
            ELSE
              Drag(3) = 0.0010D0
            ENDIF
          ENDIF
          !.. Apply a weighted average.
          PowellDrag = Weight(1)*Drag(1) + Weight(2)*Drag(2) + Weight(3)*Drag(3)
        ENDIF
      ENDIF
!   ----------------------------------------------------------------
    END FUNCTION PowellWindDrag
!   ----------------------------------------------------------------

!   ----------------------------------------------------------------
!    F U N C T I O N   W I N D  I C E  D R A G
!   ----------------------------------------------------------------
!
!   tcm v49.64.01 -- Added a function to modify wind drag coefficents
!        based on ice effects from a percentage of ice coverage
!        value (ice concentration limits 0.0 to 100.0)
!        If the percentage of ice is less than 1.0% then no
!        modifications are done to the incoming wind drag
!        coefficient values.
!
!   ----------------------------------------------------------------
    REAL(8) FUNCTION WindIceDrag(WindDrag,PercentIce)
#ifdef CMPI
      use messenger, only: msg_fini
#endif
      implicit none
      REAL(8) :: PercentIce
      REAL(8) :: pic, WindDrag, IceDrag, Cform
      REAL(8) :: Cskin = 1.5d-3, Cform_max = 3.67d-3, beta = 0.6d0

      SELECT CASE(TRIM(DragLawString))

      CASE("RaysIce")
!........If there is virtually no ice, then
!........just use the non-ice drag values
        if (PercentIce .lt. 1.d0) then
          WindIceDrag = WindDrag
        else
          pic = PercentIce*0.01d0
          IceDrag = (0.125d0 + 0.5d0*pic*(1.0d0-pic))*0.01d0
          WindIceDrag = max(IceDrag,WindDrag)
        endif

      CASE("IceCube")
!........This formula is a cubic function such that
!....... When pic = 0 the drag coefficient is equal to
!....... the minimum value for Garratt (.000075), when pic = 50
!....... it has a value of 0.0025 and a zero gradient, and
!....... when pic = 100 it has a value of 0.00125.
!....... This represents a smoother transition between Garratt
!....... and Ray's Ice for low ice values and low winds
        IF (PercentIce<0.d0) then
          WindIceDrag = WindDrag
        ELSE
          pic = PercentIce*0.01d0
          IceDrag = (0.075d0 + 0.75d0*pic - 0.9d0*pic*pic+ 0.2d0*pic*pic*pic)*0.01d0
          WindIceDrag = max(IceDrag,WindDrag)
        ENDIF

      CASE("Lupkes","default","Powell")
!........This formula is a function that explicitly
!....... combines air-sea drag (Garratt) with sea-ice skin drag and
!....... sea-ice form drag components with area fraction weights
!....... see: Lüpkes, et al. (2012), JGR: doi:10.1029/2012JD017630. &
!....... Joyce et al. (2019), OM: doi:10.1016/j.ocemod.2019.101421
        IF (PercentIce.lt.0.d0) then
          WindIceDrag = WindDrag
        ELSE
          pic = min(PercentIce*0.01d0,1d0)
          Cform = Cform_max*(1d0-pic)**beta
          WindIceDrag = WindDrag*(1d0 - pic) + (Cskin + Cform)*pic
        ENDIF

      CASE DEFAULT
        WRITE(16,*) 'ERROR: Ice drag law not recognized:'
        WRITE(16,'(A10)') DragLawString
        WRITE(16,*) 'Execution will now be terminated.'
#ifdef CMPI
         call msg_fini()
#endif
        CALL EXIT(1)

      END SELECT
      IF(WindIceDrag>WindDragLimit) WindIceDrag=WindDragLimit

      RETURN
!   ----------------------------------------------------------------
    END FUNCTION WindIceDrag
!   ----------------------------------------------------------------

#ifdef CSWAN
!     ----------------------------------------------------------------
!>    @brief Use the wind drag coefficients calculated inside SWAN
!>    for wind drag inside of adcirc
!>    @param[in] WindSpeed wind magnitude
!>    @param[in] NodeNumber adcirc mesh node index
!>    @result wind drag coefficient
!    ----------------------------------------------------------------
     FUNCTION SwanWindDrag(WindSpeed,NodeNumber) RESULT(Drag)
!    ----------------------------------------------------------------
           USE GLOBAL,ONLY: Swan_WDragCo
           IMPLICIT NONE
           REAL(8),INTENT(IN)  :: WindSpeed
           INTEGER,INTENT(IN)  :: NodeNumber
           REAL(8)             :: Drag
           Drag = SWAN_WDragCo(NodeNumber)
!    ----------------------------------------------------------------
       END FUNCTION SwanWindDrag
!    ----------------------------------------------------------------
#endif

end module wind_drag_module