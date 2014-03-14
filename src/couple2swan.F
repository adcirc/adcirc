C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      MODULE Couple2Swan
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      USE SIZES,  ONLY: SZ
      USE GLOBAL, ONLY : DEBUG, ECHO, INFO, WARNING, ERROR,
     &          setMessageSource, unsetMessageSource, allMessage
#ifdef CSWAN

      USE SWCOMM1, ONLY: NMOVAR

      IMPLICIT NONE

Casey 090302: The coupling interval controls how ADCIRC and SWAN
C             take turns during the simulation.  It is set
C             at the end of PADCSWAN_INIT.
      INTEGER :: CouplingInterval

Casey 090302: These logical variables will be reset to .TRUE.
C             if the water levels in SWAN are coupled to ADCIRC.
      LOGICAL :: COUPCUR
      LOGICAL :: COUPWIND
      LOGICAL :: COUPWLV
Casey 091216: Added another coupling variable.
      LOGICAL :: COUPFRIC

Casey 090302: These arrays contain the radiation stresses.
Casey 090820: Be explicit about the size of these REAL variables.
      REAL(4) ,ALLOCATABLE :: ADCIRC_SXX(:,:)
      REAL(4) ,ALLOCATABLE :: ADCIRC_SXY(:,:)
      REAL(4) ,ALLOCATABLE :: ADCIRC_SYY(:,:)

Casey 090302: The interpolation weight controls which value
C             is taken when information is passed between ADCIRC
C             and SWAN.  If InterpoWeight = 0, then the value
C             is taken from the beginning of the coupling interval.
C             If InterpoWeight = 1, then the value is taken from
C             the end of the coupling interval.
      REAL(SZ) :: InterpoWeight

Casey 090302: The following arrays will contain the ADCIRC values
C             that may be passed to SWAN.
      REAL(SZ),ALLOCATABLE :: SWAN_ETA2(:,:)
      REAL(SZ),ALLOCATABLE :: SWAN_UU2(:,:)
      REAL(SZ),ALLOCATABLE :: SWAN_VV2(:,:)
      REAL(SZ),ALLOCATABLE :: SWAN_WX2(:,:)
      REAL(SZ),ALLOCATABLE :: SWAN_WY2(:,:)
Casey 091216: Added another array.
      REAL(SZ),ALLOCATABLE :: SWAN_Z0(:,:)

Casey 090302: These variables control SWAN output.
      INTEGER, ALLOCATABLE :: IGSW(:)
      INTEGER              :: NumSwanOutput
      LOGICAL              :: Processed26 = .FALSE.
      LOGICAL              :: SWAN_OQPROC(NMOVAR)
      INTEGER              :: SWAN_VOQR(NMOVAR)

Casey 090303: Add variable for maximum SWAN output.
      REAL(SZ),ALLOCATABLE        :: SWAN_MAX(:,:)
      REAL(SZ),ALLOCATABLE,TARGET :: SWAN_MAXTEMP(:)
      REAL(SZ),ALLOCATABLE,TARGET :: SWAN_MAXTEMP_G(:) !...Cobell added

#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
      REAL(SZ),ALLOCATABLE :: TKXX(:)
      REAL(SZ),ALLOCATABLE :: TKXY(:)
      REAL(SZ),ALLOCATABLE :: TKYY(:)
      REAL(SZ),ALLOCATABLE :: WAVE_A(:)
      REAL(SZ),ALLOCATABLE :: WAVE_A1(:)
      REAL(SZ),ALLOCATABLE :: WAVE_A2(:)
      REAL(SZ),ALLOCATABLE :: WAVE_H(:)
      REAL(SZ),ALLOCATABLE :: WAVE_H1(:)
      REAL(SZ),ALLOCATABLE :: WAVE_H2(:)
      REAL(SZ),ALLOCATABLE :: WAVE_T(:)
      REAL(SZ),ALLOCATABLE :: WAVE_T1(:)
      REAL(SZ),ALLOCATABLE :: WAVE_T2(:)
#endif

Casey 100205: Add variables for writing of SWAN hot-start files.
      INTEGER              :: SwanHotStartUnit
      LOGICAL              :: WriteSwanHotStart = .FALSE.

Casey 120302: Merging ice changes from Taylor Asher.
CTGA 110912:  Implementing SWAN ice handling.  This variable is set to
C             0 wherever ice reaches/exceeds a threshold value
      INTEGER,ALLOCATABLE :: ICY(:)
      REAL(SZ),ALLOCATABLE :: SWAN_CICE2(:,:)

Cobell 20120510:  Moved variables for SWAN Output here
C          Added the rest of the local SWAN arrays
      REAL(SZ),ALLOCATABLE, TARGET :: SW_HSOut(:), SW_DIROut(:)
      REAL(SZ),ALLOCATABLE, TARGET :: SW_TM01Out(:), SW_TPSOut(:)
      REAL(SZ),ALLOCATABLE, TARGET :: SW_WindXOut(:), SW_WindYOut(:)
      REAL(SZ),ALLOCATABLE, TARGET :: SW_TM02Out(:), SW_TMM10Out(:)
      REAL(SZ),ALLOCATABLE, TARGET :: SW_HSMaxOut(:), SW_DIRMaxOut(:)
      REAL(SZ),ALLOCATABLE, TARGET :: SW_TM01MaxOut(:), SW_TPSMaxOut(:)
      REAL(SZ),ALLOCATABLE, TARGET :: SW_WindMaxOut(:)
      REAL(SZ),ALLOCATABLE, TARGET :: SW_TM02MaxOut(:),SW_TMM10MaxOut(:)


C-----------------------------------------------------------------------
      CONTAINS
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E
C       C O M P U T E  M O D I F I E D  W A V E  F R I C T I O N
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
#ifdef CSWANFRIC
!Casey 091020: Adopt Ethan's/Joannes's modified friction.
      SUBROUTINE ComputeModifiedWaveFriction(TK)

      USE GLOBAL, ONLY         : DP,
     &                           ETA2,
     &                           G,
     &                           H0,
     &                           IFNLFA,
     &                           NP,
     &                           PI,
     &                           UU1,
     &                           VV1
      USE NodalAttributes, ONLY: FGAMMA,
     &                           FRIC,
     &                           FTHETA,
     &                           HBREAK,
     &                           IFHYBF,
     &                           IFLINBF,
     &                           IFNLBF
      USE SIZES, ONLY          : SZ

      IMPLICIT NONE

      INTRINSIC              :: ABS
      INTRINSIC              :: COS
      INTRINSIC              :: SIN
      INTRINSIC              :: SINH
      INTRINSIC              :: SQRT

      REAL(SZ),INTENT(IN)    :: TK(:)

      INTEGER                :: IP

      REAL(SZ)               :: COSA1
      REAL(SZ)               :: CHYBR
      REAL(SZ)               :: F_WAVE
      REAL(SZ)               :: H1
      REAL(SZ)               :: L_WAVE
      REAL(SZ)               :: SINA1
      REAL(SZ)               :: UVW1
      REAL(SZ)               :: UVW2
      REAL(SZ)               :: WB
      REAL(SZ)               :: WD
      REAL(SZ)               :: WD1
      REAL(SZ)               :: WDXX
      REAL(SZ)               :: WDXY
      REAL(SZ)               :: WDYY
      REAL(SZ)               :: Y_WAVE

      call setMessageSource("ComputeModifiedWaveFriction")
#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif


      DO IP=1,NP
         H1 = DP(IP) + IFNLFA*ETA2(IP)
         TKXX(IP) = TK(IP)
         TKYY(IP) = TK(IP)
         TKXY(IP) = 0.D0
C... Added friction due to waves.
         IF((WAVE_H(IP).GT.0.D0).AND.(WAVE_T(IP).GT.0.D0).AND.(H1.GT.H0))THEN
C... Compute the wave length using a Pade approximation.
            Y_WAVE = (2.D0*PI/WAVE_T(IP))**2.D0*H1/G
            F_WAVE = Y_WAVE + 1.D0/(1.D0 + Y_WAVE*(0.66667D0 +
     &               Y_WAVE*(0.35550D0   + Y_WAVE*(0.16084D0 +
     &               Y_WAVE*(0.06320D0   + Y_WAVE*(0.02174D0 +
     &               Y_WAVE*(0.00654D0   + Y_WAVE*(0.00171D0 +
     &               Y_WAVE*(0.00039D0   + Y_WAVE*0.000111D0)))))))))
            L_WAVE = 2.D0*PI/SQRT(ABS(Y_WAVE*F_WAVE/(H1**2.D0)))
C... Compute the orbital wave velocity.
            WB  = 2.D0*PI/WAVE_T(IP)*WAVE_H(IP)/(PI*SINH(2.D0*PI/
     &            L_WAVE*H1))
C... Compute the modified friction factors due to waves.
            COSA1 = COS(WAVE_A(IP))
            SINA1 = SIN(WAVE_A(IP))
            UVW1  = UU1(IP)*UU1(IP) + VV1(IP)*VV1(IP) + WB*WB
            UVW2  = 2.D0*(UU1(IP)*COSA1 + VV1(IP)*SINA1)*WB
            WD1   = 0.5D0*((ABS(UVW1 + UVW2))**0.5D0
     &            + (ABS(UVW1 - UVW2))**0.5D0)
            WD    = 0.D0
            IF(ABS(WD1).GE.0.001D0)THEN
               WD = WB*WB/WD1
            ENDIF
            WDXX = WD*COSA1*COSA1
            WDYY = WD*SINA1*SINA1
            WDXY = WD*SINA1*COSA1
            CHYBR = IFNLBF + IFHYBF*
     &              (1+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)
            TKXX(IP) = FRIC(IP)*(IFLINBF + ((WD1 + WDXX)/H1)*CHYBR)
            TKYY(IP) = FRIC(IP)*(IFLINBF + ((WD1 + WDYY)/H1)*CHYBR)
            TKXY(IP) = FRIC(IP)*(IFLINBF + WDXY/H1*CHYBR)
         ENDIF
      ENDDO
#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE ComputeModifiedWaveFriction
C-----------------------------------------------------------------------
#endif


C-----------------------------------------------------------------------
C     S U B R O U T I N E
C       C O M P U T E  R A D I A T I O N  S T R E S S E S
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE ComputeRadiationStresses(AC2,SPCDIR,SPCSIG)

!Casey 090616: On Jade, there was a problem with getting water depths
!              from the COMPDA array.  Use ETA2 and DEPTH instead.
      USE GLOBAL, ONLY        : ETA2
      USE M_GENARR, ONLY      : DEPTH
      USE SwanGridData, ONLY  : nverts
      USE SWCOMM3, ONLY       : DDIR, DEPMIN, FRINTF, GRAV, MDC, MSC, RHO

      IMPLICIT NONE

!Casey 090616: These REAL variables must be size(4) to match SWAN.
      REAL(4) :: AC2(MDC,MSC,nverts)
      INTEGER :: IT
      REAL    :: JunkR
      REAL(4) :: SPCDIR(MDC,6)
      REAL(4) :: SPCSIG(MSC)

      INTEGER :: I, ID, IERR, IS
      REAL(4) :: DEPLOC
      REAL(4) :: WK(MSC), CG(MSC), NE(MSC), NED(MSC)

      call setMessageSource("ComputeRadiationStresses")
#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

!... (Re)allocate the arrays.
      IF(.NOT.ALLOCATED(ADCIRC_SXX))THEN
         ALLOCATE(ADCIRC_SXX(1:nverts,1:2),STAT=IERR)
         DO I=1,nverts
            ADCIRC_SXX(I,1) = 0.0
            ADCIRC_SXX(I,2) = 0.0
         ENDDO
      ENDIF
      IF(.NOT.ALLOCATED(ADCIRC_SXY))THEN
         ALLOCATE(ADCIRC_SXY(1:nverts,1:2))
         DO I=1,nverts
            ADCIRC_SXY(I,1) = 0.0
            ADCIRC_SXY(I,2) = 0.0
         ENDDO
      ENDIF
      IF(.NOT.ALLOCATED(ADCIRC_SYY))THEN
         ALLOCATE(ADCIRC_SYY(1:nverts,1:2))
         DO I=1,nverts
            ADCIRC_SYY(I,1) = 0.0
            ADCIRC_SYY(I,2) = 0.0
         ENDDO
      ENDIF
!... Loop over the nodes and compute radiation stresses.
      DO I=1,nverts
!... Transfer the stresses from the new time level to the old time level.
         ADCIRC_SXX(I,1) = ADCIRC_SXX(I,2)
         ADCIRC_SXY(I,1) = ADCIRC_SXY(I,2)
         ADCIRC_SYY(I,1) = ADCIRC_SYY(I,2)
!... Initialize the radiation stresses as zero at all of the nodes.
         ADCIRC_SXX(I,2) = 0.0
         ADCIRC_SXY(I,2) = 0.0
         ADCIRC_SYY(I,2) = 0.0
!... Compute ratio of group and phase velocity.
         DEPLOC = DEPTH(I) + ETA2(I)
         IF(DEPLOC.LE.DEPMIN)THEN
            CYCLE
         ENDIF
         CALL KSCIP1(MSC, SPCSIG, DEPLOC, WK, CG, NE, NED)
!... Loop over all sigma and theta.
         DO ID=1,MDC
            DO IS=1,MSC
!... Sum contributions to radiation stresses.
               ADCIRC_SXX(I,2) = ADCIRC_SXX(I,2)
     &                   + (NE(IS) * SPCDIR(ID,4) + NE(IS) - 0.5)
     &                   * SPCSIG(IS) * SPCSIG(IS) * AC2(ID,IS,I)
               ADCIRC_SXY(I,2) = ADCIRC_SXY(I,2)
     &                   + NE(IS) * SPCDIR(ID,5)
     &                   * SPCSIG(IS) * SPCSIG(IS) * AC2(ID,IS,I)
               ADCIRC_SYY(I,2) = ADCIRC_SYY(I,2)
     &                   + (NE(IS) * SPCDIR(ID,6) + NE(IS) - 0.5)
     &                   * SPCSIG(IS) * SPCSIG(IS) * AC2(ID,IS,I)
            ENDDO
         ENDDO
!... Multiply summed radiation stresses by the stuff outside of the sums.
         ADCIRC_SXX(I,2) = RHO * GRAV * ADCIRC_SXX(I,2) * DDIR * FRINTF
         ADCIRC_SXY(I,2) = RHO * GRAV * ADCIRC_SXY(I,2) * DDIR * FRINTF
         ADCIRC_SYY(I,2) = RHO * GRAV * ADCIRC_SYY(I,2) * DDIR * FRINTF
!Casey 080602: ADCIRC accepts wave-driven stresses "in units of velocity squared
!    (consistent with the units of gravity).  Stress in these units is obtained
!    by dividing stress in units of force/area by the reference density of water."
!    SO WE MUST DIVIDE BY RHO!
         ADCIRC_SXX(I,2) = ADCIRC_SXX(I,2) / RHO
         ADCIRC_SXY(I,2) = ADCIRC_SXY(I,2) / RHO
         ADCIRC_SYY(I,2) = ADCIRC_SYY(I,2) / RHO
!... End loop over wet nodes.  The radiation stresses are ready for ADCIRC.
      ENDDO

#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE ComputeRadiationStresses
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E
C       C O M P U T E  W A V E  D R I V E N  F O R C E S
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE ComputeWaveDrivenForces

      USE GLOBAL, ONLY: AREAS, NBDV, NBOU, NBVV, NE, NEITABELE,
     &                  NM, NODECODE, NODELE, NOFF, NOPE,
     &                  NP, NVDLL, NVELL, RSNX2, RSNY2, X, Y
      USE SIZES, ONLY: SZ

      IMPLICIT NONE

      INTEGER :: I
      INTEGER :: IE
      INTEGER :: IP
      INTEGER :: K
      INTEGER :: Node1
      INTEGER :: Node2
      INTEGER :: Node3
      INTEGER :: NUMFOUND

      LOGICAL :: Marcel = .FALSE.

      REAL(SZ),ALLOCATABLE :: DSXXDX(:)
      REAL(SZ),ALLOCATABLE :: DSXYDY(:)
      REAL(SZ),ALLOCATABLE :: DSXYDX(:)
      REAL(SZ),ALLOCATABLE :: DSYYDY(:)

      REAL(SZ)             :: NCELE

      REAL(SZ),ALLOCATABLE :: TEMP_SXX(:)
      REAL(SZ),ALLOCATABLE :: TEMP_SXY(:)
      REAL(SZ),ALLOCATABLE :: TEMP_SYY(:)

      REAL(SZ) :: TOTALAREA

      call setMessageSource("ComputeWaveDrivenForces")
#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

!... Check whether radiation stresses have already been computed.
!... If not, then apply forces of zero.
      IF( .FALSE. )THEN

         IF(.NOT.ALLOCATED(RSNX2))THEN
            ALLOCATE(RSNX2(1:NP))
            DO IP=1,NP
               RSNX2(IP) = 0.D0
            ENDDO
         ENDIF
         IF(.NOT.ALLOCATED(RSNY2))THEN
            ALLOCATE(RSNY2(1:NP))
            DO IP=1,NP
               RSNY2(IP) = 0.D0
            ENDDO
         ENDIF

!... If so, then continue to compute wave-driven forces.
      ELSE

!... Allocate arrays for radiation stresses.
         IF(.NOT.ALLOCATED(TEMP_SXX)) ALLOCATE(TEMP_SXX(1:NP))
         IF(.NOT.ALLOCATED(TEMP_SXY)) ALLOCATE(TEMP_SXY(1:NP))
         IF(.NOT.ALLOCATED(TEMP_SYY)) ALLOCATE(TEMP_SYY(1:NP))

!... Loop over all nodes and interpolate the radiation stress for this time step.
         DO IP=1,NP
           TEMP_SXX(IP) = (1.0 - InterpoWeight) * DBLE(ADCIRC_SXX(IP,1))
     &                 + InterpoWeight * DBLE(ADCIRC_SXX(IP,2))
           TEMP_SXY(IP) = (1.0 - InterpoWeight) * DBLE(ADCIRC_SXY(IP,1))
     &                 + InterpoWeight * DBLE(ADCIRC_SXY(IP,2))
           TEMP_SYY(IP) = (1.0 - InterpoWeight) * DBLE(ADCIRC_SYY(IP,1))
     &                 + InterpoWeight * DBLE(ADCIRC_SYY(IP,2))
         ENDDO

!... Allocate arrays for radiation stress gradients.
         IF(.NOT.ALLOCATED(DSXXDX)) ALLOCATE(DSXXDX(1:NE))
         IF(.NOT.ALLOCATED(DSXYDY)) ALLOCATE(DSXYDY(1:NE))
         IF(.NOT.ALLOCATED(DSXYDX)) ALLOCATE(DSXYDX(1:NE))
         IF(.NOT.ALLOCATED(DSYYDY)) ALLOCATE(DSYYDY(1:NE))

!... Loop over all elements and compute the derivatives of Sxx, Sxy and Syy.
!... These derivatives are constant on an element.  Note that the AREAS array
!... actually contains twice the area of each element.
         DO IE=1,NE

!Casey 090707: When using the serial adcswan on Zas, I received memory errors
!... when these calls were nested into the logic below.  Break them out and
!... use these variables to save on the number of calls to memory.
           Node1 = NM(IE,1)
           Node2 = NM(IE,2)
           Node3 = NM(IE,3)

           DSXXDX(IE) = (1.D0/AREAS(IE)) *
     &                ( TEMP_SXX(Node1) * (Y(Node2) - Y(Node3))
     &                + TEMP_SXX(Node2) * (Y(Node3) - Y(Node1))
     &                + TEMP_SXX(Node3) * (Y(Node1) - Y(Node2)) )

           DSXYDY(IE) = (1.D0/AREAS(IE)) *
     &                ( TEMP_SXY(Node1) * (X(Node3) - X(Node2))
     &                + TEMP_SXY(Node2) * (X(Node1) - X(Node3))
     &                + TEMP_SXY(Node3) * (X(Node2) - X(Node1)) )

           DSXYDX(IE) = (1.D0/AREAS(IE)) *
     &                ( TEMP_SXY(Node1) * (Y(Node2) - Y(Node3))
     &                + TEMP_SXY(Node2) * (Y(Node3) - Y(Node1))
     &                + TEMP_SXY(Node3) * (Y(Node1) - Y(Node2)) )

           DSYYDY(IE) = (1.D0/AREAS(IE)) *
     &                ( TEMP_SYY(Node1) * (X(Node3) - X(Node2))
     &                + TEMP_SYY(Node2) * (X(Node1) - X(Node3))
     &                + TEMP_SYY(Node3) * (X(Node2) - X(Node1)) )

         ENDDO

!... Allocate arrays for wave-driven forces.
         IF(.NOT.ALLOCATED(RSNX2)) ALLOCATE(RSNX2(1:NP))
         IF(.NOT.ALLOCATED(RSNY2)) ALLOCATE(RSNY2(1:NP))

!... Loop over all nodes and compute the wave-driven forces:
!...
!...       Fx = - DSxx/Dx - DSxy/Dy
!...
!...       Fy = - DSxy/Dx - DSyy/Dy
!...
!... We project the element-based radiation stress gradients onto the nodes
!... by taking a weighted average of the gradients in the elements connected
!... to a node.
         outer: DO IP=1,NP

           RSNX2(IP) = 0.D0
           RSNY2(IP) = 0.D0

           TOTALAREA = 0.D0

           IE = 0
           NUMFOUND = 0

           inner: DO

             IE = IE + 1

             IF(NEITABELE(IP,IE).EQ.0)THEN

               CONTINUE

             ELSE

!... Try Marcel's method of zero-ing out the forces at nodes connected to dry nodes/elements.

               NCELE = NODECODE(NM(NEITABELE(IP,IE),1))
     &               * NODECODE(NM(NEITABELE(IP,IE),2))
     &               * NODECODE(NM(NEITABELE(IP,IE),3))
     &               * NOFF(       NEITABELE(IP,IE)   )

               IF(Marcel.AND.(NCELE.EQ.0))THEN

                 RSNX2(IP) = 0.0
                 RSNY2(IP) = 0.0
                 CYCLE outer

               ELSE

                 NUMFOUND = NUMFOUND + 1

                 RSNX2(IP) = RSNX2(IP) + 0.5*AREAS(NEITABELE(IP,IE))
     &                     * ( - DSXXDX(NEITABELE(IP,IE))
     &                         - DSXYDY(NEITABELE(IP,IE)) )
                 RSNY2(IP) = RSNY2(IP) + 0.5*AREAS(NEITABELE(IP,IE))
     &                     * ( - DSXYDX(NEITABELE(IP,IE))
     &                         - DSYYDY(NEITABELE(IP,IE)) )

                 TOTALAREA = TOTALAREA + 0.5*AREAS(NEITABELE(IP,IE))

               ENDIF

             ENDIF

             IF(NUMFOUND.EQ.NODELE(IP))THEN

               EXIT inner

             ENDIF

           ENDDO inner

           RSNX2(IP) = RSNX2(IP) / TOTALAREA
           RSNY2(IP) = RSNY2(IP) / TOTALAREA

         ENDDO outer

!... Try Marcel's method of zero-ing the forces at the boundary nodes.
         IF(Marcel)THEN
           DO K=1,NOPE
             DO I=1,NVDLL(K)
               RSNX2(NBDV(K,I)) = 0.0
               RSNY2(NBDV(K,I)) = 0.0
             ENDDO
           ENDDO
           DO K=1,NBOU
             DO I=1,NVELL(K)
               RSNX2(NBVV(K,I)) = 0.0
               RSNY2(NBVV(K,I)) = 0.0
             ENDDO
           ENDDO
         ENDIF

!... Deallocate the radiation stress gradients.
         IF(ALLOCATED(DSXXDX)) DEALLOCATE(DSXXDX)
         IF(ALLOCATED(DSXYDY)) DEALLOCATE(DSXYDY)
         IF(ALLOCATED(DSXYDX)) DEALLOCATE(DSXYDX)
         IF(ALLOCATED(DSYYDY)) DEALLOCATE(DSYYDY)
         IF(ALLOCATED(TEMP_SXX)) DEALLOCATE(TEMP_SXX)
         IF(ALLOCATED(TEMP_SXY)) DEALLOCATE(TEMP_SXY)
         IF(ALLOCATED(TEMP_SYY)) DEALLOCATE(TEMP_SYY)

      ENDIF

#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE ComputeWaveDrivenForces
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E
C       C O M P U T E  W A V E  F R I C T I O N  P R O P E R T I E S
C-----------------------------------------------------------------------
Casey 091020: Adopt Ethan's/Joannes's modified friction.
Casey 091021: Copy the format of this routine from the SwanOutput
C             routine earlier in this file.  We want to interact
C             with the SWAN output routines so that they will
C             compute the wave heights, periods and directions.
C             we need for the modified friction.
C-----------------------------------------------------------------------
#ifdef CSWANFRIC
      SUBROUTINE ComputeWaveFrictionProperties

      USE Couple2Adcirc,ONLY: COMPDA
      USE GLOBAL,       ONLY: DP,
     &                        ETA2,
     &                        NP,
     &                        PI,
     &                        UU2,
     &                        VV2
      USE M_GENARR,     ONLY: AC2,
     &                        KGRPNT,
     &                        SPCDIR,
     &                        SPCSIG
      USE SwanGriddata, ONLY: xcugrd,
     &                        ycugrd
      USE SWCOMM1,      ONLY: COSCQ,
     &                        OVEXCV,
     &                        SINCQ
      USE SWCOMM2,      ONLY: XOFFS,
     &                        YOFFS
      USE SWCOMM3,      ONLY: MCGRD,
     &                        MDC,
     &                        MSC,
     &                        MTC,
     &                        MXC,
     &                        MYC

      IMPLICIT NONE

      INTRINSIC                     :: ALLOCATED
      INTRINSIC                     :: INDEX
      INTRINSIC                     :: TRIM

      CHARACTER(LEN=30)             :: TempC

      INTEGER                       :: I
      INTEGER                       :: IO
      INTEGER                       :: IONOD(NP)
      INTEGER                       :: IP
      INTEGER                       :: IS
      INTEGER                       :: IVTYPE
      INTEGER                       :: SWAN_BKC
      INTEGER,SAVE                  :: SWAN_MTC

      LOGICAL                       :: CROSS(4,NP)

      REAL(4)                       :: ACLOC(MDC,MSC)
      REAL(4)                       :: DEPXY(NP)
      REAL(4)                       :: FORCE(NP,2)
      REAL(4)                       :: SWAN_CG(MSC)
      REAL(4)                       :: SWAN_NE(MSC)
      REAL(4)                       :: SWAN_NED(MSC)
      REAL(4),ALLOCATABLE           :: SWAN_VOQ(:,:)
      REAL(4)                       :: SWAN_WK(MSC)
      REAL(4)                       :: XC(NP)
      REAL(4)                       :: YC(NP)

      call setMessageSource("ComputeWaveFrictionProperties")
#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

!... Initialize the Swan output variables.
      DO IO=1,NMOVAR
         SWAN_OQPROC(IO) = .FALSE.
      ENDDO
      NumSwanOutput = 0
!... For significant wave heights (HS).
      NumSwanOutput = NumSwanOutput + 1
      OVEXCV(10) = 0.
      SWAN_OQPROC(10) = .TRUE.
      SWAN_VOQR(10) = 7 + NumSwanOutput
!... For mean wave directions (DIR).
      NumSwanOutput = NumSwanOutput + 1
      OVEXCV(13) = 0.
      SWAN_OQPROC(13) = .TRUE.
      SWAN_VOQR(13) = 7 + NumSwanOutput
!... For mean wave periods (TM01).
      NumSwanOutput = NumSwanOutput + 1
      OVEXCV(11) = 0.
      SWAN_OQPROC(11) = .TRUE.
      SWAN_VOQR(11) = 7 + NumSwanOutput
      COSCQ = COS(0.)
      SINCQ = SIN(0.)
      IF(.NOT.ALLOCATED(IGSW))  ALLOCATE(IGSW(1:NumSwanOutput))

!... I don't know what BKC does.  Everything seems to work okay
!... when it is set to two, though.
      SWAN_BKC = 2
!... If Swan needs these values for wave numbers, etc.,
!... then maybe it solves for them inside SWOEXA.
      DO IS=1,MSC
         SWAN_CG(IS) = 0.
         SWAN_NE(IS) = 0.
         SWAN_NED(IS) = 0.
         SWAN_WK(IS) = 0.
      ENDDO
!... Allocate VOQ and set up the first seven entries inside it.
      IF(.NOT.ALLOCATED(SWAN_VOQ)) ALLOCATE(SWAN_VOQ(1:NP,1:(7+NumSwanOutput+1)))
      DO IP=1,NP
!... X,Y in problem coordinate system.
         SWAN_OQPROC(1) = .TRUE.
         SWAN_OQPROC(2) = .TRUE.
         SWAN_VOQR(1) = 1
         SWAN_VOQR(2) = 2
         SWAN_VOQ(IP,SWAN_VOQR(1)) = xcugrd(IP)
         SWAN_VOQ(IP,SWAN_VOQR(2)) = ycugrd(IP)
!... X,Y in adjusted coordinate system.
         SWAN_OQPROC(24) = .TRUE.
         SWAN_OQPROC(25) = .TRUE.
         SWAN_VOQR(24) = 3
         SWAN_VOQR(25) = 4
         SWAN_VOQ(IP,SWAN_VOQR(24)) = 0.
         SWAN_VOQ(IP,SWAN_VOQR(25)) = 0.
!... Depth of water.
         SWAN_OQPROC(4) = .TRUE.
         SWAN_VOQR(4) = 5
         SWAN_VOQ(IP,SWAN_VOQR(4)) = REAL(DP(IP)) + REAL(ETA2(IP))
!... Water current velocities.
         SWAN_OQPROC(5) = .TRUE.
         SWAN_VOQR(5) = 6
         SWAN_VOQ(IP,SWAN_VOQR(5))   = REAL(UU2(IP))
         SWAN_VOQ(IP,SWAN_VOQR(5)+1) = REAL(VV2(IP))
!... Assemble other arrays?
         DEPXY(IP) = REAL(DP(IP)) + REAL(ETA2(IP))
         XC(IP) = xcugrd(IP)
         YC(IP) = ycugrd(IP)
         FORCE(IP,1) = 0.
         FORCE(IP,2) = 0.
         IONOD(IP) = 0
      ENDDO

!... Call the Swan subroutine to interpolate internal Swan quantities.
      CALL SWOEXD(SWAN_OQPROC, NP, XC, YC, SWAN_VOQR,
     &            SWAN_VOQ, COMPDA, KGRPNT, FORCE, CROSS, IONOD, -999)

!... Call the Swan subroutine to compute the output quantities.
      CALL SWOEXA(SWAN_OQPROC, SWAN_BKC, NP, XC, YC,
     &            SWAN_VOQR, SWAN_VOQ, AC2, ACLOC, SPCSIG,
     &            SWAN_WK, SWAN_CG, SPCDIR, SWAN_NE, SWAN_NED,
     &            KGRPNT, DEPXY, CROSS)

!... Significant wave heights (HS).
      IVTYPE = 10
      DO IP=1,NP
         WAVE_H2(IP) = DBLE(SWAN_VOQ(IP,SWAN_VOQR(IVTYPE)))
         WRITE(UNIT=TempC,FMT=*) WAVE_H2(IP)
      ENDDO
!... Mean wave directions (DIR).
      IVTYPE = 13
      DO IP=1,NP
         WAVE_A2(IP) = DBLE(SWAN_VOQ(IP,SWAN_VOQR(IVTYPE)))
         WRITE(UNIT=TempC,FMT=*) WAVE_A2(IP)
         WAVE_A2(IP) = WAVE_A2(IP) * PI/180.D0
      ENDDO
!... Mean wave periods (TM01).
      IVTYPE = 11
      DO IP=1,NP
         WAVE_T2(IP) = DBLE(SWAN_VOQ(IP,SWAN_VOQR(IVTYPE)))
         WRITE(UNIT=TempC,FMT=*) WAVE_T2(IP)
      ENDDO

      IF(ALLOCATED(IGSW))     DEALLOCATE(IGSW)
      IF(ALLOCATED(SWAN_VOQ)) DEALLOCATE(SWAN_VOQ)

#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE ComputeWaveFrictionProperties
C-----------------------------------------------------------------------
#endif


C-----------------------------------------------------------------------
C     S U B R O U T I N E  M A N N I N G  2  M A D S E N
C-----------------------------------------------------------------------
Casey 091216: This routine will convert the ADCIRC Manning's n values
C             into roughness lengths that can be used with the Madsen
C             friction formulation inside SWAN.
C-----------------------------------------------------------------------
      SUBROUTINE Manning2Madsen

      USE GLOBAL,         ONLY: DP,
     &                          G,
     &                          NP
      USE NodalAttributes,ONLY: ManningsN

      IMPLICIT NONE

      INTEGER  :: IN

      REAL(SZ) :: H
      REAL(SZ) :: K = 0.4D0
      REAL(SZ) :: N
      REAL(SZ) :: Z0

      call setMessageSource("Manning2Madsen")
#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      DO IN=1,NP

         H  = SWAN_ETA2(IN,2) + DP(IN)
         N  = ManningsN(IN)

Casey 110518: Enforce a lower limit on the Manning's n seen by SWAN.
         IF(N.LT.0.02D0) N = 0.02D0

         Z0 = ( H ) * EXP( -1.D0 * ( 1.D0 + K * H**(1.D0/6.D0)
     &      / ( N * SQRT(G) ) ) )

Casey 091216: If we get a junk number, then use the default value.
         IF(Z0.LE.0.D0)THEN
            Z0 = 0.05D0
         ENDIF

         SWAN_Z0(IN,2) = Z0

      ENDDO

#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE Manning2Madsen
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E    S W A N  O U T P U T
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SwanOutput(ITIME, IT)

      USE Couple2Adcirc,ONLY: COMPDA
      USE GLOBAL,       ONLY: AID4,
     &                       DP,
     &                        DT,
     &                        ETA2,
     &                        NDSETSW,
     &                        NODECODE,
     &                        NODES_LG,
     &                        NOUTGW,
     &                        NP,
     &                        NP_G,
     &                        NSCOUGW,
     &                        NSPOOLGW,
     &                        NTCYSGW,
     &                        NTCYFGW,
     &                        NWS,
     &                        OutputDataDescript_t,
     &                        RDES4,
     &                        RID4,
     &                        UU2,
     &                        VV2
#ifdef CSWAN
     &                        ,SWAN_OutputHS,
     &                        SWAN_OutputDIR,
     &                        SWAN_OutputTM01,
     &                        SWAN_OutputTPS,
     &                        SWAN_OutputWIND,
     &                        SWAN_OutputTM02,
     &                        SWAN_OutputTMM10,
     &                        SWAN_OutputAgg
#endif

      USE GLOBAL_IO,    ONLY: Header73,
     &                        Header74,
     &                        HEADER_MAX,
     &                        OPEN_GBL_FILE,
     &                        OPEN_MINMAX_FILE,
     &                        PackOne,
     &                        PackTwo,
     &                        StoreOne,
     &                        StoreTwo,
     &                        UnPackOne,
     &                        UnPackTwo,
     &                        WRITE_GBL_FILE,
     &                        WRITE_GBL_FILE_SKIP_DEFAULT
      USE M_GENARR,     ONLY: AC2,
     &                        KGRPNT,
     &                        SPCDIR,
     &                        SPCSIG
      USE SIZES,        ONLY: GLOBALDIR,
     &                        LOCALDIR,
     &                        MNPROC,
     &                        MNWPROC,
     &                        NBYTE,
     &                        SZ,
     &                        MYPROC
      USE SwanGriddata, ONLY: nverts,
     &                        xcugrd,
     &                        ycugrd
      USE SWCOMM1,      ONLY: COSCQ,
     &                        OVEXCV,
     &                        SINCQ
      USE SWCOMM2,      ONLY: XOFFS,
     &                        YOFFS
      USE SWCOMM3,      ONLY: MCGRD,
     &                        MDC,
     &                        MSC,
     &                        MTC,
     &                        MXC,
     &                        MYC
#ifdef CMPI
      USE WRITER,       ONLY: FLUSH_WRITERS,
     &                        NUM_BUF_MAX,  !st3 100708: check writer buffer
     &                        WRITE_GBL_FILE_THROUGH_WRITER,
     &                        Writer_StoreOne,
     &                        Writer_StoreTwo
      USE MESSENGER,   ONLY : MSG_FINI,       !st3 100708: check writer buffer
     &                        subDomainFatalError
#endif

#ifdef ADCNETCDF
      USE NETCDFIO,    ONLY : initNetCDFOutputFile
#endif
      IMPLICIT NONE

      INTRINSIC                     :: ALLOCATED
      INTRINSIC                     :: INDEX
      INTRINSIC                     :: TRIM

      CHARACTER(LEN=15),ALLOCATABLE :: FileName(:)
      CHARACTER(LEN=20)             :: FileNameMax
      CHARACTER(LEN=10),ALLOCATABLE :: Names(:)
      CHARACTER(LEN=30)             :: TempC

      INTEGER                       :: I
      INTEGER                       :: IO
      INTEGER                       :: IONOD(NP)
      INTEGER                       :: IP
      INTEGER                       :: IS
      INTEGER                       :: IT
      INTEGER                       :: ITIME
      INTEGER                       :: IVTYPE
      INTEGER                       :: IW
      INTEGER                       :: SWAN_BKC
      INTEGER,SAVE                  :: SWAN_MTC
      INTEGER                       :: UnitNumber
      INTEGER                       :: IFileCounter   !st3 100708:
      INTEGER                       :: UpdateMax(NP)

      LOGICAL                       :: CROSS(4,NP)
#ifdef ADCNETCDF
      LOGICAL                       :: NETCDF_ERROR(14)
#endif

!Casey 090820: Sapphire doesn't like it if these variables
!              are declared as only REAL.  They must be declared
!              as REAL(4) to interface correctly with the
!              SWAN output subroutines.
      REAL(4)                       :: ACLOC(MDC,MSC)
      REAL(4)                       :: DEPXY(NP)
      REAL(4)                       :: FORCE(NP,2)
      REAL(4)                       :: SWAN_CG(MSC)
      REAL(4)                       :: SWAN_NE(MSC)
      REAL(4)                       :: SWAN_NED(MSC)
      REAL(4),ALLOCATABLE           :: SWAN_VOQ(:,:)
      REAL(4)                       :: SWAN_WK(MSC)
      REAL(SZ),POINTER              :: SwanOut(:)
      REAL(SZ),POINTER              :: SwanOut2(:)
C      REAL(SZ),TARGET               :: SwanOut_g(NP_G)
C      REAL(SZ),TARGET               :: SwanOut2_g(NP_G)
      REAL(4)                       :: XC(NP)
      REAL(4)                       :: YC(NP)

      TYPE(OutputDataDescript_t),POINTER    :: SwanDescript
      TYPE(OutputDataDescript_t),POINTER    :: SwanDescriptMax

Cobell 20120510: Added complete list of SWAN Output arrays
      type(OutputDataDescript_t), SAVE, TARGET :: SwanHSDescript
      type(OutputDataDescript_t), SAVE, TARGET :: SwanDIRDescript
      type(OutputDataDescript_t), SAVE, TARGET :: SwanTM01Descript
      type(OutputDataDescript_t), SAVE, TARGET :: SwanTPSDescript
      type(OutputDataDescript_t), SAVE, TARGET :: SwanWindDescript
      type(OutputDataDescript_t), SAVE, TARGET :: SwanTM02Descript
      type(OutputDataDescript_t), SAVE, TARGET :: SwanTMM10Descript
      type(OutputDataDescript_t), SAVE, TARGET :: SwanHSMaxDescript
      type(OutputDataDescript_t), SAVE, TARGET :: SwanDIRMaxDescript
      type(OutputDataDescript_t), SAVE, TARGET :: SwanTM01MaxDescript
      type(OutputDataDescript_t), SAVE, TARGET :: SwanTPSMaxDescript
      type(OutputDataDescript_t), SAVE, TARGET :: SwanWindMaxDescript
      type(OutputDataDescript_t), SAVE, TARGET :: SwanTM02MaxDescript
      type(OutputDataDescript_t), SAVE, TARGET :: SwanTMM10MaxDescript

Cobell 20120510: Added the global arrays used for SWAN globalio
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: SW_HSOut_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: SW_DIROut_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: SW_TM01Out_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: SW_TPSOut_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: SW_WINDXOut_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: SW_WINDYOut_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: SW_TM02Out_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: SW_TMM10Out_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: SW_HSMaxOut_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: SW_DIRMaxOut_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: SW_TM01MaxOut_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: SW_TPSMaxOut_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: SW_WINDMaxOut_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: SW_TM02MaxOut_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: SW_TMM10MaxOut_g(:)

      call setMessageSource("swanoutput")
#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

!... Initialize the Swan output variables.
      DO IO=1,NMOVAR
         SWAN_OQPROC(IO) = .FALSE.
      ENDDO
      NumSwanOutput = 0
!... For significant wave heights (HS).
      NumSwanOutput = NumSwanOutput + 1
      OVEXCV(10) = -99999.
      SWAN_OQPROC(10) = .TRUE.
      SWAN_VOQR(10) = 7 + NumSwanOutput
!... For mean wave directions (DIR).
      NumSwanOutput = NumSwanOutput + 1
      OVEXCV(13) = -99999.
      SWAN_OQPROC(13) = .TRUE.
      SWAN_VOQR(13) = 7 + NumSwanOutput
!... For mean wave periods (TM01).
      NumSwanOutput = NumSwanOutput + 1
      OVEXCV(11) = -99999.
      SWAN_OQPROC(11) = .TRUE.
      SWAN_VOQR(11) = 7 + NumSwanOutput
!... For peak wave periods (TPS).
      NumSwanOutput = NumSwanOutput + 1
      OVEXCV(53) = -99999.
      SWAN_OQPROC(53) = .TRUE.
      SWAN_VOQR(53) = 7 + NumSwanOutput
!... For wind speeds (WX2 and WY2).
      NumSwanOutput = NumSwanOutput + 1
      OVEXCV(26) = 0.
      SWAN_OQPROC(26) = .TRUE.
      SWAN_VOQR(26) = 7 + NumSwanOutput
!... For mean wave periods (TM02).
      NumSwanOutput = NumSwanOutput + 1
      OVEXCV(32) = -99999.
      SWAN_OQPROC(32) = .TRUE.
      SWAN_VOQR(32) = 7 + NumSwanOutput + 1
!... For mean wave periods (TMM10).
      NumSwanOutput = NumSwanOutput + 1
      OVEXCV(47) = -99999.
      SWAN_OQPROC(47) = .TRUE.
      SWAN_VOQR(47) = 7 + NumSwanOutput + 1
      COSCQ = COS(0.)
      SINCQ = SIN(0.)
      IF(.NOT.ALLOCATED(IGSW))  ALLOCATE(IGSW(1:NumSwanOutput))
      IF(.NOT.ALLOCATED(Names)) ALLOCATE(Names(1:NumSwanOutput))
      Names(1) = "HS"
      Names(2) = "DIR"
      Names(3) = "TM01"
      Names(4) = "TPS"
      Names(5) = "WIND"
      Names(6) = "TM02"
      Names(7) = "TMM10"
      IF(.NOT.ALLOCATED(FileName)) ALLOCATE(FileName(1:NumSwanOutput))
      FileName(1) = "swan_"//TRIM(Names(1))//".63"//"       "
      FileName(2) = "swan_"//TRIM(Names(2))//".63"//"       "
      FileName(3) = "swan_"//TRIM(Names(3))//".63"//"       "
      FileName(4) = "swan_"//TRIM(Names(4))//".63"//"       "
      FileName(5) = "swan_"//TRIM(Names(5))//".64"//"       "
      FileName(6) = "swan_"//TRIM(Names(6))//".63"//"       "
      FileName(7) = "swan_"//TRIM(Names(7))//".63"//"       "
      IF(.NOT.Processed26)THEN


Cobell 20120510: At first pass, we set up the output data arrays
C.....Arrange on/off for output into array
        SWAN_OutputAgg(1) = SWAN_OutputHS
        SWAN_OutputAgg(2) = SWAN_OutputDIR
        SWAN_OutputAgg(3) = SWAN_OutputTM01
        SWAN_OutputAgg(4) = SWAN_OutputTPS
        SWAN_OutputAgg(5) = SWAN_OutputWIND
        SWAN_OutputAgg(6) = SWAN_OutputTM02
        SWAN_OutputAgg(7) = SWAN_OutputTMM10

C.......Allocate what we need
        IF(SWAN_OutputHS)THEN
            ALLOCATE(SW_HSOut(NP))
            ALLOCATE(SW_HSMaxOut(NP))
        ENDIF    
        IF(SWAN_OutputDIR)THEN
            ALLOCATE(SW_DIROut(NP))
            ALLOCATE(SW_DIRMaxOut(NP))
        ENDIF    
        IF(SWAN_OutputTM01)THEN
            ALLOCATE(SW_TM01Out(NP))
            ALLOCATE(SW_TM01MaxOut(NP))
        ENDIF    
        IF(SWAN_OutputTPS)THEN
            ALLOCATE(SW_TPSOut(NP))
            ALLOCATE(SW_TPSMaxOut(NP))
        ENDIF    
        IF(SWAN_OutputWIND)THEN
            ALLOCATE(SW_WindXOut(NP), SW_WindYOut(NP))
            ALLOCATE(SW_WindMaxOut(NP))
        ENDIF    
        IF(SWAN_OutputTM02)THEN
            ALLOCATE(SW_TM02Out(NP))
            ALLOCATE(SW_TM02MaxOut(NP))
        ENDIF    
        IF(SWAN_OutputTMM10)THEN
            ALLOCATE(SW_TMM10Out(NP))
            ALLOCATE(SW_TMM10MaxOut(NP))
        ENDIF    

C.......Processor 0 gets the global arrays
        IF( (MNPROC.GT.1).AND.(MYPROC.EQ.0))THEN
            IF(SWAN_OutputHS)THEN
                ALLOCATE(SW_HSOut_g(NP_G))
                ALLOCATE(SW_HSMaxOut_g(NP_G))
            ENDIF    
            IF(SWAN_OutputDIR)THEN
                ALLOCATE(SW_DIROut_g(NP_G))
                ALLOCATE(SW_DIRMaxOut_g(NP_G))
            ENDIF    
            IF(SWAN_OutputTPS)THEN
                ALLOCATE(SW_TPSOut_g(NP_G))
                ALLOCATE(SW_TPSMaxOut_g(NP_G))
            ENDIF
            IF(SWAN_OutputTM01)THEN
                ALLOCATE(SW_TM01Out_g(NP_G))
                ALLOCATE(SW_TM01MaxOut_g(NP_G))
            ENDIF    
            IF(SWAN_OutputWIND)THEN
                ALLOCATE(SW_WindXOut_g(NP_G), SW_WindYOut_g(NP_G))
                ALLOCATE(SW_WindMaxOut_g(NP_G))
            ENDIF    
            IF(SWAN_OutputTM02)THEN
                ALLOCATE(SW_TM02Out_g(NP_G))
                ALLOCATE(SW_TM02MaxOut_g(NP_G))
            ENDIF    
            IF(SWAN_OutputTMM10)THEN
                ALLOCATE(SW_TMM10Out_g(NP_G))
                ALLOCATE(SW_TMM10MaxOut_g(NP_G))
            ENDIF    
        ENDIF

C.......Now prepare the descript arrays on all procs

        IF(SWAN_OutputHS)THEN
C........Significant Wave Height
            SwanHSDescript % specifier            = NOUTGW
            SwanHSDescript % initial_value        = 0.0
            SwanHSDescript % num_items_per_record = 1
            SwanHSDescript % num_fd_records       = NP_G
            SwanHSDescript % num_records_this     = NP
            SwanHSDescript % imap                 => NODES_LG
            SwanHSDescript % array                => SW_HSOut
            SwanHSDescript % array_g              => SW_HSOut_g
            SwanHSDescript % ConsiderWetDry       = .FALSE.
            SwanHSDescript % alternate_value      = -99999.0
            SwanHSDescript % field_name           = "swan_HS"
            SwanHSDescript % file_name            = "swan_HS.63"
            
            SwanHSMaxDescript % specifier            = NOUTGW
            SwanHSMaxDescript % initial_value        = 0.0
            SwanHSMaxDescript % num_items_per_record = 1
            SwanHSMaxDescript % num_fd_records       = NP_G
            SwanHSMaxDescript % num_records_this     = NP
            SwanHSMaxDescript % imap                 => NODES_LG
            SwanHSMaxDescript % array                => SW_HSMaxOut
            SwanHSMaxDescript % array_g              => SW_HSMaxOut_g
            SwanHSMaxDescript % ConsiderWetDry       = .FALSE.
            SwanHSMaxDescript % alternate_value      = -99999.0
            SwanHSMaxDescript % field_name           = "swan_HS_max"
            SwanHSMaxDescript % file_name            = "swan_HS_max.63"
         ENDIF

         IF(SWAN_OutputDIR)THEN
C.........Mean Wave Direction
            SwanDIRDescript % specifier            = NOUTGW
            SwanDIRDescript % initial_value        = 0.0
            SwanDIRDescript % num_items_per_record = 1
            SwanDIRDescript % num_fd_records       = NP_G
            SwanDIRDescript % num_records_this     = NP
            SwanDIRDescript % imap                 => NODES_LG
            SwanDIRDescript % array                => SW_DIROut
            SwanDIRDescript % array_g              => SW_DIROut_g
            SwanDIRDescript % ConsiderWetDry       = .FALSE.
            SwanDIRDescript % alternate_value      = -99999.0
            SwanDIRDescript % field_name           = "swan_DIR"
            SwanDIRDescript % file_name            = "swan_DIR.63"
            
            SwanDIRMaxDescript % specifier            = NOUTGW
            SwanDIRMaxDescript % initial_value        = 0.0
            SwanDIRMaxDescript % num_items_per_record = 1
            SwanDIRMaxDescript % num_fd_records       = NP_G
            SwanDIRMaxDescript % num_records_this     = NP
            SwanDIRMaxDescript % imap                 => NODES_LG
            SwanDIRMaxDescript % array                => SW_DIRMaxOut
            SwanDIRMaxDescript % array_g              => SW_DIRMaxOut_g
            SwanDIRMaxDescript % ConsiderWetDry       = .FALSE.
            SwanDIRMaxDescript % alternate_value      = -99999.0
            SwanDIRMaxDescript % field_name           = "swan_DIR_max"
            SwanDIRMaxDescript % file_name           = "swan_DIR_max.63"
        ENDIF

        IF(SWAN_OutputTM01)THEN
C.........Mean Wave Period (TM01)
           SwanTM01Descript % specifier            = NOUTGW
           SwanTM01Descript % initial_value        = 0.0
           SwanTM01Descript % num_items_per_record = 1
           SwanTM01Descript % num_fd_records       = NP_G
           SwanTM01Descript % num_records_this     = NP
           SwanTM01Descript % imap                 => NODES_LG
           SwanTM01Descript % array                => SW_TM01Out
           SwanTM01Descript % array_g              => SW_TM01Out_g
           SwanTM01Descript % ConsiderWetDry       = .FALSE.
           SwanTM01Descript % alternate_value      = -99999.0
           SwanTM01Descript % field_name           = "swan_TM01"
           SwanTM01Descript % file_name            = "swan_TM01.63"
           
           SwanTM01MaxDescript % specifier            = NOUTGW
           SwanTM01MaxDescript % initial_value        = 0.0
           SwanTM01MaxDescript % num_items_per_record = 1
           SwanTM01MaxDescript % num_fd_records       = NP_G
           SwanTM01MaxDescript % num_records_this     = NP
           SwanTM01MaxDescript % imap                 => NODES_LG
           SwanTM01MaxDescript % array                => SW_TM01MaxOut
           SwanTM01MaxDescript % array_g              => SW_TM01MaxOut_g
           SwanTM01MaxDescript % ConsiderWetDry       = .FALSE.
           SwanTM01MaxDescript % alternate_value      = -99999.0
           SwanTM01MaxDescript % field_name           = "swan_TM01_max"
           SwanTM01MaxDescript % file_name          = "swan_TM01_max.63"
         ENDIF

         IF(SWAN_OutputTPS)THEN
C.........Peak Wave Period (TPS)
            SwanTPSDescript % specifier            = NOUTGW
            SwanTPSDescript % initial_value        = 0.0
            SwanTPSDescript % num_items_per_record = 1
            SwanTPSDescript % num_fd_records       = NP_G
            SwanTPSDescript % num_records_this     = NP
            SwanTPSDescript % imap                 => NODES_LG
            SwanTPSDescript % array                => SW_TPSOut
            SwanTPSDescript % array_g              => SW_TPSOut_g
            SwanTPSDescript % ConsiderWetDry       = .FALSE.
            SwanTPSDescript % alternate_value      = -99999.0
            SwanTPSDescript % field_name           = "swan_TPS"
            SwanTPSDescript % file_name            = "swan_TPS.63"
            
            SwanTPSMaxDescript % specifier            = NOUTGW
            SwanTPSMaxDescript % initial_value        = 0.0
            SwanTPSMaxDescript % num_items_per_record = 1
            SwanTPSMaxDescript % num_fd_records       = NP_G
            SwanTPSMaxDescript % num_records_this     = NP
            SwanTPSMaxDescript % imap                 => NODES_LG
            SwanTPSMaxDescript % array                => SW_TPSMaxOut
            SwanTPSMaxDescript % array_g              => SW_TPSMaxOut_g
            SwanTPSMaxDescript % ConsiderWetDry       = .FALSE.
            SwanTPSMaxDescript % alternate_value      = -99999.0
            SwanTPSMaxDescript % field_name           = "swan_TPS_max"
            SwanTPSMaxDescript % file_name           = "swan_TPS_max.63"
         ENDIF

         IF(SWAN_OutputWIND)THEN
C.........SWAN Wind Values (WINDX,WINDY)
            SwanWindDescript % specifier            = NOUTGW
            SwanWindDescript % initial_value        = 0.0
            SwanWindDescript % num_items_per_record = 2
            SwanWindDescript % num_fd_records       = NP_G
            SwanWindDescript % num_records_this     = NP
            SwanWindDescript % imap                 => NODES_LG
            SwanWindDescript % array                => SW_WindXOut
            SwanWindDescript % array_g              => SW_WindXOut_g
            SwanWindDescript % array2               => SW_WindYOut
            SwanWindDescript % array2_g             => SW_WindYOut_g
            SwanWindDescript % ConsiderWetDry       = .FALSE.
            SwanWindDescript % alternate_value      = -99999.0
            SwanWindDescript % field_name           = "swan_WIND"
            SwanWindDescript % file_name            = "swan_WIND.64"
            
            SwanWindMaxDescript % specifier            = NOUTGW
            SwanWindMaxDescript % initial_value        = 0.0
            SwanWindMaxDescript % num_items_per_record = 1
            SwanWindMaxDescript % num_fd_records       = NP_G
            SwanWindMaxDescript % num_records_this     = NP
            SwanWindMaxDescript % imap                 => NODES_LG
            SwanWindMaxDescript % array                => SW_WindMaxOut
            SwanWindMaxDescript % array_g             => SW_WindMaxOut_g
            SwanWindMaxDescript % ConsiderWetDry       = .FALSE.
            SwanWindMaxDescript % alternate_value      = -99999.0
            SwanWindMaxDescript % field_name           = "swan_WIND_max"
            SwanWindMaxDescript % file_name         = "swan_WIND_max.63"
         ENDIF

         IF(SWAN_OutputTM02)THEN
C.........Mean Wave Period (TM02)
            SwanTM02Descript % specifier            = NOUTGW
            SwanTM02Descript % initial_value        = 0.0
            SwanTM02Descript % num_items_per_record = 1
            SwanTM02Descript % num_fd_records       = NP_G
            SwanTM02Descript % num_records_this     = NP
            SwanTM02Descript % imap                 => NODES_LG
            SwanTM02Descript % array                => SW_TM02Out
            SwanTM02Descript % array_g              => SW_TM02Out_g
            SwanTM02Descript % ConsiderWetDry       = .FALSE.
            SwanTM02Descript % alternate_value      = -99999.0
            SwanTM02Descript % field_name           = "swan_TM02"
            SwanTM02Descript % file_name            = "swan_TM02.63"
            
            SwanTM02MaxDescript % specifier            = NOUTGW
            SwanTM02MaxDescript % initial_value        = 0.0
            SwanTM02MaxDescript % num_items_per_record = 1
            SwanTM02MaxDescript % num_fd_records       = NP_G
            SwanTM02MaxDescript % num_records_this     = NP
            SwanTM02MaxDescript % imap                 => NODES_LG
            SwanTM02MaxDescript % array                => SW_TM02MaxOut
            SwanTM02MaxDescript % array_g             => SW_TM02MaxOut_g
            SwanTM02MaxDescript % ConsiderWetDry       = .FALSE.
            SwanTM02MaxDescript % alternate_value      = -99999.0
            SwanTM02MaxDescript % field_name           = "swan_TM02_max"
            SwanTM02MaxDescript % file_name         = "swan_TM02.63_max"
         ENDIF

         IF(SWAN_OutputTMM10)THEN
C.........Mean Wave Period (TMM10)
            SwanTMM10Descript % specifier            = NOUTGW
            SwanTMM10Descript % initial_value        = 0.0
            SwanTMM10Descript % num_items_per_record = 1
            SwanTMM10Descript % num_fd_records       = NP_G
            SwanTMM10Descript % num_records_this     = NP
            SwanTMM10Descript % imap                 => NODES_LG
            SwanTMM10Descript % array                => SW_TMM10Out
            SwanTMM10Descript % array_g              => SW_TMM10Out_g
            SwanTMM10Descript % ConsiderWetDry       = .FALSE.
            SwanTMM10Descript % alternate_value      = -99999.0
            SwanTMM10Descript % field_name           = "swan_TMM10"
            SwanTMM10Descript % file_name            = "swan_TMM10.63"
            
            SwanTMM10MaxDescript % specifier            = NOUTGW
            SwanTMM10MaxDescript % initial_value        = 0.0
            SwanTMM10MaxDescript % num_items_per_record = 1
            SwanTMM10MaxDescript % num_fd_records       = NP_G
            SwanTMM10MaxDescript % num_records_this     = NP
            SwanTMM10MaxDescript % imap                 => NODES_LG
            SwanTMM10MaxDescript % array               => SW_TMM10MaxOut
            SwanTMM10MaxDescript % array_g           => SW_TMM10MaxOut_g
            SwanTMM10MaxDescript % ConsiderWetDry       = .FALSE.
            SwanTMM10MaxDescript % alternate_value      = -99999.0
            SwanTMM10MaxDescript % field_name         = "swan_TMM10_max"
            SwanTMM10MaxDescript % file_name       = "swan_TMM10_max.63"
         ENDIF

         DO IW=1,NumSwanOutput
            IF(.NOT.SWAN_OutputAgg(IW))CYCLE
            UnitNumber = 300 + IW
!... Copy code from cstart.F about initializing the output file.
            IF((ABS(NOUTGW).EQ.1).OR.(ABS(NOUTGW).EQ.4))THEN
               IF(IW.NE.5)THEN
                   CALL OPEN_GBL_FILE(UnitNumber,TRIM(GLOBALDIR)//'/'//TRIM(FileName(IW)),
     &                                NP_G,NP,Header73)
                   IGSW(IW)=2
               ELSE
                   CALL OPEN_GBL_FILE(UnitNumber,TRIM(GLOBALDIR)//'/'//TRIM(FileName(IW)),
     &                                NP_G,NP,Header74)
                   IGSW(IW)=2
               ENDIF
            ENDIF
            IF(ABS(NOUTGW).EQ.2)THEN
               IGSW(IW) = 0
               OPEN(UnitNumber,FILE=TRIM(LOCALDIR)//'/'//TRIM(FileName(IW)),
     &              ACCESS='DIRECT',RECL=NByte)
               IF(NByte.EQ.4)THEN
                  DO I=1,8
                     WRITE(UnitNumber,REC=IGSW(IW)+I) RDES4(I)
                  ENDDO
                  IGSW(IW) = IGSW(IW) + 8
                  DO I=1,6
                     WRITE(UnitNumber,REC=IGSW(IW)+I) RID4(I)
                  ENDDO
                  IGSW(IW) = IGSW(IW) + 6
                  DO I=1,6
                     WRITE(UnitNumber,REC=IGSW(IW)+I) AID4(I)
                  ENDDO
                  IGSW(IW) = IGSW(IW) + 6
               ENDIF
               IF(NByte.EQ.8)THEN
                  DO I=1,4
                     WRITE(UnitNumber,REC=IGSW(IW)+I) RDES4(I)
                  ENDDO
                  IGSW(IW) = IGSW(IW) + 4
                  DO I=1,3
                     WRITE(UnitNumber,REC=IGSW(IW)+I) RID4(I)
                  ENDDO
                  IGSW(IW) = IGSW(IW) + 3
                  DO I=1,3
                     WRITE(UnitNumber,REC=IGSW(IW)+I) AID4(I)
                  ENDDO
                  IGSW(IW) = IGSW(IW) + 3
               ENDIF
               WRITE(UnitNumber,REC=IGSW(IW)+1) NDSETSW
               WRITE(UnitNumber,REC=IGSW(IW)+2) NP
               WRITE(UnitNumber,REC=IGSW(IW)+3) DT*NSPOOLGW
               WRITE(UnitNumber,REC=IGSW(IW)+4) NSPOOLGW
               WRITE(UnitNumber,REC=IGSW(IW)+5) 2
               IGSW(IW) = IGSW(IW) + 5
               CLOSE(UnitNumber)
            ENDIF
Cobell 20120510: Added for NetCDF SWAN output processing
#ifdef ADCNETCDF
            IF (myProc.eq.0) THEN ! jgf50.60.10: Only proc0 in parallel
               IF((ABS(NOUTGW).EQ.3).OR.(ABS(NOUTGW).EQ.5))THEN
                   NETCDF_ERROR(:) = .FALSE.
                   SELECT CASE(IW)
                   CASE(1)
                       CALL initNetCDFOutputFile(UnitNumber,
     &                  SwanHSDescript,NETCDF_ERROR(1))
                       CALL initNetCDFOutputFile(UnitNumber+15,
     &                  SwanHSMaxDescript,NETCDF_ERROR(2))
                   CASE(2)
                       CALL initNetCDFOutputFile(UnitNumber,
     &                  SwanDIRDescript,NETCDF_ERROR(3))
                       CALL initNetCDFOutputFile(UnitNumber+15,
     &                  SwanDIRMaxDescript,NETCDF_ERROR(4))
                   CASE(3)
                       CALL initNetCDFOutputFile(UnitNumber,
     &                  SwanTM01Descript,NETCDF_ERROR(5))
                       CALL initNetCDFOutputFile(UnitNumber+15,
     &                   SwanTM01MaxDescript,NETCDF_ERROR(6))
                   CASE(4)
                       CALL initNetCDFOutputFile(UnitNumber,
     &                  SwanTPSDescript,NETCDF_ERROR(7))
                       CALL initNetCDFOutputFile(UnitNumber+15,
     &                  SwanTPSMaxDescript,NETCDF_ERROR(8))
                   CASE(5)
                       CALL initNetCDFOutputFile(UnitNumber,
     &                  SwanWINDDescript,NETCDF_ERROR(9))
                       CALL initNetCDFOutputFile(UnitNumber+15,
     &                  SwanWINDMaxDescript,NETCDF_ERROR(10))
                   CASE(6)
                       CALL initNetCDFOutputFile(UnitNumber,
     &                  SwanTM02Descript,NETCDF_ERROR(11))
                       CALL initNetCDFOutputFile(UnitNumber+15,
     &                  SwanTM02MaxDescript,NETCDF_ERROR(12))
                   CASE(7)
                       CALL initNetCDFOutputFile(UnitNumber,
     &                  SwanTMM10Descript,NETCDF_ERROR(13))
                       CALL initNetCDFOutputFile(UnitNumber+15,
     &                  SwanTMM10MaxDescript,NETCDF_ERROR(14))
                   CASE DEFAULT
                       !...Should not be reachable
                   END SELECT
                   DO I=1,7
                       IF(NETCDF_ERROR(I))THEN
#ifdef CMPI
                           subDomainFatalError=.TRUE.
                           CALL MSG_FINI()
#endif                           
                           STOP
                       ENDIF
                   ENDDO
               ENDIF
            ENDIF
#endif
         ENDDO
         SWAN_MTC = MTC
         Processed26 = .TRUE.
!... Ready the array for the maximum SWAN output.
         IF(.NOT.ALLOCATED(SWAN_MAX)) ALLOCATE(SWAN_MAX(1:NP,1:NumSwanOutput))
         DO IP=1,nverts
            DO IW=1,NumSwanOutput
!Casey 110518: Force the max files to take the correct default values.
               IF(IW.EQ.1) IVTYPE=10
               IF(IW.EQ.2) IVTYPE=13
               IF(IW.EQ.3) IVTYPE=11
               IF(IW.EQ.4) IVTYPE=53
               IF(IW.EQ.5) IVTYPE=26
               IF(IW.EQ.6) IVTYPE=32
               IF(IW.EQ.7) IVTYPE=47
               SWAN_MAX(IP,IW) = OVEXCV(IVTYPE)
            ENDDO
         ENDDO
      ENDIF

!... I don't know what BKC does.  Everything seems to work okay
!... when it is set to two, though.
      SWAN_BKC = 2
!... If Swan needs these values for wave numbers, etc.,
!... then maybe it solves for them inside SWOEXA.
      DO IS=1,MSC
         SWAN_CG(IS) = 0.
         SWAN_NE(IS) = 0.
         SWAN_NED(IS) = 0.
         SWAN_WK(IS) = 0.
      ENDDO
!... Allocate VOQ and set up the first seven entries inside it.
      IF(.NOT.ALLOCATED(SWAN_VOQ)) ALLOCATE(SWAN_VOQ(1:NP,1:(7+NumSwanOutput+1)))
      DO IP=1,nverts
!... X,Y in problem coordinate system.
         SWAN_OQPROC(1) = .TRUE.
         SWAN_OQPROC(2) = .TRUE.
         SWAN_VOQR(1) = 1
         SWAN_VOQR(2) = 2
         SWAN_VOQ(IP,SWAN_VOQR(1)) = xcugrd(IP)
         SWAN_VOQ(IP,SWAN_VOQR(2)) = ycugrd(IP)
!... X,Y in adjusted coordinate system.
         SWAN_OQPROC(24) = .TRUE.
         SWAN_OQPROC(25) = .TRUE.
         SWAN_VOQR(24) = 3
         SWAN_VOQR(25) = 4
         SWAN_VOQ(IP,SWAN_VOQR(24)) = 0.
         SWAN_VOQ(IP,SWAN_VOQR(25)) = 0.
!... Depth of water.
         SWAN_OQPROC(4) = .TRUE.
         SWAN_VOQR(4) = 5
         SWAN_VOQ(IP,SWAN_VOQR(4)) = REAL(DP(IP)) + REAL(ETA2(IP))
!... Water current velocities.
         SWAN_OQPROC(5) = .TRUE.
         SWAN_VOQR(5) = 6
         SWAN_VOQ(IP,SWAN_VOQR(5))   = REAL(UU2(IP))
         SWAN_VOQ(IP,SWAN_VOQR(5)+1) = REAL(VV2(IP))
!... Assemble other arrays?
         DEPXY(IP) = REAL(DP(IP)) + REAL(ETA2(IP))
         XC(IP) = xcugrd(IP)
         YC(IP) = ycugrd(IP)
         FORCE(IP,1) = 0.
         FORCE(IP,2) = 0.
         IONOD(IP) = -999
      ENDDO

!... Call the Swan subroutine to interpolate internal Swan quantities.
      CALL SWOEXD(SWAN_OQPROC, nverts, XC, YC, SWAN_VOQR,
     &            SWAN_VOQ, COMPDA, KGRPNT, FORCE, CROSS, IONOD, -999)

!... Call the Swan subroutine to compute the output quantities.
      CALL SWOEXA(SWAN_OQPROC, SWAN_BKC, nverts, XC, YC,
     &            SWAN_VOQR, SWAN_VOQ, AC2, ACLOC, SPCSIG,
     &            SWAN_WK, SWAN_CG, SPCDIR, SWAN_NE, SWAN_NED,
     &            KGRPNT, DEPXY, CROSS)

      UpdateMax = 0
      IFileCounter = 0        !st3 100708: counter of file number
      DO IW=1,NumSwanOutput
         IF(.NOT.SWAN_OutputAgg(IW))CYCLE
         UnitNumber = 300 + IW
!... Copy code from write_output.F to initialize data type.
Cobell 20120510: Removed this section, since it only needs to be done once
C        Now done for all SWAN arrays found in write_output.F
C
C         SwanDescript % specifier            = NOUTGW
C         SwanDescript % initial_value        = 0.
C         IF(IW.NE.5)THEN
C            SwanDescript % num_items_per_record = 1
C         ELSE
C            SwanDescript % num_items_per_record = 2
C         ENDIF
C         SwanDescript % num_fd_records       = NP_G
C         SwanDescript % num_records_this     = NP
C         SwanDescript % imap                 => NODES_LG
C         SwanDescript % array                => SwanOut
C         SwanDescript % array_g              => SwanOut_g
C         SwanDescript % array2               => SwanOut2
C         SwanDescript % array2_g             => SwanOut2_g
!... Assign the output variable to our data structure.
         IF(IW.EQ.1) IVTYPE=10
         IF(IW.EQ.2) IVTYPE=13
         IF(IW.EQ.3) IVTYPE=11
         IF(IW.EQ.4) IVTYPE=53
         IF(IW.EQ.5) IVTYPE=26
         IF(IW.EQ.6) IVTYPE=32
         IF(IW.EQ.7) IVTYPE=47
         DO IP=1,NP
!Casey 090615: On Ranger, this information was not being assigned correctly.
!              Force it to assign by writing it to a temporary string.
Cobell 20120510: Remove the use of SWANOUT and use the newly developed
C     SWAN specific arrays
C     To make the rest of the previously used code function
C     smoothly, SwanDescript is set as a pointer. This pointer
C     will point to the correct OutputDataDescript_t array
C     for the swan variable we are in the loop for. Done to avoid
C     making the code complex for no reason.
C
C
C            SwanOut(IP) = DBLE(SWAN_VOQ(IP,SWAN_VOQR(IVTYPE)))
C            WRITE(UNIT=TempC,FMT=*) SwanOut(IP)
C            IF(IW.EQ.5)THEN
C               SwanOut2(IP) = DBLE(SWAN_VOQ(IP,SWAN_VOQR(IVTYPE)+1))
C               WRITE(UNIT=TempC,FMT=*) SwanOut2(IP)
C            ENDIF

            SELECT CASE(IW)
                CASE(1)
                    SW_HSOut(IP) = DBLE(SWAN_VOQ(IP,SWAN_VOQR(IVTYPE)))
                    SwanOut => SW_HSOut
                    SwanDescript => SwanHSDescript
                CASE(2)
                    SW_DIROut(IP) = DBLE(SWAN_VOQ(IP,SWAN_VOQR(IVTYPE)))
                    SwanOut => SW_DIROut
                    SwanDescript => SwanDIRDescript
                CASE(3)
                    SW_TM01Out(IP) =
     &                  DBLE(SWAN_VOQ(IP,SWAN_VOQR(IVTYPE)))
                    SwanOut => SW_TM01Out
                    SwanDescript => SwanTM01Descript
                CASE(4)
                    SW_TPSOut(IP) = DBLE(SWAN_VOQ(IP,SWAN_VOQR(IVTYPE)))
                    SwanOut => SW_TPSOut
                    SwanDescript => SwanTPSDescript
                CASE(5)
                    SW_WindXOut(IP) =
     &                  DBLE(SWAN_VOQ(IP,SWAN_VOQR(IVTYPE)))
                    SW_WindYOut(IP) =
     &                  DBLE(SWAN_VOQ(IP,SWAN_VOQR(IVTYPE)+1)) !jgfdebug
                    SwanOut => SW_WindXOut
                    SwanOut2 => SW_WindYOut
                    SwanDescript => SwanWindDescript
                CASE(6)
                    SW_TM02Out(IP) =
     &                  DBLE(SWAN_VOQ(IP,SWAN_VOQR(IVTYPE)))
                    SwanOut => SW_TM02Out
                    SwanDescript => SwanTM02Descript
                CASE(7)
                    SW_TMM10Out(IP) =
     &                  DBLE(SWAN_VOQ(IP,SWAN_VOQR(IVTYPE)))
                    SwanOut => SW_TMM10Out
                    SwanDescript => SwanTMM10Descript
                CASE DEFAULT
                    !...Huh? I shouldn't be here...
            END SELECT

Cobell 20120510: Preserve Ranger behavior problem fix if desired

#if 1
            WRITE(UNIT=TempC,FMT=*) SwanOut(IP)
            IF(IW.EQ.5)WRITE(UNIT=TempC,FMT=*) SwanOut2(IP)
#endif

Cobell 20120510: To make the rest of the previously used code function
C          smoothly, SwanDescript is set as a pointer. This pointer
C          will point to the correct OutputDataDescript_t array
C          for the swan variable we are in the loop for.


            IF(IW.NE.5)THEN
!Casey 110518: Ensure that dry nodes are written with default values.
!Casey 120522: Correct logic for the files that don't contain significant wave heights.
               IF(NODECODE(IP).EQ.1)THEN
                  IF(IW.EQ.1)THEN
                     IF(SwanOut(IP).GT.SWAN_MAX(IP,IW))THEN
                        UpdateMax(IP) = 1
                        SWAN_MAX(IP,IW) = SwanOut(IP)
                     ENDIF
                  ELSEIF(UpdateMax(IP).EQ.1)THEN
                     SWAN_MAX(IP,IW) = SwanOut(IP)
                  ENDIF
               ENDIF
            ELSE
               IF((NODECODE(IP).EQ.1).AND.
     &            (SQRT(SwanOut(IP)*SwanOut(IP)+SwanOut2(IP)*SwanOut2(IP)).GT.SWAN_MAX(IP,IW)))THEN
                  IF(UpdateMax(IP).EQ.1)THEN
                     SWAN_MAX(IP,IW) = SQRT(SwanOut(IP)*SwanOut(IP)+SwanOut2(IP)*SwanOut2(IP))
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
!... Copy code from write_output.F to write the output file.
         IF((NWS.NE.0).AND.(NOUTGW.NE.0))THEN
            IF((ITIME.GT.NTCYSGW).AND.(ITIME.LE.NTCYFGW))THEN
               IF((NSCOUGW.EQ.0 .OR. IT.EQ.SWAN_MTC).AND.IT.NE.0)THEN
                  IF(IW.NE.5)THEN
#ifdef CMPI
                    IF(MNWPROC.GT.0) THEN
                       IFileCounter = IFileCounter + 1               !st3 100708: writer buffer check
                       IF( IFileCounter > MNWPROC*NUM_BUF_MAX ) THEN !st3 100708: writer buffer check
                         WRITE(6,*) ' ***** WRITER BUFFER ERROR ***' !st3 100708: writer buffer check
                         WRITE(6,*) '     You need more Writer Core' !st3 100708: writer buffer check
                         subDomainFatalError=.TRUE.
                         CALL MSG_FINI()
                         STOP
                       ENDIF
                       !...Writer processor will make determination on format
                       !   (ASCII,Binary,NetCDF,NetCDF4)
                       CALL WRITE_GBL_FILE_THROUGH_WRITER
     &                    (UnitNumber,TRIM(GLOBALDIR)//'/'//
     &                     TRIM(FileName(IW)),SwanDescript,
     &                     DT*ITIME,ITIME,
     &                     Writer_StoreOne, DBLE(OVEXCV(IVTYPE)))
                    ELSE
                      SELECT CASE(ABS(NOUTGW))
                       CASE(4)  !...ASCII Sparse (4)
                          CALL WRITE_GBL_FILE_SKIP_DEFAULT
     &                      (UnitNumber,TRIM(GLOBALDIR)//'/'//TRIM(FileName(IW)),
     &                      SwanDescript,DT*ITIME,ITIME,StoreOne,DBLE(OVEXCV(IVTYPE)))
                       CASE DEFAULT  !...Full/Binary/NetCDF/NetCDF4 (1,2,3,5)
                          CALL WriteOutArray(UnitNumber,DT*ITIME,ITIME,SwanDescript,
     &                      PackOne,UnPackOne,IGSW(IW),FileName(IW))
                       END SELECT  
                    ENDIF
#else
                    SELECT CASE(ABS(NOUTGW))
                    CASE(4) !...ASCII Sparse
                        CALL WRITE_GBL_FILE_SKIP_DEFAULT
     &                      (UnitNumber,TRIM(GLOBALDIR)//'/'//TRIM(FileName(IW)),
     &                      SwanDescript,DT*ITIME,ITIME,StoreOne,DBLE(OVEXCV(IVTYPE)))
                    CASE DEFAULT !...Full/Binary/NetCDF/NetCDF4
                        CALL WriteOutArray(UnitNumber,DT*ITIME,ITIME,SwanDescript,
     &                      PackOne,UnPackOne,IGSW(IW),FileName(IW))
                    END SELECT
#endif
                  ELSE
#ifdef CMPI
                    IF(MNWPROC.GT.0) THEN
                       IFileCounter = IFileCounter + 1               !st3 100708: writer buffer check
                       IF( IFileCounter > MNWPROC*NUM_BUF_MAX ) THEN !st3 100708: writer buffer check
                         WRITE(6,*) ' ***** WRITER BUFFER ERROR ***' !st3 100708: writer buffer check
                         WRITE(6,*) '     You need more Writer Core' !st3 100708: writer buffer check
                         subDomainFatalError=.TRUE.
                         CALL MSG_FINI()
                         STOP
                       ENDIF  
                       !...Let writer processor determinate format
                       CALL WRITE_GBL_FILE_THROUGH_WRITER
     &                    (UnitNumber,TRIM(GLOBALDIR)//'/'//
     &                     TRIM(FileName(IW)),SwanDescript,
     &                     DT*ITIME,ITIME,
     &                     Writer_StoreTwo, DBLE(OVEXCV(IVTYPE)))
                    ELSE
                       SELECT CASE(ABS(NOUTGW))
                       CASE(4) !...ASCII Sparse    
                            CALL WRITE_GBL_FILE_SKIP_DEFAULT
     &                          (UnitNumber,TRIM(GLOBALDIR)//'/'//TRIM(FileName(IW)),
     &                          SwanDescript,DT*ITIME,ITIME,StoreTwo,DBLE(OVEXCV(IVTYPE)))
                       CASE DEFAULT !...ASCII Full, Binary, NetCDF, NetCDF4
                            CALL WriteOutArray(UnitNumber,DT*ITIME,ITIME,SwanDescript,
     &                       PackTwo,UnPackTwo,IGSW(IW),FileName(IW))
                       END SELECT
                    ENDIF
#else
                    SELECT CASE(ABS(NOUTGW))
                    CASE(4) !...ASCII Sparse 
                        CALL WRITE_GBL_FILE_SKIP_DEFAULT
     &                      (UnitNumber,TRIM(GLOBALDIR)//'/'//TRIM(FileName(IW)),
     &                      SwanDescript,DT*ITIME,ITIME,StoreTwo,DBLE(OVEXCV(IVTYPE)))
                    CASE DEFAULT !...ASCII full, Binary, NetCDF, NetCDF4
                        CALL WriteOutArray(UnitNumber,DT*ITIME,ITIME,SwanDescript,
     &                      PackTwo,UnPackTwo,IGSW(IW),FileName(IW))
                    END SELECT
#endif
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO

#ifdef CMPI
      CALL FLUSH_WRITERS
#endif

      IF(IT.EQ.SWAN_MTC)THEN
         IF(.NOT.ALLOCATED(SWAN_MAXTEMP)) ALLOCATE(SWAN_MAXTEMP(1:NP))
         IF(.NOT.ALLOCATED(SWAN_MAXTEMP_G).AND.(MYPROC.EQ.0)) 
     &       ALLOCATE(SWAN_MAXTEMP_G(1:NP_G))
!... Copy code from write_output.F to write the global maximum files.
         DO IW=1,NumSwanOutput
            IF(.NOT.Swan_OutputAgg(IW))CYCLE
!            DO IP=1,NP
!                SWAN_MAXTEMP(IP) = SWAN_MAX(IP,IW)
!            ENDDO
            UnitNumber = 315 + IW
            SELECT CASE(IW)
                CASE(1)
                    SwanDescriptMax => SwanHSMaxDescript
                    SW_HSMaxOut(:) = SWAN_MAX(:,IW)
                CASE(2)
                    SwanDescriptMax => SwanDirMaxDescript
                    SW_DIRMaxOut(:) = SWAN_MAX(:,IW)
                CASE(3)
                    SwanDescriptMax => SwanTM01MaxDescript
                    SW_TM01MaxOut(:) = SWAN_MAX(:,IW)
                CASE(4)
                    SwanDescriptMax => SwanTPSMaxDescript
                    SW_TPSMaxOut(:) = SWAN_MAX(:,IW)
                CASE(5)
                    SwanDescriptMax => SwanWindMaxDescript
                    SW_WINDMaxOut(:) = SWAN_MAX(:,IW)
                CASE(6)
                    SwanDescriptMax => SwanTM02MaxDescript
                    SW_TM02MaxOut(:) = SWAN_MAX(:,IW)
                CASE(7)
                    SwanDescriptMax => SwanTMM10MaxDescript
                    SW_TMM10MaxOut(:) = SWAN_MAX(:,IW)
                CASE DEFAULT
                    !...Can't get here
            END SELECT        
            SELECT CASE(ABS(NOUTGW)) 
            CASE(3,5)
                SwanDescriptMax%Specifier = NOUTGW
                CALL writeOutArrayMinMax(UnitNumber,DT*ITIME,ITIME,1,
     &              SwanDescriptMax,PackOne,UnpackOne)
            CASE DEFAULT
                !tcm v51.20.01 changed unit number from 963 to 763 to avoid conflicts with STWAVE
                FileNameMax = "swan_"//TRIM(Names(IW))//"_max.63"
                CALL OPEN_MINMAX_FILE(763, TRIM(GLOBALDIR)//'/'//TRIM(FileNameMax),
     &              NP_G, NP, 1, HEADER_MAX)
                CALL WRITE_GBL_FILE(763, TRIM(GLOBALDIR)//'/'//TRIM(FileNameMax),
     &              SwanDescriptMax, DT*ITIME, ITIME, StoreOne)
            END SELECT    
         ENDDO
         IF(ALLOCATED(SWAN_MAXTEMP)) DEALLOCATE(SWAN_MAXTEMP)
      ENDIF

!Casey 090620: More deallocations.
      IF(ALLOCATED(FileName)) DEALLOCATE(FileName)
      IF(ALLOCATED(IGSW))     DEALLOCATE(IGSW)
      IF(ALLOCATED(Names))    DEALLOCATE(Names)
      IF(ALLOCATED(SWAN_VOQ)) DEALLOCATE(SWAN_VOQ)

#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE SwanOutput
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E    P A D C S W A N _ I N I T
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE PADCSWAN_INIT

      USE GLOBAL,   ONLY: DT,
     &                    ETA2,
     &                    NP,
     &                    WVNX1,
     &                    WVNX2,
     &                    WVNY1,
     &                    WVNY2,
Casey 120305: Merging ice changes from Taylor Asher.
CTGA 110912:  Full SWAN ice implementation.  Use global variables
     &                    DP,
     &                    NCICE,
     &                    CICE1,
     &                    CICE2
CTGA 111003:  Adding ability to read in ice info from fort.26 file
      USE SWCOMM3,  ONLY: SWAN_WBICETH => WBICETH,
     &                    SWAN_IICE => IICE
      USE TIMECOMM, ONLY: SWAN_DT => DT

      IMPLICIT NONE

      INTRINSIC :: ALLOCATED

      INTEGER   :: IN  ! Node counter.

      call setMessageSource("padcswan_init")
#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

Casey 090302: Allocate memory for the ADCIRC values that may be
C             passed to SWAN.
      IF(.NOT.ALLOCATED(SWAN_ETA2)) ALLOCATE(SWAN_ETA2(1:NP,1:2))
      IF(.NOT.ALLOCATED(SWAN_UU2))  ALLOCATE(SWAN_UU2(1:NP,1:2))
      IF(.NOT.ALLOCATED(SWAN_VV2))  ALLOCATE(SWAN_VV2(1:NP,1:2))
      IF(.NOT.ALLOCATED(SWAN_WX2))  ALLOCATE(SWAN_WX2(1:NP,1:2))
      IF(.NOT.ALLOCATED(SWAN_WY2))  ALLOCATE(SWAN_WY2(1:NP,1:2))
Casey 091216: Added allocation for friction coupling.
      IF(.NOT.ALLOCATED(SWAN_Z0))   ALLOCATE(SWAN_Z0(1:NP,1:2))

#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
      IF(.NOT.ALLOCATED(TKXX))     ALLOCATE(TKXX(1:NP))
      IF(.NOT.ALLOCATED(TKXY))     ALLOCATE(TKXY(1:NP))
      IF(.NOT.ALLOCATED(TKYY))     ALLOCATE(TKYY(1:NP))
      IF(.NOT.ALLOCATED(WAVE_A))   ALLOCATE(WAVE_A(1:NP))
      IF(.NOT.ALLOCATED(WAVE_A1))  ALLOCATE(WAVE_A1(1:NP))
      IF(.NOT.ALLOCATED(WAVE_A2))  ALLOCATE(WAVE_A2(1:NP))
      IF(.NOT.ALLOCATED(WAVE_H))   ALLOCATE(WAVE_H(1:NP))
      IF(.NOT.ALLOCATED(WAVE_H1))  ALLOCATE(WAVE_H1(1:NP))
      IF(.NOT.ALLOCATED(WAVE_H2))  ALLOCATE(WAVE_H2(1:NP))
      IF(.NOT.ALLOCATED(WAVE_T))   ALLOCATE(WAVE_T(1:NP))
      IF(.NOT.ALLOCATED(WAVE_T1))  ALLOCATE(WAVE_T1(1:NP))
      IF(.NOT.ALLOCATED(WAVE_T2))  ALLOCATE(WAVE_T2(1:NP))
      TKXX    = 0.D0
      TKXY    = 0.D0
      TKYY    = 0.D0
      WAVE_A  = 0.D0
      WAVE_A1 = 0.D0
      WAVE_A2 = 0.D0
      WAVE_H  = 0.D0
      WAVE_H1 = 0.D0
      WAVE_H2 = 0.D0
      WAVE_T  = 0.D0
      WAVE_T1 = 0.D0
      WAVE_T2 = 0.D0
#endif

Casey 120302: Merging ice changes from Taylor Asher.
CTGA 110912:  Full SWAN ice implementation.  Allocating ice arrays.
      IF(NCICE.NE.0) THEN
         IF(.NOT.ALLOCATED(ICY)) ALLOCATE(ICY(1:NP))
         IF(.NOT.ALLOCATED(SWAN_CICE2)) ALLOCATE(SWAN_CICE2(1:NP,1:2))
         ICY        = 0.D0
         SWAN_CICE2 = 0.D0
      ENDIF

Casey 090302: Initialize the arrays.
      DO IN=1,NP
         SWAN_ETA2(IN,1) = ETA2(IN)
         SWAN_ETA2(IN,2) = ETA2(IN)
         SWAN_UU2 (IN,1) = 0.D0
         SWAN_UU2 (IN,2) = 0.D0
         SWAN_VV2 (IN,1) = 0.D0
         SWAN_VV2 (IN,2) = 0.D0
      ENDDO

Casey 090302: Need to pass correct first wind snaps to SWAN.
      DO IN=1,NP
         SWAN_WX2(IN,1) = WVNX1(IN)
         SWAN_WY2(IN,1) = WVNY1(IN)
         SWAN_WX2(IN,2) = WVNX2(IN)
         SWAN_WY2(IN,2) = WVNY2(IN)
      ENDDO

Casey 091216: Initialize arrays for friction coupling and compute first set of values.
      IF(COUPFRIC)THEN
         CALL Manning2Madsen
         DO IN=1,NP
            SWAN_Z0(IN,1) = SWAN_Z0(IN,2)
         ENDDO
      ENDIF

Casey 120302: Merging ice changes from Taylor Asher.
CTGA 110912:  Full SWAN ice implementation.  Applying ice coverage mask
C             wherever ice concentration reaches/exceeds threshold.
C             Mask is applied by making total water depth (water
C             elevation plus bathymetry) equal to zero.  Same as what
C             is done for dry nodes.
C             The array ICY is used so it can be output in the future.
CTGA111003:  Adding ability to read in ice parameters from fort.26 file
CTGA      IF(NCICE.NE.0) THEN
      IF((NCICE.NE.0).AND.(SWAN_IICE.EQ.2)) THEN
         ICY=0
         DO IN=1,NP
            SWAN_CICE2(IN,1) = CICE1(IN)
            SWAN_CICE2(IN,2) = CICE2(IN)
CTGA 111003:  Adding ability to read in ice parameters from fort.26 file
CTGA            IF (SWAN_CICE2(IN,1).GE.70) ICY(IN) = 1
            IF (SWAN_CICE2(IN,1).GE.SWAN_WBICETH) ICY(IN) = 1
            IF(ICY(IN).EQ.1) THEN
CTGA 110915:  Ensure SWAN considers conditions to be too icy for wave
C             generation for the current timestep by setting both
C             values of SWAN_*(IN,:) to the "dry" values.
               SWAN_ETA2(IN,:) = - DP(IN)
               SWAN_UU2 (IN,:) = 0.D0
               SWAN_VV2 (IN,:) = 0.D0
CTGA110915               SWAN_ETA2(IN,2) = - DP(IN)
CTGA110915               SWAN_UU2 (IN,2) = 0.D0
CTGA110915               SWAN_VV2 (IN,2) = 0.D0
            ENDIF
         ENDDO
      ENDIF

Casey 090302: Allow SWAN to initialize.
      CALL SWINITMPI
      CALL SWMAIN(0,0)

Casey 090302: Advance the wind speeds now so that they do not
C             get written over during the time step loop.
      DO IN=1,NP
         SWAN_WX2(IN,1) = SWAN_WX2(IN,2)
         SWAN_WY2(IN,1) = SWAN_WY2(IN,2)
      ENDDO

Casey 090302: Compute the coupling interval as equivalent
C             to the SWAN time step.
      CouplingInterval = SWAN_DT / DT

#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE PADCSWAN_INIT
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E    P A D C S W A N _ R U N
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE PADCSWAN_RUN(ITIME)

      USE GLOBAL,   ONLY: DP,
     &                    DT,
     &                    ETA2,
     &                    NODECODE,
     &                    NP,
     &                    UU2,
     &                    VV2,
Casey 120302: Merging ice changes from Taylor Asher.
CTGA 110912:  Full SWAN ice implementation.  Use global variables
     &                    STATIM,
     &                    NCICE,
     &                    CICE1,
     &                    CICE2,
     &                    CICE_TIME1,
     &                    CICE_TIMINC
CTGA 111003:  Adding ability to read in ice info from fort.26 file
      USE SWCOMM3,  ONLY: SWAN_IICE => IICE,
     &                    SWAN_WBICETH => WBICETH
      USE TIMECOMM, ONLY: SWAN_DT => DT

      IMPLICIT NONE

      INTRINSIC    :: REAL

      INTEGER      :: IN           ! Node counter.
      INTEGER      :: IT           ! Time step counter.
      INTEGER      :: ITIME        ! Time step from ADCIRC.
      INTEGER,SAVE :: SwanTimeStep ! Counter for SWAN time steps.

Casey 120302: Merging ice changes from Taylor Asher.
CTGA 110912:  Full SWAN ice implementation.  Adding variables.
      REAL(SZ)             :: SWAN_PIC
      REAL(SZ)             :: ICETIMEFRAC

      call setMessageSource("padcswan_run")
#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      DO IN=1,NP
Casey 090302: Move the 'old' values into the first columns.
         SWAN_ETA2(IN,1) = SWAN_ETA2(IN,2)
         SWAN_UU2 (IN,1) = SWAN_UU2 (IN,2)
         SWAN_VV2 (IN,1) = SWAN_VV2 (IN,2)
Casey 090302: Find the 'new' values for the second columns.
         IF(NODECODE(IN).EQ.1)THEN
            SWAN_ETA2(IN,2) = ETA2(IN)
            SWAN_UU2 (IN,2) = UU2 (IN)
            SWAN_VV2 (IN,2) = VV2 (IN)
         ELSE
            SWAN_ETA2(IN,2) = - DP(IN)
            SWAN_UU2 (IN,2) = 0.D0
            SWAN_VV2 (IN,2) = 0.D0
         ENDIF
      ENDDO

Casey 120302: Merging ice changes from Taylor Asher.
CTGA 110912:  Full SWAN ice implementation.  Applying ice coverage mask
C             wherever ice concentration reaches/exceeds threshold.
C             Mask is applied by making total water depth (water
C             elevation plus bathymetry) equal to zero.  Same as what
C             is done for dry nodes.
C             The array ICY is used so it can be output in the future.
CTGA111003:  Adding ability to read in ice parameters from fort.26 file
CTGA      IF(NCICE.NE.0) THEN
      IF((NCICE.NE.0).AND.(SWAN_IICE.EQ.2)) THEN
         ICY=0
         ICETIMEFRAC = ((STATIM*86400+ITIME*DT)-CICE_TIME1)/CICE_TIMINC
         DO IN=1,NP
            SWAN_CICE2(IN,1) = CICE1(IN)
            SWAN_CICE2(IN,2) = CICE2(IN)
            SWAN_PIC = ICETIMEFRAC*(SWAN_CICE2(IN,2) - SWAN_CICE2(IN,1))
     &               + SWAN_CICE2(IN,1)
CTGA 111003:  Adding ability to read in ice parameters from fort.26 file
CTGA            IF (SWAN_PIC.GE.70) ICY(IN) = 1
            IF (SWAN_PIC.GE.SWAN_WBICETH) ICY(IN) = 1
            IF(ICY(IN).EQ.1) THEN
CTGA 110915:  Ensure SWAN considers conditions to be too icy for wave
C             generation for the current timestep by setting both
C             values of SWAN_*(IN,:) to the "dry" values.
               SWAN_ETA2(IN,:) = - DP(IN)
               SWAN_UU2 (IN,:) = 0.D0
               SWAN_VV2 (IN,:) = 0.D0
CTGA110915               SWAN_ETA2(IN,1) = - DP(IN)
CTGA110915               SWAN_UU2 (IN,1) = 0.D0
CTGA110915               SWAN_VV2 (IN,1) = 0.D0
CTGA110915               SWAN_ETA2(IN,2) = - DP(IN)
CTGA110915               SWAN_UU2 (IN,2) = 0.D0
CTGA110915               SWAN_VV2 (IN,2) = 0.D0
            ENDIF
         ENDDO
      ENDIF

Casey 091216: Move the 'old' values into the first column,
C             and then compute the 'new' values.
      IF(COUPFRIC)THEN
         DO IN=1,NP
            SWAN_Z0(IN,1) = SWAN_Z0(IN,2)
         ENDDO
         CALL Manning2Madsen
      ENDIF

      !jgf50.28: Removed DO ... when SWAN coupling was under development,
      ! this was a loop so that SWAN could do multiple time steps per coupling
      ! interval, but it became clear that there was not a need for that and
      ! that new information should be used at the start of each SWAN step.

Casey 090729: Average values over the time step.
      InterpoWeight = 0.5D0
      SwanTimeStep = SwanTimeStep + 1
      CALL SWMAIN(ITIME,SwanTimeStep)

Casey 090302: Move the 'old' values into the first columns.
C             We do this for the winds after the SWAN time steps
C             because the 'new' values are written over in TIMESTEP.
C             We need to move them to the 'old' values now,
C             so that we do not lose them during the ADCIRC time stepping.
      DO IN=1,NP
         SWAN_WX2(IN,1) = SWAN_WX2(IN,2)
         SWAN_WY2(IN,1) = SWAN_WY2(IN,2)
      ENDDO

#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE PADCSWAN_RUN
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E    P A D C S W A N _ F I N A L
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE PADCSWAN_FINAL
      IMPLICIT NONE

      call setMessageSource("padcswan_final")
#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      CALL SWEXITMPI()

      IF(ALLOCATED(SWAN_ETA2)) DEALLOCATE(SWAN_ETA2)
      IF(ALLOCATED(SWAN_UU2))  DEALLOCATE(SWAN_UU2)
      IF(ALLOCATED(SWAN_VV2))  DEALLOCATE(SWAN_VV2)
      IF(ALLOCATED(SWAN_WX2))  DEALLOCATE(SWAN_WX2)
      IF(ALLOCATED(SWAN_WY2))  DEALLOCATE(SWAN_WY2)

Casey 091216: Add deallocation for friction coupling.
      IF(ALLOCATED(SWAN_Z0)) DEALLOCATE(SWAN_Z0)

Casey 120302: Merging ice changes from Taylor Asher.
CTGA 110912:
      IF(ALLOCATED(ICY)) DEALLOCATE(ICY)
      IF(ALLOCATED(SWAN_CICE2)) DEALLOCATE(SWAN_CICE2)

#ifdef CSWANFRIC
Casey 091020: Adopt Ethan's/Joannes's modified friction.
      IF(ALLOCATED(TKXX))    DEALLOCATE(TKXX)
      IF(ALLOCATED(TKXY))    DEALLOCATE(TKXY)
      IF(ALLOCATED(TKYY))    DEALLOCATE(TKYY)
      IF(ALLOCATED(WAVE_A))  DEALLOCATE(WAVE_A)
      IF(ALLOCATED(WAVE_A1)) DEALLOCATE(WAVE_A1)
      IF(ALLOCATED(WAVE_A2)) DEALLOCATE(WAVE_A2)
      IF(ALLOCATED(WAVE_H))  DEALLOCATE(WAVE_H)
      IF(ALLOCATED(WAVE_H1)) DEALLOCATE(WAVE_H1)
      IF(ALLOCATED(WAVE_H2)) DEALLOCATE(WAVE_H2)
      IF(ALLOCATED(WAVE_T))  DEALLOCATE(WAVE_T)
      IF(ALLOCATED(WAVE_T1)) DEALLOCATE(WAVE_T1)
      IF(ALLOCATED(WAVE_T2)) DEALLOCATE(WAVE_T2)
#endif

#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE PADCSWAN_FINAL
C-----------------------------------------------------------------------
#endif

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      END MODULE Couple2Swan
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
