
      MODULE Couple2Swan
#ifdef CSWAN

      USE SIZES,   ONLY: SZ
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

Casey 090302: These arrays contain the radiation stresses.
      REAL    ,ALLOCATABLE :: ADCIRC_SXX(:,:)
      REAL    ,ALLOCATABLE :: ADCIRC_SXY(:,:)
      REAL    ,ALLOCATABLE :: ADCIRC_SYY(:,:)

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

Casey 090302: These variables control SWAN output.
      INTEGER, ALLOCATABLE :: IGSW(:)
      INTEGER              :: NumSwanOutput
      LOGICAL              :: Processed26 = .FALSE.
      LOGICAL              :: SWAN_OQPROC(NMOVAR)
      INTEGER              :: SWAN_VOQR(NMOVAR)

Casey 090303: Add variable for maximum SWAN output.
      REAL(SZ),ALLOCATABLE        :: SWAN_MAX(:,:)
      REAL(SZ),ALLOCATABLE,TARGET :: SWAN_MAXTEMP(:)


      CONTAINS



      SUBROUTINE ComputeRadiationStresses

      USE Couple2Adcirc, ONLY : COMPDA
      USE M_GENARR, ONLY      : AC2, SPCDIR, SPCSIG
      USE SwanGridData, ONLY  : nverts
      USE SWCOMM3, ONLY       : DDIR, DEPMIN, FRINTF, GRAV, JDP2, MDC, MSC, RHO

      IMPLICIT NONE

      INTEGER :: I, ID, IERR, IS
      REAL    :: DEPLOC
      REAL    :: WK(MSC), CG(MSC), NE(MSC), NED(MSC)

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
         DEPLOC = COMPDA(I,JDP2)
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

      END SUBROUTINE ComputeRadiationStresses



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
           TEMP_SXX(IP) = (1.0 - InterpoWeight) * ADCIRC_SXX(IP,1)
     &                 + InterpoWeight * ADCIRC_SXX(IP,2)
           TEMP_SXY(IP) = (1.0 - InterpoWeight) * ADCIRC_SXY(IP,1)
     &                 + InterpoWeight * ADCIRC_SXY(IP,2)
           TEMP_SYY(IP) = (1.0 - InterpoWeight) * ADCIRC_SYY(IP,1)
     &                 + InterpoWeight * ADCIRC_SYY(IP,2)

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

           DSXXDX(IE) = (1.D0/AREAS(IE)) *
     &                ( TEMP_SXX(NM(IE,1)) * (Y(NM(IE,2)) - Y(NM(IE,3)))
     &                + TEMP_SXX(NM(IE,2)) * (Y(NM(IE,3)) - Y(NM(IE,1)))
     &                + TEMP_SXX(NM(IE,3)) * (Y(NM(IE,1)) - Y(NM(IE,2))) )

           DSXYDY(IE) = (1.D0/AREAS(IE)) *
     &                ( TEMP_SXY(NM(IE,1)) * (X(NM(IE,3)) - X(NM(IE,2)))
     &                + TEMP_SXY(NM(IE,2)) * (X(NM(IE,1)) - X(NM(IE,3)))
     &                + TEMP_SXY(NM(IE,3)) * (X(NM(IE,2)) - X(NM(IE,1))) )

           DSXYDX(IE) = (1.D0/AREAS(IE)) *
     &                ( TEMP_SXY(NM(IE,1)) * (Y(NM(IE,2)) - Y(NM(IE,3)))
     &                + TEMP_SXY(NM(IE,2)) * (Y(NM(IE,3)) - Y(NM(IE,1)))
     &                + TEMP_SXY(NM(IE,3)) * (Y(NM(IE,1)) - Y(NM(IE,2))) )

           DSYYDY(IE) = (1.D0/AREAS(IE)) *
     &                ( TEMP_SYY(NM(IE,1)) * (X(NM(IE,3)) - X(NM(IE,2)))
     &                + TEMP_SYY(NM(IE,2)) * (X(NM(IE,1)) - X(NM(IE,3)))
     &                + TEMP_SYY(NM(IE,3)) * (X(NM(IE,2)) - X(NM(IE,1))) )

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

!Casey 080324:
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

!Casey 080324:
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

      ENDIF

      RETURN

      END SUBROUTINE ComputeWaveDrivenForces



      SUBROUTINE SwanOutput(ITIME, IT)

      USE Couple2Adcirc,ONLY: COMPDA
      USE GLOBAL,       ONLY: AID4,
     &                        DP,
     &                        DT,
     &                        ETA2,
     &                        NDSETSW,
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
      USE GLOBAL_IO,    ONLY: Header73,
     &                        Header74,
     &                        HEADER_MAX,
     &                        OPEN_GBL_FILE,
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
     &                        MYPROC,
     &                        NBYTE
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
      USE WRITER,       ONLY: WRITE_GBL_FILE_THROUGH_WRITER,
     &                        Writer_StoreOne,
     &                        Writer_StoreTwo
#endif

      IMPLICIT NONE

      INTRINSIC                     :: ALLOCATED
      INTRINSIC                     :: INDEX
      INTRINSIC                     :: TRIM

      CHARACTER(LEN=15),ALLOCATABLE :: FileName(:)
      CHARACTER(LEN=20)             :: FileNameMax
      CHARACTER(LEN=10),ALLOCATABLE :: Names(:)

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

      LOGICAL                       :: CROSS(4,NP)

      REAL                          :: ACLOC(MDC,MSC)
      REAL                          :: DEPXY(NP)
      REAL                          :: FORCE(NP,2)
      REAL                          :: SWAN_CG(MSC)
      REAL                          :: SWAN_NE(MSC)
      REAL                          :: SWAN_NED(MSC)
      REAL,ALLOCATABLE              :: SWAN_VOQ(:,:)
      REAL                          :: SWAN_WK(MSC)
      REAL(8),TARGET                :: SwanOut(NP)
      REAL(8),TARGET                :: SwanOut2(NP)
      REAL(8),TARGET                :: SwanOut_g(NP_G)
      REAL(8),TARGET                :: SwanOut2_g(NP_G)
      REAL                          :: XC(NP)
      REAL                          :: YC(NP)

      TYPE(OutputDataDescript_t)    :: SwanDescript
      TYPE(OutputDataDescript_t)    :: SwanDescriptMax

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
      COSCQ = COS(0.)
      SINCQ = SIN(0.)
      ALLOCATE(IGSW(1:NumSwanOutput))
      ALLOCATE(Names(1:NumSwanOutput))
      Names(1) = "HS"
      Names(2) = "DIR"
      Names(3) = "TM01"
      Names(4) = "TPS"
      Names(5) = "WIND"
      Names(6) = "TM02"
      ALLOCATE(FileName(1:NumSwanOutput))
      FileName(1) = "swan_"//TRIM(Names(1))//".63"//"       "
      FileName(2) = "swan_"//TRIM(Names(2))//".63"//"       "
      FileName(3) = "swan_"//TRIM(Names(3))//".63"//"       "
      FileName(4) = "swan_"//TRIM(Names(4))//".63"//"       "
      FileName(5) = "swan_"//TRIM(Names(5))//".64"//"       "
      FileName(6) = "swan_"//TRIM(Names(6))//".63"//"       "
      IF(.NOT.Processed26)THEN
         DO IW=1,NumSwanOutput
            IF((IW.NE.3).AND.(IW.NE.5))THEN
               UnitNumber = 300 + IW
!... Copy code from cstart.F about initializing the output file.
               IF((ABS(NOUTGW).EQ.1).OR.(ABS(NOUTGW).EQ.4))THEN
                  IF(IW.NE.5)THEN
                      CALL OPEN_GBL_FILE(UnitNumber,TRIM(GLOBALDIR)//'/'//TRIM(FileName(IW)),
     &                                   NP_G,NP,Header73)
                      IGSW(IW)=2
                  ELSE
                      CALL OPEN_GBL_FILE(UnitNumber,TRIM(GLOBALDIR)//'/'//TRIM(FileName(IW)),
     &                                   NP_G,NP,Header74)
                      IGSW(IW)=2
                  ENDIF
               ENDIF
               IF(ABS(NOUTGW).EQ.2)THEN
                  IGSW(IW) = 0
                  OPEN(UnitNumber,FILE=TRIM(LOCALDIR)//'/'//TRIM(FileName(IW)),
     &                 ACCESS='DIRECT',RECL=NByte)
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
            ENDIF
         ENDDO
         SWAN_MTC = MTC
         Processed26 = .TRUE.
!... Ready the array for the maximum SWAN output.
         IF(.NOT.ALLOCATED(SWAN_MAX)) ALLOCATE(SWAN_MAX(1:NP,1:NumSwanOutput))
         DO IP=1,nverts
            DO IW=1,NumSwanOutput
               SWAN_MAX(IP,IW) = 0.D0
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
      ALLOCATE(SWAN_VOQ(1:NP,1:(7+NumSwanOutput+1)))
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
         SWAN_VOQ(IP,SWAN_VOQR(4)) = DP(IP) + ETA2(IP)
!... Water current velocities.
         SWAN_OQPROC(5) = .TRUE.
         SWAN_VOQR(5) = 6
         SWAN_VOQ(IP,SWAN_VOQR(5))   = UU2(IP)
         SWAN_VOQ(IP,SWAN_VOQR(5)+1) = VV2(IP)
!... Assemble other arrays?
         DEPXY(IP) = DP(IP) + ETA2(IP)
         XC(IP) = xcugrd(IP)
         YC(IP) = ycugrd(IP)
         FORCE(IP,1) = 0.
         FORCE(IP,2) = 0.
         IONOD(IP) = -999
      ENDDO

!... Call the Swan subroutine to interpolate internal Swan quantities.
      CALL SWOEXD(SWAN_OQPROC, nverts, XC, YC, SWAN_VOQR,
     &            SWAN_VOQ, COMPDA, KGRPNT, FORCE, CROSS, IONOD)

!... Call the Swan subroutine to compute the output quantities.
      CALL SWOEXA(SWAN_OQPROC, SWAN_BKC, nverts, XC, YC,
     &            SWAN_VOQR, SWAN_VOQ, AC2, ACLOC, SPCSIG,
     &            SWAN_WK, SWAN_CG, SPCDIR, SWAN_NE, SWAN_NED,
     &            KGRPNT, DEPXY, CROSS)

      DO IW=1,NumSwanOutput
         IF((IW.NE.3).AND.(IW.NE.5))THEN
            UnitNumber = 300 + NumSwanOutput
!... Copy code from write_output.F to initialize data type.
            SwanDescript % specifier            = NOUTGW
            SwanDescript % initial_value        = 0.
            IF(IW.NE.5)THEN
               SwanDescript % num_items_per_record = 1
            ELSE
               SwanDescript % num_items_per_record = 2
            ENDIF
            SwanDescript % num_fd_records       = NP_G
            SwanDescript % num_records_this     = NP
            SwanDescript % imap                 => NODES_LG
            SwanDescript % array                => SwanOut
            SwanDescript % array_g              => SwanOut_g
            SwanDescript % array2               => SwanOut2
            SwanDescript % array2_g             => SwanOut2_g
!... Assign the output variable to our data structure.
            IF(IW.EQ.1) IVTYPE=10
            IF(IW.EQ.2) IVTYPE=13
            IF(IW.EQ.3) IVTYPE=11
            IF(IW.EQ.4) IVTYPE=53
            IF(IW.EQ.5) IVTYPE=26
            IF(IW.EQ.6) IVTYPE=32
            DO IP=1,NP
               SwanOut(IP) = SWAN_VOQ(IP,SWAN_VOQR(IVTYPE))
               IF(IW.EQ.5)THEN
                  SwanOut2(IP) = SWAN_VOQ(IP,SWAN_VOQR(IVTYPE)+1)
               ENDIF
               IF(IW.NE.5)THEN
                  IF(SwanOut(IP).GT.SWAN_MAX(IP,IW))THEN
                     SWAN_MAX(IP,IW) = SwanOut(IP)
                  ENDIF
               ELSE
                  IF(SQRT(SwanOut(IP)*SwanOut(IP)+SwanOut2(IP)*SwanOut2(IP)).GT.SWAN_MAX(IP,IW))THEN
                     SWAN_MAX(IP,IW) = SQRT(SwanOut(IP)*SwanOut(IP)+SwanOut2(IP)*SwanOut2(IP))
                  ENDIF
               ENDIF
            ENDDO
!... Copy code from write_output.F to write the output file.
            IF((NWS.NE.0).AND.(NOUTGW.NE.0))THEN
               IF((ITIME.GT.NTCYSGW).AND.(ITIME.LE.NTCYFGW))THEN
                  IF((NSCOUGW.EQ.0 .OR. IT.EQ.SWAN_MTC).AND.IT.NE.0)THEN
                     IF((IW.NE.3).AND.(IW.NE.5))THEN
                        IF(ABS(NOUTGW).EQ.1.OR.ABS(NOUTGW).EQ.2)THEN
                           CALL WriteOutArray(UnitNumber,DT*ITIME,ITIME,SwanDescript,
     &                                         PackOne,UnPackOne,IGSW(IW),FileName(IW))
                        ELSEIF(ABS(NOUTGW).EQ.4)THEN
#ifdef CMPI
                           IF(MNWPROC.GT.0) THEN
                             CALL WRITE_GBL_FILE_THROUGH_WRITER
     &                         (TRIM(GLOBALDIR)//'/'//TRIM(FileName(IW)),SwanDescript,
     &                          DT*ITIME,ITIME,
     &                          Writer_StoreOne, DBLE(OVEXCV(IVTYPE)))
                           ELSE
                             CALL WRITE_GBL_FILE_SKIP_DEFAULT
     &                         (UnitNumber,TRIM(GLOBALDIR)//'/'//TRIM(FileName(IW)),
     &                          SwanDescript,DT*ITIME,ITIME,StoreOne,DBLE(OVEXCV(IVTYPE)))
                           ENDIF
#else
                           CALL WRITE_GBL_FILE_SKIP_DEFAULT
     &                        (UnitNumber,TRIM(GLOBALDIR)//'/'//TRIM(FileName(IW)),
     &                           SwanDescript,DT*ITIME,ITIME,StoreOne,DBLE(OVEXCV(IVTYPE)))
#endif
                        ENDIF
                     ELSE
                        IF(ABS(NOUTGW).EQ.1.OR.ABS(NOUTGW).EQ.2)THEN
                           CALL WriteOutArray(UnitNumber,DT*ITIME,ITIME,SwanDescript,
     &                                         PackTwo,UnPackTwo,IGSW(IW),FileName(IW))
                        ELSEIF(ABS(NOUTGW).EQ.4)THEN
#ifdef CMPI
                           IF(MNWPROC.GT.0) THEN
                             CALL WRITE_GBL_FILE_THROUGH_WRITER
     &                         (TRIM(GLOBALDIR)//'/'//TRIM(FileName(IW)),SwanDescript,
     &                          DT*ITIME,ITIME,
     &                          Writer_StoreTwo, DBLE(OVEXCV(IVTYPE)))
                           ELSE
                             CALL WRITE_GBL_FILE_SKIP_DEFAULT
     &                         (UnitNumber,TRIM(GLOBALDIR)//'/'//TRIM(FileName(IW)),
     &                          SwanDescript,DT*ITIME,ITIME,StoreTwo,DBLE(OVEXCV(IVTYPE)))
                           ENDIF
#else
                           CALL WRITE_GBL_FILE_SKIP_DEFAULT
     &                        (UnitNumber,TRIM(GLOBALDIR)//'/'//TRIM(FileName(IW)),
     &                           SwanDescript,DT*ITIME,ITIME,StoreTwo,DBLE(OVEXCV(IVTYPE)))
#endif
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO

      IF(IT.EQ.SWAN_MTC)THEN
         ALLOCATE(SWAN_MAXTEMP(1:NP))
!... Copy code from write_output.F to write the global maximum files.
         DO IW=1,NumSwanOutput
            IF((IW.NE.3).AND.(IW.NE.5))THEN
               FileNameMax = "swan_"//TRIM(Names(IW))//"_max.63"
               DO IP=1,NP
                  SWAN_MAXTEMP(IP) = SWAN_MAX(IP,IW)
               ENDDO
               CALL OPEN_GBL_FILE(963, TRIM(GLOBALDIR)//'/'//TRIM(FileNameMax),
     &              NP_G, NP, HEADER_MAX)
               SwanDescriptMax % initial_value        = 0.
               SwanDescriptMax % num_items_per_record = 1
               SwanDescriptMax % num_fd_records       = NP_G
               SwanDescriptMax % num_records_this     = NP
               SwanDescriptMax % imap                 => NODES_LG
               SwanDescriptMax % array                => SWAN_MAXTEMP
               CALL WRITE_GBL_FILE(963, TRIM(GLOBALDIR)//'/'//TRIM(FileNameMax),
     &              SwanDescriptMax, DT*ITIME, ITIME, StoreOne)
            ENDIF
         ENDDO
         DEALLOCATE(SWAN_MAXTEMP)
      ENDIF

      IF(ALLOCATED(Names))    DEALLOCATE(Names)
      IF(ALLOCATED(SWAN_VOQ)) DEALLOCATE(SWAN_VOQ)

      RETURN

      END SUBROUTINE SwanOutput



      SUBROUTINE PADCSWAN_INIT

      USE GLOBAL,   ONLY: DT,
     &                    ETA2,
     &                    NP,
     &                    WVNX1,
     &                    WVNX2,
     &                    WVNY1,
     &                    WVNY2
      USE TIMECOMM, ONLY: SWAN_DT => DT

      IMPLICIT NONE

      INTRINSIC :: ALLOCATED

      INTEGER   :: IN  ! Node counter.

Casey 090302: Allocate memory for the ADCIRC values that may be
C             passed to SWAN.
      IF(.NOT.ALLOCATED(SWAN_ETA2))THEN
         ALLOCATE(SWAN_ETA2(1:NP,1:2))
      ENDIF
      IF(.NOT.ALLOCATED(SWAN_UU2))THEN
         ALLOCATE(SWAN_UU2(1:NP,1:2))
      ENDIF
      IF(.NOT.ALLOCATED(SWAN_VV2))THEN
         ALLOCATE(SWAN_VV2(1:NP,1:2))
      ENDIF
      IF(.NOT.ALLOCATED(SWAN_WX2))THEN
         ALLOCATE(SWAN_WX2(1:NP,1:2))
      ENDIF
      IF(.NOT.ALLOCATED(SWAN_WY2))THEN
         ALLOCATE(SWAN_WY2(1:NP,1:2))
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

      END SUBROUTINE



      SUBROUTINE PADCSWAN_RUN(ITIME)

      USE GLOBAL,   ONLY: DP,
     &                    DT,
     &                    ETA2,
     &                    NODECODE,
     &                    NP,
     &                    UU2,
     &                    VV2
      USE TIMECOMM, ONLY: SWAN_DT => DT

      IMPLICIT NONE

      INTRINSIC    :: REAL

      INTEGER      :: IN           ! Node counter.
      INTEGER      :: IT           ! Time step counter.
      INTEGER      :: ITIME        ! Time step from ADCIRC.
      INTEGER,SAVE :: SwanTimeStep ! Counter for SWAN time steps.

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

      DO IT=1,(CouplingInterval*DT)/SWAN_DT

Casey 090302: If the coupling interval is not the same as the SWAN
C             time step, then the following logic would find a
C             correct interpolation weight.  As it is now,
C             we just use a weight of zero.
C        InterpoWeight = REAL(IT-1)*SWAN_DT/REAL(CouplingInterval*DT)
         InterpoWeight = 0.D0
         SwanTimeStep = SwanTimeStep + 1
         CALL SWMAIN(ITIME,SwanTimeStep)

      ENDDO

Casey 090302: Move the 'old' values into the first columns.
C             We do this for the winds after the SWAN time steps
C             because the 'new' values are written over in TIMESTEP.
C             We need to move them to the 'old' values now,
C             so that we do not lose them during the ADCIRC time stepping.
      DO IN=1,NP
         SWAN_WX2(IN,1) = SWAN_WX2(IN,2)
         SWAN_WY2(IN,1) = SWAN_WY2(IN,2)
      ENDDO

      END SUBROUTINE



      SUBROUTINE PADCSWAN_FINAL

      IMPLICIT NONE

      CALL SWEXITMPI

      DEALLOCATE(SWAN_ETA2)
      DEALLOCATE(SWAN_UU2)
      DEALLOCATE(SWAN_VV2)
      DEALLOCATE(SWAN_WX2)
      DEALLOCATE(SWAN_WY2)

      END SUBROUTINE



#endif
      END MODULE Couple2Swan

