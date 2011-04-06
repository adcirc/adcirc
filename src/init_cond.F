C******************************************************************************
C PADCIRC VERSION 48.33 04/09/2009                                            *
C  last changes in this file VERSION 48.33                                    *
C                                                                             *
C  mod history                                                                *
C  v48.33     - 04/09/09 - kmd- from 45.06 - developed a initial condition    *
C                                          file for the reading in of initial *
C                                          baroclinic values. Format for this *
C                                          file is developed from the hot     *
C                                          start subroutine file.             *
C                                                                             *
C******************************************************************************
C 
      SUBROUTINE INIT_COND()
C     
C**************************************************************************
C     
C     INITIAL CONDITION PROGRAM SETUP ROUTINE 
C     
C**************************************************************************
C     
      USE GLOBAL
      USE GLOBAL_IO
! kmd - added this in for heat flux boundary conditions
      USE GLOBAL_3DVS, ONLY : qsurf
! kmd - end of additions
      USE SIZES
      USE HARM
      USE WIND
      USE OWIWIND,ONLY : NWS12INIT,NWS12GET   ! sb46.28sb01 added 09/xx/2006
      USE OWI_ICE, ONLY : NCICE1_INIT, NCICE1_GET  !tcm v49.64.02
      USE RS2,ONLY : RS2INIT,RS2GET           ! sb46.28sb01 added 09/xx/2006
      USE NodalAttributes, ONLY : 
     &     LoadDirEffRLen,
     &     ApplyDirectionalWindReduction,
     &     OutputTau0 !jgf47.06
      IMPLICIT NONE

      INTEGER IT
      INTEGER I, J              !local loop counters
      INTEGER NCyc, NA, InputFileFmtVn
      INTEGER NTHS

      REAL(SZ) HollandTime
      REAL(SZ) ArgT, ArgTP, ArgSAlt
      REAL(SZ) H2     
      REAL(SZ) CCSFEA
      REAL(SZ) QTRatio
      REAL(SZ) SAltMul, S2SFEA
      REAL(SZ) TPMul
      REAL(SZ) WTRatio, WindX, WindY, WindMag, WDragCo
      REAL(SZ) WindDragLimit    !jgf46.01 Cap on max wind drag coefficient
      REAL(8) Time
C     jgf46.08 Fine grained ramp functions
      REAL(SZ) RampExtFlux1   ! Ramp for external flux b.c.s @ITHS-1
      REAL(SZ) RampExtFlux2   ! Ramp for external flux b.c.s @ITHS
      REAL(SZ) RampIntFlux1   ! Ramp for internal flux b.c.s @ITHS-1
      REAL(SZ) RampIntFlux2   ! Ramp for internal flux b.c.s @ITHS
      REAL(SZ) RampTip2    ! Ramp for tidal potential @ITHS
      REAL(SZ) RampMete2   ! Ramp for wind and atmospheric pressure @ITHS 
      REAL(SZ) RampWRad2   ! Ramp for wave radiation stress @ITHS
C     jgf46.14 Check for file existence before attempting to open.
      LOGICAL FileFound    ! .True. if the the file exists
      INTEGER ErrorIO      ! zero if file opened successfully
      INTEGER NP_G_IN, NE_G_IN  ! Global and
      INTEGER NP_A_IN, NE_A_IN  ! active number nodes, elements.

      CHARACTER(80) :: gridinfo

! kmd - added in for heat flux boundary condition

      INTEGER :: NOD
      REAL(SZ), ALLOCATABLE :: TMP(:,:)
      REAL(SZ) :: CD, CDQ, QWIND
      CHARACTER(80) :: CDUM80
      REAL(SZ), PARAMETER :: RHOAIR=1.225
      REAL(SZ) RNDAYHS
! kmd - end of additions

C     tcm v49.64.02 -- ADDED variables for Ice concentration fields
      REAL(SZ) PIC,CICE_TRatio  !ICE VARIABLES


C...  
C......READ IN THE INITIAL CONDITION OVER WHOLE DOMAIN
C......THIS FILE IS SET UP FOR INFORMATION FROM HYCOM FOR THE 
C......SET UP OF A 3D SIMULATION
C...  

C     Determine if the initial condition file exists.
      FileFound=.False.
      IF (IHOT.eq.17) THEN
         INQUIRE(FILE=TRIM(INPUTDIR)//'/'//'fort.17',EXIST=FileFound)
      ENDIF
      IF (.not.FileFound) THEN
         WRITE(16,1001)         ! hot start file  
         WRITE(16,1011)         ! was not found.
         WRITE(16,9973)         ! execution terminated
         IF (NScreen.ne.0.and.MyProc.eq.0) THEN
            WRITE(ScreenUnit,1001)       
            WRITE(ScreenUnit,1011)       
            WRITE(16,9973)      ! execution terminated
         ENDIF
         STOP
      ENDIF
C
C     Now open the initial condition (unit 17) file.
      WRITE(16,240) 
 240  FORMAT(/,9X,'Initial condition file was found.',
     &     ' Opening file.') 
      IF(IHOT.EQ.17) OPEN(17,FILE=TRIM(INPUTDIR)//'/'//'fort.17')
      IF (ErrorIO.GT.0) THEN 
         WRITE(16,1001)         ! hot start file
         WRITE(16,1005)         ! exists but can't be opened
         WRITE(16,9973)         ! execution terminated
         IF (NScreen.ne.0.and.MyProc.eq.0) THEN
            WRITE(ScreenUnit,1001) 
            WRITE(ScreenUnit,1005)
            WRITE(ScreenUnit,9973) 
         ENDIF
         STOP                   ! We're toast.
      ENDIF
C
 1001 FORMAT('ERROR: The initial conditions file')
 1011 FORMAT('was not found.') 
 1012 FORMAT('was a nonmatching version') 
 1005 FORMAT('exists but cannot be opened.')
 9973 FORMAT(/,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
C
C     Now read the initial condition file.
! This file will not include the version number information as we won't 
! know this until the program has run, so it has been replaced with the 
! grid information.  
      READ(IHOT,*) gridinfo 

      READ(IHOT,*) IMHS
      READ(IHOT,*) TIME
      READ(IHOT,*) ITHS
! kmd - added in these lines for 3D and baroclinic 
      IF ((IDEN.GE.1).AND.(C3D).AND.(CBAROCLINIC)) THEN
         CHOTHS=.TRUE.
      END IF
      IF (CHOTHS) THEN    
        DTDPHS=TIME/ITHS ! Set up time step information for hotstart file
        IF (DTDPHS.NE.DTDP) THEN 
          RNDAYHS=TIME/86400.D0
          NTHS=INT((RNDAY-RNDAYHS)*(86400.D0/DTDP)+0.5d0)
          NT=ITHS+NTHS
          WRITE(16,*) "The number of time steps has been reset",
     &          " due to a change in the timestep interval"
          WRITE(16,*) "It has been changed to ", NT
        END IF
       ELSE 
        DTDPHS=DTDP
       END IF 
! kmd - end of additions
      READ(IHOT,*) NP_G_IN
      READ(IHOT,*) NE_G_IN
      READ(IHOT,*) NP_A_IN
      READ(IHOT,*) NE_A_IN

      DO I=1, NP
         READ(IHOT,*) ETA1(I)
      END DO

      DO I=1, NP
         READ(IHOT,*) ETA2(I)
      END DO

      DO I=1, NP
         READ(IHOT,*) UU2(I)
      END DO

      DO I=1, NP
         READ(IHOT,*) VV2(I)
      END DO

      IF(IMHS.EQ.10) THEN
        DO I=1, NP
           READ(IHOT,*) CH1(I)
        END DO
      ENDIF

      DO I=1, NP
         READ(IHOT,*) NNODECODE(I) 
      END DO

      DO I=1,NE
         READ(IHOT,*) NOFF(I)    
      END DO

      DO I=1, NP
         ETAS(I)=ETA2(I)-ETA1(I)
         H2=DP(I)+IFNLFA*ETA2(I)
         QX2(I)=UU2(I)*H2
         QY2(I)=VV2(I)*H2
         EtaDisc(I)=ETA2(I)
      END DO

      IF (MYPROC.EQ.23) THEN 
         PRINT *, ETAS(15785)
         PRINT *, ETA1(15785)
         PRINT *, ETA2(15785)
         PRINT *, UU2(15785)
         PRINT *, VV2(15785)
      END IF 

C     jgf46.08 Fine grained ramp functions
C     jgf46.21 Split flux into internal and external, added support
C     for IBTYPE=52.
      IF (NRamp.eq.0) THEN
         RampExtFlux1=1.0d0
         RampExtFlux2=1.0d0
         RampIntFlux1=1.0d0
         RampIntFlux2=1.0d0
         RampTip2=1.0d0
         RampMete2=1.0d0
         RampWRad2=1.0d0
      ELSE
         RampExtFlux1=TANH((2.D0*(ITHS-1)*DTDPHS/86400.D0)/DRampExtFlux)
         RampExtFlux2=TANH((2.D0*(ITHS)*DTDPHS/86400.D0)/DRampExtFlux)
         RampIntFlux1=TANH((2.D0*(ITHS-1)*DTDPHS/86400.D0)/DRampIntFlux)
         RampIntFlux2=TANH((2.D0*(ITHS)*DTDPHS/86400.D0)/DRampIntFlux)
         RampTip2=TANH((2.D0*(ITHS)*DTDPHS/86400.D0)/DRampTip)
         RampMete2=TANH((2.D0*(ITHS)*DTDPHS/86400.D0)/DRampMete)
         RampWRad2=TANH((2.D0*(ITHS)*DTDPHS/86400.D0)/DRampWRad)
      ENDIF
#ifdef IBM
      FluxSettlingIT=INT(FluxSettlingTime*86400.d0/DTDPHS,KIND(0.0D0))
#else
      FluxSettlingIT=INT(FluxSettlingTime*86400.d0/DTDPHS)
#endif
      IF(ITHS.LT.(FluxSettlingIT+10)) THEN
         RampIntFlux1=0.0d0
         RampIntFlux2=0.0d0
         RampTip2=0.0d0
         RampMete2=0.0d0
         RampWRad2=0.0d0
      ENDIF
     
Cjgf46.02 comm.out for Katrina  IF(NRAMP.EQ.1) THEN
Cjgf46.02 comm.out for Katrina    RAMP1=TANH((2.D0*(ITHS-1)*DTDP/86400.D0)/DRAMP)
Cjgf46.02 comm.out for Katrina    RAMP2=TANH((2.D0*ITHS*DTDP/86400.D0)/DRAMP)
Cjgf46.02 comm.out for Katrina  ENDIF

C     
C.... SET POSITIONS IN BOUNDARY CONDITION, WIND AND OUTPUT FILES
C     
      WRITE(16,1112)
      WRITE(16,1794)
 1794 FORMAT(//,' INFORMATION ABOUT RE-STARTING THE TIME SERIES',
     &     ' OUTPUT FILES (UNITS 61-64,71-74,81,83),',
     &     /,' WIND/PRESSURE FILE (UNIT 22) AND FLOW BOUNDARY',
     &     ' CONDITION FILE (UNIT 20)',//)

C......INITIALLY, ZERO OUT THE NORMAL FLOW ON ALL BOUNDARIES

      DO I=1,NVEL
         QN2(I)=0.D0
         QN1(I)=0.D0
         QN0(I)=0.D0
      END DO

C.... FIND THE PROPER PLACE IN THE APERIODIC ELEVATION SPECIFIED
C.... BOUNDARY CONDITION FILE IF IT IS REQURIED.

      IF((NOPE.GT.0).AND.(NBFR.EQ.0)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
         WRITE(16,1112)
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1976)
         WRITE(16,1976)
 1976    FORMAT(/,1X,'LOCATING ELEVATION SPECIFIED INFORMATION IN ',
     &        'UNIT 19',/)
         OPEN(19,FILE=TRIM(INPUTDIR)//'/'//'fort.19')
         READ(19,*) ETIMINC
         ETIME1=STATIM*86400.D0
         ETIME2=ETIME1+ETIMINC
         DO J=1,NETA
            READ(19,*) ESBIN1(J)
         END DO
         DO J=1,NETA
            READ(19,*) ESBIN2(J)
         END DO
         DO IT=1,ITHS-1
            TIMEIT=IT*DTDPHS + STATIM*86400.D0
            IF(TIMEIT.GT.ETIME2) THEN
               ETIME1=ETIME2
               ETIME2=ETIME1+ETIMINC
               DO J=1,NETA
                  ESBIN1(J)=ESBIN2(J)
                  READ(19,*) ESBIN2(J)
               END DO
            ENDIF
         END DO
         IF(TIME.GT.ETIME2) THEN
            ETIME1=ETIME2
            ETIME2=ETIME1+ETIMINC
            DO J=1,NETA
               ESBIN1(J)=ESBIN2(J)
               READ(19,*) ESBIN2(J)
            END DO
         ENDIF
         ETRATIO=(TIMEIT-ETIME1)/ETIMINC
      ENDIF

C......FIND PROPER PLACE IN THE APERIODIC NORMAL FLOW BOUNDARY CONDITION
C......FILE IF IT IS REQUIRED

      IF((NFLUXF.EQ.1).AND.(NFFR.EQ.0)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
         WRITE(16,1112)
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1978)
         WRITE(16,1978)
 1978    FORMAT(/,1X,'LOCATING NORMAL FLOW INFORMATION IN UNIT 20',/)
         OPEN(20,FILE=TRIM(INPUTDIR)//'/'//'fort.20')
         READ(20,*) FTIMINC
         QTIME1=STATIM*86400.D0
         QTIME2=QTIME1+FTIMINC
         DO J=1,NVEL
            QNIN1(J)=0.D0
            IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &           .OR.(LBCODEI(J).EQ.22))
     &           READ(20,*) QNIN1(J)
         END DO
         DO J=1,NVEL
            QNIN2(J)=0.D0
            IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &           .OR.(LBCODEI(J).EQ.22))
     &           READ(20,*) QNIN2(J)
         END DO
         DO IT=1,ITHS-1
            TIMEIT=IT*DTDPHS + STATIM*86400.D0
            IF(TIMEIT.GT.QTIME2) THEN
               QTIME1=QTIME2
               QTIME2=QTIME2+FTIMINC
               DO J=1,NVEL
                  IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &                 .OR.(LBCODEI(J).EQ.22)) THEN
                     QNIN1(J)=QNIN2(J)
                     READ(20,*) QNIN2(J)
                  ENDIF
               END DO
            ENDIF
         END DO
         QTRATIO=(TIMEIT-QTIME1)/FTIMINC
         DO I=1,NVEL
            QN1(I)=RampExtFlux1*(QNIN1(I)+QTRATIO*(QNIN2(I)-QNIN1(I)))
         END DO
         IF(TIME.GT.QTIME2) THEN
            QTIME1=QTIME2
            QTIME2=QTIME1+FTIMINC
            DO J=1,NVEL
               IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &              .OR.(LBCODEI(J).EQ.22)) THEN
                  QNIN1(J)=QNIN2(J)
                  READ(20,*) QNIN2(J)
               ENDIF
            END DO
         ENDIF
         QTRATIO=(TIME-QTIME1)/FTIMINC
         DO I=1,NVEL
            QN2(I)=RampExtFlux2*(QNIN1(I)+QTRATIO*(QNIN2(I)-QNIN1(I)))
         END DO
      ENDIF

C......RESTART THE PERIODIC NORMAL FLOW BOUNDARY CONDITION

      IF((NFLUXF.EQ.1).AND.(NFFR.GT.0)) THEN
         DO J=1,NFFR
            IF(FPER(J).EQ.0.) THEN
               NCYC=0.
            ELSE
#ifdef IBM
               NCYC=INT(TIME/FPER(J),KIND(0.0d0))
#else
               NCYC=INT(TIME/FPER(J))
#endif
            ENDIF
            ARGJ1=FAMIG(J)*(TIME-DTDPHS-NCYC*FPER(J))+FFACE(J)
            ARGJ2=FAMIG(J)*(TIME-NCYC*FPER(J))+FFACE(J)
            RFF1=FFF(J)*RampExtFlux1
            RFF2=FFF(J)*RampExtFlux2
            DO I=1,NVEL
               ARG1=ARGJ1-QNPH(J,I)
               ARG2=ARGJ2-QNPH(J,I)
               QN1(I)=QN1(I)+QNAM(J,I)*RFF1*COS(ARG1)
               QN2(I)=QN2(I)+QNAM(J,I)*RFF2*COS(ARG2)
            END DO
         END DO
      ENDIF

C...  
C...  RESTART SUPERCRITICAL OUTWARD NORMAL FLOW OVER SPECIFIED
C.... EXTERNAL BARRIER BOUNDARY NODES
C...  
      IF(NFLUXB.EQ.1) THEN
         DO I=1,NVEL
            IF((LBCODEI(I).EQ.3).OR.(LBCODEI(I).EQ.13)
     &           .OR.(LBCODEI(I).EQ.23)) THEN
               NNBB=NBV(I)
               RBARWL=2.D0*(ETA1(NNBB)-BARLANHT(I))/3.D0
               IF(RBARWL.GT.0.0D0) THEN
                  QN1(I)=-RampIntFlux1
     &                 *BARLANCFSP(I)*RBARWL*(RBARWL*G)**0.5D0
               ENDIF
               RBARWL=2.D0*(ETA2(NNBB)-BARLANHT(I))/3.D0
               IF(RBARWL.GT.0.0D0) THEN
                  QN2(I)=-RampIntFlux2
     &                 *BARLANCFSP(I)*RBARWL*(RBARWL*G)**0.5D0
               ENDIF
            ENDIF
         END DO
      ENDIF

C...  
C...  RESTART INWARD/OUTWARD NORMAL FLOW OVER SPECIFIED
cjj   wm001 - modified/added the following 3 lines
C.... INTERNAL BARRIERS AND FOR INTERNAL BARRIER BOUNDARIES
C.... WITH CROSS BARRIER PIPES
C.... THIS SECTION ONLY RESTARTS THE OVER BARRIER FLOW COMPONENT
C...  
cjj   wm001 - modified following line
      IF(NFLUXIB.EQ.1) THEN
         DO I=1,NVEL
cjj   wm001 - modified following 2 lines
            IF((LBCODEI(I).EQ.4).OR.(LBCODEI(I).EQ.24) 
     &           .OR.(LBCODEI(I).EQ.5).OR.(LBCODEI(I).EQ.25)) THEN
               NNBB1=NBV(I)     ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
               NNBB2=IBCONN(I)  ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
C.............RESET INFORMATION FOR K-1 TIME LEVEL
               RBARWL1=ETA1(NNBB1)-BARINHT(I)
               RBARWL2=ETA1(NNBB2)-BARINHT(I)
               RBARWL1F=2.0D0*RBARWL1/3.0D0
               RBARWL2F=2.0D0*RBARWL2/3.0D0
               IF((RBARWL1.LT.0.0).AND.(RBARWL2.LT.0.0)) THEN ! WATER LEVEL BELOW BARRIER
                  QN1(I)=0.0D0  ! NO FLOW
                  GOTO 1998
               ENDIF
               IF(RBARWL1.EQ.RBARWL2) THEN ! WATER LEVEL EQUAL ON BOTH SIDES OF BARRIER
                  QN1(I)=0.0D0  ! NO FLOW
                  GOTO 1998
               ENDIF
               IF(RBARWL1.GT.RBARWL2) THEN ! WATER LEVEL GREATER ON THIS SIDE OF BARRIER
                  IF(RBARWL2.GT.RBARWL1F) THEN ! OUTWARD SUBCRITICAL FLOW
                     QN1(I)=-RampIntFlux1*BARINCFSB(I)*RBARWL2*
     &                    (2.d0*G*(RBARWL1-RBARWL2))**0.5D0
                     GOTO 1998
                  ELSE          ! OUTWARD SUPERCRITICAL FLOW
                     QN1(I)=-RampIntFlux1
     &                    *BARINCFSP(I)*RBARWL1F*(RBARWL1F*G)**0.5D0
                     GOTO 1998
                  ENDIF
               ENDIF
               IF(RBARWL2.GT.RBARWL1) THEN ! WATER LEVEL LOWER ON THIS SIDE OF BARRIER
                  IF(RBARWL1.GT.RBARWL2F) THEN ! INWARD SUBCRITICAL FLOW
                     QN1(I)=RampIntFlux1*BARINCFSB(I)*RBARWL1*
     &                    (2.d0*G*(RBARWL2-RBARWL1))**0.5D0
                     GOTO 1998
                  ELSE          ! INWARD SUPERCRITICAL FLOW
                     QN1(I)=RampIntFlux1
     &                    *BARINCFSP(I)*RBARWL2F*(RBARWL2F*G)**0.5D0
                     GOTO 1998
                  ENDIF
               ENDIF
 1998          CONTINUE
C.............RESET INFORMATION FOR K TIME LEVEL
               RBARWL1=ETA2(NNBB1)-BARINHT(I)
               RBARWL2=ETA2(NNBB2)-BARINHT(I)
               RBARWL1F=2.0D0*RBARWL1/3.0D0
               RBARWL2F=2.0D0*RBARWL2/3.0D0
               IF((RBARWL1.LT.0.0).AND.(RBARWL2.LT.0.0)) THEN ! WATER LEVEL BELOW BARRIER
                  QN2(I)=0.0D0  ! NO FLOW
                  GOTO 1999
               ENDIF
               IF(RBARWL1.EQ.RBARWL2) THEN ! WATER LEVEL EQUAL ON BOTH SIDES OF BARRIER
                  QN2(I)=0.0D0  ! NO FLOW
                  GOTO 1999
               ENDIF
               IF(RBARWL1.GT.RBARWL2) THEN ! WATER LEVEL GREATER ON THIS SIDE OF BARRIER
                  IF(RBARWL2.GT.RBARWL1F) THEN ! OUTWARD SUBCRITICAL FLOW
                     QN2(I)=-RampIntFlux2*BARINCFSB(I)*RBARWL2*
     &                    (2.d0*G*(RBARWL1-RBARWL2))**0.5D0
                     GOTO 1999
                  ELSE          ! OUTWARD SUPERCRITICAL FLOW
                     QN2(I)=-RampIntFlux2
     &                    *BARINCFSP(I)*RBARWL1F*(RBARWL1F*G)**0.5D0
                     GOTO 1999
                  ENDIF
               ENDIF
               IF(RBARWL2.GT.RBARWL1) THEN !WATER LEVEL LOWER ON THIS SIDE OF BARRIER
                  IF(RBARWL1.GT.RBARWL2F) THEN ! INWARD SUBCRITICAL FLOW
                     QN2(I)=RampIntFlux2*BARINCFSB(I)*RBARWL1*
     &                    (2.d0*G*(RBARWL2-RBARWL1))**0.5D0
                     GOTO 1999
                  ELSE          ! INWARD SUPERCRITICAL FLOW
                     QN2(I)=RampIntFlux2
     &                    *BARINCFSP(I)*RBARWL2F*(RBARWL2F*G)**0.5D0
                     GOTO 1999
                  ENDIF
               ENDIF
 1999          CONTINUE
            ENDIF
         END DO
      ENDIF

cjj   wm001 - start add
C...  
C...  RESTART INWARD/OUTWARD NORMAL FLOW OVER SPECIFIED
C.... INTERNAL BARRIERS WITH CROSS BARRIER PIPES
C.... THIS SECTION RESTARTS THE PIPE FLOW COMPONENT
C.... NOTE THAT PIPE FLOW COMPONENT IS ADDED INTO BARRIER FLOW COMPONENT
C.... THAT WAS PREVIOUSLY SET
C...  
      IF(NFLUXIBP.EQ.1) THEN
         DO I=1,NVEL
            IF((LBCODEI(I).EQ.5).OR.(LBCODEI(I).EQ.25)) THEN
               NNBB1=NBV(I)     ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
               NNBB2=IBCONN(I)  ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
C.............RESET INFORMATION FOR K-1 TIME LEVEL
               RBARWL1=ETA1(NNBB1)-PIPEHT(I)
               RBARWL2=ETA1(NNBB2)-PIPEHT(I)
               IF((RBARWL1.LT.0.0).AND.(RBARWL2.LT.0.0)) THEN ! WATER LEVEL BELOW PIPE
                  QN1(I)=QN1(I)+0.0D0 ! NO FLOW
                  GOTO 2002
               ENDIF
               IF(RBARWL1.EQ.RBARWL2) THEN ! WATER LEVEL EQUAL ON BOTH SIDES OF PIPE
                  QN1(I)=QN1(I)+0.0D0 ! NO FLOW
                  GOTO 2002
               ENDIF
               IF(RBARWL1.GT.RBARWL2) THEN ! WATER LEVEL GREATER ON THIS SIDE OF PIPE
                  IF(RBARWL2.LE.0) THEN ! OUTWARD FREE DISCHARGE 
                     QN1(I)=QN1(I)-RampIntFlux1
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*RBARWL1/(1+PIPECOEF(I)))**0.5D0
                     GOTO 2002
                  ELSE          ! OUTWARD SUBMERGED DISCHARGE
                     QN1(I)=QN1(I)-RampIntFlux1
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*(RBARWL1-RBARWL2)/PIPECOEF(I))**0.5D0
                     GOTO 2002
                  ENDIF
               ENDIF
               IF(RBARWL2.GT.RBARWL1) THEN ! WATER LEVEL LOWER ON THIS SIDE OF PIPE
                  IF(RBARWL1.LE.0) THEN ! INWARD FREE DISCHARGE 
                     QN1(I)=QN1(I)+RampIntFlux1
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*RBARWL2/(1+PIPECOEF(I)))**0.5D0
                     GOTO 2002
                  ELSE          ! INWARD SUBMERGED DISCHARGE
                     QN1(I)=QN1(I)+RampIntFlux1
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*(RBARWL2-RBARWL1)/PIPECOEF(I))**0.5D0
                     GOTO 2002
                  ENDIF
               ENDIF
 2002          CONTINUE
C.............RESET INFORMATION FOR K TIME LEVEL
               RBARWL1=ETA2(NNBB1)-PIPEHT(I)
               RBARWL2=ETA2(NNBB2)-PIPEHT(I)
               IF((RBARWL1.LT.0.0).AND.(RBARWL2.LT.0.0)) THEN ! WATER LEVEL BELOW PIPE
                  QN2(I)=QN2(I)+0.0D0 ! NO FLOW
                  GOTO 2003
               ENDIF
               IF(RBARWL1.EQ.RBARWL2) THEN ! WATER LEVEL EQUAL ON BOTH SIDES OF PIPE
                  QN2(I)=QN2(I)+0.0D0 ! NO FLOW
                  GOTO 2003
               ENDIF
               IF(RBARWL1.GT.RBARWL2) THEN ! WATER LEVEL GREATER ON THIS SIDE OF PIPE
                  IF(RBARWL2.LE.0) THEN ! OUTWARD FREE DISCHARGE
                     QN2(I)=QN2(I)-RampIntFlux2
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*RBARWL1/(1+PIPECOEF(I)))**0.5D0
                     GOTO 2003
                  ELSE          ! OUTWARD SUBMERGED DISCHARGE
                     QN2(I)=QN2(I)-RampIntFlux2
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*(RBARWL1-RBARWL2)/PIPECOEF(I))**0.5D0
                     GOTO 2003
                  ENDIF
               ENDIF
               IF(RBARWL2.GT.RBARWL1) THEN !WATER LEVEL LOWER ON THIS SIDE OF PIPE
                  IF(RBARWL1.LE.0) THEN ! INWARD FREE DISCHARGE
                     QN2(I)=QN2(I)+RampIntFlux2
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*RBARWL2/(1+PIPECOEF(I)))**0.5D0
                     GOTO 2003
                  ELSE          ! INWARD SUBMERGED DISCHARGE
                     QN2(I)=QN2(I)+RampIntFlux2
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*(RBARWL2-RBARWL1)/PIPECOEF(I))**0.5D0
                     GOTO 2003
                  ENDIF
               ENDIF
 2003          CONTINUE
            ENDIF
         END DO
      ENDIF
cjjwm001 - end add    

C------------------------ICE FIELDS----------------------------------------
C...  UPDATE THE ICE CONCENTRATION FIELDS FROM UNIT 25, 225, 227
C...  TCM V49.64.02 ADDED THE ICE FIELDS SECTION
C
      IF (NCICE.EQ.12) THEN 
          CICE_TIME1 = STATIM*86400.D0
          CICE_TIME2 = CICE_TIME1 + CICE_TIMINC
          ! This just initializes some variables and might set CICE,etc...
          ! if skipping ahead is used
          CALL NCICE1_INIT(CICE1,NP)
          CALL NCICE1_GET(CICE1,NP)
          CALL NCICE1_GET(CICE2,NP)
          DO IT=1,ITHS
             TIMEIT=IT*DTDP + STATIM*86400.D0
             IF(TIMEIT.GT.CICE_TIME2) THEN
                CICE_TIME1 = CICE_TIME2
                CICE_TIME2 = CICE_TIME2 + CICE_TIMINC
                DO I=1,NP
                   CICE1(I)=CICE2(I)
                END DO
                CALL NCICE1_GET(CICE2,NP)
             ENDIF
          END DO
      ENDIF
C
C

C------------------------MET FORCING---------------------------------------
C
C......RESTART WIND AND PRESSURE INFORMATION
C

      !kmd - added in the zeroing out of the wind and pressure arrays
      DO I=1,NP
         WSX1(I)=0.D0
         WSY1(I)=0.D0
         WVNXOUT(I)=0.D0
         WVNYOUT(I)=0.D0
         PR1(I)=0.D0
         WSX2(I)=0.D0
         WSY2(I)=0.D0
         PR2(I)=0.D0
         RSNX2(I)=0.D0
         RSNY2(I)=0.D0
      END DO 

      WindDragLimit = 0.0035
C     No wind, radiation stress or atmospheric pressure forcings are used.
      IF(NWS.EQ.1) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         DO J=1,ITHS
            DO I=1,NP
               READ(22,*) NHG,WSX2(I),WSY2(I),PR2(I)
            END DO
         END DO
         DO I=1,NP
            WSX2(I)=RampMete2*WSX2(I)
            WSY2(I)=RampMete2*WSY2(I)
            PR2(I)=RampMete2*PR2(I)
         END DO
      ENDIF
      
C     Wind stress and atmospheric pressure are read in at all grid nodes
C     at a time interval that does not equal the model time
C     step. Interpolation in time is used to synchronize the wind and
C     pressure information with the model time step.
      IF(NWS.EQ.2) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         DO IT=1,ITHS
            TIMEIT=IT*DTDPHS + STATIM*86400.D0
            IF(TIMEIT.GT.WTIME2) THEN
               WTIME1=WTIME2
               WTIME2=WTIME2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
                  PRN1(I)=PRN2(I)
                  READ(22,*) NHG,WVNX2(I),WVNY2(I),PRN2(I)
               END DO
            ENDIF
         END DO
         WTRATIO=(TIME-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WSX2(I) = RampMete2*WINDX
            WSY2(I) = RampMete2*WINDY
            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
         END DO
      ENDIF

      IF(NWS.EQ.-2) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = TIME
         WTIME2 = WTIME1 + WTIMINC
         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         DO I=1,NP
            WSX2(I) = RampMete2*WVNX1(I)
            WSY2(I) = RampMete2*WVNY1(I)
            PR2(I) = RampMete2*PRN1(I)
         END DO
      ENDIF

C     Wind velocity in US Navy Fleet Numeric format interpolated in
C     space onto the ADCIRC grid and in time to synchronize the wind and
C     pressure information with the model time step. Garratt's formula
C     is used to compute wind stress from the wind velocity.
      IF(NWS.EQ.3) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
 2223    CALL NWS3GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,IWTIME,IWYR,WTIMED,NP,
     &        NWLON,NWLAT,WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,
     &        NScreen,ScreenUnit)
         IF(IWYR.NE.IREFYR) THEN
            IWTIMEP=IWTIME
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
            END DO
            GOTO 2223
         ENDIF
         IF(WTIMED.LE.WREFTIM) THEN
            IWTIMEP=IWTIME
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
            END DO
            GOTO 2223
         ENDIF
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) 
     &        WRITE(ScreenUnit,*)'FOUND WIND DATA AT TIME= ',IWTIMEP
         WRITE(16,*) 'FOUND WIND DATA AT TIME =',IWTIMEP
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) 
     &        WRITE(ScreenUnit,*)'FOUND WIND DATA AT TIME= ',IWTIME
         WRITE(16,*) 'FOUND WIND DATA AT TIME =',IWTIME
         WTIME2=WTIMED-WREFTIM  !CAST INTO MODEL TIME REFERENCE
         WTIME1=WTIME2-WTIMINC
         DO IT=1,ITHS
            TIMEIT=IT*DTDPHS + STATIM*86400.D0
            IF(TIMEIT.GT.WTIME2) THEN
               WTIME1=WTIME2
               WTIME2=WTIME2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
               END DO
               CALL NWS3GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,IWTIME,
     &              IWYR,WTIMED,NP,NWLON,NWLAT,WLATMAX,WLONMIN,
     &              WLATINC,WLONINC,ICS,NScreen,ScreenUnit)
               IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) 
     &              WRITE(ScreenUnit,*) 
     &              'WIND FILE ADVANCED TO TIME',' = ', IWTIME
               WRITE(16,*) 'WIND FILE ADVANCED TO TIME = ',IWTIME
            ENDIF
         END DO
         WTRATIO=(TIME-WTIME1)/WTIMINC
         DO I=1,NP              !INTERPOLATE IN TIME
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG=SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
C    TCM V49.64.02 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TIME-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
            ENDIF
            WSX2(I)=RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I)=RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
         END DO
      ENDIF

C     Wind velocity and atmospheric pressure are read in (PBL/JAG
C     format) at selected ADCIRC grid nodes. Interpolation in time is
C     used to synchronize the wind and pressure information with the
C     model time step. Garratt's formula is used to compute wind stress
C     from wind velocity.
      IF(NWS.EQ.4) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
         CALL NWS4GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
         DO IT=1,ITHS
            TIMEIT=IT*DTDPHS + STATIM*86400.D0
            IF(TIMEIT.GT.WTIME2) THEN
               WTIME1=WTIME2
               WTIME2=WTIME2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
                  PRN1(I)=PRN2(I)
               END DO
               CALL NWS4GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
            ENDIF
         END DO
         WTRATIO=(TIME-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
C    TCM V49.64.02 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TIME-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I) = RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
         END DO
      ENDIF

      IF(NWS.EQ.-4) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = TIME
         WTIME2 = WTIME1 + WTIMINC
         CALL NWS4GET(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
         CALL NWS4GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
         DO I=1,NP
            WINDX = WVNX1(I)
            WINDY = WVNY1(I)
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
C    TCM V49.64.02 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TIME-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I) = RampMete2*PRN1(I)
         END DO
      ENDIF

C     Wind velocity and atmospheric pressure are read in at all grid
C     nodes. Interpolation in time is used to synchronize the wind and
C     pressure information with the model time step. Garratt's formula
C     is used to compute wind stress from wind velocity.
      IF(NWS.EQ.5) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         DO IT=1,ITHS
            TIMEIT=IT*DTDPHS + STATIM*86400.D0
            IF(TIMEIT.GT.WTIME2) THEN
               WTIME1=WTIME2
               WTIME2=WTIME2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
                  PRN1(I)=PRN2(I)
                  READ(22,*) NHG,WVNX2(I),WVNY2(I),PRN2(I)
               END DO
            ENDIF
         END DO
         WTRATIO=(TIME-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
C    TCM V49.64.02 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TIME-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I) = RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I))) 
C     jjw-42.06j wrote, jgf46.00 added following two lines:
            WVNXOUT(I)=RampMete2*WINDX
            WVNYOUT(I)=RampMete2*WINDY
         END DO
      ENDIF

      IF(NWS.EQ.-5) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = TIME
         WTIME2 = WTIME1 + WTIMINC
         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         WTRATIO=(TIME-WTIME1)/WTIMINC !jgf46.00 added
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I)) 
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I)) 
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
C    TCM V49.64.02 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TIME-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I) = RampMete2*PRN1(I)
            WVNXOUT(I)=RampMete*WINDX !jgf46.00 added 
            WVNYOUT(I)=RampMete*WINDY !jgf46.00 added
         END DO
      ENDIF

C     Wind velocity and atmospheric pressure are read in for a
C     rectangular grid (either in Longitude, Latitude or Cartesian
C     coordinates, consistent with the grid coordinates) and
C     interpolated in space onto the ADCIRC grid and in time to
C     synchronize the wind and pressure information with the model time
C     step. Garratt's formula is used to compute wind stress from the
C     wind velocity.
      IF(NWS.EQ.6) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
C     The following 3 lines are a hardwire to allow a non standard met
C     file to be read in at time zero in a hot start.  They should be
C     eliminated or commented out for normal operation
c     OPEN(199,FILE=TRIM(INPUTDIR)//'/'//'fort.199')
c     READ(199,*) (NHG,PRN1(I),WVNX1(I),WVNY1(I),I=1,NP)
c     CLOSE(199)
C     The following CALL statement should be uncommented for normal operation
         CALL NWS6GET(X,Y,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         CALL NWS6GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         WTIME1=TIME
         WTIME2=WTIME1+WTIMINC
         WTRATIO=(TIME-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
C    TCM V49.64.02 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TIME-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I) = RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
         END DO
      ENDIF

C     jgf46.01 New option to read in surface wind stress and atmospheric
C     pressure for a rectangular grid (either in Longitude, Latitude or
C     Cartesian coordinates, consistent with the grid coordinates) and
C     interpolate in space onto the ADCIRC grid. Interpolation in time
C     is used to synchronize the wind and pressure information with the
C     model time step.
C
      IF(NWS.EQ.7) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         CALL NWS7GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         CALL NWS7GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         WTIME1=TIME
         WTIME2=WTIME1+WTIMINC
         WTRATIO=(TIME-WTIME1)/WTIMINC
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WSX2(I) = RampMete2*WindX !apply ramp
            WSY2(I) = RampMete2*WindY
            PR2(I)=RampMete2*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            wvnxout(i)=WSX2(i) !for met recording sta. output
            wvnyout(i)=WSY2(i)             
         END DO
      ENDIF

      IF(NWS.EQ.-7) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         CALL NWS7GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         CALL NWS7GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         WTIME1=TIME
         WTIME2=WTIME1+WTIMINC
         WTRATIO=(TIME-WTIME1)/WTIMINC
C     jgfdebug46.01 How do we convert these marine wind stresses into
C     directional land surface stresses?
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WSX2(I) = RampMete2*WindX !apply ramp
            WSY2(I) = RampMete2*WindY
            PR2(I)=RampMete2*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            wvnxout(i)=WSX2(i) !for met recording sta. output
            wvnyout(i)=WSY2(i)             
         END DO
      ENDIF

C     jgf46.02 New option to read in hurricane locations and generate 
C     generate hurricane winds from the Holland Wind Model.
      IF(ABS(NWS).EQ.8) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         HollandTime = TIME
         CALL HollandGet(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,
     &        ICS,RHOWAT0,G,HollandTime,NSCREEN,ScreenUnit)
         DO I=1,NP
            WindX = WVNX2(I)
            WindY = WVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
C    TCM V49.64.02 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TIME-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete2*0.001293d0*WDragCo*WindY*WindMag
            PR2(I) = RampMete2*PRN2(I)
            WVNXOUT(I)=RampMete2*WindX
            WVNYOUT(I)=RampMete2*WindY
         ENDDO
      ENDIF 

!     cf & cm added nws = 9: asymmetric hurricane winds

      IF(ABS(NWS).EQ.9) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
             CALL NWS9GET(SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,TIME, ICS)
!         IF(NWS.EQ.9) THEN
!             CALL NWS9GET(SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,TIME, ICS)
!         ENDIF
!         IF(NWS.EQ.-9) THEN
!             CALL NWS9GET(SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,TIME, ICS)
!         ENDIF
         DO I=1,NP
            WindX = WVNX2(I)
            WindY = WVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
C    TCM V49.64.02 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TIME-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete2*0.001293d0*WDragCo*WindY*WindMag
            PR2(I) = RampMete2*PRN2(I)
            WVNXOUT(I)=RampMete2*WindX
            WVNYOUT(I)=RampMete2*WindY
         ENDDO
      ENDIF
C
C     Wind velocity (10 m) and atmospheric pressure are read in from a
C     sequence of National Weather Service (NWS) Aviation (AVN) model
C     output files. Each AVN file is assumed to contain data on a
C     Gaussian longitude, latitude grid at a single time. Consecutive
C     files in the sequence are separated by N hours in time. Garratt's
C     formula is used to compute wind stress from the wind velocity.
      IF(NWS.EQ.10) THEN
         WTIME1=TIME
         WTIME2=WTIME1+WTIMINC
         NWSGGWI=-1
         CALL NWS10GET(NWSGGWI,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,RHOWAT0,G,
     &        NWLON,NWLAT,WTIMINC)
         NWSGGWI=0
         CALL NWS10GET(NWSGGWI,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,RHOWAT0,G,
     &        NWLON,NWLAT,WTIMINC)
         NWSGGWI=1
         CALL NWS10GET(NWSGGWI,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,RHOWAT0,G,
     &        NWLON,NWLAT,WTIMINC)
         WTRATIO=(TIME-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
C    TCM V49.64.02 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TIME-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I) = RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
         END DO
      ENDIF

C     Wind velocity (10 m) and atmospheric pressure are read in from a
C     sequence of stripped down National Weather Service (NWS) ETA 29km
C     model output files. Each ETA file is assumed to contain data on an
C     E grid for a single day (8 data sets, one every 3 hours, beginning
C     @ 03:00 and continuing through 24:00 of the given day). The wind
C     data is converted to an east-west, north-south coordinate system
C     inside ADCIRC. Garratt's formula is used to compute wind stress
C     from the wind velocity.
      IF(NWS.EQ.11) THEN
         WTIME1=TIME
         WTIME2=WTIME1+WTIMINC
         NWSEGWI=0
         IDSETFLG=0
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1197)
         WRITE(16,1197)
         CALL NWS11GET(NWSEGWI,IDSETFLG,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,
     &        RHOWAT0,G)        !JUST COMPUTE INTERPOLATING FACTORS
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1198)
         WRITE(16,1198)
         NWSEGWI=0
         IDSETFLG=8
         CALL NWS11GET(NWSEGWI,IDSETFLG,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,
     &        RHOWAT0,G)        !NOW COMPUTE HOTSTART WIND FILED
         NWSEGWI=1
         IDSETFLG=1
         CALL NWS11GET(NWSEGWI,IDSETFLG,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,
     &        RHOWAT0,G)        !NOW COMPUTE NEXT WIND FIELD
         WTRATIO=(TIME-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
C    TCM V49.64.02 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TIME-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PR2(I) = RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
         END DO
      ENDIF

C.....sb46.28sb01 NWS=12,-12 were added to deal with raw OWI files.  09/xx/2006
      IF(NWS.EQ.12) THEN
          WTIME1 = STATIM*86400.D0
          WTIME2 = WTIME1 + WTIMINC
          CALL NWS12INIT(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
          CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
          DO IT=1,ITHS
             TIMEIT=IT*DTDPHS + STATIM*86400.D0
             IF(TIMEIT.GT.WTIME2) THEN
                WTIME1=WTIME2
                WTIME2=WTIME2+WTIMINC
                DO I=1,NP
                   WVNX1(I)=WVNX2(I)
                   WVNY1(I)=WVNY2(I)
                   PRN1(I)=PRN2(I)
                 END DO
                 CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
              ENDIF
          END DO
          WTRATIO=(TIME-WTIME1)/WTIMINC
          DO I=1,NP
             WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
             WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
             WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
             WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
             IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY) 
                WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
             ENDIF
C    TCM V49.64.02 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TIME-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
            ENDIF
             WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
             WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
             PR2(I) = RampMete2*PRN1(I)
             WVNXOUT(I)=RampMete2*WINDX
             WVNYOUT(I)=RampMete2*WINDY
          END DO
       ENDIF

       IF(NWS.EQ.-12) THEN
          WTIME1 = TIME
          WTIME2 = WTIME1 + WTIMINC
          CALL NWS12INIT(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
          CALL NWS12GET(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
          CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
          DO I=1,NP
             WINDX = WVNX1(I)
             WINDY = WVNY1(I)
             WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
             WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
             IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY) 
                WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
             ENDIF
C    TCM V49.64.02 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               CICE_TRatio = (TIME-CICE_TIME1)/CICE_TIMINC
               PIC = CICE1(I) + CICE_TRatio*(CICE2(I)-CICE1(I))
               WDragCo = WindIceDrag(WDragCo,PIC,WindDragLimit,"IceCube   ")
            ENDIF
             WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
             WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
             PR2(I) = RampMete2*PRN1(I)
             WVNXOUT(I)=RampMete*WINDX
             WVNYOUT(I)=RampMete*WINDY
          END DO
       ENDIF

!kmd - sponge layer for winds (add in sponge)
        IF ((C3D).AND.(RES_BC_FLAG.GT.0).AND.(NWS.NE.0)) THEN 
         DO I=1,NP
           WSX2(I)=sponge(I)*WSX2(I)
           WSY2(I)=sponge(I)*WSY2(I)
         END DO 
        END IF

!kmd - end of additions

C
C--------------------END MET FORCING---------------------------------------

C......RESTART THE WAVE RADIATION STRESS

        IF(NRS.GE.1) THEN ! sb46.28sb03
          RSTIME1 = TIME
          RSTIME2 = RSTIME1 + RSTIMINC
          ! NRS=2 was added. sb46.28sb03 09/xx/2006
          IF(NRS.EQ.1) THEN
            OPEN(23,FILE=TRIM(INPUTDIR)//'/'//'fort.23')
            CALL RSGET(RSNX1,RSNY1,NP)
            CALL RSGET(RSNX2,RSNY2,NP)
          ENDIF
          IF(NRS.EQ.2) THEN
            CALL RS2INIT(RSNX1,RSNY1,NP)
            CALL RS2GET(RSNX1,RSNY1,NP)
            CALL RS2GET(RSNX2,RSNY2,NP)
          ENDIF
          DO I=1,NP
            WSX2(I) = WSX2(I)+RampWRad2*RSNX1(I)
            WSY2(I) = WSY2(I)+RampWRad2*RSNY1(I)
            END DO
          ENDIF

       if (CTIP) then
Cjromo 11-01-00  Initialize TIP2 for HOTSTART
          DO I=1,NP
             TIP2(I)=0.0
          END DO
CTIP  LINES TO USE TIDAL POTENTIAL FORCING
       IF(NTIP.GE.1) THEN
         DO J=1,NTIF
           IF(PERT(J).EQ.0.) THEN
             NCYC=0
             ELSE
#ifdef IBM                
             NCYC=INT(TIME/PERT(J),KIND(0.0d0))
#else
             NCYC=INT(TIME/PERT(J))
#endif
             ENDIF
           ARGT=AMIGT(J)*(TIME-NCYC*PERT(J))+FACET(J)
           TPMUL=RampTip2*ETRF(J)*TPK(J)*FFT(J)
           SALTMUL=RampTip2*FFT(J)
#ifdef IBM
           NA=NINT(0.00014/AMIGT(J),KIND(0.0d0))
#else
           NA=NINT(0.00014/AMIGT(J))
#endif
           IF(NA.EQ.1) THEN                        !SEMI-DIURNAL SPECIES
             DO I=1,NP
               ARGTP=ARGT+2.*SLAM(I)
               ARGSALT=ARGT-SALTPHA(J,I)
               CCSFEA=COS(SFEA(I))
               CCSFEA=CCSFEA*CCSFEA
               TIP2(I)=TIP2(I)+TPMUL*CCSFEA*COS(ARGTP)
     &                 +SALTMUL*SALTAMP(J,I)*COS(ARGSALT)
               END DO
             ENDIF
           IF(NA.EQ.2) THEN
             DO I=1,NP
               ARGTP=ARGT+SLAM(I)
               ARGSALT=ARGT-SALTPHA(J,I)
cjjw/vjpm002 - modified/added the following 5 lines
#ifdef REAL8
               S2SFEA=SIN(2.d0*SFEA(I))
#else
               S2SFEA=SIN(2.e0*SFEA(I))
#endif
               TIP2(I)=TIP2(I)+TPMUL*S2SFEA*COS(ARGTP)
     &                +SALTMUL*SALTAMP(J,I)*COS(ARGSALT)
               END DO
             ENDIF
           END DO
         ENDIF
      endif     !   CTIP


C...
C....SET UP TO RESTART TIMESERIES OUTPUT FILES
C....
C...
      IF(NBYTE.EQ.4) ITEMPSTP=20
      IF(NBYTE.EQ.8) ITEMPSTP=10

C...
C....IF RESTARTING THE ELEVATION STATION OUTPUT FILE, GO TO THE PROPER PLACE
C....IN THE FILE.  OTHERWISE ZERO OUT NSCOUE.
C...
      IESTP=0
      NSCOUE=0
      WRITE(16,1040) IESTP,NSCOUE
 1040 FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN ELEVATION ',
     &                'STATION FILE BY THE TIME OF THE HOT START',
     &           /,8X,'SPOOL COUNTER = ',I6)
      IF(NOUTE.LT.0) THEN
        IESTP=0
        NSCOUE=0
        IF((NTCYSE.LT.ITHS).AND.(NSPOOLE.GT.0)) THEN
          NTCYSE=NTCYSE+((ITHS-NTCYSE)/NSPOOLE)*NSPOOLE
          IF(NTCYSE.LT.ITHS) NTCYSE=NTCYSE+NSPOOLE
          IF(NSPOOLE.NE.0) NTRSPE=(NTCYFE-NTCYSE)/NSPOOLE
        ENDIF
        WRITE(16,1041)
 1041   FORMAT(//,' A NEW ELEVATION STATION FILE WILL BE STARTED')
      ENDIF
      IF(NOUTE.EQ.-2) THEN
        OPEN(61,FILE=TRIM(LOCALDIR)//'/'//'fort.61',
     &                      ACCESS='DIRECT',RECL=NBYTE)
        IF(NBYTE.EQ.4) THEN
          DO I=1,8
            WRITE(61,REC=IESTP+I) RDES4(I)
          ENDDO
          IESTP=IESTP+8
          DO I=1,6
            WRITE(61,REC=IESTP+I) RID4(I)
          ENDDO
          IESTP=IESTP+6
          DO I=1,6
            WRITE(61,REC=IESTP+I) AID4(I)
          ENDDO
          IESTP=IESTP+6
        ENDIF
        IF(NBYTE.EQ.8) THEN
          DO I=1,4
            WRITE(61,REC=IESTP+I) RDES8(I)
          ENDDO
          IESTP=IESTP+4
          DO I=1,3
            WRITE(61,REC=IESTP+I) RID8(I)
          ENDDO
          IESTP=IESTP+3
          DO I=1,3
            WRITE(61,REC=IESTP+I) AID8(I)
          ENDDO
          IESTP=IESTP+3
        ENDIF
        WRITE(61,REC=IESTP+1) NTRSPE
        WRITE(61,REC=IESTP+2) NSTAE
        WRITE(61,REC=IESTP+3) DT*NSPOOLE
        WRITE(61,REC=IESTP+4) NSPOOLE
        WRITE(61,REC=IESTP+5) 1
        IESTP=IESTP+5
        CLOSE(61)                    ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF
      IF(NOUTE.EQ.-1) THEN
        CALL OPEN_GBL_FILE(61, TRIM(GLOBALDIR)//'/'//'fort.61',
     $    NSTAE_G, NSTAE, HEADER61) 
        IESTP=2
      ENDIF

      IF(NOUTE.EQ.2) THEN
        OPEN(61,FILE=TRIM(LOCALDIR)//'/'//'fort.61',
     &                      ACCESS='DIRECT',RECL=NBYTE)
c       WRITE(61,REC=ITEMPSTP+1) NTRSPE   ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
        CLOSE(61)                         ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF

C...
C....GO TO THE PROPER PLACE IN THE VELOCITY STATION OUTPUT FILE
C...
      IVSTP=0
      NSCOUV=0
      WRITE(16,1042) IVSTP,NSCOUV
 1042 FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN VELOCITY ',
     &                'STATION FILE BY THE TIME OF THE HOT START',
     &           /,8X,'SPOOL COUNTER =',I6)
      IF(NOUTV.LT.0) THEN
        IVSTP=0
        NSCOUV=0
        IF((NTCYSV.LT.ITHS).AND.(NSPOOLV.GT.0)) THEN
          NTCYSV=NTCYSV+((ITHS-NTCYSV)/NSPOOLV)*NSPOOLV
          IF(NTCYSV.LT.ITHS) NTCYSV=NTCYSV+NSPOOLV
          NTRSPV=(NTCYFV-NTCYSV)/NSPOOLV
        ENDIF
        WRITE(16,1043)
 1043   FORMAT(//,' A NEW VELOCITY STATION FILE WILL BE STARTED')
      ENDIF

      IF(NOUTV.EQ.-2) THEN
        OPEN(62,FILE=TRIM(LOCALDIR)//'/'//'fort.62',
     &                      ACCESS='DIRECT',RECL=NBYTE)
        IF(NBYTE.EQ.4) THEN
          DO I=1,8
            WRITE(62,REC=IVSTP+I) RDES4(I)
          ENDDO
          IVSTP=IVSTP+8
          DO I=1,6
            WRITE(62,REC=IVSTP+I) RID4(I)
          ENDDO
          IVSTP=IVSTP+6
          DO I=1,6
            WRITE(62,REC=IVSTP+I) AID4(I)
          ENDDO
          IVSTP=IVSTP+6
          ENDIF
        IF(NBYTE.EQ.8) THEN
          DO I=1,4
            WRITE(62,REC=IVSTP+I) RDES8(I)
          ENDDO
          IVSTP=IVSTP+4
          DO I=1,3
            WRITE(62,REC=IVSTP+I) RID8(I)
          ENDDO
          IVSTP=IVSTP+3
          DO I=1,3
            WRITE(62,REC=IVSTP+I) AID8(I)
          ENDDO
          IVSTP=IVSTP+3
        ENDIF
        WRITE(62,REC=IVSTP+1) NTRSPV
        WRITE(62,REC=IVSTP+2) NSTAV
        WRITE(62,REC=IVSTP+3) DT*NSPOOLV
        WRITE(62,REC=IVSTP+4) NSPOOLV
        WRITE(62,REC=IVSTP+5) 2
        IVSTP=IVSTP+5
        CLOSE(62)                    ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF
      IF(NOUTV.EQ.-1) THEN
        CALL OPEN_GBL_FILE(62, TRIM(GLOBALDIR)//'/'//'fort.62',
     $    NSTAV_G, NSTAV, HEADER62)
        IVSTP=2
      ENDIF
      IF(NOUTV.EQ.2) THEN
        OPEN(62,FILE=TRIM(LOCALDIR)//'/'//'fort.62',
     &                      ACCESS='DIRECT',RECL=NBYTE)
c       WRITE(62,REC=ITEMPSTP+1) NTRSPV   ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
        CLOSE(62)                         ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF

C...
C....GO TO THE PROPER PLACE IN THE CONCENTRATION STATION OUTPUT FILE
C...
      ICSTP=0
      NSCOUC=0
      WRITE(16,1044) ICSTP,NSCOUC
 1044 FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN CONCENTRATION ',
     &                'STATION FILE BY THE TIME OF THE HOT START',
     &           /,8X,'SPOOL COUNTER = ',I6)
      IF(NOUTC.LT.0) THEN
        ICSTP=0
        NSCOUC=0
        IF((NTCYSC.LT.ITHS).AND.(NSPOOLC.GT.0)) THEN
          NTCYSC=NTCYSC+((ITHS-NTCYSC)/NSPOOLC)*NSPOOLC
          IF(NTCYSC.LT.ITHS) NTCYSC=NTCYSC+NSPOOLC
          NTRSPC=(NTCYFC-NTCYSC)/NSPOOLC
        ENDIF
        WRITE(16,1045)
 1045   FORMAT(//,' A NEW CONCENTRATION STATION FILE WILL BE STARTED')
      ENDIF

      IF(NOUTC.EQ.-2) THEN
        OPEN(81,FILE=TRIM(LOCALDIR)//'/'//'fort.81',
     &                      ACCESS='DIRECT',RECL=NBYTE)
        IF(NBYTE.EQ.4) THEN
          DO I=1,8
            WRITE(81,REC=ICSTP+I) RDES4(I)
          ENDDO
          ICSTP=ICSTP+8
          DO I=1,6
            WRITE(81,REC=ICSTP+I) RID4(I)
          ENDDO
          ICSTP=ICSTP+6
          DO I=1,6
            WRITE(81,REC=ICSTP+I) AID4(I)
          ENDDO
          ICSTP=ICSTP+6
        ENDIF
        IF(NBYTE.EQ.8) THEN
          DO I=1,4
            WRITE(81,REC=ICSTP+I) RDES8(I)
          ENDDO
          ICSTP=ICSTP+4
          DO I=1,3
            WRITE(81,REC=ICSTP+I) RID8(I)
          ENDDO
          ICSTP=ICSTP+3
          DO I=1,3
            WRITE(81,REC=ICSTP+I) AID8(I)
          ENDDO
          ICSTP=ICSTP+3
        ENDIF
        WRITE(81,REC=ICSTP+1) NTRSPC
        WRITE(81,REC=ICSTP+2) NSTAC
        WRITE(81,REC=ICSTP+3) DT*NSPOOLC
        WRITE(81,REC=ICSTP+4) NSPOOLC
        WRITE(81,REC=ICSTP+5) 1
        ICSTP=ICSTP+5
        CLOSE(81)                    ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF
      IF(NOUTC.EQ.-1) THEN
        CALL OPEN_GBL_FILE(81, TRIM(GLOBALDIR)//'/'//'fort.81',
     $    NSTAC_G, NSTAC, HEADER81) 
        ICSTP=2
        CLOSE(81)
      ENDIF
      IF(NOUTC.EQ.2) THEN
        OPEN(81,FILE=TRIM(LOCALDIR)//'/'//'fort.81',
     &                      ACCESS='DIRECT',RECL=NBYTE)
c       WRITE(81,REC=ITEMPSTP+1) NTRSPC   ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
        CLOSE(81)                         ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF

C...
C....GO TO THE PROPER PLACE IN THE METEOROLOGICAL STATION OUTPUT FILE
C...
      IPSTP=0
      IWSTP=0
      NSCOUM=0
      IICESTP = 0
c... v49.64.01 tcm -added support for ice stations
      IF(NCICE.EQ.0)  THEN
         WRITE(16,1038) IWSTP,IPSTP,NSCOUM
      ELSE
         IICESTP = IPSTP  !SET TO THE SAME PLACE AS THE PRESSURE FILE
         WRITE(16,2038) IWSTP,IPSTP,IICESTP,NSCOUM
      ENDIF
      WRITE(16,1038) IWSTP,IPSTP,NSCOUM
 1038 FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE WIND STATION',
     &                ' FILE BY THE TIME OF THE HOT START',
     &        /,1X,I6,' LINES OR RECORDS WRITTIN IN THE PRES STATION',
     &                ' FILE BY THE TIME OF THE HOT START',
     &           /,8X,'SPOOL COUNTER = ',I6)
 2038 FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE WIND STATION',
     &                ' FILE BY THE TIME OF THE HOT START',
     &        /,1X,I6,' LINES OR RECORDS WRITTEN IN THE PRES STATION',
     &                ' FILE BY THE TIME OF THE HOT START',
     &        /,1X,I6,' LINES OR RECORDS WRITTEN IN THE ICE STATION',
     &                ' FILE BY THE TIME OF THE HOT START',
     &           /,8X,'SPOOL COUNTER = ',I6)
C
      IF(NOUTM.LT.0) THEN
        IPSTP=0
        IWSTP=0
        IICESTP = 0  !v49.64.01 tcm added ice station
        NSCOUM=0
        IF((NTCYSM.LT.ITHS).AND.(NSPOOLM.GT.0)) THEN
          NTCYSM=NTCYSM+((REAL(ITHS)-REAL(NTCYSM))/
     &           REAL(NSPOOLM))*REAL(NSPOOLM)
          IF(NTCYSM.LT.ITHS) NTCYSM=NTCYSM+NSPOOLM
          NTRSPM=(NTCYFM-NTCYSM)/NSPOOLM
        ENDIF
        WRITE(16,1039)
 1039   FORMAT(//,' A NEW METEOROLOGICAL STATION FILE WILL BE STARTED')
      ENDIF

      IF(NOUTM.EQ.-2) THEN
        OPEN(71,FILE=TRIM(LOCALDIR)//'/'//'fort.71',
     &                      ACCESS='DIRECT',RECL=NBYTE)
        OPEN(72,FILE=TRIM(LOCALDIR)//'/'//'fort.72',
     &                      ACCESS='DIRECT',RECL=NBYTE)
!       v49.64.01 tcm added ice station
        IF(NCICE.NE.0) THEN
           OPEN(91,FILE=TRIM(LOCALDIR)//'/'//'fort.91',
     &                      ACCESS='DIRECT',RECL=NBYTE)
        ENDIF
        IF(NBYTE.EQ.4) THEN
          DO I=1,8
            WRITE(71,REC=IPSTP+I) RDES4(I)
            WRITE(72,REC=IWSTP+I) RDES4(I)
            IF (NCICE.NE.0) WRITE(91,REC=IICESTP+I) RDES4(I)
          ENDDO
          IPSTP=IPSTP+8
          IWSTP=IWSTP+8
          IF (NCICE.NE.0) IICESTP = IICESTP+8
          DO I=1,6
            WRITE(71,REC=IPSTP+I) RID4(I)
            WRITE(72,REC=IWSTP+I) RID4(I)
            IF (NCICE.NE.0) WRITE(91,REC=IICESTP+I) RID4(I)
          ENDDO
          IPSTP=IPSTP+6
          IWSTP=IWSTP+6
          IF (NCICE.NE.0) IICESTP = IICESTP+6
          DO I=1,6
            WRITE(71,REC=IPSTP+I) AID4(I)
            WRITE(72,REC=IWSTP+I) AID4(I)
            IF (NCICE.NE.0) WRITE(91,REC=IICESTP+I) AID4(I)
          ENDDO
          IPSTP=IPSTP+6
          IWSTP=IWSTP+6
          IF (NCICE.NE.0) IICESTP = IICESTP+6
        ENDIF
        IF(NBYTE.EQ.8) THEN
          DO I=1,4
            WRITE(71,REC=IPSTP+I) RDES8(I)
            WRITE(72,REC=IWSTP+I) RDES8(I)
            IF (NCICE.NE.0) WRITE(91,REC=IICESTP+I) RDES8(I)
          ENDDO
          IPSTP=IPSTP+4
          IWSTP=IWSTP+4
          IF (NCICE.NE.0) IICESTP = IICESTP+4
          DO I=1,3
            WRITE(71,REC=IPSTP+I) RID8(I)
            WRITE(72,REC=IWSTP+I) RID8(I)
            IF (NCICE.NE.0) WRITE(91,REC=IICESTP+I) RID8(I)
          ENDDO
          IPSTP=IPSTP+3
          IWSTP=IWSTP+3
          IF (NCICE.NE.0) IICESTP = IICESTP+3
          DO I=1,3
            WRITE(71,REC=IPSTP+I) AID8(I)
            WRITE(72,REC=IWSTP+I) AID8(I)
            IF (NCICE.NE.0) WRITE(91,REC=IICESTP+I) AID8(I)
          ENDDO
          IPSTP=IPSTP+3
          IWSTP=IWSTP+3
          IF (NCICE.NE.0) IICESTP = IICESTP+3
        ENDIF
        WRITE(71,REC=IPSTP+1) NTRSPM
        WRITE(71,REC=IPSTP+2) NSTAM
        WRITE(71,REC=IPSTP+3) DT*NSPOOLM
        WRITE(71,REC=IPSTP+4) NSPOOLM
        WRITE(71,REC=IPSTP+5) 1
        WRITE(72,REC=IWSTP+1) NTRSPM
        WRITE(72,REC=IWSTP+2) NSTAM
        WRITE(72,REC=IWSTP+3) DT*NSPOOLM
        WRITE(72,REC=IWSTP+4) NSPOOLM
        WRITE(72,REC=IWSTP+5) 2
        IPSTP=IPSTP+5
        IWSTP=IWSTP+5
        CLOSE(71)                    ! DO THIS TO FLUSH THE WRITE BUFFER
        CLOSE(72)
        ! v49.64.01 tcm -- added ice stations
        IF(NCICE.NE.0) THEN
           WRITE(91,REC=IICESTP+1) NTRSPM
           WRITE(91,REC=IICESTP+2) NSTAM
           WRITE(91,REC=IICESTP+3) DT*NSPOOLM
           WRITE(91,REC=IICESTP+4) NSPOOLM
           WRITE(91,REC=IICESTP+5) 1
           IICESTP = IICESTP + 5
           CLOSE(91)
        ENDIF
      ENDIF

      IF(NOUTM.EQ.-1) THEN
        CALL OPEN_GBL_FILE(71, TRIM(GLOBALDIR)//'/'//'fort.71',
     $    NSTAM_G, NSTAM, HEADER71)
        IPSTP=2
        CALL OPEN_GBL_FILE(72, TRIM(GLOBALDIR)//'/'//'fort.72',
     $    NSTAM_G, NSTAM, HEADER72)
        IWSTP=2
!       v49.64.01 tcm added ice station support
        IF(NCICE.NE.0) THEN
           CALL OPEN_GBL_FILE(91, TRIM(GLOBALDIR)//'/'//'fort.91',
     $       NSTAM_G, NSTAM, HEADER91)
           IICESTP=2
        ENDIF
      ENDIF
      IF(NOUTM.EQ.2) THEN
        OPEN(71,FILE=TRIM(LOCALDIR)//'/'//'fort.71',
     &                      ACCESS='DIRECT',RECL=NBYTE)
        OPEN(72,FILE=TRIM(LOCALDIR)//'/'//'fort.72',
     &                      ACCESS='DIRECT',RECL=NBYTE)
c         WRITE(71,REC=ITEMPSTP+1) NTRSPM   ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
c         WRITE(72,REC=ITMEPSTP+1) NTRSPM
        CLOSE(71)                         ! DO THIS TO FLUSH THE WRITE BUFFER
        CLOSE(72)
!       v49.64.01 tcm added ice station support
        IF(NCICE.NE.0) THEN
           OPEN(91,FILE=TRIM(LOCALDIR)//'/'//'fort.91',
     &                      ACCESS='DIRECT',RECL=NBYTE)
           CLOSE(91)
        ENDIF
      ENDIF

C...
C....GO TO THE PROPER PLACE IN THE GLOBAL ELEVATION OUTPUT FILE
C...
        IGEP=0
        NSCOUGE=0
        WRITE(16,1046) IGEP,NSCOUGE
 1046   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'ELEVATION FILE BY THE TIME OF THE HOT START',
     &          /,8X,'SPOOL COUNTER =',I6)
        IF(NOUTGE.LT.0) THEN
          IGEP=0
          NSCOUGE=0
          IF((NTCYSGE.LT.ITHS).AND.(NSPOOLGE.GT.0)) THEN
            NTCYSGE=NTCYSGE+((ITHS-NTCYSGE)/NSPOOLGE)*NSPOOLGE
            IF(NTCYSGE.LT.ITHS) NTCYSGE=NTCYSGE+NSPOOLGE
            NDSETSE=(NTCYFGE-NTCYSGE)/NSPOOLGE
            ENDIF
          WRITE(16,1047)
 1047     FORMAT(//,' A NEW GLOBAL ELEVATION FILE WILL BE STARTED')
          ENDIF

        IF(NOUTGE.EQ.-2) THEN
          OPEN(63,FILE=TRIM(LOCALDIR)//'/'//'fort.63',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(63,REC=IGEP+I) RDES4(I)
              ENDDO
            IGEP=IGEP+8
            DO I=1,6
              WRITE(63,REC=IGEP+I) RID4(I)
              ENDDO
            IGEP=IGEP+6
            DO I=1,6
              WRITE(63,REC=IGEP+I) AID4(I)
              ENDDO
            IGEP=IGEP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(63,REC=IGEP+I) RDES8(I)
              ENDDO
            IGEP=IGEP+4
            DO I=1,3
              WRITE(63,REC=IGEP+I) RID8(I)
              ENDDO
            IGEP=IGEP+3
            DO I=1,3
              WRITE(63,REC=IGEP+I) AID8(I)
              ENDDO
            IGEP=IGEP+3
            ENDIF
          WRITE(63,REC=IGEP+1) NDSETSE
          WRITE(63,REC=IGEP+2) NP
          WRITE(63,REC=IGEP+3) DT*NSPOOLGE
          WRITE(63,REC=IGEP+4) NSPOOLGE
          WRITE(63,REC=IGEP+5) 1
          IGEP=IGEP+5
          CLOSE(63)                    ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF
        IF((NOUTGE.EQ.-1).OR.(NOUTGE.EQ.-4)) THEN
           CALL OPEN_GBL_FILE(63, TRIM(GLOBALDIR)//'/'//'fort.63',
     $         NP_G, NP, HEADER63)
C   kmd48.33 add the ability to output the sponge layer
          IF (OUTPUTSPONGE) THEN 
           CALL writeDomainHeader(92, 
     &          TRIM(GLOBALDIR)//'/'//'fort.92',
     $           NP_G, NP, 'sponge     ')
          END IF 
C
C        jgf47.06 Tau0 output is produced on the same schedule as
C        global elevation
           IF (OUTPUTTAU0) THEN
               CALL writeDomainHeader(10,
     &              TRIM(GLOBALDIR)//'/'//'fort.tau0',
     $               NP_G, NP, 'Tau0      ')
           ENDIF
           IGEP=2
        ENDIF
        IF(NOUTGE.EQ.2) THEN
          OPEN(63,FILE=TRIM(LOCALDIR)//'/'//'fort.63',
     &          ACCESS='DIRECT',RECL=NBYTE)
c         WRITE(63,REC=ITEMPSTP+1) NDSETSE  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          CLOSE(63)                         ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF

C...
C....GO TO THE PROPER PLACE IN THE GLOBAL VELOCITY OUTPUT FILE
C...
        IGVP=0
        NSCOUGV=0
        WRITE(16,1048) IGVP,NSCOUGV
 1048   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'VELOCITY FILE BY THE TIME OF THE HOT START',
     &          /,8X,'SPOOL COUNTER =',I6)
        IF(NOUTGV.LT.0) THEN
          IGVP=0
          NSCOUGV=0
          IF((NTCYSGV.LT.ITHS).AND.(NSPOOLGV.GT.0)) THEN
            NTCYSGV=NTCYSGV+((ITHS-NTCYSGV)/NSPOOLGV)*NSPOOLGV
            IF(NTCYSGV.LT.ITHS) NTCYSGV=NTCYSGV+NSPOOLGV
            NDSETSV=(NTCYFGV-NTCYSGV)/NSPOOLGV
            ENDIF
          WRITE(16,1049)
 1049     FORMAT(//,' A NEW GLOBAL VELOCITY FILE WILL BE STARTED')
          ENDIF

        IF(NOUTGV.EQ.-2) THEN
          OPEN(64,FILE=TRIM(LOCALDIR)//'/'//'fort.64',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(64,REC=IGVP+I) RDES4(I)
              ENDDO
            IGVP=IGVP+8
            DO I=1,6
              WRITE(64,REC=IGVP+I) RID4(I)
              ENDDO
            IGVP=IGVP+6
            DO I=1,6
              WRITE(64,REC=IGVP+I) AID4(I)
              ENDDO
            IGVP=IGVP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(64,REC=IGVP+I) RDES8(I)
              ENDDO
            IGVP=IGVP+4
            DO I=1,3
              WRITE(64,REC=IGVP+I) RID8(I)
              ENDDO
            IGVP=IGVP+3
            DO I=1,3
              WRITE(64,REC=IGVP+I) AID8(I)
              ENDDO
            IGVP=IGVP+3
            ENDIF
          WRITE(64,REC=IGVP+1) NDSETSV
          WRITE(64,REC=IGVP+2) NP
          WRITE(64,REC=IGVP+3) DT*NSPOOLGV
          WRITE(64,REC=IGVP+4) NSPOOLGV
          WRITE(64,REC=IGVP+5) 2
          IGVP=IGVP+5
          CLOSE(64)                    ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF
        IF((NOUTGV.EQ.-1).OR.(NOUTGV.EQ.-4)) THEN
          CALL OPEN_GBL_FILE(64, TRIM(GLOBALDIR)//'/'//'fort.64',
     $      NP_G, NP, HEADER64)
          IGVP=2
        ENDIF
        IF(NOUTGV.EQ.2) THEN
          OPEN(64,FILE=TRIM(LOCALDIR)//'/'//'fort.64',
     &          ACCESS='DIRECT',RECL=NBYTE)
c         WRITE(64,REC=ITEMPSTP+1) NDSETSV  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          CLOSE(64)                         ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF

C...
C....GO TO THE PROPER PLACE IN THE GLOBAL CONCENTRATION OUTPUT FILE
C...
        IGCP=0
        NSCOUGC=0
        WRITE(16,1053) IGCP,NSCOUGC
 1053   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                'CONCENTRATION FILE BY THE TIME OF THE HOT START',
     &          /,8X,'SPOOL COUNTER =',I6)
        IF(NOUTGC.LT.0) THEN
          IGCP=0
          NSCOUGC=0
          IF((NTCYSGC.LT.ITHS).AND.(NSPOOLGC.GT.0)) THEN
            NTCYSGC=NTCYSGC+((ITHS-NTCYSGC)/NSPOOLGC)*NSPOOLGC
            IF(NTCYSGC.LT.ITHS) NTCYSGC=NTCYSGC+NSPOOLGC
            NDSETSC=(NTCYFGC-NTCYSGC)/NSPOOLGC
            ENDIF
          WRITE(16,1054)
 1054     FORMAT(//,' A NEW GLOBAL CONCENTRATION FILE WILL BE STARTED')
        ENDIF

        IF(NOUTGC.EQ.-2) THEN
          OPEN(83,FILE=TRIM(LOCALDIR)//'/'//'fort.83',
     &           ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(83,REC=IGCP+I) RDES4(I)
              ENDDO
            IGCP=IGCP+8
            DO I=1,6
              WRITE(83,REC=IGCP+I) RID4(I)
              ENDDO
            IGCP=IGCP+6
            DO I=1,6
              WRITE(83,REC=IGCP+I) AID4(I)
              ENDDO
            IGCP=IGCP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(83,REC=IGCP+I) RDES8(I)
              ENDDO
            IGCP=IGCP+4
            DO I=1,3
              WRITE(83,REC=IGCP+I) RID8(I)
              ENDDO
            IGCP=IGCP+3
            DO I=1,3
              WRITE(83,REC=IGCP+I) AID8(I)
              ENDDO
            IGCP=IGCP+3
            ENDIF
          WRITE(83,REC=IGCP+1) NDSETSC
          WRITE(83,REC=IGCP+2) NP
          WRITE(83,REC=IGCP+3) DT*NSPOOLGC
          WRITE(83,REC=IGCP+4) NSPOOLGC
          WRITE(83,REC=IGCP+5) 1
          IGCP=IGCP+5
          CLOSE(83)                    ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF
        IF(NOUTGC.EQ.-1) THEN
          OPEN(83,FILE=TRIM(LOCALDIR)//'/'//'fort.83')
          WRITE(83,3220) RUNDES,RUNID,AGRID
          WRITE(83,3645) NDSETSC,NP,DTDP*NSPOOLGC,NSPOOLGC,1
          IGCP=2
          CLOSE(83)
        ENDIF
        IF(NOUTGC.EQ.1) THEN
          OPEN(83,FILE=TRIM(LOCALDIR)//'/'//'fort.83')
          DO I=1,IGCP           !I DON'T KNOW OF A PRACTICAL WAY TO CHANGE NDSETSC
            READ(83,1050)
 1050       FORMAT(1X)
          ENDDO
          ENDFILE(83)
          CLOSE(83)
        ENDIF
        IF(NOUTGC.EQ.2) THEN
          OPEN(83,FILE=TRIM(LOCALDIR)//'/'//'fort.83',
     &          ACCESS='DIRECT',RECL=NBYTE)
c         WRITE(83,REC=ITEMPSTP+1) NDSETSC  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          CLOSE(83)                         ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF

C...
C....GO TO THE PROPER PLACE IN THE GLOBAL METEOROLOGICAL OUTPUT FILES
C...
        IGPP=0
        IGWP=0
        NSCOUGW=0
        IGIP=0  !tcm v49.64.01 added for ice
!       V49.64.01 TCM ADDED GLOBAL ICE CONCENTRATION FILES
        IF(NCICE.NE.0) THEN
           WRITE(16,2055) IGWP,IGPP,IGIP,NSCOUGW
        ELSE
           WRITE(16,1055) IGWP,IGPP,NSCOUGW
        ENDIF
 1055   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'WIND FILE BY THE TIME OF THE HOT START',
     &          /,1X,I6,'LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'PRESSURE FILE BY THE TIME OF THE HOT START',
     &          /,8X,'SPOOL COUNTER =',I6)
 2055   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'WIND FILE BY THE TIME OF THE HOT START',
     &          /,1X,I6,'LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'PRESSURE FILE BY THE TIME OF THE HOT START',   
     &          /,1X,I6,'LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'ICE FILE BY THE TIME OF THE HOT START',
     &          /,8X,'SPOOL COUNTER =',I6)
c
        IF(NOUTGW.LT.0) THEN
          IGPP=0
          IGWP=0
          IGIP=0 !tcm v49.64.01 added for ice
          NSCOUGW=0
          IF((NTCYSGW.LT.ITHS).AND.(NSPOOLGW.GT.0)) THEN
            NTCYSGW=NTCYSGW+((ITHS-NTCYSGW)/NSPOOLGW)*NSPOOLGW
            IF(NTCYSGW.LT.ITHS) NTCYSGW=NTCYSGW+NSPOOLGW
            NDSETSW=(NTCYFGW-NTCYSGW)/NSPOOLGW
          ENDIF
! TCM V49.64.01 -- ADDED ICE SUPPORT
          IF(NCICE.NE.0) THEN
             WRITE(16,2056)
          ELSE
             WRITE(16,1056)
          ENDIF
          WRITE(16,1056)
 1056     FORMAT(//,' NEW GLOBAL WIND & pressure FILEs WILL BE STARTED')
 2056   FORMAT(//,' NEW GLOBAL WIND & pressure & Ice FILEs WILL BE STARTED')
        ENDIF

        IF(NOUTGW.EQ.-2) THEN
          OPEN(73,file=trim(LOCALDIR)//'/'//'fort.73',
     &          access='DIRECT',recl=nbyte)
          OPEN(74,FILE=TRIM(LOCALDIR)//'/'//'fort.74',
     &          ACCESS='DIRECT',RECL=NBYTE)
! TCM V49.64.01 -- ADDED ICE SUPPORT
          IF(NCICE.NE.0) THEN
             OPEN(93,FILE=TRIM(LOCALDIR)//'/'//'fort.93',
     &          ACCESS='DIRECT',RECL=NBYTE)
          ENDIF
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              write(73,rec=igpp+i) rdes4(i)
              WRITE(74,REC=IGWP+I) RDES4(I)
              IF(NCICE.NE.0) write(93,rec=igip+i) rdes4(i)
              ENDDO
            igpp=igpp+8
            IGWP=IGWP+8
            if(ncice.ne.0) igip = igip+8
            DO I=1,6
              write(73,rec=igpp+i) rid4(i)
              WRITE(74,REC=IGWP+I) RID4(I)
              IF(NCICE.NE.0) write(93,rec=igip+i) rid4(i)
            ENDDO
            igpp=igpp+6
            IGWP=IGWP+6
            if(ncice.ne.0) igip = igip+6
            DO I=1,6
              write(73,rec=igpp+i) aid4(i)
              WRITE(74,REC=IGWP+I) AID4(I)
              IF(NCICE.NE.0) write(93,rec=igip+i) aid4(i)
            ENDDO
            igpp=igpp+6
            IGWP=IGWP+6
            if(ncice.ne.0) igip = igip+6
          ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              write(73,rec=igpp+i) rdes8(i)
              WRITE(74,REC=IGWP+I) RDES8(I)
              IF(NCICE.NE.0) write(93,rec=igip+i) rdes8(i)
            ENDDO
            igpp=igpp+4
            IGWP=IGWP+4
            if(ncice.ne.0) igip = igip+4
            DO I=1,3
              write(73,rec=igpp+i) rid8(i)
              WRITE(74,REC=IGWP+I) RID8(I)
              IF(NCICE.NE.0) write(93,rec=igip+i) rid8(i)
            ENDDO
            igpp=igpp+3
            IGWP=IGWP+3
            if(ncice.ne.0) igip = igip+3
            DO I=1,3
              write(73,rec=igpp+i) aid8(i)
              WRITE(74,REC=IGWP+I) AID8(I)
              IF(NCICE.NE.0) write(93,rec=igip+i) aid8(i)
            ENDDO
            igpp=igpp+3
            IGWP=IGWP+3
            if(ncice.ne.0) igip = igip+3
          ENDIF
          WRITE(73,rec=igpp+1) ndsetsw
          WRITE(73,rec=igpp+2) np
          WRITE(73,rec=igpp+3) dt*nspoolgw
          WRITE(73,rec=igpp+4) nspoolgw
          WRITE(73,rec=igpp+5) 2
          IGPP=IGPP+5
          CLOSE(73)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          WRITE(74,REC=IGWP+1) NDSETSW
          WRITE(74,REC=IGWP+2) NP
          WRITE(74,REC=IGWP+3) DT*NSPOOLGW
          WRITE(74,REC=IGWP+4) NSPOOLGW
          WRITE(74,REC=IGWP+5) 2
          IGWP=IGWP+5
          CLOSE(74)                    ! DO THIS TO FLUSH THE WRITE BUFFER
! TCM V49.64.01 -- ADDED ICE SUPPORT
          IF(NCICE.NE.0) THEN
             WRITE(93,REC=IGIP+1) NDSETSW
             WRITE(93,REC=IGIP+2) NP
             WRITE(93,REC=IGIP+3) DT*NSPOOLGW
             WRITE(93,REC=IGIP+4) NSPOOLGW
             WRITE(93,REC=IGIP+5) 2
             IGIP=IGIP+5
             CLOSE(93)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          ENDIF
        ENDIF
        IF((NOUTGW.EQ.-1).OR.(NOUTGW.EQ.-4)) THEN
          CALL OPEN_GBL_FILE(73, TRIM(GLOBALDIR)//'/'//'fort.73',
     $      NP_G, NP, HEADER73)
          IGPP=2
          CALL OPEN_GBL_FILE(74, TRIM(GLOBALDIR)//'/'//'fort.74',
     $      NP_G, NP, HEADER74)
          IGWP=2
! TCM V49.64.01 -- ADDED ICE SUPPORT
          IF (NCICE.NE.0) THEN
            CALL OPEN_GBL_FILE(93, TRIM(GLOBALDIR)//'/'//'fort.93',
     $        NP_G, NP, HEADER93)
            IGIP=2
          ENDIF
        ENDIF
        IF(NOUTGW.EQ.2) THEN
          open(73,file=TRIM(LOCALDIR)//'/'//'fort.73',
     &          access='DIRECT',recl=nbyte)
c         WRITE(73,REC=ITempStp+1) ndsetsw  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          close(73)                         ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(74,FILE=TRIM(LOCALDIR)//'/'//'fort.74',
     &          ACCESS='DIRECT',RECL=NBYTE)
c         WRITE(74,REC=ITEMPSTP+1) NDSETSW  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          CLOSE(74)                         ! DO THIS TO FLUSH THE WRITE BUFFER
! TCM V49.64.01 -- ADDED ICE SUPPORT
          IF(NCICE.NE.0) THEN
             OPEN(93,FILE=TRIM(LOCALDIR)//'/'//'fort.93',
     &          ACCESS='DIRECT',RECL=NBYTE)
c         WRITE(93,REC=ITempStp+1) ndsetsw  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
             CLOSE(93)                      ! DO THIS TO FLUSH THE WRITE BUFFER
          ENDIF
        ENDIF


C...
C...  HOT START INFORMATION FOR A 3D RUN
C...

!kmd - added time and the time step of the hotstart file
      IF(C3D) THEN
        CALL INIT_COND_3D(TIME,ITHS)
      ENDIF


C...
C......HOT START INFORMATION FOR HARMONIC ANALYSIS
C...
C...
C........IF HARMONIC ANALYSIS HAS ALREADY BEGUN, READ IN HOT START
C........HARMONIC ANALYSIS, MEAN AND SQUARE INFO
C...
      IF (IHARIND.EQ.1) THEN
         IF (ITHS.GT.ITHAS) THEN
            IHOTSTP=IHOTSTP+1
            READ(IHOT,REC=IHOTSTP) ICHA
         ENDIF
         IF(ITHS.GE.IHABEG) THEN
            CALL readBinaryHAHotstart(IHOT, IHOTSTP)
            CALL checkHAHotstart()
         ENDIF
      ENDIF  
      CLOSE(IHOT)
C
 1112 FORMAT(/,1X,79('_'))
 1197 FORMAT(/,1X,'THE E29 MET GRID INTERPOLATING FACTORS ARE ',
     &                'BEING COMPUTED ')
 1198     FORMAT(1X,'FINISHED COMPUTING E29 INTERPOLATING FACTORS',/)
 3220 FORMAT(1X,A32,2X,A24,2X,A24)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C
      RETURN
      END

C******************************************************************************
C   Subroutine to initialize the 3D routines for a hot start including        *
C   reading the 3D portion of the HOTSTART file and finding proper places in  *
C   3D output files.                                                          *
C                                                                             *
C                                                                             *
c******************************************************************************

      SUBROUTINE INIT_COND_3D(TIME,ITHS)

      USE GLOBAL_3DVS
!kmd - added these for the boundary conditions
      USE GLOBAL, ONLY : RES_BC_FLAG, RBCTIME1, RBCTIME2,
     &                    RBCTIMEINC, SBCTIME1, SBCTIME2,
     &                    SBCTIMEINC, TBCTIME1, TBCTIME2,
     &                    TBCTIMEINC, BCSTATIM, RBCRATIO,
     &                    LNM_BC1, LNM_BC2, TIMEIT, STATIM,
     &                    SBCSTATIM, TBCSTATIM, SBCRATIO,
     &                    TBCRATIO, DTDPHS, LNM_BC, UU2, 
     &                    VV2, NBD, BCFLAG_LNM, BCFLAG_TEMP,
     &                    TTBCTIME1, TTBCTIME2, TTBCTIMEINC,
     &                    TTBCRATIO, q_heat1, q_heat2, HFLUX,
     &                    TTBCSTATIM
 
!kmd - end of additions


      IMPLICIT NONE

      INTEGER :: IRType                  !number of fields in output file
      INTEGER :: NH, N                   !horizontal & vertical loop counters
      REAL(SZ) :: ImagPartOfQ
!kmd - added these 
      INTEGER :: NumofBCNodes, J, NOD, IT
      INTEGER, INTENT(IN) :: ITHS
      REAL(SZ), INTENT(IN) :: TIME
      INTEGER :: IDenHS, IINDX
      CHARACTER(80) :: CDUM80
      REAL(SZ), ALLOCATABLE :: RealQ_Tmp(:,:)
      REAL(SZ), ALLOCATABLE :: TMP(:,:)
      INTEGER :: NumofNodes, K
!kmd - end of additions

C...
C...  Define format statements used to initialize 3D output files
C...
 499  FORMAT(1X,A32,2X,A24,2X,A24)
 498  FORMAT(1X,I10,1X,I10,1X,E15.7,I10,1X,I10,1X,I3)
 497  FORMAT(5X,'UNIT ',I2,' FORMAT WILL BE ASCII')
 496  FORMAT(5X,'UNIT ',I2,' FORMAT WILL BE BINARY')

      ALLOCATE(RealQ_Tmp(NP,NFEN))
C...
C...  Read in 3D portion of hotstart file
C...
      READ(IHOT,*) IDenHS
      N3DSD=0
      I3DSDRec=0
      N3DSV=0
      I3DSVRec=0
      N3DST=0
      I3DSTRec=0
      N3DGD=0
      I3DGDRec=0
      N3DGV=0
      I3DGVRec=0
      N3DGT=0
      I3DGTRec=0


!kmd - need to determine if IDEN parameter in Hotstart is different than IDEN from
!      fort.15 file. If so, determine which one to use in the simulation.
      IF (IDenHS.NE.IDEN) THEN
         IF ((CBAROCLINIC).AND.(RES_BC_FLAG.GE.1)) THEN
            ! IDEN value from fort.15 should be used
         ELSE IF ((CBAROCLINIC).AND.(RES_BC_FLAG.LE.0)) THEN
            IDEN=IDenHS
         END IF
      ENDIF
! kmd - end of additions
! kmd - changed the format for the global IO 
      DO NH=1,NP
         DUU(NH)=0.d0
         DUV(NH)=0.d0
         DVV(NH)=0.d0
         UU(NH)=UU2(NH)
         VV(NH)=VV2(NH)
      END DO 

      DO NH=1,NP
         READ(IHOT,*) BSX(NH)
      END DO 

      DO NH=1,NP
         READ(IHOT,*) BSY(NH)
      END DO 

      DO N=1, NFEN 
        DO NH=1,NP
          READ(IHOT,*) RealQ_Tmp(NH,N)
        END DO 
      END DO 

      DO N=1, NFEN 
        DO NH=1,NP
          READ(IHOT,*) ImagPartofQ
          Q(NH,N) = RealQ_Tmp(NH,N) + iy*ImagPartofQ 
        END DO 
      END DO 

      DO N=1, NFEN 
        DO NH=1,NP
          READ(IHOT,*) WZ(NH,N) 
        END DO 
      END DO 
      DO N=1, NFEN 
        DO NH=1,NP
          READ(IHOT,*) q20(NH,N)  
        END DO 
      END DO 
      DO N=1, NFEN 
        DO NH=1,NP
          READ(IHOT,*) l(NH,N) 
        END DO 
      END DO 

      IF (ABS(IDEN).EQ.1) THEN 
        DO N=1, NFEN 
          DO NH=1,NP
            READ(IHOT,*) SIGT(NH,N) 
          END DO 
        END DO 
      END IF 
      IF (ABS(IDEN).EQ.2) THEN 
        DO N=1, NFEN 
          DO NH=1,NP
            READ(IHOT,*) SAL(NH,N) 
          END DO 
        END DO 
      END IF 
      IF (ABS(IDEN).EQ.3) THEN 
        DO N=1, NFEN 
          DO NH=1,NP
            READ(IHOT,*) TEMP(NH,N) 
          END DO 
        END DO 
      END IF 
      IF (ABS(IDEN).EQ.4) THEN 
        DO N=1, NFEN 
          DO NH=1,NP
            READ(IHOT,*) SAL(NH,N) 
          END DO 
        END DO 
        DO N=1, NFEN 
          DO NH=1,NP
            READ(IHOT,*) TEMP(NH,N) 
          END DO 
        END DO 
      END IF 

      DEALLOCATE(RealQ_Tmp)

! kmd - Need to restart the boundary conditions for the baroclinic simulations
      IF((RES_BC_FLAG.GT.0).AND.(CBAROCLINIC)) THEN
         IF (ABS(RES_BC_FLAG).GE.1) THEN 
            IF (BCFLAG_LNM.EQ.1) THEN
              OPEN(35,FILE=TRIM(INPUTDIR)//'/'//'fort.35')
              RBCTIME1=BCSTATIM*86400.d0
              RBCTIME2=RBCTIME1+RBCTIMEINC
              READ(35,'(A)') CDUM80
              DO NumofBCNodes=1,NETA
                 READ(35,*) NOD,LNM_BC1(NumofBCNodes)
              END DO
              READ(35,'(A)') CDUM80
              DO NumofBCNodes=1,NETA
                 READ(35,*) NOD,LNM_BC2(NumofBCNodes)
              END DO
              DO IT=1,ITHS
                 TIMEIT=IT*DTDPHS + STATIM*86400.D0
                 IF(TIMEIT.GT.RBCTIME2) THEN
                    RBCTIME1=RBCTIME2
                    RBCTIME2=RBCTIME1+RBCTIMEINC
                    READ(35,'(A)') CDUM80
                    DO NumofBCNodes=1,NETA
                       LNM_BC1(NumofBCNodes)=LNM_BC2(NumofBCNodes)
                       READ(35,*) NOD,LNM_BC2(NumofBCNodes)
                    END DO
                 END IF
              END DO
              IF(TIME.GT.RBCTIME2) THEN
                 RBCTIME1=RBCTIME2
                 RBCTIME2=RBCTIME2+RBCTIMEINC
                 READ(35,'(A)') CDUM80
                 DO NumofBCNodes=1,NETA
                   LNM_BC1(NumofBCNodes)=LNM_BC2(NumofBCNodes)
                   READ(35,*) NOD,LNM_BC2(NumofBCNodes)
                 END DO
              END IF
              RBCRATIO=(TIMEIT-RBCTIME1)/RBCTIMEINC
            ELSE IF (BCFLAG_LNM.EQ.2) THEN
              ! to be added later but currently information 
              ! is taken in two ways from the HYCOM results
            ELSE IF (BCFLAG_LNM.EQ.3) THEN 
              DO NumofBCNodes=1,NETA
                 IINDX=NBD(NumofBCNodes)
                 LNM_BC1(NumofBCNodes)=ETA2(IINDX)
                 LNM_BC2(NumofBCNodes)=LNM_BC1(NumofBCNodes)
              END DO 
            END IF 
         END IF 
         IF ((ABS(RES_BC_FLAG).EQ.2).OR.(ABS(RES_BC_FLAG).EQ.4)) THEN
            OPEN(36,FILE=TRIM(INPUTDIR)//'/'//'fort.36')
            SBCTIME1=SBCSTATIM*86400.D0
            SBCTIME2=SBCTIME1+SBCTIMEINC
            READ(36,'(A)') CDUM80
            DO NumofBCNodes=1,NETA
               READ(36,*) NOD,(RESSAL1(NumofBCNodes,J),J=1,NFEN)
            END DO
            READ(36,'(A)') CDUM80
            DO NumofBCNodes=1,NETA
               READ(36,*) NOD,(RESSAL2(NumofBCNodes,J),J=1,NFEN)
            END DO
            DO IT=1,ITHS
               TIMEIT=IT*DTDPHS + STATIM*86400.D0
               IF(TIMEIT.GT.SBCTIME2) THEN
                  SBCTIME1=SBCTIME2
                  SBCTIME2=SBCTIME1+SBCTIMEINC
                  READ(36,'(A)') CDUM80
                  DO NumofBCNodes=1,NETA
                     DO J=1,NFEN
                       RESSAL1(NumofBCNodes,J)=RESSAL2(NumofBCNodes,J)
                     END DO
                     READ(36,*) NOD,(RESSAL2(NumofBCNodes,J),J=1,NFEN)
                  END DO
               END IF
            END DO
            IF(TIME.GT.SBCTIME2) THEN
               SBCTIME1=SBCTIME2
               SBCTIME2=SBCTIME2+SBCTIMEINC
               READ(36,'(A)') CDUM80
               DO NumofBCNodes=1,NETA
                 DO J=1,NFEN
                   RESSAL1(NumofBCNodes,J)=RESSAL2(NumofBCNodes,J)
                 END DO
                 READ(36,*) NOD,(RESSAL2(NumofBCNodes,J),J=1,NFEN)
               END DO
            END IF
            SBCRATIO=(TIMEIT-SBCTIME1)/SBCTIMEINC
         END IF
         IF ((ABS(RES_BC_FLAG).EQ.3).OR.
     &            (ABS(RES_BC_FLAG).EQ.4)) THEN
            OPEN(37,FILE=TRIM(INPUTDIR)//'/'//'fort.37')
            TBCTIME1=TBCSTATIM*86400.D0
            TBCTIME2=TBCTIME1+TBCTIMEINC
            READ(37,'(A)') CDUM80
            DO NumofBCNodes=1,NETA
               READ(37,*) NOD,(RESTEMP1(NumofBCNodes,J),J=1,NFEN)
            END DO
            READ(37,'(A)') CDUM80
            DO NumofBCNodes=1,NETA
               READ(37,*) NOD,(RESTEMP2(NumofBCNodes,J),J=1,NFEN)
            END DO
            DO IT=1,ITHS
               TIMEIT=IT*DTDPHS + STATIM*86400.D0
               IF(TIMEIT.GT.TBCTIME2) THEN
                  TBCTIME1=TBCTIME2
                  TBCTIME2=TBCTIME1+TBCTIMEINC
                  READ(37,'(A)') CDUM80
                  DO NumofBCNodes=1,NETA
                     DO J=1,NFEN
                       RESTEMP1(NumofBCNodes,J)=RESTEMP2(NumofBCNodes,J)
                     END DO
                     READ(37,*) NOD,(RESTEMP2(NumofBCNodes,J),J=1,NFEN)
                  END DO
               END IF
            END DO
            IF(TIME.GT.TBCTIME2) THEN
               TBCTIME1=TBCTIME2
               TBCTIME2=TBCTIME2+TBCTIMEINC
               READ(37,'(A)') CDUM80
               DO NumofBCNodes=1,NETA
                 DO J=1,NFEN
                   RESTEMP1(NumofBCNodes,J)=RESTEMP2(NumofBCNodes,J)
                 END DO
                 READ(37,*) NOD,(RESTEMP2(NumofBCNodes,J),J=1,NFEN)
               END DO
            END IF
            TBCRATIO=(TIMEIT-TBCTIME1)/TBCTIMEINC
            IF (BCFLAG_TEMP.EQ.1) THEN  ! read in one value for surface heat flux
               OPEN(38,FILE=TRIM(INPUTDIR)//'/'//'fort.38')
               TTBCTIME1=TTBCSTATIM*86400.D0
               TTBCTIME2=TTBCTIME1+TTBCTIMEINC
               READ(38,'(A)') CDUM80
               DO NumofNodes=1,NH
                  READ(38,*) NOD, q_heat1(Numofnodes)
               END DO
               READ(38,'(A)') CDUM80
               DO NumofNodes=1,NH
                  READ(38,*) NOD, q_heat2(NumofNodes)
               END DO
               DO IT=1,ITHS
                  TIMEIT=IT*DTDPHS + STATIM*86400.D0
                  IF(TIMEIT.GT.TTBCTIME2) THEN
                     TTBCTIME1=TTBCTIME2
                     TTBCTIME2=TTBCTIME1+TTBCTIMEINC
                     READ(38,'(A)') CDUM80
                     DO NumofNodes=1,NH
                       q_heat1(NumofNodes)=q_heat2(NumofNodes)
                       READ(38,*) NOD, q_heat2(NumofNodes)
                     END DO
                  END IF
               END DO
               TTBCRATIO=(TIMEIT-TTBCTIME1)/TTBCTIMEINC
               DO NumofNodes=1,NH
                 HFLUX(NumofNodes) = q_heat1(NumofNodes)+TTBCRATIO*
     &                (q_heat2(NumofNodes)-q_heat1(NumofNodes))
                  qsurf(NumofNodes)=-HFLUX(NumofNodes)
               END DO
            ELSE IF (BCFLAG_TEMP.EQ.2) THEN ! read in file with multiple values
               ALLOCATE (TMP(NP,6))
               OPEN(38,FILE=TRIM(INPUTDIR)//'/'//'fort.38')
               TTBCTIME1=TTBCSTATIM*86400.D0
               TTBCTIME2=TTBCTIME1+TTBCTIMEINC
               READ(38,'(A)') CDUM80
               READ(38,*) (NOD,(TMP(K,J),J=1,6),K=1,NP)
               DO NumofNodes=1,NP
                  q_heat1(NumofNodes)=-TMP(NumofNodes,1)-
     &               TMP(NumofNodes,2)+TMP(NumofNodes,3)-
     &               TMP(NumofNodes,5)+TMP(NumofNodes,4)-
     &               TMP(NumofNodes,6)
               END DO
               READ(38,'(A)') CDUM80
               READ(38,*) (NOD,(TMP(K,J),J=1,6),K=1,NP)
               DO NumofNodes=1,NP
                  q_heat2(NumofNodes)=-TMP(NumofNodes,1)-
     &               TMP(NumofNodes,2)+TMP(NumofNodes,3)-
     &               TMP(NumofNodes,5)+TMP(NumofNodes,4)-
     &               TMP(NumofNodes,6)
               END DO
               DO IT=1,ITHS
                  TIMEIT=IT*DTDPHS + STATIM*86400.D0
                  IF(TIMEIT.GT.TTBCTIME2) THEN
                     TTBCTIME1=TTBCTIME2
                     TTBCTIME2=TTBCTIME1+TTBCTIMEINC
                     DO NumofNodes=1,NP
                        q_heat1(NumofNodes)=q_heat2(NumofNodes)
                     END DO
                     READ(38,'(A)') CDUM80
                     READ(38,*) (NOD,(TMP(K,J),J=1,6),K=1,NP)
                     DO NumofNodes=1,NP
                        q_heat2(NumofNodes)=-TMP(NumofNodes,1)-
     &                     TMP(NumofNodes,2)+TMP(NumofNodes,3)-
     &                     TMP(NumofNodes,5)+TMP(NumofNodes,4)-
     &                     TMP(NumofNodes,6)
                     END DO
                  END IF
               END DO
               TTBCRATIO=(TIMEIT-TTBCTIME1)/TTBCTIMEINC
               DO NumofNodes=1,NP
                  HFLUX(NumofNodes) = q_heat1(NumofNodes)+TTBCRATIO*
     &                (q_heat2(NumofNodes)-q_heat1(NumofNodes))
                  qsurf(NumofNodes)=-HFLUX(NumofNodes)
               END DO
               DEALLOCATE(TMP)
            ELSE IF (BCFLAG_TEMP.EQ.3) THEN ! read in from file with 4 components 
               ALLOCATE (TMP(NP,4))
               OPEN(38,FILE=TRIM(INPUTDIR)//'/'//'fort.38')
               TTBCTIME1=TTBCSTATIM*86400.D0
               TTBCTIME2=TTBCTIME1+TTBCTIMEINC
               READ(38,'(A)') CDUM80
               READ(38,*) (NOD,(TMP(K,J),J=1,4),K=1,NP)
               DO NumofNodes=1,NP
                  q_heat1(NumofNodes)=-TMP(NumofNodes,4)-
     &               TMP(NumofNodes,3)+TMP(NumofNodes,1)+
     &               TMP(NumofNodes,2)
               END DO
               READ(38,'(A)') CDUM80
               READ(38,*) (NOD,(TMP(K,J),J=1,4),K=1,NP)
               DO NumofNodes=1,NP
                  q_heat2(NumofNodes)=-TMP(NumofNodes,4)-
     &               TMP(NumofNodes,3)+TMP(NumofNodes,1)+
     &               TMP(NumofNodes,2)
               END DO
               DO IT=1,ITHS
                  TIMEIT=IT*DTDPHS + STATIM*86400.D0
                  IF(TIMEIT.GT.TTBCTIME2) THEN
                     TTBCTIME1=TTBCTIME2
                     TTBCTIME2=TTBCTIME1+TTBCTIMEINC
                     DO NumofNodes=1,NP
                        q_heat1(NumofNodes)=q_heat2(NumofNodes)
                     END DO
                     READ(38,'(A)') CDUM80
                     READ(38,*) (NOD,(TMP(K,J),J=1,4),K=1,NP)
                     DO NumofNodes=1,NP
                        q_heat2(NumofNodes)=-TMP(NumofNodes,4)-
     &                     TMP(NumofNodes,3)+TMP(NumofNodes,1)+
     &                     TMP(NumofNodes,2)
                     END DO
                  END IF
               END DO
               TTBCRATIO=(TIMEIT-TTBCTIME1)/TTBCTIMEINC
               DO NumofNodes=1,NP
                  HFLUX(NumofNodes) = q_heat1(NumofNodes)+TTBCRATIO*
     &                (q_heat2(NumofNodes)-q_heat1(NumofNodes))
                  qsurf(NumofNodes)=-HFLUX(NumofNodes)
               END DO
               DEALLOCATE(TMP)
            END IF ! end if for heat flux file
         END IF ! end if for temperature boundary and heat flux
      ELSE IF ((RES_BC_FLAG.LT.0).AND.(CBAROCLINIC)) THEN 
         IF (BCFLAG_LNM.EQ.1) THEN 
           OPEN(35,FILE=TRIM(INPUTDIR)//'/'//'fort.35')
           READ(35,'(A)') CDUM80
           DO NumofBCNodes=1,NETA
             READ(35,*) NOD,LNM_BC(NumofBCNodes)
           END DO
         ELSE IF (BCFLAG_LNM.EQ.2) THEN 
           ! to be added later but currently information
           ! is taken in two ways from the HYCOM results
         ELSE IF (BCFLAG_LNM.EQ.3) THEN 
            DO NumofBCNodes=1,NETA
               IINDX=NBD(NumofBCNodes)
               LNM_BC(NumofBCNodes)=ETA2(IINDX)
            END DO 
         END IF 
      END IF

!kmd - end of additions

C.... Initialize the station 3D density output file 

      IF(IDen.EQ.1) IRType=1
      IF((IDen.EQ.2).OR.(IDen.EQ.3)) IRType=2
      IF(IDen.EQ.4) IRType=3

      IF(I3DSD.EQ.-1) THEN             !start a new ASCII file
        WRITE(16,497) 41
        OPEN(41,FILE=TRIM(LOCALDIR)//'/'//'fort.41')
        WRITE(41,499) RUNDes,RunID,AGrid
        WRITE(41,498) NDSet3DSD,NSta3DD,DTDP*NSpo3DSD,NSpo3DSD,NFEN,
     &                                                            IRType
        CLOSE(41)
        I3DSDRec=2
      ENDIF

      IF(I3DSD.EQ.-2) THEN             !start a new BINARY file
        WRITE(16,496) 41
        OPEN(41,FILE=TRIM(LOCALDIR)//'/'//'fort.41',ACCESS='DIRECT',
     &       RECL=NByte)
        I3DSDRec=0
        IF(NByte.EQ.4) THEN
          DO N=1,8
            I3DSDRec=I3DSDRec+1
            WRITE(41,REC=I3DSDRec) RDES4(N)
          ENDDO
          DO N=1,6
            I3DSDREC=I3DSDREC+1
            WRITE(41,REC=I3DSDRec) RID4(N)
          ENDDO
          DO N=1,6
            I3DSDRec=I3DSDRec+1
            WRITE(41,REC=I3DSDRec) AID4(N)
          ENDDO
        ENDIF
        IF(NByte.EQ.8) THEN
          DO N=1,4
            I3DSDRec=I3DSDRec+1
            WRITE(41,REC=I3DSDRec) RDES8(N)
          ENDDO
          DO N=1,3
            I3DSDRec=I3DSDRec+1
            WRITE(41,REC=I3DSDRec) RID8(N)
          ENDDO
          DO N=1,3
            I3DSDRec=I3DSDRec+1
            WRITE(41,REC=I3DSDRec) AID8(N)
          ENDDO
        ENDIF
        WRITE(41,REC=I3DSDRec+1) NDSet3DSD
        WRITE(41,REC=I3DSDRec+2) NSta3DD
        WRITE(41,REC=I3DSDRec+3) DTDP*NSpo3DSD
        WRITE(41,REC=I3DSDRec+4) NSpo3DSD
        WRITE(41,REC=I3DSDRec+5) NFEN   
        WRITE(41,REC=I3DSDRec+6) IRType
        I3DSDRec=I3DSDRec+6
        CLOSE(41)
      ENDIF

      IF(I3DSD.EQ.1) THEN          !find spot in existing ASCII file
        WRITE(16,497) 41
        OPEN(41,FILE=TRIM(INPUTDIR)//'/'//'fort.41')
        DO N=1,I3DSDRec
          READ(41,*) 
        ENDDO
        ENDFILE(41)
        CLOSE(41)
      ENDIF

      IF(I3DSD.EQ.2) THEN          !find spot in existing BINARY file
        WRITE(16,496) 41
        OPEN(41,FILE=TRIM(LOCALDIR)//'/'//'fort.41',ACCESS='DIRECT',
     &       RECL=NByte)
        CLOSE(41)
      ENDIF

C.... Initialize the 3D velocity station output file (Unit 42)

      IRType=3

      IF(I3DSV.EQ.-1) THEN             !start a new ASCII file
        WRITE(16,497) 42
        OPEN(42,FILE=TRIM(INPUTDIR)//'/'//'fort.42')
        WRITE(42,499) RUNDES,RUNID,AGRID
        WRITE(42,498) NDSet3DSV,NSta3DV,DTDP*NSpo3DSV,NSpo3DSV,NFEN,
     &                                                            IRType
        CLOSE(42)
        I3DSVRec=2
      ENDIF

      IF(ABS(I3DSV).EQ.-2) THEN        !start a new BINARY file
        WRITE(16,496) 42
        OPEN(42,FILE=TRIM(INPUTDIR)//'/'//'fort.42',ACCESS='DIRECT',
     &       RECL=NBYTE)
        I3DSVRec=0
        IF(NByte.EQ.4) THEN
          DO N=1,8
            I3DSVRec=I3DSVRec+1
            WRITE(42,REC=I3DSVRec) RDES4(N)
          ENDDO
          DO N=1,6
            I3DSVRec=I3DSVRec+1
            WRITE(42,REC=I3DSVRec) RID4(N)
          ENDDO
          DO N=1,6
            I3DSVRec=I3DSVRec+1
            WRITE(42,REC=I3DSVRec) AID4(N)
          ENDDO
        ENDIF
        IF(NByte.EQ.8) THEN
          DO N=1,4
            I3DSVRec=I3DSVRec+1
            WRITE(42,REC=I3DSVRec) RDES8(N)
          ENDDO
          DO N=1,3
            I3DSVRec=I3DSVRec+1
            WRITE(42,REC=I3DSVRec) RID8(N)
          ENDDO
          DO N=1,3
            I3DSVRec=I3DSVRec+1
            WRITE(42,REC=I3DSVRec) AID8(N)
          ENDDO
        ENDIF
        WRITE(42,REC=I3DSVRec+1) NDSet3DSV
        WRITE(42,REC=I3DSVRec+2) NSta3DV
        WRITE(42,REC=I3DSVRec+3) DTDP*NSpo3DSV
        WRITE(42,REC=I3DSVRec+4) NSpo3DSV
        WRITE(42,REC=I3DSVRec+5) NFEN
        WRITE(42,REC=I3DSVRec+6) IRType
        I3DSVRec=I3DSVRec+6
        CLOSE(42)
      ENDIF

      IF(I3DSV.EQ.1) THEN              !find spot in existing ASCII file
        WRITE(16,497) 42
        OPEN(42,FILE=TRIM(INPUTDIR)//'/'//'fort.42')
        DO N=1,I3DSVRec
          READ(42,*) 
        ENDDO
        ENDFILE(42)
        CLOSE(42)
      ENDIF

      IF(I3DSV.EQ.2) THEN              !find spot in existing BINARY file
        WRITE(16,496) 42
        OPEN(42,FILE=TRIM(INPUTDIR)//'/'//'fort.42',ACCESS='DIRECT',
     &       RECL=NByte)
        CLOSE(42)
      ENDIF

C.... Initialize the 3D turbulence station output file (Unit 43)

      IRType=3

      IF(I3DST.EQ.-1) THEN             !start a new ASCII file
        WRITE(16,497) 43
        OPEN(43,FILE=TRIM(LOCALDIR)//'/'//'fort.43')
        WRITE(43,499) RUNDES,RUNID,AGRID
        WRITE(43,498) NDSet3DST,NSta3DT,DTDP*NSpo3DST,NSpo3DST,NFEN,
     &                                                            IRType
        CLOSE(43)
        I3DSTRec=2
      ENDIF

      IF(ABS(I3DST).EQ.2) THEN         !start a new BINARY file
        WRITE(16,496) 43
        OPEN(43,FILE=TRIM(LOCALDIR)//'/'//'fort.43',ACCESS='DIRECT',
     &       RECL=NByte)
        I3DSTRec=0
        IF(NByte.EQ.4) THEN
          DO N=1,8
            I3DSTRec=I3DSTRec+1
            WRITE(43,REC=I3DSTRec) RDES4(N)
          ENDDO
          DO N=1,6
            I3DSTRec=I3DSTRec+1
            WRITE(43,REC=I3DSTRec) RID4(N)
          ENDDO
          DO N=1,6
            I3DSTRec=I3DSTRec+1
            WRITE(43,REC=I3DSTRec) AID4(N)
          ENDDO
        ENDIF
        IF(NByte.EQ.8) THEN
          DO N=1,4
            I3DSTRec=I3DSTRec+1
            WRITE(43,REC=I3DSTRec) RDES8(N)
          ENDDO
          DO N=1,3
            I3DSTRec=I3DSTRec+1
            WRITE(43,REC=I3DSTRec) RID8(N)
          ENDDO
          DO N=1,3
            I3DSTRec=I3DSTRec+1
            WRITE(43,REC=I3DSTRec) AID8(N)
          ENDDO
        ENDIF
        WRITE(43,REC=I3DSTRec+1) NDSet3DST
        WRITE(43,REC=I3DSTRec+2) NSta3DT
        WRITE(43,REC=I3DSTRec+3) DTDP*NSpo3DST
        WRITE(43,REC=I3DSTRec+4) NSpo3DST
        WRITE(43,REC=I3DSTRec+5) NFEN
        WRITE(43,REC=I3DSTRec+6) IRType
        I3DSTRec=I3DSTRec+6
        CLOSE(43)
      ENDIF

      IF(I3DST.EQ.1) THEN          !find spot in existing ASCII file
        WRITE(16,497) 43
        OPEN(43,FILE=TRIM(INPUTDIR)//'/'//'fort.43')
        DO N=1,I3DSTRec
          READ(43,*) 
        ENDDO
        ENDFILE(43)
        CLOSE(43)
      ENDIF

      IF(I3DST.EQ.2) THEN          !find spot in existing BINARY file
        WRITE(16,496) 43
        OPEN(43,FILE=TRIM(INPUTDIR)//'/'//'fort.43',ACCESS='DIRECT',
     &       RECL=NByte)
        CLOSE(43)
      ENDIF

C.... Initialize the global 3D density output file (Unit 44)
      IF(I3DGD.LT.0) THEN 
         N3DGD=0
         I3DGDRec=0
         IF((NTO3DGDS.LT.ITHS).AND.(NSpo3DGD.GT.0)) THEN 
           NTO3DGDS=NTO3DGDS+((ITHS-NTO3DGDS)/NSpo3DGD)*NSpo3DGD
           IF(NTO3DGDS.LT.ITHS) NTO3DGDS=NTO3DGDS+NSpo3DGD
           NDSet3DGD=(NTO3DGDF-NTO3DGDS)/NSpo3DGD
         END IF 
      END IF 

      IF(IDen.EQ.1) IRType=1
      IF((IDen.EQ.2).OR.(IDen.EQ.3)) IRType=2
      IF(IDen.EQ.4) IRType=3
C     jgf46.27 Replaced IRType with IDen
!kmd48.33bc - added in the information to make a file with all the
!             top temperature boundary condition information in
      IF(I3DGD.EQ.-1) THEN         !start a new ASCII file
        WRITE(16,497) 44
        OPEN(44,FILE=TRIM(LOCALDIR)//'/'//'fort.44')
        WRITE(44,499) RUNDES,RUNID,AGRID
        WRITE(44,498) NDSet3DGD,NP,DTDP*NSpo3DGD,NSpo3DGD,NFEN,IDen
        CLOSE(44)
        IF((IDEN.EQ.3).OR.(IDEN.EQ.4)) THEN
          IF (BCFLAG_TEMP.NE.0) THEN
            OPEN(47,FILE=TRIM(LOCALDIR)//'/'//'fort.47')
            WRITE(47,499) RUNDES,RUNID,AGRID
            WRITE(47,498) NDSet3DGD,NP,DTDP*NSpo3DGD,NSpo3DGD,1
            CLOSE(47)
          END IF
        END IF
        I3DGDRec=2
      ENDIF

      IF(I3DGD.EQ.-2) THEN         !start a new BINARY file
        WRITE(16,496) 44
        OPEN(44,FILE=TRIM(LOCALDIR)//'/'//'fort.44',ACCESS='DIRECT',
     &       RECL=NByte)
        I3DGDRec=0
        IF(NByte.EQ.4) THEN
          DO N=1,8
            I3DGDRec=I3DGDRec+1
            WRITE(44,REC=I3DGDRec) RDES4(N)
          ENDDO
          DO N=1,6
            I3DGDRec=I3DGDRec+1
            WRITE(44,REC=I3DGDRec) RID4(N)
          ENDDO
          DO N=1,6
            I3DGDRec=I3DGDRec+1
            WRITE(44,REC=I3DGDRec) AID4(N)
          ENDDO
        ENDIF
        IF(NByte.EQ.8) THEN
          DO N=1,4
            I3DGDRec=I3DGDRec+1
            WRITE(44,REC=I3DGDRec) RDES8(N)
          ENDDO
          DO N=1,3
            I3DGDRec=I3DGDRec+1
            WRITE(44,REC=I3DGDRec) RID8(N)
          ENDDO
          DO N=1,3
            I3DGDRec=I3DGDRec+1
            WRITE(44,REC=I3DGDRec) AID8(N)
          ENDDO
        ENDIF
        WRITE(44,REC=I3DGDRec+1) NDSet3DGD
        WRITE(44,REC=I3DGDRec+2) NP
        WRITE(44,REC=I3DGDRec+3) DTDP*NSpo3DGD
        WRITE(44,REC=I3DGDRec+4) NSpo3DGD
        WRITE(44,REC=I3DGDRec+5) NFEN
        WRITE(44,REC=I3DGDRec+6) IDEN
        I3DGDRec=I3DGDRec+6
        CLOSE(44)
      ENDIF

      IF(I3DGD.EQ.1) THEN          !find spot in existing ASCII file
        WRITE(16,497) 44
        OPEN(44,FILE=TRIM(INPUTDIR)//'/'//'fort.44')
        DO N=1,I3DGDRec
          READ(44,*) 
        ENDDO
        ENDFILE(44)
        CLOSE(44)
      ENDIF

!kmd - changed to fort.44 from fort.41
      IF(I3DGD.EQ.2) THEN          !find spot in existing BINARY file
        WRITE(16,496) 44
        OPEN(44,FILE=TRIM(INPUTDIR)//'/'//'fort.44',ACCESS='DIRECT',
     &       RECL=NByte)
        CLOSE(44)
      ENDIF
      
c.... Initialize the global 3D velocity output file (Unit 45)
      IF(I3DGV.LT.0) THEN
        N3DGV=0
        I3DGVRec=0
        IF((NTO3DGVS.LT.ITHS).AND.(NSpo3DGV.GT.0)) THEN
          NTO3DGVS=NTO3DGVS+((ITHS-NTO3DGVS)/NSpo3DGV)*NSpo3DGV
          IF(NTO3DGVS.LT.ITHS) NTO3DGVS=NTO3DGVS+NSpo3DGV
          NDSet3DGV=(NTO3DGVF-NTO3DGVS)/NSpo3DGV
        ENDIF
      ENDIF

      IF(I3DGV.EQ.-1) THEN         !start a new ASCII file
        WRITE(16,497) 45
        OPEN(45,FILE=TRIM(LOCALDIR)//'/'//'fort.45')
        WRITE(45,499) RUNDES,RUNID,AGRID
        WRITE(45,498) NDSet3DGV,NP,DTDP*NSpo3DGV,NSpo3DGV,NFEN,IRType
        CLOSE(45)
        I3DGVRec=2
      ENDIF

      IF(I3DGV.EQ.-2) THEN        !start a new BINARY file
        WRITE(16,496) 45
        OPEN(45,FILE=TRIM(LOCALDIR)//'/'//'fort.45',ACCESS='DIRECT',
     &       RECL=NBYTE)
        I3DGVRec=0
        IF(NByte.EQ.4) THEN
          DO N=1,8
            I3DGVRec=I3DGVRec+1
            WRITE(45,REC=I3DGVRec) RDES4(N)
          ENDDO
          DO N=1,6
            I3DGVRec=I3DGVRec+1
            WRITE(45,REC=I3DGVRec) RID4(N)
          ENDDO
          DO N=1,6
            I3DGVRec=I3DGVRec+1
            WRITE(45,REC=I3DGVRec) AID4(N)
          ENDDO
        ENDIF
        IF(NByte.EQ.8) THEN
          DO N=1,4
            I3DGVRec=I3DGVRec+1
            WRITE(45,REC=I3DGVRec) RDES8(N)
          ENDDO
          DO N=1,3
            I3DGVRec=I3DGVRec+1
            WRITE(45,REC=I3DGVRec) RID8(N)
          ENDDO
          DO N=1,3
            I3DGVRec=I3DGVRec+1
            WRITE(45,REC=I3DGVRec) AID8(N)
          ENDDO
        ENDIF
        WRITE(45,REC=I3DGVRec+1) NDSet3DGV
        WRITE(45,REC=I3DGVRec+2) NP
        WRITE(45,REC=I3DGVRec+3) DTDP*NSpo3DGV
        WRITE(45,REC=I3DGVRec+4) NSpo3DGV
        WRITE(45,REC=I3DGVRec+5) NFEN
        WRITE(45,REC=I3DGVRec+6) IRType
        I3DGVRec=I3DGVRec+6
        CLOSE(45)
      ENDIF

      IF(I3DGV.EQ.1) THEN          !find spot in existing ASCII file
        WRITE(16,497) 45
        OPEN(45,FILE=TRIM(INPUTDIR)//'/'//'fort.45')
        DO N=1,I3DGVRec
          READ(45,*) 
        ENDDO
        ENDFILE(45)
        CLOSE(45)
      ENDIF

      IF(I3DGV.EQ.2) THEN          !find spot in existing BINARY file
        WRITE(16,496) 45
        OPEN(45,FILE=TRIM(INPUTDIR)//'/'//'fort.45',ACCESS='DIRECT',
     &       RECL=NByte)
        CLOSE(45)
      ENDIF

c.... Initialize the global 3D turbulence output file (Unit 46)
      IF(I3DGT.LT.0) THEN
        N3DGT=0
        I3DGTRec=0
        IF((NTO3DGTS.LT.ITHS).AND.(NSpo3DGT.GT.0)) THEN
          NTO3DGTS=NTO3DGTS+((ITHS-NTO3DGTS)/NSpo3DGT)*NSpo3DGT
          IF(NTO3DGTS.LT.ITHS) NTO3DGTS=NTO3DGTS+NSpo3DGT
          NDSet3DGT=(NTO3DGTF-NTO3DGTS)/NSpo3DGT
        ENDIF
      ENDIF

      IF(I3DGT.EQ.-1) THEN         !start a new ASCII file
        WRITE(16,497) 46
        OPEN(46,FILE=TRIM(LOCALDIR)//'/'//'fort.46')
        WRITE(46,499) RUNDES,RUNID,AGRID
        WRITE(46,498) NDSet3DGT,NP,DTDP*NSpo3DGT,NSpo3DGT,NFEN,IRType
        CLOSE(46)
        I3DGTRec=2
      ENDIF

      IF(I3DGT.EQ.-2) THEN         !start a new BINARY file
        WRITE(16,496) 46
        OPEN(46,FILE=TRIM(LOCALDIR)//'/'//'fort.46',ACCESS='DIRECT',
     &       RECL=NByte)
        I3DGTRec=0
        IF(NByte.EQ.4) THEN
          DO N=1,8
            I3DGTRec=I3DGTRec+1
            WRITE(46,REC=I3DGTRec) RDES4(N)
          ENDDO
          DO N=1,6
            I3DGTRec=I3DGTRec+1
            WRITE(46,REC=I3DGTRec) RID4(N)
          ENDDO
          DO N=1,6
            I3DGTRec=I3DGTRec+1
            WRITE(46,REC=I3DGTRec) AID4(N)
          ENDDO
        ENDIF
        IF(NByte.EQ.8) THEN
          DO N=1,4
            I3DGTRec=I3DGTRec+1
            WRITE(46,REC=I3DGTRec) RDES8(N)
          ENDDO
          DO N=1,3
            I3DGTRec=I3DGTRec+1
            WRITE(46,REC=I3DGTRec) RID8(N)
          ENDDO
          DO N=1,3
            I3DGTRec=I3DGTRec+1
            WRITE(46,REC=I3DGTRec) AID8(N)
          ENDDO
        ENDIF
        WRITE(46,REC=I3DGTRec+1) NDSet3DGT
        WRITE(46,REC=I3DGTRec+2) NP
        WRITE(46,REC=I3DGTRec+3) DTDP*NSpo3DGT
        WRITE(46,REC=I3DGTRec+4) NSpo3DGT
        WRITE(46,REC=I3DGTRec+5) NFEN
        WRITE(46,REC=I3DGTRec+6) IRType
        I3DGTRec=I3DGTRec+6
        CLOSE(46)
      ENDIF

      IF(I3DGT.EQ.1) THEN          !find spot in existing ASCII file
        WRITE(16,497) 46
        OPEN(46,FILE=TRIM(INPUTDIR)//'/'//'fort.46')
        DO N=1,I3DGTRec
          READ(46,*) 
        ENDDO
        ENDFILE(46)
        CLOSE(46)
      ENDIF

      IF(I3DGT.EQ.2) THEN          !find spot in existing BINARY file
        WRITE(16,496) 46
        OPEN(46,FILE=TRIM(INPUTDIR)//'/'//'fort.46',ACCESS='DIRECT',
     &       RECL=NByte)
        CLOSE(46)
      ENDIF

C.... Set up a few final odds and ends for a 3D run

      CALL VSSTUP ()
        
      RETURN
      END
