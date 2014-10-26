C******************************************************************************
C PADCIRC VERSION 45.12 03/17/2006                                            *
C  last changes in this file VERSION 45.12                                    *
C                                                                             *
C                                                                             *
C This module handles the model initialization for a cold start.  The primary *
C 2D coldstart initialization is handled in SUBROUTINE COLDSTART.  The primary*
C 3D initialization is handled in SUBROUTINE COLDSTART_3D.                    *
C                                                                             *
C******************************************************************************

C
      SUBROUTINE COLDSTART()
C
      USE SIZES
      USE GLOBAL
      USE GLOBAL_3DVS, ONLY: I3DSD, I3DSV, I3DST, I3DGD, I3DGV, I3DGT
      USE MESH, ONLY : NE, NP, X, Y, SLAM, SFEA, ICS, DP, NM, MJU, AID4
      USE BOUNDARIES, ONLY : NETA, NFLUXF, NOPE, NVEL, LBCODEI
      USE HARM
      USE WIND
      USE GLOBAL_IO

      USE ADCIRC_MOD, ONLY : ADCIRC_Terminate
      USE OWIWIND,ONLY : NWS12INIT,NWS12GET    ! sb added 09/xx/2006
      USE OWI_ICE, ONLY : NCICE1_INIT, NCICE1_GET  !tcm v49.64.01
      USE RS2,ONLY : RS2INIT,RS2GET            ! sb added 09/xx/2006
      USE NodalAttributes, ONLY :
     &         STARTDRY, GeoidOffset, LoadGeoidOffset,
     &         OutputTau0, LoadRiver_et_WSE, River_et_WSE !jgf47.06

      IMPLICIT NONE

      INTEGER I, J,IOS, IJK
      INTEGER NM1, NM2, NM3
      INTEGER NC1, NC2, NC3, NCELE

      REAL(SZ) HTOT
      REAL(SZ) HollandTime, AsymmetricTime, GeneralAsymTime
      REAL(SZ) :: PRBCKGRND_MH2O  !tcm v49.16 20100617
!kmd - Evan's updates for rivers above MSL
!      REAL(SZ),ALLOCATABLE :: et_WSE(:)
      LOGICAL FOUND

C     jgf48.4628 Add capability to turn off solution if only met output
C     was requested.
      LOGICAL metOutputOn     ! .true. if met output was requested
      LOGICAL nonMetOutputOff ! .true. no other output was requested

c.... tcm v50.66.02 additions for time varying bathymetry
      CHARACTER(80) :: CDUM80
C
      call setMessageSource("coldstart")
#if defined(COLDSTART_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      metOutputOn = .FALSE.
      nonMetOutputOff = .FALSE.

C------------------------------------------------------------------------------
C... tcm v50.66.01 Time Varying Bathymetry
C------------------------------------------------------------------------------
C
      IF (NDDT.NE.0) then
C...    Set the first two bathymetry values, these arrays are
C...    allocated in alloc_main13 called during read_input
         DP1(:) = DP(:)  !set dp1 to initial depth from fort.14 grid
         DP2(:) = DP(:)  !set dp2 to the same values

         BTIME1 = STATIM*86400.D0
         BTIME2=BTIME2+BTIMINC  !new ending time for this record
         BTIME_END = BTIME1 + BCHGTIMINC  

C...    read in the second depth dp2 for btime2
         IOS = 0
         CALL openFileForRead(141,TRIM(INPUTDIR)//'/'//'fort.141',IOS)
         IF (IOS.GT.0) THEN
            CALL ADCIRC_Terminate()
         ENDIF

         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
         WRITE(16,1112)
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1976)
         WRITE(16,1976)
 1976    FORMAT(/,1X,'TIME VARYING BATHYMETRY INFORMATION ',
     &        'READ FROM UNIT 141',/)

         IF (ABS(NDDT).EQ.1) THEN
            DO I=1,NP
               READ(141,*) IJK,DP2(IJK)
            ENDDO
         ENDIF
         IF (ABS(NDDT).EQ.2) THEN
!!!        go get new record for only some nodes, all
!!!        other nodes keep their current value
            CALL NDDT2GET(141,DP2(:),-99999.d0)
         ENDIF

C...     IF WETTING AND DRYING WILL NOT BE USED, MAKE SURE ALL BATHYMETRIC
C...     DEPTHS ARE > OR = TO H0.
         IF ((NOLIFA.EQ.0).OR.(NOLIFA.EQ.1)) THEN
            DO I=1,NP
               IF (DP2(I).LT.H0) DP2(I) = H0
            ENDDO
         ENDIF
      ENDIF  !end test for time varying bathymetry
C
C------------------------------------------------------------------------------
C

! kmd - Evan's updates for rivers above MSL
! open the fort.88 files if needed
!      INQUIRE(FILE=TRIM(INPUTDIR)//'/'//'fort.88',EXIST=FOUND)
!      IF (FOUND) THEN
!         ALLOCATE (et_WSE(NP) )
!         OPEN (88, FILE=TRIM(INPUTDIR)//'/'//'fort.88')
!         DO I=1,NP
!            READ(88,*) et_WSE(i)
!         END DO
!         River_above_MSL=.TRUE.
!      END IF

! tcm v51.27 20140502 Changed from using a fort.88 file
! to a nodal attribute
      IF (LoadRiver_et_WSE) then
         River_above_MSL=.TRUE.
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
         WRITE(16,1112)
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,2788)
         WRITE(16,2788)
 2788    FORMAT(/,1X,'INITIAL RIVER ELEVATIONS WERE SPECIFIED AS A NODAL',
     &        /,1x,'ATTRIBUTE AND WILL BE USED.',/)         
      ENDIF

      ! tcm v51.27 Check to see if a fort.88 file is supplied
      ! if so issue a warning.      
      INQUIRE(FILE=TRIM(INPUTDIR)//'/'//'fort.88',EXIST=FOUND)
      IF (FOUND) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
         WRITE(16,1112)
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,2789)
         WRITE(16,2789)
 2789    FORMAT(/,1X,'--------WARNING--WARNING--WARNING---------',
     &         /,1x,'A FORT.88 FILE WAS FOUND BUT IS NO LONGER SUPPORTED.',
     &        /,1x,'INITIAL RIVER ELEVATIONS MUST BE SET AS A NODAL ATTRIBUTE.',/)
      ENDIF

C
C.... tcm v49.16 20100617
C.... convert background pressure from millibars to meters of water
      PRBCKGRND_MH2O = 100.0D0*PRBCKGRND/(RHOWAT0*G)

C
C...  SET AT REST INITIAL CONDITION OVER WHOLE DOMAIN
C...  IF BOTTOM IS ABOVE THE GEIOD -> DRY NODE
C...  IF BOTTOM IS INITIALLY BELOW THE GEOID AND STARTDRY=1 -> DRY NODE
C...
      DO I=1,NP
         UU1(I) =0.D0
         VV1(I) =0.D0
         UU2(I) =0.D0
         VV2(I) =0.D0
         QX1(I) =0.D0
         QY1(I) =0.D0
         QX2(I) =0.D0
         QY2(I) =0.D0
C     jgf46.01 Added the ability to include steric effects.
         IF (LoadGeoidOffset) THEN
            ETA2(I)=GeoidOffset(I)
         ELSE
            ETA2(I)=0.D0
         ENDIF
! kmd - Evan's updates for rivers above MSL
! It is only necessary to read in and reset the ETA2 values
! in this case for the cold start conditions.
!         IF (River_above_MSL) THEN
!            IF (LoadGeoidOffset) THEN
!               IF (et_WSE(i).GT.GeoidOffset(I)) THEN
!                  ETA2(I)=et_WSE(i)
!               END IF
!            ELSE
!               IF (et_WSE(i).GT.0.d0) THEN
!                  ETA2(I)=et_WSE(i)
!               END IF
!            END IF
!         END IF
! tcm v51.27 20140502 Made the GeoidOffset and Initial
!            River Elevation an additive quantitiy 
!            instead of the max of the two         
         IF (River_above_MSL) then
           IF (LoadGeoidOffSet) then
              IF (River_et_WSE(I).GT.0.d0)
     &                ETA2(I) = River_et_WSE(I) + GeoidOffset(I)
           ELSE
              IF (River_et_WSE(I).GT.0.d0) ETA2(I) = River_et_WSE(I)
           ENDIF
         ENDIF


         NNODECODE(I)=1
         IF(NOLIFA.EQ.2) THEN
            HTOT=DP(I)+ETA2(I)
            IF(HTOT.LE.H0) THEN
               NNODECODE(I)=0
               ETA2(I)=H0-DP(I)
            ELSE
               IF(STARTDRY(I).EQ.1) THEN
                  NNODECODE(I)=0
                  ETA2(I)=H0-DP(I)
               ENDIF
            ENDIF
         ENDIF
         ETA1(I)=ETA2(I)
         ETAS(I)=0.D0
         CH1(I)=0.d0
      END DO
C...
C...  INITIALIZE DRY ELEMENT CODE
C...
      DO I=1,NE
         NOFF(I)=1
      END DO
C...
C...  INITIALIZE THE ELEVATION SPECIFIED BOUNDARY CONDITION IF IT
C...  REQUIRES THE USE OF THE UNIT 19 FILE.
C...

      IF((NOPE.GT.0).AND.(NBFR.EQ.0)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
         WRITE(16,1112)
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1977)
         WRITE(16,1977)
 1977    FORMAT(/,1X,'ELEVATION SPECIFIED INFORMATION READ FROM UNIT ',
     &        '19',/)
         OPEN(19,FILE=TRIM(INPUTDIR)//'/'//'fort.19')
         READ(19,*) ETIMINC
         DO J=1,NETA
            READ(19,*) ESBIN1(J)
         END DO
         DO J=1,NETA
            READ(19,*) ESBIN2(J)
         END DO
         ETIME1 = STATIM*86400.D0
         ETIME2 = ETIME1 + ETIMINC
      ENDIF
C
C....INITIALIZE THE NORMAL FLOW BOUNDARY CONDITION
C
      DO I=1,NVEL
         QN2(I)=0.D0
         QN1(I)=0.D0
         QN0(I)=0.D0
      END DO

      IF((NFLUXF.EQ.1).AND.(NFFR.EQ.0)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
         WRITE(16,1112)
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1979)
         WRITE(16,1979)
 1979    FORMAT(/,1X,'NORMAL FLOW INFORMATION READ FROM UNIT 20',/)
         OPEN(20,FILE=TRIM(INPUTDIR)//'/'//'fort.20')
         READ(20,*) FTIMINC
         DO J=1,NVEL
            QNIN1(J)=0.D0
            IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &          .OR.(LBCODEI(J).EQ.22).OR.(LBCODEI(J).EQ.32)) THEN
               READ(20,*) QNIN1(J)
            ENDIF
         END DO
         DO J=1,NVEL
            QNIN2(J)=0.D0
            IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &           .OR.(LBCODEI(J).EQ.22).OR.(LBCODEI(J).EQ.32)) THEN
               READ(20,*) QNIN2(J)
            ENDIF
         END DO
         QTIME1 = STATIM*86400.D0
         QTIME2 = QTIME1 + FTIMINC
      ENDIF

C...INPUT METEOROLOGICAL INFORMATION FROM UNIT 22 OR UNIT 200 SERIES
C....IF FLEET NUMERIC WIND DATA IS USED, FIND BEGINNING TIME IN FILE,
C....NOTE: CAN'T DEAL WITH WIND THAT STARTS AFTER WREFTIM!!!!!!!!!!!!
C....READ IN AND INTERPOLATE IN SPACE ONTO THE ADCIRC GRID THE
C....TIME LEVEL 1 AND LEVEL 2 WIND FIELDS

      DO I=1,NP
         WSX1(I)=0.D0
         WSY1(I)=0.D0
!        PR1(I)=0.D0   !tcm v49.16 20100617
         PR1(I) = PRBCKGRND_MH2O  !tcm v49.16 20100617 added
         WSX2(I)=0.D0
         WSY2(I)=0.D0
!        PR2(I)=0.D0   !tcm v49.16 20100617
         PR2(I) = PRBCKGRND_MH2O  !tcm v49.16 20100617 added
      ENDDO

      IF(NWS.NE.0) THEN
C        jgf48.4627 If user only requests meteorological output, then
C        set a flag that turns off the calculation of the GWCE and
C        momentum equations
         IF ((NOUTE.EQ.0) .AND. (NOUTV.EQ.0) .AND. (NOUTC.EQ.0) .AND.
     &       (NOUTGE.EQ.0) .AND. (NOUTGV.EQ.0) .AND. (NOUTGC.EQ.0 ).AND.
     &       (I3DSD.EQ.0) .AND. (I3DSV.EQ.0) .AND. (I3DST.EQ.0) .AND.
     &       (I3DGD.EQ.0) .AND. (I3DGV.EQ.0) .AND. (I3DGT.EQ.0) .AND.
     &       (NHASE.EQ.0) .AND. (NHASV.EQ.0) .AND.
     &       (NHAGE.EQ.0) .AND. (NHAGV.EQ.0) )
     &       THEN
            nonMetOutputOff = .TRUE.
         ENDIF
         IF ( (NOUTM.NE.0) .OR. (NOUTGW.NE.0) ) THEN
            metOutputOn = .TRUE.
         ENDIF
         IF ( metOutputOn .AND. nonMetOutputOff ) THEN
            IF ((NSCREEN.NE.0).AND.MYPROC.EQ.0) THEN
               WRITE(ScreenUnit,*)
     &      'INFO: Only meterological output was requested. ADCIRC'
               WRITE(ScreenUnit,*)
     &      'will not solve the GWCE or momentum equations.'
             ENDIF
             WRITE(16,*)
     &      'INFO: Only meterological output was requested. ADCIRC'
             WRITE(16,*)
     &      'will not solve the GWCE or momentum equations.'
            METONLY = .TRUE. ! set flag
         ENDIF

         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
         WRITE(16,1112)
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1980)
         WRITE(16,1980)
 1980    FORMAT(/,1X,'WIND (AND PRESSURE) INFORMATION READ.',/)
      ENDIF

      IF(NWS.EQ.1) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
      ENDIF

      IF(ABS(NWS).EQ.2) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
      ENDIF

      IF(NWS.EQ.3) THEN
C....    tcm_v49.04 20091124 -- Changed from local inputdir to global inputdir
C....                           to correspond with having only a global wind file.
!         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22')

 2222    CALL NWS3GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,IWTIME,IWYR,WTIMED,NP,
     &        NWLON,NWLAT,WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,
     &        NScreen, ScreenUnit)
         IF(IWYR.NE.IREFYR) THEN
            IWTIMEP=IWTIME
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
            END DO
            GOTO 2222
         ENDIF
         IF(WTIMED.LE.WREFTIM) THEN
            IWTIMEP=IWTIME
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
            END DO
            GOTO 2222
         ENDIF
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0)
     &        WRITE(ScreenUnit,*)'FOUND WIND DATA AT TIME= ',IWTIMEP
         WRITE(16,*) 'FOUND WIND DATA AT TIME= ',IWTIMEP
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0)
     &        WRITE(ScreenUnit,*)'FOUND WIND DATA AT TIME= ',IWTIME
         WRITE(16,*) 'FOUND WIND DATA AT TIME= ',IWTIME
         WTIME2=WTIMED-WREFTIM  !CAST INTO MODEL TIME REFERENCE
         WTIME1=WTIME2-WTIMINC
      ENDIF

      IF(ABS(NWS).EQ.4) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         WTIME2=WTIME1+WTIMINC
         CALL NWS4GET(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
         CALL NWS4GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
      ENDIF

      IF(ABS(NWS).EQ.5) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
      ENDIF

      IF(NWS.EQ.6) THEN
C....    tcm_v49.04 20091124 -- Changed from local inputdir to global inputdir
C....                           to correspond with having only a global wind file.
!         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22')
         CALL NWS6GET(X,Y,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         CALL NWS6GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
      ENDIF

C     jgf46.00 Added option to directly apply surface stress without any
C     other correction factors.
      IF(ABS(NWS).EQ.7) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
      ENDIF

C     jgf46.02 New option to read in hurricane locations and generate
C     generate hurricane winds from the Holland Wind Model.
      IF(ABS(NWS).EQ.8) THEN
         HollandTime = STATIM*86400.D0
         CALL HollandGet(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,
     &     ICS,RHOWAT0,G,HollandTime,NSCREEN,ScreenUnit)
      ENDIF

C     RJW Merged:
!    rjw added nws = 19: asymmetric hurricane winds v2.0
      IF(NWS.EQ.19) THEN
         AsymmetricTime = STATIM*86400.D0
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22',STATUS='OLD')
        CALL NWS19GET(SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,AsymmetricTime,ICS)
      ENDIF
      
!    jie added nws = 20: generalized asymmetric vortex model
      IF(NWS.EQ.20) THEN
         GeneralAsymTime = STATIM*86400.D0
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22',STATUS='OLD')
        CALL NWS20GET(SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,GeneralAsymTime,ICS)
      ENDIF      

      ! jgf49.1001 Added NWS=29, asymmetric vortex wind (NWS19) embedded in
      ! OWI (NWS12) basin scale wind field.
      IF(NWS.EQ.29) THEN
         CALL NWS12INIT(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
         CALL NWS12GET(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
         CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
         ! now blend in the vortex met field
         AsymmetricTime = STATIM*86400.D0
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'NWS_19_fort.22',
     &        STATUS='OLD')
         ALLOCATE(vortexWVNX2(NP),vortexWVNY2(NP),vortexPRN2(NP))
         CALL NWS19GET(SLAM,SFEA,vortexWVNX2,vortexWVNY2,vortexPRN2,
     &                 NP,AsymmetricTime,ICS)
      ENDIF

      IF(NWS.EQ.10) THEN
         WTIME1=STATIM*86400.D0
         WTIME2=WTIME1+WTIMINC
         NWSGGWI=-1
         CALL NWS10GET(NWSGGWI,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,RHOWAT0,G,
     &        NWLON,NWLAT,WTIMINC) !JUST COMPUTE INTERPOLATING FACTORS
         NWSGGWI=1
         CALL NWS10GET(NWSGGWI,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,RHOWAT0,G,
     &        NWLON,NWLAT,WTIMINC) !NOW INTERPOLATE 1st WIND FIELD
      ENDIF

      IF(NWS.EQ.11) THEN
         WTIME1=STATIM*86400.D0
         WTIME2=WTIME1+WTIMINC
         NWSEGWI=0
         IDSETFLG=0
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1197)
         WRITE(16,1197)
 1197    FORMAT(/,1X,'THE E29 MET GRID INTERPOLATING FACTORS ARE ',
     &        'BEING COMPUTED ')
         CALL NWS11GET(NWSEGWI,IDSETFLG,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,
     &        RHOWAT0,G)        !JUST COMPUTE INTERPOLATING FACTORS
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1198)
         WRITE(16,1198)
 1198    FORMAT(1X,'FINISHED COMPUTING E29 INTERPOLATING FACTORS',/)
         NWSEGWI=1
         IDSETFLG=1
         CALL NWS11GET(NWSEGWI,IDSETFLG,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,
     &        RHOWAT0,G)        !NOW INTERPOLATE 1st WIND FIELD
      ENDIF

C.....sb46.28sb01 NWS=-12 and 12 was added to deal with raw OWI files.  09/xx/2006
      IF(ABS(NWS).EQ.12) THEN
        CALL NWS12INIT(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
        CALL NWS12GET(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
        CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
        WTIME1 = STATIM*86400.D0
        WTIME2 = WTIME1 + WTIMINC
      ENDIF

      !
      ! jgf50.38.05: Added NWS=15,-15 for reading HWind data
      IF (ABS(NWS).EQ.15) THEN
         CALL NWS15INIT(STATIM*86400.0)
      ENDIF

C     
C....tcm v51.06.02: Added NWS=16,-16 for reading ASCII GFDL Met data
      IF (ABS(NWS).EQ.16) THEN
         CALL INIT_GFDL(STATIM*86400.0)
      ENDIF


C...INPUT RADIATION STRESS INFORMATION FROM UNIT 23
C....READ IN THE TIME LEVEL 1 AND LEVEL 2 FIELDS

      ! NRS=2 was added. 09/xx/2006 sb
      IF(NRS.GE.1) THEN ! sb46.28sb03
         IF(NWS.EQ.0) THEN
            DO I=1,NP
               WSX1(I)=0.D0
               WSY1(I)=0.D0
               WSX2(I)=0.D0
               WSY2(I)=0.D0
               PRN1(I)=0.D0     !need to be initialized
               PRN2(I)=0.D0     !even if not used
            ENDDO
         ENDIF
         RSTIME1 = STATIM*86400.D0
         RSTIME2 = RSTIME1+RSTIMINC
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
#ifdef CSWAN
Casey 090302: Added this call to initialize the radiation stress gradients.
         IF(NRS.EQ.3) THEN
           IF(.NOT.ALLOCATED(RSNX1)) ALLOCATE(RSNX1(1:NP))
           IF(.NOT.ALLOCATED(RSNX2)) ALLOCATE(RSNX2(1:NP))
           IF(.NOT.ALLOCATED(RSNY1)) ALLOCATE(RSNY1(1:NP))
           IF(.NOT.ALLOCATED(RSNY2)) ALLOCATE(RSNY2(1:NP))
           DO I=1,NP
             RSNX1(I) = 0.D0
             RSNX2(I) = 0.D0
             RSNY1(I) = 0.D0
             RSNY2(I) = 0.D0
           ENDDO
         ENDIF
#endif
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
         WRITE(16,1112)
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1981)
         WRITE(16,1981)
 1981    FORMAT(/,1X,'RADIATION STRESS INFORMATION READ.',/)
      ENDIF

C... v49.64.01 - tcm added
C...INPUT ICE CONCENTRATION INFORMATION FROM UNIT 25
C....READ IN THE TIME LEVEL 1 AND LEVEL 2 FIELDS
      IF (NCICE.GE.1) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
         WRITE(16,1112)
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,2980)
         WRITE(16,2980)
 2980    FORMAT(/,1X,'ICE CONCENTRATION INFORMATION READ.',/)

         IF(NCICE.EQ.12) THEN
            CALL NCICE1_INIT(CICE1,NP)
            CALL NCICE1_GET(CICE1,NP)
            CALL NCICE1_GET(CICE2,NP)
            CICE_TIME1 = STATIM*86400.D0
            CICE_TIME2 = CICE_TIME1 + CICE_TIMINC
         ENDIF
      ENDIF



C...
C...LINES TO USE TIDAL POTENTIAL FORCING
C...
      if (CTIP) then
         DO I=1,NP
            TIP2(I)=0.0
         END DO
      endif
CWET...
CWET...THE FOLLOWING LINES ARE FOR WETTING AND DRYING
CWET...Dry any landlocked nodes by checking that they are connected to at
CWET...least 1 functioning element.
CWET...
      IF(NOLIFA.EQ.2) THEN
         DO I=1,NP
            MJU(I)=0
         ENDDO
         DO I=1,NE
            NM1=NM(I,1)
            NM2=NM(I,2)
            NM3=NM(I,3)
            NC1=NNODECODE(NM1)
            NC2=NNODECODE(NM2)
            NC3=NNODECODE(NM3)
            NCELE=NC1*NC2*NC3*NOFF(I)
            MJU(NM1)=MJU(NM1)+NCELE
            MJU(NM2)=MJU(NM2)+NCELE
            MJU(NM3)=MJU(NM3)+NCELE
         ENDDO
         DO I=1,NP
            IF((NNODECODE(I).EQ.1).AND.(MJU(I).EQ.0)) THEN
               NNODECODE(I)=0
               IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9883) I
               WRITE(16,9883) I
            ENDIF
         ENDDO
      ENDIF

C....INITIALIZE ELEVATION STATION SPOOL COUNTER
C....OPEN ELEVATION STATION OUTPUT FILE
C....WRITE OUT HEADER INFORMATION INCLUDING NTRSPE (NO. OF DATA PTS. AT EACH
C....ELEVATION STATION), NSTAE, DT*NSPOOLE, NSPOOLE, IRTYPE
C...
      NSCOUE=0
      IESTP=0

 3220 FORMAT(1X,A32,2X,A24,2X,A24)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)

      IF(ABS(NOUTE).EQ.1) THEN
         CALL OPEN_GBL_FILE(61, TRIM(GLOBALDIR)//'/'//'fort.61',
     &    NSTAE_G, NSTAE, HEADER61)
         IESTP=2
      ENDIF

      IF(ABS(NOUTE).EQ.2) THEN
         OPEN(61,FILE=TRIM(LOCALDIR)//'/'//'fort.61',
     &           ACCESS='DIRECT',RECL=NBYTE)
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
         CLOSE(61)
         ENDIF

C...
C....INITIALIZE VELOCITY STATION SPOOL COUNTER
C....OPEN VELOCITY STATION OUTPUT FILE
C....WRITE OUT HEADER INFORMATION INCLUDING NTRSPV (NO. OF DATA PTS. AT EACH
C....VELOCITY STATION), NSTAV, DT*NSPOOLV, NSPOOLV, IRTYPE
C...
      NSCOUV=0
      IVSTP=0

      IF(ABS(NOUTV).EQ.1) THEN
         CALL OPEN_GBL_FILE(62, TRIM(GLOBALDIR)//'/'//'fort.62',
     $    NSTAV_G, NSTAV, HEADER62)
         IVSTP=2
      ENDIF

      IF(ABS(NOUTV).EQ.2) THEN
         OPEN(62,FILE=TRIM(LOCALDIR)//'/'//'fort.62',
     &           ACCESS='DIRECT',RECL=NBYTE)
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
         CLOSE(62)
         ENDIF

C...
C....INITIALIZE CONCENTRATION STATION SPOOL COUNTER
C....OPEN ELEVATION STATION OUTPUT FILE
C....WRITE OUT HEADER INFORMATION INCLUDING NTRSPC (NO. OF DATA PTS. AT EACH
C....CONCENTRATION STATION), NSTAC, DT*NSPOOLC, NSPOOLC, IRTYPE
C...
      NSCOUC=0
      ICSTP=0

      IF(ABS(NOUTC).EQ.1) THEN
         CALL OPEN_GBL_FILE(81, TRIM(GLOBALDIR)//'/'//'fort.81',
     $    NSTAC_G, NSTAC, HEADER81)
         ICSTP=2
      ENDIF

      IF(ABS(NOUTC).EQ.2) THEN
         OPEN(81,FILE=TRIM(LOCALDIR)//'/'//'fort.81',
     &           ACCESS='DIRECT',RECL=NBYTE)
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
         CLOSE(81)
         ENDIF

C...
C....INITIALIZE METEOROLOGICAL STATION SPOOL COUNTERS
C....OPEN METEOROLOGICAL STATION OUTPUT FILES
C....WRITE OUT HEADER INFORMATION INCLUDING NTRSPM (NO. OF DATA PTS. AT EACH
C....METEOROLOGICAL STATION), NSTAM, DT*NSPOOLM, NSPOOLM, IRTYPE
C...
      NSCOUM=0
      IPSTP=0
      IWSTP=0
      IICESTP = 0  ! v49.64.01 tcm -added ice stations

      IF(ABS(NOUTM).EQ.1) THEN
         CALL OPEN_GBL_FILE(71, TRIM(GLOBALDIR)//'/'//'fort.71',
     $    NSTAM_G, NSTAM, HEADER71)
         IPSTP=2
         CALL OPEN_GBL_FILE(72, TRIM(GLOBALDIR)//'/'//'fort.72',
     $     NSTAM_G, NSTAM, HEADER72)
         IWSTP=2
! V49.64.01 TCM -ADDED ICE STATIONS
         IF(NCICE.NE.0) THEN
            CALL OPEN_GBL_FILE(91, TRIM(GLOBALDIR)//'/'//'fort.91',
     $     NSTAM_G, NSTAM, HEADER91)
            IICESTP=2
         ENDIF
      ENDIF

      IF(ABS(NOUTM).EQ.2) THEN
         OPEN(71,FILE=TRIM(LOCALDIR)//'/'//'fort.71',
     &           ACCESS='DIRECT',RECL=NBYTE)
         OPEN(72,FILE=TRIM(LOCALDIR)//'/'//'fort.72',
     &           ACCESS='DIRECT',RECL=NBYTE)
! V49.64.01 TCM -ADDED ICE STATIONS
         IF(NCICE.NE.0) THEN
            OPEN(91,FILE=TRIM(LOCALDIR)//'/'//'fort.91',
     &           ACCESS='DIRECT',RECL=NBYTE)
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
         CLOSE(71)
         CLOSE(72)
!  v49.64.01 tcm - added ice station support
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

C....tcm v50.66.01 Addition for time varying bathymetry
C...
C....INITIALIZE BATHYMETRY STATION SPOOL COUNTER
C....THESE WILL BE SAVED AT THE SAME LOCATIONS AS ELEVATION STATIONS 
C....AND ARE CONTROLLED BY THE ELEVATION STATION CONTROLS
C....OPEN BATHYMETRY STATION OUTPUT FILE
C....WRITE OUT HEADER INFORMATION INCLUDING NTRSPE (NO. OF DATA PTS. AT EACH
C....BATHYMETRY STATION), NSTAE, DT*NSPOOLE, NSPOOLE, IRTYPE
C...
      NSCOUB=0
      IBSTP=0

      IF ((NDDT.NE.0).AND.(ABS(NOUTE).EQ.1)) THEN
         CALL OPEN_GBL_FILE(75, TRIM(GLOBALDIR)//'/'//'fort.75',
     &    NSTAE_G, NSTAE, HEADER61)  !no need to create a header75
         IBSTP=2
      ENDIF

      IF ((NDDT.NE.0).AND.(ABS(NOUTE).EQ.2)) THEN
         OPEN(75,FILE=TRIM(LOCALDIR)//'/'//'fort.75',
     &           ACCESS='DIRECT',RECL=NBYTE)
         IF(NBYTE.EQ.4) THEN
            DO I=1,8
               WRITE(75,REC=IBSTP+I) RDES4(I)
            ENDDO
            IBSTP=IBSTP+8
            DO I=1,6
               WRITE(75,REC=IBSTP+I) RID4(I)
            ENDDO
            IBSTP=IBSTP+6
            DO I=1,6
               WRITE(75,REC=IBSTP+I) AID4(I)
            ENDDO
            IBSTP=IBSTP+6
         ENDIF
         IF(NBYTE.EQ.8) THEN
            DO I=1,4
               WRITE(75,REC=IBSTP+I) RDES8(I)
            ENDDO
            IBSTP=IBSTP+4
            DO I=1,3
               WRITE(75,REC=IBSTP+I) RID8(I)
            ENDDO
            IBSTP=IBSTP+3
            DO I=1,3
               WRITE(75,REC=IBSTP+I) AID8(I)
            ENDDO
            IBSTP=IBSTP+3
         ENDIF
         WRITE(75,REC=IBSTP+1) NTRSPE
         WRITE(75,REC=IBSTP+2) NSTAE
         WRITE(75,REC=IBSTP+3) DT*NSPOOLE
         WRITE(75,REC=IBSTP+4) NSPOOLE
         WRITE(75,REC=IBSTP+5) 1
         IBSTP=IBSTP+5
         CLOSE(75)
      ENDIF

C...
C....INITIALIZE GLOBAL ELEVATION SPOOL COUNTER
C....OPEN GLOBAL ELEVATION OUTPUT FILE
C....WRITE OUT HEADER INFORMATION INCLUDING NDSETSE
C....(NO. OF GLOBAL ELEVATION DATA SETS TO BE SPOOLED),
C....NP, DT*NSPOOLGE, NSPOOLGE, IRTYPE
C...
      NSCOUGE=0
      IGEP=0

      IF((ABS(NOUTGE).EQ.1).OR.(ABS(NOUTGE).EQ.4)) THEN
         CALL OPEN_GBL_FILE(63, TRIM(GLOBALDIR)//'/'//'fort.63',
     $    NP_G, NP, HEADER63)

C   kmd48.33 add the ability to output the sponge layer
          IF (OUTPUTSPONGE) THEN
           CALL writeDomainHeader(92,
     &          TRIM(GLOBALDIR)//'/'//'fort.92',
     $           NP_G, NP, 'sponge     ')
          END IF
C
C        jgf47.06 Tau0 output is produced on the same schedule as
C        global elevation
         IF (OutputTau0) THEN
Casey 090302: Corrected the unit number from 10 to 90.
            CALL writeDomainHeader(90,
     &           TRIM(GLOBALDIR)//'/'//'fort.tau0',
     $           NP_G, NP, 'Tau0      ')
         ENDIF
         IGEP=2
      ENDIF

      IF(ABS(NOUTGE).EQ.2) THEN
         OPEN(63,FILE=TRIM(LOCALDIR)//'/'//'fort.63',
     &           ACCESS='DIRECT',RECL=NBYTE)
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
         CLOSE(63)
         ENDIF
C...
C....INITIALIZE GLOBAL VELOCITY SPOOL COUNTER
C....OPEN GLOBAL VELOCITY OUTPUT FILE
C....WRITE OUT HEADER INFORMATION INCLUDING NDSETSV
C....(NO. OF GLOBAL VELOCITY DATA SETS TO BE SPOOLED),
C....NP, DT*NSPOOLGV, NSPOOLGV, IRTYPE
C...
      NSCOUGV=0
      IGVP=0

!kmd - changed for -4 option      IF((ABS(NOUTGV).EQ.1)) THEN
!Casey 090828: Added the possibility of NOUTGV = +/- 4.
      IF((ABS(NOUTGV).EQ.1).OR.(ABS(NOUTGV).EQ.4)) THEN
         CALL OPEN_GBL_FILE(64, TRIM(GLOBALDIR)//'/'//'fort.64',
     $    NP_G, NP, HEADER64)
         IGVP=2
         ENDIF

      IF(ABS(NOUTGV).EQ.2) THEN
         OPEN(64,FILE=TRIM(LOCALDIR)//'/'//'fort.64',
     &        ACCESS='DIRECT',RECL=NBYTE)
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
         CLOSE(64)
         ENDIF
C...
C....INITIALIZE GLOBAL WIND and pressure SPOOL COUNTER
C....OPEN GLOBAL WIND and pressure OUTPUT FILEs
C....WRITE OUT HEADER INFORMATION INCLUDING NDSETSW
C....(NO. OF GLOBAL WIND DATA SETS TO BE SPOOLED),
C....NP, DT*NSPOOLGW, NSPOOLGW, IRTYPE
C...
      NSCOUGW=0
      IGWP=0
      IGPP=0
      IGIP=0  !tcm v49.64.01 added for global ice concentration field

!Casey 090302: Added the possibility of NOUTGW = +/- 4.
!kmd - changed for -4 option      IF((ABS(NOUTGW).EQ.1)) THEN
      IF((ABS(NOUTGW).EQ.1).OR.(ABS(NOUTGW).EQ.4)) THEN
         CALL OPEN_GBL_FILE(73, TRIM(GLOBALDIR)//'/'//'fort.73',
     $    NP_G, NP, HEADER73)
         IGPP=2
         CALL OPEN_GBL_FILE(74, TRIM(GLOBALDIR)//'/'//'fort.74',
     $     NP_G, NP, HEADER74)
         IGWP=2
! tcm v49.64.01 added for ice
         IF(NCICE.NE.0) then
            CALL OPEN_GBL_FILE(93, TRIM(GLOBALDIR)//'/'//'fort.93',
     $        NP_G, NP, HEADER93)
            IGIP=2
         ENDIF
      ENDIF

      IF(ABS(NOUTGW).EQ.2) THEN
         OPEN(73,file=trim(LOCALDIR)//'/'//'fort.73',
     &           ACCESS='DIRECT',RECL=NByte)
         OPEN(74,FILE=TRIM(LOCALDIR)//'/'//'fort.74',
     &           ACCESS='DIRECT',RECL=NByte)
! tcm v49.64.01 added for ice
         IF(NCICE.NE.0) THEN
            OPEN(93,FILE=TRIM(LOCALDIR)//'/'//'fort.93',
     &           ACCESS='DIRECT',RECL=NByte)
         ENDIF
         IF(NBYTE.EQ.4) THEN
            DO I=1,8
               write(73,rec=igpp+i) rdes4(i)
               WRITE(74,REC=IGWP+I) RDES4(I)
               if(ncice.ne.0) write(93,rec=igip+i) rdes4(i)
            ENDDO
            igpp=igpp+8
            IGWP=IGWP+8
            if(ncice.ne.0) igip = igip+8
            DO I=1,6
               write(73,rec=igpp+i) rid4(i)
               WRITE(74,REC=IGWP+I) RID4(I)
               if(ncice.ne.0) write(93,rec=igip+i) rid4(i)
            ENDDO
            igpp=igpp+6
            IGWP=IGWP+6
            if(ncice.ne.0) igip = igip+6
            DO I=1,6
               write(73,rec=igpp+i) aid4(i)
               WRITE(74,REC=IGWP+I) AID4(I)
               if(ncice.ne.0) write(93,rec=igip+i) aid4(i)
            ENDDO
            igpp=igpp+6
            IGWP=IGWP+6
            if(ncice.ne.0) igip = igip+6
         ENDIF
         IF(NBYTE.EQ.8) THEN
            DO I=1,4
               write(73,rec=igpp+i) rdes8(i)
               WRITE(74,REC=IGWP+I) RDES8(I)
               if(ncice.ne.0) write(93,rec=igip+i) rdes8(i)
            ENDDO
            igpp=igpp+4
            IGWP=IGWP+4
            if(ncice.ne.0) igip = igip+4
            DO I=1,3
               write(73,rec=igpp+i) rid8(i)
               WRITE(74,REC=IGWP+I) RID8(I)
               if(ncice.ne.0) write(93,rec=igip+i) rid8(i)
            ENDDO
            igpp=igpp+3
            IGWP=IGWP+3
            if(ncice.ne.0) igip = igip+3
            DO I=1,3
               write(73,rec=igpp+i) aid8(i)
               WRITE(74,REC=IGWP+I) AID8(I)
               if(ncice.ne.0) write(93,rec=igip+i) aid8(i)
            ENDDO
            igpp=igpp+3
            IGWP=IGWP+3
            if(ncice.ne.0) igip = igip+3
         ENDIF
         write(73,rec=igpp+1) ndsetsw
         write(73,rec=igpp+2) np
         write(73,rec=igpp+3) dt*nspoolgw
         write(73,rec=igpp+4) nspoolgw
         write(73,rec=igpp+5) 2
         igpp=igpp+5
         close(73)
         WRITE(74,REC=IGWP+1) NDSETSW
         WRITE(74,REC=IGWP+2) NP
         WRITE(74,REC=IGWP+3) DT*NSPOOLGW
         WRITE(74,REC=IGWP+4) NSPOOLGW
         WRITE(74,REC=IGWP+5) 2
         IGWP=IGWP+5
         CLOSE(74)
         IF(NCICE.NE.0) THEN
            WRITE(93,REC=IGIP+1) NDSETSW
            WRITE(93,REC=IGIP+2) NP
            WRITE(93,REC=IGIP+3) DT*NSPOOLGW
            WRITE(93,REC=IGIP+4) NSPOOLGW
            WRITE(93,REC=IGIP+5) 2
            IGIP=IGIP+5
            CLOSE(93)
         ENDIF
      ENDIF

!  tcm v50.75 changed ifdef cwan to output for nrs=3 or 4
      IF ((ABS(NRS).EQ.3).or.(ABS(NRS).eq.4)) then
!#ifdef CSWAN  tcm temp
Casey 090302: Added the output of radiation stress gradients.
         IGRadS=0
         IF((ABS(NOUTGW).EQ.1).OR.(ABS(NOUTGW).EQ.4)) THEN
            CALL OPEN_GBL_FILE(164, TRIM(GLOBALDIR)//'/'//'rads.64',
     &          NP_G, NP, HEADER74)
            IGRadS=2
         ENDIF
         IF(ABS(NOUTGW).EQ.2) THEN
            OPEN(164,FILE=TRIM(LOCALDIR)//'/'//'rads.64',
     &              ACCESS='DIRECT',RECL=NByte)
            IF(NBYTE.EQ.4) THEN
               DO I=1,8
                  WRITE(164,REC=IGRadS+I) RDES4(I)
                  ENDDO
               IGRadS=IGRadS+8
               DO I=1,6
                  WRITE(164,REC=IGRadS+I) RID4(I)
                  ENDDO
               IGRadS=IGRadS+6
               DO I=1,6
                  WRITE(164,REC=IGRadS+I) AID4(I)
               ENDDO
               IGRadS=IGRadS+6
            ENDIF
            IF(NBYTE.EQ.8) THEN
               DO I=1,4
                  WRITE(164,REC=IGRadS+I) RDES8(I)
               ENDDO
               IGRadS=IGRadS+4
               DO I=1,3
                  WRITE(164,REC=IGRadS+I) RID8(I)
               ENDDO
               IGRadS=IGRadS+3
               DO I=1,3
                  WRITE(164,REC=IGRadS+I) AID8(I)
               ENDDO
               IGRadS=IGRadS+3
            ENDIF
            WRITE(164,REC=IGRadS+1) NDSETSW
            WRITE(164,REC=IGRadS+2) NP
            WRITE(164,REC=IGRadS+3) DT*NSPOOLGW
            WRITE(164,REC=IGRadS+4) NSPOOLGW
            WRITE(164,REC=IGRadS+5) 2
            IGRadS=IGRadS+5
            CLOSE(164)
         ENDIF
      endif
!#endif

C...
C....INITIALIZE GLOBAL CONCENTRATION SPOOL COUNTER
C....OPEN GLOBAL CONCENTRATION OUTPUT FILE
C....WRITE OUT HEADER INFORMATION INCLUDING NDSETSC
C....(NO. OF GLOBAL CONCENTRATION DATA SETS TO BE SPOOLED),
C....NP, DT*NSPOOLGC, NSPOOLGC, IRTYPE
C...
      NSCOUGC=0
      IGCP=0

      IF(ABS(NOUTGC).EQ.1) THEN
         OPEN(83,FILE=TRIM(INPUTDIR)//'/'//'fort.83')
         WRITE(83,3220) RUNDES,RUNID,AGRID
         WRITE(83,3645) NDSETSC,NP,DTDP*NSPOOLGC,NSPOOLGC,1
         CLOSE(83)
         IGCP=2
         ENDIF

      IF(ABS(NOUTGC).EQ.2) THEN
         OPEN(83,FILE=TRIM(INPUTDIR)//'/'//'fort.83',
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
         CLOSE(83)
      ENDIF
C
 1112 FORMAT(/,1X,79('_'))
 9883 FORMAT(' !!! NODE ',I6,' DRIED (LANDLOCKING)')


C... TCM V50.66.01 ADDITIONS FOR TIME VARYING BATHYMETRY
C... THESE FILES ARE CONTROLED BY THE GLOBAL ELEVATION
C...
C....INITIALIZE GLOBAL BATHYMETRY SPOOL COUNTER
C....OPEN GLOBAL BATHYMETRY OUTPUT FILE
C....WRITE OUT HEADER INFORMATION INCLUDING NDSETSE
C....(NO. OF GLOBAL BATHYMETRY DATA SETS TO BE SPOOLED),
C....NP, DT*NSPOOLGE, NSPOOLGE, IRTYPE
C...
      NSCOUGB=0
      IGBP=0
      
      IF ((NDDT.NE.0).AND.
     &       ((ABS(NOUTGE).EQ.1).OR.(ABS(NOUTGE).EQ.4))) THEN
         CALL OPEN_GBL_FILE(76, TRIM(GLOBALDIR)//'/'//'fort.76',
     &    NP_G, NP, HEADER63)  !NO NEED TO CREATE HEADER76
         IGBP=2
      ENDIF

      IF ((NDDT.NE.0).AND.(ABS(NOUTGE).EQ.2)) THEN
         OPEN(76,FILE=TRIM(LOCALDIR)//'/'//'fort.76',
     &           ACCESS='DIRECT',RECL=NBYTE)
         IF(NBYTE.EQ.4) THEN
            DO I=1,8
               WRITE(76,REC=IGBP+I) RDES4(I)
               ENDDO
            IGBP=IGBP+8
            DO I=1,6
               WRITE(76,REC=IGBP+I) RID4(I)
               ENDDO
            IGBP=IGBP+6
            DO I=1,6
               WRITE(76,REC=IGBP+I) AID4(I)
               ENDDO
            IGBP=IGBP+6
         ENDIF
         IF(NBYTE.EQ.8) THEN
            DO I=1,4
               WRITE(76,REC=IGBP+I) RDES8(I)
               ENDDO
            IGBP=IGBP+4
            DO I=1,3
               WRITE(76,REC=IGBP+I) RID8(I)
               ENDDO
            IGBP=IGBP+3
            DO I=1,3
               WRITE(76,REC=IGBP+I) AID8(I)
               ENDDO
            IGBP=IGBP+3
         ENDIF
         WRITE(76,REC=IGBP+1) NDSETSE
         WRITE(76,REC=IGBP+2) NP
         WRITE(76,REC=IGBP+3) DT*NSPOOLGE
         WRITE(76,REC=IGBP+4) NSPOOLGE
         WRITE(76,REC=IGBP+5) 1
         IGBP=IGBP+5
         CLOSE(76)
      ENDIF


C...
C......INITIALIZE 3D SOLUTION
C...

C...  LINES TO RUN THE CODE IN 3D VS MODE.

      if (C3D) then
         CALL COLDSTART_3D()
      endif

C...LINES TO RUN THE CODE IN 3D DSS MODE

c     if (C3DDSS) then
c       CALL DSSSTUP(DT,NT)
c     endif
C...
      ! tcmv51.27 20140502 Added to remove this nodal attribute
      ! as it is only used at Cold Start
      IF (LoadRiver_et_WSE) THEN
         IF (ALLOCATED (River_et_WSE) ) THEN
            DEALLOCATE (River_et_WSE)
         ENDIF
      ENDIF

#if defined(COLDSTART_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN
      END SUBROUTINE COLDSTART


C******************************************************************************
C   Subroutine to initialize the 3D routines for a cold start including       *
C   reading in an initial density field.                                      *
C                                                                             *
C                                                                             *
c******************************************************************************

      SUBROUTINE COLDSTART_3D()
      USE SIZES, ONLY : WRITE_LOCAL_FILES, MNPROC, GLOBALDIR
      USE WRITE_OUTPUT, ONLY : initOutput3D
      USE GLOBAL_3DVS
      USE GLOBAL, ONLY : BCFLAG_TEMP, RES_BC_FLAG, BCFLAG_LNM,
     &                   RBCTIME1, BCSTATIM, RBCTIME2,
     &                   RBCTIMEINC, SBCTIME1,
     &                   SBCTIME2, SBCSTATIM, TTBCTIME1,
     &                   TTBCTIME2, TBCSTATIM, SBCTIMEINC,
     &                   TBCTIMEINC, TTBCTIMEINC, TTBCSTATIM,
     &                   TBCTIME1, TBCTIME2, q_heat1, q_heat2, 
     &                   LNM_BC, LNM_BC1, LNM_BC2, RIVBCTIMINC,
     &                   RIVBCTIME1, RIVBCTIME2, 
     &                   RIVBCSTATIM, DEBUG, setMessageSource,
     &                   unsetMessageSource, allMessage, NP_G
      USE MESH, ONLY : NP, AGRID
      USE BOUNDARIES, ONLY : NETA, NOPE, NBD,
     &  totalbcrivernodes, bndbcriver
      USE ADCIRC_MOD, ONLY : ADCIRC_Terminate
      IMPLICIT NONE

      INTEGER :: IRType
      INTEGER :: N                       !loop counter
      INTEGER :: NH, IHN, IVN            !horizontal & vertical loop counters
      INTEGER :: NHNN, NVNN              !horizontal and vertical node numbers
      INTEGER :: NVN
      INTEGER NVP                        ! number of horizontal nodes

ckmd  additional variable for baroclinic boundary conditions
      CHARACTER(80) :: CDUM80
      INTEGER :: NumofBCNodes
      INTEGER :: NOD, J, NumofNodes, K
      INTEGER :: IINDX
      REAL(SZ), ALLOCATABLE :: TMP(:,:)
C
      call setMessageSource("coldstart_3D")
#if defined(COLDSTART_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C...
C...  Define format statements used to initialize 3D output files
C...
 499  FORMAT(1X,A32,2X,A24,2X,A24)

C...  IF A BAROCLINIC RUN, READ IN INITIAL DENSITY FIELD

      IF(CBaroclinic.and.(.not.C3DVS)) THEN
         IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
            WRITE(16,*) "CBaroclinic is ",CBaroclinic
            WRITE(16,*) "CB3DVS is ",C3DVS
            WRITE(16,*)
     &           "ERROR: 2DDI baroclinic runs not currently supported."
         ENDIF
         CALL ADCIRC_Terminate()
      ENDIF
C

      IF(CBaroclinic .and. C3DVS) THEN
         WRITE(16,424)
 424     FORMAT(/,5X,'INITIAL DENSITY FIELD READ IN FROM UNIT 11',/)
         OPEN(11,FILE=TRIM(INPUTDIR)//'/'//'fort.11')
         READ(11,*)             !skip over header line
         READ(11,*)             !skip over header line
         READ(11,*) NVN, NVP
         IF(NVN.NE.NFEN) THEN
            IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
               WRITE(ScreenUnit,351) NVN,NFEN
            ENDIF
            WRITE(16,351) NVN,NFEN
 351        FORMAT(/,2X,'***** INVALID INPUT IN THE DENSITY INITIAL ',
     &           'CONDITION FILE (UNIT 11) *****',
     &           /,2X,'***** NVN = ',I4,' MUST MATCH NFEN = ',I4,
     &           ' *****',
     &           /,10X,'****** RUN TERMINATED ******')
            CALL ADCIRC_Terminate()
         ENDIF
         DO IHN=1,NP
            DO IVN=1,NFEN
               IF(ABS(IDen).EQ.1)
     &              READ(11,*) NHNN,NVNN,SigT(NHNN,NVNN)
               IF(ABS(IDen).EQ.2)
     &              READ(11,*) NHNN,NVNN,Sal(NHNN,NVNN)
               IF(ABS(IDen).EQ.3)
     &              READ(11,*) NHNN,NVNN,Temp(NHNN,NVNN)
               IF(ABS(IDen).EQ.4)
     &              READ(11,*) NHNN,NVNN,Temp(NHNN,NVNN),Sal(NHNN,NVNN)
            ENDDO
         ENDDO
         CLOSE(11)
      ENDIF

      IF (CBaroclinic) THEN
!     Kendra45.12: Read in the temperature boundary condition for the
!     temperature field. Note, we will add the salinity condition
!     later. Currently, we are adopting the sign convention that the
!     flux into the domain is consider positive and the flux out of the
!     domain is negative. However, remember the constitutive law changes
!     the direction of the flux.
         IF ((IDEN.EQ.3).OR.(IDEN.EQ.-3)) THEN
            IF (NTF.EQ.0) THEN
               DO IHN=1,NP
                  qsurfkp1(IHN)=0.d0
                  qsurf(IHN)=0.d0
               END DO
            ELSE IF (NTF.EQ.1) THEN
               OPEN(111,FILE=TRIM(INPUTDIR)//'/'//'fort.111')
               DO IHN=1,NP
                  READ(111,*) NHNN,qsurfkp1(IHN),qsurf(IHN)
               END DO
               CLOSE(111)
            END IF
         END IF
      END IF

ckmd  Read in boundary condition information from files if
C...  we are using the transport portion of the code (IDEN.GE.1)
C...  and baroclinic.

      IF ((CBAROCLINIC).AND.(IDEN.GE.1)) THEN
         IF ((ABS(RES_BC_FLAG).GE.1).AND.(NOPE.GT.0)) THEN
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
           IF (NOPE.GT.0) THEN
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
           END IF
         END IF
         IF ((ABS(RES_BC_FLAG).EQ.3).OR.
     &            (ABS(RES_BC_FLAG).EQ.4)) THEN
           IF (NOPE.GT.0) THEN
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
           END IF
           IF (BCFLAG_TEMP.EQ.1) THEN ! read in file with multiple values but no calcs
              OPEN(38,FILE=TRIM(INPUTDIR)//'/'//'fort.38')
              TTBCTIME1=TTBCSTATIM*86400.D0
              TTBCTIME2=TTBCTIME1+TTBCTIMEINC
              DO NumofNodes=1,NP
                 READ(38,*) NOD, q_heat1(Numofnodes)
              END DO
              DO NumofNodes=1,NP
                 READ(38,*) NOD, q_heat2(NumofNodes)
              END DO
           ELSE IF (BCFLAG_TEMP.EQ.2) THEN ! read in from file with 6 components
              ALLOCATE (TMP(NP,6))
              OPEN(38,FILE=TRIM(INPUTDIR)//'/'//'fort.38')
              TTBCTIME1=TTBCSTATIM*86400.D0
              TTBCTIME2=TTBCTIME1+TTBCTIMEINC
              READ(38,*) (NOD,(TMP(K,J),J=1,6),K=1,NP)
              DO NumofNodes=1,NP
                 q_heat1(NumofNodes)=-TMP(NumofNodes,1)-
     &              TMP(NumofNodes,2)+TMP(NumofNodes,3)-
     &              TMP(NumofNodes,5)+TMP(NumofNodes,4)-
     &              TMP(NumofNodes,6)
              END DO
              READ(38,*) (NOD,(TMP(K,J),J=1,6),K=1,NP)
              DO NumofNodes=1,NP
                 q_heat2(NumofNodes)=-TMP(NumofNodes,1)-
     &              TMP(NumofNodes,2)+TMP(NumofNodes,3)-
     &              TMP(NumofNodes,5)+TMP(NumofNodes,4)-
     &              TMP(NumofNodes,6)
              END DO
              DEALLOCATE(TMP)
           ELSE IF (BCFLAG_TEMP.EQ.3) THEN ! read in from file with 4 components
              ALLOCATE (TMP(NP,4))
              OPEN(38,FILE=TRIM(INPUTDIR)//'/'//'fort.38')
              TTBCTIME1=TTBCSTATIM*86400.D0
              TTBCTIME2=TTBCTIME1+TTBCTIMEINC
              READ(38,*) (NOD,(TMP(K,J),J=1,4),K=1,NP)
              DO NumofNodes=1,NP
                 q_heat1(NumofNodes)=TMP(NumofNodes,4)+
     &              TMP(NumofNodes,3)-TMP(NumofNodes,1)+
     &              TMP(NumofNodes,2)
              END DO
              READ(38,*) (NOD,(TMP(K,J),J=1,4),K=1,NP)
              DO NumofNodes=1,NP
                 q_heat2(NumofNodes)=TMP(NumofNodes,4)+
     &              TMP(NumofNodes,3)-TMP(NumofNodes,1)+
     &              TMP(NumofNodes,2)
              END DO
              DEALLOCATE(TMP)
           END IF ! end if for heat flux file
         END IF ! end if for temperature boundary and heat flux
      ELSE IF ((CBAROCLINIC).AND.(IDEN.LE.-1)) THEN
         IF ((RES_BC_FLAG.LT.0).AND.(NOPE.GT.0)) THEN
          IF (BCFLAG_LNM.EQ.1) THEN
            OPEN(35,FILE=TRIM(INPUTDIR)//'/'//'fort.35')
            READ(35,'(A)') CDUM80
            DO NumofBCNodes=1,NETA
                READ(35,*) NOD,LNM_BC(NumofBCNodes)
            END DO
          ELSE IF (BCFLAG_LNM.EQ.2) THEN
            ! add later
          ELSE IF (BCFLAG_LNM.EQ.3) THEN
            DO NumofBCNodes=1,NETA
               IINDX=NBD(NumofBCNodes)
               LNM_BC(NumofBCNodes)=ETA2(IINDX)
            END DO
          END IF
         END IF
      END IF

C kmd - River information for boundary conditions in baroclinic simulation
      IF (BndBCRiver) THEN
        OPEN(39,FILE=TRIM(INPUTDIR)//'/'//'fort.39')
        READ(39,*) RIVBCTIMINC, RIVBCSTATIM
        RIVBCTIME1=RIVBCSTATIM*86400.d0
        RIVBCTIME2=RIVBCTIME1+RIVBCTIMINC
        DO J=1,totalbcrivernodes
         IF (IDEN.EQ.2) THEN
           READ(39,*) (BCRivSalN1(J,K),K=1,NFEN)
         ELSE IF (IDEN.EQ.3) THEN
           READ(39,*) (BCRivTempN1(J,K),K=1,NFEN)
         ELSE IF (IDEN.EQ.4) THEN
           READ(39,*) (BCRivSalN1(J,K), BCRivTempN1(J,K),K=1,NFEN)
         END IF
        END DO
        DO J=1,totalbcrivernodes
         IF (IDEN.EQ.2) THEN
           READ(39,*) (BCRivSalN2(J,K),K=1,NFEN)
         ELSE IF (IDEN.EQ.3) THEN
           READ(39,*) (BCRivTempN2(J,K),K=1,NFEN)
         ELSE IF (IDEN.EQ.4) THEN
           READ(39,*) (BCRivSalN2(J,K), BCRivTempN2(J,K),K=1,NFEN)
         END IF
        END DO
      END IF

C...  ZERO OUT STUFF PASSED FROM 3D SOLUTION TO EXTERNAL MODE

      DO NH=1,NP
         DUU(NH)=0.d0
         DUV(NH)=0.d0
         DVV(NH)=0.d0
         UU(NH)=0.d0
         VV(NH)=0.d0
         DAFluxX(NH)=0.d0
         DAFluxY(NH)=0.d0
         BSX(NH)=0.d0
         BSY(NH)=0.d0
      ENDDO

C...  INITIALIZE 3D VELOCITY AND TURBULENCE SOLUTION

!     kmd45.12 Need to initialize the new qkp1 and wzkp1 arrays

      DO NH=1,NP
         DO N=1,NFEN
            Q(NH,N)=(0.d0,0.d0)
            Qkp1(NH,N)=(0.d0,0.d0)
C. same problem as with L
            q20(NH,N)=-9999d0
C.RJW bug fix cannot initialize L to zero, must initialize to lmin or non-zero
C.value, so when elements wet, the first turb length scale is nonzero
            l(NH,N)=-9999d0
            wz(NH,N)=0.d0
            wzkp1(NH,N)=0.d0
         ENDDO
      ENDDO
      bpg(:,:)=(0.d0,0.d0) ! jgf50.44: was not initialized anywhere else

C...  INITIALIZE 3D OUTPUT FILES
C
C.... Initialize the 3D density station output file (Unit 41)
      IF((IDen.EQ.1).OR.(IDen.EQ.0)) IRType=1
      IF((IDen.EQ.2).OR.(IDen.EQ.3)) IRType=2
      IF(IDen.EQ.4) IRType=3
      CALL initOutput3D(41, I3DSD, NDSET3DSD, NSTA3DD, NSTA3DD_G,
     &                  NSPO3DSD, IRType, I3DSDRec)
C
C.... Initialize the 3D velocity station output file (Unit 42)
      IRType=3
      CALL initOutput3D(42, I3DSV, NDSET3DSV, NSTA3DV, NSTA3DV_G,
     &                  NSPO3DSV, IRType, I3DSVRec)
C
C.... Initialize the 3D turbulence station output file (Unit 43)
      IRType=3
      CALL initOutput3D(43, I3DST, NDSET3DST, NSTA3DT, NSTA3DT_G,
     &                  NSPO3DST, IRType, I3DSTRec)
C.... Initialize the global 3D density output file (Unit 44)

      IF((IDen.EQ.1).OR.(IDen.EQ.0)) IRType=1
      IF((IDen.EQ.2).OR.(IDen.EQ.3)) IRType=2
      IF(IDen.EQ.4) IRType=3
C     jgf46.27 Replaced IRType with IDen in header.
!kmd48.33bc - add in infomration for file containing the
!            top temperature boundary condition information
      CALL initOutput3D(44, I3DGD, NDSET3DGD, NP, NP_G,
     &                  NSPO3DGD, IRType, I3DGDRec)
      IF (ABS(I3DGD).EQ.1) THEN
         IF ((MNPROC.EQ.1.).OR.(WRITE_LOCAL_FILES.eqv..true.)) THEN
            IF((IDEN.EQ.3).OR.(IDEN.EQ.4)) THEN
               IF (BCFLAG_TEMP.NE.0) THEN
                  OPEN(47,FILE=TRIM(LOCALDIR)//'/'//'fort.47')
                  WRITE(47,499) RUNDES,RUNID,AGRID
                  WRITE(47,498) NDSet3DGD,NP,DTDP*NSpo3DGD,NSpo3DGD,1
                  CLOSE(47)
               END IF
            END IF
         ENDIF
         IF ((MNPROC.GT.1.).AND.(WRITE_LOCAL_FILES.eqv..false.).and.
     &       (myProc.eq.0)) THEN
            IF((IDEN.EQ.3).OR.(IDEN.EQ.4)) THEN
               IF (BCFLAG_TEMP.NE.0) THEN
                  OPEN(47,FILE=TRIM(GLOBALDIR)//'/'//'fort.47')
                  WRITE(47,499) RUNDES,RUNID,AGRID
                  WRITE(47,498) NDSet3DGD,NP,DTDP*NSpo3DGD,NSpo3DGD,1
                  CLOSE(47)
               END IF
            END IF
         ENDIF
      ENDIF
C..RJW bug fix in 498 (kendra found this)
 498  FORMAT(1X,I10,1X,I10,1X,E15.7,I10,1X,I10,1X,I3)

C
c.... Initialize the global 3D velocity output file (Unit 45)
      CALL initOutput3D(45, I3DGV, NDSET3DGV, NP, NP_G,
     &                  NSPO3DGV, IRType, I3DGVRec)
C
c.... Initialize the global 3D turbulence output file (Unit 46)
      CALL initOutput3D(46, I3DGT, NDSET3DGT, NP, NP_G,
     &                  NSPO3DGT, IRType, I3DGTRec)
C
C.... Set up a few final odds and ends for a 3D run
      CALL VSSTUP ()

#if defined(COLDSTART_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE COLDSTART_3D
C----------------------------------------------------------------------

