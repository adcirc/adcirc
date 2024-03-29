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
      USE GLOBAL, only: DP1, DP2, eta1, qx2, qy2, tk, etas, eta2, ch1,
     &    en2, nddt, btime1, statim, btime2, btiminc, btime_end, qn1,
     &    qn0, nffr, ftiminc, qnin1, enin1, qnin2, enin2, qtime1, nws,
     &    qtime2, noute, noutv, noutc, noutge, noutgv, noutgc, nrs, qn2,
     &    noute, noutm, noutgw, metonly, rstime1, rstime2, rstiminc,
     &    rsnx1, rsny1, rsnx2, rsny2, ncice, cice1, cice2, cice_time1,
     &    cice_time2, cice_timinc, ctip, tip2, nscoue, iestp, nstae_g,
     &    nstae, rdes4, rid4, rdes8, rid8, aid8, esbin1, vv2, nolifa,
     &    uu2, inundationoutput, esbin2, nspoolv, nscouc, icstp, nstac,
     &    nstac_g, ntrspc, nspoolc, nscoum, iwstp, bchgtiminc, nscreen,
     &    screenunit, h0, ntrspe, dt, nspoole, nscouv, ivstp, np_g,
     &    outputsponge, ne_g, usingDynamicWaterLevelCorrection, ndsetse,
     &    nspoolge,
     &    nscougv, igvp, ndsetsv, nspoolgv, use_tvw, nscou_tvw, ig_tvw,
     &    nout_tvw, ndsets_tvw, nspool_tvw, nscougw, igwp, igpp, igip,
     &    igep, ndsetsw, nspoolgw, igrads, nscougc, igcp, c2ddi, iden, 
     &    rundes, runid, dtdp, ntrspm, nspoolm, nscoub, ibstp, nscouge, 
     &    ndsetsc, nspoolgc, nscougb, igbp, c3d, ipstp, nstav_g, nstav,
     &    ntrspv, iicestp, nstam_g, nstam, openfileforread, cbaroclinic
      USE GLOBAL_3DVS, ONLY: I3DSD, I3DSV, I3DST, I3DGD, I3DGV, I3DGT
      USE MESH, ONLY : NE, NP, X, Y, SLAM, SFEA, ICS, DP, NM, MJU, AID4,
     &    agrid 
      USE BOUNDARIES, ONLY : NETA, NFLUXF, NOPE, NVEL, LBCODEI, NPEBC
      USE HARM
      USE WETDRY, ONLY : initializeWettingAndDrying
      USE WIND, ONLY : coldStartMeteorologicalForcing, rsget
      USE GLOBAL_IO, only: header63, header64, header77, header73,
     &    header74, header93, header91, header61, header62, header71,
     &    header72, header81, writedomainheader, open_gbl_file, nddt2get
      USE SUBDOMAIN, ONLY : subdomainOn, readFort015
      USE ADCIRC_MOD, ONLY : ADCIRC_Terminate
      USE OWI_ICE, ONLY : NCICE1_INIT, NCICE1_GET  !tcm v49.64.01
      USE RS2,ONLY : RS2INIT,RS2GET            ! sb added 09/xx/2006
      USE NODALATTRIBUTES, ONLY : outputTau0
      USE WRITE_OUTPUT, only : outputNodeCode, outputNOFF,
     &    initInundationOutput
      USE TIME_VARYING_WEIR_BOUNDARY,ONLY:WEIR_SETUP
      USE GWCE, ONLY : etime1, etime2, etiminc
#ifdef ADCNETCDF
      use Couple2BC3D, only: Initial_BC3D_NetCDF
#endif
      IMPLICIT NONE

      INTEGER I, J,IOS, IJK
      INTEGER NM1, NM2, NM3
      INTEGER NC1, NC2, NC3, NCELE

      REAL(8) HTOT
!kmd - Evan's updates for rivers above MSL
!      REAL(8),ALLOCATABLE :: et_WSE(:)
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
      !
      call initializeWettingAndDrying()

C SET UP TIME VARYING WEIRS IF NECESSARY
      CALL WEIR_SETUP()

      if (subdomainOn) then
         call readFort015() ! NCSU Subdomain Modeling
      endif
      
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
C...  SET AT REST INITIAL CONDITION OVER WHOLE DOMAIN
C...  IF BOTTOM IS ABOVE THE GEIOD -> DRY NODE
C...  IF BOTTOM IS INITIALLY BELOW THE GEOID AND STARTDRY=1 -> DRY NODE
C...
      UU2(:) = 0.D0
      VV2(:) = 0.D0
      QX2(:) = 0.D0
      QY2(:) = 0.D0
      TK(:) = 0.D0     
      ETA1(:) = 0.D0
      ETA2(:) = 0.D0
      ETAS(:) = 0.D0
      CH1(:) = 0.d0
      EN2(:) = 0.d0
      !
      ! Set the initial wet/dry state and water surface elevation
      call setColdWetDryStateAndWaterSurfaceElevation()  
      !
      ! jgf52.08.01: Initialize inundation arrays if necessary. This
      ! can only be called after nnodecode has been set via 
      ! setColdWetDryState() in coldstart() or hstart().
      if (inundationOutput.eqv..true.) then
         call initInundationOutput()
      endif
C...
C...  INITIALIZE THE ELEVATION SPECIFIED BOUNDARY CONDITION IF IT
C...  REQUIRES THE USE OF THE UNIT 19 FILE.
C...
      IF (NPEBC) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
         WRITE(16,1112)
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1977)
         WRITE(16,1977)
 1977    FORMAT(/,1X,'ELEVATION SPECIFIED INFORMATION READ FROM UNIT ',
     &        '19',/)
         ! NCSU Subdomain Modeling
         if (subdomainOn.eqv..false.) then
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
         endif
      ENDIF
C
C....INITIALIZE THE NORMAL FLOW BOUNDARY CONDITION
C
      DO I=1,NVEL
         QN2(I)=0.D0
         QN1(I)=0.D0
         QN0(I)=0.D0
      END DO

      IF(NFLUXF.EQ.1.AND.NFFR.LE.0)THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
         WRITE(16,1112)
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1979)
         WRITE(16,1979)
 1979    FORMAT(/,1X,'NORMAL FLOW INFORMATION READ FROM UNIT 20',/)
         OPEN(20,FILE=TRIM(INPUTDIR)//'/'//'fort.20')
         READ(20,*) FTIMINC

         DO J=1,NVEL
            QNIN1(J)=0.D0
            IF ((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &          .OR.(LBCODEI(J).EQ.22)) THEN
               READ(20,*) QNIN1(J)
            ELSEIF (LBCODEI(J).EQ.32) THEN
               READ(20,*) QNIN1(J),ENIN1(J)
            ENDIF
         END DO
         DO J=1,NVEL
            QNIN2(J)=0.D0
            IF ((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &           .OR.(LBCODEI(J).EQ.22)) THEN
               READ(20,*) QNIN2(J)
            ELSEIF (LBCODEI(J).EQ.32) THEN
               READ(20,*) QNIN2(J),ENIN2(J)
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
     &      'INFO: Only meteorological output was requested. ADCIRC'
               WRITE(ScreenUnit,*)
     &      'will not solve the GWCE or momentum equations.'
             ENDIF
             WRITE(16,*)
     &      'INFO: Only meteorological output was requested. ADCIRC'
             WRITE(16,*)
     &      'will not solve the GWCE or momentum equations.'
            METONLY = .TRUE. ! set flag
         ENDIF

         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
         WRITE(16,1112)
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1980)
         WRITE(16,1980)
 1980    FORMAT(/,1X,'WIND (AND PRESSURE) INFORMATION READ.',/)

         call coldstartMeteorologicalForcing(nws)

      ENDIF

C...INPUT RADIATION STRESS INFORMATION FROM UNIT 23
C....READ IN THE TIME LEVEL 1 AND LEVEL 2 FIELDS

      ! NRS=2 was added. 09/xx/2006 sb
      IF(NRS.GE.1) THEN ! sb46.28sb03
         RSTIME1 = STATIM*86400.D0
         RSTIME2 = RSTIME1+RSTIMINC
         IF(NRS.EQ.1) THEN
           OPEN(23,FILE=TRIM(INPUTDIR)//'/'//'fort.23')
           CALL RSGET(RSNX1,RSNY1)
           CALL RSGET(RSNX2,RSNY2)
         ENDIF
         IF(NRS.EQ.2) THEN
           CALL RS2INIT(RSNX1,RSNY1,NP)
           CALL RS2GET(RSNX1,RSNY1,NP)
           CALL RS2GET(RSNX2,RSNY2,NP)
         ENDIF
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

#ifdef ADCNETCDF     
      ! For the Coupling to 3D Baroclinic mode 
      IF (C2DDI.and.CBaroclinic.and.abs(IDEN).ge.5) THEN
         call Initial_BC3D_NetCDF(STATIM*86400.d0)
      endif
#endif       

C...
C...LINES TO USE TIDAL POTENTIAL FORCING
C...
      if (CTIP) then
         DO I=1,NP
            TIP2(I)=0.0
         END DO
      endif

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
         if (write_local_files.eqv..false.) then
            CALL OPEN_GBL_FILE(61, TRIM(GLOBALDIR)//'/'//'fort.61',
     &       NSTAE_G, NSTAE, HEADER61)
         else 
            CALL OPEN_GBL_FILE(61, TRIM(LOCALDIR)//'/'//'fort.61',
     &       NSTAE_G, NSTAE, HEADER61)
         endif
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
         if (write_local_files.eqv..false.) then
            CALL OPEN_GBL_FILE(62, TRIM(GLOBALDIR)//'/'//'fort.62',
     &          NSTAV_G, NSTAV, HEADER62)
         else
            CALL OPEN_GBL_FILE(62, TRIM(LOCALDIR)//'/'//'fort.62',
     &          NSTAV_G, NSTAV, HEADER62)         
         endif
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
         if (write_local_files.eqv..false.) then
            CALL OPEN_GBL_FILE(81, TRIM(GLOBALDIR)//'/'//'fort.81',
     &          NSTAC_G, NSTAC, HEADER81)
         else
            CALL OPEN_GBL_FILE(81, TRIM(LOCALDIR)//'/'//'fort.81',
     &          NSTAC_G, NSTAC, HEADER81)
         endif 
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
         if (write_local_files.eqv..false.) then
            CALL OPEN_GBL_FILE(71, TRIM(GLOBALDIR)//'/'//'fort.71',
     &          NSTAM_G, NSTAM, HEADER71)
            IPSTP=2
            CALL OPEN_GBL_FILE(72, TRIM(GLOBALDIR)//'/'//'fort.72',
     &         NSTAM_G, NSTAM, HEADER72)
            IWSTP=2
! V49.64.01 TCM -ADDED ICE STATIONS
            IF(NCICE.NE.0) THEN
               CALL OPEN_GBL_FILE(91, TRIM(GLOBALDIR)//'/'//'fort.91',
     &             NSTAM_G, NSTAM, HEADER91)
                IICESTP=2
            ENDIF
         else
            CALL OPEN_GBL_FILE(71, TRIM(LOCALDIR)//'/'//'fort.71',
     &          NSTAM_G, NSTAM, HEADER71)
            IPSTP=2
            CALL OPEN_GBL_FILE(72, TRIM(LOCALDIR)//'/'//'fort.72',
     &         NSTAM_G, NSTAM, HEADER72)
            IWSTP=2
! V49.64.01 TCM -ADDED ICE STATIONS
            IF(NCICE.NE.0) THEN
               CALL OPEN_GBL_FILE(91, TRIM(LOCALDIR)//'/'//'fort.91',
     &             NSTAM_G, NSTAM, HEADER91)
                IICESTP=2
            ENDIF
         endif 
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
         if (write_local_files.eqv..false.) then
            CALL OPEN_GBL_FILE(75, TRIM(GLOBALDIR)//'/'//'fort.75',
     &         NSTAE_G, NSTAE, HEADER61)  !no need to create a header75
         else
            CALL OPEN_GBL_FILE(75, TRIM(LOCALDIR)//'/'//'fort.75',
     &         NSTAE_G, NSTAE, HEADER61)  !no need to create a header75
         endif       
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
         if (write_local_files.eqv..false.) then
            CALL OPEN_GBL_FILE(63, TRIM(GLOBALDIR)//'/'//'fort.63',
     &       NP_G, NP, HEADER63)
C            kmd48.33 add the ability to output the sponge layer
             IF (OUTPUTSPONGE) THEN
                 CALL writeDomainHeader(92,
     &             TRIM(GLOBALDIR)//'/'//'fort.92',
     &             NP_G, NP, 'sponge     ')
             END IF
C
C           jgf47.06 Tau0 output is produced on the same schedule as
C           global elevation
            IF (OutputTau0) THEN
Casey 090302: Corrected the unit number from 10 to 90.
               CALL writeDomainHeader(90,
     &              TRIM(GLOBALDIR)//'/'//'fort.90',
     &              NP_G, NP, 'Tau0      ')
            ENDIF
         else

            CALL OPEN_GBL_FILE(63, TRIM(LOCALDIR)//'/'//'fort.63',
     &         NP_G, NP, HEADER63)
C              kmd48.33 add the ability to output the sponge layer
            IF (OUTPUTSPONGE) THEN
               CALL writeDomainHeader(92,
     &           TRIM(LOCALDIR)//'/'//'fort.92',
     &           NP_G, NP, 'sponge     ')
            END IF
C
C           jgf47.06 Tau0 output is produced on the same schedule as
C           global elevation
            IF (OutputTau0.eqv..true.) THEN
Casey 090302: Corrected the unit number from 10 to 90.
               CALL writeDomainHeader(90,
     &           TRIM(GLOBALDIR)//'/'//'fort.90',
     $           NP_G, NP, 'Tau0      ')
            ENDIF
         endif
         IGEP=2
      ENDIF

      ! jgf52.04.03: Create files for nodecode and noff output 
      ! (wet/dry state of nodes and wet/dry state of elements)      
      IF (ABS(NOUTGE).NE.0) then
         if (outputNodeCode.eqv..true.) then
            if (write_local_files.eqv..true.) then
               call writeDomainHeader(101,trim(localdir)//'/'//
     &            'nodecode.63',NP_G, NP, 'nodecode    ')
            else
               call writeDomainHeader(101,trim(globaldir)//'/'//
     &            'nodecode.63',NP_G, NP, 'nodecode    ')
            endif
         endif
         if (outputNOFF.eqv..true.) then
            if (write_local_files.eqv..true.) then
               call writeDomainHeader(102,trim(localdir)//'/'//
     &            'noff.100',NE_G, NE, 'noff       ')
            else
               call writeDomainHeader(102,trim(globaldir)//'/'//
     &            'noff.100',NE_G, NE, 'noff       ')
            endif
         endif
      endif

      ! jgf: Added init for full domain dynamic waterlevel correction output file 
      if ( (ABS(NOUTGE).eq.1.or.ABS(NOUTGE).eq.4).and.
     &      (usingDynamicWaterLevelCorrection.eqv..true.) ) then
         if (write_local_files.eqv..true.) then
            call writeDomainHeader(108,trim(localdir)//'/'//
     &           'dynamicWaterlevelCorrection.63',NP_G, NP, 'dynamicWaterlevelCorrection')
           else
              call writeDomainHeader(108,trim(globaldir)//'/'//
     &           'dynamicWaterlevelCorrection.63',NP_G, NP, 'dynamicWaterlevelCorrection')
           endif
      endif
      ! jgf: Added init for station dynamic waterlevel correction output file 
      IF( (ABS(NOUTE).EQ.1).and.(usingDynamicWaterLevelCorrection.eqv..true.)) THEN
         if (write_local_files.eqv..false.) then
            CALL OPEN_GBL_FILE(109, TRIM(GLOBALDIR)//'/'//'dynamicWaterlevelCorrection.61',
     &       NSTAE_G, NSTAE, HEADER61)
         else 
            CALL OPEN_GBL_FILE(109, TRIM(LOCALDIR)//'/'//'dynamicWaterlevelCorrection.61',
     &       NSTAE_G, NSTAE, HEADER61)
         endif
         IESTP=2
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
         if (write_local_files.eqv..false.) then
            CALL OPEN_GBL_FILE(64, TRIM(GLOBALDIR)//'/'//'fort.64',
     &        NP_G, NP, HEADER64)
         else
            CALL OPEN_GBL_FILE(64, TRIM(LOCALDIR)//'/'//'fort.64',
     &        NP_G, NP, HEADER64)
         endif
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


Cobell - Open the fort.77 (weir elevation change) file
         IF(USE_TVW)THEN
            NSCOU_TVW = 0
            IG_TVW = 0
            IF((ABS(NOUT_TVW).EQ.1).OR.(ABS(NOUT_TVW).EQ.4))THEN
                CALL OPEN_GBL_FILE(77,TRIM(GLOBALDIR)//'/'//'fort.77',
     &              NP_G,NP,HEADER77)
                IG_TVW = 2
            ENDIF   
            IF(ABS(NOUT_TVW).EQ.2)THEN
                OPEN(77,FILE=TRIM(LOCALDIR)//'/'//'fort.77',
     &              ACCESS='DIRECT',RECL=NBYTE)
                IF(NBYTE.EQ.4) THEN
                    DO I=1,8
                        WRITE(77,REC=IG_TVW+I) RDES4(I)
                    ENDDO
                    IG_TVW=IG_TVW+8
                    DO I=1,6
                        WRITE(77,REC=IG_TVW+I) RID4(I)
                    ENDDO
                    IG_TVW=IG_TVW+6
                    DO I=1,6
                        WRITE(77,REC=IG_TVW+I) AID4(I)
                    ENDDO
                    IG_TVW=IG_TVW+6
                ENDIF
                IF(NBYTE.EQ.8) THEN
                    DO I=1,4
                        WRITE(77,REC=IG_TVW+I) RDES8(I)
                    ENDDO
                    IG_TVW=IG_TVW+4
                    DO I=1,3
                        WRITE(77,REC=IG_TVW+I) RID8(I)
                    ENDDO
                    IG_TVW=IG_TVW+3
                    DO I=1,3
                        WRITE(77,REC=IG_TVW+I) AID8(I)
                    ENDDO
                    IG_TVW=IG_TVW+3
                ENDIF
                WRITE(77,REC=IG_TVW+1) NDSETS_TVW
                WRITE(77,REC=IG_TVW+2) NP
                WRITE(77,REC=IG_TVW+3) DT*NSPOOL_TVW
                WRITE(77,REC=IG_TVW+4) NSPOOL_TVW
                WRITE(77,REC=IG_TVW+5) 2
                IG_TVW=IG_TVW+5
                CLOSE(77)
            ENDIF    
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
         if (write_local_files.eqv..false.) then
            CALL OPEN_GBL_FILE(73, TRIM(GLOBALDIR)//'/'//'fort.73',
     &         NP_G, NP, HEADER73)
            IGPP=2
            CALL OPEN_GBL_FILE(74, TRIM(GLOBALDIR)//'/'//'fort.74',
     $        NP_G, NP, HEADER74)
            IGWP=2
! tcm v49.64.01 added for ice
            IF(NCICE.NE.0) then
               CALL OPEN_GBL_FILE(93, TRIM(GLOBALDIR)//'/'//'fort.93',
     $           NP_G, NP, HEADER93)
               IGIP=2
            ENDIF
         else
         
            CALL OPEN_GBL_FILE(73, TRIM(LOCALDIR)//'/'//'fort.73',
     &            NP_G, NP, HEADER73)
            IGPP=2
            CALL OPEN_GBL_FILE(74, TRIM(LOCALDIR)//'/'//'fort.74',
     $        NP_G, NP, HEADER74)
            IGWP=2
! tcm v49.64.01 added for ice
            IF(NCICE.NE.0) then
               CALL OPEN_GBL_FILE(93, TRIM(LOCALDIR)//'/'//'fort.93',
     $           NP_G, NP, HEADER93)
               IGIP=2
            ENDIF
         
         endif
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
!      IF ((ABS(NRS).EQ.3).or.(ABS(NRS).eq.4)) then
! tcm v56 2024 added NRS = 1 and NRS = 5
      IF ( (ABS(NRS).EQ.1).or.(ABS(NRS).EQ.3).or.(ABS(NRS).eq.4).or.
     &     (ABS(NRS).EQ.5) ) then
!#ifdef CSWAN  tcm temp
Casey 090302: Added the output of radiation stress gradients.
         IGRadS=0
         IF((ABS(NOUTGW).EQ.1).OR.(ABS(NOUTGW).EQ.4)) THEN
            if (write_local_files.eqv..false.) then
               CALL OPEN_GBL_FILE(164, TRIM(GLOBALDIR)//'/'//'rads.64',
     &             NP_G, NP, HEADER74)
            else
               CALL OPEN_GBL_FILE(164, TRIM(LOCALDIR)//'/'//'rads.64',
     &             NP_G, NP, HEADER74)
            endif
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
         if (write_local_files.eqv..false.) then
            CALL OPEN_GBL_FILE(76, TRIM(GLOBALDIR)//'/'//'fort.76',
     &       NP_G, NP, HEADER63)  !NO NEED TO CREATE HEADER76
         else
            CALL OPEN_GBL_FILE(76, TRIM(LOCALDIR)//'/'//'fort.76',
     &       NP_G, NP, HEADER63)  !NO NEED TO CREATE HEADER76 
         endif 
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
      REAL(8), ALLOCATABLE :: TMP(:,:)
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

      IF (CBaroclinic.and..not.C3DVS) THEN
         IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
            WRITE(16,*) "CBaroclinic is ",CBaroclinic
            WRITE(16,*) "CB3DVS is ",C3DVS
            WRITE(16,*)
     &           "ERROR: 2DDI baroclinic runs not ", 
     &           "currently supported."
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

!        write(*,*) 'file fort.39 opened.'

        READ(39,*) RIVBCTIMINC, RIVBCSTATIM
        RIVBCTIME1=RIVBCSTATIM*86400.d0
        RIVBCTIME2=RIVBCTIME1+RIVBCTIMINC
        DO J=1,totalbcrivernodes

!         write(*,*) 'inside the DO J=1,totalbcrivernodes loop.'

         IF (IDEN.EQ.2) THEN
           READ(39,*) (BCRivSalN1(J,K),K=1,NFEN)
         ELSE IF (IDEN.EQ.3) THEN
           READ(39,*) (BCRivTempN1(J,K),K=1,NFEN)
         ELSE IF (IDEN.EQ.4) THEN
           READ(39,*) (BCRivSalN1(J,K), BCRivTempN1(J,K),K=1,NFEN)



           ! arash debug
           BCRivSalN1(:,:) = 0.0
           BCRivTempN1(:,:) = 0.0



         END IF
        END DO
        DO J=1,totalbcrivernodes
         IF (IDEN.EQ.2) THEN
           READ(39,*) (BCRivSalN2(J,K),K=1,NFEN)
         ELSE IF (IDEN.EQ.3) THEN
           READ(39,*) (BCRivTempN2(J,K),K=1,NFEN)
         ELSE IF (IDEN.EQ.4) THEN
           READ(39,*) (BCRivSalN2(J,K), BCRivTempN2(J,K),K=1,NFEN)


           ! arash debug
           BCRivSalN2(:,:) = 0.0
           BCRivTempN2(:,:) = 0.0


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
!      CALL initOutput3D(44, I3DGD, NDSET3DGD, NP, NP_G,
!     &                  NSPO3DGD, IRType, I3DGDRec)
! arash
      CALL initOutput3D(44, I3DGD, NDSET3DGD, NP, NP_G,
     &                  NSPO3DGD, IDen, I3DGDRec)

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

C----------------------------------------------------------------------
C                      S U B R O U T I N E     
C           S E T   C O L D   W E T   D R Y   S T A T E
C       A N D   W A T E R   S U R F A C E   E L E V A T I O N 
C----------------------------------------------------------------------
C     jgf52.08.01: Sets value of nnodecode at cold start. The values
C     of nnodecode are copied to the nodecode array in totalAreaCalc().
C     This is called from subroutine hotstart() for use in producing 
C     detailed inundation output only. When called from hotstart(),
C     the hotstart wet/dry state and water surface elevation ends up
C     overwriting the values computed here. 
C----------------------------------------------------------------------
      subroutine setColdWetDryStateAndWaterSurfaceElevation()
      use sizes, only : inputdir
      use ADCIRC_MOD, only : ADCIRC_Terminate
      use mesh, only : dp, np, mju, ne, nm
      use boundaries, only : nvel, lbcodei, nbv
      use global, only : eta2, nnodecode, noff, logMessage, allMessage,
     &   scratchMessage, INFO, h0, setMessageSource, unsetMessageSource,
     &   DEBUG, WARNING, INFO, ERROR, scratchMessage, eta1, nolifa,
     &   EtaDisc, River_above_MSL
      use nodalattributes, only : loadStartDry, STARTDRY, GeoidOffset, 
     &   LoadGeoidOffset, OutputTau0, LoadRiver_et_WSE, River_et_WSE,
     &   BFCdlLimit, fric, fgamma, ftheta, manningsn, hbreak,
     &   loadManningsN, iflinbf, ifnlbf, ifhybf
!JLW: adding subgrid
      use subgrid, only: level0, readsubgridlookup
     &                   , getvertlookup
     &                   , initVar
     &                   , gridDepthVertETA2
     &                   , subgridVertList
     &                   , wetFracVertTab
     &                   , numPhi
     &                   , wetDepthVertETA2
     &                   , gridDepthVertTab
      implicit none
      integer :: nm1, nm2, nm3  ! node numbers around the element
      integer :: nc1, nc2, nc3  ! nodecodes around the element
      integer :: ncele          ! wet dry state of the element
      integer :: i
      real(8) :: H1
      real(8) :: HABSMIN
      integer :: numGreater
      integer :: j
      ! jgf52.08.24: Investigate reasons why a node was dried by the
      ! setColdWetDryState() subroutine.  
      character(len=1000), allocatable :: dryReason(:)
      logical, allocatable :: dryReasonInitialized(:)
      logical :: found ! .true. if obsolete init river elevation file was found
      !
      
      call setMessageSource(
     &   "setColdWetDryStateAndWaterSurfaceElevation")
#if defined(COLDSTART_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      HABSMIN = h0*0.8d0

      !
      ! tcm v51.27 20140502 Changed from using a fort.88 file
      ! to a nodal attribute
      IF (LoadRiver_et_WSE.eqv..true.) then
         River_above_MSL=.true.
         call allMessage(INFO,
     &      'Initial river elevations were specified as a nodal '
     &       //'attribute and will be used.')
      ENDIF
      !
      ! tcm v51.27 Check to see if a fort.88 file is supplied
      ! if so issue a warning.      
      inquire(file=trim(inputdir)//'/'//'fort.88',exist=found)
      if (found.eqv..true.) then
         call allMessage(WARNING,
     &    'A fort.88 file was found but is no longer supported. '
     &   //'Initial river elevations must be set as a nodal attribute.')
      endif
      !
      ! jgf46.01 Added the ability to include steric effects.
      if (LoadGeoidOffset.eqv..true.) then
         ETA2(:)=GeoidOffset(:)
      endif
      !
      ! tcm v51.27 20140502 Made the GeoidOffset and Initial River Elevation
      ! an additive quantity instead of the max of the two.
      !
      ! jgf52.08.01: Simplified to just adding the initial river elevation to 
      ! eta2 wherever the initial river elevation is greater than zero.
      if (LoadRiver_et_WSE.eqv..true.) then
         where (River_et_WSE.GT.0.d0) 
            eta2 = eta2 + River_et_WSE 
         end where
         ! tcmv51.27 20140502 Added to remove this nodal attribute
         ! as it is only used at Cold Start
         deallocate (River_et_WSE)
      endif
      !...Set the radiation bc information 
      DO I=1,NVEL
         IF(LBCODEI(I).EQ.52) THEN
            EtaDisc(NBV(I))=Eta2(NBV(I))
         ENDIF
      END DO

!JLW: adding subgrid lookup table read in
      IF(level0.EQV..true.)THEN
#ifdef ADCNETCDF
         ! Set all nodes to wet by default.
         nnodecode(:)=1
         ! set all elements to wet by default.
         noff(:)=1
!JLW: call subgrid lookup table read in
         call readSubgridLookup()
         if(nolifa.eq.2)then
!JLW: initialize subgrid variables
            call initVar()
!JLW: call the vertex lookup subroutine
            call getVertLookup()
            ! jgf52.08.23: Collect, store, and report the reasons why each
            ! node was dried.
            allocate(dryReason(np))
            allocate(dryReasonInitialized(np))
            dryReasonInitialized = .false.
            do i=1,np
               H1 = gridDepthVertETA2(i)
               if (H1.le.h0) then
                  nnodecode(i) = 0
                  write(dryReason(i),'("node=",i0)') i
                  write(scratchMessage,
     &               '(" dp+eta2<=h0: ",f10.3," + ",f10.3," <= ",f6.3)')
     &               dp(i), eta2(i), h0
                  dryReason(i) = trim(dryReason(i)) // trim(scratchMessage)
                  dryReasonInitialized(i) = .true.
               endif
               if ((loadStartDry.eqv..true.).and.(startDry(i).eq.1)) then
                  nnodecode(i) = 0
                  if (dryReasonInitialized(i).eqv..false.) then
                     write(dryReason(i),'("node=",i0)') i
                     dryReasonInitialized(i) = .true.
                  else
                     dryReason(i) = trim(dryReason(i)) // ", "
                  endif
                  dryReason(i) = trim(dryReason(i)) // "startDry=1"
               endif
            end do
         !
            ! Dry any landlocked nodes by checking that they are connected to at
            ! least 1 functioning element.
            mju(:)=0
            do i=1,ne
               nm1=nm(i,1)
               nm2=nm(i,2)
               nm3=nm(i,3)
               nc1=nnodecode(nm1)
               nc2=nnodecode(nm2)
               nc3=nnodecode(nm3)
               ncele=nc1*nc2*nc3*noff(i)
               mju(nm1)=mju(nm1)+ncele
               mju(nm2)=mju(nm2)+ncele
               mju(nm3)=mju(nm3)+ncele
            enddo
            do i=1,np
               if((nnodecode(i).eq.1).and.(mju(i).eq.0)) then
                  nnodecode(i)=0
                  write(scratchMessage,9883) I
                  call allMessage(INFO,scratchMessage)
                  if (dryReasonInitialized(i).eqv..false.) then
                     write(dryReason(i),'("node=",i0)') i
                     dryReasonInitialized(i) = .true.
                  else
                     dryReason(i) = trim(dryReason(i)) // ", "
                  endif
                  dryReason(i) = trim(dryReason(i)) // " landlocked"
               endif
            enddo
! 9883       format('Node ',i0,' dried (landlocking).')
         !
            !jgf52.08.23: Report the reasons why a node was initially dried.
            call logMessage(DEBUG,
     &         "Nodes were set to an initially dry state "
     &         // "for the following reasons:")
            do i=1,np
               if (nnodecode(i).eq.0) then
                  if (dryReasonInitialized(i).eqv..false.) then
                     write(dryReason(i),
     &            '("node=",i0," This node is dry for no good reason.")') i
                  dryReasonInitialized(i) = .true.
                  endif
                  call logMessage(DEBUG,trim(dryReason(i)))
               endif
            end do
            ! release memory
            deallocate(dryReason)
            deallocate(dryReasonInitialized)
            !
            ! Finish initializing water surface elevation (eta2 and eta1)
            do i=1,np
               if(nnodecode(i).eq.0)then
                  if(gridDepthVertETA2(i).lt.HABSMIN)then
                     gridDepthVertETA2(i) = HABSMIN
                     wetDepthVertETA2(i) = HABSMIN
                     if(subgridVertList(i).eq.1)then
                        numGreater = 0
                        do j=1,numPhi
                           if(gridDepthVertETA2(i).lt.
     &                        gridDepthVertTab(i,j))then
                              numGreater = numGreater + 1
                           endif
                        enddo
                        if(numGreater.eq.numPhi)then
                           eta2(i) = wetFracVertTab(i,1)
                        elseif(numGreater.eq.0)then
                           eta2(i) = wetFracVertTab(i,numPhi)
     &                        + (HABSMIN-gridDepthVertTab(i,numPhi))
                        else
                        !JLW: now back calculate eta2 based on this change
                           eta2(i) = ((gridDepthVertETA2(i)
     &                     - gridDepthVertTab(i,numPhi-numGreater))
     &                     / (gridDepthVertTab(i,numPhi-numGreater+1)
     &                     - gridDepthVertTab(i,numPhi-numGreater))
     &                     * (wetFracVertTab(i,numPhi-numGreater+1)
     &                     - wetFracVertTab(i,numPhi-numGreater))
     &                     + wetFracVertTab(i,numPhi-numGreater))
                        endif
                     else
                           eta2(i)=HABSMIN-dp(i)
                     endif
                  endif
               endif
            enddo
         endif    
#else
         call allMessage(ERROR, "Subgrid calculation requires building"
     &                          //"with netcdf enabled")
         call ADCIRC_Terminate()
#endif
      ELSE
      !
      ! Set all nodes to wet by default.
      nnodecode(:)=1
      ! set all elements to wet by default.
      noff(:)=1        
      !
      ! If wetting and drying is active, set the nodal wet dry state
      if (nolifa.eq.2) then
         !
         ! jgf52.08.23: Collect, store, and report the reasons why each
         ! node was dried. 
         allocate(dryReason(np))
         allocate(dryReasonInitialized(np))
         dryReasonInitialized = .false.
         do i=1,np      
            ! set dry areas based on eta2
            if (dp(i)+eta2(i).le.h0) then
               nnodecode(i) = 0
               write(dryReason(i),'("node=",i0)') i
               write(scratchMessage,
     &            '(" dp+eta2<=h0: ",f10.3," + ",f10.3," <= ",f6.3)') 
     &            dp(i), eta2(i), h0
               dryReason(i) = trim(dryReason(i)) // trim(scratchMessage)
               dryReasonInitialized(i) = .true.
            endif
            if ((loadStartDry.eqv..true.).and.(startDry(i).eq.1)) then
               nnodecode(i) = 0
               if (dryReasonInitialized(i).eqv..false.) then
                  write(dryReason(i),'("node=",i0)') i
                  dryReasonInitialized(i) = .true.
               else
                  dryReason(i) = trim(dryReason(i)) // ", "
               endif
               dryReason(i) = trim(dryReason(i)) // "startDry=1"
            endif       
         end do
         !
         ! Dry any landlocked nodes by checking that they are connected to at
         ! least 1 functioning element.
         mju(:)=0
         do i=1,ne
            nm1=nm(i,1)
            nm2=nm(i,2)
            nm3=nm(i,3)
            nc1=nnodecode(nm1)
            nc2=nnodecode(nm2)
            nc3=nnodecode(nm3)
            ncele=nc1*nc2*nc3*noff(i)
            mju(nm1)=mju(nm1)+ncele
            mju(nm2)=mju(nm2)+ncele
            mju(nm3)=mju(nm3)+ncele
         enddo
         do i=1,np
            if((nnodecode(i).eq.1).and.(mju(i).eq.0)) then
               nnodecode(i)=0
               write(scratchMessage,9883) I
               call allMessage(INFO,scratchMessage)
               if (dryReasonInitialized(i).eqv..false.) then
                  write(dryReason(i),'("node=",i0)') i
                  dryReasonInitialized(i) = .true.
               else
                  dryReason(i) = trim(dryReason(i)) // ", "
               endif
               dryReason(i) = trim(dryReason(i)) // " landlocked"
            endif
         enddo
 9883    format('Node ',i0,' dried (landlocking).')
         !
         !jgf52.08.23: Report the reasons why a node was initially dried.
         call logMessage(DEBUG,
     &      "Nodes were set to an initially dry state "
     &      // "for the following reasons:")
         do i=1,np
            if (nnodecode(i).eq.0) then
               if (dryReasonInitialized(i).eqv..false.) then
                  write(dryReason(i),
     &            '("node=",i0," This node is dry for no good reason.")') i
                  dryReasonInitialized(i) = .true.
               endif
               call logMessage(DEBUG,trim(dryReason(i)))
            endif
         end do 
         ! release memory 
         deallocate(dryReason)
         deallocate(dryReasonInitialized)
         !
         ! Finish initializing water surface elevation (eta2 and eta1)
         where (nnodecode.eq.0) 
            eta2=h0-dp
         end where
      endif
!JLW: end subgrid if statement
      ENDIF
      eta1(:) = eta2(:) 

#if defined(COLDSTART_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C----------------------------------------------------------------------
      end subroutine setColdWetDryStateAndWaterSurfaceElevation
C----------------------------------------------------------------------

