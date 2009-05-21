C----------------------------------------------------------------------
C                     W R I T E   O U T P U T . F  
C----------------------------------------------------------------------
C jgf48.03 This file contains all the subroutines that write output to 
C files. It was created to consolidate all the output routines that 
C were located in other parts of the code, to make it easier to incorporate
C them into globalio, and to make it easier to add new output modes or 
C formats in the future.
C----------------------------------------------------------------------


C----------------------------------------------------------------------
C       S U B R O U T I N E   W R I T E _ O U T P U T _ 2 D 
C----------------------------------------------------------------------
C
C     R.L. 8/22/05 Subroutine to write primary 2D model output not
C     including hotstart and harmonic analysis.
C
C----------------------------------------------------------------------
      SUBROUTINE WRITE_OUTPUT_2D (IT,Time)
      USE SIZES, ONLY : SZ, INPUTDIR, NBYTE
      USE GLOBAL
      USE GLOBAL_IO, ONLY: storeOne, storeTwo, store63, HEADER_MAX,
     &    open_gbl_file, write_gbl_file, write_gbl_63_skip_drynode,
     &    write_gbl_file_skip_default, collectFullDomainArray,
     &    packOne, unpackOne, packTwo, unpackTwo, pack63, pack83 
      USE NodalAttributes, ONLY : OutputTau0, Tau0Var
#ifdef CMPI
      USE WRITER, ONLY : WRITE_GBL_FILE_THROUGH_WRITER, 
     $     FLUSH_WRITERS, WRITER_STOREONE, 
     $     WRITER_STORETWO, WRITER_STORE63, WRITER_INIT
#endif
#ifdef NETCDF
      USE NETCDFIO, ONLY :
     &     timenc, timenc_len, fileunit, irtype, 
     &     xenc, yenc, xvnc, yvnc, xpnc, ypnc, xnc, ync, xwnc, ywnc,
     &     et00nc, UU00nc, VV00nc, RMP00nc, RMU00nc, RMV00nc, nbvvnc,
     &     eta22, uu22, vv22, pr22, WNCX22, WNCY22, 
     &     nbounc,nvelnc,nvellnc,
     &     nopenc,netanc,nvdllnc,
     &     writenc61, writenc62, writenc63, writenc64,
     &     writenc71, writenc72, writenc73, writenc74 
       USE NodalAttributes, ONLY : 
     &     nolibf, nwp, tau0, cf, eslm
#endif

      IMPLICIT NONE 
      INTEGER, intent(in) :: IT
      REAL(8), intent(in) :: Time

      INTEGER I                       !local loop counters
      INTEGER NC1, NC2, NC3, NCEle
      INTEGER NM1, NM2, NM3

      REAL(SZ) C1, C2, C3
      REAL(SZ) EE1, EE2, EE3
      REAL(SZ) H2, H2N1, H2N2, H2N3
      REAL(SZ) U11, U22, U33
      REAL(SZ) V11, V22, V33

      type(OutputDataDescript_t) :: descript
      type(OutputDataDescript_t), SAVE :: ElevStaDescript
      type(OutputDataDescript_t), SAVE :: VelStaDescript
      type(OutputDataDescript_t), SAVE :: ElevDescript
      type(OutputDataDescript_t), SAVE :: Tau0Descript
      type(OutputDataDescript_t), SAVE :: VelDescript
      type(OutputDataDescript_t), SAVE :: PrStaDescript
      type(OutputDataDescript_t), SAVE :: WindVelStaDescript
      type(OutputDataDescript_t), SAVE :: PrDescript
      type(OutputDataDescript_t), SAVE :: WindVelDescript
      type(OutputDataDescript_t), SAVE :: ConcStaDescript
      type(OutputDataDescript_t), SAVE :: ConcDescript

C     jgf48.03 Full domain arrays used when writing globalio output.
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: ET00_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: UU00_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: VV00_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: ETA2_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: TAU0VAR_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: UU2_g(:)     
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: VV2_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: RMP00_g(:)     
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: RMU00_g(:)     
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: RMV00_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: Pr2_g(:)     
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: WVNXOut_g(:)     
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: WVNYOut_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: CC00_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: CH1_g(:)

      LOGICAL, SAVE :: FirstCall = .true.
      
      IF (FirstCall) THEN

         ! fort.61
         write(16,*) 'Allocating ElevStaDescript' !jgfdebug48.03
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(ET00_g(NSTAE_G)) 
         ENDIF
         ElevStaDescript % specifier            = NOUTE
         ElevStaDescript % initial_value        = 0.0
         ElevStaDescript % num_items_per_record = 1
         ElevStaDescript % num_fd_records       = NSTAE_G
         ElevStaDescript % num_records_this     = NSTAE
         ElevStaDescript % imap                 => IMAP_STAE_LG
         ElevStaDescript % array                => ET00
         ElevStaDescript % array_g              => ET00_g

         ! fort.62
         write(16,*) 'Allocating VelStaDescript' !jgfdebug48.03
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(UU00_g(NSTAV_G))
            ALLOCATE(VV00_g(NSTAV_G)) 
         ENDIF
         VelStaDescript % specifier            = NOUTV
         VelStaDescript % initial_value        = 0.0
         VelStaDescript % num_items_per_record = 2
         VelStaDescript % num_fd_records       = NSTAV_G
         VelStaDescript % num_records_this     = NSTAV
         VelStaDescript % imap                 => IMAP_STAV_LG
         VelStaDescript % array                => UU00
         VelStaDescript % array2               => VV00
         VelStaDescript % array_g              => UU00_g
         VelStaDescript % array2_g             => VV00_g

         ! fort.63
         write(16,*) 'Allocating ElevDescript' !jgfdebug48.03
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(ETA2_g(NP_G))
         ENDIF
         ElevDescript % specifier            = NOUTGE
         ElevDescript % initial_value        = 0.0
         ElevDescript % num_items_per_record = 1
         ElevDescript % num_fd_records       = NP_G
         ElevDescript % num_records_this     = NP
         ElevDescript % imap                 => NODES_LG
         ElevDescript % array                => ETA2
         ElevDescript % array_g              => ETA2_g

         ! tau0
         write(16,*) 'Allocating Tau0Descript' !jgfdebug48.03
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(TAU0VAR_g(NP_G))
         ENDIF
         Tau0Descript % specifier            = 1  ! hard-code to ascii text
         Tau0Descript % initial_value        = 0.0
         Tau0Descript % num_items_per_record = 1
         Tau0Descript % num_fd_records       = NP_G
         Tau0Descript % num_records_this     = NP
         Tau0Descript % imap                 => NODES_LG
         Tau0Descript % array                => TAU0VAR
         Tau0Descript % array_g              => TAU0VAR_g

         ! fort.64
         write(16,*) 'Allocating VelDescript' !jgfdebug48.03
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(UU2_g(NP_G))
            ALLOCATE(VV2_g(NP_G)) 
         ENDIF
         VelDescript % specifier            = NOUTGV
         VelDescript % initial_value        = 0.0
         VelDescript % num_items_per_record = 2
         VelDescript % num_fd_records       = NP_G
         VelDescript % num_records_this     = NP
         VelDescript % imap                 => NODES_LG
         VelDescript % array                => UU2
         VelDescript % array2               => VV2
         VelDescript % array_g              => UU2_g
         VelDescript % array2_g             => VV2_g

         ! fort.71
         write(16,*) 'Allocating PrStaDescript' !jgfdebug48.03
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(RMP00_g(NSTAM_G))
         ENDIF
         PrStaDescript % specifier            = NOUTM
         PrStaDescript % initial_value        = 0.0
         PrStaDescript % num_items_per_record = 1
         PrStaDescript % num_fd_records       = NSTAM_G
         PrStaDescript % num_records_this     = NSTAM
         PrStaDescript % imap                 => IMAP_STAM_LG
         PrStaDescript % array                => RMP00
         PrStaDescript % array_g              => RMP00_g
     
         ! fort.72
         write(16,*) 'Allocating WindVelStaDescript' !jgfdebug48.03
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(RMU00_g(NSTAM_G))
            ALLOCATE(RMV00_g(NSTAM_G))
         ENDIF
         WindVelStaDescript % specifier            = NOUTM
         WindVelStaDescript % initial_value        = 0.0
         WindVelStaDescript % num_items_per_record = 2
         WindVelStaDescript % num_fd_records       = NSTAM_G
         WindVelStaDescript % num_records_this     = NSTAM
         WindVelStaDescript % imap                 => IMAP_STAM_LG
         WindVelStaDescript % array                => RMU00
         WindVelStaDescript % array_g              => RMU00_g
         WindVelStaDescript % array2               => RMV00
         WindVelStaDescript % array2_g             => RMV00_g

         ! fort.73
         write(16,*) 'Allocating PrDescript' !jgfdebug48.03
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(Pr2_g(NP_G))
         ENDIF
         PrDescript % specifier            = NOUTGW
         PrDescript % initial_value        = 0.0
         PrDescript % num_items_per_record = 1
         PrDescript % num_fd_records       = NP_G
         PrDescript % num_records_this     = NP
         PrDescript % imap                 => NODES_LG
         PrDescript % array                => Pr2
         PrDescript % array_g              => Pr2_g

         ! fort.74
         write(16,*) 'Allocating WindVelDescript' !jgfdebug48.03
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(WVNXOut_g(NP_G))
            ALLOCATE(WVNYOut_g(NP_G)) 
         ENDIF
         WindVelDescript % specifier            = NOUTGW
         WindVelDescript % initial_value        = 0.0
         WindVelDescript % num_items_per_record = 2
         WindVelDescript % num_fd_records       = NP_G
         WindVelDescript % num_records_this     = NP
         WindVelDescript % imap                 => NODES_LG
         WindVelDescript % array                => WVNXOut
         WindVelDescript % array_g              => WVNXOut_g
         WindVelDescript % array2               => WVNYOut
         WindVelDescript % array2_g             => WVNYOut_g

         ! fort.81
         write(16,*) 'Allocating ConcStaDescript' !jgfdebug48.03
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(CC00_g(NP_G))
         ENDIF
         ConcStaDescript % specifier            = NOUTC
         ConcStaDescript % initial_value        = 0.0
         ConcStaDescript % num_items_per_record = 1
         ConcStaDescript % num_fd_records       = NSTAC_G
         ConcStaDescript % num_records_this     = NSTAC
         ConcStaDescript % imap                 => IMAP_STAC_LG
         ConcStaDescript % array                => CC00
         ConcStaDescript % array_g              => CC00_g        

         ! fort.83
         write(16,*) 'Allocating ConcDescript' !jgfdebug48.03
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(CH1_g(NP_G))
         ENDIF
         ConcDescript % specifier            = NOUTGC
         ConcDescript % initial_value        = 0.0
         ConcDescript % num_items_per_record = 1
         ConcDescript % num_fd_records       = NP_G
         ConcDescript % num_records_this     = NP
         ConcDescript % imap                 => NODES_LG
         ConcDescript % array                => CH1
         ConcDescript % array_g              => CH1_g        

         FirstCall = .false.
      ENDIF

C...  Output elevation recording station information if noute<>0 and the
C...  time step falls within the specified window calculate elevation
C...  solutions at stations using interpolation

      IF(NOUTE.NE.0) THEN
         IF((IT.GT.NTCYSE).AND.(IT.LE.NTCYFE)) THEN
             NSCOUE=NSCOUE+1
             NCCOUE=NSCOUE
         ENDIF
         IF(NSCOUE.EQ.NSPOOLE) THEN
            DO I=1,NSTAE
               EE1=ETA2(NM(NNE(I),1))
               EE2=ETA2(NM(NNE(I),2))
               EE3=ETA2(NM(NNE(I),3))
               NC1=NODECODE(NM(NNE(I),1))
               NC2=NODECODE(NM(NNE(I),2))
               NC3=NODECODE(NM(NNE(I),3))
               NCELE=NC1*NC2*NC3*NOFF(NNE(I))
               IF(NCELE.EQ.1) ET00(I)=EE1*STAIE1(I)+EE2*STAIE2(I)
     &                                             +EE3*STAIE3(I)
               IF(NCELE.EQ.0) ET00(I)=-99999.
            END DO
C
C           save elevation recording station data to file
            CALL writeOutArray(61, TIME, IT, ElevStaDescript,
     &                         packOne, unpackOne, IESTP) 

            write(16,*) 'Finished with fort.61 output routine.' !jgfdebug48.03
            NSCOUE=0
            NCCOUE=NSCOUE
         ENDIF
      ENDIF

C...  OUTPUT VELOCITY RECORDING STATION TIME SERIES INFORMATION IF
C...  NOUTV<>0 AND THE TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
C...  CALCULATE VELOCITY SOLUTIONS AT STATIONS USING INTERPOLATION
C...  
      IF(NOUTV.NE.0) THEN
         IF((IT.GT.NTCYSV).AND.(IT.LE.NTCYFV)) THEN
             NSCOUV=NSCOUV+1
             NCCOUV=NSCOUV
         ENDIF
         IF(NSCOUV.EQ.NSPOOLV) THEN
            write(16,*) 'About to write fort.62 file.' !jgfdebug48.03
            DO I=1,NSTAV
               U11=UU2(NM(NNV(I),1))
               U22=UU2(NM(NNV(I),2))
               U33=UU2(NM(NNV(I),3))
               V11=VV2(NM(NNV(I),1))
               V22=VV2(NM(NNV(I),2))
               V33=VV2(NM(NNV(I),3))
               UU00(I)=U11*STAIV1(I)+U22*STAIV2(I)+U33*STAIV3(I)
               VV00(I)=V11*STAIV1(I)+V22*STAIV2(I)+V33*STAIV3(I)
            END DO
C
C           save velocity recording station data to file
            CALL writeOutArray(62, TIME, IT, VelStaDescript,
     &                         packTwo, unpackTwo, IVSTP) 

            write(16,*) 'Finished with fort.62 output routine.' !jgfdebug48.03
            NSCOUV=0
            NCCOUV=NSCOUV
         ENDIF
      ENDIF

C...  OUTPUT CONCENTRATION RECORDING STATION INFORMATION IF NOUTC<>0
C...  AND THE TIME STEP FALLS WITHIN THE SPECIFIED WINDOW CALCULATE
C...  CONCENTRATION SOLUTIONS AT STATIONS USING INTERPOLATION
C...  
      IF(NOUTC.NE.0) THEN
         IF((IT.GT.NTCYSC).AND.(IT.LE.NTCYFC)) NSCOUC=NSCOUC+1
         IF(NSCOUC.EQ.NSPOOLC) THEN
            DO I=1,NSTAC
               NM1=NM(NNC(I),1)
               NM2=NM(NNC(I),2)
               NM3=NM(NNC(I),3)
               H2N1=DP(NM1)+IFNLFA*ETA2(NM1)
               H2N2=DP(NM2)+IFNLFA*ETA2(NM2)
               H2N3=DP(NM3)+IFNLFA*ETA2(NM3)
               C1=CH1(NM1)/H2N1
               C2=CH1(NM2)/H2N2
               C3=CH1(NM3)/H2N3
               NC1=NODECODE(NM1)
               NC2=NODECODE(NM2)
               NC3=NODECODE(NM3)
               NCELE=NC1*NC2*NC3*NOFF(NNC(I))
               IF(NCELE.EQ.1) CC00(I)=C1*STAIC1(I)+C2*STAIC2(I)
     &                                            +C3*STAIC3(I)
               IF(NCELE.EQ.0) CC00(I)=-99999.
            END DO
C
C           save concentration recording station data to file
            CALL writeOutArray(81, TIME, IT, ConcStaDescript,
     &                         packOne, unpackOne, ICSTP)
 
            write(16,*) 'Finished with fort.81 output routine.' !jgfdebug48.03
            NSCOUC=0
         ENDIF
      ENDIF

C...  OUTPUT METEOROLOGICAL RECORDING STATION INFORMATION IF NWS>0 AND
C...  THE TIME STEP FALLS WITHIN THE SPECIFIED WINDOW CALCULATE
C...  METEOROLOGICAL SOLUTIONS AT STATIONS USING INTERPOLATION
C...  
      IF((NWS.NE.0).AND.(NOUTM.NE.0)) THEN
         IF((IT.GT.NTCYSM).AND.(IT.LE.NTCYFM)) THEN
            NSCOUM=NSCOUM+1
            NCCOUM=NSCOUM
         ENDIF
         IF(NSCOUM.EQ.NSPOOLM) THEN
            write(16,*) 'About to write fort.71,72 files.' !jgfdebug48.03
            DO I=1,NSTAM
               NM1=NM(NNM(I),1)
               NM2=NM(NNM(I),2)
               NM3=NM(NNM(I),3)
               U11=wvnxout(NM1)
               U22=wvnxout(NM2)
               U33=wvnxout(NM3)
               V11=wvnyout(NM1)
               V22=wvnyout(NM2)
               V33=wvnyout(NM3)
               P11=PR2(NM1)
               P22=PR2(NM2)
               P33=PR2(NM3)
               RMU00(I)=U11*STAIM1(I)+U22*STAIM2(I)+U33*STAIM3(I)
               RMV00(I)=V11*STAIM1(I)+V22*STAIM2(I)+V33*STAIM3(I)
               RMP00(I)=P11*STAIM1(I)+P22*STAIM2(I)+P33*STAIM3(I)
            END DO
C
C           save atmospheric pressure recording station data to file
            CALL writeOutArray(71, TIME, IT, PrStaDescript,
     &                         packOne, unpackOne, IPSTP) 
C
C           save wind velocity recording station data to file
            CALL writeOutArray(72, TIME, IT, WindVelStaDescript,
     &                         packTwo, unpackTwo, IWSTP) 
            NSCOUM=0
            NCCOUM=NSCOUM
         ENDIF
      ENDIF
C...   OUTPUT GLOBAL ELEVATION DATA IF NOUTGE<>0 AND THE
C.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
C...  
      IF(NOUTGE.NE.0) THEN
!         IF(MYPROC.EQ.0) PRINT *, "W: NSCOUGE = ", NSCOUGE
!         IF(MYPROC.EQ.0) PRINT *, "W: NSPOOLGE = ", NSPOOLGE
!         IF(MYPROC.EQ.0) PRINT *, "W: NTCYSGE = ", NTCYSGE
!         IF(MYPROC.EQ.0) PRINT *, "W: NTCYFGE = ", NTCYFGE
         IF((IT.GT.NTCYSGE).AND.(IT.LE.NTCYFGE)) THEN
            NSCOUGE=NSCOUGE+1
            NCCOUGE=NSCOUGE
         ENDIF
         IF(NSCOUGE.EQ.NSPOOLGE) THEN
!         IF(MYPROC.EQ.0) PRINT *, "W: B-write ", NSCOUGE
C
C           save fulldomain elevation data to file
            CALL writeOutArray(63, TIME, IT, ElevDescript,
     &                         pack63, unpackOne, IGEP) 
C
C           jgf47.06: Write out tau0 if it is time varying and the user has 
C           requested it with tau0 = -4 in fort.15.
            IF (OUTPUTTAU0) THEN
               CALL writeOutArray(90, TIME, IT, Tau0Descript,
     &                            packOne, unpackOne, IGEP) 
            ENDIF

C          jgf47.08 Made sparse output an option (if NOUTGE=-4 or 4 in fort.15)
           IF(ABS(NOUTGE).EQ.4) THEN
C              jgf47.08 Create compact output 
C              jgf48.03 TODO: Make this globalio compliant.
#ifdef CMPI
              IF(MNWPROC.GT.0) THEN
                 CALL WRITE_GBL_FILE_THROUGH_WRITER
     $              (TRIM(GLOBALDIR)//'/'//'fort.63',descript,time,it,
     $                 writer_store63, -99999.D0)
              ELSE
                 CALL WRITE_GBL_63_SKIP_DRYNODE  ! sb 11/10/2006
     $                (63, TRIM(GLOBALDIR)//'/'//'fort.63',
     $                 descript, time,it, store63)
              ENDIF
#else
              CALL WRITE_GBL_63_SKIP_DRYNODE  ! sb 11/10/2006
     $            (63, TRIM(GLOBALDIR)//'/'//'fort.63',
     &               descript, time,it, store63)
#endif
              IGEP=IGEP+1+NP
           ENDIF
           write(16,*) 'Finished with fort.63 output routine.' !jgfdebug48.03
           NSCOUGE=0
           NCCOUGE=NSCOUGE
!         IF(MYPROC.EQ.0) PRINT *, "W: B-write after ", NSCOUGE
         ENDIF
      ENDIF

C...  OUTPUT GLOBAL VELOCITY DATA IF NOUTGV<>0 AND THE
C.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
C...  
      IF(NOUTGV.NE.0) THEN
         IF((IT.GT.NTCYSGV).AND.(IT.LE.NTCYFGV)) THEN
            NSCOUGV=NSCOUGV+1
            NCCOUGV=NSCOUGV
         ENDIF
         IF(NSCOUGV.EQ.NSPOOLGV) THEN

            CALL writeOutArray(64, TIME, IT, VelDescript,
     &                      packTwo, unpackTwo, IGVP) 

            IF (ABS(NOUTGV).EQ.4) THEN
            ! jgf48.03 TODO: Make this globalio compliant.
#ifdef CMPI
                  IF(MNWPROC.GT.0) THEN
                    CALL WRITE_GBL_FILE_THROUGH_WRITER
     $                (TRIM(GLOBALDIR)//'/'//'fort.64',descript,time,it,
     $                 writer_storeTwo, 0.D0)
                  ELSE
                    CALL WRITE_GBL_FILE_SKIP_DEFAULT ! sb 11/10/2006
     $                   (64, TRIM(GLOBALDIR)//'/'//'fort.64',
     $                   descript, time, it, storeTwo,0.d0)
                  ENDIF
#else
                  CALL WRITE_GBL_FILE_SKIP_DEFAULT ! sb 11/10/2006
     $                 (64, TRIM(GLOBALDIR)//'/'//'fort.64',
     $                 descript, time, it, storeTwo,0.d0)
#endif
            ENDIF
            NSCOUGV=0
            NCCOUGV=NSCOUGV
         ENDIF
      ENDIF

C...
C...  OUTPUT GLOBAL WIND STRESS and atmospheric pressure data IF
C.... NOUTGW<>0 AND THE TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
C...
      IF((NWS.NE.0).AND.(NOUTGW.NE.0)) THEN
         IF((IT.GT.NTCYSGW).AND.(IT.LE.NTCYFGW)) THEN
             NSCOUGW=NSCOUGW+1
             NCCOUGW=NSCOUGW
         ENDIF
         IF(NSCOUGW.EQ.NSPOOLGW) THEN
C
C           save atmospheric pressure recording station data to file
            CALL writeOutArray(73, TIME, IT, PrDescript,
     &                         packOne, unpackOne, IGPP) 
C
C           save wind velocity recording station data to file
            CALL writeOutArray(74, TIME, IT, WindVelDescript,
     &                         packTwo, unpackTwo, IGWP) 

C           jgf48.03 TODO: Make this globalio compliant.
            IF (ABS(NOUTGW).EQ.4) THEN
#ifdef CMPI
               IF(MNWPROC.GT.0) THEN
                 CALL WRITE_GBL_FILE_THROUGH_WRITER
     $             (TRIM(GLOBALDIR)//'/'//'fort.73',descript,time,it,
     $              writer_storeOne, 0.D0)
               ELSE
                 CALL WRITE_GBL_FILE_SKIP_DEFAULT
     $             (73, TRIM(GLOBALDIR)//'/'//'fort.73',
     $              descript, time, it, storeOne,0.d0)
               ENDIF
#else
               CALL WRITE_GBL_FILE_SKIP_DEFAULT
     $            (73, TRIM(GLOBALDIR)//'/'//'fort.73',
     $               descript, time, it, storeOne,0.d0)
#endif
            ENDIF
C           jgf48.03 TODO: Make this globalio compliant.
            IF (ABS(NOUTGW).EQ.4) THEN
#ifdef CMPI
                  IF(MNWPROC.GT.0) THEN
                    CALL WRITE_GBL_FILE_THROUGH_WRITER
     $                   (TRIM(GLOBALDIR)//'/'//'fort.74',
     &                   descript,time,it,
     $                   writer_storeTwo, 0.D0)
                  ELSE
                    CALL WRITE_GBL_FILE_SKIP_DEFAULT
     $                   (74, TRIM(GLOBALDIR)//'/'//'fort.74',
     $                   descript, time, it, storeTwo,0.d0)
                  ENDIF
#else
                  CALL WRITE_GBL_FILE_SKIP_DEFAULT
     $              (74, TRIM(GLOBALDIR)//'/'//'fort.74',
     $              descript, time, it, storeTwo,0.d0)
#endif
            ENDIF
            NSCOUGW=0
            NCCOUGW=NSCOUGW
         ENDIF
      ENDIF

C...  OUTPUT GLOBAL CONCENTRATION DATA IF NOUTGC<>0 AND THE
C.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
C...  
      IF(NOUTGC.NE.0) THEN
         IF((IT.GT.NTCYSGC).AND.(IT.LE.NTCYFGC)) NSCOUGC=NSCOUGC+1
         IF(NSCOUGC.EQ.NSPOOLGC) THEN
C
C           save fulldomain concentration data to file
C           jgf48.03 TODO: add the proper processingi (see below).
            CALL writeOutArray(83, TIME, IT, ConcDescript,
     &                         pack83, unpackOne, IGCP) 
C               DO I=1,NP
C                  H2=DP(I)+IFNLFA*ETA2(I)
C                  C1=CH1(I)/H2
C                  IF(NODECODE(I).EQ.1) WRITE(83,2453) I,C1
C                  IF(NODECODE(I).EQ.0) WRITE(83,2453) I,-99999.
C               ENDDO
            NSCOUGC=0
         ENDIF
      ENDIF

#ifdef CMPI
Csb 02/09/2007
C...  Communicate writer processors and let them write out what
C...  they have gotton.
      CALL FLUSH_WRITERS()
#endif
      

Csb 11/11/2006
C...   OUTPUT MAXIMUM/MINIMUM GLOBAL DATA AT THE LAST TIME STEP.
C...  

      IF(IT.EQ.NT) THEN
         CALL OPEN_GBL_FILE(963, TRIM(GLOBALDIR)//'/'//'maxele.63',
     $        NP_G, NP, HEADER_MAX)
         descript % initial_value        = 0.0
         descript % num_items_per_record = 1
         descript % num_fd_records       = NP_G
         descript % num_records_this     = NP
         descript % imap                 => NODES_LG
         descript % array                => ETAMAX
         CALL WRITE_GBL_FILE(963, TRIM(GLOBALDIR)//'/'//'maxele.63',
     $        descript, time, it, storeOne)
         CALL OPEN_GBL_FILE(963, TRIM(GLOBALDIR)//'/'//'maxvel.63',
     $        NP_G, NP, HEADER_MAX)
         descript % initial_value        = 0.0
         descript % num_items_per_record = 1
         descript % num_fd_records       = NP_G
         descript % num_records_this     = NP
         descript % imap                 => NODES_LG
         descript % array                => UMAX
         CALL WRITE_GBL_FILE(963, TRIM(GLOBALDIR)//'/'//'maxvel.63',
     $        descript, time, it, storeOne)
         CALL OPEN_GBL_FILE(963, TRIM(GLOBALDIR)//'/'//'minpr.63',
     $        NP_G, NP, HEADER_MAX)
         descript % initial_value        = 0.0
         descript % num_items_per_record = 1
         descript % num_fd_records       = NP_G
         descript % num_records_this     = NP
         descript % imap                 => NODES_LG
         descript % array                => PRMIN
         CALL WRITE_GBL_FILE(963, TRIM(GLOBALDIR)//'/'//'minpr.63',
     $        descript, time, it, storeOne)
         CALL OPEN_GBL_FILE(963, TRIM(GLOBALDIR)//'/'//'maxwvel.63',
     $        NP_G, NP, HEADER_MAX)
         descript % initial_value        = 0.0
         descript % num_items_per_record = 1
         descript % num_fd_records       = NP_G
         descript % num_records_this     = NP
         descript % imap                 => NODES_LG
         descript % array                => WVNOUTMAX
         CALL WRITE_GBL_FILE(963, TRIM(GLOBALDIR)//'/'//'maxwvel.63',
     $        descript, time, it, storeOne)
         CALL OPEN_GBL_FILE(963, TRIM(GLOBALDIR)//'/'//'maxrs.63',
     $        NP_G, NP, HEADER_MAX)
         descript % initial_value        = 0.0
         descript % num_items_per_record = 1
         descript % num_fd_records       = NP_G
         descript % num_records_this     = NP
         descript % imap                 => NODES_LG
         descript % array                => RSNMAX
         CALL WRITE_GBL_FILE(963, TRIM(GLOBALDIR)//'/'//'maxrs.63',
     $        descript, time, it, storeOne)
      ENDIF

      RETURN
 2120 FORMAT(2X,E20.10,5X,I10)
 2453 FORMAT(2x, i8, 2x, E20.10, E20.10, E20.10, E20.10) 
 2454 FORMAT(2X,I8,2(2X,E15.8))
C-----------------------------------------------------------------------
      END SUBROUTINE WRITE_OUTPUT_2D
C-----------------------------------------------------------------------


C********************************************************
C  Subroutine to generate 3D Model Output 
C                                        
C                          
C  Note, the FORM='FORMATTED' field in the binary open statements 
C  means that the record length (RECL) is specified in bytes.  
C  If this field is omitted, the record length must be specified 
C  in 4 byte words. 
C                   
C  r.l.  8/23/2005      
C*********************************************************
C-----------------------------------------------------------------------
C     S U B R O U T I N E   W R I T E  O U T P U T  3 D  
C-----------------------------------------------------------------------
C     jgf48.11 Moved from vsmy.F to write_output.F. 
C-----------------------------------------------------------------------
      SUBROUTINE Write_Output_3D (Time,IT)
      USE SIZES, ONLY : SZ, INPUTDIR, NBYTE
      USE GLOBAL, ONLY : ScreenUnit
      USE GLOBAL_3DVS 
      USE GLOBAL_IO, ONLY: storeOne, storeTwo, HEADER_MAX,
     &    open_gbl_file, write_gbl_file, 
     &    write_gbl_file_skip_default, collectFullDomainArray,
     &    packOne, unpackOne, packTwo, unpackTwo, pack63, pack83 
#ifdef NETCDF
      USE NETCDFIO, ONLY :
     &     timenc, timenc_len, fileunit, irtype, 
     &     xenc, yenc, xvnc, yvnc, xpnc, ypnc, xnc, ync, xwnc, ywnc,
     &     et00nc, UU00nc, VV00nc, RMP00nc, RMU00nc, RMV00nc, nbvvnc,
     &     eta22, uu22, vv22, pr22, WNCX22, WNCY22, 
     &     nbounc,nvelnc,nvellnc,
     &     nopenc,netanc,nvdllnc
!     &     writenc41, writenc42, writenc43, writenc44,
!     &     writenc45, writenc46
       USE NodalAttributes, ONLY : 
     &     nolibf, nwp, tau0, cf, eslm
#endif


      

      IMPLICIT NONE 
      INTEGER, intent(in) :: IT
      REAL(8), intent(in) :: Time

      INTEGER NC1, NC2, NC3, NCEle
      INTEGER NM1, NM2, NM3
      INTEGER ::     NN, NH, NEle, N1, N2, N3, k
      REAL(SZ) C1, C2, C3
      REAL(SZ) EE1, EE2, EE3
      REAL(SZ) H2, H2N1, H2N2, H2N3
      REAL(SZ) U11, U22, U33
      REAL(SZ) V11, V22, V33

!      type(OutputDataDescript_t) :: descript
!      type(OutputDataDescript_t), SAVE :: ElevStaDescript
!      type(OutputDataDescript_t), SAVE :: VelStaDescript
!      type(OutputDataDescript_t), SAVE :: ElevDescript
!      type(OutputDataDescript_t), SAVE :: Tau0Descript
!      type(OutputDataDescript_t), SAVE :: VelDescript
!      type(OutputDataDescript_t), SAVE :: PrStaDescript
!      type(OutputDataDescript_t), SAVE :: WindVelStaDescript
!      type(OutputDataDescript_t), SAVE :: PrDescript
!      type(OutputDataDescript_t), SAVE :: WindVelDescript
!      type(OutputDataDescript_t), SAVE :: ConcStaDescript
!      type(OutputDataDescript_t), SAVE :: ConcDescript

C     jgf48.03 Full domain arrays used when writing globalio output.
!      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: ET00_g(:)
!      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: UU00_g(:)
!      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: VV00_g(:)
!      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: ETA2_g(:)
!      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: TAU0VAR_g(:)
!      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: UU2_g(:)     
!      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: VV2_g(:)
!      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: RMP00_g(:)     
!      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: RMU00_g(:)     
!      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: RMV00_g(:)
!      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: Pr2_g(:)     
!      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: WVNXOut_g(:)     
!      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: WVNYOut_g(:)
!      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: CC00_g(:)
!      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: CH1_g(:)

!     LOGICAL, SAVE :: FirstCall = .true.
      
!      IF (FirstCall) THEN

!         ! fort.61
!         write(16,*) 'Allocating ElevStaDescript' !jgfdebug48.03
!         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
!            ALLOCATE(ET00_g(NSTAE_G)) 
!         ENDIF
!         ElevStaDescript % specifier            = NOUTE
!         ElevStaDescript % initial_value        = 0.0
!         ElevStaDescript % num_items_per_record = 1
!         ElevStaDescript % num_fd_records       = NSTAE_G
!         ElevStaDescript % num_records_this     = NSTAE
!         ElevStaDescript % imap                 => IMAP_STAE_LG
!         ElevStaDescript % array                => ET00
!         ElevStaDescript % array_g              => ET00_g

         ! fort.62
!         write(16,*) 'Allocating VelStaDescript' !jgfdebug48.03
!         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
!            ALLOCATE(UU00_g(NSTAV_G))
!            ALLOCATE(VV00_g(NSTAV_G)) 
!         ENDIF
!         VelStaDescript % specifier            = NOUTV
!         VelStaDescript % initial_value        = 0.0
!         VelStaDescript % num_items_per_record = 2
!         VelStaDescript % num_fd_records       = NSTAV_G
!         VelStaDescript % num_records_this     = NSTAV
!         VelStaDescript % imap                 => IMAP_STAV_LG
!         VelStaDescript % array                => UU00
!         VelStaDescript % array2               => VV00
!         VelStaDescript % array_g              => UU00_g
!         VelStaDescript % array2_g             => VV00_g
!
         ! fort.63
!         write(16,*) 'Allocating ElevDescript' !jgfdebug48.03
!         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
!            ALLOCATE(ETA2_g(NP_G))
!         ENDIF
!         ElevDescript % specifier            = NOUTGE
!         ElevDescript % initial_value        = 0.0
!         ElevDescript % num_items_per_record = 1
!         ElevDescript % num_fd_records       = NP_G
!         ElevDescript % num_records_this     = NP
!         ElevDescript % imap                 => NODES_LG
!         ElevDescript % array                => ETA2
!         ElevDescript % array_g              => ETA2_g

         ! tau0
!        write(16,*) 'Allocating Tau0Descript' !jgfdebug48.03
!         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
!            ALLOCATE(TAU0VAR_g(NP_G))
!         ENDIF
!         Tau0Descript % specifier            = 1  ! hard-code to ascii text
!         Tau0Descript % initial_value        = 0.0
!         Tau0Descript % num_items_per_record = 1
!         Tau0Descript % num_fd_records       = NP_G
!         Tau0Descript % num_records_this     = NP
!         Tau0Descript % imap                 => NODES_LG
!         Tau0Descript % array                => TAU0VAR
!         Tau0Descript % array_g              => TAU0VAR_g

         ! fort.64
!         write(16,*) 'Allocating VelDescript' !jgfdebug48.03
!         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
!            ALLOCATE(UU2_g(NP_G))
!            ALLOCATE(VV2_g(NP_G)) 
!         ENDIF
!         VelDescript % specifier            = NOUTGV
!         VelDescript % initial_value        = 0.0
!         VelDescript % num_items_per_record = 2
!         VelDescript % num_fd_records       = NP_G
!         VelDescript % num_records_this     = NP
!         VelDescript % imap                 => NODES_LG
!         VelDescript % array                => UU2
!         VelDescript % array2               => VV2
!         VelDescript % array_g              => UU2_g
!         VelDescript % array2_g             => VV2_g
!
1100  FORMAT(1X,E16.10,1X,I10,32000(2X,E12.6))
1104  FORMAT(9X,I6,4X,32000(E12.6,2X))

C
C.... 3D Density, Temperature, Salinity Station Output (Unit 41)
C
      IF(I3DSD.NE.0) THEN 

C       Check to see if it is time to generate 3D station density output.

        IF((IT.GT.NTO3DSDS).AND.(IT.LE.NTO3DSDF)) N3DSD=N3DSD+1
        IF(N3DSD.EQ.NSpo3DSD) THEN                                           !if yes

C         Write time header into file

          IF(ABS(I3DSD).EQ.1) THEN                                          !ASCII FORMAT
            OPEN(41,FILE=TRIM(LOCALDIR)//'/'//'fort.41',
     &                            ACCESS='SEQUENTIAL',POSITION='APPEND')
            IF(IDen.EQ.1) WRITE(41,1100) Time,IT,(Sigma(k),k=1,NFEN-1)                           
            IF(IDen.EQ.2) WRITE(41,1100) Time,IT,(Sigma(k),Sigma(k),
     &                               k=1,NFEN-1),Sigma(NFEN)            
            IF(IDen.EQ.3) WRITE(41,1100) Time,IT,(Sigma(k),Sigma(k),
     &                               k=1,NFEN-1),Sigma(NFEN)             
            IF(IDen.EQ.4) WRITE(41,1100) Time,IT,(Sigma(k),Sigma(k),
     &                      Sigma(k),k=1,NFEN-1),Sigma(NFEN),Sigma(NFEN)
            I3DSDRec=I3DSDRec+1
          ELSEIF(ABS(I3DSD).EQ.2) THEN                                      !BINARY FORMAT
            OPEN(41,FILE=TRIM(LOCALDIR)//'/'//'fort.41',ACCESS='DIRECT',
     &            RECL=NByte)
            WRITE(41,REC=I3DSDRec+1) Time
            WRITE(41,REC=I3DSDRec+2) IT
            I3DSDRec=I3DSDRec+2
          ENDIF

C         Process each station

          DO NN=1,NSta3DD

C           Interpolate to station locations

            NEle=NE3DD(NN)
            N1=NM(NEle,1)
            N2=NM(NEle,2)
            N3=NM(NEle,3)		               
            DO k=1,NFEN
              SigTSta(k)=SigT(N1,k)*StaI3DD1(NN)+SigT(N2,k)*StaI3DD2(NN)
     &                                          +SigT(N3,k)*StaI3DD3(NN)
              IF((ABS(IDen).EQ.2).OR.(ABS(IDen).EQ.4)) SalSta(k) =
     &            Sal(N1,k) *StaI3DD1(NN)+Sal(N2,k) *StaI3DD2(NN)
     &           +Sal(N3,k) *StaI3DD3(NN)
              IF((ABS(IDen).EQ.3).OR.(ABS(IDen).EQ.4)) TempSta(k)=
     &            Temp(N1,k)*StaI3DD1(NN)+Temp(N2,k)*StaI3DD2(NN)
     &           +Temp(N3,k)*StaI3DD3(NN)
            END DO

C           Write station output

            IF(ABS(I3DSD).EQ.1) THEN                                        !ASCII FORMAT
              IF(ABS(IDen).EQ.1) WRITE(41,1104) NN,(SigTSta(k),k=1,NFEN)  
              IF(ABS(IDen).EQ.2) WRITE(41,1104) NN,(SigTSta(k),
     &                                               SalSta(k),k=1,NFEN)			
              IF(ABS(IDen).EQ.3) WRITE(41,1104) NN,(SigTSta(k),
     &                                              TempSta(k),k=1,NFEN)
              IF(ABS(IDen).EQ.4) WRITE(41,1104) NN,(SigTSta(k),
     &                                    SalSta(k),TempSta(k),k=1,NFEN)	
              I3DSDRec=I3DSDRec+1
            ELSEIF(ABS(I3DSD).EQ.2) THEN                                    !BINARY FORMAT
              WRITE(41,REC=I3DSDRec+1) NN
              I3DSDRec=I3DSDRec+1
              DO k=1,NFEN            
                IF(ABS(IDen).EQ.1) THEN
                  WRITE(41,REC=I3DSDRec+1) SigTSta(k)
                  I3DSDRec=I3DSDRec+1
                ELSEIF(ABS(IDen).EQ.2) THEN
                  WRITE(41,REC=I3DSDRec+1) SigTSta(k)
                  WRITE(41,REC=I3DSDRec+2) SalSta(k)
                  I3DSDRec=I3DSDRec+2
                ELSEIF(ABS(IDen).EQ.3) THEN
                  WRITE(41,REC=I3DSDRec+1) SigTSta(k)
                  WRITE(41,REC=I3DSDRec+2) TempSta(k)		
                  I3DSDRec=I3DSDRec+2		
                ELSEIF(ABS(IDen).EQ.4) THEN
                  WRITE(41,REC=I3DSDRec+1) SigTSta(k)
                  WRITE(41,REC=I3DSDRec+2) SalSta(k)
                  WRITE(41,REC=I3DSDRec+3) TempSta(k)
                  I3DSDRec=I3DSDRec+3
                ENDIF
              ENDDO
            ENDIF
          ENDDO
          N3DSD=0
          CLOSE(41)
        ENDIF
      ENDIF

C
C.... 3D Velocity Station Output (Unit 42)
C
      IF(I3DSV.NE.0) THEN

C       Check to see if it is time to generate 3D station velocity output.

        IF((IT.GT.NTO3DSVS).AND.(IT.LE.NTO3DSVF)) N3DSV=N3DSV+1
        IF(N3DSV.EQ.NSpo3DSV) THEN                                          !if yes

C         Write time header into file

          IF(ABS(I3DSV).EQ.1) THEN                                          !ASCII FORMAT
            OPEN(42,FILE=TRIM(LOCALDIR)//'/'//'fort.42',
     &                            ACCESS='SEQUENTIAL',POSITION='APPEND')
            WRITE(42,1100) Time,IT,(Sigma(k),Sigma(k),Sigma(k),
     &                       k=1,NFEN-1),Sigma(NFEN),Sigma(NFEN)
            I3DSVRec=I3DSVRec+1
          ELSEIF(ABS(I3DSV).EQ.2) THEN                                      !BINARY FORMAT
            OPEN(42,FILE=TRIM(LOCALDIR)//'/'//'fort.42',ACCESS='DIRECT',
     &            RECL=NByte)
            WRITE(42,REC=I3DSVRec+1) Time
            WRITE(42,REC=I3DSVRec+2) IT
            I3DSVRec=I3DSVRec+2
          ENDIF

C         Process each station

          DO NN=1,NSta3DV

C           Interpolate to station locations

            NEle=NE3DV(NN)
            N1=NM(NEle,1)
            N2=NM(NEle,2)
            N3=NM(NEle,3)
            DO k=1,NFEN
              qSta(k) =q(N1,k) *StaI3DV1(NN)+q(N2,k) *StaI3DV2(NN)
     &                                      +q(N3,k) *StaI3DV3(NN)
              WZSta(k)=WZ(N1,k)*StaI3DV1(NN)+WZ(N2,k)*StaI3DV2(NN)
     &                                      +WZ(N3,k)*StaI3DV3(NN)
            END DO

C           Write station output
          
            IF(ABS(I3DSV).EQ.1) THEN                                        !ASCII FORMAT
              WRITE(42,1104) NN,(REAL(qSta(k)),AIMAG(qSta(k)),
     &                                                WZSta(k),k=1,NFEN)
              I3DSVRec=I3DSVRec+1
            ELSEIF(ABS(I3DSV).EQ.2) THEN                                    !BINARY FORMAT
              WRITE(42,REC=I3DSVRec+1) NN
              I3DSVRec=I3DSVRec+1
              DO k=1,NFEN
                WRITE(42,REC=I3DSVRec+1) REAL(qSta(k))
                WRITE(42,REC=I3DSVRec+2) AIMAG(qSta(k))
                WRITE(42,REC=I3DSVRec+3) WZSta(k)
                I3DSVRec=I3DSVRec+3
              END DO
            ENDIF
          ENDDO
          N3DSV=0
          CLOSE(42)
        ENDIF
      ENDIF

C
C.... 3D Turbulence Station Output (Unit 43)
C
      IF(I3DST.NE.0) THEN

C       Check to see if it is time to generate 3D station turbulence output.

        IF((IT.GT.NTO3DSTS).AND.(IT.LE.NTO3DSTF)) N3DST=N3DST+1
        IF(N3DST.EQ.NSpo3DST) THEN                                          !if yes

C         Write time header into file

          IF(ABS(I3DST).EQ.1) THEN                                          !ASCII FORMAT
            OPEN(43,FILE=TRIM(LOCALDIR)//'/'//'fort.43',
     &                            ACCESS='SEQUENTIAL',POSITION='APPEND')        
            WRITE(43,1100) Time,IT,(Sigma(k),Sigma(k),Sigma(k),
     &                       k=1,NFEN-1),Sigma(NFEN),Sigma(NFEN)
            I3DSTRec=I3DSTRec+1
          ELSEIF(ABS(I3DST).EQ.2) THEN                                      !BINARY FORMAT
            OPEN(43,FILE=TRIM(LOCALDIR)//'/'//'fort.43',ACCESS='DIRECT',
     &            RECL=NByte)
            WRITE(43,REC=I3DSTRec+1) Time
            WRITE(43,REC=I3DSTRec+2) IT
            I3DSTRec=I3DSTRec+2
          ENDIF

C         Process each station

          DO NN=1,NSta3DT     

C           Interpolate to station locations

            NEle=NE3DT(NN)
            N1=NM(NEle,1)
            N2=NM(NEle,2)
            N3=NM(NEle,3)
            DO k=1,NFEN
              q20Sta(k)=q20(N1,k)*StaI3DT1(NN)+q20(N2,k)*StaI3DT2(NN)
     &                                        +q20(N3,k)*StaI3DT3(NN)
              lSta(k)  =l(N1,k)  *StaI3DT1(NN)+l(N2,k)  *StaI3DT2(NN)
     &                                        +l(N3,k)  *StaI3DT3(NN)
              EVSta(k) =EV(N1,k) *StaI3DT1(NN)+EV(N2,k) *StaI3DT2(NN)
     &                                        +EV(N3,k) *StaI3DT3(NN)
            END DO

C           Write station output

            IF(ABS(I3DST).EQ.1) THEN                                        !ASCII FORMAT
              WRITE(43,1104) NN,(q20Sta(k),lSta(k),EVSta(k),k=1,NFEN)
              I3DSTRec=I3DSTRec+1
            ELSEIF(ABS(I3DST).EQ.2) THEN                                    !BINARY FORMAT
              WRITE(43,REC=I3DSTRec+1) NN
              I3DSTRec=I3DSTRec+1
              DO k=1,NFEN
                WRITE(43,REC=I3DSTRec+1) q20Sta(k)
                WRITE(43,REC=I3DSTRec+2) lSta(k)
                WRITE(43,REC=I3DSTRec+3) EVSta(k)
                I3DSTRec=I3DSTRec+3
              END DO
            ENDIF
          ENDDO        
          N3DST=0
          CLOSE(43)
	  ENDIF
      ENDIF

C
C.... 3D Density, Temperature, Salinity Global Output (Unit 44)
C
      IF(I3DGD.NE.0) THEN

C     Check to see if it is time to generate 3D global density output.

        IF((IT.GT.NTO3DGDS).AND.(IT.LE.NTO3DGDF)) N3DGD=N3DGD+1
        IF(N3DGD.EQ.NSpo3DGD) THEN                                           !if yes

C         Write time header into file

          IF(ABS(I3DGD).EQ.1) THEN                                          !ASCII FORMAT
            OPEN(44,FILE=TRIM(LOCALDIR)//'/'//'fort.44',
     &                            ACCESS='SEQUENTIAL',POSITION='APPEND')
            IF(IDen.EQ.1) WRITE(44,1100) Time,IT,(Sigma(k),k=1,NFEN-1)                           
            IF(IDen.EQ.2) WRITE(44,1100) Time,IT,(Sigma(k),Sigma(k),
     &                               k=1,NFEN-1),Sigma(NFEN)            
            IF(IDen.EQ.3) WRITE(44,1100) Time,IT,(Sigma(k),Sigma(k),
     &                               k=1,NFEN-1),Sigma(NFEN)             
            IF(IDen.EQ.4) WRITE(44,1100) Time,IT,(Sigma(k),Sigma(k),
     &                      Sigma(k),k=1,NFEN-1),Sigma(NFEN),Sigma(NFEN)
            I3DGDRec=I3DGDRec+1
          ELSEIF(ABS(I3DGD).EQ.2) THEN                                      !BINARY FORMAT
            OPEN(44,FILE=TRIM(LOCALDIR)//'/'//'fort.44',ACCESS='DIRECT',
     &            RECL=NByte)
            WRITE(44,REC=I3DGDRec+1) Time
            WRITE(44,REC=I3DGDRec+2) IT
            I3DGDRec=I3DGDRec+2
          ENDIF

C         Write out at each node

          DO NH=1,NP

            IF(ABS(IDen).EQ.1) THEN
              IF(ABS(I3DGD).EQ.1) THEN
	          WRITE(44,1104) NH,(SigT(NH,k),k=1,NFEN)                      !ASCII
                I3DGDRec=I3DGDRec+1
              ELSEIF(ABS(I3DGD).EQ.2) THEN                                   !BINARY
                WRITE(44,REC=I3DGDRec+1) NH
                I3DGDRec=I3DGDRec+1
                DO k=1,NFEN
                  WRITE(44,REC=I3DGDRec+1) SigT(NH,k)
                  I3DGDRec=I3DGDRec+1
                ENDDO
              ENDIF
            ENDIF

            IF(ABS(IDen).EQ.2) THEN
              IF(ABS(I3DGD).EQ.1) THEN
	          WRITE(44,1104) NH,(SigT(NH,k),Sal(NH,k),k=1,NFEN)            !ASCII
                I3DGDRec=I3DGDRec+1
              ELSEIF(ABS(I3DGD).EQ.2) THEN                                   !BINARY
                WRITE(44,REC=I3DGDRec+1) NH
                I3DGDRec=I3DGDRec+1
                DO k=1,NFEN
                  WRITE(44,REC=I3DGDRec+1) SigT(NH,k)
                  WRITE(44,REC=I3DGDRec+2) Sal(NH,k)
                  I3DGDRec=I3DGDRec+2
                END DO
              ENDIF
            ENDIF

            IF(ABS(IDen).EQ.3) THEN
              IF(ABS(I3DGD).EQ.1) THEN
                WRITE(44,1104) NH,(SigT(NH,k),Temp(NH,k),k=1,NFEN)           !ASCII
                I3DGDRec=I3DGDRec+1
              ELSEIF(ABS(I3DGD).EQ.2) THEN                                   !BINARY
                WRITE(44,REC=I3DGDRec+1) NH
                I3DGDRec=I3DGDRec+1
                DO k=1,NFEN
                  WRITE(44,REC=I3DGDRec+1) SigT(NH,k)
                  WRITE(44,REC=I3DGDRec+2) Temp(NH,k)
                  I3DGDRec=I3DGDRec+2
                ENDDO
              ENDIF
            ENDIF

            IF(ABS(IDen).EQ.4) THEN
              IF(ABS(I3DGD).EQ.1) THEN
                WRITE(44,1104) NH,(SigT(NH,k),Temp(NH,k),Sal(NH,k)
     &                                                        ,k=1,NFEN)     !ASCII
                I3DGDRec=I3DGDRec+1
              ELSEIF(ABS(I3DGD).EQ.2) THEN                                   !BINARY
                WRITE(44,REC=I3DGDRec+1) NH
                I3DGDRec=I3DGDRec+1
                DO k=1,NFEN
                  WRITE(44,REC=I3DGDRec+1) SigT(NH,k)
                  WRITE(44,REC=I3DGDRec+2) Temp(NH,k)
                  WRITE(44,REC=I3DGDRec+3) Sal(NH,k)
                  I3DGDRec=I3DGDRec+3
                ENDDO
              ENDIF
            ENDIF
          ENDDO
          N3DGD=0
          CLOSE(44)
        ENDIF
      ENDIF

C
C.... 3D Velocity Global Output (Unit 45)
C
      IF(I3DGV.GT.0) THEN

C     Check to see if it is time to generate 3D global velocity output.
         
         IF((IT.GT.NTO3DGVS).AND.(IT.LE.NTO3DGVF)) N3DGV=N3DGV+1

         IF(N3DGV.EQ.NSpo3DGV) THEN !if yes

C         Write time header into file

            IF(ABS(I3DGV).EQ.1) THEN !ASCII FORMAT
               OPEN(45,FILE=TRIM(LOCALDIR)//'/'//'fort.45',
     &              ACCESS='SEQUENTIAL',POSITION='APPEND')
               WRITE(45,1100) Time,IT,(Sigma(k),Sigma(k),Sigma(k),
     &              k=1,NFEN-1),Sigma(NFEN),Sigma(NFEN)
               I3DGVRec=I3DGVRec+1
            ELSEIF(ABS(I3DGV).EQ.2) THEN !BINARY FORMAT
               OPEN(45,FILE=TRIM(LOCALDIR)//'/'//'fort.45',
     &              ACCESS='DIRECT',RECL=NByte)
               WRITE(45,REC=I3DGVRec+1) Time
               WRITE(45,REC=I3DGVRec+2) IT
               I3DGVRec=I3DGVRec+2
            ENDIF
            
C     Write out at each node

            DO NH=1,NP
               IF(ABS(I3DGV).EQ.1) THEN
                  WRITE(45,1104) NH,(REAL(q(NH,k)),AIMAG(q(NH,k)),
     &                 WZ(NH,k),k=1,NFEN) !ASCII
                  I3DGVRec=I3DGVRec+1
               ELSEIF(ABS(I3DGV).EQ.2) THEN !BINARY
                  WRITE(45,REC=I3DGVRec+1) NH
                  I3DGVRec=I3DGVRec+1
                  DO k=1,NFEN
                     WRITE(45,REC=I3DGVRec+1) REAL(q(NH,k))
                     WRITE(45,REC=I3DGVRec+2) AIMAG(q(NH,k))
                     WRITE(45,REC=I3DGVRec+3) WZ(NH,k)
                     I3DGVRec=I3DGVRec+3
                  END DO
               ENDIF
            ENDDO
            N3DGV=0
            CLOSE(45)
         ENDIF
      ENDIF
	
C
C.... 3D Turbulence Global Output (Unit 46)
C
      IF(I3DGT.GT.0) THEN

C     Check to see if it is time to generate 3D global turbulence output.

        IF((IT.GT.NTO3DGTS).AND.(IT.LE.NTO3DGTF)) N3DGT=N3DGT+1
        IF(N3DGT.EQ.NSpo3DGT) THEN                                           !if yes

C         Write time header into file

          IF(ABS(I3DGT).EQ.1) THEN                                          !ASCII FORMAT
            OPEN(46,FILE=TRIM(LOCALDIR)//'/'//'fort.46',
     &                            ACCESS='SEQUENTIAL',POSITION='APPEND')        
            WRITE(46,1100) Time,IT,(Sigma(k),Sigma(k),Sigma(k),
     &                       k=1,NFEN-1),Sigma(NFEN),Sigma(NFEN)
            I3DGTRec=I3DGTRec+1
          ELSEIF(ABS(I3DGT).EQ.2) THEN                                      !BINARY FORMAT
            OPEN(46,FILE=TRIM(LOCALDIR)//'/'//'fort.46',ACCESS='DIRECT',
     &            RECL=NByte)
            WRITE(46,REC=I3DGTRec+1) Time
            WRITE(46,REC=I3DGTRec+2) IT
            I3DGTRec=I3DGTRec+2
          ENDIF

C         Write out at each node

          DO NH=1,NP
            IF(ABS(I3DGT).EQ.1) THEN
              WRITE(46,1104) NH,(q20(NH,k),l(NH,k),EV(NH,k),k=1,NFEN)        !ASCII
              I3DGTRec=I3DGTRec+1
            ELSEIF(ABS(I3DGT).EQ.2) THEN                                     !BINARY
              WRITE(46,REC=I3DGTRec+1) NH
              I3DGTRec=I3DGTRec+1
              DO k=1,NFEN
                WRITE(46,REC=I3DGTRec+1) q20(NH,k)
                WRITE(46,REC=I3DGTRec+2) l(NH,k)
                WRITE(46,REC=I3DGTRec+3) EV(NH,k)
                I3DGTREC=I3DGTREC+3
              ENDDO
            ENDIF
          ENDDO
          N3DGT=0
          CLOSE(46)
        ENDIF
      ENDIF

 
      RETURN
      END

C-----------------------------------------------------------------------
C     S U B R O U T I N E   W R I T E  O U T  A R R A Y  
C-----------------------------------------------------------------------
C     jgf48.03 This subroutine was created to write out a column 
C     vector (i.e., nodal data such as water surface elevation or 
C     pressure) to a file.
C-----------------------------------------------------------------------
      SUBROUTINE writeOutArray(lun, time, it, descript, pack_cmd,
     &                         unpack_cmd, filepos)

      USE SIZES
      USE GLOBAL
      USE GLOBAL_IO, ONLY : collectFullDomainArray
      IMPLICIT NONE
C     args
      INTEGER, intent(in) :: lun ! logical unit number of file to write to
      REAL(8), intent(in) :: time ! seconds since cold start
      INTEGER, intent(in) :: it   ! number of time steps since cold start
      type(OutputDataDescript_t), intent(in) :: descript !describes output data
      EXTERNAL :: pack_cmd   ! subroutine used to pack data on subdomain
      EXTERNAL :: unpack_cmd ! subroutine used to unpack data on proc 0
      INTEGER, intent(inout) :: filepos  ! current position in the output file
C     local vars
      CHARACTER(7) :: fn     ! outfile name (valid for lun between 10 and 99)
      INTEGER :: I           ! loop counter

C     initialize output file name
      fn(1:5) = 'fort.'
      WRITE(fn(6:7),2) lun

C     collect up the data from subdomains if running in parallel
      IF ( MNPROC.gt.1 ) THEN
         write(16,*) 'About to collectFullDomainArray' !jgfdebug48.03
         CALL collectFullDomainArray(descript, pack_cmd, unpack_cmd)
      ENDIF

C     write data according to format specifier from fort.15 (e.g., NOUTE)
      SELECT CASE (ABS(descript % specifier))
 
      CASE(1) !ascii text

         write(16,*) 'About to open globalio text file.' !jgfdebug48.03
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            OPEN(lun,FILE=TRIM(GLOBALDIR)//'/'//fn,
     &           ACCESS='SEQUENTIAL',POSITION='APPEND')
            WRITE(lun,2120) time,IT
            IF (descript % num_items_per_record .eq. 1) THEN
               DO I=1, descript % num_fd_records
                  WRITE(16,*) I ! jgfdebug48.03
                  WRITE(lun,2453) I, descript % array_g(I)
               ENDDO
            ENDIF
            IF (descript % num_items_per_record .eq. 2) THEN
               DO I=1, descript % num_fd_records
                  WRITE(16,*) I ! jgfdebug48.03
                  WRITE(lun,2454) I, descript % array_g(I),
     &                               descript % array2_g(I)  
               ENDDO
            ENDIF
            WRITE(16,*) 'About to close globalio text file.'!jgfdebug48.03
            CLOSE(lun)
         ENDIF

         IF ((MNPROC.eq.1).or.(WRITE_LOCAL_FILES)) THEN
            OPEN(lun,FILE=TRIM(LOCALDIR)//'/'//fn,
     &           ACCESS='SEQUENTIAL',POSITION='APPEND')
            WRITE(lun,2120) time,IT
            IF (descript % num_items_per_record .eq. 1) THEN           
               DO I=1, descript % num_records_this
                  WRITE(lun,2453) I, descript % array(I)
               END DO
            ENDIF
            IF (descript % num_items_per_record .eq. 2) THEN 
               DO I=1, descript % num_records_this
                  WRITE(lun,2454) I, descript % array(I),
     &                               descript % array2(I)
               END DO
            ENDIF
 
         ENDIF
         filepos = filepos+1+descript % num_records_this

      CASE(2) !binary (nonportable)

         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            OPEN(lun,FILE=TRIM(GLOBALDIR)//'/'//fn,
     &           ACCESS='DIRECT',RECL=NBYTE)
            WRITE(lun,REC=filepos+1) time
            WRITE(lun,REC=filepos+2) IT
            filepos = filepos + 2
            IF ( descript % num_items_per_record .eq. 1 ) THEN
               DO I=1, descript % num_fd_records
                  WRITE(lun,REC=filepos+I) descript % array_g(I)
               END DO
            ENDIF
            IF ( descript % num_items_per_record .eq. 2 ) THEN
               DO I=1, descript % num_fd_records
                  WRITE(lun,REC=filepos+2*I-1) descript % array_g(I)
                  WRITE(lun,REC=filepos+2*I)   descript % array2_g(I)
               END DO
            ENDIF
            CLOSE(lun)
         ENDIF

         IF ((MNPROC.eq.1).or.(WRITE_LOCAL_FILES)) THEN
            OPEN(lun,FILE=TRIM(LOCALDIR)//'/'//fn,
     &           ACCESS='DIRECT',RECL=NBYTE)
            WRITE(lun,REC=filepos+1) time
            WRITE(lun,REC=filepos+2) IT
            filepos = filepos + 2
            IF ( descript % num_items_per_record .eq. 1 ) THEN
               DO I=1, descript % num_records_this
                  WRITE(lun,REC=filepos+I) descript % array(I)
               END DO
            ENDIF
            IF ( descript % num_items_per_record .eq. 2 ) THEN
               DO I=1, descript % num_records_this
                  WRITE(lun,REC=filepos+2*I-1) descript % array_g(I)
                  WRITE(lun,REC=filepos+2*I)   descript % array2_g(I)
               END DO
            ENDIF
            CLOSE(lun)
         ENDIF
         filepos = filepos + descript % num_records_this

      CASE(3) !netcdf (portable)
#ifdef NETCDF
!       IF(MYPROC.EQ.0) PRINT *,"W: BEFORE WRITING OUT ARRAY "
         CALL writeOutArrayNetCDF(lun, time, it, descript)
!       IF(MYPROC.EQ.0) PRINT *,"W: AFTER WRITING OUT ARRAY "
#else
         WRITE(ScreenUnit,*) 'ERROR: NetCDF is not available.'
         WRITE(16,*) 'ERROR: NetCDF is not available.'
#endif

      CASE DEFAULT
         WRITE(ScreenUnit,*) 'ERROR: Invalid output specifier.' 
      END SELECT

 2    FORMAT(I2)
 2120 FORMAT(2X,E20.10,5X,I10)
 2453 FORMAT(2x, i8, 2x, E20.10, E20.10, E20.10, E20.10) 
 2454 FORMAT(2X,I8,2(2X,E15.8))
C-----------------------------------------------------------------------
       END SUBROUTINE writeOutArray
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E   W R I T E  O U T  A R R A Y  N E T C D F 
C-----------------------------------------------------------------------
C     jgf48.03 This subroutine was created from cf's code in timestep.F to
C     write output files in NetCDF format.
C-----------------------------------------------------------------------
      SUBROUTINE writeOutArrayNetCDF(lun, time, it, descript)
      USE SIZES
      USE GLOBAL
#ifdef NETCDF
      USE NETCDFIO, ONLY :
     &     timenc, timenc_len, fileunit, irtype, 
     &     xenc, yenc, xvnc, yvnc, xpnc, ypnc, xnc, ync, xwnc, ywnc,
     &     et00nc, UU00nc, VV00nc, RMP00nc, RMU00nc, RMV00nc, nbvvnc,
     &     eta22, uu22, vv22, pr22, WNCX22, WNCY22, 
     &     nbounc,nvelnc,nvellnc,
     &     nopenc,netanc,nvdllnc,
     &     writenc61, writenc62, writenc63, writenc64,
     &     writenc71, writenc72, writenc73, writenc74 
       USE NodalAttributes, ONLY : 
     &     nolibf, nwp, tau0, cf, eslm
#endif
      IMPLICIT NONE
C     args
      INTEGER, intent(in) :: lun ! logical unit number of file to write to
      REAL(8), intent(in) :: time ! seconds since cold start
      INTEGER, intent(in) :: it   ! number of time steps since cold start
      type(OutputDataDescript_t), intent(in) :: descript !describes output data
C     local vars
      INTEGER :: I,J,K           ! loop counter
      INTEGER :: num_records ! to save in the netcdf file
      INTEGER ncmode
!      INTEGER, ALLOCATABLE :: ibtypenc(:)

#ifdef NETCDF
!            WRITE(6,*)"================================================"
!            WRITE(6,*)"NOPE,NETA,NBOU,NVEL",NOPE,NETA,NBOU,NVEL

!            WRITE(6,*)"WRITE INPUT NBVV(K,I) ***** ===", NBVV

      IF (MNPROC.gt.1) ncmode=2
      IF (MNPROC.eq.1) ncmode=1
      IF (MNPROC.gt.1) num_records = descript % num_fd_records   ! parallel
      IF (MNPROC.eq.1) num_records = descript % num_records_this ! serial
      IF (MyProc.eq.0) THEN
         fileunit = lun
C     ==================================================================
C     START NETCDF I/O OPTION - MCF
C
!         WRITE(6,*)"fileunit = ", fileunit
         timenc_len=1
!         WRITE(6,*)"timenc_len ", timenc_len
!         WRITE(6,*)"time", time

         IF(.NOT. ALLOCATED(timenc)) 
     &        ALLOCATE(timenc(timenc_len))
         timenc(timenc_len)=time
!         WRITE(6,*)"timenc(timenc_len) =  ", timenc(timenc_len)
C
C     END NETCDF I/O OPTION - MCF
C     ==================================================================

         SELECT CASE(lun)
 
         CASE(61) ! netcdf fort.61 file

C     ==================================================================
C     START NETCDF I/O OPTION - MCF
C
            IF(.NOT. ALLOCATED(et00nc))
     &           ALLOCATE(et00nc(num_records,timenc_len))
            IF(.NOT. ALLOCATED(xenc)) 
     &           ALLOCATE(xenc(num_records))
            IF(.NOT. ALLOCATED(yenc)) 
     &           ALLOCATE(yenc(num_records))
            IF(.NOT. ALLOCATED(nbvvnc)) 
     &           ALLOCATE(nbvvnc(nbou,nvel))

            DO I=1, num_records
               IF (MNPROC.gt.1) THEN      ! parallel
                  et00nc (i,timenc_len) = descript % array_g(I)
               ELSE IF (MNPROC.eq.1) THEN ! serial
                  et00nc (i,timenc_len) = descript % array(I) 
               ENDIF
!               WRITE(6,*) "LOOP =  1 ", et00nc (i,timenc_len)
            END DO
!            WRITE(6,*) "BEFORE ICS =2 "
            
            if (ncmode.eq.1) then
            IF(ICS.EQ.2) THEN
               DO I=1, num_records
                   xenc(i)=SLEL(I)*RAD2DEG
                   yenc(i)=SFEL(I)*RAD2DEG
!                   WRITE(6,*) "LOOP =  ", I, xenc(i),yenc(i)
               END DO
            ELSE
               DO I=1,num_records
!                  IF(MYPROC.EQ.0)PRINT *,"wo: num_records xel", 
!     .             xel(I),num_records
             	  xenc(i)=XEL(I)
                  yenc(i)=YEL(I)
!                  IF(MYPROC.EQ.0)PRINT *,"wa: num_records xel", 
!     .             xel(I),num_records, xenc(i), yenc(i)
               END DO
            ENDIF

               DO I=1,nbou
                  DO J=1,nvel
             	  nbvvnc(i,j)=nbvv(i,j)
                  END DO
               END DO
            endif   

            CALL writenc61 (descript % specifier, 
     &             fileunit, timenc_len, num_records, timenc,
     &             et00nc, nt, NTRSPE, rundes, runid, agrid, irtype,
     &             dtdp, ihot, ics, nolibf, nolifa,
     &             nolica, nolicat, nwp,
     &             ncor, ntip, nws, nramp, tau0,
     &             statim,reftim, rnday, dramp,
     &             a00, b00, c00, h0,
     &             slam0, sfea0, cf,eslm, cori,
     &             ntif, nbfr, ne, nhy, nm,
     &             xenc, yenc, dp,
     &             nope, neta, nbdv,
     &             nbou, nvel, nbvvnc,
     &             ncmode,
     &             title, institution, source, history,
     &             references, comments, host, convention, contact,
     &             base_date)
     
             NCIESTP=NCIESTP+num_records
C
C     END NETCDF I/O OPTION - MCF
C     ==================================================================
C              
         CASE(62) ! netcdf fort.62 file
C
C     ==================================================================
C     START NETCDF I/O OPTION - MCF
C
C
            IF(.NOT.ALLOCATED(UU00nc))
     &          ALLOCATE(UU00nc(num_records,timenc_len))
            IF(.NOT.ALLOCATED(VV00nc))
     &          ALLOCATE(VV00nc(num_records,timenc_len))
            IF(.NOT.ALLOCATED(xvnc)) ALLOCATE(xvnc(num_records))
            IF(.NOT.ALLOCATED(yvnc)) ALLOCATE(yvnc(num_records))
            IF(.NOT. ALLOCATED(nbvvnc)) ALLOCATE(nbvvnc(nbou,nvel))

            DO I=1, num_records
               IF (MNPROC.gt.1) THEN      ! parallel
                  UU00nc(i,timenc_len)= descript % array_g(i)
                  VV00nc(i,timenc_len)= descript % array2_g(i)
               ELSE IF (MNPROC.eq.1) THEN ! serial
                  UU00nc(i,timenc_len)= descript % array(i)
                  VV00nc(i,timenc_len)= descript % array2(i)
               ENDIF
            END DO
            if (ncmode.eq.1) then
            IF(ICS.EQ.2) THEN
               DO I=1, num_records
                  xvnc(i)=SLEV(I)*RAD2DEG
                  yvnc(i)=SFEV(I)*RAD2DEG
               END DO
            ELSE
               DO I=1,num_records
             	  xvnc(i)=XEV(I)
                  yvnc(i)=YEV(I)
               END DO
            ENDIF
            
               DO I=1,nbou
                  DO J=1,nvel
             	  nbvvnc(i,j)=nbvv(i,j)
                  END DO
               END DO
           endif

            CALL writenc62 (descript % specifier, fileunit, 
     &             timenc_len, num_records, timenc,
     &             nt, NTRSPV, UU00nc, VV00nc,
     &             rundes, runid, agrid, irtype,
     &             dtdp,ihot, ics, nolibf, nolifa,
     &             nolica, nolicat, nwp,
     &             ncor, ntip, nws, nramp, tau0,
     &             statim,reftim, rnday, dramp,
     &             a00, b00, c00, h0,
     &             slam0, sfea0, cf,eslm, cori,
     &             ntif, nbfr, ne, nhy, nm,
     &             xvnc, yvnc, dp,
     &             nope, neta, nbdv,
     &             nbou, nvel, nbvvnc,
     &             ncmode,
     &             title, institution, source, history,
     &             references, comments, host, convention, contact,
     &             base_date)
            NCIVSTP=NCIVSTP+num_records

C
C     END NETCDF I/O OPTION - MCF
C     ==================================================================
C
         CASE(63) ! netcdf fort.63 file
C
C     ==================================================================
C     START NETCDF I/O OPTION - MCF
C

            IF(.NOT. ALLOCATED(eta22)) 
     &               ALLOCATE (eta22(num_records,timenc_len))
            IF(.NOT. ALLOCATED(xnc)) ALLOCATE(xnc(num_records))
            IF(.NOT. ALLOCATED(ync)) ALLOCATE(ync(num_records))
            IF(.NOT. ALLOCATED(nbvvnc)) ALLOCATE(nbvvnc(nbou,nvel))
            IF(.NOT. ALLOCATED(nvellnc)) ALLOCATE(nvellnc(nbou))
            IF(.NOT. ALLOCATED(nvdllnc)) ALLOCATE(nvdllnc(nope))
            IF(.NOT. ALLOCATED(ibtypenc)) ALLOCATE(ibtypenc(nbou))

            DO I=1, num_records
               IF (MNPROC.gt.1) THEN
                  eta22(i,timenc_len) = descript % array_g(i)
               ELSE IF (MNPROC.eq.1) THEN
                  eta22(i,timenc_len) = descript % array(i)
               ENDIF
            END DO

            if (ncmode.eq.1) then
            IF(ICS.EQ.2) THEN
               DO I=1,num_records
                     xnc(i)=slam(i)*RAD2DEG
                     ync(i)=sfea(i)*RAD2DEG
               END DO
            ELSE
               DO I=1,num_records
             	     xnc(i)=X(i)
                     ync(i)=Y(i)
       	       END DO
            ENDIF
             endif
!
!            elevation specified boundary forcing segments
!
               netanc=0
               DO i=1,nope
                  nvdllnc(i)=nvdll(i)
                  netanc=netanc+nvdllnc(i)
!                 DO j=1,nvel
!             	     nbvvnc(i,j)=nbvv(i,j)
!                 END DO
              END DO

            IF(.NOT. ALLOCATED(nopenc)) ALLOCATE(nopenc(netanc))
               
               k=0
                  DO i=1,nope
                        DO j=1,nvdllnc(i)
                        k=k+1
             	        nopenc(k)=nbdv(i,j)
                     END DO
                  END DO
!
!              normal flow (discharge) specified boundary segments
!
               nvelnc=0
               DO i=1,nbou
                  nvellnc(i)=nvell(i)
                  nvelnc=nvelnc+nvellnc(i)
                  DO j=1,nvel
             	     nbvvnc(i,j)=nbvv(i,j)
                 END DO
               END DO

            IF(.NOT. ALLOCATED(nbounc)) ALLOCATE(nbounc(nvelnc))
               
               k=0
                  DO i=1,nbou
                        DO j=1,nvellnc(i)
                        k=k+1
             	        nbounc(k)=nbvv(i,j)
                     END DO
                  END DO
             
!            WRITE(6,*)"================================================"
!            WRITE(6,*)"num_records =====================",num_records
!           WRITE(6,*)"WRITE INPUT NBVV(K,I) ===", NBVV
!            WRITE(633,*)"**noutge, fileunit,timenc_len, np, timenc",
!     &             noutge, fileunit,timenc_len, np, timenc,"**eta22",
!     &             eta22, "** nt = ", nt, "** ndsetse = ", ndsetse,
!     &             "** rundes = ", rundes, runid, agrid,
!     &             "** irtype = ", irtype,
!     &             "** dtdp = ", dtdp,  "** ihot = ", ihot,
!     &             "** ics, nolibf, nolifa = ", ics, nolibf, nolifa,
!     &             "** nolica, nolicat, nwp = ",nolica, nolicat, nwp,
!     &             "** ncor, ntip, nws, nramp, tau0 = ", ncor, ntip,
!     &             nws, nramp, tau0,
!     &             "** statim,reftim, rnday, dramp = ",statim,
!     &             reftim, rnday, dramp,
!     &             "** a00, b00, c00, h0 = ",a00, b00, c00, h0,
!     &             "** slam0, sfea0, cf,eslm, cori = ", slam0, sfea0,
!     &             cf,eslm, cori,
!     &             "**  ntif, nbfr, ne, nhy, nm = ", ntif, nbfr, ne,
!     &             nhy, nm,
!     &             "**  x = ", xnc, "**  y = ", ync, "**  dp = ", dp,
!     &             "nope, neta, nbdv = ", nope, neta, nbdv,
!     &             "nbou, nvel, nbvv = ", nbou, nvel, nbvvnc, 
!     &              "nbounc", nbounc,"nvellnc", nvellnc,
!     &             "ibtypenc",ibtypenc,"netanc",netanc,"nopenc",nopenc,
!     &             "nvdllnc",nvdllnc,
!     &             title, institution, source, history,
!     &             references, comments, host, convention, contact,
!     &             base_date

!            CALL writenc63 (noutge, fileunit,timenc_len, np, timenc,
!            WRITE(6,*)"==================63========================="
!             NCIGEP=NCIGEP+num_records
            CALL writenc63 (noutge, fileunit,timenc_len, 
     &             num_records, timenc,
     &             eta22, nt, ndsetse, rundes, runid, agrid, irtype,
     &             dtdp, ihot, ics, nolibf, nolifa,
     &             nolica, nolicat, nwp,
     &             ncor, ntip, nws, nramp, tau0,
     &             statim,reftim, rnday, dramp,
     &             a00, b00, c00, h0,
     &             slam0, sfea0, cf,eslm, cori,
     &             ntif, nbfr, ne, nhy, nm,
     &             xnc, ync, dp,
     &             nope, neta, nbdv,nopenc,netanc,
     &             nvdllnc,
     &             nbou, nvel, nbvvnc,nbounc,nvelnc,
     &             nvellnc,ibtypenc,
     &             ncmode,
     &             title, institution, source, history,
     &             references, comments, host, convention, contact,
     &             base_date)
             NCIGEP=NCIGEP+num_records
!             WRITE(634,*)"**noutge, fileunit,timenc_len, np, timenc",
!     &             noutge, fileunit,timenc_len, np, timenc,"**eta22",
!     &             eta22, "** nt = ", nt, "** ndsetse = ", ndsetse,
!     &             "** rundes = ", rundes, runid, agrid,
!     &             "** irtype = ", irtype,
!     &             "** dtdp = ", dtdp,  "** ihot = ", ihot,
!     &             "** ics, nolibf, nolifa = ", ics, nolibf, nolifa,
!     &             "** nolica, nolicat, nwp = ",nolica, nolicat, nwp,
!     &             "** ncor, ntip, nws, nramp, tau0 = ", ncor, ntip,
!     &             nws, nramp, tau0,
!     &             "** statim,reftim, rnday, dramp = ",statim,
!     &             reftim, rnday, dramp,
!     &             "** a00, b00, c00, h0 = ",a00, b00, c00, h0,
!     &             "** slam0, sfea0, cf,eslm, cori = ", slam0, sfea0,
!     &             cf,eslm, cori,
!     &             "**  ntif, nbfr, ne, nhy, nm = ", ntif, nbfr, ne,
!     &             nhy, nm,
!     &             "**  x = ", xnc, "**  y = ", ync, "**  dp = ", dp,
!     &             "nope, neta, nbdv = ", nope, neta, nbdv,
!     &             "nbou, nvel, nbvv = ", nbou, nvel, nbvv,
!     &             title, institution, source, history,
!     &             references, comments, host, convention, contact,
!     &             base_date

!             IF((myproc.eq.0.and.ncmode.eq.2).or.ncmode.eq.1)
!     &       WRITE(6,*)"NCCOUGE = ",NCCOUGE," NCIGEP = ",NCIGEP

C
C     END NETCDF I/O OPTION - MCF
C     ==================================================================
C
         CASE(64) ! netcdf fort.64 file
C
C     ==================================================================
C     START NETCDF I/O OPTION - MCF
C
            IF(.NOT. ALLOCATED(uu22)) 
     &               ALLOCATE(uu22(num_records,timenc_len))
            IF(.NOT. ALLOCATED(vv22)) 
     &               ALLOCATE(vv22(num_records,timenc_len))
            IF(.NOT. ALLOCATED(xnc)) ALLOCATE(xnc(num_records))
            IF(.NOT. ALLOCATED(ync)) ALLOCATE(ync(num_records))
            IF(.NOT. ALLOCATED(nbvvnc)) ALLOCATE(nbvvnc(nbou,nvel))
            IF(.NOT. ALLOCATED(nvellnc)) ALLOCATE(nvellnc(nbou))
            IF(.NOT. ALLOCATED(nvdllnc)) ALLOCATE(nvdllnc(nope))
            IF(.NOT. ALLOCATED(ibtypenc)) ALLOCATE(ibtypenc(nbou))

            DO I=1, num_records
                IF (MNPROC.eq.1) THEN
                   uu22(i,timenc_len)= descript % array(i)
                   vv22(i,timenc_len)= descript % array2(i)
                ELSE IF (MNPROC.gt.1) THEN
                   uu22(i,timenc_len)= descript % array_g(i)
                   vv22(i,timenc_len)= descript % array2_g(i)
                ENDIF
            END DO

            if (ncmode.eq.1) then
            IF(ICS.EQ.2) THEN
               DO I=1, num_records
                  xnc(i)=slam(i)*RAD2DEG
                  ync(i)=sfea(i)*RAD2DEG
               END DO
            ELSE
               DO I=1,num_records
             	  xnc(i)=X(i)
                  ync(i)=Y(i)
               END DO
            ENDIF
            
           endif
 !               DO I=1,nbou
!                  DO J=1,nvel
!             	  nbvvnc(i,j)=nbvv(i,j)
!                  END DO
!               END DO
!
!            elevation specified boundary forcing segments
!
               netanc=0
               DO i=1,nope
                  nvdllnc(i)=nvdll(i)
                  netanc=netanc+nvdllnc(i)
!                 DO j=1,nvel
!             	     nbvvnc(i,j)=nbvv(i,j)
!                 END DO
              END DO

            
            IF(.NOT. ALLOCATED(nopenc)) ALLOCATE(nopenc(netanc))
               
               k=0
                  DO i=1,nope
                        DO j=1,nvdllnc(i)
                        k=k+1
             	        nopenc(k)=nbdv(i,j)
                     END DO
                  END DO
!
!              normal flow (discharge) specified boundary segments
!
               nvelnc=0
               DO i=1,nbou
                  nvellnc(i)=nvell(i)
                  nvelnc=nvelnc+nvellnc(i)
                  DO j=1,nvel
             	     nbvvnc(i,j)=nbvv(i,j)
                 END DO
               END DO

            IF(.NOT. ALLOCATED(nbounc)) ALLOCATE(nbounc(nvelnc))
               
               k=0
                  DO i=1,nbou
                        DO j=1,nvellnc(i)
                        k=k+1
             	        nbounc(k)=nbvv(i,j)
                     END DO
                  END DO
                  
            CALL writenc64 (descript % specifier,
     &             fileunit, timenc_len,
     &             num_records, timenc, 
     &             nt, ndsetsv, uu22, vv22,
     &             rundes, runid, agrid, irtype,
     &             dtdp,ihot, ics, nolibf, nolifa,
     &             nolica, nolicat, nwp,
     &             ncor, ntip, nws, nramp, tau0,
     &             statim,reftim, rnday, dramp,
     &             a00, b00, c00, h0,
     &             slam0, sfea0, cf,eslm, cori,
     &             ntif, nbfr, ne, nhy, nm,
     &             xnc, ync, dp,
     &             nope, neta, nbdv,nopenc,netanc,
     &             nvdllnc,
     &             nbou, nvel, nbvvnc,nbounc,nvelnc,
     &             nvellnc,ibtypenc,
     &             ncmode,
     &             title, institution, source, history,
     &             references, comments, host, convention, contact,
     &             base_date)

             NCIGVP=NCIGVP+num_records

C
C     END NETCDF I/O OPTION - MCF
C     ==================================================================
C
         CASE(71) ! netcdf fort.71 file
C
C     ==================================================================
C     START NETCDF I/O OPTION - MCF
C
            IF(.NOT. ALLOCATED(RMP00nc))
     &               ALLOCATE( RMP00nc (num_records, timenc_len) )
            IF(.NOT. ALLOCATED(xpnc)) ALLOCATE( xpnc(num_records) )
            IF(.NOT. ALLOCATED(ypnc)) ALLOCATE( ypnc(num_records) )
            IF(.NOT. ALLOCATED(nbvvnc)) ALLOCATE(nbvvnc(nbou,nvel))

            DO I=1, num_records
               IF (MNPROC.eq.1) THEN      ! serial
                  RMP00nc (i,timenc_len) = descript % array(I)
               ELSE IF (MNPROC.gt.1) THEN ! parallel
                  RMP00nc (i,timenc_len) = descript % array_g(I)
               ENDIF
            END DO
            if(ncmode.eq.1) then
            IF(ICS.EQ.2) THEN
               DO I=1,num_records
                  xpnc(i)=SLEM(I)*RAD2DEG
                  ypnc(i)=SFEM(I)*RAD2DEG
               END DO
            ELSE
                DO I=1,num_records
          	   xpnc(i)=XEM(i)
                   ypnc(i)=YEM(i)
               END DO
            ENDIF
!            endif
               DO I=1,nope
                  DO J=1,nvel
             	  nbvvnc(i,j)=nbvv(i,j)
                  END DO
               END DO
               
            endif 
            CALL writenc71 (descript % specifier, 
     &             fileunit, timenc_len, num_records, timenc,
     &             RMP00nc, nt,NTRSPM, rundes, runid, agrid, irtype,
     &             dtdp, ihot, ics, nolibf, nolifa,
     &             nolica, nolicat, nwp,
     &             ncor, ntip, nws, nramp, tau0,
     &             statim,reftim, rnday, dramp,
     &             a00, b00, c00, h0,
     &             slam0, sfea0, cf,eslm, cori,
     &             ntif, nbfr, ne, nhy, nm,
     &             xpnc, ypnc, dp,
     &             nope, neta, nbdv,
     &             nbou, nvel, nbvvnc,
     &             ncmode,
     &             title, institution, source, history,
     &             references, comments, host, convention, contact,
     &             base_date)

             NCIPSTP=NCIPSTP+num_records
C
C     END NETCDF I/O OPTION - MCF
C     ==================================================================
C
          CASE(72) ! netcdf fort.72 file
C
C     ==================================================================
C     START NETCDF I/O OPTION - MCF
C
!            WRITE(6,*)"fileunit = ", fileunit

            IF(.NOT.ALLOCATED(RMU00nc))
     &          ALLOCATE(RMU00nc(num_records,timenc_len))
            IF(.NOT.ALLOCATED(RMV00nc))
     &          ALLOCATE(RMV00nc(num_records,timenc_len))
            IF(.NOT.ALLOCATED(xwnc)) ALLOCATE(xwnc(num_records))
            IF(.NOT.ALLOCATED(ywnc)) ALLOCATE(ywnc(num_records))
            IF(.NOT. ALLOCATED(nbvvnc)) ALLOCATE(nbvvnc(nbou,nvel))

            DO I=1, num_records
               IF (MNPROC.eq.1) THEN ! serial
                  RMU00nc(i,timenc_len)= descript % array(i)
                  RMV00nc(i,timenc_len)= descript % array2(i)
               ELSE IF (MNPROC.gt.1) THEN ! parallel
                  RMU00nc(i,timenc_len)= descript % array_g(i)
                  RMV00nc(i,timenc_len)= descript % array2_g(i)
               ENDIF
            END DO
            
            if (ncmode.eq.1) then
            IF(ICS.EQ.2) THEN
               DO I=1, num_records
                  xwnc(i)=SLEM(I)*RAD2DEG
                  ywnc(i)=SFEM(I)*RAD2DEG
               END DO
            ELSE
               DO I=1,num_records
             	  xwnc(i)=XEM(i)
                  ywnc(i)=YEM(i)
               END DO
            ENDIF
!            endif
            
               DO I=1,nbou
                  DO J=1,nvel
             	  nbvvnc(i,j)=nbvv(i,j)
                  END DO
               END DO

            endif
            
            CALL writenc72 (descript % specifier, 
     &             fileunit, timenc_len, num_records, timenc,
     &             nt, NTRSPM, RMU00nc, RMV00nc,
     &             rundes, runid, agrid, irtype,
     &             dtdp,ihot, ics, nolibf, nolifa,
     &             nolica, nolicat, nwp,
     &             ncor, ntip, nws, nramp, tau0,
     &             statim,reftim, rnday, dramp,
     &             a00, b00, c00, h0,
     &             slam0, sfea0, cf,eslm, cori,
     &             ntif, nbfr, ne, nhy, nm,
     &             xwnc, ywnc, dp,
     &             nope, neta, nbdv,
     &             nbou, nvel, nbvvnc,
     &             ncmode,
     &             title, institution, source, history,
     &             references, comments, host, convention, contact,
     &             base_date)

             NCIWSTP=NCIWSTP+num_records


C
C     END NETCDF I/O OPTION - MCF
C     ==================================================================
C
          CASE(73) ! netcdf fort.73 file       
C
C     ==================================================================
C     START NETCDF I/O OPTION - MCF
C
            IF(.NOT. ALLOCATED(pr22))
     &               ALLOCATE (pr22(num_records,timenc_len))
            IF(.NOT. ALLOCATED(xnc)) ALLOCATE(xnc(num_records))
            IF(.NOT. ALLOCATED(ync)) ALLOCATE(ync(num_records))
            IF(.NOT. ALLOCATED(nbvvnc)) ALLOCATE(nbvvnc(nbou,nvel))
            IF(.NOT. ALLOCATED(nvellnc)) ALLOCATE(nvellnc(nbou))
            IF(.NOT. ALLOCATED(nvdllnc)) ALLOCATE(nvdllnc(nope))
            IF(.NOT. ALLOCATED(ibtypenc)) ALLOCATE(ibtypenc(nbou))

            DO I=1, num_records
               IF (MNPROC.eq.1) THEN ! serial
                   pr22 (i,timenc_len) = descript % array(i)
               ELSE IF (MNPROC.gt.1) THEN ! parallel
                   pr22 (i,timenc_len) = descript % array_g(i)
               ENDIF
            END DO
               
            if (ncmode.eq.1) then

            IF(ICS.EQ.2) THEN
               DO I=1, num_records
                        xnc(i)=slam(i)*RAD2DEG
                        ync(i)=sfea(i)*RAD2DEG
               END DO
            ELSE
               DO I=1,num_records
             	        xnc(i)=X(i)
                	ync(i)=Y(i)
               END DO
            ENDIF
            
           endif
!
!            elevation specified boundary forcing segments
!
               netanc=0
               DO i=1,nope
                  nvdllnc(i)=nvdll(i)
                  netanc=netanc+nvdllnc(i)
              END DO

            IF(.NOT. ALLOCATED(nopenc)) ALLOCATE(nopenc(netanc))
               
               k=0
                  DO i=1,nope
                        DO j=1,nvdllnc(i)
                        k=k+1
             	        nopenc(k)=nbdv(i,j)
                     END DO
                  END DO
!
!              normal flow (discharge) specified boundary segments
!
               nvelnc=0
               DO i=1,nbou
                  nvellnc(i)=nvell(i)
                  nvelnc=nvelnc+nvellnc(i)
                  DO j=1,nvel
             	     nbvvnc(i,j)=nbvv(i,j)
                 END DO
               END DO

            IF(.NOT. ALLOCATED(nbounc)) ALLOCATE(nbounc(nvelnc))
               
               k=0
                  DO i=1,nbou
                        DO j=1,nvellnc(i)
                        k=k+1
             	        nbounc(k)=nbvv(i,j)
                     END DO
                  END DO
             
               CALL writenc73 (descript % specifier, 
     &                fileunit, timenc_len, num_records, timenc,
     &                pr22, nt, NDSETSW, rundes, runid, agrid, irtype,
     &                dtdp, ihot, ics, nolibf, nolifa,
     &                nolica, nolicat, nwp,
     &                ncor, ntip, nws, nramp, tau0,
     &                statim,reftim, rnday, dramp,
     &                a00, b00, c00, h0,
     &                slam0, sfea0, cf,eslm, cori,
     &                ntif, nbfr, ne, nhy, nm,
     &                xnc, ync, dp,
     &                nope, neta, nbdv,nopenc,netanc,
     &                nvdllnc,
     &                nbou, nvel, nbvvnc,nbounc,nvelnc,
     &                nvellnc,ibtypenc,
     &                ncmode,
     &                title, institution, source, history,
     &                references, comments, host, convention, contact,
     &                base_date)

             NCIGPP=NCIGPP+num_records

C
C     END NETCDF I/O OPTION - MCF
C     ==================================================================
C
          CASE(74) ! netcdf fort.74 file 
C
C     ==================================================================
C     START NETCDF I/O OPTION - MCF
C
               IF(.NOT.ALLOCATED(WNCX22))
     &                 ALLOCATE(WNCX22(num_records,timenc_len))
               IF(.NOT.ALLOCATED(WNCY22))
     &                 ALLOCATE(WNCY22(num_records,timenc_len))
            IF(.NOT. ALLOCATED(nbvvnc)) ALLOCATE(nbvvnc(nbou,nvel))
            IF(.NOT. ALLOCATED(nvellnc)) ALLOCATE(nvellnc(nbou))
            IF(.NOT. ALLOCATED(nvdllnc)) ALLOCATE(nvdllnc(nope))
            IF(.NOT. ALLOCATED(ibtypenc)) ALLOCATE(ibtypenc(nbou))

               DO I=1, num_records
                   IF (MNPROC.eq.1) THEN ! serial
                      WNCX22(i,timenc_len)= descript % array(I)
                      WNCY22(i,timenc_len)= descript % array2(I)
                   ELSE IF (MNPROC.gt.1) THEN ! parallel
                      WNCX22(i,timenc_len)= descript % array_g(I)
                      WNCY22(i,timenc_len)= descript % array2_g(I)
                   ENDIF
               END DO
               
            if (ncmode.eq.1) then
               IF(ICS.EQ.2) THEN
                  DO I=1, num_records
                     xnc(i)=slam(i)*RAD2DEG
                     ync(i)=sfea(i)*RAD2DEG
                  END DO
               ELSE
                  DO I=1,num_records
             	     xnc(i)=X(i)
                     ync(i)=Y(i)
                  END DO
               ENDIF

            endif
!
!            elevation specified boundary forcing segments
!
               netanc=0
               DO i=1,nope
                  nvdllnc(i)=nvdll(i)
                  netanc=netanc+nvdllnc(i)
!                 DO j=1,nvel
!             	     nbvvnc(i,j)=nbvv(i,j)
!                 END DO
              END DO

            IF(.NOT. ALLOCATED(nopenc)) ALLOCATE(nopenc(netanc))
               
               k=0
                  DO i=1,nope
                        DO j=1,nvdllnc(i)
                        k=k+1
             	        nopenc(k)=nbdv(i,j)
                     END DO
                  END DO
!
!              normal flow (discharge) specified boundary segments
!
               nvelnc=0
               DO i=1,nbou
                  nvellnc(i)=nvell(i)
                  nvelnc=nvelnc+nvellnc(i)
                  DO j=1,nvel
             	     nbvvnc(i,j)=nbvv(i,j)
                 END DO
               END DO

            IF(.NOT. ALLOCATED(nbounc)) ALLOCATE(nbounc(nvelnc))
               
               k=0
                  DO i=1,nbou
                        DO j=1,nvellnc(i)
                        k=k+1
             	        nbounc(k)=nbvv(i,j)
                     END DO
                  END DO

               CALL writenc74 (descript % specifier, 
     &                fileunit, timenc_len, num_records, timenc,
     &                nt, NDSETSW, WNCX22, WNCY22,
     &                rundes, runid, agrid, irtype,
     &                dtdp,ihot, ics, nolibf, nolifa,
     &                nolica, nolicat, nwp,
     &                ncor, ntip, nws, nramp, tau0,
     &                statim,reftim, rnday, dramp,
     &                a00, b00, c00, h0,
     &                slam0, sfea0, cf,eslm, cori,
     &                ntif, nbfr, ne, nhy, nm,
     &                xnc, ync, dp,
     &                nope, neta, nbdv,nopenc,netanc,
     &                nvdllnc,
     &                nbou, nvel, nbvvnc,nbounc,nvelnc,
     &                nvellnc,ibtypenc,
     &                ncmode,
     &                title, institution, source, history,
     &                references, comments, host, convention, contact,
     &                base_date)
             NCIGWP=NCIGWP+num_records
C
C     END NETCDF I/O OPTION - MCF
C     ==================================================================
C
         CASE DEFAULT
            WRITE(ScreenUnit,*) 'ERROR: No netCDF for this filetype.'
         END SELECT
      ENDIF ! if MyProc.eq.0

#endif 
C-----------------------------------------------------------------------
       END SUBROUTINE writeOutArrayNetCDF
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   W R I T E  H O T S T A R T
C-----------------------------------------------------------------------
C     jgf48.03 This subroutine was created from code in timestep.F to
C     write hotstart files.
C-----------------------------------------------------------------------
      SUBROUTINE writeHotstart(Time, IT)

      USE SIZES
      USE GLOBAL
      USE GLOBAL_IO, ONLY: packOne, unpackOne, packOneInt, unpackOneInt,
     &    packTwo, unpackTwo, HEADER_MAX,
     &    collectFullDomainArray, collectFullDomainIntArray
#ifdef CMPI
      USE MESSENGER
#endif
#ifdef NETCDF
      USE NETCDFIO, ONLY :
     &     timenc, timenc_len, fileunit, irtype,
     &     xenc, yenc, xvnc, yvnc, xpnc, ypnc, xnc, ync, xwnc, ywnc,
     &     et00nc, UU00nc, VV00nc, RMP00nc, RMU00nc, RMV00nc, nbvvnc,
     &     eta11,eta22,uu22, vv22, pr22, WNCX22, WNCY22,
     &     nbounc,nvelnc,nvellnc,
     &     nopenc,netanc,nvdllnc, nodecodenc,noffnc, etadd,
     &     writenc67
       USE NodalAttributes, ONLY :
     &     nolibf, nwp, tau0, cf, eslm

#endif

      IMPLICIT NONE
      REAL(8), intent(in) :: Time
      INTEGER, intent(in) :: IT
C
C     START ADDITIONS FOR NETCDF by MCF
C
      INTEGER ncmode
      INTEGER num_records, num_recordse ! to save in the netcdf file
      INTEGER I,J,K
C
C     FINISH ADDITIONS FOR NETCDF by MCF
C

      type(OutputDataDescript_t) :: descript 
      LOGICAL TEST_HOTSTART     !jgf45.07 used for testing hot start capability
      INTEGER npx, nex

      type(OutputDataDescript_t), SAVE :: Elev1Descript
      type(OutputDataDescript_t), SAVE :: Elev2Descript 
      type(OutputDataDescript_t), SAVE :: VelDescript
      type(OutputDataDescript_t), SAVE :: CH1Descript  
      type(OutputDataDescript_t), SAVE :: EtaDiscDescript
      type(OutputDataDescript_t), SAVE :: NodeCodeDescript
      type(OutputDataDescript_t), SAVE :: NOFFDescript
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: ETA1_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: ETA2_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: EtaDisc_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: UU2_g(:)     
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: VV2_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: CH1_g(:)
      INTEGER, SAVE, ALLOCATABLE, TARGET :: NodeCode_g(:)
      INTEGER, SAVE, ALLOCATABLE, TARGET :: NOFF_g(:)
      LOGICAL, SAVE :: FirstCall = .true.
      
      WRITE(16,*) "***********WRITE HOTSTART*********************"
      WRITE(16,*) "NHSTAR = ", NHSTAR
      WRITE(16,*) "FirstCall = ", FirstCall
      IF (FirstCall) THEN
         write(16,*) 'Allocating Elev1Descript' 
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(ETA1_g(NP_G))
         ENDIF
         Elev1Descript % specifier            = NHSTAR
         Elev1Descript % initial_value        = 0.0
         Elev1Descript % num_items_per_record = 1
         Elev1Descript % num_fd_records       = NP_G
         Elev1Descript % num_records_this     = NP
         Elev1Descript % imap                 => NODES_LG
         Elev1Descript % array                => ETA1
         Elev1Descript % array_g              => ETA1_g

         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            IF (.not.ALLOCATED(ETA2_g)) THEN
               write(16,*) 'Allocating Elev2Descript' 
               ALLOCATE(ETA2_g(NP_G))
            ENDIF
         ENDIF
         Elev2Descript % specifier            = NHSTAR
         Elev2Descript % initial_value        = 0.0
         Elev2Descript % num_items_per_record = 1
         Elev2Descript % num_fd_records       = NP_G
         Elev2Descript % num_records_this     = NP
         Elev2Descript % imap                 => NODES_LG
         Elev2Descript % array                => ETA2
         Elev2Descript % array_g              => ETA2_g
         
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            IF (.not.ALLOCATED(UU2_g)) THEN
               write(16,*) 'Allocating VelDescript' 
               ALLOCATE(UU2_g(NP_G))
               ALLOCATE(VV2_g(NP_G)) 
            ENDIF
         ENDIF
         VelDescript % specifier            = NHSTAR
         VelDescript % initial_value        = 0.0
         VelDescript % num_items_per_record = 2
         VelDescript % num_fd_records       = NP_G
         VelDescript % num_records_this     = NP
         VelDescript % imap                 => NODES_LG
         VelDescript % array                => UU2
         VelDescript % array2               => VV2
         VelDescript % array_g              => UU2_g
         VelDescript % array2_g             => VV2_g

         IF (IM.eq.10) THEN
            write(16,*) 'Allocating CH1Descript' 
            IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
               ALLOCATE(EtaDisc_g(NP_G))
            ENDIF
            CH1Descript % specifier            = NHSTAR
            CH1Descript % initial_value        = 0.0
            CH1Descript % num_items_per_record = 1
            CH1Descript % num_fd_records       = NP_G
            CH1Descript % num_records_this     = NP
            CH1Descript % imap                 => NODES_LG
            CH1Descript % array                => CH1
            CH1Descript % array_g              => CH1_g
         ENDIF        
 
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            write(16,*) 'Allocating EtaDiscDescript' 
            ALLOCATE(EtaDisc_g(NP_G))
         ENDIF
         EtaDiscDescript % specifier            = NHSTAR
         EtaDiscDescript % initial_value        = 0.0
         EtaDiscDescript % num_items_per_record = 1
         EtaDiscDescript % num_fd_records       = NP_G
         EtaDiscDescript % num_records_this     = NP
         EtaDiscDescript % imap                 => NODES_LG
         EtaDiscDescript % array                => EtaDisc
         EtaDiscDescript % array_g              => EtaDisc_g
         
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            write(16,*) 'Allocating NodeCodeDescript' 
            ALLOCATE(NodeCode_g(NP_G))
         ENDIF
         NodeCodeDescript % specifier            = NHSTAR
         NodeCodeDescript % int_initial_value        = 0
         NodeCodeDescript % num_items_per_record = 1
         NodeCodeDescript % num_fd_records       = NP_G
         NodeCodeDescript % num_records_this     = NP
         NodeCodeDescript % imap                 => NODES_LG
         NodeCodeDescript % iarray               => NODECODE
         NodeCodeDescript % iarray_g             => NODECODE_g

         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            write(16,*) 'Allocating NOFFDescript' 
            ALLOCATE(NOFF_g(NE_G))
         ENDIF
         NOFFDescript % specifier            = NHSTAR
         NOFFDescript % int_initial_value        = 0
         NOFFDescript % num_items_per_record = 1
         NOFFDescript % num_fd_records       = NE_G
         NOFFDescript % num_records_this     = NE
         NOFFDescript % imap                 => IMAP_EL_LG
         NOFFDescript % iarray               => NOFF
         NOFFDescript % iarray_g             => NOFF_g

         FirstCall = .false.
      ENDIF

C     collect up the data from subdomains if running in parallel
#ifdef CMPI
         WRITE(16,*) "START collectFullDomainArray"
         WRITE(16,*) 
     &    'About to collectFullDomainArrays for hotstart file.'
         WRITE(16,*) "collectFullDomainArray Elev1Descript"
         CALL collectFullDomainArray(Elev1Descript, packOne, unpackOne)
         WRITE(16,*) "collectFullDomainArray Elev2Descript"
         CALL collectFullDomainArray(Elev2Descript, packOne, unpackOne)
         WRITE(16,*) "collectFullDomainArray VelDescript"
         CALL collectFullDomainArray(VelDescript, packTwo, unpackTwo)
         IF (IM.eq.10) THEN
            CALL collectFullDomainArray(CH1Descript, packOne, unpackOne)
         ENDIF
         WRITE(16,*) "collectFullDomainArray EtaDiscDescript"
         CALL collectFullDomainArray(EtaDiscDescript, 
     &             packOne, unpackOne)
         WRITE(16,*) "collectFullDomainIntArray NodeCodeDescript"
         CALL collectFullDomainIntArray(NodeCodeDescript, 
     &             packOneInt, unpackOneInt)
         WRITE(16,*) "collectFullDomainIntArray NOFFDescript"
         CALL collectFullDomainIntArray(NOFFDescript, 
     &             packOneInt, unpackOneInt)
         WRITE(16,*) "FINISH collectFullDomainArray"
         IF ( (MyProc.eq.0) ) THEN
           DO I=1, NP_G
!             PRINT *,"before case 1 ",I, Elev1Descript % array_g(I)
!             PRINT *,"before case 1 ",I, Elev1Descript % array_g(I)
          ENDDO
         ENDIF

#endif
      SELECT CASE (NHSTAR)

      CASE(1) ! nonportable binary (1 for backwards compatibility)

         NPX = NP_G
         NEX = NE_G
!         PRINT *, "write_output: NPX - 1", NPX
!         PRINT *, "write_output: NEX - 1", NEX

         if (MNPROC == 1 .or. WRITE_LOCAL_HOT_START_FILES) THEN
            NPX = NP
            NEX = NE
!            PRINT *, "write_output: NPX - 2", NPX
!            PRINT *, "write_output: NEX - 2", NEX
         endif

      IF (MYPROC == 0 .or. WRITE_LOCAL_HOT_START_FILES) THEN 
          OPEN(hss % lun ,FILE=TRIM(HOTSTARTDIR)//'/'//hss % filename,
     $      ACCESS='DIRECT',RECL=8)
          IHOTSTP=1
!          print *, FileFmtVersion
          WRITE(hss % lun,REC=IHOTSTP) FileFmtVersion 
          IHOTSTP = IHOTSTP + 1
!          print *, IM
          WRITE(hss % lun,REC=IHOTSTP) IM   ; IHOTSTP = IHOTSTP + 1
!          print *, TIME
          WRITE(hss % lun,REC=IHOTSTP) TIME ; IHOTSTP = IHOTSTP + 1
!          print *, IT
          WRITE(hss % lun,REC=IHOTSTP) IT   ; IHOTSTP = IHOTSTP + 1
!          print *, NPX
          WRITE(hss % lun,REC=IHOTSTP) NPX  ; IHOTSTP = IHOTSTP + 1
!          print *, NEX
          WRITE(hss % lun,REC=IHOTSTP) NEX  ; IHOTSTP = IHOTSTP + 1
!          print *, NPX
          WRITE(hss % lun,REC=IHOTSTP) NPX  ; IHOTSTP = IHOTSTP + 1
!          print *, NEX
          WRITE(hss % lun,REC=IHOTSTP) NEX  ; IHOTSTP = IHOTSTP + 1
          
!          PRINT *,"========================before bin=================="
!          PRINT *, "==== write_output: IM = ", IM
!          PRINT *, "==== write_output: TIME = ", TIME
!          PRINT *, "==== write_output: IM = ", IT
!          PRINT *, "==== write_output: NPX = ", NPX
!          PRINT *, "==== write_output: NEX = ", NEX
!          PRINT *, "==== write_output: IESTP = ", IESTP
!          PRINT *, "==== write_output: NSCOUE = ", NSCOUE
!          PRINT *, "==== write_output: IVSTP = ", IVSTP
!          PRINT *, "==== write_output: NSCOUV = ", NSCOUV
!          PRINT *, "==== write_output: ICSTP = ", ICSTP
!          PRINT *, "==== write_output: NSCOUC = ", NSCOUC
!          PRINT *, "==== write_output: IPSTP = ", IPSTP
!
!          PRINT *, "==== write_output: IWSTP = ", IWSTP
!          PRINT *, "==== write_output: NSCOUM = ", NSCOUM
!          PRINT *, "==== write_output: IGEP = ", IGEP
!          PRINT *, "==== write_output: NSCOUGE = ", NSCOUGE
!          PRINT *, "==== write_output: IGVP = ", IGVP
!
!          PRINT *, "==== write_output: NSCOUGV = ", NSCOUGV
!          PRINT *, "==== write_output: IGCP = ", IGCP
!          PRINT *, "==== write_output: NSCOUGC = ", NSCOUGC
!          PRINT *, "==== write_output: IGPP = ", IGPP
!          PRINT *, "==== write_output: IGWP = ", IGWP
!          PRINT *, "==== write_output: NSCOUGW = ", NSCOUGW
!          PRINT *,"========================before bin=================="
      
#ifdef CMPI
          DO I=1, NP_G
             WRITE(hss % lun,REC=IHOTSTP) Elev1Descript % array_g(I)
             IHOTSTP=IHOTSTP+1
          ENDDO
          DO I=1, NP_G
             WRITE(hss % lun,REC=IHOTSTP) Elev2Descript % array_g(I)
             IHOTSTP=IHOTSTP+1
          ENDDO
          DO I=1, NP_G
             WRITE(hss % lun,REC=IHOTSTP) EtaDiscDescript % array_g(I)
             IHOTSTP=IHOTSTP+1
          ENDDO
          DO I=1, NP_G
             WRITE(hss % lun,REC=IHOTSTP) VelDescript % array_g(I)
             IHOTSTP=IHOTSTP+1
          ENDDO
          DO I=1, NP_G
             WRITE(hss % lun,REC=IHOTSTP) VelDescript % array2_g(I)
             IHOTSTP=IHOTSTP+1
          ENDDO
          IF (IM.eq.10) THEN
             DO I=1, NP_G
                WRITE(hss % lun,REC=IHOTSTP) CH1Descript % array_g(I)
                IHOTSTP=IHOTSTP+1
             ENDDO
          ENDIF
          DO I=1, NP_G
             WRITE(hss % lun,REC=IHOTSTP) NodeCodeDescript % iarray_g(I)
             IHOTSTP=IHOTSTP+1
          ENDDO
          DO I=1, NE_G
             WRITE(hss % lun,REC=IHOTSTP) NOFFDescript % iarray_g(I)
             IHOTSTP=IHOTSTP+1
          ENDDO
#else
          WRITE(16,*) '=================ETA1:===================='
          DO I=1, NP_G
             WRITE(16,*) I, ETA1(I)
             WRITE(hss % lun,REC=IHOTSTP) ETA1(I)
             IHOTSTP=IHOTSTP+1
          ENDDO
          WRITE(16,*) '=================ETA2:===================='
          DO I=1, NP_G
             WRITE(16,*) I, ETA2(I)
             WRITE(hss % lun,REC=IHOTSTP) ETA2(I)
             IHOTSTP=IHOTSTP+1
          ENDDO
          WRITE(16,*) '=================ETADisc:===================='
          DO I=1, NP_G
             WRITE(16,*) I, ETADisc(I)
             WRITE(hss % lun,REC=IHOTSTP) EtaDisc(I)
             IHOTSTP=IHOTSTP+1
          ENDDO
          WRITE(16,*) '=================UU2:===================='
          DO I=1, NP_G
             WRITE(16,*) I, UU2(I)
             WRITE(hss % lun,REC=IHOTSTP) UU2(I)
             IHOTSTP=IHOTSTP+1
          ENDDO
          WRITE(16,*) '=================VV2:===================='
          DO I=1, NP_G
             WRITE(16,*) I, VV2(I) 
             WRITE(hss % lun,REC=IHOTSTP) VV2(I)
             IHOTSTP=IHOTSTP+1
          ENDDO
          IF (IM.eq.10) THEN
             DO I=1, NP_G
                WRITE(hss % lun,REC=IHOTSTP) CH1(I)
                IHOTSTP=IHOTSTP+1
             ENDDO
          ENDIF
          WRITE(16,*) '=================NodeCode:===================='
          DO I=1, NP_G
             WRITE(16,*) I, NodeCode(I)
             WRITE(hss % lun,REC=IHOTSTP) NodeCode(I)
             IHOTSTP=IHOTSTP+1
          ENDDO
          WRITE(16,*) '=================NOFF:===================='
          DO I=1, NE_G
             WRITE(16,*) I, NOFF(I)
             WRITE(hss % lun,REC=IHOTSTP) NOFF(I)
             IHOTSTP=IHOTSTP+1
          ENDDO
#endif

         WRITE(hss % lun,REC=IHOTSTP) IESTP   ; IHOTSTP = IHOTSTP + 1
         WRITE(hss % lun,REC=IHOTSTP) NSCOUE  ; IHOTSTP = IHOTSTP + 1

         WRITE(hss % lun,REC=IHOTSTP) IVSTP   ; IHOTSTP = IHOTSTP + 1
         WRITE(hss % lun,REC=IHOTSTP) NSCOUV  ; IHOTSTP = IHOTSTP + 1

         WRITE(hss % lun,REC=IHOTSTP) ICSTP   ; IHOTSTP = IHOTSTP + 1
         WRITE(hss % lun,REC=IHOTSTP) NSCOUC  ; IHOTSTP = IHOTSTP + 1

         WRITE(hss % lun,REC=IHOTSTP) IPSTP   ; IHOTSTP = IHOTSTP + 1
         WRITE(hss % lun,REC=IHOTSTP) IWSTP   ; IHOTSTP = IHOTSTP + 1
         WRITE(hss % lun,REC=IHOTSTP) NSCOUM  ; IHOTSTP = IHOTSTP + 1

         WRITE(hss % lun,REC=IHOTSTP) IGEP    ; IHOTSTP = IHOTSTP + 1
         WRITE(hss % lun,REC=IHOTSTP) NSCOUGE ; IHOTSTP = IHOTSTP + 1

         WRITE(hss % lun,REC=IHOTSTP) IGVP    ; IHOTSTP = IHOTSTP + 1
         WRITE(hss % lun,REC=IHOTSTP) NSCOUGV ; IHOTSTP = IHOTSTP + 1

         WRITE(hss % lun,REC=IHOTSTP) IGCP    ; IHOTSTP = IHOTSTP + 1
         WRITE(hss % lun,REC=IHOTSTP) NSCOUGC ; IHOTSTP = IHOTSTP + 1

         WRITE(hss % lun,REC=IHOTSTP) IGPP    ; IHOTSTP = IHOTSTP + 1
         WRITE(hss % lun,REC=IHOTSTP) IGWP    ; IHOTSTP = IHOTSTP + 1
         WRITE(hss % lun,REC=IHOTSTP) NSCOUGW ; IHOTSTP = IHOTSTP + 1
       
!          PRINT *,"========================before bin=================="
!          PRINT *, "==== write_output: IM = ", IM
!          PRINT *, "==== write_output: TIME = ", TIME
!          PRINT *, "==== write_output: IM = ", IT
!          PRINT *, "==== write_output: NPX = ", NPX
!          PRINT *, "==== write_output: NEX = ", NEX
!          PRINT *, "==== write_output: IESTP = ", IESTP
!          PRINT *, "==== write_output: NSCOUE = ", NSCOUE
!          PRINT *, "==== write_output: IVSTP = ", IVSTP
!          PRINT *, "==== write_output: NSCOUV = ", NSCOUV
!          PRINT *, "==== write_output: ICSTP = ", ICSTP
!          PRINT *, "==== write_output: NSCOUC = ", NSCOUC
!          PRINT *, "==== write_output: IPSTP = ", IPSTP
!
!          PRINT *, "==== write_output: IWSTP = ", IWSTP
!          PRINT *, "==== write_output: NSCOUM = ", NSCOUM
!          PRINT *, "==== write_output: IGEP = ", IGEP
!          PRINT *, "==== write_output: NSCOUGE = ", NSCOUGE
!          PRINT *, "==== write_output: IGVP = ", IGVP

!          PRINT *, "==== write_output: NSCOUGV = ", NSCOUGV
!          PRINT *, "==== write_output: IGCP = ", IGCP
!          PRINT *, "==== write_output: NSCOUGC = ", NSCOUGC
!          PRINT *, "==== write_output: IGPP = ", IGPP
!          PRINT *, "==== write_output: IGWP = ", IGWP
!          PRINT *, "==== write_output: NSCOUGW = ", NSCOUGW
!          PRINT *,"========================after bin==================="
          
       ENDIF

!This code is not supported yet.
#if 0
C...
C...  ADD IN 3D STUFF IF APPROPRIATE
C...
      IF (C3D) THEN
         CALL HSTART3D_OUT()
      ENDIF
C...  
C...  IF APPROPRIATE ADD HARMONIC ANALYSIS INFORMATION TO HOT START FILE
C...  
      IF((IHARIND.EQ.1).AND.(IT.GT.ITHAS)) THEN
         WRITE(hss % lun,REC=IHOTSTP+1) ICHA
         IHOTSTP = IHOTSTP + 1
         CALL HAHOUT(NP,NSTAE,NSTAV,NHASE,NHASV,NHAGE,NHAGV,
     &      hss % lun,IHOTSTP)
C         
         IF(NHASE.EQ.1) CALL HAHOUTES(NSTAE,hss % lun,IHOTSTP)
         IF(NHASV.EQ.1) CALL HAHOUTVS(NSTAV,hss % lun,IHOTSTP)
         IF(NHAGE.EQ.1) CALL HAHOUTEG(NP,hss % lun,IHOTSTP)
         IF(NHAGV.EQ.1) CALL HAHOUTVG(NP,hss % lun,IHOTSTP)
      ENDIF
        
      if( CHARMV) then
         IF((IHARIND.EQ.1).AND.(IT.GT.ITMV)) THEN
            WRITE(hss % lun,REC=IHOTSTP) NTSTEPS
            IHOTSTP=IHOTSTP+1
            IF(NHAGE.EQ.1) THEN
               DO I=1,NP
                  WRITE(hss % lun,REC=IHOTSTP) ELAV(I)
                  IHOTSTP = IHOTSTP + 1
                  WRITE(hss % lun,REC=IHOTSTP) ELVA(I)
                  IHOTSTP = IHOTSTP + 1
               END DO
            ENDIF
            IF(NHAGV.EQ.1) THEN
              DO I=1,NP
                WRITE(hss % lun,REC=IHOTSTP) XVELAV(I)
                IHOTSTP=IHOTSTP+1
                WRITE(hss % lun,REC=IHOTSTP) YVELAV(I)
                IHOTSTP=IHOTSTP+1
                WRITE(hss % lun,REC=IHOTSTP) XVELVA(I)
                IHOTSTP=IHOTSTP+1
                WRITE(hss % lun,REC=IHOTSTP) YVELVA(I)
                IHOTSTP=IHOTSTP+1
              END DO
            ENDIF
          ENDIF
      endif                   !  charmv
        
#endif
        
      CASE(2) ! ascii?
         write(16,*) 'ASCII Hotstart file not implemented yet.'

      CASE(3) ! netcdf 
         write(16,*) 'NetCDF Hotstart undergoing testing.'

#ifdef NETCDF
C
C     START ADDITIONS FOR NETCDF by MCF 4/8/08
C

      IF (MNPROC.gt.1) ncmode=2
      IF (MNPROC.eq.1) ncmode=1
      
      IF (ncmode.eq.1.or. (ncmode.eq.2. and.MYPROC.eq. 0)) THEN
      
         NPX = NP_G
         NEX = NE_G

#ifdef NETCDF_DEBUG
        IF ( (MyProc.eq.0) )PRINT *, "ncmode = ", ncmode
        IF ( (MyProc.eq.0) )PRINT *, "NP_G = ", NP_G,"NE_G = ", NE_G
#endif 
      
         timenc_len=1
         IF(.NOT. ALLOCATED(timenc)) ALLOCATE(timenc(timenc_len))
         timenc(timenc_len)=time

#ifdef NETCDF_DEBUG
        IF ( (MyProc.eq.0) )PRINT *, "time = ", time
#endif 
C
C     FINISH ADDITIONS FOR NETCDF by MCF
C
!#endif
!      IF (MNPROC.gt.1) num_records = descript % num_fd_records   ! parallel
!      IF (MNPROC.eq.1) num_records = descript % num_records_this ! serial
       IF (MNPROC.gt.1) num_records = NP_G   ! parallel
       IF (MNPROC.eq.1) num_records = NP ! serial

#ifdef NETCDF_DEBUG
      IF ( (MyProc.eq.0) )PRINT *, "num_records = ", num_records
#endif 

      IF(.NOT.ALLOCATED(eta11))ALLOCATE (eta11(num_records,timenc_len))

!      PRINT *, " Start Do Loop "
            DO I=1, num_records
             IF (MNPROC.gt.1.and.MyProc.eq.0) THEN                                     ! parallel
!                   eta11(i,timenc_len) = descript % array_g(i)
                   eta11(i,timenc_len) = Elev1Descript % array_g(i)
!                   PRINT *,i,"  Parallel ETA11 = ", ETA11(i, timenc_len)
               ELSE IF (MNPROC.eq.1) THEN                               ! serial
                   eta11(i,timenc_len) = Elev1Descript % array(i)
!                   eta11(i,timenc_len) = descript % array(i)
!                   PRINT *,i,"  Serial ETA11 = ", ETA11(i, timenc_len)
               ENDIF
            END DO
#ifdef NETCDF_DEBUG
            IF ( MNPROC.eq.1.or.(MNPROC.gt.1.and.MyProc.eq.0) ) 
     &           PRINT *,"  ETA11 = ", ETA11(11, timenc_len)
#endif 
#ifdef NETCDF_DEBUG
            IF ( MyProc.eq.0 ) 
     &           PRINT *,"  ETA11 = ", ETA11(11, timenc_len)
#endif 
C   End Addition for netCDF by MCF 4/8/08
#endif

#ifdef NETCDF
C   Start Addition for netCDF by MCF 4/8/08
      IF(.NOT.ALLOCATED(eta22))ALLOCATE (eta22(num_records,timenc_len))

            DO I=1, num_records
              IF (MNPROC.gt.1.and.MyProc.eq.0) THEN                                     ! parallel
                   eta22(i,timenc_len) = Elev2Descript % array_g(i)
              ELSE IF (MNPROC.eq.1) THEN                                ! serial
                   eta22(i,timenc_len) = Elev2Descript % array(i)
!                   PRINT *,i,"  ETA22 = ", ETA22(i, timenc_len)
              ENDIF
            END DO
#ifdef NETCDF_DEBUG
        IF ( (MyProc.eq.0) ) PRINT *,"  ETA22 = ", ETA22(11, timenc_len)
#endif 
C   End Addition for netCDF by MCF 4/8/08
#endif

#ifdef NETCDF
C   Start Addition for netCDF by MCF 4/8/08
      IF(.NOT.ALLOCATED(etadd))ALLOCATE (etadd(num_records,timenc_len))

            DO I=1, num_records
              IF (MNPROC.gt.1.and.MyProc.eq.0) THEN                                     ! parallel
                   etadd(i,timenc_len) = EtaDiscDescript % array_g(i)
              ELSE IF (MNPROC.eq.1) THEN                                ! serial
                   etadd(i,timenc_len) = EtaDiscDescript % array(i)
!                   PRINT *,i,"  ETADD = ", ETADD(i, timenc_len)
              ENDIF
            END DO
#ifdef NETCDF_DEBUG
        IF ( (MyProc.eq.0) ) PRINT *,"  ETADD = ", ETADD(11, timenc_len)
#endif 

C   End Addition for netCDF by MCF 4/8/08
#endif

#ifdef NETCDF
C   Start Addition for netCDF by MCF 4/8/08
      IF(.NOT.ALLOCATED(uu22))ALLOCATE (uu22(num_records,timenc_len))

            DO I=1, num_records
              IF (MNPROC.gt.1.and.MyProc.eq.0) THEN                                     ! parallel
                   uu22(i,timenc_len) = VelDescript % array_g(i)
              ELSE IF (MNPROC.eq.1) THEN                                ! serial
                   uu22(i,timenc_len) = VelDescript % array(i)
!                   PRINT *,i,"  uu22 = ", uu22(i, timenc_len)
              ENDIF
            END DO
#ifdef NETCDF_DEBUG
        IF ( (MyProc.eq.0) ) PRINT *,"  uu22 = ", uu22(11, timenc_len)
#endif         

C   End Addition for netCDF by MCF 4/8/08
#endif         
#ifdef NETCDF
C   Start Addition for netCDF by MCF 4/8/08
      IF(.NOT.ALLOCATED(vv22))ALLOCATE (vv22(num_records,timenc_len))

            DO I=1, num_records
               IF (MNPROC.gt.1.and.MyProc.eq.0) THEN                                     ! parallel
                  vv22(i,timenc_len) = VelDescript % array2_g(i)
              ELSE IF (MNPROC.eq.1) THEN                                ! serial
                  vv22(i,timenc_len) = VelDescript % array2(i)
!                   PRINT *,i,"  vv22 = ", vv22(i, timenc_len)
              ENDIF
            END DO
#ifdef NETCDF_DEBUG
        IF ( (MyProc.eq.0) ) PRINT *,"  vv22 = ", vv22(11, timenc_len)
#endif         
C   End Addition for netCDF by MCF 4/8/08
#endif 
#ifdef NETCDF     
C   Start Addition for netCDF by MCF 4/8/08
      num_records=np_g
      IF(.NOT.ALLOCATED(nodecodenc))ALLOCATE (nodecodenc(num_records))
!      IF (MNPROC.gt.1.and.MyProc.eq.0 )PRINT *,'num_records',num_records
            DO I=1, num_records
               IF (MNPROC.gt.1.and.MyProc.eq.0 ) THEN                                     ! parallel
                   nodecodenc(i) = NodeCodeDescript % iarray_g(i)     
!                      PRINT *,i,"  nodecodenc = ", nodecodenc(i)
              ELSE IF (MNPROC.eq.1) THEN                                ! serial
                   nodecodenc(i) = NodeCodeDescript % iarray(i)
!                   PRINT *,i,"  nodecodenc = ", nodecodenc(i)
              ENDIF
            END DO
#ifdef NETCDF_DEBUG
        IF ( (MyProc.eq.0) ) PRINT *,"  nodecodenc = ", nodecodenc(11)
#endif         
C   End Addition for netCDF by MCF 4/8/08
#endif
        
#ifdef NETCDF
C   Start Addition for netCDF by MCF 4/8/08
!      num_recordse=ne_g
!      IF(.NOT.ALLOCATED(noffnc))ALLOCATE (noffnc(num_recordse))

       IF (MNPROC.gt.1) num_recordse = NE_G   ! parallel
       IF (MNPROC.eq.1) num_recordse = NE ! serial

      IF(.NOT.ALLOCATED(noffnc))ALLOCATE (noffnc(num_recordse))
!      IF (MNPROC.gt.1.and.MyProc.eq.0 )
!     . PRINT *,'wo: num_recordse',num_recordse

            DO I=1, num_recordse
              IF (MNPROC.gt.1.and.MyProc.eq.0) THEN                     ! parallel
                    noffnc(i) = NOFFDescript % iarray_g(i)
                    noffnc(i)=min(noffnc(i),1)                          ! noff ranges from 0 to 1
!                    PRINT *,i,"  noffnc = ", noffnc(i)
              ELSE IF (MNPROC.eq.1) THEN                                ! serial
                    noffnc(i) = NOFFDescript % iarray(i)
              ENDIF
            END DO
#ifdef NETCDF_DEBUG
        IF ( (MyProc.eq.0) ) PRINT *,"  noffnc = ", noffnc(11)
#endif         
C   End Addition for netCDF by MCF 4/8/08
#endif
C
C    ADDITION OF netCDF Hostart file 67.nc and 68.nc by MCF 4/8/08
C
#ifdef NETCDF
!        IF ( (MyProc.eq.0) ) PRINT *,"  write output ncmode = ", ncmode
            IF(.NOT. ALLOCATED(xnc)) ALLOCATE(xnc(num_records))
            IF(.NOT. ALLOCATED(ync)) ALLOCATE(ync(num_records))
            IF(.NOT. ALLOCATED(nbvvnc)) ALLOCATE(nbvvnc(nbou,nvel))
            IF(.NOT. ALLOCATED(nvellnc)) ALLOCATE(nvellnc(nbou))
            IF(.NOT. ALLOCATED(nvdllnc)) ALLOCATE(nvdllnc(nope))
            IF(.NOT. ALLOCATED(ibtypenc)) ALLOCATE(ibtypenc(nbou))
         if (ncmode.eq.2) then
            if(myproc.eq.0) PRINT *,"wo: allocate nopenc"
!            IF(.NOT. ALLOCATED(nopenc)) ALLOCATE(nopenc(netanc))
!            IF(.NOT. ALLOCATED(nbounc)) ALLOCATE(nbounc(nvelnc))
         endif
        IF (ncmode.eq.1) THEN                                           ! serial
        IF ( (MyProc.eq.0) ) PRINT *,"  arrays allocated = "


            IF(ICS.EQ.2) THEN
               DO I=1, num_records
                  xnc(i)=slam(i)*RAD2DEG
                  ync(i)=sfea(i)*RAD2DEG
               END DO
            ELSE
               DO I=1,num_records
                  xnc(i)=X(i)
                  ync(i)=Y(i)
               END DO
            ENDIF

          ENDIF                                                         ! end serial
!
!            elevation specified boundary forcing segments
!
        IF ((MyProc.eq.0))PRINT *," write output: elevation specified "
               netanc=0
               DO i=1,nope
                  nvdllnc(i)=nvdll(i)
                  netanc=netanc+nvdllnc(i)
!                 DO j=1,nvel
!                    nbvvnc(i,j)=nbvv(i,j)
!                 END DO
              END DO

             IF(.NOT. ALLOCATED(nopenc)) ALLOCATE(nopenc(netanc))
#ifdef NETCDF_DEBUG
      IF ( (MyProc.eq.0) ) PRINT *," write output: nopenc "
#endif         

               k=0
                  DO i=1,nope
                        DO j=1,nvdllnc(i)
                        k=k+1
                        nopenc(k)=nbdv(i,j)
                     END DO
                  END DO
!
!              normal flow (discharge) specified boundary segments
!
#ifdef NETCDF_DEBUG
        IF ((MyProc.eq.0))PRINT *," write output: normal flow  "
#endif         

               nvelnc=0
               DO i=1,nbou
                  nvellnc(i)=nvell(i)
                  nvelnc=nvelnc+nvellnc(i)
                  DO j=1,nvel
                     nbvvnc(i,j)=nbvv(i,j)
                 END DO
               END DO

#ifdef NETCDF_DEBUG
        IF ((MyProc.eq.0))PRINT *," write output: nbounc  "
#endif         
       
            IF(.NOT. ALLOCATED(nbounc)) ALLOCATE(nbounc(nvelnc))

               k=0
                  DO i=1,nbou
                        DO j=1,nvellnc(i)
                        k=k+1
                        nbounc(k)=nbvv(i,j)
                     END DO
                  END DO

!       IF (ncmode.eq.1. or. MYPROC == 0) THEN
!       
!          PRINT *,"========================net===================="
!          PRINT *,"FileFmtMajor, FileFmtMinor, FileFmtRev",
!     &             FileFmtMajor, FileFmtMinor, FileFmtRev
!          PRINT *,"FileFmtVersion = ",FileFmtVersion
!          PRINT *, "==== write_output: IM = ", IM
!          PRINT *, "==== write_output: TIME = ", TIME
!          PRINT *, "==== write_output: IM = ", IT
!          PRINT *, "==== write_output: NPX = ", NPX
!          PRINT *, "==== write_output: NEX = ", NEX
!          PRINT *, "==== write_output: IESTP = ", IESTP
!          PRINT *, "==== write_output: NSCOUE = ", NSCOUE
!          PRINT *, "==== write_output: IVSTP = ", IVSTP
!          PRINT *, "==== write_output: NSCOUV = ", NSCOUV
!          PRINT *, "==== write_output: ICSTP = ", ICSTP
!          PRINT *, "==== write_output: NSCOUC = ", NSCOUC
!          PRINT *, "==== write_output: IPSTP = ", IPSTP
!
!          PRINT *, "==== write_output: IWSTP = ", IWSTP
!          PRINT *, "==== write_output: NSCOUM = ", NSCOUM
!          PRINT *, "==== write_output: IGEP = ", IGEP
!          PRINT *, "==== write_output: NSCOUGE = ", NSCOUGE
!          PRINT *, "==== write_output: IGVP = ", IGVP
!
!          PRINT *, "==== write_output: NSCOUGV = ", NSCOUGV
!          PRINT *, "==== write_output: IGCP = ", IGCP
!          PRINT *, "==== write_output: NSCOUGC = ", NSCOUGC
!          PRINT *, "==== write_output: IGPP = ", IGPP
!          PRINT *, "==== write_output: IGWP = ", IGWP
!          PRINT *, "==== write_output: NSCOUGW = ", NSCOUGW
!          PRINT *,"========================net===================="
!          
!       ENDIF

!     &                    iestp,nscoue,ivstp,nscouv, icstp,
!     &                    nscouc, ipstp, iwstp, nscoum, 
!     &                    igep, nscouge, igvp, nscougv, 
!     &                    igcp, nscougc, igpp, igwp,
!     &                    nscougw,

!      Note that: igcp,nscougc,nscouc,icstp are not yet supported in netCDF
!      IF((myproc.eq.0.and.ncmode.eq.2).or.ncmode.eq.1) THEN
!         WRITE(6,*)"wo: - BEF 67 NPX = ",NPX," NEX = ",NEX
!         do i=1, 10
!             PRINT *, i, "wo: z1 z2 = ", eta11(i,1),eta22(i,1)
!         end do
!         do i=1, 10
!             PRINT *, i, "wo: uu22 = ", uu22(i,1)
!         end do
!         do i=1, 10
!             PRINT *, i, "wo: vv22 = ", vv22(i,1)
!         end do
!         do i=1, 10
!             PRINT *, i, "wo: nodecodenc = ", nodecodenc(i)
!         end do
!         do i=1, 10
!            PRINT *, i, " wo: noffnc = ", noffnc(i)
!         end do
!      ENDIF
!      IF(ncmode.ne. 2. or . (ncmode.eq.2.and. MYPROC.EQ.0)) THEN 
!         WRITE(678,*)"nt = ",nt
!         DO I=1,npx
!            WRITE(678,*) "Store zeta1 ", i, eta11(i,1)
!         END DO
!         DO I=1,npx
!            WRITE(678,*) "Store zeta2 ", i, eta22(i,1)
!         END DO
!         DO I=1,npx
!            WRITE(678,*) "Store u ",i, uu22(i,1)
!         END DO
!         DO I=1,npx
!            WRITE(678,*) "Store v ",i, vv22(i,1)
!         END DO
!         DO I=1,npx
!            WRITE(678,*) "Store nodecodenc ", i, nodecodenc(i)
!         END DO
!         DO I=1,nex
!            WRITE(678,*) "Store noffnc ", i, noffnc(i)
!         END DO
!      ENDIF

       IF ( (MyProc.eq.0) ) PRINT *," call writenc67 "
          call writenc67 (nhstar,abs(hss % lun),timenc_len,npx,timenc,
     &                    nt, nhsinc, uu22,vv22,
     &                    eta11, eta22, etadd, nodecodenc,noffnc,
     &                    imhs, it,
     &                    nciestp,nccoue,ncivstp,nccouv, icstp,
     &                    nscouc, ncipstp, nciwstp, nccoum, 
     &                    ncigep, nccouge, ncigvp, nccougv, 
     &                    igcp, nscougc, ncigpp, ncigwp,
     &                    nccougw,
     &                    rundes, runid, agrid, irtype,
     &                    dtdp, ihot, ics, nolibf, nolifa,
     &                    nolica, nolicat, nwp,
     &                    ncor, ntip, nws, nramp, tau0,
     &                    statim, reftim, rnday, dramp,
     &                    a00, b00, c00, h0,
     &                    slam0, sfea0, cf ,eslm, cori,
     &                    ntif, nbfr, nex, nhy, nm,
     &                    xnc, ync, dp,
     &                    nope, neta, nbdv,nopenc,netanc,
     &                    nvdllnc,
     &                    nbou, nvel, nbvvnc,nbounc,nvelnc,
     &                    nvellnc,ibtypenc,
     &                    ncmode,
     &                    FileFmtMajor, FileFmtMinor, FileFmtRev,
     &                    title, institution, source, history,
     &                    references, comments, host, convention, 
     &                    contact, base_date)

       IF ( (MyProc.eq.0) ) PRINT *," finish call writenc67 "

!      IF(ncmode.ne. 2. or . (ncmode.eq.2.and. MYPROC.EQ.0)) THEN 
!         WRITE(677,*)"nt = ",nt
!         DO I=1,npx
!            WRITE(677,*) "Store zeta1 ", i, eta11(i,1)
!         END DO
!         DO I=1,npx
!            WRITE(677,*) "Store zeta2 ", i, eta22(i,1)
!         END DO
!         DO I=1,npx
!            WRITE(677,*) "Store u ",i, uu22(i,1)
!         END DO
!         DO I=1,npx
!            WRITE(677,*) "Store v ",i, vv22(i,1)
!         END DO
!         DO I=1,npx
!            WRITE(677,*) "Store nodecodenc ", i, nodecodenc(i)
!         END DO
!         DO I=1,nex
!            WRITE(677,*) "Store noffnc ", i, noffnc(i)
!         END DO
!      ENDIF

!      IF((myproc.eq.0.and.ncmode.eq.2).or.ncmode.eq.1) THEN
!         WRITE(6,*)"wo: - AFTER 67 NPX = ",NPX," NEX = ",NEX
!         do i=1, 5
!            PRINT *, i, "wo: z1 z2 = ", eta11(i,1),eta22(i,1)
!         end do
!         do i=1, 5
!            PRINT *, i, "wo: uu22 = ", uu22(i,1)
!         end do
!         do i=1, 5
!            PRINT *, i, "wo: vv22 = ", vv22(i,1)
!         end do
!         do i=1, 5
!            PRINT *, i, "wo: nodecodenc = ", nodecodenc(i)
!         end do
!         do i=1, 5
!            PRINT *, i, " wo: noffnc = ", noffnc(i)
!         end do
!      ENDIF

       ENDIF

#endif

      CASE DEFAULT
         write(ScreenUnit,*) 'The value of NHSTAR=',NHSTAR
         write(ScreenUnit,*) 'is not supported.'
         write(ScreenUnit,*) 'Hotstart file not written.'

      END SELECT

C...    
C...    CLOSE THE HOT START OUTPUT FILE
C...    
      IF(NHSTAR.EQ.1) THEN                                              ! added by mcf 7/14/08
      if (myproc == 0 .or. WRITE_LOCAL_HOT_START_FILES) CLOSE(hss % lun)
      ENDIF                                                             ! end NHSTAR=1 
        
      IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
         WRITE(ScreenUnit,24541) hss % lun,IT,TIME
      ENDIF
      WRITE(16,24541) hss % lun,IT,TIME
24541 FORMAT(1X,'HOT START OUTPUT WRITTEN TO UNIT ',I2,
     &    ' AT TIME STEP = ',I9,' TIME = ',E15.8)
      TEST_HOTSTART = .FALSE.
      IF(hss % lun.EQ.67) THEN
C        jgf45.07 added option to stop ADCIRC after writing hot start file. 
         IF (NHSTAR.EQ.67) THEN
            WRITE(ScreenUnit,*) 'NHSTAR.EQ.67; ADCIRC stopping.'
            WRITE(16,*) 'NHSTAR.EQ.67; ADCIRC stopping.'
            TEST_HOTSTART = .TRUE.
         ENDIF
         hss % lun      = 68
         IF(NHSTAR.EQ.1) hss % filename = 'fort.68'                     ! added by mcf 7/14/08
      ELSE
         IF (NHSTAR.EQ.68) THEN
            WRITE(ScreenUnit,*) 'NHSTAR.EQ.68; ADCIRC stopping.'
            WRITE(16,*) 'NHSTAR.EQ.68; ADCIRC stopping.'
            TEST_HOTSTART = .TRUE.
         ENDIF
         hss % lun=67
         IF(NHSTAR.EQ.1) hss % filename = 'fort.67'                     ! added by mcf 7/14/08
      ENDIF
!      PRINT *, "****TEST_HOTSTART CALL MSG_FINI ******* ", TEST_HOTSTART
      IF (TEST_HOTSTART) THEN
#ifdef CMPI
          CALL MSG_FINI()
#endif               
         STOP
      ENDIF
C-----------------------------------------------------------------------
      END SUBROUTINE writeHotstart
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   W R I T E  W A R N  E L E V
C-----------------------------------------------------------------------
C
C     jgf46.10 Subroutine write an elevation file to fort.69 if any
C     elevation exceeds the WarnElev and the user has set WarnElevDump
C     to .True. (1 in the fort.15 input file). Terminate the run
C     if the number of dumps to the fort.69 file exceeds the user's 
C     limit of WarnElevDumpLimit.
C
C-----------------------------------------------------------------------
      SUBROUTINE WriteWarnElev(time, it)
      USE SIZES
      USE GLOBAL
#ifdef CMPI
      USE MESSENGER
#endif
      IMPLICIT NONE
      REAL(8), intent(in) :: time
      INTEGER, intent(in) :: it
C      
      INTEGER I ! node loop counter
C
C     Open the file and write the header on the first warning. Format
C     the actual header data to make it look more or less like a fort.63
C     file.
      IF (WarnElevDumpCounter.eq.0) THEN
         OPEN(69,FILE=TRIM(LOCALDIR)//'/'//'fort.69')
         WRITE(69,'(A80)') TRIM(RUNDES)
         WRITE(69,6000) WarnElevDumpLimit,NP,DTDP*1.d0,1,1
         CLOSE(69)
      ENDIF
 6000 FORMAT(1x,I10,1x,I10,1x,E15.7,1x,I5,1x,I5)
C
C     Open the file and append the new elevation data.
      OPEN(69,FILE=TRIM(LOCALDIR)//'/'//'fort.69',
     &     ACCESS='SEQUENTIAL',POSITION='APPEND')
      WRITE(69,2120) time,IT
 2120 FORMAT(2X,1pE20.10,5X,I10)
      DO I=1,NP
         IF(NODECODE(I).EQ.1) WRITE(69,2453) I,ETA2(I)
 2453    FORMAT(2X,I8,2X,E15.8)
      ENDDO
      CLOSE(69)
C
      WarnElevDumpCounter = WarnElevDumpCounter+1
C
C     If we have exceeded the user's limit, terminate the run.
      IF (WarnElevDumpCounter.gt.WarnElevDumpLimit) THEN
         WRITE(16,*) 'ERROR: WarnElevDumpLimit Exceeded.'
         WRITE(ScreenUnit,5000) MyProc
 5000    FORMAT('ERROR: WarnElevDumpLimit Exceeded on MyPROC=',i4)
#ifdef CMPI
         CALL MSG_FINI()
#endif               
         STOP
      ENDIF
C
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE WriteWarnElev
C-----------------------------------------------------------------------

