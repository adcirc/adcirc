C----------------------------------------------------------------------
C                     W R I T E   O U T P U T . F
C----------------------------------------------------------------------
C jgf48.03 This file contains all the subroutines that write output to
C files. It was created to consolidate all the output routines that
C were located in other parts of the code, to make it easier to incorporate
C them into globalio, and to make it easier to add new output modes or
C formats in the future.
C
C tcm48.4618 Fixed Bug in writing out global elevation values that are
C dry.
C----------------------------------------------------------------------


C----------------------------------------------------------------------
C       S U B R O U T I N E   W R I T E _ O U T P U T _ 2 D
C----------------------------------------------------------------------
C
C     R.L. 8/22/05 Subroutine to write primary 2D model output not
C     including hotstart and harmonic analysis.
C
C----------------------------------------------------------------------
      SUBROUTINE WRITE_OUTPUT_2D (IT,TimeLoc)
      USE ADCIRC_MOD,ONLY: ADCIRC_TERMINATE
      USE SIZES, ONLY : SZ, INPUTDIR, NBYTE, MNWPROC, MYPROC, MNPROC
C RJW added MYPROC, MNPROC
      USE GLOBAL
      USE GLOBAL_IO, ONLY: storeOne, storeTwo, store63, HEADER_MAX,
     &    open_gbl_file, write_gbl_file, write_gbl_63_skip_drynode,
     &    write_gbl_file_skip_default, collectFullDomainArray,
     &    packOne, unpackOne, packTwo, unpackTwo, pack63, pack83,
     &    open_minmax_file
      USE NodalAttributes, ONLY : OutputTau0, Tau0Var, LoadEleSlopeLim
#ifdef CMPI
      USE WRITER, ONLY : WRITE_GBL_FILE_THROUGH_WRITER,
     $     FLUSH_WRITERS, WRITER_STOREONE,
     $     WRITER_STORETWO, WRITER_STORE63, WRITER_INIT
      USE MESSENGER, ONLY : msg_fini, msg_lbcast
#endif
#ifdef ADCNETCDF
      USE NETCDFIO, ONLY : initNetCDFOutputFile
#endif

      IMPLICIT NONE
      INTEGER, intent(in) :: IT
      REAL(8), intent(in) :: TimeLoc

      INTEGER I                       !local loop counters
      INTEGER NC1, NC2, NC3, NCEle
      INTEGER NM1, NM2, NM3
      INTEGER,SAVE :: countsponge=1

      REAL(SZ) C1, C2, C3
      REAL(SZ) EE1, EE2, EE3
      REAL(SZ) H2, H2N1, H2N2, H2N3
      REAL(SZ) U11, U22, U33
      REAL(SZ) V11, V22, V33
      REAL(SZ) CIC11,CIC22,CIC33  ! v49.64.01 tcm -- added for ice
C
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
      type(OutputDataDescript_t), SAVE :: IceDescript    !tcm v49.64.01 added for ice
      type(OutputDataDescript_t), SAVE :: IceStaDescript !tcm v49.64.01 added for ice stations
      type(OutputDataDescript_t), SAVE :: EtaMaxDescript
      type(OutputDataDescript_t), SAVE :: UMaxDescript
      type(OutputDataDescript_t), SAVE :: WVMaxDescript
      type(OutputDataDescript_t), SAVE :: PrMinDescript
      type(OutputDataDescript_t), SAVE :: RSMaxDescript
      type(OutputDataDescript_t), SAVE :: SpongeDescript
      type(OutputDataDescript_t), SAVE :: ESLDescript
      type(OutputDataDescript_t), SAVE :: BathyDescript   !tcm v50.66.01 added time varying bathy
      type(OutputDataDescript_t), SAVE :: BathyStaDescript !tcm v50.66.01 added time varying bathy stations
!  tcm v50.75 removed ifdef cswan to allow for use whenever nrs=3 or nrs=4
!#ifdef CSWAN
Casey 090302: Added this type for output of radiation stress gradients.
      type(OutputDataDescript_t), SAVE :: RSDescript
!#endif

C     jgf48.03 Full domain arrays used when writing globalio output.
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: ET00_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: UU00_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: VV00_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: TAU0VAR_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: RMP00_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: RMU00_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: RMV00_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: Pr2_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: WVNXOut_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: WVNYOut_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: CC00_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: CICEOut_g(:)  !tcm v49.64.01 added for ice
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: RMICE00_g(:)  !tcm v49.64.01 added for ice
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: ETAMAX_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: UMAX_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: PRMIN_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: WVNOUTMAX_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: RSNMAX_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: sponge_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: ESLONOFF_G(:) !zc added for ESL On/Off output
      REAL(SZ), SAVE, Allocatable, Target :: DPOut_g(:)    !tcm v50.66.01 added time varying bathy
      REAL(SZ), SAVE, Allocatable, Target :: DP00_g(:)  !tcm v50.66.01 added time varying bathy
!  tcm v50.75 removed ifdef cswan to allow for use whenever nrs=3 or nrs=4
!#ifdef CSWAN
Casey 090302: Added these arrays for output of radiation stress gradients.
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: RSNXOut_g(:)
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: RSNYOut_g(:)
!#endif
      LOGICAL, dimension(14) :: netcdf_error
      LOGICAL, SAVE :: FirstCall = .true.

      IF (FirstCall) THEN

         ! fort.61
C         write(16,*) 'Allocating ElevStaDescript' !jgfdebug48.03
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
         ElevStaDescript % ConsiderWetDry       = .TRUE.
         ElevStaDescript % alternate_value      = -99999.0
         ElevStaDescript % field_name           = 'ElevSta'
         IF (ICS.eq.1) THEN
            ElevStaDescript % x_coord              => XEL ! use orig coord
            ElevStaDescript % y_coord              => YEL
         ELSE
            ElevStaDescript % x_coord              => SLEL ! radians
            ElevStaDescript % y_coord              => SFEL
         ENDIF

         ! fort.62
C         write(16,*) 'Allocating VelStaDescript' !jgfdebug48.03
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
         VelStaDescript % ConsiderWetDry       = .FALSE.
         VelStaDescript % alternate_value      = 0.0
         VelStaDescript % field_name           = 'VelSta'
         IF (ICS.eq.1) THEN
            VelStaDescript % x_coord              => XEV ! use orig coord
            VelStaDescript % y_coord              => YEV
         ELSE
            VelStaDescript % x_coord              => SLEV ! radians
            VelStaDescript % y_coord              => SFEV
         ENDIF

         ! fort.63
         ElevDescript % specifier            = NOUTGE
         ElevDescript % initial_value        = 0.0
         ElevDescript % num_items_per_record = 1
         ElevDescript % num_fd_records       = NP_G
         ElevDescript % num_records_this     = NP
         ElevDescript % imap                 => NODES_LG
         ElevDescript % array                => ETA2
         ElevDescript % array_g              => ETA2_g
         ElevDescript % ConsiderWetDry       = .TRUE.
         ElevDescript % alternate_value      = -99999.0
         ElevDescript % field_name           = 'Elev'

         ! tau0
C         write(16,*) 'Allocating Tau0Descript' !jgfdebug48.03
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
         Tau0Descript % ConsiderWetDry       = .FALSE.
         Tau0Descript % alternate_value      = 0.0
         Tau0Descript % field_name           = 'Tau0'

         ! sponge layer
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(sponge_g(NP_G))
         ENDIF
         SpongeDescript % specifier            = 1  ! hard-code to ascii text
         SpongeDescript % initial_value        = 0.0
         SpongeDescript % num_items_per_record = 1
         SpongeDescript % num_fd_records       = NP_G
         SpongeDescript % num_records_this     = NP
         SpongeDescript % imap                 => NODES_LG
         SpongeDescript % array                => sponge
         SpongeDescript % array_g              => sponge_g

         ! fort.64
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
         VelDescript % ConsiderWetDry       = .FALSE.
         VelDescript % alternate_value      = 0.0
         VelDescript % field_name           = 'Vel'

         ! fort.71
C         write(16,*) 'Allocating PrStaDescript' !jgfdebug48.03
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
         PrStaDescript % ConsiderWetDry       = .FALSE.
         PrStaDescript % alternate_value      = 0.0
         PrStaDescript % field_name           = 'PrSta'
         IF (ICS.eq.1) THEN
            PrStaDescript % x_coord           => XEM ! use orig coord
            PrStaDescript % y_coord           => YEM ! use orig coord
         ELSE
            PrStaDescript % x_coord           => SLEM ! radians
            PrStaDescript % y_coord           => SFEM 
         ENDIF

         ! fort.72
C         write(16,*) 'Allocating WindVelStaDescript' !jgfdebug48.03
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
         WindVelStaDescript % ConsiderWetDry       = .FALSE.
         WindVelStaDescript % alternate_value      = 0.0
         WindVelStaDescript % field_name           = 'WindVelSta'

         ! fort.73
C         write(16,*) 'Allocating PrDescript' !jgfdebug48.03
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
         PrDescript % ConsiderWetDry       = .FALSE.
         PrDescript % alternate_value      = 0.0
         PrDescript % field_name           = 'Pr'
         IF (ICS.eq.1) THEN
            WindVelStaDescript % x_coord           => XEM ! use orig coord
            WindVelStaDescript % y_coord           => YEM ! use orig coord
         ELSE
            WindVelStaDescript % x_coord           => SLEM ! use orig coord
            WindVelStaDescript % y_coord           => SFEM ! use orig coord
         ENDIF

         ! fort.74
C         write(16,*) 'Allocating WindVelDescript' !jgfdebug48.03
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
         WindVelDescript % ConsiderWetDry       = .FALSE.
         WindVelDescript % alternate_value      = 0.0
         WindVelDescript % field_name           = 'WindVel'

         ! fort.75  !tcm v50.66.01 Added for Time Varying Bathymetry
C         write(16,*) 'Allocating BathyStaDescript'
         IF (NDDT.NE.0 ) then
            IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
               ALLOCATE(DP00_g(NSTAE_G))
            ENDIF
            BathyStaDescript % specifier            = NOUTE
            BathyStaDescript % initial_value        = 0.0
            BathyStaDescript % num_items_per_record = 1
            BathyStaDescript % num_fd_records       = NSTAE_G
            BathyStaDescript % num_records_this     = NSTAE
            BathyStaDescript % imap                 => IMAP_STAE_LG
            BathyStaDescript % array                => DP00
            BathyStaDescript % array_g              => DP00_g
            BathyStaDescript % ConsiderWetDry       = .FALSE.
            BathyStaDescript % alternate_value      = 0.d0
            BathyStaDescript % field_name           = 'BathySta'
            IF (ICS.eq.2) THEN
               BathyStaDescript % x_coord           => SLEL
               BathyStaDescript % y_coord           => SFEL
            ELSE
               BathyStaDescript % x_coord           => XEL
               BathyStaDescript % y_coord           => YEL
            ENDIF
         ENDIF

         ! fort.76  !tcm v50.66.01  Added for Time Varying Bathymetry
C         write(16,*) 'Allocating BathyDescript'
         IF (NDDT.NE.0) then
            IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
               ALLOCATE(DPOUT_g(NP_G))
            ENDIF
            BathyDescript % specifier            = NOUTGE
            BathyDescript % initial_value        = 0.0
            BathyDescript % num_items_per_record = 1
            BathyDescript % num_fd_records       = NP_G
            BathyDescript % num_records_this     = NP
            BathyDescript % imap                 => NODES_LG
            BathyDescript % array                => DP
            BathyDescript % array_g              => DPOut_g
            BathyDescript % ConsiderWetDry       = .FALSE.
            BathyDescript % alternate_value      = 0.0
            BathyDescript % field_name           = 'Bathy'
         ENDIF

         ! fort.93
C         write(16,*) 'Allocating IceDescript' !tcm v49.64.01 added for ice
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(CICEOUT_g(NP_G))
         ENDIF
         IceDescript % specifier            = NOUTGW
         IceDescript % initial_value        = 0.0
         IceDescript % num_items_per_record = 1
         IceDescript % num_fd_records       = NP_G
         IceDescript % num_records_this     = NP
         IceDescript % imap                 => NODES_LG
         IceDescript % array                => CICEOUT
         IceDescript % array_g              => CICEOUT_g
         IceDescript % ConsiderWetDry       = .FALSE.
         IceDescript % alternate_value      = 0.0
         IceDescript % field_name           = 'Ice'

C... v49.64.01 tcm -- added for ice stations
         ! fort.91
C         write(16,*) 'Allocating IceStaDescript'
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(RMICE00_g(NSTAM_G))
         ENDIF
         IceStaDescript % specifier            = NOUTM
         IceStaDescript % initial_value        = 0.0
         IceStaDescript % num_items_per_record = 1
         IceStaDescript % num_fd_records       = NSTAM_G
         IceStaDescript % num_records_this     = NSTAM
         IceStaDescript % imap                 => IMAP_STAM_LG
         IceStaDescript % array                => RMICE00
         IceStaDescript % array_g              => RMICE00_g
         IceStaDescript % ConsiderWetDry       = .FALSE.
         IceStaDescript % alternate_value      = 0.0
         IceStaDescript % field_name           = 'IceSta'
         IF (ICS.eq.1) THEN
            IceStaDescript % x_coord           => XEM ! use orig coord
            IceStaDescript % y_coord           => YEM
         ELSE
            IceStaDescript % x_coord           => SLEM ! radians
            IceStaDescript % y_coord           => SFEM
         ENDIF

!  tcm v50.75 removed ifdef cswan to allow for use whenever nrs=3 or nrs=4
         IF ((ABS(NRS).EQ.3).OR.(ABS(NRS).EQ.4)) THEN
!#ifdef CSWAN   !tcm
Casey 090302: Copied these lines for the rads.64 file.
         ! rads.64
C         write(16,*) 'Allocating RSDescript' !jgfdebug48.03
            IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
               ALLOCATE(RSNXOUT_g(NP_G))
               ALLOCATE(RSNYOUT_g(NP_G))
            ENDIF
C........Radiation Stress
            RSDescript % specifier            = NOUTGW
            RSDescript % initial_value        = 0.0
            RSDescript % num_items_per_record = 2
            RSDescript % num_fd_records       = NP_G
            RSDescript % num_records_this     = NP
            RSDescript % imap                 => NODES_LG
            RSDescript % array                => RSNXOUT
            RSDescript % array_g              => RSNXOUT_g
            RSDescript % array2               => RSNYOUT
            RSDescript % array2_g             => RSNYOUT_g
            RSDescript % alternate_value      = -99999.0
            RSDescript % field_name           = 'RadStress'
            RSDescript % file_name            = "rads.64"
         ENDIF
!#endif

         ! fort.81
C         write(16,*) 'Allocating ConcStaDescript' !jgfdebug48.03
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
         ConcStaDescript % ConsiderWetDry       = .FALSE.
         ConcStaDescript % alternate_value      = 0.0
         ConcStaDescript % field_name           = 'ConcSta'
         IF (ICS.eq.1) THEN
            ConcStaDescript % x_coord           => XEC ! use orig coord
            ConcStaDescript % y_coord           => YEC
         ELSE
            ConcStaDescript % x_coord           => SLEC ! radians
            ConcStaDescript % y_coord           => SFEC
         ENDIF

         ! fort.83
C         write(16,*) 'Allocating ConcDescript' !jgfdebug48.03
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
         ConcDescript % ConsiderWetDry       = .FALSE.
         ConcDescript % alternate_value      = 0.0
         ConcDescript % field_name           = 'Conc'

C        maxele.63
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(ETAMAX_g(NP_G))
         ENDIF
         EtaMaxDescript % specifier            = NOUTGE
         EtaMaxDescript % initial_value        = 0.0
         EtaMaxDescript % num_items_per_record = 1
         EtaMaxDescript % num_fd_records       = NP_G
         EtaMaxDescript % num_records_this     = NP
         EtaMaxDescript % imap                 => NODES_LG
         EtaMaxDescript % array                => ETAMAX
         EtaMaxDescript % array_g              => ETAMAX_G
         EtaMaxDescript % hotstart             => HOT_ETAMAX
         EtaMaxDescript % hotstart_g           => HOT_ETAMAX_G
         EtaMaxDescript % ConsiderWetDry       = .FALSE.
         EtaMaxDescript % alternate_value      = 0.0
         EtaMaxDescript % field_name           = 'EtaMax'
         EtaMaxDescript % file_name            = 'maxele.63'
C        maxvel.63
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(UMAX_g(NP_G))
         ENDIF
         UMaxDescript % specifier              = NOUTGV
         UMaxDescript % initial_value          = 0.0
         UMaxDescript % num_items_per_record   = 1
         UMaxDescript % num_fd_records         = NP_G
         UMaxDescript % num_records_this       = NP
         UMaxDescript % imap                   => NODES_LG
         UMaxDescript % array                  => UMAX
         UMaxDescript % array_g                => UMAX_G
         UMaxDescript % hotstart               => HOT_UMAX
         UMaxDescript % hotstart_g             => HOT_UMAX_G
         UMaxDescript % ConsiderWetDry         = .FALSE.
         UMaxDescript % alternate_value        = 0.0
         UMaxDescript % field_name             = 'UMax'
         UMaxDescript % file_name              = 'maxvel.63'
C        prmin.63
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(PRMIN_g(NP_G))
         ENDIF
         PrMinDescript % specifier             = NOUTGW
         PrMinDescript % initial_value         = 0.0
         PrMinDescript % num_items_per_record  = 1
         PrMinDescript % num_fd_records        = NP_G
         PrMinDescript % num_records_this      = NP
         PrMinDescript % imap                  => NODES_LG
         PrMinDescript % array                 => PRMIN
         PrMinDescript % array_g               => PRMIN_G
         PrMinDescript % hotstart              => HOT_PRMIN
         PrMinDescript % hotstart_g            => HOT_PRMIN_G
         PrMinDescript % ConsiderWetDry        = .FALSE.
         PrMinDescript % alternate_value       = 0.0
         PrMinDescript % field_name            = 'PrMin'
         PrMinDescript % file_name             = 'minpr.63'
C        maxwvel.63
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(WVNOUTMAX_g(NP_G))
         ENDIF
         WVMaxDescript % specifier             = NOUTGW
         WVMaxDescript % initial_value         = 0.0
         WVMaxDescript % num_items_per_record  = 1
         WVMaxDescript % num_fd_records        = NP_G
         WVMaxDescript % num_records_this      = NP
         WVMaxDescript % imap                  => NODES_LG
         WVMaxDescript % array                 => WVNOUTMAX
         WVMaxDescript % array_g               => WVNOUTMAX_G
         WVMaxDescript % hotstart              => HOT_WVNOUTMAX
         WVMaxDescript % hotstart_g            => HOT_WVNOUTMAX_G
         WVMaxDescript % ConsiderWetDry        = .FALSE.
         WVMaxDescript % alternate_value       = 0.0
         WVMaxDescript % field_name            = 'WVMax'
         WVMaxDescript % file_name             = 'maxwvel.63'
C        maxrs.63
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(RSNMAX_g(NP_G))
         ENDIF
         RSMaxDescript % specifier             = NOUTGW
         RSMaxDescript % initial_value         = 0.0
         RSMaxDescript % num_items_per_record  = 1
         RSMaxDescript % num_fd_records        = NP_G
         RSMaxDescript % num_records_this      = NP
         RSMaxDescript % imap                  => NODES_LG
         RSMaxDescript % array                 => RSNMAX
         RSMaxDescript % array_g               => RSNMAX_G
         RSMaxDescript % hotstart              => HOT_RSNMAX
         RSMaxDescript % hotstart_g            => HOT_RSNMAX_G
         RSMaxDescript % ConsiderWetDry        = .FALSE.
         RSMaxDescript % alternate_value       = 0.0
         RSMaxDescript % field_name            = 'RSMax'
         RSMaxDescript % file_name             = 'maxrs.63'
C        ESLNodes.63
         IF(LoadEleSlopeLim)THEN
            IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
                ALLOCATE(ESLONOFF_G(NP_G))
            ENDIF
            ESLDescript % specifier            = 1
            ESLDescript % initial_value        = 0.0
            ESLDescript % num_items_per_record = 1
            ESLDescript % num_fd_records       = NP_G
            ESLDescript % num_records_this     = NP
            ESLDescript % imap                 => NODES_LG
            ESLDescript % array                => ESLONOFF
            ESLDescript % array_g              => ESLONOFF_G
            ESLDescript % hotstart             => HOT_ESLONOFF
            ESLDescript % hotstart_g           => HOT_ESLONOFF_G
            ESLDescript % ConsiderWetDry       = .FALSE.
            ESLDescript % alternate_value      = 0.0
            ESLDescript % field_name           = 'ESLActive'
            ESLDescript % file_name            = 'ESLNodes.63'
         ENDIF
C
#ifdef ADCNETCDF
         netcdf_error = .false.
         ! jgf49.43.11: if netcdf is not required for a particular output file,
         ! the associated subroutine call will do nothing and just return
         IF (myProc.eq.0) THEN
            CALL initNetCDFOutputFile(61,ElevStaDescript, netcdf_error(1))
            CALL initNetCDFOutputFile(62,VelStaDescript, netcdf_error(2))
            CALL initNetCDFOutputFile(63,ElevDescript, netcdf_error(3))
            CALL initNetCDFOutputFile(64,VelDescript, netcdf_error(4))
            CALL initNetCDFOutputFile(71,PrStaDescript, netcdf_error(5))
            CALL initNetCDFOutputFile(72,WindVelStaDescript,
     &                                netcdf_error(6))
            CALL initNetCDFOutputFile(73,PrDescript, netcdf_error(7))
            CALL initNetCDFOutputFile(74,WindVelDescript, netcdf_error(8))
            CALL initNetCDFOutputFile(311,EtaMaxDescript,netcdf_error(10))
            CALL initNetCDFOutputFile(312,UMaxDescript,netcdf_error(11))
            CALL initNetCDFOutputFile(313,PrMinDescript,netcdf_error(12))
            CALL initNetCDFOutputFile(314,WVMaxDescript,netcdf_error(13))
            CALL initNetCDFOutputFile(315,RSMaxDescript,netcdf_error(14))
            IF(ABS(NRS).EQ.3.OR.ABS(NRS).EQ.4)THEN
                CALL initNetCDFOutputFile(164,RSDescript, netcdf_error(9))
            ENDIF    
         ENDIF
Cobell - We no longer need this bcast if we set subdomainFatalError first
C        CALL msg_lbcast(netcdf_error,14)
         DO i=1,14
            IF (netcdf_error(i).eqv..true.) THEN
               CALL ADCIRC_TERMINATE()               
               STOP
            ENDIF
         ENDDO
#endif
C
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
Casey 090302: Added file name.
            CALL writeOutArray(61, TimeLoc, IT, ElevStaDescript,
     &                         packOne, unpackOne, IESTP,
     &                         'fort.61        ')

c.....     tcm v50.66.01 added for time varying bathymetry
c.....     save bathymetry recording station data to file
            IF (NDDT.NE.0) THEN
               DO I=1,NSTAE
                  EE1=DP(NM(NNE(I),1))
                  EE2=DP(NM(NNE(I),2))
                  EE3=DP(NM(NNE(I),3))
                  DP00(I)=EE1*STAIE1(I)+EE2*STAIE2(I)
     &                                 +EE3*STAIE3(I)
               ENDDO
               CALL writeOutArray(75,TimeLoc,IT,BathyStaDescript,
     &                         packOne, unpackOne, IBSTP,
     &                         'fort.75        ')
            ENDIF

C            write(16,*) 'Finished with fort.61 output routine.' !jgfdebug48.03
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
C            write(16,*) 'About to write fort.62 file.' !jgfdebug48.03
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
Casey 090302: Added file name.
            CALL writeOutArray(62, TimeLoc, IT, VelStaDescript,
     &                         packTwo, unpackTwo, IVSTP,
     &                         'fort.62        ')

C            write(16,*) 'Finished with fort.62 output routine.' !jgfdebug48.03
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
               IF(NCELE.EQ.0) CC00(I)=-99999.0
            END DO
C
C           save concentration recording station data to file
Casey 090302: Added file name.
            CALL writeOutArray(81, TimeLoc, IT, ConcStaDescript,
     &                         packOne, unpackOne, ICSTP,
     &                         'fort.81        ')

C            write(16,*) 'Finished with fort.81 output routine.' !jgfdebug48.03
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
C            write(16,*) 'About to write fort.71,72 files.' !jgfdebug48.03
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
C... v49.64.01 tcm -- added support for ice station output
C           IF using ICE fields, collect ice station data fort.91
            IF (NCICE.NE.0) THEN
               DO I = 1,NSTAM
                  NM1=NM(NNM(I),1)
                  NM2=NM(NNM(I),2)
                  NM3=NM(NNM(I),3)
                  CIC11 = CICEOUT(NM1)
                  CIC22 = CICEOUT(NM2)
                  CIC33 = CICEOUT(NM3)
                  RMICE00(I) = CIC11*STAIM1(I) +
     &                         CIC22*STAIM2(I)+CIC33*STAIM3(I)
               END DO
            ENDIF
C
C           save atmospheric pressure recording station data to file
Casey 090302: Added file name.
            CALL writeOutArray(71, TimeLoc, IT, PrStaDescript,
     &                         packOne, unpackOne, IPSTP,
     &                         'fort.71        ')
C
C           save wind velocity recording station data to file
Casey 090302: Added file name.
            CALL writeOutArray(72, TimeLoc, IT, WindVelStaDescript,
     &                         packTwo, unpackTwo, IWSTP,
     &                         'fort.72        ')
c  V49.64.01 TCM -- ADDED ICE STATIONS
C           save ice field recording station data to file
            IF(NCICE.NE.0) THEN
               CALL writeOutArray(91, TimeLoc, IT, IceStaDescript,
     &                         packOne, unpackOne, IICESTP,
     &                         'fort.91        ')
            ENDIF

            NSCOUM=0
            NCCOUM=NSCOUM
         ENDIF
      ENDIF
C...   OUTPUT GLOBAL ELEVATION DATA IF NOUTGE<>0 AND THE
C.... TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
C...
      IF(NOUTGE.NE.0) THEN
         IF((IT.GT.NTCYSGE).AND.(IT.LE.NTCYFGE)) THEN
            NSCOUGE=NSCOUGE+1
            NCCOUGE=NSCOUGE
         ENDIF
         IF(NSCOUGE.EQ.NSPOOLGE) THEN
!         IF(MYPROC.EQ.0) PRINT *, "W: B-write ", NSCOUGE
C
Casey 090302: Added outer IF statement.
            IF((ABS(NOUTGE).EQ.1).OR.
     &         (ABS(NOUTGE).EQ.2).OR.
     &         (ABS(NOUTGE).EQ.3).OR.
     &         (ABS(NOUTGE).EQ.5)) THEN
C           save fulldomain elevation data to file

C....tcm v50.27 added using writer processors for full domain output
               if(ABS(NOUTGE).eq.1.OR.ABS(NOUTGE).EQ.3.OR.
     &            ABS(NOUTGE).EQ.5) then
#ifdef CMPI
Casey 090302: Changed descript to ElevDescript.
                  IF(MNWPROC.GT.0) THEN
                     CALL WRITE_GBL_FILE_THROUGH_WRITER
     $              (63,TRIM(GLOBALDIR)//'/'//'fort.63',ElevDescript,TimeLoc,it,
     $                 writer_store63, -99999.D0)

                     !tcm v50.66.01 added for time varying bathymetry
                   if (nddt.ne.0) call WRITE_GBL_FILE_THROUGH_WRITER
     $              (76,TRIM(GLOBALDIR)//'/'//'fort.76',BathyDescript,TimeLoc,it,
     $                 writer_storeOne, 0.D0)

                  ELSE
                     CALL writeOutArray(63, TimeLoc, IT, ElevDescript,
     &                            pack63, unpackOne, IGEP,
     &                            'fort.63        ')
                     !tcm v50.66.01 added for time varying bathymetry
                     if (nddt.ne.0) CALL writeOutArray(76, TimeLoc, IT,
     &                      BathyDescript,packOne, unpackOne, IGBP,
     &                            'fort.76        ')

                  ENDIF
#else
                  CALL writeOutArray(63, TimeLoc, IT, ElevDescript,
     &                            pack63, unpackOne, IGEP,
     &                            'fort.63        ')

                  !tcm v50.66.01 added for time varying bathymetry
                  if (nddt.ne.0) CALL writeOutArray(76, TimeLoc, IT,
     &                   BathyDescript,packOne, unpackOne, IGBP,
     &                            'fort.76        ')
#endif
               else
Casey 090302: Added file name.
                  CALL writeOutArray(63, TimeLoc, IT, ElevDescript,
     &                            pack63, unpackOne, IGEP,
     &                            'fort.63        ')

                  !tcm v50.66.01 added for time varying bathymetry
                  if (nddt.ne.0) CALL writeOutArray(76, TimeLoc, IT,
     &                   BathyDescript,packOne, unpackOne, IGBP,
     &                            'fort.76        ')

               endif
C
C              jgf47.06: Write out tau0 if it is time varying and the user has
C              requested it with tau0 = -4 in fort.15.
               IF (OUTPUTTAU0) THEN
Casey 090302: Added file name.
                  CALL writeOutArray(90, TimeLoc, IT, Tau0Descript,
     &                               packOne, unpackOne, IGEP,
     &                               'fort.90        ')
               ENDIF
C
C           kmd48.33bc: Write out sponge layer so user can see which nodes
C           are included in the sponge layer. Note only outputs if spongedist
C           given in fort.15 is not equal to zero.
               IF ((OUTPUTSPONGE).AND.(countsponge.eq.1)) THEN
                  CALL writeOutArray(92, TimeLoc, IT, SpongeDescript,
     &                               packOne, unpackOne, IGEP,
     &                               'fort.92        ')
                  countsponge=countsponge + 1 ! only need to output once as it does not change with time
               ENDIF

C          jgf47.08 Made sparse output an option (if NOUTGE=-4 or 4 in fort.15)
            ELSEIF(ABS(NOUTGE).EQ.4) THEN
C              jgf47.08 Create compact output
C              jgf48.03 TODO: Make this globalio compliant.
C              kmd - updated
#ifdef CMPI
Casey 090302: Changed descript to ElevDescript.
              IF(MNWPROC.GT.0) THEN
                 CALL WRITE_GBL_FILE_THROUGH_WRITER
     $              (63,TRIM(GLOBALDIR)//'/'//'fort.63',ElevDescript,TimeLoc,it,
     $                 writer_store63, -99999.D0)

                  !tcm v50.66.01 added for time varying bathymetry
                  !this will be a full record not compact
                 if (nddt.ne.0) CALL WRITE_GBL_FILE_THROUGH_WRITER
     $              (63,TRIM(GLOBALDIR)//'/'//'fort.76',BathyDescript,TimeLoc,it,
     $                 writer_storeOne, 0.D0)

              ELSE
                 CALL WRITE_GBL_63_SKIP_DRYNODE  ! sb 11/10/2006
     $                (63, TRIM(GLOBALDIR)//'/'//'fort.63',
     $                 ElevDescript, TimeLoc,it, store63)

                  !tcm v50.66.01 added for time varying bathymetry
                  ! note that this will be a full record not compact
                  IF (NDDT.NE.0) CALL writeOutArray(76, TimeLoc, IT,
     &                   BathyDescript,packOne, unpackOne, IGBP,
     &                         'fort.76        ')

              ENDIF
#else
              CALL WRITE_GBL_63_SKIP_DRYNODE  ! sb 11/10/2006
     $            (63, TRIM(GLOBALDIR)//'/'//'fort.63',
     &               ElevDescript, TimeLoc,it, store63)

              !tcm v50.66.01 added for time varying bathymetry
              ! note that this will be a full record not compact
              IF (NDDT.NE.0) CALL writeOutArray(76, TimeLoc, IT,
     &                BathyDescript,packOne, unpackOne, IGBP,
     &                      'fort.76        ')

#endif
              IGEP=IGEP+1+NP
              IGBP=IGBP+1+NP
           ENDIF
C           write(16,*) 'Finished with fort.63 output routine.' !jgfdebug48.03
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
Casey 090302: Added outer IF statement.
            IF((ABS(NOUTGV).EQ.1).OR.
     &         (ABS(NOUTGV).EQ.2).OR.
     &         (ABS(NOUTGV).EQ.5).OR.
     &         (ABS(NOUTGV).EQ.3)) THEN

C....tcm v50.27 added using writer processors for full domain output
               if(ABS(NOUTGV).eq.1.OR.ABS(NOUTGV).eq.3.OR.
     &            ABS(NOUTGV).EQ.5) then
#ifdef CMPI
Casey 090302: Changed descript to VelDescript.
                  IF(MNWPROC.GT.0) THEN
                     CALL WRITE_GBL_FILE_THROUGH_WRITER
     $                (64,TRIM(GLOBALDIR)//'/'//'fort.64',VelDescript,TimeLoc,it,
     $                 writer_storeTwo, 0.D0)
                  ELSE
                     CALL writeOutArray(64, TimeLoc, IT, VelDescript,
     &                         packTwo, unpackTwo, IGVP,
     &                         'fort.64        ')
                  ENDIF
#else
                  CALL writeOutArray(64, TimeLoc, IT, VelDescript,
     &                         packTwo, unpackTwo, IGVP,
     &                         'fort.64        ')
#endif
               else
Casey 090302: Added file name.
                  CALL writeOutArray(64, TimeLoc, IT, VelDescript,
     &                         packTwo, unpackTwo, IGVP,
     &                         'fort.64        ')
               endif
            ELSEIF (ABS(NOUTGV).EQ.4) THEN
            ! jgf48.03 TODO: Make this globalio compliant.
            ! kmd - updated
#ifdef CMPI
Casey 090302: Changed descript to VelDescript.
                  IF(MNWPROC.GT.0) THEN
                    CALL WRITE_GBL_FILE_THROUGH_WRITER
     $                (64,TRIM(GLOBALDIR)//'/'//'fort.64',VelDescript,TimeLoc,it,
     $                 writer_storeTwo, 0.D0)
                  ELSE
                    CALL WRITE_GBL_FILE_SKIP_DEFAULT ! sb 11/10/2006
     $                   (64, TRIM(GLOBALDIR)//'/'//'fort.64',
     $                   VelDescript, TimeLoc, it, storeTwo,0.d0)
                  ENDIF
#else
                  CALL WRITE_GBL_FILE_SKIP_DEFAULT ! sb 11/10/2006
     $                 (64, TRIM(GLOBALDIR)//'/'//'fort.64',
     $                 VelDescript, TimeLoc, it, storeTwo,0.d0)
#endif
            ENDIF
            NSCOUGV=0  ! kmd49 - added this back into the code
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
Casey 090302: Added outer IF statement.
            IF((ABS(NOUTGW).EQ.1).OR.
     &         (ABS(NOUTGW).EQ.2).OR.
     &         (ABS(NOUTGW).EQ.5).OR.
     &         (ABS(NOUTGW).EQ.3)) THEN
C           save atmospheric pressure recording station data to file
C....tcm v50.27 added using writer processors for full domain output
               IF(ABS(NOUTGW).eq.1.OR.ABS(NOUTGW).EQ.3.OR.
     &            ABS(NOUTGW).EQ.5) then
#ifdef CMPI
Casey 090302: Changed descript to VelDescript.
                  IF(MNWPROC.GT.0) THEN
                     CALL WRITE_GBL_FILE_THROUGH_WRITER
     $                  (73,TRIM(GLOBALDIR)//'/'//'fort.73',PrDescript,TimeLoc,
     $                   it,writer_storeOne, 0.D0)
                     CALL WRITE_GBL_FILE_THROUGH_WRITER
     $                   (74,TRIM(GLOBALDIR)//'/'//'fort.74',
     &                   WindVelDescript,TimeLoc,it,writer_storeTwo, 0.D0)
                     IF(NCICE.NE.0) then
                        CALL WRITE_GBL_FILE_THROUGH_WRITER
     $                       (93,TRIM(GLOBALDIR)//'/'//'fort.93',
     $                       IceDescript,TimeLoc,it,writer_storeOne, 0.D0)
                     ENDIF
!  tcm v50.75 removed ifdef cswan to allow for use whenever nrs=3 or nrs=4
                     IF ((ABS(NRS).EQ.3).OR.(ABS(NRS).EQ.4)) THEN
!#ifdef CSWAN
                        CALL WRITE_GBL_FILE_THROUGH_WRITER
     &                      (164,TRIM(GLOBALDIR)//'/'//'rads.64',RSDescript,
     &                       TimeLoc,it,writer_storeTwo, 0.D0)
                     ENDIF
!#endif
                  ELSE  !no writer processors

                     CALL writeOutArray(73, TimeLoc, IT, PrDescript,
     &                            packOne, unpackOne, IGPP,
     &                            'fort.73        ')
C           save wind velocity global data to file
Casey 090302: Added file name.
                     CALL writeOutArray(74, TimeLoc, IT, WindVelDescript,
     &                            packTwo, unpackTwo, IGWP,
     &                            'fort.74        ')
                     IF (NCICE.NE.0) THEN
C TCM V49.64.01 -- ADDED ICE SUPPORT
C           save global ICE data to file
                        CALL writeOutArray(93, TimeLoc, IT, IceDescript,
     &                         packOne, unpackOne, IGIP,
     &                         'fort.93        ')
                     ENDIF
!  tcm v50.75 removed ifdef cswan to allow for use whenever nrs=3 or nrs=4
                     IF ((ABS(NRS).EQ.3).OR.(ABS(NRS).EQ.4)) THEN
!#ifdef CSWAN
Casey 090302: Added these lines for output of radiation stresses.
                        CALL writeOutArray(164,TimeLoc, IT, RSDescript,
     &                             packTwo, unpackTwo, IGRadS,
     &                            'rads.64        ')
                     ENDIF
!#endif
                  ENDIF  !test for writer processors

#else
       !Not parallel
                  CALL writeOutArray(73, TimeLoc, IT, PrDescript,
     &                            packOne, unpackOne, IGPP,
     &                            'fort.73        ')
C           save wind velocity global data to file
Casey 090302: Added file name.
                  CALL writeOutArray(74, TimeLoc, IT, WindVelDescript,
     &                            packTwo, unpackTwo, IGWP,
     &                            'fort.74        ')
                  IF (NCICE.NE.0) THEN
C TCM V49.64.01 -- ADDED ICE SUPPORT
C           save global ICE data to file
                     CALL writeOutArray(93, TimeLoc, IT, IceDescript,
     &                         packOne, unpackOne, IGIP,
     &                         'fort.93        ')
                  ENDIF
!  tcm v50.75 removed ifdef cswan to allow for use whenever nrs=3 or nrs=4
                  IF ((ABS(NRS).EQ.3).OR.(ABS(NRS).EQ.4)) THEN
!#ifdef CSWAN
Casey 090302: Added these lines for output of radiation stresses.
                     CALL writeOutArray(164,TimeLoc, IT, RSDescript,
     &                               packTwo, unpackTwo, IGRadS,
     &                               'rads.64        ')
                  ENDIF
!#endif

#endif
               ELSE  ! abs(noutgw) == 2 or 3

Casey 090302: Added file name.
                  CALL writeOutArray(73, TimeLoc, IT, PrDescript,
     &                            packOne, unpackOne, IGPP,
     &                            'fort.73        ')
C           save wind velocity global data to file
Casey 090302: Added file name.
                  CALL writeOutArray(74, TimeLoc, IT, WindVelDescript,
     &                            packTwo, unpackTwo, IGWP,
     &                            'fort.74        ')
                  IF (NCICE.NE.0) THEN
C TCM V49.64.01 -- ADDED ICE SUPPORT
C           save global ICE data to file
                     CALL writeOutArray(93, TimeLoc, IT, IceDescript,
     &                         packOne, unpackOne, IGIP,
     &                         'fort.93        ')
                  ENDIF
!  tcm v50.75 removed ifdef cswan to allow for use whenever nrs=3 or nrs=4
                  IF ((ABS(NRS).EQ.3).OR.(ABS(NRS).EQ.4)) THEN
!#ifdef CSWAN
Casey 090302: Added these lines for output of radiation stresses.
                     CALL writeOutArray(164,TimeLoc, IT, RSDescript,
     &                            packTwo, unpackTwo, IGRadS,
     &                            'rads.64        ')
                  ENDIF
!#endif
               ENDIF  !test for abs(noutgw) == 1

C           jgf48.03 TODO: Make this globalio compliant.
            ELSEIF (ABS(NOUTGW).EQ.4) THEN
#ifdef CMPI
Casey 090302: Changed descript to PrDescript
               IF(MNWPROC.GT.0) THEN
                 CALL WRITE_GBL_FILE_THROUGH_WRITER
     $             (73,TRIM(GLOBALDIR)//'/'//'fort.73',PrDescript,TimeLoc,it,
     $              writer_storeOne, 0.D0)
               ELSE
                 CALL WRITE_GBL_FILE_SKIP_DEFAULT
     $             (73, TRIM(GLOBALDIR)//'/'//'fort.73',
     $              PrDescript, TimeLoc, it, storeOne,0.d0)
               ENDIF
#else
               CALL WRITE_GBL_FILE_SKIP_DEFAULT
     $            (73, TRIM(GLOBALDIR)//'/'//'fort.73',
     $               PrDescript, TimeLoc, it, storeOne,0.d0)
#endif
C           ENDIF
C           jgf48.03 TODO: Make this globalio compliant.
C           IF (ABS(NOUTGW).EQ.4) THEN
#ifdef CMPI
Casey 090302: Changed descript to WindVelDescript
               IF(MNWPROC.GT.0) THEN
                 CALL WRITE_GBL_FILE_THROUGH_WRITER
     $                   (74,TRIM(GLOBALDIR)//'/'//'fort.74',
     &                   WindVelDescript,TimeLoc,it,
     $                   writer_storeTwo, 0.D0)
               ELSE
                 CALL WRITE_GBL_FILE_SKIP_DEFAULT
     $                   (74, TRIM(GLOBALDIR)//'/'//'fort.74',
     $                   WindVelDescript, TimeLoc, it, storeTwo,0.d0)
               ENDIF
#else
               CALL WRITE_GBL_FILE_SKIP_DEFAULT
     $              (74, TRIM(GLOBALDIR)//'/'//'fort.74',
     $              WindVelDescript, TimeLoc, it, storeTwo,0.d0)
#endif

C  ! tcm v49.64.01 ADDED for ICE
C           jgf48.03 TODO: Make this globalio compliant.
               IF (NCICE.NE.0) THEN
#ifdef CMPI
                  IF(MNWPROC.GT.0) THEN
                    CALL WRITE_GBL_FILE_THROUGH_WRITER
     $                (93,TRIM(GLOBALDIR)//'/'//'fort.93',
     $                IceDescript,TimeLoc,it,writer_storeOne, 0.D0)
                  ELSE
                     CALL WRITE_GBL_FILE_SKIP_DEFAULT
     $                 (93, TRIM(GLOBALDIR)//'/'//'fort.93',
     $                 IceDescript, TimeLoc, it, storeOne,0.d0)
                  ENDIF
#else
                  CALL WRITE_GBL_FILE_SKIP_DEFAULT
     $               (93, TRIM(GLOBALDIR)//'/'//'fort.93',
     $               IceDescript, TimeLoc, it, storeOne,0.d0)
#endif
               ENDIF
!  tcm v50.75 removed ifdef cswan to allow for use whenever nrs=3 or nrs=4
               IF ((ABS(NRS).EQ.3).OR.(ABS(NRS).EQ.4)) THEN
!#ifdef CSWAN
Casey 090302: Added these lines for output of radiation stresses.
#ifdef CMPI
                  IF(MNWPROC.GT.0) THEN
                    CALL WRITE_GBL_FILE_THROUGH_WRITER
     &                   (164,TRIM(GLOBALDIR)//'/'//'rads.64',
     &                   RSDescript,TimeLoc,it,
     &                   writer_storeTwo, 0.D0)
                  ELSE
                     CALL WRITE_GBL_FILE_SKIP_DEFAULT
     &                 (164, TRIM(GLOBALDIR)//'/'//'rads.64',
     &                 RSDescript, TimeLoc, it, storeTwo, 0.D0)
                  ENDIF
#else
                  CALL WRITE_GBL_FILE_SKIP_DEFAULT
     &              (164, TRIM(GLOBALDIR)//'/'//'rads.64',
     &              RSDescript, TimeLoc, it, storeTwo, 0.D0)
#endif
!#endif           
               ENDIF
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
Casey 090302: Added file name.
            CALL writeOutArray(83, TimeLoc, IT, ConcDescript,
     &                         pack83, unpackOne, IGCP,
     &                         'fort.83        ')
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
C...  they have.
      CALL FLUSH_WRITERS()
#endif


Csb 11/11/2006
C...   OUTPUT MAXIMUM/MINIMUM GLOBAL DATA AT THE LAST TIME STEP.
C...
C     jgf48.4636 Updated to newer data structure, so that we can
C     compare with values from hot start, so that mins and maxes
C     can be preserved across hot starts.
      IF(IT.EQ.NT) THEN
         CALL writeOutArrayMinMax(311, TimeLoc, IT, EtaMaxDescript,
     &        packOne, unpackOne)
         CALL writeOutArrayMinMax(312, TimeLoc, IT, UMaxDescript,
     &        packOne, unpackOne)
         IF (NWS.NE.0) THEN
            CALL writeOutArrayMinMax(313, TimeLoc, IT, PrMinDescript,
     &           packOne, unpackOne)
            CALL writeOutArrayMinMax(314, TimeLoc, IT, WVMaxDescript,
     &           packOne, unpackOne)
            IF (NRS.NE.0) THEN
               CALL writeOutArrayMinMax(315, TimeLoc, IT, RSMaxDescript,
     &              packOne, unpackOne)
            ENDIF
         ENDIF
C       OUTPUT THE NODES WHERE THE ELEMENTAL SLOPE LIMITER HAS BEEN
C       TURNED ON DURING THE RUN
         IF(LoadEleSlopeLim)THEN
            CALL writeOutArrayMinMax(323, TimeLoc, IT, ESLDescript,
     &               packOne, unpackOne)
         ENDIF
      ENDIF

      RETURN
 2120 FORMAT(2X,1pE20.10E3,5X,I10)
 2453 FORMAT(2x, i8, 2x, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3)
 2454 FORMAT(2X,I8,2(2X,1pE20.10E3))
C-----------------------------------------------------------------------
      END SUBROUTINE WRITE_OUTPUT_2D
C-----------------------------------------------------------------------

C----------------------------------------------------------------------
C       S U B R O U T I N E
C          W R I T E  H A R M O N I C   A N A L Y S I S   O U T P U T
C----------------------------------------------------------------------
C     jgf49.44: Subroutine to write harmonic analysis output files
C     in serial or in parallel using globalio.
C----------------------------------------------------------------------
      SUBROUTINE writeHarmonicAnalysisOutput(ITIME)
      USE SIZES, ONLY : SZ, MNHARF, LOCALDIR, GLOBALDIR,
     &                  WRITE_LOCAL_FILES, MNPROC, MYPROC
C RJW added MNPROC and MYPROC
C removed from global  myProc, MNPROC,
      USE GLOBAL, ONLY : OutputDataDescript_t, NSTAE_G,
     &                   NOUTE, NSTAE, IMAP_STAE_LG, ICS, SLEL,
     &                   XEL, SFEL, YEL, NSTAV_G, NSTAV,  MNP,
     &                   IMAP_STAV_LG, SLEV, XEV, SFEV,
     &                   YEV, NP_G, NP, NODES_LG, screenMessage, DEBUG,
     &                   setMessageSource, unsetMessageSource
      USE ADCIRC_MOD, ONLY : ADCIRC_Terminate
      USE GLOBAL_IO, ONLY:
     &    open_gbl_file, write_gbl_file, collectFullDomainArray,
     &    packOne, unpackOne, packTwo, unpackTwo, packMbyNP, unpackMbyNP
      USE HARM, ONLY :
     &    emag, phasede, umag, vmag, phasedu, phasedv, emagt, phaseden,
     &    umagt, vmagt, phasedut, phasedvt, elav_g, elva_g,
     &    xvelav_g, xvelva_g, yvelav_g, yvelva_g,
     &    elav, elva, xvelav, uav, uavdif, nhase, nhasv, nhage, nhagv,
     &    charmv, fmv, eav, esq, eavdif, evadif, iharind, ITHAS,
     &    xvelva, usq, uvadif, yvelav, vav, vavdif, yvelva, vsq, vvadif

#ifdef CMPI
      USE MESSENGER, ONLY : msg_fini, msg_lbcast
#endif

      IMPLICIT NONE
      INTEGER I, J, N                      !local loop counters
      INTEGER, intent(in) :: ITIME  !tcm v49.64.02 added
C
      type(OutputDataDescript_t), SAVE :: HAElevStaMagDescript
      type(OutputDataDescript_t), SAVE :: HAElevStaPhaseDescript
C
      type(OutputDataDescript_t), SAVE :: HAVelStaUMagDescript
      type(OutputDataDescript_t), SAVE :: HAVelStaVMagDescript
      type(OutputDataDescript_t), SAVE :: HAVelStaUPhaseDescript
      type(OutputDataDescript_t), SAVE :: HAVelStaVPhaseDescript
C
      type(OutputDataDescript_t), SAVE :: HAElevMagDescript
      type(OutputDataDescript_t), SAVE :: HAElevPhaseDescript
C
      type(OutputDataDescript_t), SAVE :: HAVelUMagDescript
      type(OutputDataDescript_t), SAVE :: HAVelVMagDescript
      type(OutputDataDescript_t), SAVE :: HAVelUPhaseDescript
      type(OutputDataDescript_t), SAVE :: HAVelVPhaseDescript
C
      type(OutputDataDescript_t), SAVE :: MVDescript
C
      type(OutputDataDescript_t), SAVE :: dummyDescript
C
C     Fulldomain Stations in parallel:
      REAL(SZ), ALLOCATABLE, TARGET :: EMAG_g(:,:)    ! elevation magnitudes
      REAL(SZ), ALLOCATABLE, TARGET :: PHASEDE_g(:,:) ! elevation phases
      REAL(SZ), ALLOCATABLE, TARGET :: UMAG_g(:,:)    ! u velocity magnitudes
      REAL(SZ), ALLOCATABLE, TARGET :: VMAG_g(:,:)    ! v velocity magnitudes
      REAL(SZ), ALLOCATABLE, TARGET :: PHASEDU_g(:,:) ! u velocity phases
      REAL(SZ), ALLOCATABLE, TARGET :: PHASEDV_g(:,:) ! v velocity phases
C     Fulldomain Nodes in parallel:
      REAL(SZ), ALLOCATABLE, TARGET :: EMAGT_g(:,:)    ! elevation magnitudes
      REAL(SZ), ALLOCATABLE, TARGET :: PHASEDEN_g(:,:) ! elevation phases
      REAL(SZ), ALLOCATABLE, TARGET :: UMAGT_g(:,:)    ! u velocity magnitudes
      REAL(SZ), ALLOCATABLE, TARGET :: VMAGT_g(:,:)    ! v velocity magnitudes
      REAL(SZ), ALLOCATABLE, TARGET :: PHASEDUT_g(:,:) ! u velocity phases
      REAL(SZ), ALLOCATABLE, TARGET :: PHASEDVT_g(:,:) ! v velocity phases
C     Fulldomain Time means and variance calculation in parallel:
      REAL(SZ), ALLOCATABLE, TARGET :: EAV_g(:)
      REAL(SZ), ALLOCATABLE, TARGET :: ESQ_g(:)
      REAL(SZ), ALLOCATABLE, TARGET :: EAVDIF_g(:)
      REAL(SZ), ALLOCATABLE, TARGET :: EVADIF_g(:)
C
      REAL(SZ), ALLOCATABLE, TARGET :: UAV_g(:)
      REAL(SZ), ALLOCATABLE, TARGET :: USQ_g(:)
      REAL(SZ), ALLOCATABLE, TARGET :: UAVDIF_g(:)
      REAL(SZ), ALLOCATABLE, TARGET :: UVADIF_g(:)
C
      REAL(SZ), ALLOCATABLE, TARGET :: VAV_g(:)
      REAL(SZ), ALLOCATABLE, TARGET :: VSQ_g(:)
      REAL(SZ), ALLOCATABLE, TARGET :: VAVDIF_g(:)
      REAL(SZ), ALLOCATABLE, TARGET :: VVADIF_g(:)
C
      LOGICAL, SAVE :: FirstCall = .true.
C
      call setMessageSource("writeHarmonicAnalysisOutput")
#ifdef WRITE_OUTPUT_TRACE
      call screen(DEBUG,"Enter.")
#endif
      ! JUST RETURN if harmonic analysis was not specified
      ! tcm v49.64.02 -added itime le ithas to be
      ! consistent with what is in solveHarmonicAnalysis
      IF ((IHARIND.eq.0).or.(ITIME.LE.ITHAS)) THEN
#ifdef WRITE_OUTPUT_TRACE
      call screen(DEBUG,"Return.")
#endif
         RETURN
      ENDIF
C
C     This subroutine should only be called once, at the end of a run, so
C     a check of FirstCall is not necessary. However, we'll still do it, in
C     case the subroutine is ever called more than once in some future
C     scenario.
      IF (FirstCall) THEN
         IF (NHASE.ne.0) THEN
            ! fort.51 magnitude
            IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
              ALLOCATE(EMAG_g(MNHARF,NSTAE_G))
            ENDIF
            HAElevStaMagDescript % specifier            = NHASE
            HAElevStaMagDescript % initial_value        = 0.0
            HAElevStaMagDescript % num_items_per_record = MNHARF
            HAElevStaMagDescript % num_fd_records       = NSTAE_G
            HAElevStaMagDescript % num_records_this     = NSTAE
            HAElevStaMagDescript % imap                 => IMAP_STAE_LG
            HAElevStaMagDescript % array2D              => EMAG
            HAElevStaMagDescript % array2D_g            => EMAG_g
            HAElevStaMagDescript % ConsiderWetDry       = .FALSE.
            HAElevStaMagDescript % field_name           = 'HAElevStaMag'
            IF (ICS.eq.2) THEN
               HAElevStaMagDescript % x_coord           => SLEL
               HAElevStaMagDescript % y_coord           => SFEL
            ELSE
               HAElevStaMagDescript % x_coord           => XEL
               HAElevStaMagDescript % y_coord           => YEL
            ENDIF
            ! fort.51 phase
            IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
               ALLOCATE(PHASEDE_g(MNHARF,NSTAE_G))
            ENDIF
            HAElevStaPhaseDescript % specifier            = NHASE
            HAElevStaPhaseDescript % initial_value        = 0.0
            HAElevStaPhaseDescript % num_items_per_record = MNHARF
            HAElevStaPhaseDescript % num_fd_records       = NSTAE_G
            HAElevStaPhaseDescript % num_records_this     = NSTAE
            HAElevStaPhaseDescript % imap                => IMAP_STAE_LG
            HAElevStaPhaseDescript % array2D              => PHASEDE
            HAElevStaPhaseDescript % array2D_g            => PHASEDE_g
            HAElevStaPhaseDescript % ConsiderWetDry       = .FALSE.
            HAElevStaPhaseDescript % field_name     = 'HAElevStaPhase'
            IF (ICS.eq.2) THEN
               HAElevStaPhaseDescript % x_coord           => SLEL
               HAElevStaPhaseDescript % y_coord           => SFEL
            ELSE
               HAElevStaPhaseDescript % x_coord           => XEL
               HAElevStaPhaseDescript % y_coord           => YEL
            ENDIF
         ENDIF
C
C
         IF (NHASV.ne.0) THEN
            ! fort.52 u velocity magnitude
            IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
               ALLOCATE(UMAG_g(MNHARF,NSTAV_G))
            ENDIF
            HAVelStaUMagDescript % specifier            = NHASV
            HAVelStaUMagDescript % initial_value        = 0.0
            HAVelStaUMagDescript % num_items_per_record = MNHARF
            HAVelStaUMagDescript % num_fd_records       = NSTAV_G
            HAVelStaUMagDescript % num_records_this     = NSTAV
            HAVelStaUMagDescript % imap                 => IMAP_STAV_LG
            HAVelStaUMagDescript % array2D              => UMAG
            HAVelStaUMagDescript % array2D_g            => UMAG_g
            HAVelStaUMagDescript % ConsiderWetDry       = .FALSE.
            HAVelStaUMagDescript % alternate_value      = 0.0
            HAVelStaUMagDescript % field_name           = 'HAVelStaUMag'
            IF (ICS.eq.2) THEN
               HAVelStaUMagDescript % x_coord           => SLEV
               HAVelStaUMagDescript % y_coord           => SFEV
            ELSE
               HAVelStaUMagDescript % x_coord           => XEV
               HAVelStaUMagDescript % y_coord           => YEV
            ENDIF
            ! fort.52 v velocity magnitude
            IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
               ALLOCATE(VMAG_g(MNHARF,NSTAV_G))
            ENDIF
            HAVelStaVMagDescript % specifier            = NHASV
            HAVelStaVMagDescript % initial_value        = 0.0
            HAVelStaVMagDescript % num_items_per_record = MNHARF
            HAVelStaVMagDescript % num_fd_records       = NSTAV_G
            HAVelStaVMagDescript % num_records_this     = NSTAV
            HAVelStaVMagDescript % imap                 => IMAP_STAV_LG
            HAVelStaVMagDescript % array2D              => VMAG
            HAVelStaVMagDescript % array2D_g            => VMAG_g
            HAVelStaVMagDescript % ConsiderWetDry       = .FALSE.
            HAVelStaVMagDescript % alternate_value      = 0.0
            HAVelStaVMagDescript % field_name           = 'HAVelStaVMag'
            IF (ICS.eq.2) THEN
               HAVelStaVMagDescript % x_coord           => SLEV
               HAVelStaVMagDescript % y_coord           => SFEV
            ELSE
               HAVelStaVMagDescript % x_coord           => XEV
               HAVelStaVMagDescript % y_coord           => YEV
            ENDIF
            ! fort.52 u velocity phase
            IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
               ALLOCATE(PHASEDU_g(MNHARF,NSTAV_G))
            ENDIF
            HAVelStaUPhaseDescript % specifier            = NHASV
            HAVelStaUPhaseDescript % initial_value        = 0.0
            HAVelStaUPhaseDescript % num_items_per_record = MNHARF
            HAVelStaUPhaseDescript % num_fd_records       = NSTAV_G
            HAVelStaUPhaseDescript % num_records_this     = NSTAV
            HAVelStaUPhaseDescript % imap         => IMAP_STAV_LG
            HAVelStaUPhaseDescript % array2D              => PHASEDU
            HAVelStaUPhaseDescript % array2D_g            => PHASEDU_g
            HAVelStaUPhaseDescript % ConsiderWetDry       = .FALSE.
            HAVelStaUPhaseDescript % alternate_value      = 0.0
            HAVelStaUPhaseDescript % field_name   = 'HAVelStaUPhase'
            IF (ICS.eq.2) THEN
               HAVelStaUPhaseDescript % x_coord           => SLEV
               HAVelStaUPhaseDescript % y_coord           => SFEV
            ELSE
               HAVelStaUPhaseDescript % x_coord           => XEV
               HAVelStaUPhaseDescript % y_coord           => YEV
            ENDIF
            ! fort.52 v velocity phase
            IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
               ALLOCATE(PHASEDV_g(MNHARF,NSTAV_G))
            ENDIF
            HAVelStaVPhaseDescript % specifier            = NHASV
            HAVelStaVPhaseDescript % initial_value        = 0.0
            HAVelStaVPhaseDescript % num_items_per_record = MNHARF
            HAVelStaVPhaseDescript % num_fd_records       = NSTAV_G
            HAVelStaVPhaseDescript % num_records_this     = NSTAV
            HAVelStaVPhaseDescript % imap            => IMAP_STAV_LG
            HAVelStaVPhaseDescript % array2D              => PHASEDV
            HAVelStaVPhaseDescript % array2D_g            => PHASEDV_g
            HAVelStaVPhaseDescript % ConsiderWetDry       = .FALSE.
            HAVelStaVPhaseDescript % alternate_value      = 0.0
            HAVelStaVPhaseDescript % field_name      = 'HAVelStaVPhase'
            IF (ICS.eq.2) THEN
               HAVelStaVPhaseDescript % x_coord           => SLEV
               HAVelStaVPhaseDescript % y_coord           => SFEV
            ELSE
               HAVelStaVPhaseDescript % x_coord           => XEV
               HAVelStaVPhaseDescript % y_coord           => YEV
            ENDIF

         ENDIF
C
C
         IF (NHAGE.ne.0) THEN
            ! fort.53 Full Domain Elevation Magnitude
            IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
               ALLOCATE(EMAGT_g(MNHARF,NP_G))
            ENDIF
            HAElevMagDescript % specifier            = NHAGE
            HAElevMagDescript % initial_value        = 0.0
            HAElevMagDescript % num_items_per_record = MNHARF
            HAElevMagDescript % num_fd_records       = NP_G
            HAElevMagDescript % num_records_this     = NP
            HAElevMagDescript % imap                 => NODES_LG
            HAElevMagDescript % array2D              => EMAGT
            HAElevMagDescript % array2D_g            => EMAGT_g
            HAElevMagDescript % ConsiderWetDry       = .FALSE.
            HAElevMagDescript % alternate_value      = -99999.0
            HAElevMagDescript % field_name           = 'HAElevMag'
            ! fort.53 Full Domain Elevation Phase
            IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
               ALLOCATE(PHASEDEN_g(MNHARF,NP_G))
            ENDIF
            HAElevPhaseDescript % specifier            = NHAGE
            HAElevPhaseDescript % initial_value        = 0.0
            HAElevPhaseDescript % num_items_per_record = MNHARF
            HAElevPhaseDescript % num_fd_records       = NP_G
            HAElevPhaseDescript % num_records_this     = NP
            HAElevPhaseDescript % imap                 => NODES_LG
            HAElevPhaseDescript % array2D              => PHASEDEN
            HAElevPhaseDescript % array2D_g            => PHASEDEN_g
            HAElevPhaseDescript % ConsiderWetDry       = .FALSE.
            HAElevPhaseDescript % alternate_value      = -99999.0
            HAElevPhaseDescript % field_name           = 'HAElevPhase'

         ENDIF
C
C
         IF (NHAGV.ne.0) THEN
            ! fort.54 Full Domain U Velocity Magnitude
            IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
               ALLOCATE(UMAGT_g(MNHARF,NP_G))
            ENDIF
            HAVelUMagDescript % specifier            = NHAGV
            HAVelUMagDescript % initial_value        = 0.0
            HAVelUMagDescript % num_items_per_record = MNHARF
            HAVelUMagDescript % num_fd_records       = NP_G
            HAVelUMagDescript % num_records_this     = NP
            HAVelUMagDescript % imap                 => NODES_LG
            HAVelUMagDescript % array2D              => UMAGT
            HAVelUMagDescript % array2D_g            => UMAGT_g
            HAVelUMagDescript % ConsiderWetDry       = .FALSE.
            HAVelUMagDescript % alternate_value      = -99999.0
            HAVelUMagDescript % field_name           = 'HAVelUMag'
            ! fort.54 Full Domain V Velocity Magnitude
            IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
               ALLOCATE(VMAGT_g(MNHARF,NP_G))
            ENDIF
            HAVelVMagDescript % specifier            = NHAGV
            HAVelVMagDescript % initial_value        = 0.0
            HAVelVMagDescript % num_items_per_record = MNHARF
            HAVelVMagDescript % num_fd_records       = NP_G
            HAVelVMagDescript % num_records_this     = NP
            HAVelVMagDescript % imap                 => NODES_LG
            HAVelVMagDescript % array2D              => VMAGT
            HAVelVMagDescript % array2D_g            => VMAGT_g
            HAVelVMagDescript % ConsiderWetDry       = .FALSE.
            HAVelVMagDescript % alternate_value      = -99999.0
            HAVelVMagDescript % field_name           = 'HAVelVMag'
            ! fort.54 Full Domain U Velocity Phase
            IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
               ALLOCATE(PHASEDUT_g(MNHARF,NP_G))
            ENDIF
            HAVelUPhaseDescript % specifier            = NHAGV
            HAVelUPhaseDescript % initial_value        = 0.0
            HAVelUPhaseDescript % num_items_per_record = MNHARF
            HAVelUPhaseDescript % num_fd_records       = NP_G
            HAVelUPhaseDescript % num_records_this     = NP
            HAVelUPhaseDescript % imap                 => NODES_LG
            HAVelUPhaseDescript % array2D              => PHASEDUT
            HAVelUPhaseDescript % array2D_g            => PHASEDUT_g
            HAVelUPhaseDescript % ConsiderWetDry       = .FALSE.
            HAVelUPhaseDescript % alternate_value      = -99999.0
            HAVelUPhaseDescript % field_name           = 'HAVelUPhase'
            ! fort.54 Full Domain V Velocity Phase
            IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
               ALLOCATE(PHASEDVT_g(MNHARF,NP_G))
            ENDIF
            HAVelVPhaseDescript % specifier            = NHAGV
            HAVelVPhaseDescript % initial_value        = 0.0
            HAVelVPhaseDescript % num_items_per_record = MNHARF
            HAVelVPhaseDescript % num_fd_records       = NP_G
            HAVelVPhaseDescript % num_records_this     = NP
            HAVelVPhaseDescript % imap                 => NODES_LG
            HAVelVPhaseDescript % array2D              => PHASEDVT
            HAVelVPhaseDescript % array2D_g            => PHASEDVT_g
            HAVelVPhaseDescript % ConsiderWetDry       = .FALSE.
            HAVelVPhaseDescript % alternate_value      = -99999.0
            HAVelVPhaseDescript % field_name           = 'HAVelVPhase'
         ENDIF
C
C        ! fort.55 means and variance calculations

C        ! fort.55 will use the same output format (ascii, netcdf, etc)
C        ! as the full domain elevation harmonic analysis file (fort.53)
         IF (CHARMV.eqv..true.) THEN
            IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
               ALLOCATE(eav_g(NP_g),eavdif_g(NP_G),
     &                  esq_g(NP_G),evadif_g(NP_G))
               ALLOCATE(uav_g(NP_g),uavdif_g(NP_G),
     &                  usq_g(NP_G),uvadif_g(NP_G))
               ALLOCATE(vav_g(NP_g),vavdif_g(NP_G),
     &                  vsq_g(NP_G),vvadif_g(NP_G))
            ENDIF
            ! MVDescript is used over and over for different arrays
            MVDescript % specifier            = NHAGE
            MVDescript % initial_value        = 0.0
            MVDescript % num_items_per_record = 1
            MVDescript % num_fd_records       = NP_G
            MVDescript % num_records_this     = NP
            MVDescript % imap                 => NODES_LG
            MVDescript % array                => ELAV
            MVDescript % array_g              => ELAV_g
            MVDescript % ConsiderWetDry       = .false.
            MVDescript % alternate_value      = -99999.0
            MVDescript % field_name           = 'HAMeansVariances'
         ENDIF
      ENDIF  ! firstcall
C
C     Collect fulldomain data if running in parallel.
#ifdef CMPI
      IF (NHASE.ne.0) THEN
         CALL collectFullDomainArray(HAElevStaMagDescript,
     &             packMbyNP, unpackMbyNP)
         CALL collectFullDomainArray(HAElevStaPhaseDescript,
     &             packMbyNP, unpackMbyNP)
      ENDIF
      IF (NHASV.ne.0) THEN
         CALL collectFullDomainArray(HAVelStaUMagDescript,
     &             packMbyNP, unpackMbyNP)
         CALL collectFullDomainArray(HAVelStaVMagDescript,
     &             packMbyNP, unpackMbyNP)
         CALL collectFullDomainArray(HAVelStaUPhaseDescript,
     &             packMbyNP, unpackMbyNP)
         CALL collectFullDomainArray(HAVelStaVPhaseDescript,
     &             packMbyNP, unpackMbyNP)
      ENDIF
      IF (NHAGE.ne.0) THEN
         CALL collectFullDomainArray(HAElevMagDescript,
     &             packMbyNP, unpackMbyNP)
         CALL collectFullDomainArray(HAElevPhaseDescript,
     &             packMbyNP, unpackMbyNP)
      ENDIF
      IF (NHAGV.ne.0) THEN
         CALL collectFullDomainArray(HAVelUMagDescript,
     &             packMbyNP, unpackMbyNP)
         CALL collectFullDomainArray(HAVelVMagDescript,
     &             packMbyNP, unpackMbyNP)
         CALL collectFullDomainArray(HAVelUPhaseDescript,
     &             packMbyNP, unpackMbyNP)
         CALL collectFullDomainArray(HAVelVPhaseDescript,
     &             packMbyNP, unpackMbyNP)
      ENDIF
      IF (CHARMV.eqv..true.) THEN
         CALL collectFullDomainArray(MVDescript, packOne, unpackOne)
         MVDescript%array => eav
         MVDescript%array_g => eav_g
         CALL collectFullDomainArray(MVDescript, packOne, unpackOne)
         MVDescript%array => eavdif
         MVDescript%array_g => eavdif_g
         CALL collectFullDomainArray(MVDescript, packOne, unpackOne)
         MVDescript%array => elva
         MVDescript%array_g => elva_g
         CALL collectFullDomainArray(MVDescript, packOne, unpackOne)
         MVDescript%array => esq
         MVDescript%array_g => esq_g
         CALL collectFullDomainArray(MVDescript, packOne, unpackOne)
         MVDescript%array => evadif
         MVDescript%array_g => evadif_g
         CALL collectFullDomainArray(MVDescript, packOne, unpackOne)
         MVDescript%array => xvelav
         MVDescript%array_g => xvelav_g
         CALL collectFullDomainArray(MVDescript, packOne, unpackOne)
         MVDescript%array => uav
         MVDescript%array_g => uav_g
         CALL collectFullDomainArray(MVDescript, packOne, unpackOne)
         MVDescript%array => uavdif
         MVDescript%array_g => uavdif_g
         CALL collectFullDomainArray(MVDescript, packOne, unpackOne)
         MVDescript%array => xvelva
         MVDescript%array_g => xvelva_g
         CALL collectFullDomainArray(MVDescript, packOne, unpackOne)
         MVDescript%array => usq
         MVDescript%array_g => usq_g
         CALL collectFullDomainArray(MVDescript, packOne, unpackOne)
         MVDescript%array => uvadif
         MVDescript%array_g => uvadif_g
         CALL collectFullDomainArray(MVDescript, packOne, unpackOne)
         MVDescript%array => yvelav
         MVDescript%array_g => yvelav_g
         CALL collectFullDomainArray(MVDescript, packOne, unpackOne)
         MVDescript%array => vav
         MVDescript%array_g => vav_g
         CALL collectFullDomainArray(MVDescript, packOne, unpackOne)
         MVDescript%array => vavdif
         MVDescript%array_g => vavdif_g
         CALL collectFullDomainArray(MVDescript, packOne, unpackOne)
         MVDescript%array => yvelva
         MVDescript%array_g => yvelva_g
         CALL collectFullDomainArray(MVDescript, packOne, unpackOne)
         MVDescript%array => vsq
         MVDescript%array_g => vsq_g
         CALL collectFullDomainArray(MVDescript, packOne, unpackOne)
         MVDescript%array => vvadif
         MVDescript%array_g => vvadif_g
         CALL collectFullDomainArray(MVDescript, packOne, unpackOne)
      ENDIF
#endif
c
c     Now write data to file; format according to the file type specifier.
C
c     Harmonic analysis at elevation stations (fort.51).
      call writeOutHarmonicArrays(51, 1, HAElevStaMagDescript,
     &        HAElevStaPhaseDescript, dummyDescript, dummyDescript)
C     Harmonic analysis at velocity stations (fort.52).
      call writeOutHarmonicArrays(52, 2, HAVelStaUMagDescript,
     &        HAVelStaUPhaseDescript, HAVelStaVMagDescript,
     &        HAVelStaVPhaseDescript)
C     Harmonic analysis at all nodes for elevation (fort.53).
      call writeOutHarmonicArrays(53, 1, HAElevMagDescript,
     &        HAElevPhaseDescript, dummyDescript, dummyDescript)
C     Harmonic analysis at all nodes for velocity (fort.54).
      call writeOutHarmonicArrays(54, 2,HAVelUMagDescript,
     &        HAVelUPhaseDescript, HAVelVMagDescript,
     &        HAVelVPhaseDescript)
C
C     Means and variance calculations, if specified (fort.55).
      if ((CHARMV.eqv..true.).and.(FMV.gt.1.0d-3)) then
         IF ((MNPROC.gt.1).and.(MyProc.eq.0)
     &        .and.(WRITE_LOCAL_FILES.eqv..false.)) THEN
            OPEN(55,FILE=TRIM(GLOBALDIR)//'/'//'fort.55',
     &         STATUS='REPLACE',ACCESS='SEQUENTIAL',ACTION='WRITE')
            WRITE(55,*) NP_G
            DO n=1,NP_G
               write(55,*) n
               write(55,7637) elav_g(n),eav_g(n),eavdif_g(n),
     &                        elva_g(n),esq_g(n),evadif_g(n)
            END DO
            DO n=1,NP_G
               write(55,*) n
               write(55,7637) xvelav_g(n),uav_g(n),uavdif_g(n),
     &                        xvelva_g(n),usq_g(n),uvadif_g(n)
               write(55,7637) yvelav_g(n),vav_g(n),vavdif_g(n),
     &                        yvelva_g(n),vsq_g(n),vvadif_g(n)
            END DO
            CLOSE(55)
         ENDIF
         IF ((MNPROC.eq.1).or.(WRITE_LOCAL_FILES.eqv..true.)) THEN
            OPEN(55,FILE=TRIM(LOCALDIR)//'/'//'fort.55',
     &         STATUS='REPLACE',ACCESS='SEQUENTIAL',ACTION='WRITE')
            WRITE(55,*) NP
            DO n=1,NP
               write(55,*) n
               write(55,7637) elav(n),eav(n),eavdif(n),
     &                        elva(n),esq(n),evadif(n)
            END DO
            DO n=1,NP
               write(55,*) n
               write(55,7637) xvelav(n),uav(n),uavdif(n),
     &                        xvelva(n),usq(n),uvadif(n)
               write(55,7637) yvelav(n),vav(n),vavdif(n),
     &                        yvelva(n),vsq(n),vvadif(n)
            END DO
            CLOSE(55)
         ENDIF
      ENDIF

 7637 format(2x,3(1pE16.8E3,1x),2x,3(1pE16.8E3,1x))
#ifdef WRITE_OUTPUT_TRACE
      call screen(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE writeHarmonicAnalysisOutput
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E
C            W R I T E   O U T   H A R M O N I C  A R R A Y S
C-----------------------------------------------------------------------
C     jgf49.44 Writes the harmonic arrays in output files of the desired
C     format. Only ascii is supported at this point.
C-----------------------------------------------------------------------
      SUBROUTINE writeOutHarmonicArrays(lun, sets, mag1, ph1, mag2, ph2)
      USE SIZES, ONLY : WRITE_LOCAL_FILES, MNPROC, LOCALDIR, GLOBALDIR,
     & MYPROC
      USE HARM, ONLY : nf, nfreq, hafreq, haff, haface, namefr
      USE GLOBAL, ONLY : OutputDataDescript_t, screenMessage,
     &                   DEBUG, setMessageSource, unsetMessageSource
c     & ,myproc
      IMPLICIT NONE
      INTEGER, intent(in) :: lun  ! adcirc logical unit number for i/o
      INTEGER, intent(in) :: sets ! how many sets of mags/phases to write
      type(OutputDataDescript_t), intent(in) :: mag1
      type(OutputDataDescript_t), intent(in) :: ph1
      type(OutputDataDescript_t), intent(in) :: mag2
      type(OutputDataDescript_t), intent(in) :: ph2
C
      INTEGER :: i,j,n           ! loop counters
C     local vars
      CHARACTER(7) :: fn     ! outfile name (valid for lun between 10 and 99)
C
      call setMessageSource("writeOutHarmonicArrays")
#ifdef WRITE_OUTPUT_TRACE
      call screen(DEBUG,"Enter.")
#endif
C
C     initialize output file name
      fn(1:5) = 'fort.'
      WRITE(fn(6:7),2) lun

C     write data according to format specifier from fort.15 (e.g., NHASE)
      SELECT CASE (mag1%specifier)
      CASE(0) ! no output; do nothing
      CASE(1) ! ascii text
         IF ( (MNPROC.gt.1).and.(MyProc.eq.0)
     &        .and.(WRITE_LOCAL_FILES.eqv..false.)) THEN
            OPEN(lun,FILE=TRIM(GLOBALDIR)//'/'//fn, STATUS='REPLACE',
     &           ACCESS='SEQUENTIAL',ACTION='WRITE')
            write(lun,*) nfreq+nf
            do j=1,nfreq+nf
               write(lun,3679) hafreq(j),HAFF(j),HAFACE(j),namefr(j)
            end do
            write(lun,*) mag1%num_fd_records
            DO N=1,mag1%num_fd_records
               write(lun,*) N
               do i=1,nfreq+nf
                  if (sets.eq.2) then
                     write(lun,6636) mag1%array2D_g(i,n),
     &                               ph1%array2D_g(i,n),
     &                               mag2%array2D_g(i,n),
     &                               ph2%array2D_g(i,n)
                  else
                     write(lun,6635) mag1%array2D_g(i,n),
     &                               ph1%array2D_g(i,n)
                  endif
               enddo
            ENDDO
            CLOSE(lun)
         ENDIF
C
         IF ((MNPROC.eq.1).or.(WRITE_LOCAL_FILES.eqv..true.)) THEN
            OPEN(lun,FILE=TRIM(LOCALDIR)//'/'//fn,STATUS='REPLACE',
     &           ACCESS='SEQUENTIAL',ACTION='WRITE')
            write(lun,*) nfreq+nf
            do j=1,nfreq+nf
               write(lun,3679) hafreq(j),HAFF(j),HAFACE(j),namefr(j)
            end do
            write(lun,*) mag1%num_records_this
            DO N=1,mag1%num_records_this
               write(lun,*) N
               do i=1,nfreq+nf
                  if (sets.eq.2) then
                     write(lun,6636) mag1%array2D(i,n),
     &                               ph1%array2D(i,n),
     &                               mag2%array2D(i,n),
     &                               ph2%array2D(i,n)
                  else
                     write(lun,6635) mag1%array2D(i,n),
     &                               ph1%array2D(i,n)
                  endif
               enddo
            ENDDO
            CLOSE(lun)
         ENDIF
      CASE(3,5)
         write(*,*)
     & "ERROR: Harmonic analysis output is not yet available in NetCDF."
      CASE DEFAULT
        write(*,*)
     &     "ERROR: Invalid harmonic analysis format specifier: ",
     &     mag1%specifier
      END SELECT

 2    FORMAT(I2)
 6635 format(2x,1pE16.8E3,1x,0pf11.4)
 6636 format(2x,1pE16.8E3,1x,0pf11.4,2x,1pE16.8E3,1x,0pf11.4)
 3679 format(1x,e20.10,1x,f10.7,1x,f12.8,1x,a10)
#ifdef WRITE_OUTPUT_TRACE
      call screen(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE writeOutHarmonicArrays
C-----------------------------------------------------------------------

C----------------------------------------------------------------------
C       S U B R O U T I N E   I N I T   O U T P U T   3 D
C----------------------------------------------------------------------
C     jgf49.43.19:  Subroutine to initialize ascii or nonportable binary
C     output for a 3D run.
C----------------------------------------------------------------------
      SUBROUTINE initOutput3D(lun, specifier, nset, nPerSet, nPerSet_g,
     &                        tsPeriod, outputShape, recNum)
      USE SIZES, ONLY : SZ, LOCALDIR, GLOBALDIR, MNPROC,
     &   WRITE_LOCAL_FILES, NBYTE, Myproc
C RJW added myproc
C MYPROC,
      USE GLOBAL, ONLY : RUNDES, RUNID, AGRID, RDES4, RID4,
     &   AID4, RDES8, RID8, AID8, DTDP, IHOT, DEBUG, scratchMessage,
     &   logMessage, ERROR, allMessage
      USE GLOBAL_3DVS, ONLY : NFEN
      IMPLICIT NONE
C
      INTEGER, intent(in) :: lun       ! ADCIRC logical unit number
      INTEGER, intent(in) :: specifier ! format and append behavior
      INTEGER, intent(in) :: nset      ! number of output data sets in the file
      INTEGER, intent(in) :: nPerSet   ! num of stations/nodes in each data set
      INTEGER, intent(in) :: nPerSet_g ! full domain num of stations/nodes
      INTEGER, intent(in) :: tsPeriod  ! sim time steps between output writes
      INTEGER, intent(in) :: outputShape ! number of data columns in output file
      INTEGER, intent(inout) :: recNum ! number of lines written to output file
C
      CHARACTER(len=7) :: fn ! file name of 3D output file
      INTEGER :: N ! loop counter
      INTEGER :: ios ! i/o status, 0 indicates success
C
      ios=0
C
C     form file name from lun
      fn(1:5) = 'fort.'
      WRITE(fn(6:7),'(I2)') lun
C
      IF ((specifier.EQ.-1).OR.((specifier.eq.1).AND.(IHOT.EQ.0))) THEN
         !start a new ASCII file
         IF ((MNPROC.eq.1).or.(WRITE_LOCAL_FILES.eqv..true.)) THEN
            WRITE(16,497) lun
            OPEN(lun,FILE=TRIM(LOCALDIR)//'/'//fn,STATUS="REPLACE")
            WRITE(lun,499) RUNDes,RunID,AGrid
            WRITE(lun,498)
     &         nset,nPerSet,(tsPeriod*DTDP),tsPeriod,NFEN,outputShape

            CLOSE(lun)
            recNum=2
         ENDIF
         IF ((MNPROC.GT.1).AND.(WRITE_LOCAL_FILES.eqv..false.).AND.
     &        (myProc.eq.0)) THEN
            WRITE(16,497) lun
            OPEN(lun,FILE=TRIM(GLOBALDIR)//'/'//fn,STATUS="REPLACE")
            WRITE(lun,499) RUNDes,RunID,AGrid
            WRITE(lun,498) nset,nPerSet_g,tsPeriod*DTDP,tsPeriod,
     &                     NFEN,outputShape
            CLOSE(lun)
            recNum=2
         ENDIF
      ENDIF
      IF ((specifier.eq.1).and.(IHOT.NE.0)) THEN
         ! find spot in existing ASCII file
 123     FORMAT('Skipping down ',I12,' lines in the output file.')
         WRITE(scratchMessage,123) recNum
         CALL logMessage(DEBUG,scratchMessage)
         IF ((MNPROC.eq.1).or.(WRITE_LOCAL_FILES.eqv..true.)) THEN
            WRITE(16,497) lun
            OPEN(lun,FILE=TRIM(LOCALDIR)//'/'//fn)
            DO N=1,recNum
               READ(lun,*,END=456,ERR=456)
            ENDDO
            ! jgf49.59 Tell user there were not enough records.
456         IF (N.lt.recNum) THEN
               CALL allMessage(ERROR,
     &            "There are not enough records in the output file.")
               CALL allMessage(ERROR,
     &            "There will be a gap in the output data.")
            ENDIF
            ENDFILE(lun,IOSTAT=ios,ERR=789)
            ! jgf49.59 The gfortran compiler complains when we write an
            ! end-of-file record at the end of the file. Sheesh.
789         IF (ios.lt.0) THEN
               CALL logMessage(DEBUG,"End of output file reached.")! no big deal
            ELSE IF (ios.gt.0) THEN
               CALL allMessage(ERROR,"I/O error ending output file.")
            ENDIF
            CLOSE(lun)
         ENDIF
         IF ((MNPROC.GT.1).AND.(WRITE_LOCAL_FILES.eqv..false.).AND.
     &        (myProc.eq.0)) THEN
            WRITE(16,497) lun
            OPEN(lun,FILE=TRIM(GLOBALDIR)//'/'//fn)
            DO N=1,recNum
               READ(lun,*,END=654,ERR=654)
            ENDDO
            ! jgf49.59 Tell user there were not enough records.
654         IF (N.lt.recNum) THEN
               CALL allMessage(ERROR,
     &            "There are not enough records in the output file.")
               CALL allMessage(ERROR,
     &            "There will be a gap in the output data.")
            ENDIF
            ENDFILE(lun,IOSTAT=ios,ERR=987)
            ! jgf49.59 The gfortran compiler complains when we write an
            ! end-of-file record at the end of the file. Sheesh.
987         IF (ios.lt.0) THEN
               CALL logMessage(DEBUG,"End of output file reached.")! no big deal
            ELSE IF (ios.gt.0) THEN
               CALL allMessage(ERROR,"I/O error ending output file.")
            ENDIF
            CLOSE(lun)
         ENDIF
      ENDIF
      IF ((specifier.EQ.-2).OR.((specifier.eq.2).AND.(IHOT.EQ.0))) THEN
        !start a new BINARY file
        WRITE(16,496) lun
        OPEN(lun,FILE=TRIM(LOCALDIR)//'/'//fn,ACCESS='DIRECT',
     &       RECL=NByte,STATUS="REPLACE")
        recNum=0
        IF(NByte.EQ.4) THEN
          DO N=1,8
            recNum=recNum+1
            WRITE(lun,REC=recNum) RDES4(N)
          ENDDO
          DO N=1,6
            recNum=recNum+1
            WRITE(lun,REC=recNum) RID4(N)
          ENDDO
          DO N=1,6
            recNum=recNum+1
            WRITE(lun,REC=recNum) AID4(N)
          ENDDO
        ENDIF
        IF(NByte.EQ.8) THEN
          DO N=1,4
            recNum=recNum+1
            WRITE(lun,REC=recNum) RDES8(N)
          ENDDO
          DO N=1,3
            recNum=recNum+1
            WRITE(lun,REC=recNum) RID8(N)
          ENDDO
          DO N=1,3
            recNum=recNum+1
            WRITE(lun,REC=recNum) AID8(N)
          ENDDO
        ENDIF
        WRITE(lun,REC=recNum+1) nset
        WRITE(lun,REC=recNum+2) nPerSet
        WRITE(lun,REC=recNum+3) tsPeriod*DTDP
        WRITE(lun,REC=recNum+4) tsPeriod
        WRITE(lun,REC=recNum+5) NFEN
        WRITE(lun,REC=recNum+6) outputShape
        recNum=recNum+6
        CLOSE(lun)
      ENDIF
      IF ((specifier.eq.2).and.(IHOT.NE.0)) THEN
         !find spot in existing BINARY file
         WRITE(16,496) lun
         OPEN(lun,FILE=TRIM(LOCALDIR)//'/'//fn,ACCESS='DIRECT',
     &       RECL=NByte)
         CLOSE(lun)
      ENDIF
C

 499  FORMAT(1X,A32,2X,A24,2X,A24)
C..RJW bug fix in 498 (kendra found this)
 498  FORMAT(1X,I10,1X,I10,1X,1pE15.7E3,1X,I10,1X,I10,1X,I3)
 497  FORMAT(5X,'UNIT ',I2,' FORMAT WILL BE ASCII')
 496  FORMAT(5X,'UNIT ',I2,' FORMAT WILL BE BINARY')
C----------------------------------------------------------------------
      END SUBROUTINE initOutput3D
C----------------------------------------------------------------------


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
C     jgf49.43 Rewrote for globalio.
C-----------------------------------------------------------------------
      SUBROUTINE writeOutput3D(TimeLoc,IT)
      USE SIZES, ONLY : SZ, INPUTDIR, NBYTE, MNWPROC,
     &    WRITE_LOCAL_FILES, MNPROC, GLOBALDIR, MYPROC
C RJW added MYPROC
      USE GLOBAL, ONLY : ScreenUnit, BCFLAG_TEMP, scratchMessage,
     &    DEBUG, INFO, WARNING, ERROR, screenMessage, setMessageSource,
     &    unsetMessageSource, allMessage, OutputDataDescript_t,
     &    NP_G, NODES_LG
      USE GLOBAL_3DVS
      USE GLOBAL_IO, ONLY: packNPbyM, unpackNPbyM,
     &    collectFullDomainArray, packOne, unpackOne
#ifdef ADCNETCDF
       USE NodalAttributes, ONLY :
     &     nolibf, nwp, tau0, cf, eslm
       USE NetCDFIO, ONLY : initNetCDFOutputFile,
     &     writeOutArrayNetCDF
#endif
#ifdef CMPI
       USE MESSENGER, ONLY : msg_lbcast, msg_fini
#endif
      IMPLICIT NONE
      INTEGER, intent(in) :: IT
      REAL(8), intent(in) :: TimeLoc

      INTEGER NC1, NC2, NC3, NCEle
      INTEGER NM1, NM2, NM3
      INTEGER :: NN, NH, NEle, N1, N2, N3, k, i
      REAL(SZ) C1, C2, C3
      REAL(SZ) EE1, EE2, EE3
      REAL(SZ) H2, H2N1, H2N2, H2N3
      REAL(SZ) U11, U22, U33
      REAL(SZ) V11, V22, V33
      LOGICAL, SAVE :: FirstCall = .true.
      LOGICAL, dimension(14) :: netcdf_error
C
      type(OutputDataDescript_t), save :: SigTStaDescript
      type(OutputDataDescript_t), save :: SalStaDescript
      type(OutputDataDescript_t), save :: TempStaDescript
      type(OutputDataDescript_t), save :: QSurfKp1Descript
      type(OutputDataDescript_t), save :: RealQStaDescript
      type(OutputDataDescript_t), save :: ImaginaryQStaDescript
      type(OutputDataDescript_t), save :: WZStaDescript
      type(OutputDataDescript_t), save :: Q20StaDescript
      type(OutputDataDescript_t), save :: LStaDescript
      type(OutputDataDescript_t), save :: EVStaDescript
      type(OutputDataDescript_t), save :: SigTDescript
      type(OutputDataDescript_t), save :: SalDescript
      type(OutputDataDescript_t), save :: TempDescript
      type(OutputDataDescript_t), save :: RealQDescript
      type(OutputDataDescript_t), save :: ImaginaryQDescript
      type(OutputDataDescript_t), save :: WZDescript
      type(OutputDataDescript_t), save :: Q20Descript
      type(OutputDataDescript_t), save :: LDescript
      type(OutputDataDescript_t), save :: EVDescript
C
C     parts of complex variables
      REAL(SZ), ALLOCATABLE, TARGET, SAVE :: rp(:,:) !real part subdomain data
      REAL(SZ), ALLOCATABLE, TARGET, SAVE :: ip(:,:) !imaginary part subdom dat
      REAL(SZ), ALLOCATABLE, TARGET, SAVE :: rp_g(:,:) !real part fulldomain
      REAL(SZ), ALLOCATABLE, TARGET, SAVE :: ip_g(:,:) !imaginary part fulldom
      REAL(SZ), ALLOCATABLE, TARGET, SAVE :: rpSta(:,:) !real part subd sta
      REAL(SZ), ALLOCATABLE, TARGET, SAVE :: ipSta(:,:) !im part subdom sta
      REAL(SZ), ALLOCATABLE, TARGET, SAVE :: rpSta_g(:,:) !real part fd sta
      REAL(SZ), ALLOCATABLE, TARGET, SAVE :: ipSta_g(:,:) !im part fd sta
C
      REAL(SZ), ALLOCATABLE, TARGET, SAVE :: qsurfkp1_g(:)
C
      REAL(SZ), ALLOCATABLE, TARGET, SAVE :: q20Sta_g(:,:)
      REAL(SZ), ALLOCATABLE, TARGET, SAVE :: lSta_g(:,:)
      REAL(SZ), ALLOCATABLE, TARGET, SAVE :: EVSta_g(:,:)
      REAL(SZ), ALLOCATABLE, TARGET, SAVE :: SigTSta_g(:,:)
      REAL(SZ), ALLOCATABLE, TARGET, SAVE :: SalSta_g(:,:)
      REAL(SZ), ALLOCATABLE, TARGET, SAVE :: TempSta_g(:,:)
      REAL(SZ), ALLOCATABLE, TARGET, SAVE :: WZSta_g(:,:)
C
      call setMessageSource("writeoutput3D")
#ifdef WRITE_OUTPUT_TRACE
      call allMessage(DEBUG,"Enter.")
#endif
C
1099  FORMAT(2X,1pE20.10E3,5X,I10)
1100  FORMAT(1X,E17.10,1X,I10,32000(2X,E13.6))
1104  FORMAT(9X,I6,4X,32000(E13.6,2X))
1105  FORMAT(2X,I8,2X,1pE20.10E3)
C 
      IF (FirstCall.eqv..true.) THEN
         ALLOCATE(rp(NP,NFEN),ip(NP,NFEN))
         ALLOCATE(rpSta(MNSTA3DV,NFEN),ipSta(MNSTA3DV,NFEN))
C
         SigTStaDescript % specifier            =  I3DSD
         SigTStaDescript % initial_value        =  0.0
         SigTStaDescript % num_items_per_record =  NFEN
         SigTStaDescript % num_fd_records       =  NSTA3DD_G
         SigTStaDescript % num_records_this     =  NSTA3DD
         SigTStaDescript % imap                 => IMAP_STA3DD_LG
         SigTStaDescript % array2D              => SigTSta
         SigTStaDescript % array2D_g            => SigTSta_g
         SigTStaDescript % considerWetDry       = .FALSE.
         SigTStaDescript % alternate_value      = -99999.0
         SigTStaDescript % field_name           = 'SigmaTStations'
         IF (ICS.eq.1) THEN
            SigTStaDescript % x_coord           => XED  ! use orig coord
            SigTStaDescript % y_coord           => YED  ! use orig coord
         ELSE
            SigTStaDescript % x_coord           => SLED  ! radians
            SigTStaDescript % y_coord           => SFED           
         ENDIF

         SalStaDescript % specifier            =  I3DSD
         SalStaDescript % initial_value        =  0.0
         SalStaDescript % num_items_per_record =  NFEN
         SalStaDescript % num_fd_records       =  NSTA3DD_G
         SalStaDescript % num_records_this     =  NSTA3DD
         SalStaDescript % imap                 => IMAP_STA3DD_LG
         SalStaDescript % array2D              => SalSta
         SalStaDescript % array2D_g            => SalSta_g
         SalStaDescript % considerWetDry       = .FALSE.
         SalStaDescript % alternate_value      = -99999.0
         SalStaDescript % field_name           = 'SalinityStations'
         IF (ICS.eq.1) THEN
            SalStaDescript % x_coord           => XED ! use orig coord
            SalStaDescript % y_coord           => YED ! use orig coord
         ELSE
            SalStaDescript % x_coord           => SLED  ! radians
            SalStaDescript % y_coord           => SFED           
         ENDIF

         TempStaDescript % specifier            =  I3DSD
         TempStaDescript % initial_value        =  0.0
         TempStaDescript % num_items_per_record =  NFEN
         TempStaDescript % num_fd_records       =  NSTA3DD_G
         TempStaDescript % num_records_this     =  NSTA3DD
         TempStaDescript % imap                 => IMAP_STA3DD_LG
         TempStaDescript % array2D              => TempSta
         TempStaDescript % array2D_g            => TempSta_g
         TempStaDescript % considerWetDry       = .FALSE.
         TempStaDescript % alternate_value      = -99999.0
         TempStaDescript % field_name           = 'TemperatureStations'
         IF (ICS.eq.1) THEN
            TempStaDescript % x_coord           => XED ! use orig coord
            TempStaDescript % y_coord           => YED
         ELSE
            TempStaDescript % x_coord           => SLED ! radians
            TempStaDescript % y_coord           => SFED
         ENDIF

         RealQStaDescript % specifier            =  I3DSV
         RealQStaDescript % initial_value        =  0.0
         RealQStaDescript % num_items_per_record =  NFEN
         RealQStaDescript % num_fd_records       =  NSta3DV_G
         RealQStaDescript % num_records_this     =  NSta3DV
         RealQStaDescript % imap                 => IMAP_STA3DV_LG
         RealQStaDescript % array2D              => rpSta
         RealQStaDescript % array2D_g            => rpSta_g
         RealQStaDescript % considerWetDry       = .FALSE.
         RealQStaDescript % alternate_value      = -99999.0
         RealQStaDescript % field_name           = 'RealQStations'
         IF (ICS.eq.1) THEN
            RealQStaDescript % x_coord           => XE3DV ! use orig coord
            RealQStaDescript % y_coord           => YE3DV
         ELSE
            RealQStaDescript % x_coord           => SLE3DV ! use orig coord
            RealQStaDescript % y_coord           => SFE3DV
         ENDIF

         ImaginaryQStaDescript % specifier            =  I3DSV
         ImaginaryQStaDescript % initial_value        =  0.0
         ImaginaryQStaDescript % num_items_per_record =  NFEN
         ImaginaryQStaDescript % num_fd_records       =  NSTA3DV_G
         ImaginaryQStaDescript % num_records_this     =  NSTA3DV
         ImaginaryQStaDescript % imap                 => IMAP_STA3DV_LG
         ImaginaryQStaDescript % array2D              => ipSta
         ImaginaryQStaDescript % array2D_g            => ipSta_g
         ImaginaryQStaDescript % considerWetDry       = .FALSE.
         ImaginaryQStaDescript % alternate_value      = -99999.0
         ImaginaryQStaDescript % field_name      ='ImaginaryQStations'
         IF (ICS.eq.1) THEN
            ImaginaryQStaDescript % x_coord           => XE3DV ! use orig coord
            ImaginaryQStaDescript % y_coord           => YE3DV
         ELSE
            ImaginaryQStaDescript % x_coord           => SLE3DV ! radians
            ImaginaryQStaDescript % y_coord           => SFE3DV
         ENDIF

         WZStaDescript % specifier            =  I3DSV
         WZStaDescript % initial_value        =  0.0
         WZStaDescript % num_items_per_record =  NFEN
         WZStaDescript % num_fd_records       =  NSTA3DV_G
         WZStaDescript % num_records_this     =  NSTA3DV
         WZStaDescript % imap                 => IMAP_STA3DV_LG
         WZStaDescript % array2D              => WZSta
         WZStaDescript % array2D_g            => WZSta_g
         WZStaDescript % considerWetDry       = .FALSE.
         WZStaDescript % alternate_value      = -99999.0
         WZStaDescript % field_name           = 'WZStations'
         IF (ICS.eq.1) THEN
            WZStaDescript % x_coord           => XE3DV  ! use orig coord
            WZStaDescript % y_coord           => YE3DV
         ELSE
            WZStaDescript % x_coord           => SLE3DV  ! radians
            WZStaDescript % y_coord           => SFE3DV         
         ENDIF

         Q20StaDescript % specifier            =  I3DST
         Q20StaDescript % initial_value        =  0.0
         Q20StaDescript % num_items_per_record =  NFEN
         Q20StaDescript % num_fd_records       =  NSta3DT_G
         Q20StaDescript % num_records_this     =  NSta3DT
         Q20StaDescript % imap                 => IMAP_STA3DT_LG
         Q20StaDescript % array2D              => q20Sta
         Q20StaDescript % array2D_g            => q20Sta_g
         Q20StaDescript % considerWetDry       = .FALSE.
         Q20StaDescript % alternate_value      = -99999.0
         Q20StaDescript % field_name           = 'q20Stations'
         IF (ICS.eq.1) THEN
            Q20StaDescript % x_coord           => XET  ! use orig coord
            Q20StaDescript % y_coord           => YET  ! use orig coord
         ELSE
            Q20StaDescript % x_coord           => SLET  ! radians
            Q20StaDescript % y_coord           => SFET  ! 
         ENDIF

         LStaDescript % specifier        =  I3DST
         LStaDescript % initial_value        =  0.0
         LStaDescript % num_items_per_record =  NFEN
         LStaDescript % num_fd_records       =  NSTA3DT_G
         LStaDescript % num_records_this     =  NSTA3DT
         LStaDescript % imap                 => IMAP_STA3DT_LG
         LStaDescript % array2D              => lSta
         LStaDescript % array2D_g            => lSta_g
         LStaDescript % considerWetDry       = .FALSE.
         LStaDescript % alternate_value      = -99999.0
         LStaDescript % field_name           = 'LStations'
         IF (ICS.eq.1) THEN
            LStaDescript % x_coord           => XET  ! use orig coord
            LStaDescript % y_coord           => YET
         ELSE
            LStaDescript % x_coord           => SLET  ! radians
            LStaDescript % y_coord           => SFET
         ENDIF

         EVStaDescript % specifier            =  I3DST
         EVStaDescript % initial_value        =  0.0
         EVStaDescript % num_items_per_record =  NFEN
         EVStaDescript % num_fd_records       =  NSTA3DT_G
         EVStaDescript % num_records_this     =  NSTA3DT
         EVStaDescript % imap                 => IMAP_STA3DT_LG
         EVStaDescript % array2D              => EVSta
         EVStaDescript % array2D_g            => EVSta_g
         EVStaDescript % considerWetDry       = .FALSE.
         EVStaDescript % alternate_value      = -99999.0
         EVStaDescript % field_name           = 'EVStations'
         IF (ICS.eq.1) THEN
            EVStaDescript % x_coord           => XET ! use orig coord
            EVStaDescript % y_coord           => YET
         ELSE
            EVStaDescript % x_coord           => SLET ! use orig coord
            EVStaDescript % y_coord           => SFET
         ENDIF

         QSurfKp1Descript % specifier            =  I3DGD
         QSurfKp1Descript % initial_value        =  0.0
         QSurfKp1Descript % num_items_per_record =  1
         QSurfKp1Descript % num_fd_records     =  NP_G
         QSurfKp1Descript % num_records_this   =  NP
         QSurfKp1Descript % imap               => NODES_LG
         QSurfKp1Descript % array              => qsurfkp1
         QSurfKp1Descript % array_g            => qsurfkp1_g
         QSurfKp1Descript % considerWetDry       = .FALSE.
         QSurfKp1Descript % alternate_value      = -99999.0
         QSurfKp1Descript % field_name           = 'qsurfkp1'

         SigTDescript % specifier            =  I3DGD
         SigTDescript % initial_value        =  0.0
         SigTDescript % num_items_per_record =  NFEN
         SigTDescript % num_fd_records        =  NP_G
         SigTDescript % num_records_this     =  NP
         SigTDescript % imap                 => nodes_lg
         SigTDescript % array2D              => SigT
         SigTDescript % array2D_g            => SigT_g
         SigTDescript % considerWetDry       = .FALSE.
         SigTDescript % alternate_value      = -99999.0
         SigTDescript % field_name           = 'SigmaT'

         SalDescript % specifier            =  I3DGD
         SalDescript % initial_value        =  0.0
         SalDescript % num_items_per_record =  NFEN
         SalDescript % num_fd_records       =  NP_G
         SalDescript % num_records_this     =  NP
         SalDescript % imap                 => nodes_lg
         SalDescript % array2D              => Sal
         SalDescript % array2D_g            => Sal_g
         SalDescript % considerWetDry       = .FALSE.
         SalDescript % alternate_value      = -99999.0
         SalDescript % field_name           = 'Salinity'

         TempDescript % specifier            =  I3DGD
         TempDescript % initial_value        =  0.0
         TempDescript % num_items_per_record =  NFEN
         TempDescript % num_fd_records       =  NP_G
         TempDescript % num_records_this     =  NP
         TempDescript % imap                 => nodes_lg
         TempDescript % array2D              => Temp
         TempDescript % array2D_g            => Temp_g
         TempDescript % considerWetDry       = .FALSE.
         TempDescript % alternate_value      = -99999.0
         TempDescript % field_name           = 'Temperature'

         RealQDescript % specifier            =  I3DGV
         RealQdescript % initial_value        =  0.0
         RealQDescript % num_items_per_record =  NFEN
         RealQDescript % num_fd_records       =  NP_G
         RealQDescript % num_records_this     =  NP
         RealQDescript % imap                 => nodes_lg
         RealQDescript % array2D              => rp
         RealQDescript % array2D_g            => rp_g
         RealQDescript % considerWetDry       = .FALSE.
         RealQDescript % alternate_value      = -99999.0
         RealQDescript % field_name           = 'RealQ'

         ImaginaryQDescript % specifier            =  I3DGV
         ImaginaryQDescript % initial_value        =  0.0
         ImaginaryQDescript % num_items_per_record =  NFEN
         ImaginaryQDescript % num_fd_records       =  NP_G
         ImaginaryQDescript % num_records_this     =  NP
         ImaginaryQDescript % imap                 => nodes_lg
         ImaginaryQDescript % array2D              => ip
         ImaginaryQDescript % array2D_g            => ip_g
         ImaginaryQDescript % considerWetDry       = .FALSE.
         ImaginaryQDescript % alternate_value      = -99999.0
         ImaginaryQDescript % field_name           = 'ImaginaryQ'

         WZDescript % specifier            =  I3DGV
         WZDescript % initial_value        =  0.0
         WZDescript % num_items_per_record =  NFEN
         WZDescript % num_fd_records       =  NP_G
         WZDescript % num_records_this     =  NP
         WZDescript % imap                 => nodes_lg
         WZDescript % array2D              => WZ
         WZDescript % array2D_g            => WZ_g
         WZDescript % considerWetDry       = .FALSE.
         WZDescript % alternate_value      = -99999.0
         WZDescript % field_name           = 'WZ'

         Q20Descript % specifier            =  I3DGT
         Q20Descript % initial_value        =  0.0
         Q20Descript % num_items_per_record =  NFEN
         Q20Descript % num_fd_records       =  NP_G
         Q20Descript % num_records_this     =  NP
         Q20Descript % imap                 => nodes_lg
         Q20Descript % array2D              => q20
         Q20Descript % array2D_g            => q20_g
         Q20Descript % considerWetDry       = .FALSE.
         Q20Descript % alternate_value      = -99999.0
         Q20Descript % field_name           = 'q20'

         LDescript % specifier            =  I3DGT
         LDescript % initial_value        =  0.0
         LDescript % num_items_per_record =  NFEN
         LDescript % num_fd_records       =  NP_G
         LDescript % num_records_this     =  NP
         LDescript % imap                 => nodes_lg
         LDescript % array2D              => l
         LDescript % array2D_g            => l_g
         LDescript % considerWetDry       = .FALSE.
         LDescript % alternate_value      = -99999.0
         LDescript % field_name           = 'L'

         EVDescript % specifier            =  I3DGT
         EVDescript % initial_value        =  0.0
         EVDescript % num_items_per_record =  NFEN
         EVDescript % num_fd_records       =  NP_G
         EVDescript % num_records_this     =  NP
         EVDescript % imap                 => nodes_lg
         EVDescript % array2D              => EV
         EVDescript % array2D_g            => EV_g
         EVDescript % considerWetDry       = .FALSE.
         EVDescript % alternate_value      = -99999.0
         EVDescript % field_name           = 'EV'

#ifdef ADCNETCDF
         netcdf_error = .false.
         ! jgf49.48.01: if netcdf is not required for a particular output file,
         ! the associated subroutine call will do nothing and just return
         IF (myProc.eq.0) THEN
            CALL initNetCDFOutputFile(41, SigTStaDescript,
     &           netcdf_error(1), SalStaDescript, TempStaDescript)
            CALL initNetCDFOutputFile(42, RealQStaDescript,
     &           netcdf_error(2), ImaginaryQStaDescript, WZStaDescript)
            CALL initNetCDFOutputFile(43, Q20StaDescript,
     &           netcdf_error(3), LStaDescript, EVStaDescript)
            CALL initNetCDFOutputFile(44, SigTDescript,
     &           netcdf_error(4), SalDescript, TempDescript)
            CALL initNetCDFOutputFile(45, RealQDescript,
     &           netcdf_error(5), ImaginaryQDescript, WZDescript)
            CALL initNetCDFOutputFile(46, Q20Descript,
     &           netcdf_error(6), LDescript, EVDescript)
            CALL initNetCDFOutputFile(47, QSurfKp1Descript,
     &           netcdf_error(7))
         ENDIF
#ifdef CMPI
         CALL msg_lbcast(netcdf_error,7)
         DO i=1,7
            IF (netcdf_error(i).eqv..true.) THEN
               CALL msg_fini()
               STOP
            ENDIF
         ENDDO
#endif
#endif
         FirstCall = .false.
      ENDIF
C
C.... 3D Density, Temperature, Salinity Station Output (Unit 41)
C
      IF(I3DSD.NE.0) THEN
C        Check to see if it is time to generate 3D station density output.
         IF((IT.GT.NTO3DSDS).AND.(IT.LE.NTO3DSDF)) N3DSD=N3DSD+1
         IF(N3DSD.EQ.NSpo3DSD) THEN
C           Process each station
            DO NN=1,NSta3DD
C              Interpolate to station locations
               NEle=NE3DD(NN)
               N1=NM(NEle,1)
               N2=NM(NEle,2)
               N3=NM(NEle,3)
               DO k=1,NFEN
                  SigTSta(NN,k)=SigT(N1,k)*StaI3DD1(NN)
     &              +SigT(N2,k)*StaI3DD2(NN)+SigT(N3,k)*StaI3DD3(NN)
                  IF ((ABS(IDen).EQ.2).OR.(ABS(IDen).EQ.4)) THEN
                     SalSta(NN,k) = Sal(N1,k)*StaI3DD1(NN)
     &                             +Sal(N2,k)*StaI3DD2(NN)
     &                             +Sal(N3,k)*StaI3DD3(NN)
                  ENDIF
                  IF ((ABS(IDen).EQ.3).OR.(ABS(IDen).EQ.4)) THEN
                     TempSta(NN,k)= Temp(N1,k)*StaI3DD1(NN)
     &                             +Temp(N2,k)*StaI3DD2(NN)
     &                             +Temp(N3,k)*StaI3DD3(NN)
                  ENDIF
               ENDDO
            END DO

            ! jgf49.43.19 Collect fulldomain data in parallel.
            IF ((MNPROC.GT.1).AND.(WRITE_LOCAL_FILES.eqv..false.)) THEN
               IF (myProc.eq.0)  THEN
                  ALLOCATE(SigTSta_g(NSTA3DD_G,NFEN))
                  SigTStaDescript % array2D_g => SigTSta_g
               ENDIF
               CALL collectFullDomainArray(SigTStaDescript,
     &                  packNPbyM, unpackNPbyM)

               IF ((ABS(IDEN).EQ.2).OR.(ABS(IDEN).EQ.4)) THEN
                  IF (myProc.eq.0) THEN
                     ALLOCATE(SalSta_g(NSTA3DD_G,NFEN))
                     SalStaDescript % array2D_g => SalSta_g
                  ENDIF
                  CALL collectFullDomainArray(SalStaDescript,
     &                     packNPbyM, unpackNPbyM)
               ENDIF
               IF ((ABS(IDEN).EQ.3).OR.(ABS(IDEN).EQ.4)) THEN
                  IF (myProc.eq.0) THEN
                     ALLOCATE(TempSta_g(NSTA3DD_G,NFEN))
                     !kmd : Took out these two lines as already allocated above
                     !ALLOCATE(SalSta_g(NSTA3DD_G,NFEN))
                     !SalStaDescript % array2D_g => SalSta_g
                     TempStaDescript % array2D_g => TempSta_g
                  ENDIF
                  CALL collectFullDomainArray(TempStaDescript,
     &                     packNPbyM, unpackNPbyM)
               ENDIF
            ENDIF
            SELECT CASE(ABS(I3DSD))
            CASE(1) ! ascii
               ! serial ascii
               IF ((MNPROC.EQ.1).OR.(WRITE_LOCAL_FILES.eqv..true.)) THEN
                  OPEN(41,FILE=TRIM(LOCALDIR)//'/'//'fort.41',
     &                            ACCESS='SEQUENTIAL',POSITION='APPEND')
C                 Write time header into file
                  SELECT CASE(ABS(IDEN))
                  CASE(0,1)
                     WRITE(41,1100) TimeLoc,IT,(Sigma(k),k=1,NFEN-1)
                     I3DSDRec=I3DSDRec+1
                     DO NN=1,NSta3DD
                        WRITE(41,1104) NN,(SigTSta(NN,k),k=1,NFEN)
                        I3DSDRec=I3DSDRec+1
                     END DO
                  CASE(2)
                     WRITE(41,1100) TimeLoc,IT,(Sigma(k),Sigma(k),
     &                               k=1,NFEN-1),Sigma(NFEN)
                     I3DSDRec=I3DSDRec+1
                     DO NN=1,NSta3DD
                        WRITE(41,1104) NN,
     &                     (SigTSta(NN,k),SalSta(NN,k),k=1,NFEN)
                        I3DSDRec=I3DSDRec+1
                     END DO
                  CASE(3)
                     WRITE(41,1100) TimeLoc,IT,(Sigma(k),Sigma(k),
     &                               k=1,NFEN-1),Sigma(NFEN)
                     I3DSDRec=I3DSDRec+1
                     DO NN=1,NSta3DD
                        WRITE(41,1104) NN,
     &                     (SigTSta(NN,k),TempSta(NN,k),k=1,NFEN)
                        I3DSDRec=I3DSDRec+1
                     END DO
                  CASE(4)
                     WRITE(41,1100) TimeLoc,IT,(Sigma(k),Sigma(k),
     &                      Sigma(k),k=1,NFEN-1),Sigma(NFEN),Sigma(NFEN)
                     I3DSDRec=I3DSDRec+1
                     DO NN=1,NSta3DD
                        WRITE(41,1104) NN,(SigTSta(NN,k),
     &                      SalSta(NN,k),TempSta(NN,k),k=1,NFEN)
                        I3DSDRec=I3DSDRec+1
                     END DO
                  CASE DEFAULT
                     CALL allMessage(WARNING,
     &                 "IDEN is not +/- 0 to 4 and is invalid.")
                  END SELECT
               ENDIF
               ! globalio ascii
               IF ((MNPROC.GT.1).AND.(WRITE_LOCAL_FILES.eqv..false.).AND.
     &             (myProc.eq.0)) THEN
                  OPEN(41,FILE=TRIM(GLOBALDIR)//'/'//'fort.41',
     &                            ACCESS='SEQUENTIAL',POSITION='APPEND')
C                 Write time header into file
                  SELECT CASE(ABS(IDEN))
                  CASE(0,1)
                     WRITE(41,1100) TimeLoc,IT,(Sigma(k),k=1,NFEN-1)
                     I3DSDRec=I3DSDRec+1
                     DO NN=1,NSta3DD_G
                        WRITE(41,1104) NN,(SigTSta_g(NN,k),k=1,NFEN)
                        I3DSDRec=I3DSDRec+1
                     END DO
                  CASE(2)
                     WRITE(41,1100) TimeLoc,IT,(Sigma(k),Sigma(k),
     &                               k=1,NFEN-1),Sigma(NFEN)
                     I3DSDRec=I3DSDRec+1
                     DO NN=1,NSta3DD_G
                        WRITE(41,1104) NN,
     &                     (SigTSta_g(NN,k),SalSta_g(NN,k),k=1,NFEN)
                        I3DSDRec=I3DSDRec+1
                     END DO
                  CASE(3)
                     WRITE(41,1100) TimeLoc,IT,(Sigma(k),Sigma(k),
     &                               k=1,NFEN-1),Sigma(NFEN)
                     I3DSDRec=I3DSDRec+1
                     DO NN=1,NSta3DD_G
                        WRITE(41,1104) NN,
     &                     (SigTSta_g(NN,k),TempSta_g(NN,k),k=1,NFEN)
                        I3DSDRec=I3DSDRec+1
                     END DO
                  CASE(4)
                     WRITE(41,1100) TimeLoc,IT,(Sigma(k),Sigma(k),
     &                      Sigma(k),k=1,NFEN-1),Sigma(NFEN),Sigma(NFEN)
                     I3DSDRec=I3DSDRec+1
                     DO NN=1,NSta3DD_G
                        WRITE(41,1104) NN,(SigTSta_g(NN,k),
     &                      SalSta_g(NN,k),TempSta_g(NN,k),k=1,NFEN)
                        I3DSDRec=I3DSDRec+1
                     END DO
                  CASE DEFAULT
                     CALL allMessage(WARNING,
     &                 "IDEN is not +/- 0 to 4 and is invalid.")
                  END SELECT
               ENDIF
            CASE(2) ! binary
               OPEN(41,FILE=TRIM(LOCALDIR)//'/'//'fort.41',
     &            ACCESS='DIRECT', RECL=NByte)
               WRITE(41,REC=I3DSDRec+1) TimeLoc
               WRITE(41,REC=I3DSDRec+2) IT
               I3DSDRec=I3DSDRec+2
C              Write station output
               DO NN=1,NSTA3DD
                  WRITE(41,REC=I3DSDRec+1) NN
                  I3DSDRec=I3DSDRec+1
                  DO k=1,NFEN
                     SELECT CASE(ABS(IDEN))
                     CASE(1)
                        WRITE(41,REC=I3DSDRec+1) SigTSta(NN,k)
                        I3DSDRec=I3DSDRec+1
                     CASE(2)
                        WRITE(41,REC=I3DSDRec+1) SigTSta(NN,k)
                        WRITE(41,REC=I3DSDRec+2) SalSta(NN,k)
                        I3DSDRec=I3DSDRec+2
                     CASE(3)
                        WRITE(41,REC=I3DSDRec+1) SigTSta(NN,k)
                        WRITE(41,REC=I3DSDRec+2) TempSta(NN,k)
                        I3DSDRec=I3DSDRec+2
                     CASE(4)
                        WRITE(41,REC=I3DSDRec+1) SigTSta(NN,k)
                        WRITE(41,REC=I3DSDRec+2) SalSta(NN,k)
                        WRITE(41,REC=I3DSDRec+3) TempSta(NN,k)
                        I3DSDRec=I3DSDRec+3
                     CASE DEFAULT
                        CALL allMessage(WARNING,
     &                    "IDEN is not +/- 0 to 4 and is invalid.")
                     END SELECT
                  ENDDO
               ENDDO
               CLOSE(41)
#ifdef ADCNETCDF
            CASE(3,5) ! netcdf
               IF (myProc.eq.0) THEN
                  CALL writeOutArrayNetCDF(41, TimeLoc, it,
     &             SigTStaDescript, SalStaDescript, TempStaDescript)
               ENDIF
#endif
            CASE DEFAULT
               CALL allMessage(ERROR,
     &            "ABS(I3DSD) must be either 0, 1, 2, 3, or 5.")
            END SELECT
            N3DSD=0
            ! deallocate memory used in globalio
            IF ((MNPROC.GT.1).AND.(WRITE_LOCAL_FILES.eqv..false.).AND.
     &          (myProc.eq.0)) THEN
               DEALLOCATE(SigTSta_g)
               IF ((ABS(IDEN).EQ.2).OR.(ABS(IDEN).EQ.4)) THEN
                  DEALLOCATE(SalSta_g)
               ENDIF
               IF ((ABS(IDEN).EQ.3).OR.(ABS(IDEN).EQ.4)) THEN
                  DEALLOCATE(TempSta_g)
                  !kmd : Changed Deallocate line above as Salinity is already deallocated
                  !DEALLOCATE(TempSta_g,SalSta_g)
               ENDIF
            ENDIF

         ENDIF
      ENDIF

C
C.... 3D Velocity Station Output (Unit 42)
C
      IF(I3DSV.NE.0) THEN
C        Check to see if it is time to generate 3D station velocity output.
         IF((IT.GT.NTO3DSVS).AND.(IT.LE.NTO3DSVF)) N3DSV=N3DSV+1
         IF(N3DSV.EQ.NSpo3DSV) THEN
C           Process each station
            DO NN=1,NSta3DV
C              Interpolate to station locations
               NEle=NE3DV(NN)
               N1=NM(NEle,1)
               N2=NM(NEle,2)
               N3=NM(NEle,3)
               DO k=1,NFEN
                  qSta(NN,k) =q(N1,k)*StaI3DV1(NN)
     &                       +q(N2,k)*StaI3DV2(NN)
     &                       +q(N3,k)*StaI3DV3(NN)
                  WZSta(NN,k)=WZ(N1,k)*StaI3DV1(NN)
     &                       +WZ(N2,k)*StaI3DV2(NN)
     &                       +WZ(N3,k)*StaI3DV3(NN)
               ENDDO
            END DO

            ! jgf49.43.19 Collect fulldomain data in parallel.
            IF ((MNPROC.GT.1).AND.(WRITE_LOCAL_FILES.eqv..false.)) THEN
               IF (myProc.eq.0) THEN
                  ALLOCATE(rpSta_g(NSTA3DV_G,NFEN),
     &               ipSta_g(NSTA3DV_G,NFEN))
                  ALLOCATE(WZSta_g(NSTA3DV_G,NFEN))
                  RealQStaDescript % array2D_g => rpSta_g
                  ImaginaryQStaDescript % array2D_g => ipSta_g
                  WZStaDescript % array2D_g => WZSta_g
               ENDIF
               DO NN=1,NSta3DV
                  DO k=1,NFEN
                     rpSta(nn,k) = real(qSta(nn,k))
                     ipSta(nn,k) = aimag(qSta(nn,k))
                  END DO
               END DO
               CALL collectFullDomainArray(RealQStaDescript,
     &                     packNPbyM, unpackNPbyM)
               CALL collectFullDomainArray(ImaginaryQStaDescript,
     &                     packNPbyM, unpackNPbyM)
               CALL collectFullDomainArray(WZStaDescript,
     &                     packNPbyM, unpackNPbyM)
            ENDIF
            SELECT CASE(ABS(I3DSV))
            CASE(1)  !ASCII FORMAT
               IF ((MNPROC.EQ.1).OR.(WRITE_LOCAL_FILES.eqv..true.)) THEN
                  OPEN(42,FILE=TRIM(LOCALDIR)//'/'//'fort.42',
     &                            ACCESS='SEQUENTIAL',POSITION='APPEND')
C                 Write time header into file
                  WRITE(42,1100) TimeLoc,IT,(Sigma(k),Sigma(k),Sigma(k),
     &                       k=1,NFEN-1),Sigma(NFEN),Sigma(NFEN)
                  I3DSVRec=I3DSVRec+1
                  DO NN=1,NSTA3DV
                     WRITE(42,1104) NN,
     &                  (REAL(qSta(NN,k)),AIMAG(qSta(NN,k)),
     &                  WZSta(NN,k),k=1,NFEN)
                     I3DSVRec=I3DSVRec+1
                  ENDDO
                  CLOSE(42)
               ENDIF
               IF ((MNPROC.GT.1).AND.(WRITE_LOCAL_FILES.eqv..false.).AND.
     &             (myProc.eq.0)) THEN
                  OPEN(42,FILE=TRIM(GLOBALDIR)//'/'//'fort.42',
     &                            ACCESS='SEQUENTIAL',POSITION='APPEND')
C                 Write time header into file
                  WRITE(42,1100) TimeLoc,IT,(Sigma(k),Sigma(k),Sigma(k),
     &                       k=1,NFEN-1),Sigma(NFEN),Sigma(NFEN)
                  I3DSVRec=I3DSVRec+1
                  DO NN=1,NSTA3DV_G
                     WRITE(42,1104) NN,(rpSta_g(NN,k),ipSta_g(NN,k),
     &                  WZSta_g(NN,k),k=1,NFEN)
                     I3DSVRec=I3DSVRec+1
                  ENDDO
                  CLOSE(42)
               ENDIF
            CASE(2) !BINARY FORMAT
               OPEN(42,FILE=TRIM(LOCALDIR)//'/'//'fort.42',
     &            ACCESS='DIRECT', RECL=NByte)
               WRITE(42,REC=I3DSVRec+1) TimeLoc
               WRITE(42,REC=I3DSVRec+2) IT
               I3DSVRec=I3DSVRec+2
               WRITE(42,REC=I3DSVRec+1) NN
               I3DSVRec=I3DSVRec+1
               DO NN=1,NSTA3DV
                  DO k=1,NFEN
                     WRITE(42,REC=I3DSVRec+1) REAL(qSta(NN,k))
                     WRITE(42,REC=I3DSVRec+2) AIMAG(qSta(NN,k))
                     WRITE(42,REC=I3DSVRec+3) WZSta(NN,k)
                     I3DSVRec=I3DSVRec+3
                  END DO
               END DO
               CLOSE(42)
#ifdef ADCNETCDF
            CASE(3,5) ! netcdf
               IF (myProc.eq.0) THEN
                  CALL writeOutArrayNetCDF(42, TimeLoc, it,
     &               RealQStaDescript,ImaginaryQStaDescript,
     &               WZStaDescript)
               ENDIF
#endif
            CASE DEFAULT
               CALL allMessage(ERROR,
     &            "ABS(I3DSV) must be either 0, 1, 2, 3 or 5.")
            END SELECT
            N3DSV=0
            ! deallocate memory for globalio.
            IF ((MNPROC.GT.1).AND.(WRITE_LOCAL_FILES.eqv..false.).AND.
     &           (myProc.eq.0)) THEN
               DEALLOCATE(rpSta_g,ipSta_g,WZSta_g)
            ENDIF
        ENDIF
      ENDIF

C
C.... 3D Turbulence Station Output (Unit 43)
C
      IF(I3DST.NE.0) THEN
C        Check to see if it is time to generate 3D station turbulence output.
         IF ((IT.GT.NTO3DSTS).AND.(IT.LE.NTO3DSTF)) N3DST=N3DST+1
         IF (N3DST.EQ.NSpo3DST) THEN
C           Process each station
            DO NN=1,NSta3DT
C              Interpolate to station locations
               NEle=NE3DT(NN)
               N1=NM(NEle,1)
               N2=NM(NEle,2)
               N3=NM(NEle,3)
               DO k=1,NFEN
                  q20Sta(NN,k)=q20(N1,k)*StaI3DT1(NN)
     &               +q20(N2,k)*StaI3DT2(NN)+q20(N3,k)*StaI3DT3(NN)
                  lSta(NN,k)  =l(N1,k)  *StaI3DT1(NN)
     &               +l(N2,k)  *StaI3DT2(NN)+l(N3,k)  *StaI3DT3(NN)
                  EVSta(NN,k) =EV(N1,k) *StaI3DT1(NN)
     &               +EV(N2,k) *StaI3DT2(NN)+EV(N3,k) *StaI3DT3(NN)
               END DO
            END DO

            ! jgf49.43.18 Collect fulldomain data in parallel.
            IF ((MNPROC.GT.1).AND.(WRITE_LOCAL_FILES.eqv..false.)) THEN
               IF (myProc.eq.0) THEN
                  ALLOCATE(q20Sta_g(NSTA3DT_G,NFEN))
                  ALLOCATE(lSta_g(NSTA3DT_G,NFEN))
                  ALLOCATE(EVSta_g(NSTA3DT_G,NFEN))
                  Q20StaDescript % array2D_g => q20Sta_g
                  LStaDescript % array2D_g => lSta_g
                  EVStaDescript % array2D_g => EVSta_g
               ENDIF
               CALL collectFullDomainArray(Q20StaDescript,
     &                     packNPbyM, unpackNPbyM)
               CALL collectFullDomainArray(LStaDescript,
     &                     packNPbyM, unpackNPbyM)
               CALL collectFullDomainArray(EVStaDescript,
     &                     packNPbyM, unpackNPbyM)
            ENDIF
            SELECT CASE(ABS(I3DST))
            CASE(1)  !ASCII FORMAT
               IF ((MNPROC.EQ.1).OR.(WRITE_LOCAL_FILES.eqv..true.)) THEN
                  OPEN(43,FILE=TRIM(LOCALDIR)//'/'//'fort.43',
     &                ACCESS='SEQUENTIAL',POSITION='APPEND')
C                 Write time header into file
                  WRITE(43,1100) TimeLoc,IT,(Sigma(k),Sigma(k),Sigma(k),
     &                k=1,NFEN-1),Sigma(NFEN),Sigma(NFEN)
                  I3DSTRec=I3DSTRec+1
                  DO NN=1,NSta3DT
                     WRITE(43,1104)
     &                 NN,(q20Sta(NN,k),lSta(NN,k),EVSta(NN,k),k=1,NFEN)
                     I3DSTRec=I3DSTRec+1
                  END DO
                  CLOSE(43)
               ENDIF
               IF ((MNPROC.GT.1).AND.(WRITE_LOCAL_FILES.eqv..false.).AND.
     &              (myProc.eq.0)) THEN
                  OPEN(43,FILE=TRIM(GLOBALDIR)//'/'//'fort.43',
     &                ACCESS='SEQUENTIAL',POSITION='APPEND')
C                 Write time header into file
                  WRITE(43,1100) TimeLoc,IT,(Sigma(k),Sigma(k),Sigma(k),
     &                k=1,NFEN-1),Sigma(NFEN),Sigma(NFEN)
                  I3DSTRec=I3DSTRec+1
                  DO NN=1,NSta3DT_G
                     WRITE(43,1104) NN,(q20Sta_g(NN,k),lSta_g(NN,k),
     &                  EVSta_g(NN,k),k=1,NFEN)
                     I3DSTRec=I3DSTRec+1
                  END DO
                  CLOSE(43)
               ENDIF
            CASE(2)  ! binary
               OPEN(43,FILE=TRIM(LOCALDIR)//'/'//'fort.43',
     &            ACCESS='DIRECT', RECL=NByte)
               WRITE(43,REC=I3DSTRec+1) TimeLoc
               WRITE(43,REC=I3DSTRec+2) IT
               I3DSTRec=I3DSTRec+2
               DO NN=1,NSta3DT
                  WRITE(43,REC=I3DSTRec+1) NN
                  I3DSTRec=I3DSTRec+1
                  DO k=1,NFEN
                     WRITE(43,REC=I3DSTRec+1) q20Sta(NN,k)
                     WRITE(43,REC=I3DSTRec+2) lSta(NN,k)
                     WRITE(43,REC=I3DSTRec+3) EVSta(NN,k)
                     I3DSTRec=I3DSTRec+3
                  END DO
               END DO
               CLOSE(43)
#ifdef ADCNETCDF
            CASE(3,5) ! netcdf
               IF (myProc.eq.0) THEN
                  CALL writeOutArrayNetCDF(43, TimeLoc, it,
     &               Q20StaDescript,LStaDescript,EVStaDescript)
               ENDIF
#endif
            CASE DEFAULT
               CALL allMessage(ERROR,
     &            "ABS(I3DST) must be either 0, 1, 2, 3 or 5.")
            END SELECT
            N3DST=0
            ! deallocate memory for globalio
            IF ((MNPROC.GT.1).AND.(WRITE_LOCAL_FILES.eqv..false.).AND.
     &            (myProc.eq.0)) THEN
               DEALLOCATE(q20Sta_g,lSta_g,EVSta_g)
            ENDIF

         ENDIF
      ENDIF

C
C.... 3D Density, Temperature, Salinity Global Output (Unit 44)
C
      IF(I3DGD.NE.0) THEN
         IF((IT.GT.NTO3DGDS).AND.(IT.LE.NTO3DGDF)) N3DGD=N3DGD+1
C        Check to see if it is time to generate 3D fulldomain density output.
         IF(N3DGD.EQ.NSpo3DGD) THEN
            ! collect up output data in parallel if fulldomain files are needed
            IF ((MNPROC.gt.1).and.(WRITE_LOCAL_FILES.eqv..false.)) THEN
               IF (myProc.eq.0) THEN
                  ALLOCATE(SigT_g(NP_G,NFEN))
                  SigTDescript % array2D_g => SigT_g
               ENDIF
               CALL collectFullDomainArray(SigTDescript,
     &            packNPbyM, unpackNPbyM)
               SELECT CASE(ABS(IDEN))
               CASE(0) ! barotropic
                  CALL allMessage(WARNING,
     &            "IDEN was set to 0 but density output was requested.")
               CASE(1)
                  ! already collected SigT above
               CASE(2)
                  IF (myProc.eq.0) THEN
                     ALLOCATE(Sal_g(NP_G,NFEN))
                     SalDescript % array2D_g => Sal_g
                  ENDIF
                  CALL collectFullDomainArray(SalDescript,
     &                  packNPbyM, unpackNPbyM)
               CASE(3)
                  IF (myProc.eq.0) THEN
                     ALLOCATE(Temp_g(NP_G,NFEN))
                     TempDescript % array2D_g => Temp_g
                  ENDIF
                  CALL collectFullDomainArray(TempDescript,
     &                  packNPbyM, unpackNPbyM)
               CASE(4)
                  IF (myProc.eq.0) THEN
                     ALLOCATE(Sal_g(NP_G,NFEN),Temp_g(NP_G,NFEN))
                     SalDescript % array2D_g => Sal_g
                     TempDescript % array2D_g => Temp_g
                  ENDIF
                  CALL collectFullDomainArray(SalDescript,
     &                  packNPbyM, unpackNPbyM)
                  CALL collectFullDomainArray(TempDescript,
     &                  packNPbyM, unpackNPbyM)
               CASE DEFAULT
                  CALL allMessage(WARNING,
     &            "IDEN is not +/- 0 to 4 and is invalid.")
               END SELECT
            ENDIF
!kmd48.33bc - add in the information for outputting the top temperature
!             boundary condition.
            SELECT CASE(ABS(I3DGD))
            CASE(1)  ! ASCII
               IF ((MNPROC.eq.1).or.(WRITE_LOCAL_FILES.eqv..true.)) THEN
                  OPEN(44,FILE=TRIM(LOCALDIR)//'/'//'fort.44',
     &                ACCESS='SEQUENTIAL',POSITION='APPEND')
                  IF((IDEN.EQ.3).OR.(IDEN.EQ.4)) THEN
                     IF(BCFLAG_TEMP.NE.0) THEN
                        OPEN(47,FILE=TRIM(LOCALDIR)//'/'//'fort.47',
     &                         ACCESS='SEQUENTIAL',POSITION='APPEND')
                        WRITE(47,1099) TimeLoc,IT
                        DO NH=1,NP
                           WRITE(47,1105) NH, qsurfkp1(NH)
                        END DO
                        CLOSE(47)
                     END IF
                  END IF
                  SELECT CASE(ABS(IDEN))
                  CASE(0) ! barotropic
                     CALL allMessage(WARNING,
     &            "IDEN was set to 0 but density output was requested.")
                  CASE(1) ! SigmaT density
                     WRITE(44,1100) TimeLoc,IT,(Sigma(k),k=1,NFEN-1)
                     DO NH=1,NP
                        WRITE(44,1104) NH,(SigT(NH,k),k=1,NFEN)
                        I3DGDRec=I3DGDRec+1
                     ENDDO
                  CASE(2) ! Salinity density
                     WRITE(44,1100) TimeLoc,IT,(Sigma(k),Sigma(k),
     &                                k=1,NFEN-1),Sigma(NFEN)
                     DO NH=1,NP
                        WRITE(44,1104)
     &                     NH,(SigT(NH,k),Sal(NH,k),k=1,NFEN)
                           I3DGDRec=I3DGDRec+1
                     ENDDO
                  CASE(3) ! Temperature density
                     WRITE(44,1100) TimeLoc,IT,(Sigma(k),Sigma(k),
     &                                k=1,NFEN-1),Sigma(NFEN)
                     DO NH=1,NP
                        WRITE(44,1104)
     &                     NH,(SigT(NH,k),Temp(NH,k),k=1,NFEN)
                        I3DGDRec=I3DGDRec+1
                     ENDDO
                  CASE(4) ! Salinity+Temperature density
                     WRITE(44,1100) TimeLoc,IT,(Sigma(k),Sigma(k),
     &                    Sigma(k),k=1,NFEN-1),Sigma(NFEN),Sigma(NFEN)
                     DO NH=1,NP
                        WRITE(44,1104)
     &                     NH,(SigT(NH,k),Temp(NH,k),Sal(NH,k),k=1,NFEN)
                        I3DGDRec=I3DGDRec+1
                     ENDDO
                  CASE DEFAULT
                     CALL allMessage(WARNING,
     &               "IDEN is not +/- 0 to 4 and is invalid.")
                  END SELECT
                  I3DGDRec=I3DGDRec+1 ! for the header line
                  CLOSE(44)
               ENDIF
               IF ((MNPROC.GT.1).AND.(WRITE_LOCAL_FILES.eqv..false.)) THEN
                  IF((IDEN.EQ.3).OR.(IDEN.EQ.4)) THEN
                     IF(BCFLAG_TEMP.NE.0) THEN
                        IF (myProc.eq.0) THEN
                           ALLOCATE(qsurfkp1_g(NP_G))
                        ENDIF
                        CALL collectFullDomainArray(QSurfKp1Descript,
     &                     packOne, unpackOne)
                        OPEN(47,FILE=TRIM(GLOBALDIR)//'/'//'fort.47',
     &                     ACCESS='SEQUENTIAL',POSITION='APPEND')
                        WRITE(47,1099) TimeLoc,IT
                        DO NH=1,NP_G
                           WRITE(47,1105) NH, qsurfkp1_g(NH)
                        END DO
                        DEALLOCATE(qsurfkp1_g)
                        CLOSE(47)
                     END IF
                  END IF
                  OPEN(44,FILE=TRIM(GLOBALDIR)//'/'//'fort.44',
     &                ACCESS='SEQUENTIAL',POSITION='APPEND')
                  SELECT CASE(ABS(IDEN))
                  CASE(0) ! barotropic
                     CALL allMessage(WARNING,
     &            "IDEN was set to 0 but density output was requested.")
                  CASE(1) ! SigmaT density
                     IF (myProc.eq.0) THEN
                        WRITE(44,1100) TimeLoc,IT,(Sigma(k),k=1,NFEN-1)
                        DO NH=1,NP_G
                           WRITE(44,1104) NH,(SigT_g(NH,k),k=1,NFEN)
                           I3DGDRec=I3DGDRec+1
                        ENDDO
                     ENDIF
                  CASE(2) ! Salinity density
                     IF (myProc.eq.0) THEN
                        WRITE(44,1100) TimeLoc,IT,(Sigma(k),Sigma(k),
     &                                   k=1,NFEN-1),Sigma(NFEN)
                        DO NH=1,NP_G
                           WRITE(44,1104)
     &                        NH,(SigT_g(NH,k),Sal_g(NH,k),k=1,NFEN)
                              I3DGDRec=I3DGDRec+1
                        ENDDO
                     ENDIF
                  CASE(3) ! Temperature density
                     IF (myProc.eq.0) THEN
                        WRITE(44,1100) TimeLoc,IT,(Sigma(k),Sigma(k),
     &                                   k=1,NFEN-1),Sigma(NFEN)
                        DO NH=1,NP_G
                           WRITE(44,1104)
     &                        NH,(SigT_g(NH,k),Temp_g(NH,k),k=1,NFEN)
                           I3DGDRec=I3DGDRec+1
                        ENDDO
                     ENDIF
                  CASE(4) ! Salinity+Temperature density
                     IF (myProc.eq.0) THEN
                        WRITE(44,1100) TimeLoc,IT,(Sigma(k),Sigma(k),
     &                      Sigma(k),k=1,NFEN-1),Sigma(NFEN),Sigma(NFEN)
                        DO NH=1,NP_G
                           WRITE(44,1104) NH,(SigT_g(NH,k),Temp_g(NH,k),
     &                          Sal_g(NH,k),k=1,NFEN)
                           I3DGDRec=I3DGDRec+1
                        ENDDO
                     ENDIF
                  CASE DEFAULT
                     CALL allMessage(WARNING,
     &               "IDEN is not +/- 0 to 4 and is invalid.")
                  END SELECT
                  I3DGDRec=I3DGDRec+1 ! for the header line
                  CLOSE(44)
               ENDIF
            CASE(2)  !BINARY FORMAT
               OPEN(44,FILE=TRIM(LOCALDIR)//'/'//'fort.44',
     &            ACCESS='DIRECT', RECL=NByte)
               WRITE(44,REC=I3DGDRec+1) TimeLoc
               WRITE(44,REC=I3DGDRec+2) IT
               I3DGDRec=I3DGDRec+2
               SELECT CASE(ABS(IDEN))
               CASE(0) ! barotropic
                  CALL allMessage(WARNING,
     &            "IDEN was set to 0 but density output was requested.")
               CASE(1) ! sigmaT density
                  DO NH=1,NP
                     WRITE(44,REC=I3DGDRec+1) NH
                     I3DGDRec=I3DGDRec+1
                     DO k=1,NFEN
                        WRITE(44,REC=I3DGDRec+1) SigT(NH,k)
                        I3DGDRec=I3DGDRec+1
                     ENDDO
                  END DO
               CASE(2) ! salinity density
                  DO NH=1,NP
                     WRITE(44,REC=I3DGDRec+1) NH
                     I3DGDRec=I3DGDRec+1
                     DO k=1,NFEN
                        WRITE(44,REC=I3DGDRec+1) SigT(NH,k)
                        WRITE(44,REC=I3DGDRec+2) Sal(NH,k)
                        I3DGDRec=I3DGDRec+2
                     END DO
                  END DO
               CASE(3) ! temperature density
                  DO NH=1,NP
                     WRITE(44,REC=I3DGDRec+1) NH
                     I3DGDRec=I3DGDRec+1
                     DO k=1,NFEN
                        WRITE(44,REC=I3DGDRec+1) SigT(NH,k)
                        WRITE(44,REC=I3DGDRec+2) Temp(NH,k)
                        I3DGDRec=I3DGDRec+2
                     ENDDO
                  END DO
               CASE(4) ! salinity+temperature density
                  DO NH=1,NP
                     WRITE(44,REC=I3DGDRec+1) NH
                     I3DGDRec=I3DGDRec+1
                     DO k=1,NFEN
                        WRITE(44,REC=I3DGDRec+1) SigT(NH,k)
                        WRITE(44,REC=I3DGDRec+2) Temp(NH,k)
                        WRITE(44,REC=I3DGDRec+3) Sal(NH,k)
                        I3DGDRec=I3DGDRec+3
                     ENDDO
                  END DO
               CASE DEFAULT
                  CALL allMessage(WARNING,
     &              "IDEN is not +/- 0 to 4 and is invalid.")
               END SELECT
               CLOSE(44)
#ifdef ADCNETCDF
            CASE(3,5) ! netcdf
               IF (myProc.eq.0) THEN
                  CALL writeOutArrayNetCDF(44, TimeLoc, it,
     &               SigTDescript,SalDescript,TempDescript)
               ENDIF
#endif
            CASE DEFAULT
               CALL allMessage(ERROR,
     &            "ABS(I3DGD) must be either 0, 1, 2, 3 or 5.")
            END SELECT
            N3DGD=0
            ! deallocate memory used in globalio
            IF ((MNPROC.gt.1).and.(WRITE_LOCAL_FILES.eqv..false.).AND.
     &          (myProc.eq.0)) THEN
               DEALLOCATE(SigT_g)
               SELECT CASE(ABS(IDEN))
               CASE(0) ! barotropic
                  CALL allMessage(WARNING,
     &            "IDEN was set to 0 but density output was requested.")
               CASE(1)
                  ! already deallocated SigT above
               CASE(2)
                  DEALLOCATE(Sal_g)
               CASE(3)
                  DEALLOCATE(Temp_g)
               CASE(4)
                  DEALLOCATE(Sal_g,Temp_g)
               CASE DEFAULT
                  CALL allMessage(WARNING,
     &            "IDEN is not +/- 0 to 4 and is invalid.")
               END SELECT
            ENDIF
         ENDIF
      ENDIF

C
C.... 3D Velocity Global Output (Unit 45)
C
C kmd48.33bc changed to NE from GT
      IF(I3DGV.NE.0) THEN
         IF ((IT.GT.NTO3DGVS).AND.(IT.LE.NTO3DGVF)) N3DGV=N3DGV+1
C        Check to see if it is time to generate 3D fulldomain velocity output.
         IF (N3DGV.EQ.NSpo3DGV) THEN
            rp = real(q)
            ip = aimag(q)
            IF ((MNPROC.GT.1).and.(WRITE_LOCAL_FILES.eqv..false.)) THEN ! parallel
               IF (myProc.eq.0) THEN
                  ALLOCATE(rp_g(NP_G,NFEN),ip_g(NP_G,NFEN),
     &               WZ_g(NP_G,NFEN))
                  RealQDescript % array2D_g => rp_g
                  ImaginaryQDescript % array2D_g => ip_g
                  WZDescript % array2D_g => WZ_g
               ENDIF
               CALL collectFullDomainArray(RealQDescript,
     &            packNPbyM, unpackNPbyM)
               CALL collectFullDomainArray(ImaginaryQDescript,
     &            packNPbyM, unpackNPbyM)
               CALL collectFullDomainArray(WZDescript,
     &            packNPbyM, unpackNPbyM)
            ENDIF
            SELECT CASE(ABS(I3DGV))
            CASE(1) ! ascii
               IF ((MNPROC.eq.1).or.(WRITE_LOCAL_FILES.eqv..true.)) THEN
C                 Write time header into file
                  OPEN(45,FILE=TRIM(LOCALDIR)//'/'//'fort.45',
     &                 ACCESS='SEQUENTIAL',POSITION='APPEND')
                  WRITE(45,1100) TimeLoc,IT,(Sigma(k),Sigma(k),Sigma(k),
     &                 k=1,NFEN-1),Sigma(NFEN),Sigma(NFEN)
                  I3DGVRec=I3DGVRec+1
                  DO NH=1,NP
                     WRITE(45,1104) NH,(REAL(q(NH,k)),AIMAG(q(NH,k)),
     &                  WZ(NH,k),k=1,NFEN) !ASCII
                     I3DGVRec=I3DGVRec+1
                  END DO
                  CLOSE(45)
               ENDIF
               IF ((MNPROC.GT.1).and.(WRITE_LOCAL_FILES.eqv..false.)) THEN
                  IF (myProc.eq.0) THEN
                     OPEN(45,FILE=TRIM(GLOBALDIR)//'/'//'fort.45',
     &                  ACCESS='SEQUENTIAL',POSITION='APPEND')
                     WRITE(45,1100) TimeLoc,IT,(Sigma(k),Sigma(k),Sigma(k),
     &                  k=1,NFEN-1),Sigma(NFEN),Sigma(NFEN)
                     I3DGVRec=I3DGVRec+1
                     DO NH=1,NP_G
                        WRITE(45,1104) NH,
     &                     (rp_g(NH,k),ip_g(NH,k),WZ_g(NH,k),k=1,NFEN)
                        I3DGVRec=I3DGVRec+1
                     ENDDO
                     CLOSE(45)
                  ENDIF
               ENDIF
            CASE(2) !BINARY FORMAT
               OPEN(45,FILE=TRIM(LOCALDIR)//'/'//'fort.45',
     &              ACCESS='DIRECT',RECL=NByte)
               WRITE(45,REC=I3DGVRec+1) TimeLoc
               WRITE(45,REC=I3DGVRec+2) IT
               I3DGVRec=I3DGVRec+2
               DO NH=1,NP
                  WRITE(45,REC=I3DGVRec+1) NH
                  I3DGVRec=I3DGVRec+1
                  DO k=1,NFEN
                     WRITE(45,REC=I3DGVRec+1) REAL(q(NH,k))
                     WRITE(45,REC=I3DGVRec+2) AIMAG(q(NH,k))
                     WRITE(45,REC=I3DGVRec+3) WZ(NH,k)
                     I3DGVRec=I3DGVRec+3
                  END DO
               END DO
               CLOSE(45)
#ifdef ADCNETCDF
            CASE(3,5) ! netcdf
               IF (myProc.eq.0) THEN
                  CALL writeOutArrayNetCDF(45, TimeLoc, it,
     &               RealQDescript,ImaginaryQDescript,WZDescript)
               ENDIF
#endif
            CASE DEFAULT
               write(scratchMessage,'("Invalid I3DGV: ",I2)') ABS(I3DGV)
               CALL allMessage(ERROR, scratchMessage)
            END SELECT
            N3DGV=0
            IF ((MNPROC.GT.1).and.(WRITE_LOCAL_FILES.eqv..false.).and.
     &          (myProc.eq.0)) THEN
               DEALLOCATE(rp_g,ip_g,WZ_g)
            ENDIF
         ENDIF
      ENDIF

C
C.... 3D Turbulence Global Output (Unit 46)
C
C kmd48.33bc changed to NE from GT
      IF(I3DGT.NE.0) THEN
C        Check to see if it is time to generate 3D fulldomain turbulence output.
         IF((IT.GT.NTO3DGTS).AND.(IT.LE.NTO3DGTF)) N3DGT=N3DGT+1
         IF (N3DGT.EQ.NSpo3DGT) THEN
            IF ((MNPROC.GT.1).and.(WRITE_LOCAL_FILES.eqv..false.)) THEN
               IF (myProc.eq.0) THEN
                  ALLOCATE(q20_g(NP_G,NFEN),l_g(NP_G,NFEN),
     &               EV_g(NP_G,NFEN))
                  Q20Descript % array2D_g => q20_g
                  LDescript % array2D_g => l_g
                  EVDescript % array2D_g => EV_g
               ENDIF
               CALL collectFullDomainArray(Q20Descript,
     &            packNPbyM, unpackNPbyM)
               CALL collectFullDomainArray(LDescript,
     &            packNPbyM, unpackNPbyM)
               CALL collectFullDomainArray(EVDescript,
     &            packNPbyM, unpackNPbyM)
            ENDIF
            SELECT CASE(ABS(I3DGT))
            CASE(1) ! ascii
               IF((MNPROC.EQ.1).OR.(WRITE_LOCAL_FILES.eqv..true.)) THEN
                  OPEN(46,FILE=TRIM(LOCALDIR)//'/'//'fort.46',
     &              ACCESS='SEQUENTIAL',POSITION='APPEND')
                  WRITE(46,1100) TimeLoc,IT,(Sigma(k),Sigma(k),Sigma(k),
     &                       k=1,NFEN-1),Sigma(NFEN),Sigma(NFEN)
                  I3DGTRec=I3DGTRec+1
                  DO NH=1,NP
                     WRITE(46,1104) NH,
     &               (q20(NH,k),l(NH,k),EV(NH,k),k=1,NFEN)
                     I3DGTRec=I3DGTRec+1
                  ENDDO
                  CLOSE(46)
               ENDIF
               IF ((MNPROC.GT.1).and.(WRITE_LOCAL_FILES.eqv..false.)) THEN
                  IF (myProc.eq.0) THEN
                     OPEN(46,FILE=TRIM(GLOBALDIR)//'/'//'fort.46',
     &                  ACCESS='SEQUENTIAL',POSITION='APPEND')
                     WRITE(46,1100) TimeLoc,IT,(Sigma(k),Sigma(k),Sigma(k),
     &                  k=1,NFEN-1),Sigma(NFEN),Sigma(NFEN)
                     I3DGTRec=I3DGTRec+1
                     DO NH=1,NP_G
                        WRITE(46,1104) NH,
     &                     (q20_g(NH,k),l_g(NH,k),EV_g(NH,k),k=1,NFEN)
                        I3DGTRec=I3DGTRec+1
                     ENDDO
                     CLOSE(46)
                  ENDIF
               ENDIF
            CASE(2)
               OPEN(46,FILE=TRIM(LOCALDIR)//'/'//'fort.46',
     &            ACCESS='DIRECT',RECL=NByte)
               WRITE(46,REC=I3DGTRec+1) TimeLoc
               WRITE(46,REC=I3DGTRec+2) IT
               I3DGTRec=I3DGTRec+2
               DO NH=1,NP
                  WRITE(46,REC=I3DGTRec+1) NH
                  I3DGTRec=I3DGTRec+1
                  DO k=1,NFEN
                     WRITE(46,REC=I3DGTRec+1) q20(NH,k)
                     WRITE(46,REC=I3DGTRec+2) l(NH,k)
                     WRITE(46,REC=I3DGTRec+3) EV(NH,k)
                     I3DGTREC=I3DGTREC+3
                  ENDDO
               ENDDO
               CLOSE(46)
#ifdef ADCNETCDF
            CASE(3,5) ! netcdf
               IF (myProc.eq.0) THEN
                  CALL writeOutArrayNetCDF(46, TimeLoc, it,
     &               Q20Descript,LDescript,EVDescript)
               ENDIF
#endif
            CASE DEFAULT
               write(scratchMessage,'("Invalid I3DGT: ",I2)') ABS(I3DGT)
               CALL allMessage(ERROR, scratchMessage)
            END SELECT
            N3DGT=0
            IF ((MNPROC.GT.1).and.(WRITE_LOCAL_FILES.eqv..false.).and.
     &          (myProc.eq.0)) THEN
               DEALLOCATE(q20_g,l_g,EV_g)
            ENDIF

         ENDIF
      ENDIF

#ifdef WRITE_OUTPUT_TRACE
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
      END

C-----------------------------------------------------------------------
C     S U B R O U T I N E   W R I T E  O U T  A R R A Y
C-----------------------------------------------------------------------
C     jgf48.03 This subroutine was created to write out a column
C     vector (i.e., nodal data such as water surface elevation or
C     pressure) to a file.
C-----------------------------------------------------------------------
Casey 090302: Added the filename as the last argument.
      SUBROUTINE writeOutArray(lun, TimeLoc, it, descript, pack_cmd,
     &                         unpack_cmd, filepos, fn)

      USE SIZES
      USE GLOBAL
      USE GLOBAL_IO, ONLY : collectFullDomainArray
#ifdef ADCNETCDF
      USE NETCDFIO, ONLY : writeOutArrayNetCDF
#endif
      IMPLICIT NONE
C     args
      INTEGER, intent(in) :: lun ! logical unit number of file to write to
      REAL(8), intent(in) :: TimeLoc ! seconds since cold start
      INTEGER, intent(in) :: it   ! number of time steps since cold start
      type(OutputDataDescript_t), intent(in) :: descript !describes output data
      EXTERNAL :: pack_cmd   ! subroutine used to pack data on subdomain
      EXTERNAL :: unpack_cmd ! subroutine used to unpack data on proc 0
      INTEGER, intent(inout) :: filepos  ! current position in the output file
C     local vars
Casey 090302: Increase the length of the filename from 7 to 15.
      CHARACTER(15) :: fn     ! outfile name (valid for lun between 10 and 99)
      INTEGER :: I           ! loop counter

C     initialize output file name
Casey 090302: The filename is now passed as an input argument.
C     fn(1:5) = 'fort.'
C     WRITE(fn(6:7),2) lun

C     collect up the data from subdomains if running in parallel
      IF ((MNPROC.gt.1).and.(WRITE_LOCAL_FILES.eqv..false.)) THEN
C         write(16,*) 'About to collectFullDomainArray' !jgfdebug48.03
         CALL collectFullDomainArray(descript, pack_cmd, unpack_cmd)
      ENDIF

C     write data according to format specifier from fort.15 (e.g., NOUTE)
      SELECT CASE (ABS(descript % specifier))

      CASE(1) !ascii text

C         write(16,*) 'About to open globalio text file.' !jgfdebug48.03
         IF ( (MNPROC.gt.1).and.(MyProc.eq.0)
     &        .and.(.not.WRITE_LOCAL_FILES)) THEN
            OPEN(lun,FILE=TRIM(GLOBALDIR)//'/'//fn,
     &           ACCESS='SEQUENTIAL',POSITION='APPEND')
            WRITE(lun,2120) TimeLoc,IT
            IF (descript % num_items_per_record .eq. 1) THEN
               DO I=1, descript % num_fd_records
C                  WRITE(16,*) I ! jgfdebug48.03
                  WRITE(lun,2453) I, descript % array_g(I)
               ENDDO
            ENDIF
            IF (descript % num_items_per_record .eq. 2) THEN
               DO I=1, descript % num_fd_records
C                  WRITE(16,*) I ! jgfdebug48.03
                  WRITE(lun,2454) I, descript % array_g(I),
     &                               descript % array2_g(I)
               ENDDO
            ENDIF
C            WRITE(16,*) 'About to close globalio text file.'!jgfdebug48.03
            CLOSE(lun)
         ENDIF

         IF ((MNPROC.eq.1).or.(WRITE_LOCAL_FILES)) THEN
            OPEN(lun,FILE=TRIM(LOCALDIR)//'/'//fn,
     &           ACCESS='SEQUENTIAL',POSITION='APPEND')
            WRITE(lun,2120) TimeLoc,IT
            IF (descript % num_items_per_record .eq. 1) THEN
!......TCM - v48.4618 -- Fixed dry node output for serial run case
               IF ((trim(descript % field_name) .eq. 'Elev').and.
     &             (descript % ConsiderWetDry .EQV. .TRUE.)) THEN
                  DO I=1, descript % num_records_this
                     if(NODECODE(I).EQ.1) THEN
                        WRITE(lun,2453) I, descript % array(I)
                     ELSE
                        WRITE(lun,2453) I, descript % alternate_value  !-99999.0 for dry nodes
                     ENDIF
                  END DO
               ELSE
                  DO I=1, descript % num_records_this
                     WRITE(lun,2453) I, descript % array(I)
                  END DO
               ENDIF
            ENDIF
            IF (descript % num_items_per_record .eq. 2) THEN
               DO I=1, descript % num_records_this
                  WRITE(lun,2454) I, descript % array(I),
     &                               descript % array2(I)
               END DO
            ENDIF
            CLOSE(lun)
         ENDIF
         filepos = filepos+1+descript % num_records_this

      CASE(2) !binary (nonportable)

         IF ( (MNPROC.gt.1).and.(MyProc.eq.0)
     &         .and.(.not.WRITE_LOCAL_FILES)) THEN
            OPEN(lun,FILE=TRIM(GLOBALDIR)//'/'//fn,
     &           ACCESS='DIRECT',RECL=NBYTE)
            WRITE(lun,REC=filepos+1) TimeLoc
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
            WRITE(lun,REC=filepos+1) TimeLoc
            WRITE(lun,REC=filepos+2) IT
            filepos = filepos + 2
            IF ( descript % num_items_per_record .eq. 1 ) THEN
               IF ((trim(descript % field_name) .eq. 'Elev').and.
     &             (descript % ConsiderWetDry .EQV. .TRUE.)) THEN
                  DO I=1, descript % num_records_this
                     if(NODECODE(I).EQ.1) THEN
                        WRITE(lun,REC=filepos+I) descript % array(I)
                     ELSE
                        WRITE(lun,REC=filepos+I)
     &                    descript % alternate_value !-99999.0 for dry nodes
                     ENDIF
                  END DO
               ELSE
                  DO I=1, descript % num_records_this
                     WRITE(lun,REC=filepos+I) descript % array(I)
                  END DO
               ENDIF
            ENDIF
            IF ( descript % num_items_per_record .eq. 2 ) THEN
               DO I=1, descript % num_records_this
                  !tcmv48.4618 -- changed from array_g to array
                  WRITE(lun,REC=filepos+2*I-1) descript % array(I)
                  !tcmv48.4618 -- changed from array2_g to array2
                  WRITE(lun,REC=filepos+2*I)   descript % array2(I)
               END DO
            ENDIF
            CLOSE(lun)
         ENDIF
         filepos = filepos + descript % num_records_this

      CASE(3,5) !netcdf (portable)
#ifdef ADCNETCDF
         IF (MYPROC.EQ.0) THEN
            CALL writeOutArrayNetCDF(lun, TimeLoc, it, descript)
         ENDIF
#else
         WRITE(ScreenUnit,*) 'ERROR: NetCDF is not available.'
         WRITE(16,*) 'ERROR: NetCDF is not available.'
#endif

      CASE DEFAULT
         WRITE(ScreenUnit,*) 'ERROR: Invalid output specifier.'
      END SELECT

 2    FORMAT(I2)
 2120 FORMAT(2X,1pE20.10E3,5X,I10)
 2453 FORMAT(2x, i8, 2x, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3)
 2454 FORMAT(2X,I8,2(2X,1pE20.10E3))
C-----------------------------------------------------------------------
       END SUBROUTINE writeOutArray
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   W R I T E  O U T  A R R A Y  M I N  M A X
C-----------------------------------------------------------------------
C     jgf48.4636 This subroutine was created to write out a column
C     vector (i.e., nodal data such as water surface elevation or
C     pressure) to a file, for quantities whose min or max values
C     have been stored over the course of the run.
C-----------------------------------------------------------------------
      SUBROUTINE writeOutArrayMinMax(lun, TimeLoc, it, descript, pack_cmd,
     &                         unpack_cmd)

      USE SIZES
      USE GLOBAL
      USE GLOBAL_IO, ONLY : collectFullDomainArray, HEADER_MAX,
     &                      open_minmax_file
#ifdef ADCNETCDF
      USE NetCDFIO, ONLY :  writeOutArrayNetCDF
#endif
      IMPLICIT NONE

C     args
      INTEGER, intent(in) :: lun ! logical unit number of file to write to
      REAL(8), intent(in) :: TimeLoc ! seconds since cold start
      INTEGER, intent(in) :: it   ! number of time steps since cold start
      type(OutputDataDescript_t), intent(in) :: descript !describes output data
      EXTERNAL :: pack_cmd   ! subroutine used to pack data on subdomain
      EXTERNAL :: unpack_cmd ! subroutine used to unpack data on proc 0
      CHARACTER(len=3) :: cmp ! comparison to make vs hotstart val ('min' or 'max')
      INTEGER :: I          ! loop counter

C     collect up the data from subdomains if running in parallel
      IF ((MNPROC.gt.1).and.(.not.WRITE_LOCAL_FILES)) THEN
         CALL collectFullDomainArray(descript, pack_cmd, unpack_cmd)
      ENDIF

C     write data in ascii text only (other output formats are not supported)
      SELECT CASE (ABS(descript % specifier))

      CASE(0,1,2,4) !ascii text
         CALL OPEN_MINMAX_FILE(lun, TRIM(GLOBALDIR)//'/'
     &           // descript%file_name, NP_G, NP, HEADER_MAX)
         IF ((MNPROC.gt.1).and.(MyProc.eq.0).and.
     &        (.not.WRITE_LOCAL_FILES)) THEN
            OPEN(lun,
     &           FILE=TRIM(GLOBALDIR)//'/'//TRIM(descript%file_name),
     &           ACCESS='SEQUENTIAL',POSITION='APPEND')
            WRITE(lun,2120) TimeLoc,IT
            IF ((IHOT.EQ.0).OR.(IHOT.EQ.17)) THEN
               DO I=1, descript % num_fd_records
                  WRITE(lun,2453) I, descript % array_g(I)
               ENDDO
            ELSE IF ( descript%file_name.eq."minpr.63" ) THEN
                   DO I=1, descript % num_fd_records
                     WRITE(lun,2453) I, MIN(descript % array_g(I),
     &                   descript % hotstart_g(I))
                  ENDDO
            ELSE
               DO I=1, descript % num_fd_records
                  WRITE(lun,2453) I, MAX(descript % array_g(I),
     &                descript % hotstart_g(I))
               ENDDO
            ENDIF
            CLOSE(lun)
         ENDIF

         IF ((MNPROC.eq.1).or.(WRITE_LOCAL_FILES)) THEN
            CALL OPEN_MINMAX_FILE(lun, TRIM(LOCALDIR)//'/'
     &           // descript%file_name, NP_G, NP, HEADER_MAX)
            OPEN(lun,FILE=TRIM(LOCALDIR)//'/'//TRIM(descript%file_name),
     &           ACCESS='SEQUENTIAL',POSITION='APPEND')
            WRITE(lun,2120) TimeLoc,IT
            IF ((IHOT.EQ.0).OR.(IHOT.EQ.17)) THEN
               DO I=1, descript % num_records_this
                  WRITE(lun,2453) I, descript % array(I)
               ENDDO
            ELSE IF ( descript%file_name.eq."minpr.63" ) THEN
               DO I=1, descript % num_records_this
                  WRITE(lun,2453) I, MIN(descript % array(I),
     &                descript % hotstart(I))
               ENDDO
            ELSE
               DO I=1, descript % num_records_this
                  WRITE(lun,2453) I, MAX(descript % array(I),
     &                descript % hotstart(I))
               ENDDO
            ENDIF
         ENDIF
Casey 120830: Reversed the next two lines.
      CASE(3,5) !netcdf (portable)
#ifdef ADCNETCDF
         IF (MYPROC.EQ.0) THEN
             CALL writeOutArrayNetCDF(lun, TimeLoc, it, descript)
         ENDIF
#else
         WRITE(ScreenUnit,*) 'ERROR: NetCDF is not available.'
         WRITE(16,*) 'ERROR: NetCDF is not available.'
#endif
      CASE DEFAULT
         WRITE(ScreenUnit,*) 'ERROR: Invalid output specifier.'
      END SELECT

 2    FORMAT(I2)
 2120 FORMAT(2X,1pE20.10E3,5X,I10)
 2453 FORMAT(2x, i8, 2x, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3)
 2454 FORMAT(2X,I8,2(2X,1pE20.10E3))
C-----------------------------------------------------------------------
       END SUBROUTINE writeOutArrayMinMax
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   W R I T E  H O T S T A R T
C-----------------------------------------------------------------------
C     jgf48.03 This subroutine was created from code in timestep.F to
C     write hotstart files.
C-----------------------------------------------------------------------
      SUBROUTINE writeHotstart(TimeLoc, IT)

      USE SIZES
      USE GLOBAL
      USE HARM
      USE GLOBAL_IO, ONLY: packOne, unpackOne, packOneInt, unpackOneInt,
     &    packTwo, unpackTwo, packMbyNP, unpackMbyNP, HEADER_MAX,
     &    collectFullDomainArray, collectFullDomainIntArray
      USE ADCIRC_MOD, ONLY : ADCIRC_Terminate
#ifdef CMPI
      USE MESSENGER
#endif
#ifdef ADCNETCDF
       USE NETCDFIO, ONLY : initNetCDFHotstart,
     &                      writeNetCDFHotstart,
     &                      initNetCDFHotstartHarmonic,
     &                      initNetCDFHotstartHarmonicMeansVariances,
     &                      writeNetCDFHotstartHarmonic,
     &                      writeNetCDFHotstartHarmonicMeansVariances
       USE NodalAttributes, ONLY :
     &     nolibf, nwp, tau0, cf, eslm
#endif

      IMPLICIT NONE
      REAL(8), intent(in) :: TimeLoc
      INTEGER, intent(in) :: IT
      INTEGER I,J,K
      type(OutputDataDescript_t) :: descript
      LOGICAL TEST_HOTSTART     !jgf45.07 used for testing hot start capability
      INTEGER npx, nex
C
      type(OutputDataDescript_t), SAVE :: Elev1Descript
      type(OutputDataDescript_t), SAVE :: Elev2Descript
      type(OutputDataDescript_t), SAVE :: VelDescript
      type(OutputDataDescript_t), SAVE :: CH1Descript
      type(OutputDataDescript_t), SAVE :: EtaDiscDescript
      type(OutputDataDescript_t), SAVE :: NodeCodeDescript
      type(OutputDataDescript_t), SAVE :: NOFFDescript
      LOGICAL ncerror
      INTEGER numHotstartWrites ! number writes to hot start files
      INTEGER nextLun           ! next LUN to write to, after initial write
C
C     !jgf49.44: harmonic analysis vars
      INTEGER N ! number of stations
      INTEGER NSTAEX  ! num of elev stations in either fulldomain or subdomain
      INTEGER NSTAVX  ! num of vel stations in either fulldomain or subdomain
      type(OutputDataDescript_t), SAVE :: HarmElevFDLVDescript
      type(OutputDataDescript_t), SAVE :: HarmElevSLVDescript
      type(OutputDataDescript_t), SAVE :: HarmUVelFDLVDescript
      type(OutputDataDescript_t), SAVE :: HarmVVelFDLVDescript
      type(OutputDataDescript_t), SAVE :: HarmUVelSLVDescript
      type(OutputDataDescript_t), SAVE :: HarmVVelSLVDescript
      type(OutputDataDescript_t), SAVE :: ELAVDescript
      type(OutputDataDescript_t), SAVE :: ELVADescript
      type(OutputDataDescript_t), SAVE :: XVELAVDescript
      type(OutputDataDescript_t), SAVE :: YVELAVDescript
      type(OutputDataDescript_t), SAVE :: XVELVADescript
      type(OutputDataDescript_t), SAVE :: YVELVADescript
C
      LOGICAL, SAVE :: FirstCall = .true.
C
      call setMessageSource("writeHotstart")
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ncerror = .false.
C
      IF (FirstCall) THEN

         Elev1Descript % specifier            = NHSTAR
         Elev1Descript % initial_value        = 0.0
         Elev1Descript % num_items_per_record = 1
         Elev1Descript % num_fd_records       = NP_G
         Elev1Descript % num_records_this     = NP
         Elev1Descript % imap                 => NODES_LG
         Elev1Descript % array                => ETA1
         Elev1Descript % array_g              => ETA1_g

         Elev2Descript % specifier            = NHSTAR
         Elev2Descript % initial_value        = 0.0
         Elev2Descript % num_items_per_record = 1
         Elev2Descript % num_fd_records       = NP_G
         Elev2Descript % num_records_this     = NP
         Elev2Descript % imap                 => NODES_LG
         Elev2Descript % array                => ETA2
         Elev2Descript % array_g              => ETA2_g

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
            CH1Descript % specifier            = NHSTAR
            CH1Descript % initial_value        = 0.0
            CH1Descript % num_items_per_record = 1
            CH1Descript % num_fd_records       = NP_G
            CH1Descript % num_records_this     = NP
            CH1Descript % imap                 => NODES_LG
            CH1Descript % array                => CH1
            CH1Descript % array_g              => CH1_g
         ENDIF

         EtaDiscDescript % specifier            = NHSTAR
         EtaDiscDescript % initial_value        = 0.0
         EtaDiscDescript % num_items_per_record = 1
         EtaDiscDescript % num_fd_records       = NP_G
         EtaDiscDescript % num_records_this     = NP
         EtaDiscDescript % imap                 => NODES_LG
         EtaDiscDescript % array                => EtaDisc
         EtaDiscDescript % array_g              => EtaDisc_g

         NodeCodeDescript % specifier            = NHSTAR
         NodeCodeDescript % int_initial_value        = 0
         NodeCodeDescript % num_items_per_record = 1
         NodeCodeDescript % num_fd_records       = NP_G
         NodeCodeDescript % num_records_this     = NP
         NodeCodeDescript % imap                 => NODES_LG
         NodeCodeDescript % iarray               => NODECODE
         NodeCodeDescript % iarray_g             => NODECODE_g

         NOFFDescript % specifier            = NHSTAR
         NOFFDescript % int_initial_value        = 0
         NOFFDescript % num_items_per_record = 1
         NOFFDescript % num_fd_records       = NE_G
         NOFFDescript % num_records_this     = NE
         NOFFDescript % imap                 => IMAP_EL_LG
         NOFFDescript % iarray               => NOFF
         NOFFDescript % iarray_g             => NOFF_g
C
         IF (IHARIND.eq.1) THEN
            HarmElevFDLVDescript % specifier            = NHSTAR
            HarmElevFDLVDescript % initial_value        = 0.0
            HarmElevFDLVDescript % num_items_per_record = 2*MNHARF
            HarmElevFDLVDescript % num_fd_records       = NP_G
            HarmElevFDLVDescript % num_records_this     = NP
            HarmElevFDLVDescript % imap                 => NODES_LG
            HarmElevFDLVDescript % array2D              => GLOELV
            HarmElevFDLVDescript % array2D_g            => GLOELV_g

            HarmElevSLVDescript % specifier            = NHSTAR
            HarmElevSLVDescript % initial_value        = 0.0
            HarmElevSLVDescript % num_items_per_record = 2*MNHARF
            HarmElevSLVDescript % num_fd_records       = NSTAE_G
            HarmElevSLVDescript % num_records_this     = NSTAE
            HarmElevSLVDescript % imap                 => IMAP_STAE_LG
            HarmElevSLVDescript % array2D              => STAELV
            HarmElevSLVDescript % array2D_g            => STAELV_g

            HarmUVelFDLVDescript % specifier            = NHSTAR
            HarmUVelFDLVDescript % initial_value        = 0.0
            HarmUVelFDLVDescript % num_items_per_record = 2*MNHARF
            HarmUVelFDLVDescript % num_fd_records       = NP_G
            HarmUVelFDLVDescript % num_records_this     = NP
            HarmUVelFDLVDescript % imap                 => NODES_LG
            HarmUVelFDLVDescript % array2D              => GLOULV
            HarmUVelFDLVDescript % array2D_g            => GLOULV_g

            HarmVVelFDLVDescript % specifier            = NHSTAR
            HarmVVelFDLVDescript % initial_value        = 0.0
            HarmVVelFDLVDescript % num_items_per_record = 2*MNHARF
            HarmVVelFDLVDescript % num_fd_records       = NP_G
            HarmVVelFDLVDescript % num_records_this     = NP
            HarmVVelFDLVDescript % imap                 => NODES_LG
            HarmVVelFDLVDescript % array2D              => GLOVLV
            HarmVVelFDLVDescript % array2D_g            => GLOVLV_g

            HarmUvelSLVDescript % specifier            = NHSTAR
            HarmUVelSLVDescript % initial_value        = 0.0
            HarmUVelSLVDescript % num_items_per_record = 2*MNHARF
            HarmUVelSLVDescript % num_fd_records       = NSTAV_G
            HarmUVelSLVDescript % num_records_this     = NSTAV
            HarmUVelSLVDescript % imap                 => IMAP_STAV_LG
            HarmUVelSLVDescript % array2D              => STAULV
            HarmUVelSLVDescript % array2D_g            => STAULV_g

            HarmVVelSLVDescript % specifier            = NHSTAR
            HarmVVelSLVDescript % initial_value        = 0.0
            HarmVVelSLVDescript % num_items_per_record = 2*MNHARF
            HarmVVelSLVDescript % num_fd_records       = NSTAV_G
            HarmVVelSLVDescript % num_records_this     = NSTAV
            HarmVVelSLVDescript % imap                 => IMAP_STAV_LG
            HarmVVelSLVDescript % array2D              => STAVLV
            HarmVVelSLVDescript % array2D_g            => STAVLV_g
         ENDIF

         IF (CHARMV.eqv..true.) THEN
            ELAVDescript % specifier            = NHSTAR
            ELAVDescript % initial_value        = 0.0
            ELAVDescript % num_items_per_record = 1
            ELAVDescript % num_fd_records       = NP_G
            ELAVDescript % num_records_this     = NP
            ELAVDescript % imap                 => NODES_LG
            ELAVDescript % array                => ELAV
            ELAVDescript % array_g              => ELAV_g

            ELVADescript % specifier            = NHSTAR
            ELVADescript % initial_value        = 0.0
            ELVADescript % num_items_per_record = 1
            ELVADescript % num_fd_records       = NP_G
            ELVADescript % num_records_this     = NP
            ELVADescript % imap                 => NODES_LG
            ELVADescript % array                => ELVA
            ELVADescript % array_g              => ELVA_g

            XVELAVDescript % specifier            = NHSTAR
            XVELAVDescript % initial_value        = 0.0
            XVELAVDescript % num_items_per_record = 1
            XVELAVDescript % num_fd_records       = NP_G
            XVELAVDescript % num_records_this     = NP
            XVELAVDescript % imap                 => NODES_LG
            XVELAVDescript % array                => XVELAV
            XVELAVDescript % array_g              => XVELAV_g

            YVELAVDescript % specifier            = NHSTAR
            YVELAVDescript % initial_value        = 0.0
            YVELAVDescript % num_items_per_record = 1
            YVELAVDescript % num_fd_records       = NP_G
            YVELAVDescript % num_records_this     = NP
            YVELAVDescript % imap                 => NODES_LG
            YVELAVDescript % array                => YVELAV
            YVELAVDescript % array_g              => YVELAV_g

            XVELVADescript % specifier            = NHSTAR
            XVELVADescript % initial_value        = 0.0
            XVELVADescript % num_items_per_record = 1
            XVELVADescript % num_fd_records       = NP_G
            XVELVADescript % num_records_this     = NP
            XVELVADescript % imap                 => NODES_LG
            XVELVADescript % array                => XVELVA
            XVELVADescript % array_g              => XVELVA_g

            YVELVADescript % specifier            = NHSTAR
            YVELVADescript % initial_value        = 0.0
            YVELVADescript % num_items_per_record = 1
            YVELVADescript % num_fd_records       = NP_G
            YVELVADescript % num_records_this     = NP
            YVELVADescript % imap                 => NODES_LG
            YVELVADescript % array                => YVELVA
            YVELVADescript % array_g              => YVELVA_g
         ENDIF
C
#ifdef ADCNETCDF
         IF ((NHSTAR.eq.3).or.(NHSTAR.eq.367).or.(NHSTAR.eq.368).or.
     &       (NHSTAR.eq.5).or.(NHSTAR.eq.567).or.(NHSTAR.eq.568)) THEN
            ! jgf49.41: In parallel, we don't need to create the hotstart file, it was
            ! created by adcprep.
            IF (MNPROC.eq.1) THEN  ! serial
C              jgf49.35 TODO: Allocate arrays to write 3D data
               CALL initNetCDFHotstart(hss%lun, Elev1Descript,
     &         Elev2Descript, VelDescript, CH1Descript,
     &         EtaDiscDescript, NodeCodeDescript, NOFFDescript, ncerror)
               IF ( ncerror.eqv..true. ) CALL ADCIRC_Terminate()
               IF (IHARIND.eq.1) THEN
                  CALL initNetCDFHotstartHarmonic(hss%lun,
     &            HarmElevFDLVDescript, HarmElevSLVDescript,
     &            HarmUVelFDLVDescript, HarmVVelFDLVDescript,
     &            HarmUVelSLVDescript, HarmVVelSLVDescript, ncerror)
                  IF ( ncerror.eqv..true. ) CALL ADCIRC_Terminate()
                  IF (CHARMV.eqv..true.) THEN
                     CALL initNetCDFHotstartHarmonicMeansVariances(
     &                  hss%lun, ELAVDescript, ELVADescript,
     &                  XVELAVDescript, YVELAVDescript, XVELVADescript,
     &                  YVELVADescript,ncerror)
                     IF ( ncerror.eqv..true. ) CALL ADCIRC_Terminate()
                  ENDIF
               ENDIF
               ! Determine if we need to init both hotstart files, i.e., if the
               ! hotstart timestep increment is short enough that we will need
               ! to write both the 67 and 68 file during this run
               numHotstartWrites = (NT-IT)/NHSINC
               IF (numHotstartWrites.ge.1) THEN
                  IF (hss%lun.eq.67) THEN
                     nextLun = 68
                  ELSE
                     nextLun = 67
                  ENDIF
                  CALL initNetCDFHotstart(nextLun, Elev1Descript,
     &            Elev2Descript, VelDescript, CH1Descript,
     &            EtaDiscDescript, NodeCodeDescript, NOFFDescript,
     &            ncerror)
                  IF ( ncerror.eqv..true. ) CALL ADCIRC_Terminate()
                  IF (IHARIND.eq.1) THEN
                     CALL initNetCDFHotstartHarmonic(nextLun,
     &               HarmElevFDLVDescript, HarmElevSLVDescript,
     &               HarmUVelFDLVDescript, HarmVVelFDLVDescript,
     &               HarmUVelSLVDescript, HarmVVelSLVDescript, ncerror)
                     IF ( ncerror.eqv..true. ) CALL ADCIRC_Terminate()
                     IF (CHARMV.eqv..true.) THEN
                        CALL initNetCDFHotstartHarmonicMeansVariances(
     &                     nextLun, ELAVDescript, ELVADescript,
     &                     XVELAVDescript, YVELAVDescript,
     &                     XVELVADescript, YVELVADescript, ncerror)
                        IF ( ncerror.eqv..true. ) CALL ADCIRC_Terminate()
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
#endif
C
         FirstCall = .false.
      ENDIF

C     collect up the data from subdomains if running in parallel
ckmd48.33bc - changed the NOFF and NodeCode calls
#ifdef CMPI
      IF (WRITE_LOCAL_HOT_START_FILES.eqv..FALSE.) THEN
         CALL collectFullDomainArray(Elev1Descript, packOne, unpackOne)
         CALL collectFullDomainArray(Elev2Descript, packOne, unpackOne)
         CALL collectFullDomainArray(VelDescript, packTwo, unpackTwo)
         IF (IM.eq.10) THEN
            CALL collectFullDomainArray(CH1Descript, packOne, unpackOne)
         ENDIF
         CALL collectFullDomainArray(EtaDiscDescript,
     &             packOne, unpackOne)
         CALL collectFullDomainIntArray(NodeCodeDescript,
     &             packOneInt, unpackOneInt, 1)
         CALL collectFullDomainIntArray(NOFFDescript,
     &             packOneInt, unpackOneInt, 2)
!         TCM v48.4638 20090902 -- NOFF could have the wrong
!         values after globalizing because there are no rules for which
!         processor owns an element and Noff is an elemental value
!         The way around is to take the minimum of (Noff, 1) as
!         NOFF is supposed to be either 0 or 1.
!
!        jgf48.4642 Modified syntax to operate directly on NOFF array
!        instead of the pointer to the NOFF array (i.e.
!        NOFFDescript%iarray) b/c g95 choked on the original syntax.
!         WHERE (NOFF.gt.1) NOFF = 1

!         TCM49.14 20100604 -- Neither of the previous fixes for NOFF
!             worked.  The correct way is to fix NOFF_G
!             the Global list which is owned only by MyProc = 0
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            WHERE (NOFF_G.gt.1) NOFF_G = 1
         endif
C
C        Collect harmonic analysis data if needed.
         IF ((IHARIND.eq.1).and.(IT.gt.ITHAS)) THEN
            IF (NHAGE.ne.0) THEN
               CALL collectFullDomainArray(HarmElevFDLVDescript,
     &              packMbyNP, unpackMbyNP)
            ENDIF
            IF (NHASE.ne.0) THEN
               CALL collectFullDomainArray(HarmElevSLVDescript,
     &              packMbyNP, unpackMbyNP)
            ENDIF
            IF (NHAGV.ne.0) THEN
               CALL collectFullDomainArray(HarmUVelFDLVDescript,
     &              packMbyNP, unpackMbyNP)
               CALL collectFullDomainArray(HarmVVelFDLVDescript,
     &              packMbyNP, unpackMbyNP)
            ENDIF
            IF (NHASV.ne.0) THEN
               CALL collectFullDomainArray(HarmUVelSLVDescript,
     &              packMbyNP, unpackMbyNP)
               CALL collectFullDomainArray(HarmVVelSLVDescript,
     &              packMbyNP, unpackMbyNP)
            ENDIF
         ENDIF
C
C        Collect timeseries reconstruction data if needed.
         IF (CHARMV) THEN
            IF ((IHARIND.eq.1).and.(IT.gt.ITMV)) THEN
               CALL collectFullDomainArray(ELAVDescript,
     &              packOne, unpackOne)
               CALL collectFullDomainArray(ELVADescript,
     &              packOne, unpackOne)
               CALL collectFullDomainArray(XVELAVDescript,
     &              packOne, unpackOne)
               CALL collectFullDomainArray(YVELAVDescript,
     &              packOne, unpackOne)
               CALL collectFullDomainArray(XVELVADescript,
     &              packOne, unpackOne)
               CALL collectFullDomainArray(YVELVADescript,
     &              packOne, unpackOne)
            ENDIF
         ENDIF

      ENDIF
#endif
C
      SELECT CASE (NHSTAR)
C
      CASE(1,67,68) ! 1=nonportable binary
C        jgf49.35 Values of 67 and 68 are used for testing and cause the
C        program to exit immediately after writing the hotstart file.
         NPX = NP_G
         NEX = NE_G
         if ((MNPROC.eq.1) .or.
     &      (WRITE_LOCAL_HOT_START_FILES.eqv..TRUE.)) THEN
            NPX = NP
            NEX = NE
         endif
         IF ((MYPROC.eq.0).or.
     &         (WRITE_LOCAL_HOT_START_FILES.eqv..TRUE.)) THEN

            OPEN(hss % lun ,FILE=TRIM(HOTSTARTDIR)//'/'//hss % filename,
     $      ACCESS='DIRECT',STATUS='REPLACE',RECL=8)
            IHOTSTP=1
            WRITE(hss % lun,REC=IHOTSTP) FileFmtVersion
            IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) IM   ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) TimeLoc ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) IT   ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) NPX  ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) NEX  ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) NPX  ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) NEX  ; IHOTSTP = IHOTSTP + 1
         ENDIF
         IF ((MNPROC.gt.1).and.(MYPROC.eq.0).and.
     &         (WRITE_LOCAL_HOT_START_FILES.eqv..FALSE.)) THEN
            DO I=1, Elev1Descript % num_fd_records
               WRITE(hss % lun,REC=IHOTSTP) Elev1Descript % array_g(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1, Elev2Descript % num_fd_records
               WRITE(hss % lun,REC=IHOTSTP) Elev2Descript % array_g(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1, EtaDiscDescript % num_fd_records
               WRITE(hss % lun,REC=IHOTSTP) EtaDiscDescript % array_g(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1, VelDescript % num_fd_records
               WRITE(hss % lun,REC=IHOTSTP) VelDescript % array_g(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1, VelDescript % num_fd_records
               WRITE(hss % lun,REC=IHOTSTP) VelDescript % array2_g(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            IF (IM.eq.10) THEN
               DO I=1, CH1Descript % num_fd_records
                  WRITE(hss % lun,REC=IHOTSTP) CH1Descript % array_g(I)
                  IHOTSTP=IHOTSTP+1
               ENDDO
            ENDIF
            DO I=1, NodeCodeDescript % num_fd_records
              WRITE(hss % lun,REC=IHOTSTP)NodeCodeDescript % iarray_g(I)
              IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1, NOFFDescript % num_fd_records
               WRITE(hss % lun,REC=IHOTSTP) NOFFDescript % iarray_g(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
         ENDIF
         IF ((MNPROC.eq.1).or.
     &      (WRITE_LOCAL_HOT_START_FILES.eqv..TRUE.)) THEN
            DO I=1, Elev1Descript % num_records_this
               WRITE(hss % lun,REC=IHOTSTP) Elev1Descript % array(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1, Elev2Descript % num_records_this
               WRITE(hss % lun,REC=IHOTSTP) Elev2Descript % array(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1, EtaDiscDescript % num_records_this
               WRITE(hss % lun,REC=IHOTSTP) EtaDiscDescript % array(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1, VelDescript % num_records_this
               WRITE(hss % lun,REC=IHOTSTP) VelDescript % array(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
           DO I=1, VelDescript % num_records_this
               WRITE(hss % lun,REC=IHOTSTP) VelDescript % array2(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            IF (IM.eq.10) THEN
               DO I=1, CH1Descript % num_records_this
                  WRITE(hss % lun,REC=IHOTSTP) CH1Descript % array(I)
                  IHOTSTP=IHOTSTP+1
               ENDDO
            ENDIF
            DO I=1, NodeCodeDescript % num_records_this
               WRITE(hss % lun,REC=IHOTSTP) NodeCodeDescript % iarray(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1, NOFFDescript % num_records_this
               WRITE(hss % lun,REC=IHOTSTP) NOFFDescript % iarray(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
         ENDIF
C
         IF ((MYPROC.eq.0).or.
     &            (WRITE_LOCAL_HOT_START_FILES.eqv..TRUE.)) THEN
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
!         if(ncice.ne.0) then
!           WRITE(hss % lun,REC=IHOTSTP) IGIP    ; IHOTSTP = IHOTSTP + 1
!         endif
!         igip is in sync. with igpp so no need to store it
            WRITE(hss % lun,REC=IHOTSTP) NSCOUGW ; IHOTSTP = IHOTSTP + 1
         ENDIF
C...
C...  ADD IN 3D STUFF IF APPROPRIATE
C...
         IF (C3D.eqv..true.) THEN
            CALL logMessage(INFO,"Writing 3D hotstart information.")
            CALL HSTART3D_OUT(IT) ! -> see vsmy.F
         ENDIF
C...
C...  IF APPROPRIATE ADD HARMONIC ANALYSIS INFORMATION TO HOT START FILE
C...
         IF((IHARIND.EQ.1).AND.(IT.GT.ITHAS)) THEN
            IF ((MNPROC.eq.1).or.WRITE_LOCAL_HOT_START_FILES) THEN
               NSTAEX = NSTAE
               NSTAVX = NSTAV
            ELSE
               NSTAEX = NSTAE_G
               NSTAVX = NSTAV_G
            ENDIF
            IF ((MyPROC.eq.0).or.WRITE_LOCAL_HOT_START_FILES) THEN
               WRITE(hss % lun,REC=IHOTSTP) ICHA
               WRITE(hss % lun,REC=IHOTSTP+1) NZ
               WRITE(hss % lun,REC=IHOTSTP+2) NF
               WRITE(hss % lun,REC=IHOTSTP+3) MM
               WRITE(hss % lun,REC=IHOTSTP+4) NPX
               WRITE(hss % lun,REC=IHOTSTP+5) NSTAEX
               WRITE(hss % lun,REC=IHOTSTP+6) NSTAVX
               WRITE(hss % lun,REC=IHOTSTP+7) NHASE
               WRITE(hss % lun,REC=IHOTSTP+8) NHASV
               WRITE(hss % lun,REC=IHOTSTP+9) NHAGE
               WRITE(hss % lun,REC=IHOTSTP+10) NHAGV
               WRITE(hss % lun,REC=IHOTSTP+11) ICALL
               WRITE(hss % lun,REC=IHOTSTP+12) NFREQ
               IHOTSTP = IHOTSTP+12
C
               do i=1,nfreq+nf
                  FNAME=NAMEFR(I)
                  WRITE(hss % lun,REC=IHOTSTP+1) FNAM8(1)
                  WRITE(hss % lun,REC=IHOTSTP+2) FNAM8(2)
                  IHOTSTP = IHOTSTP + 2
                  WRITE(hss % lun,REC=IHOTSTP+1) HAFREQ(i)
                  WRITE(hss % lun,REC=IHOTSTP+2) HAFF(i)
                  WRITE(hss % lun,REC=IHOTSTP+3) HAFACE(i)
                  IHOTSTP=IHOTSTP+3
               end do
C
C              Write Out time of most recent H.A. update
               WRITE(hss % lun,REC=IHOTSTP+1) TIMEUD
               WRITE(hss % lun,REC=IHOTSTP+2) ITUD
               IHOTSTP=IHOTSTP+2
c
c              Write Out LHS Matrix
c
               IHOTSTP=IHOTSTP+1
               do i=1,mm
                  do j=1,mm
                     WRITE(hss % lun,REC=IHOTSTP) HA(I,J)
                     IHOTSTP = IHOTSTP + 1
                  END DO
               END DO
C
C              In serial, write load vectors to hotstart file.
#ifndef CMPI
               IF(NHASE.ne.0) THEN
                  do n=1,NSTAE
                     do i=1,mm
                        WRITE(hss % lun,REC=IHOTSTP) STAELV(I,N)
                        IHOTSTP=IHOTSTP+1
                     end do
                  end do
               ENDIF
C
               IF(NHASV.ne.0) THEN
                  do N=1,NSTAV
                     do i=1,mm
                        WRITE(hss % lun,REC=IHOTSTP) STAULV(I,N)
                        IHOTSTP=IHOTSTP+1
                        WRITE(hss % lun,REC=IHOTSTP) STAVLV(I,N)
                        IHOTSTP=IHOTSTP+1
                     end do
                  end do
               ENDIF
C
               IF(NHAGE.ne.0) THEN
                  do n=1,np
                     do i=1,mm
                        WRITE(hss % lun,REC=IHOTSTP) GLOELV(I,N)
                        IHOTSTP=IHOTSTP+1
                     end do
                  end do
               ENDIF
C
               IF(NHAGV.ne.0) THEN
                  do n=1,np
                     do i=1,mm
                        WRITE(hss % lun,REC=IHOTSTP) GLOULV(I,N)
                        IHOTSTP=IHOTSTP+1
                        WRITE(hss % lun,REC=IHOTSTP) GLOVLV(I,N)
                        IHOTSTP=IHOTSTP+1
                     end do
                  end do
               ENDIF
#else
C
C              In parallel, write the full domain load vectors to hotstart
C              file.
               IF(NHASE.ne.0) THEN
                  do n=1,NSTAE_G
                     do i=1,mm
                        WRITE(hss % lun,REC=IHOTSTP)
     &                           HarmElevSLVDescript % array2D_g(i,n)
                        IHOTSTP=IHOTSTP+1
                     end do
                  end do
               ENDIF
C
               IF(NHASV.ne.0) THEN
                  do N=1,NSTAV_G
                     do i=1,mm
                        WRITE(hss % lun,REC=IHOTSTP)
     &                              HarmUVelSLVDescript % array2D_g(i,n)
                        IHOTSTP=IHOTSTP+1
                        WRITE(hss % lun,REC=IHOTSTP)
     &                              HarmVVelSLVDescript % array2D_g(i,n)
                        IHOTSTP=IHOTSTP+1
                     end do
                  end do
               ENDIF
C
               IF(NHAGE.ne.0) THEN
                  do n=1,np_g
                     do i=1,mm
                        WRITE(hss % lun,REC=IHOTSTP)
     &                            HarmElevFDLVDescript % array2D_g(i,n)
                        IHOTSTP=IHOTSTP+1
                     end do
                  end do
               ENDIF
C
               IF(NHAGV.ne.0) THEN
                  do n=1,np_g
                     do i=1,mm
                        WRITE(hss % lun,REC=IHOTSTP)
     &                          HarmUVelFDLVDescript % array2D_g(i,n)
                        IHOTSTP=IHOTSTP+1
                        WRITE(hss % lun,REC=IHOTSTP)
     &                           HarmVVelFDLVDescript % array2D_g(i,n)
                        IHOTSTP=IHOTSTP+1
                     end do
                  end do
               ENDIF
#endif
            ENDIF ! myProc.eq.0 .or. write_local_hotstart_files
         ENDIF  ! IHARIND
C
         if (CHARMV.and.
     &           ((MNPROC.eq.1).or.WRITE_LOCAL_HOT_START_FILES)) then
            IF((IHARIND.EQ.1).AND.(IT.GT.ITMV)) THEN

               WRITE(hss % lun,REC=IHOTSTP) NTSTEPS
               IHOTSTP=IHOTSTP+1
               IF(NHAGE.ne.0) THEN
                  DO I=1,NP
                     WRITE(hss % lun,REC=IHOTSTP) ELAV(I)
                     IHOTSTP = IHOTSTP + 1
                     WRITE(hss % lun,REC=IHOTSTP) ELVA(I)
                     IHOTSTP = IHOTSTP + 1
                  END DO
               ENDIF
               IF(NHAGV.ne.0) THEN
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
C
         if (CHARMV.and.(MNPROC.gt.1).and.(MYPROC.eq.0)) then
            IF((IHARIND.EQ.1).AND.(IT.GT.ITMV)) THEN
               WRITE(hss % lun,REC=IHOTSTP) NTSTEPS
               IHOTSTP=IHOTSTP+1
               IF(NHAGE.ne.0) THEN
                  DO I=1,NP_G
                     WRITE(hss % lun,REC=IHOTSTP)
     &                  ELAVDescript % array_g(i)
                     IHOTSTP = IHOTSTP + 1
                     WRITE(hss % lun,REC=IHOTSTP)
     &                  ELVADescript % array_g(i)
                     IHOTSTP = IHOTSTP + 1
                  END DO
               ENDIF
               IF(NHAGV.ne.0) THEN
                  DO I=1,NP_G
                     WRITE(hss % lun,REC=IHOTSTP)
     &                  XVELAVDescript % array_g(i)
                     IHOTSTP=IHOTSTP+1
                     WRITE(hss % lun,REC=IHOTSTP)
     &                  YVELAVDescript % array_g(i)
                     IHOTSTP=IHOTSTP+1
                     WRITE(hss % lun,REC=IHOTSTP)
     &                  XVELVADescript % array_g(i)
                     IHOTSTP=IHOTSTP+1
                     WRITE(hss % lun,REC=IHOTSTP)
     &                  YVELVADescript % array_g(i)
                     IHOTSTP=IHOTSTP+1
                  END DO
               ENDIF
            ENDIF
         endif                   !  charmv
C
         IF ((MYPROC.eq.0).or.
     &         (WRITE_LOCAL_HOT_START_FILES.eqv..TRUE.)) THEN
            CLOSE(hss%lun)
         ENDIF

      CASE(2,267,268) ! ascii hotstart file, only used for debugging
C
C        jgf49.35 write both subdomain and fulldomain hotstart files,
C        just to have them, since this is just for debugging purposes
C        (i.e., to check the values written to hotstart files in other
C        formats) anyway.
         OPEN(hss%lun,FILE=TRIM(INPUTDIR)//'/'//hss % filename,
     &      STATUS='REPLACE')
         WRITE(hss%lun,*) "FileFmtVersion = ",FileFmtVersion
         WRITE(hss%lun,*) "IM = ",IM
         WRITE(hss%lun,*) "TIME = ",TimeLoc
         WRITE(hss%lun,*) "IT = ",IT
         WRITE(hss%lun,*) "NP = ",NP
         WRITE(hss%lun,*) "NE = ",NE
         WRITE(hss%lun,*) "Elev1:"
         DO I=1, Elev1Descript % num_records_this
            WRITE(hss%lun,*) i," ",Elev1Descript % array(I)
         ENDDO
         WRITE(hss%lun,*) "Elev2:"
         DO I=1, Elev2Descript % num_records_this
            WRITE(hss%lun,*) i," ",Elev2Descript % array(I)
         ENDDO
         WRITE(hss%lun,*) "EtaDisc:"
         DO I=1, EtaDiscDescript % num_records_this
            WRITE(hss%lun,*) i," ",EtaDiscDescript % array(I)
         ENDDO
         WRITE(hss%lun,*) "VelU:"
         DO I=1, VelDescript % num_records_this
            WRITE(hss%lun,*) i," ",VelDescript % array(I)
         ENDDO
         WRITE(hss%lun,*) "VelV:"
         DO I=1, VelDescript % num_records_this
            WRITE(hss%lun,*) i," ",VelDescript % array2(I)
         ENDDO
         IF (IM.eq.10) THEN
            WRITE(hss%lun,*) "CH1:"
            DO I=1, CH1Descript % num_records_this
               WRITE(hss%lun,*) i," ",CH1Descript % array(I)
            ENDDO
         ENDIF
         WRITE(hss%lun,*) "NodeCode:"
         DO I=1, NodeCodeDescript % num_records_this
            WRITE(hss%lun,*) i," ",NodeCodeDescript % iarray(I)
         ENDDO
         WRITE(hss%lun,*) "NOFF:"
         DO I=1, NOFFDescript % num_records_this
            WRITE(hss%lun,*) i," ",NOFFDescript % iarray(I)
         ENDDO
         WRITE(hss%lun,*) "IESTP = ",IESTP
         WRITE(hss%lun,*) "NSCOUE = ",NSCOUE
         WRITE(hss%lun,*) "IVSTP = ",IVSTP
         WRITE(hss%lun,*) "NSCOUV = ",NSCOUV
         WRITE(hss%lun,*) "ICSTP = ",ICSTP
         WRITE(hss%lun,*) "NSCOUC = ",NSCOUC
         WRITE(hss%lun,*) "IPSTP = ",IPSTP
         WRITE(hss%lun,*) "IWSTP = ",IWSTP
         WRITE(hss%lun,*) "NSCOUM = ",NSCOUM
         WRITE(hss%lun,*) "IGEP = ",IGEP
         WRITE(hss%lun,*) "NSCOUGE = ",NSCOUGE
         WRITE(hss%lun,*) "IGVP = ",IGVP
         WRITE(hss%lun,*) "NSCOUGV = ",NSCOUGV
         WRITE(hss%lun,*) "IGCP = ",IGCP
         WRITE(hss%lun,*) "NSCOUGC = ",NSCOUGC
         WRITE(hss%lun,*) "IGPP = ",IGPP
         WRITE(hss%lun,*) "IGWP = ",IGWP
         WRITE(hss%lun,*) "NSCOUGW = ",NSCOUGW
         CLOSE(hss%lun)
C
C
         IF ((MNPROC.gt.1).and.(MYPROC.eq.0)) THEN
            OPEN(hss%lun,FILE=TRIM(GBLINPUTDIR)//'/'//hss % filename,
     &         STATUS='REPLACE')
            WRITE(hss%lun,*) "FileFmtVersion = ",FileFmtVersion
            WRITE(hss%lun,*) "IM = ",IM
            WRITE(hss%lun,*) "TIME = ",TimeLoc
            WRITE(hss%lun,*) "IT = ",IT
            WRITE(hss%lun,*) "NP = ",NP_G
            WRITE(hss%lun,*) "NE = ",NE_G
            WRITE(hss%lun,*) "Elev1:"
            DO I=1, Elev1Descript % num_fd_records
               WRITE(hss%lun,*) i," ",Elev1Descript % array_g(I)
            ENDDO
            WRITE(hss%lun,*) "Elev2:"
            DO I=1, Elev2Descript % num_fd_records
               WRITE(hss%lun,*) i," ",Elev2Descript % array_g(I)
            ENDDO
            WRITE(hss%lun,*) "EtaDisc:"
            DO I=1, EtaDiscDescript % num_fd_records
               WRITE(hss%lun,*) i," ",EtaDiscDescript % array_g(I)
            ENDDO
            WRITE(hss%lun,*) "VelU:"
            DO I=1, VelDescript % num_fd_records
               WRITE(hss%lun,*) i," ",VelDescript % array_g(I)
            ENDDO
            WRITE(hss%lun,*) "VelV:"
            DO I=1, VelDescript % num_fd_records
               WRITE(hss%lun,*) i," ",VelDescript % array2_g(I)
            ENDDO
            IF (IM.eq.10) THEN
               WRITE(hss%lun,*) "CH1:"
               DO I=1, CH1Descript % num_fd_records
                  WRITE(hss%lun,*) i," ",CH1Descript % array_g(I)
               ENDDO
            ENDIF
            WRITE(hss%lun,*) "NodeCode:"
            DO I=1, NodeCodeDescript % num_fd_records
               WRITE(hss%lun,*) i," ",NodeCodeDescript % iarray_g(I)
            ENDDO
            WRITE(hss%lun,*) "NOFF:"
            DO I=1, NOFFDescript % num_fd_records
               WRITE(hss%lun,*) i," ",NOFFDescript % iarray_g(I)
            ENDDO
            WRITE(hss%lun,*) "IESTP = ",IESTP
            WRITE(hss%lun,*) "NSCOUE = ",NSCOUE
            WRITE(hss%lun,*) "IVSTP = ",IVSTP
            WRITE(hss%lun,*) "NSCOUV = ",NSCOUV
            WRITE(hss%lun,*) "ICSTP = ",ICSTP
            WRITE(hss%lun,*) "NSCOUC = ",NSCOUC
            WRITE(hss%lun,*) "IPSTP = ",IPSTP
            WRITE(hss%lun,*) "IWSTP = ",IWSTP
            WRITE(hss%lun,*) "NSCOUM = ",NSCOUM
            WRITE(hss%lun,*) "IGEP = ",IGEP
            WRITE(hss%lun,*) "NSCOUGE = ",NSCOUGE
            WRITE(hss%lun,*) "IGVP = ",IGVP
            WRITE(hss%lun,*) "NSCOUGV = ",NSCOUGV
            WRITE(hss%lun,*) "IGCP = ",IGCP
            WRITE(hss%lun,*) "NSCOUGC = ",NSCOUGC
            WRITE(hss%lun,*) "IGPP = ",IGPP
            WRITE(hss%lun,*) "IGWP = ",IGWP
            WRITE(hss%lun,*) "NSCOUGW = ",NSCOUGW
            CLOSE(hss%lun)
         ENDIF
C
C        TODO: 3D data in ascii
C        TODO: harmonic analysis data in ascii
C
      CASE(3,5,367,368,567,568) ! netcdf

#ifdef ADCNETCDF
         IF (myProc.eq.0) THEN
            CALL writeNetCDFHotstart(hss%lun, Elev1Descript,
     &         Elev2Descript, VelDescript, CH1Descript, EtaDiscDescript,
     &         NodeCodeDescript, NOFFDescript, TimeLoc, it)
            IF ((IHARIND.eq.1).and.(IT.gt.ITHAS)) THEN
               CALL writeNetCDFHotstartHarmonic(hss%lun,
     &            HarmElevFDLVDescript, HarmElevSLVDescript,
     &            HarmUVelFDLVDescript, HarmVVelFDLVDescript,
     &            HarmUVelSLVDescript, HarmVVelSLVDescript)
            ENDIF
            IF (CHARMV.eqv..true.) THEN
               CALL writeNetCDFHotstartHarmonicMeansVariances(hss%lun,
     &            ELAVDescript, ELVADescript, XVELAVDescript,
     &            YVELAVDescript, XVELVADescript, YVELVADescript)
            ENDIF
         ENDIF
         IF (C3D.eqv..true.) THEN
            CALL HSTART3D_OUT(IT) ! -> see vsmy.F
         ENDIF
#endif

      CASE DEFAULT
         write(ScreenUnit,*) 'The value of NHSTAR=',NHSTAR
         write(ScreenUnit,*) 'is not supported.'
         write(ScreenUnit,*) 'Hotstart file not written.'

      END SELECT

      write(scratchMessage,24541) hss % lun,IT,TimeLoc
      call allMessage(INFO,scratchMessage)
24541 FORMAT(1X,'HOT START OUTPUT WRITTEN TO UNIT ',I2,
     &    ' AT TIME STEP = ',I9,' TIME = ',E15.8)
      TEST_HOTSTART = .FALSE.
      IF(hss % lun.EQ.67) THEN
C        jgf45.07 added option to stop ADCIRC after writing hot start file.
         IF ((NHSTAR.EQ.67).OR.(NHSTAR.EQ.267).OR.(NHSTAR.EQ.367).OR.
     &        (NHSTAR.EQ.567)) THEN
            call allMessage(INFO,"NHSTAR.EQ.67; ADCIRC stopping.")
            TEST_HOTSTART = .TRUE.
         ENDIF
         hss % lun      = 68
         IF((NHSTAR.EQ.1).OR.(NHSTAR.EQ.2)) THEN
            hss % filename = 'fort.68'
         ENDIF
      ELSE
         IF ((NHSTAR.EQ.68).OR.(NHSTAR.EQ.268).OR.(NHSTAR.EQ.368).OR.
     &       (NHSTAR.EQ.568)) THEN
            call allMessage(INFO,"NHSTAR.EQ.68; ADCIRC stopping.")
            TEST_HOTSTART = .TRUE.
         ENDIF
         hss % lun=67
         IF((NHSTAR.EQ.1).OR.(NHSTAR.EQ.2)) THEN
            hss % filename = 'fort.67'
         ENDIF
      ENDIF
      IF (TEST_HOTSTART.eqv..true.) THEN
#ifdef CMPI
          CALL MSG_FINI()
#endif
         STOP
      ENDIF

#ifdef WRITE_OUTPUT_TRACE
      call screen(DEBUG,"Return.")
#endif
      call unsetMessageSource()
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
      SUBROUTINE WriteWarnElev(TimeLoc, it)
      USE SIZES
      USE GLOBAL
      USE GLOBAL_IO, ONLY: collectFullDomainArray, unpackOne, pack63
#ifdef CMPI
      USE MESSENGER
#endif
      IMPLICIT NONE
      REAL(8), intent(in) :: TimeLoc
      INTEGER, intent(in) :: it
C
      INTEGER I ! node loop counter
      type(OutputDataDescript_t), SAVE :: DebugElevDescript
      REAL(SZ), SAVE, ALLOCATABLE, TARGET :: ETA_DEBUG_g(:)
C
C     Open the file and write the header on the first warning. Format
C     the actual header data to make it look more or less like a fort.63
C     file.
      IF (WarnElevDumpCounter.eq.0) THEN
         OPEN(69,FILE=TRIM(LOCALDIR)//'/'//'fort.69')
         WRITE(69,'(A80)') TRIM(RUNDES)
         WRITE(69,6000) WarnElevDumpLimit,NP,DTDP*1.d0,1,1
         CLOSE(69)
#ifdef CMPI
         !jgf49.0601: if running in parallel, also initialize a fulldomain file
         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            ALLOCATE(ETA_DEBUG_g(NP_G))
         ENDIF
         DebugElevDescript % specifier            = 1 ! TODO: make configurable
         DebugElevDescript % initial_value        = 0.0
         DebugElevDescript % num_items_per_record = 1
         DebugElevDescript % num_fd_records       = NP_G
         DebugElevDescript % num_records_this     = NP
         DebugElevDescript % imap                 => NODES_LG
         DebugElevDescript % array                => ETA2
         DebugElevDescript % array_g              => ETA_DEBUG_g
         DebugElevDescript % ConsiderWetDry       = .TRUE.
         DebugElevDescript % alternate_value      = -99999.0
         DebugElevDescript % field_name           = 'DebugElev'

         OPEN(69,FILE=TRIM(GLOBALDIR)//'/'//'fort.69')
         WRITE(69,'(A80)') TRIM(RUNDES)
         WRITE(69,6000) WarnElevDumpLimit,NP_G,DTDP*1.d0,1,1
         CLOSE(69)
#endif
      ENDIF
 6000 FORMAT(1x,I10,1x,I10,1x,1pE15.7E3,1x,I5,1x,I5)
C
C     Open the file and append the new elevation data.
      OPEN(69,FILE=TRIM(LOCALDIR)//'/'//'fort.69',
     &     ACCESS='SEQUENTIAL',POSITION='APPEND')
      WRITE(69,2120) TimeLoc,IT
 2120 FORMAT(2X,1pE20.10E3,5X,I10)
      DO I=1,NP
         IF(NODECODE(I).EQ.1)THEN
            WRITE(69,2453) I,ETA2(I)
         ELSE
            WRITE(69,2453) I,-99999D0
         ENDIF
 2453    FORMAT(2X,I8,2X,1pE20.10E3)
      ENDDO
      CLOSE(69)


#ifdef CMPI
      !jgf49.0601: if running in parallel, also append to fulldomain file
      CALL writeOutArray(69, TimeLoc, IT, DebugElevDescript,
     &                         pack63, unpackOne, IGEP,
     &                         'fort.69        ')
#endif
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

