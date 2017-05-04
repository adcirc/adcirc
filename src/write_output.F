C----------------------------------------------------------------------
C      M O D U L E    W R I T E   O U T P U T
C----------------------------------------------------------------------
C jgf48.03 This module contains all the subroutines that write output to
C files. It was created to consolidate all the output routines that
C were located in other parts of the code, to make it easier to incorporate
C them into globalio, and to make it easier to add new output modes or
C formats in the future.
C
C tcm48.4618 Fixed Bug in writing out global elevation values that are
C dry.
C----------------------------------------------------------------------
      MODULE WRITE_OUTPUT
C----------------------------------------------------------------------
      USE SIZES, ONLY : SZ
      USE GLOBAL, ONLY : OutputDataDescript_t, setMessageSource,
     &  unsetMessageSource, DEBUG, allMessage, screenMessage, ERROR,
     &  netcdf_avail

      !
      !    2 D   D A T A
      !
#ifdef CSWAN
      integer, parameter :: numOutputDescript2D = 46
#else
      integer, parameter :: numOutputDescript2D = 32
#endif
      type descript2D_ptr_t
         type(OutputDataDescript_t), pointer :: descript2D
      end type descript2D_ptr_t
      type(descript2D_ptr_t), allocatable :: ptr(:)
      !
      type(OutputDataDescript_t), target :: ElevStaDescript
      type(OutputDataDescript_t), target :: VelStaDescript
      type(OutputDataDescript_t), target :: ElevDescript
      type(OutputDataDescript_t), target :: Tau0Descript
      type(OutputDataDescript_t), target :: VelDescript
      type(OutputDataDescript_t), target :: PrStaDescript
      type(OutputDataDescript_t), target :: WindVelStaDescript
      type(OutputDataDescript_t), target :: PrDescript
      type(OutputDataDescript_t), target :: WindVelDescript
      type(OutputDataDescript_t), target :: ConcStaDescript
      type(OutputDataDescript_t), target :: ConcDescript
      !tcm v49.64.01 added for ice
      type(OutputDataDescript_t), target :: IceDescript
      !tcm v49.64.01 added for ice stations
      type(OutputDataDescript_t), target :: IceStaDescript
      type(OutputDataDescript_t), target :: EtaMaxDescript
      type(OutputDataDescript_t), target :: UMaxDescript
      type(OutputDataDescript_t), target :: WVMaxDescript
      type(OutputDataDescript_t), target :: PrMinDescript
      type(OutputDataDescript_t), target :: RSMaxDescript
      type(OutputDataDescript_t), target :: SpongeDescript
      type(OutputDataDescript_t), target :: ESLDescript
      !tcm v50.66.01 added time varying bathy
      type(OutputDataDescript_t), target :: BathyDescript
      !tcm v50.66.01 added time varying bathy stations
      type(OutputDataDescript_t), target :: BathyStaDescript
      ! Casey 090302: Added this type for output of radiation stress gradients.
      ! tcm v50.75 removed ifdef cswan to allow for use whenever nrs=3 or nrs=4
      type(OutputDataDescript_t), target :: RSDescript
      !jgf52.05: Added for regular output of wet/dry quantities.
      type(OutputDataDescript_t), target :: noffOutDescript
      type(OutputDataDescript_t), target :: nodeCodeOutDescript
Cobell Added a descript for time varying weir tracking
      type(OutputDataDescript_t), target :: weirElevDescript
      ! @jasonfleming: adding output for wind drag
      type(OutputDataDescript_t), target :: windDragDescript

C     jgf48.03 Full domain arrays used when writing globalio output.
      real(sz), allocatable, target :: et00_g(:)
      real(sz), allocatable, target :: uu00_g(:)
      real(sz), allocatable, target :: vv00_g(:)
      real(sz), allocatable, target :: tau0var_g(:)
      real(sz), allocatable, target :: rmp00_g(:)
      real(sz), allocatable, target :: rmu00_g(:)
      real(sz), allocatable, target :: rmv00_g(:)
      real(sz), allocatable, target :: pr2_g(:)
      real(sz), allocatable, target :: wvnxout_g(:)
      real(sz), allocatable, target :: wvnyout_g(:)
      real(sz), allocatable, target :: cc00_g(:)
      real(sz), allocatable, target :: ciceout_g(:)  !tcm v49.64.01 added for ice
      real(sz), allocatable, target :: rmice00_g(:)  !tcm v49.64.01 added for ice
      real(sz), allocatable, target :: etamax_g(:)
      real(sz), allocatable, target :: umax_g(:)
      real(sz), allocatable, target :: prmin_g(:)
      real(sz), allocatable, target :: wvnoutmax_g(:)
      real(sz), allocatable, target :: rsnmax_g(:)
      real(sz), allocatable, target :: sponge_g(:)
      real(sz), allocatable, target :: eslonoff_g(:) !zc added for esl on/off output
      real(sz), allocatable, target :: dpout_g(:)    !tcm v50.66.01 added time varying bathy
      real(sz), allocatable, target :: dp00_g(:)  !tcm v50.66.01 added time varying bathy
      ! Casey 090302: Added these arrays for output of radiation stress gradients.
      ! tcm v50.75 removed ifdef cswan to allow for use whenever nrs=3 or nrs=4
      real(sz), allocatable, target :: rsnxout_g(:)
      real(sz), allocatable, target :: rsnyout_g(:)
      !jgf52.05: Added for normal output of wet/dry quantities (nodecode and noff)
      integer, allocatable, target :: nodecode_out_g(:)
      integer, allocatable, target :: noff_out_g(:)
      real(sz), allocatable, target :: windDragOut_g(:)
      !
      !  H A R M O N I C   A N A L Y S I S   D A T A
      !
      type(OutputDataDescript_t) :: HAElevStaMagDescript
      type(OutputDataDescript_t) :: HAElevStaPhaseDescript
      type(OutputDataDescript_t) :: HAVelStaUMagDescript
      type(OutputDataDescript_t) :: HAVelStaVMagDescript
      type(OutputDataDescript_t) :: HAVelStaUPhaseDescript
      type(OutputDataDescript_t) :: HAVelStaVPhaseDescript
      type(OutputDataDescript_t) :: HAElevMagDescript
      type(OutputDataDescript_t) :: HAElevPhaseDescript
      type(OutputDataDescript_t) :: HAVelUMagDescript
      type(OutputDataDescript_t) :: HAVelVMagDescript
      type(OutputDataDescript_t) :: HAVelUPhaseDescript
      type(OutputDataDescript_t) :: HAVelVPhaseDescript
      type(OutputDataDescript_t) :: MVDescript
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
      !
      !    3 D   D A T A
      !
      type(OutputDataDescript_t) :: SigTStaDescript
      type(OutputDataDescript_t) :: SalStaDescript
      type(OutputDataDescript_t) :: TempStaDescript
      type(OutputDataDescript_t) :: QSurfKp1Descript
      type(OutputDataDescript_t) :: RealQStaDescript
      type(OutputDataDescript_t) :: ImaginaryQStaDescript
      type(OutputDataDescript_t) :: WZStaDescript
      type(OutputDataDescript_t) :: Q20StaDescript
      type(OutputDataDescript_t) :: LStaDescript
      type(OutputDataDescript_t) :: EVStaDescript
      type(OutputDataDescript_t) :: SigTDescript
      type(OutputDataDescript_t) :: SalDescript
      type(OutputDataDescript_t) :: TempDescript
      type(OutputDataDescript_t) :: RealQDescript
      type(OutputDataDescript_t) :: ImaginaryQDescript
      type(OutputDataDescript_t) :: WZDescript
      type(OutputDataDescript_t) :: Q20Descript
      type(OutputDataDescript_t) :: LDescript
      type(OutputDataDescript_t) :: EVDescript
C
C     parts of complex variables
      REAL(SZ), ALLOCATABLE, TARGET :: rp(:,:) !real part subdomain data
      REAL(SZ), ALLOCATABLE, TARGET :: ip(:,:) !imaginary part subdom dat
      REAL(SZ), ALLOCATABLE, TARGET :: rp_g(:,:) !real part fulldomain
      REAL(SZ), ALLOCATABLE, TARGET :: ip_g(:,:) !imaginary part fulldom
      REAL(SZ), ALLOCATABLE, TARGET :: rpSta(:,:) !real part subd sta
      REAL(SZ), ALLOCATABLE, TARGET :: ipSta(:,:) !im part subdom sta
      REAL(SZ), ALLOCATABLE, TARGET :: rpSta_g(:,:) !real part fd sta
      REAL(SZ), ALLOCATABLE, TARGET :: ipSta_g(:,:) !im part fd sta
C
      REAL(SZ), ALLOCATABLE, TARGET :: qsurfkp1_g(:)
C
      REAL(SZ), ALLOCATABLE, TARGET :: q20Sta_g(:,:)
      REAL(SZ), ALLOCATABLE, TARGET :: lSta_g(:,:)
      REAL(SZ), ALLOCATABLE, TARGET :: EVSta_g(:,:)
      REAL(SZ), ALLOCATABLE, TARGET :: SigTSta_g(:,:)
      REAL(SZ), ALLOCATABLE, TARGET :: SalSta_g(:,:)
      REAL(SZ), ALLOCATABLE, TARGET :: TempSta_g(:,:)
      REAL(SZ), ALLOCATABLE, TARGET :: WZSta_g(:,:)
      REAL(SZ), ALLOCATABLE, TARGET :: TVW_g(:)
      !
      !  H O T S T A R T   D A T A
      !
      type(OutputDataDescript_t) :: Elev1Descript
      type(OutputDataDescript_t) :: Elev2Descript
      type(OutputDataDescript_t) :: HotstartVelDescript
      type(OutputDataDescript_t) :: CH1Descript
      type(OutputDataDescript_t) :: EtaDiscDescript
      type(OutputDataDescript_t) :: NodeCodeDescript
      type(OutputDataDescript_t) :: NOFFDescript
      type(OutputDataDescript_t) :: HarmElevFDLVDescript
      type(OutputDataDescript_t) :: HarmElevSLVDescript
      type(OutputDataDescript_t) :: HarmUVelFDLVDescript
      type(OutputDataDescript_t) :: HarmVVelFDLVDescript
      type(OutputDataDescript_t) :: HarmUVelSLVDescript
      type(OutputDataDescript_t) :: HarmVVelSLVDescript
      type(OutputDataDescript_t) :: ELAVDescript
      type(OutputDataDescript_t) :: ELVADescript
      type(OutputDataDescript_t) :: XVELAVDescript
      type(OutputDataDescript_t) :: YVELAVDescript
      type(OutputDataDescript_t) :: XVELVADescript
      type(OutputDataDescript_t) :: YVELVADescript
      !
      !  S W A N  O U T P U T  D A T A
      !
#ifdef CSWAN
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
      !...Spool Counters
      INTEGER, TARGET                          :: SWAN_HS_SPOOL,
     &                                            SWAN_TPS_SPOOL,
     &                                            SWAN_DIR_SPOOL,
     &                                            SWAN_TM01_SPOOL,
     &                                            SWAN_TM02_SPOOL,
     &                                            SWAN_TMM10_SPOOL,
     &                                            SWAN_WIND_SPOOL
      !...File position
      INTEGER, TARGET                          :: SWAN_HS_POS,
     &                                            SWAN_TPS_POS,
     &                                            SWAN_DIR_POS,
     &                                            SWAN_TM01_POS,
     &                                            SWAN_TM02_POS,
     &                                            SWAN_TMM10_POS,
     &                                            SWAN_WIND_POS
#endif
      !
      !   I N U N D A T I O N   O U T P U T   D A T A
      !
      ! Casey 120523: Add Nishant's changes for the Texas ASGS.
      ! jgf52.08.01: Additions for detailed inundation output products.
      !
      ! initiallyDry is 1 for a node that is dry at cold start.
      integer, allocatable, target :: initiallyDry(:)
      integer, allocatable, target :: initiallyDry_g(:)
      type(OutputDataDescript_t), target :: initiallyDryDescript
      !
      ! endRisingInun is 1 if the total water depth increases at a node
      ! that was initiallyDry, if that increase occurs between the
      ! penultimate and ultimate time steps of the simulation (i.e.,
      ! the water is still rising when the simulation ends)
      integer, allocatable, target :: endRisingInun(:)
      integer, allocatable, target :: endRisingInun_g(:)
      type(OutputDataDescript_t), target :: endRisingInunDescript
      !
      ! maxInunDepth is used for a new min/max file, where the maximum
      ! total water column depth of inundation is recorded.
      real(sz), allocatable, target :: maxInunDepth(:)
      real(sz), allocatable, target :: maxInunDepth_time(:)
      real(sz), allocatable, target :: maxInunDepth_g(:)
      real(sz), allocatable, target :: maxInunDepth_time_g(:)
      type(OutputDataDescript_t), target :: maxInunDepthDescript
      !
      ! penultimateNodeCode is used as a temp array to record the
      ! wet/dry state on the 2nd to last time step
      integer, allocatable, target :: penultimateNodeCode(:)
      !
      ! inundationTime records the total time an initiallyDry node
      ! has been wet (deeper than thresh)
      type(OutputDataDescript_t), target :: inundationTimeDescript
      real(sz), allocatable, target :: inundationTime(:)
      ! record the time of onset of inundation beyond the threshold
      real(sz), allocatable, target :: inundationTime_onset(:)
      real(sz), allocatable, target :: inundationTime_g(:)
      real(sz), allocatable, target :: inundationTime_onset_g(:)
      !
      ! everDried is used to record whether a node has ever become dry,
      ! which is useful to know when computing tidal phases and amplitudes
      real(sz), allocatable, target :: everDried(:)
      real(sz), allocatable, target :: everDried_time(:)
      real(sz), allocatable, target :: everDried_g(:)
      real(sz), allocatable, target :: everDried_time_g(:)
      type(OutputDataDescript_t), target :: everDriedDescript

      ! jgf52.05: Set via fort.15 namelist in read_input.F.
      logical :: outputNodeCode = .false.
      logical :: outputNOFF = .false.

C     ---------
      CONTAINS
C     ---------

C----------------------------------------------------------------------
C       S U B R O U T I N E   I N I T   O U T P U T   2 D
C----------------------------------------------------------------------
C     jgf51.21.24: Initialize 2D output data structures.
C----------------------------------------------------------------------
      SUBROUTINE initOutput2D(timeloc)
      USE SIZES, ONLY : SZ, INPUTDIR, NBYTE, MNWPROC, MYPROC, MNPROC,
     &                  GLOBALDIR, OFF, ASCII, NETCDF3, NETCDF4, XDMF,
     &                  numFormats, write_local_files, localdir,
     &                  controlFileName
      USE GLOBAL, ONLY : imap_stae_lg, imap_stav_lg, et00, eta2,
     &   nscoue, iestp, nne, staie1, staie2, staie3, xel, yel,
     &   slel, sfel, uu00, vv00, uu2, vv2, nscouv, ivstp, nnv,
     &   staiv1, staiv2, staiv3, xev, yev, slev, sfev, nscouge,
     &   igep, nodes_lg, eta2_g, nscougt, itau0p, nscougs, igsp,
     &   sponge, nscougv, igvp, uu2_g, vv2_g, imap_stam_lg, rmp00,
     &   pr2, nscoum, ipstp, nnm, staim1, staim2, staim3,
     &   xem, yem, slem, sfem, rmu00, rmv00,
     &   iwstp, nscougw, igpp, igwp, dp00, ibstp, igbp, igip,
     &   ciceout, rmice00, iicestp, igrads, rsnxout, rsnyout,
     &   imap_stac_lg, cc00, ch1, xec, yec, slec, sfec, ch1_g,
     &   nscougc, igcp, ch1, ch1_g, etamax_time_g, etamax,
     &   etamax_time, umax_time_g, umax, umax_time, prmin_time_g,
     &   prmin, wvnoutmax_time_g, wvnoutmax_time, rsnmax_time_g,
     &   rsnmax_time, eslonoff, prmin_time, wvnoutmax, rsnmax,
     &   nstae_g, noute, nstae, wvnxout, wvnyout,
     &   ntcyse, ntcyfe, nspoole, nstav_g, noutv, nstav, ntcysv,
     &   ntcyfv, nspoolv, noutge, ntcysge, ntcyfge, nspoolge,
     &   np_g, outputsponge, noutgv, ntcysgv, nspoolgv, nstam_g,
     &   noutm, nstam, ntcyfgv, ntcysm, ntcyfm, nspoolm, nws,
     &   noutgw, ntcysgw, ntcyfgw, nspoolgw, nddt, ncice, nrs,
     &   noutc, nstac_g, nstac, im, noutgc, ntcyfgc, ntcysgc,
     &   nspoolgc, INFO, scratchMessage, allMessage, nscougi, nscougp,
     &   nscours, nscoub, nscougb, nscoui, nscoup,
     &   use_tvw, nout_tvw, ntcys_tvw, ntcyf_tvw, ig_tvw, tvw,
     &   nspool_tvw, nscou_tvw,
     &   nodeCodeOutSpoolCounter, noffOutSpoolCounter, noffOutFilePos,
     &   nodeCodeOutFilePos, NE_G, IMAP_EL_LG, nodecode, noff,
     &   fauxFilePosition, initiallyDrySpoolCounter,
     &   inundationOutput, ihot, logMessage, windDragFilePosition,
     &   windDragOut, windDragSpoolCounter, outputWindDrag
      USE MESH, ONLY : NP, NE, DP, NM, ICS
      USE NodalAttributes, ONLY : OutputTau0, Tau0Var, LoadEleSlopeLim,
     &  elemental_slope_limiter_grad_max,elemental_slope_limiter_active
      USE GLOBAL_IO, ONLY : readAndMapToSubdomainMaxMin
#ifdef ADCNETCDF
      USE NETCDFIO, ONLY : initNetCDFOutputFile,
     &                      readAndMapToSubdomainMaxMinNetCDF
#endif
#ifdef ADCXDMF
      USE XDMFIO, ONLY : initOutputXDMF, writeControlXDMF
      USE CONTROL, ONLY : readControlFile
#ifdef CMPI
      USE WRITER, ONLY : sendInitWriterXDMF
#endif
#endif
      implicit none
C
      real(8), intent(in) :: timeloc
      character(len=20) :: extString ! string version of integer file extension
      logical :: fileFound = .false. ! .true. if the file exists
      logical :: nerr
      integer :: i
C
      call setMessageSource("initOutput2D")
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      !
      ! jgf52.21.24: Create an array of pointers to all the 2D output
      ! data structures, which enables us to iterate over them; this
      ! simplifies the task of setting default values as well as writing
      ! output
      allocate(ptr(numOutputDescript2D))
      ptr(1)%descript2D => ElevStaDescript
      ptr(2)%descript2D => VelStaDescript
      ptr(3)%descript2D => ElevDescript
      ptr(4)%descript2D => Tau0Descript
      ptr(5)%descript2D => VelDescript
      ptr(6)%descript2D => PrStaDescript
      ptr(7)%descript2D => WindVelStaDescript
      ptr(8)%descript2D => PrDescript
      ptr(9)%descript2D => WindVelDescript
      ptr(10)%descript2D => ConcStaDescript
      ptr(11)%descript2D => ConcDescript
      ptr(12)%descript2D => IceDescript
      ptr(13)%descript2D => IceStaDescript
      ptr(14)%descript2D => EtaMaxDescript
      ptr(15)%descript2D => UMaxDescript
      ptr(16)%descript2D => WVMaxDescript
      ptr(17)%descript2D => PrMinDescript
      ptr(18)%descript2D => RSMaxDescript
      ptr(19)%descript2D => SpongeDescript
      ptr(20)%descript2D => ESLDescript
      ptr(21)%descript2D => BathyDescript
      ptr(22)%descript2D => BathyStaDescript
      ptr(23)%descript2D => RSDescript
      ptr(24)%descript2D => weirElevDescript
      ptr(25)%descript2D => nodeCodeOutDescript
      ptr(26)%descript2D => noffOutDescript
      ptr(27)%descript2D => inundationTimeDescript
      ptr(28)%descript2D => maxInunDepthDescript
      ptr(29)%descript2D => endRisingInunDescript
      ptr(30)%descript2D => initiallyDryDescript
      ptr(31)%descript2D => everDriedDescript
      ptr(32)%descript2D => windDragDescript

!... For SWAN, there are additional output arrays
#ifdef CSWAN
      ptr(33)%descript2D => SwanHSDescript
      ptr(34)%descript2D => SwanHSMaxDescript
      ptr(35)%descript2D => SwanDIRDescript
      ptr(36)%descript2D => SwanDIRMaxDescript
      ptr(37)%descript2D => SwanTM01Descript
      ptr(38)%descript2D => SwanTM01MaxDescript
      ptr(39)%descript2D => SwanTPSDescript
      ptr(40)%descript2D => SwanTPSMaxDescript
      ptr(41)%descript2D => SwanWindDescript
      ptr(42)%descript2D => SwanWindMaxDescript
      ptr(43)%descript2D => SwanTM02Descript
      ptr(44)%descript2D => SwanTM02MaxDescript
      ptr(45)%descript2D => SwanTMM10Descript
      ptr(46)%descript2D => SwanTMM10MaxDescript
#endif
      !
      ! jgf51.21.24: Set some defaults
      do i=1,numOutputDescript2D
         ptr(i) % descript2D % specifier = OFF
         ptr(i) % descript2D % writeFlag = .true.
         ptr(i) % descript2D % initial_value = 0.d0
         ptr(i) % descript2D % num_items_per_record = 1
         ptr(i) % descript2D % num_fd_records = np_g
         ptr(i) % descript2D % num_records_this = np
         ptr(i) % descript2D % imap => nodes_lg
         ptr(i) % descript2D % considerWetDry = .false.
         ptr(i) % descript2D % alternate_value = -99999.0
         ptr(i) % descript2D % isStation = .false.
         ptr(i) % descript2D % divideByDepth = .false.
         ! jgf52.08.03: initiallydry.63 already initialized and written by now
         if (trim(ptr(i)%descript2D%field_name).ne.'initiallyDry') then
            allocate(ptr(i) % descript2D % writerFormats(numFormats))
         endif
         ptr(i) % descript2D % writerFormats(:) = -99999
         ptr(i) % descript2D % useWriter = .false.
         ptr(i) % descript2D % file_extension = -99999
         ptr(i) % descript2D % file_basename = 'fort'
         ptr(i) % descript2D % initialized = .false.
         ptr(i) % descript2D % minmax_timestamp = .false.
         ptr(i) % descript2D % readMaxMin = .false.
         ptr(i) % descript2D % isInteger = .false.
      end do

      ! fort.61
      IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
         ALLOCATE(ET00_g(NSTAE_G))
      ENDIF
      ElevStaDescript % lun                  = 61
      ElevStaDescript % specifier            = NOUTE
      ElevStaDescript % num_fd_records       = NSTAE_G
      ElevStaDescript % num_records_this     = NSTAE
      ElevStaDescript % imap                 => IMAP_STAE_LG
      ElevStaDescript % array                => ET00
      ElevStaDescript % array_g              => ET00_g
      ElevStaDescript % interped_array       => ETA2
      ElevStaDescript % ConsiderWetDry       = .TRUE.
      ElevStaDescript % field_name           = 'ElevSta'
      ElevStaDescript % isStation            = .true.
      ElevStaDescript % startTimeStep           = NTCYSE
      ElevStaDescript % endTimeStep             = NTCYFE
      ElevStaDescript % outputTimeStepIncrement = NSPOOLE
      ElevStaDescript % spoolCounter            => NSCOUE
      ElevStaDescript % filepos                 => IESTP
      ElevStaDescript % elements                => NNE
      ElevStaDescript % interp_fac1             => STAIE1
      ElevStaDescript % interp_fac2             => STAIE2
      ElevStaDescript % interp_fac3             => STAIE3
      IF (ICS.eq.1) THEN
         ElevStaDescript % x_coord              => XEL ! use orig coord
         ElevStaDescript % y_coord              => YEL
      ELSE
         ElevStaDescript % x_coord              => SLEL ! radians
         ElevStaDescript % y_coord              => SFEL
      ENDIF
      !
      ! fort.62
      IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
         ALLOCATE(UU00_g(NSTAV_G))
         ALLOCATE(VV00_g(NSTAV_G))
      ENDIF
      VelStaDescript % lun                  = 62
      VelStaDescript % specifier            = NOUTV
      VelStaDescript % num_items_per_record = 2
      VelStaDescript % num_fd_records       = NSTAV_G
      VelStaDescript % num_records_this     = NSTAV
      VelStaDescript % imap                 => IMAP_STAV_LG
      VelStaDescript % array                => UU00
      VelStaDescript % array2               => VV00
      VelStaDescript % array_g              => UU00_g
      VelStaDescript % array2_g             => VV00_g
      VelStaDescript % interped_array       => UU2
      VelStaDescript % interped_array2      => VV2
      VelStaDescript % field_name           = 'VelSta'
      VelStaDescript % isStation            = .true.
      VelStaDescript % startTimeStep           = NTCYSV
      VelStaDescript % endTimeStep             = NTCYFV
      VelStaDescript % outputTimeStepIncrement = NSPOOLV
      VelStaDescript % spoolCounter            => NSCOUV
      VelStaDescript % filepos                 => IVSTP
      VelStaDescript % elements                => NNV
      VelStaDescript % interp_fac1             => STAIV1
      VelStaDescript % interp_fac2             => STAIV2
      VelStaDescript % interp_fac3             => STAIV3
      IF (ICS.eq.1) THEN
         VelStaDescript % x_coord              => XEV ! use orig coord
         VelStaDescript % y_coord              => YEV
      ELSE
         VelStaDescript % x_coord              => SLEV ! radians
         VelStaDescript % y_coord              => SFEV
      ENDIF

      ! fort.63
      ElevDescript % lun                  = 63
      ElevDescript % specifier            = NOUTGE
      ElevDescript % startTimeStep           = NTCYSGE
      ElevDescript % endTimeStep             = NTCYFGE
      ElevDescript % outputTimeStepIncrement = NSPOOLGE
      ElevDescript % spoolCounter            => NSCOUGE
      ElevDescript % filepos                 => IGEP
      ElevDescript % array                => ETA2
      ElevDescript % array_g              => ETA2_g
      ElevDescript % ConsiderWetDry       = .TRUE.
      ElevDescript % field_name           = 'Elev'
      ElevDescript % writerFormats(1:5)   = (/ 1, 3, 4, 5, 7 /)

      ! tau0 fort.90
      IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
         ALLOCATE(TAU0VAR_g(NP_G))
      ENDIF
      Tau0Descript % lun                  = 90
      Tau0Descript % specifier            = NOUTGE
      Tau0Descript % startTimeStep           = NTCYSGE
      Tau0Descript % endTimeStep             = NTCYFGE
      Tau0Descript % outputTimeStepIncrement = NSPOOLGE
      Tau0Descript % spoolCounter            => NSCOUGT
      Tau0Descript % filepos                 => ITAU0P
      Tau0Descript % array                => TAU0VAR
      Tau0Descript % array_g              => TAU0VAR_g
      Tau0Descript % field_name           = 'Tau0'
      ! jgf51.52.30: Fixing the writing of tau0 so it only occurs
      ! when tau0=-x.1
      Tau0Descript % writeFlag            = OutputTau0

      ! sponge layer fort.92
      IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
         ALLOCATE(sponge_g(NP_G))
      ENDIF
      SpongeDescript % lun                  = 92
      SpongeDescript % specifier            = ASCII
      SpongeDescript % startTimeStep           = NTCYSGE
      SpongeDescript % endTimeStep             = NTCYFGE
      SpongeDescript % outputTimeStepIncrement = NSPOOLGE
      SpongeDescript % spoolCounter            => NSCOUGS
      SpongeDescript % filepos                 => IGSP
      SpongeDescript % array                => sponge
      SpongeDescript % array_g              => sponge_g
      SpongeDescript % writeFlag            = outputsponge

      ! NodeCode output
      IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
         ALLOCATE(NodeCode_out_g(NP_G))
      ENDIF
      nodeCodeOutDescript % lun                     = 101
      if ( NOUTGE.ne.0 ) then
         nodeCodeOutDescript % specifier            = ASCII
      endif
      nodeCodeOutDescript % startTimeStep           = NTCYSGE
      nodeCodeOutDescript % endTimeStep             = NTCYFGE
      nodeCodeOutDescript % outputTimeStepIncrement = NSPOOLGE
      nodeCodeOutDescript % spoolCounter    => nodeCodeOutSpoolCounter
      nodeCodeOutDescript % filepos         => nodeCodeOutFilePos
      nodeCodeOutDescript % iarray               => NODECODE
      nodeCodeOutDescript % iarray_g             => NodeCode_out_g
      nodeCodeOutDescript % field_name           = 'NodeCode'
      nodeCodeOutDescript % file_basename        = 'nodecode'
      nodeCodeOutDescript % file_extension       = 63
      nodeCodeOutDescript % writeFlag            = outputNodeCode
      nodeCodeOutDescript % isInteger            = .true.

      ! NOFF output
      IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
         ALLOCATE(NOFF_out_g(NE_G))
      ENDIF
      noffOutDescript % lun                  = 102
      if ( NOUTGE.ne.0 ) then
         noffOutDescript % specifier            = ASCII
      endif
      noffOutDescript % startTimeStep           = NTCYSGE
      noffOutDescript % endTimeStep             = NTCYFGE
      noffOutDescript % outputTimeStepIncrement = NSPOOLGE
      noffOutDescript % spoolCounter            => noffOutSpoolCounter
      noffOutDescript % filepos                 => noffOutFilePos
      noffOutDescript % num_fd_records       = NE_G
      noffOutDescript % num_records_this     = NE
      noffOutDescript % imap                 => IMAP_EL_LG
      noffOutDescript % iarray                => NOFF
      noffOutDescript % iarray_g              => noff_out_g
      noffOutDescript % field_name           = 'NOFF'
      noffOutDescript % file_basename        = 'noff'
      noffOutDescript % file_extension       = 100
      noffOutDescript % writeFlag            = outputNOFF
      noffOutDescript % isInteger            = .true.

      ! fort.64
      VelDescript % lun                  = 64
      VelDescript % specifier            = NOUTGV
      VelDescript % startTimeStep           = NTCYSGV
      VelDescript % endTimeStep             = NTCYFGV
      VelDescript % outputTimeStepIncrement = NSPOOLGV
      VelDescript % spoolCounter            => NSCOUGV
      VelDescript % filepos                 => IGVP
      VelDescript % num_items_per_record = 2
      VelDescript % array                => UU2
      VelDescript % array2               => VV2
      VelDescript % array_g              => UU2_g
      VelDescript % array2_g             => VV2_g
      VelDescript % field_name           = 'Vel'
      VelDescript % writerFormats(1:5)   = (/ 1, 3, 4, 5, 7 /)

      ! fort.71
      IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
         ALLOCATE(RMP00_g(NSTAM_G))
      ENDIF
      NSCOUP = NSCOUM ! jgf51.51: So we don't share a spool counter.
      PrStaDescript % lun                  = 71
      PrStaDescript % specifier            = NOUTM
      PrStaDescript % num_fd_records       = NSTAM_G
      PrStaDescript % num_records_this     = NSTAM
      PrStaDescript % imap                 => IMAP_STAM_LG
      PrStaDescript % array                => RMP00
      PrStaDescript % array_g              => RMP00_g
      PrStaDescript % interped_array       => PR2
      PrStaDescript % field_name           = 'PrSta'
      PrStaDescript % isStation            = .true.
      PrStaDescript % startTimeStep           = NTCYSM
      PrStaDescript % endTimeStep             = NTCYFM
      PrStaDescript % outputTimeStepIncrement = NSPOOLM
      PrStaDescript % spoolCounter            => NSCOUP
      PrStaDescript % filepos                 => IPSTP
      PrStaDescript % elements                => NNM
      PrStaDescript % interp_fac1             => STAIM1
      PrStaDescript % interp_fac2             => STAIM2
      PrStaDescript % interp_fac3             => STAIM3
      IF (ICS.eq.1) THEN
         PrStaDescript % x_coord           => XEM ! use orig coord
         PrStaDescript % y_coord           => YEM ! use orig coord
      ELSE
         PrStaDescript % x_coord           => SLEM ! radians
         PrStaDescript % y_coord           => SFEM
      ENDIF
      if (nws.eq.0) then
         PrStaDescript % writeFlag         = .false.
      endif

      ! fort.72
      IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
         ALLOCATE(RMU00_g(NSTAM_G))
         ALLOCATE(RMV00_g(NSTAM_G))
      ENDIF
      WindVelStaDescript % lun                  = 72
      WindVelStaDescript % specifier            = NOUTM
      WindVelStaDescript % num_items_per_record = 2
      WindVelStaDescript % num_fd_records       = NSTAM_G
      WindVelStaDescript % num_records_this     = NSTAM
      WindVelStaDescript % imap                 => IMAP_STAM_LG
      WindVelStaDescript % array                => RMU00
      WindVelStaDescript % array_g              => RMU00_g
      WindVelStaDescript % array2               => RMV00
      WindVelStaDescript % array2_g             => RMV00_g
      WindVelStaDescript % interped_array       => wvnxout
      WindVelStaDescript % interped_array2      => wvnyout
      WindVelStaDescript % isStation            = .true.
      WindVelStaDescript % field_name           = 'WindVelSta'
      WindVelStaDescript % startTimeStep           = NTCYSM
      WindVelStaDescript % endTimeStep             = NTCYFM
      WindVelStaDescript % outputTimeStepIncrement = NSPOOLM
      WindVelStaDescript % spoolCounter            => NSCOUM
      WindVelStaDescript % filepos                 => IWSTP
      WindVelStaDescript % elements                => NNM
      WindVelStaDescript % interp_fac1             => STAIM1
      WindVelStaDescript % interp_fac2             => STAIM2
      WindVelStaDescript % interp_fac3             => STAIM3
      IF (ICS.eq.1) THEN
         WindVelStaDescript % x_coord           => XEM ! use orig coord
         WindVelStaDescript % y_coord           => YEM ! use orig coord
      ELSE
         WindVelStaDescript % x_coord           => SLEM ! radians
         WindVelStaDescript % y_coord           => SFEM
      ENDIF
      if (nws.eq.0) then
         WindVelStaDescript % writeFlag         = .false.
      endif

      ! fort.73
      IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
         ALLOCATE(Pr2_g(NP_G))
      ENDIF
      PrDescript % lun                  = 73
      PrDescript % specifier            = NOUTGW
      PrDescript % startTimeStep           = NTCYSGW
      PrDescript % endTimeStep             = NTCYFGW
      PrDescript % outputTimeStepIncrement = NSPOOLGW
      PrDescript % spoolCounter            => NSCOUGP
      PrDescript % filepos                 => IGPP
      PrDescript % writerFormats(1:5)   = (/ 1, 3, 4, 5, 7 /)
      PrDescript % array                => Pr2
      PrDescript % array_g              => Pr2_g
      PrDescript % field_name           = 'Pr'
      if (nws.eq.0) then
         PrDescript % writeFlag         = .false.
      endif

      ! fort.74
      IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
         ALLOCATE(WVNXOut_g(NP_G))
         ALLOCATE(WVNYOut_g(NP_G))
      ENDIF
      WindVelDescript % lun                  = 74
      WindVelDescript % specifier            = NOUTGW
      WindVelDescript % startTimeStep           = NTCYSGW
      WindVelDescript % endTimeStep             = NTCYFGW
      WindVelDescript % outputTimeStepIncrement = NSPOOLGW
      WindVelDescript % spoolCounter            => NSCOUGW
      WindVelDescript % filepos                 => IGWP
      WindVelDescript % writerFormats(1:5)   = (/ 1, 3, 4, 5, 7 /)
      WindVelDescript % num_items_per_record = 2
      WindVelDescript % array                => WVNXOut
      WindVelDescript % array_g              => WVNXOut_g
      WindVelDescript % array2               => WVNYOut
      WindVelDescript % array2_g             => WVNYOut_g
      WindVelDescript % field_name           = 'WindVel'
      if (nws.eq.0) then
         WindVelDescript % writeFlag         = .false.
      endif
      !
      ! wind Drag
      if ( (mnproc.gt.1) .and. (myproc.eq.0) ) then
         allocate(windDragOut_g(np_g))
      endif
      windDragDescript % lun                  = 173
      windDragDescript % specifier            = NOUTGW
      windDragDescript % startTimeStep           = NTCYSGW
      windDragDescript % endTimeStep             = NTCYFGW
      windDragDescript % outputTimeStepIncrement = NSPOOLGW
      windDragDescript % spoolCounter            => windDragSpoolCounter
      windDragDescript % filepos                 => windDragFilePosition
      windDragDescript % array                => windDragOut
      windDragDescript % array_g              => windDragOut_g
      windDragDescript % field_name           = 'windDrag'
      windDragDescript % file_basename        = 'windDrag'
      windDragDescript % writeFlag            = outputWindDrag
      !
      ! fort.75  !tcm v50.66.01 Added for Time Varying Bathymetry
      IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
         ALLOCATE(DP00_g(NSTAE_G))
      ENDIF
      NSCOUB = NSCOUE ! jgf51.51: So we don't share a spool counter.
      BathyStaDescript % lun                  = 75
      BathyStaDescript % specifier            = NOUTE
      BathyStaDescript % num_fd_records       = NSTAE_G
      BathyStaDescript % num_records_this     = NSTAE
      BathyStaDescript % imap                 => IMAP_STAE_LG
      BathyStaDescript % array                => DP00
      BathyStaDescript % array_g              => DP00_g
      BathyStaDescript % interped_array       => DP
      BathyStaDescript % field_name           = 'BathySta'
      BathyStaDescript % isStation            = .true.
      BathyStaDescript % startTimeStep           = NTCYSE
      BathyStaDescript % endTimeStep             = NTCYFE
      BathyStaDescript % outputTimeStepIncrement = NSPOOLE
      BathyStaDescript % spoolCounter            => NSCOUB
      BathyStaDescript % filepos                 => IBSTP
      BathyStaDescript % elements                => NNE
      BathyStaDescript % interp_fac1             => STAIE1
      BathyStaDescript % interp_fac2             => STAIE2
      BathyStaDescript % interp_fac3             => STAIE3
      IF (ICS.eq.2) THEN
         BathyStaDescript % x_coord           => SLEL
         BathyStaDescript % y_coord           => SFEL
      ELSE
         BathyStaDescript % x_coord           => XEL
         BathyStaDescript % y_coord           => YEL
      ENDIF
      if (nddt.eq.0) then
         BathyStaDescript % writeFlag         = .false.
      endif

      ! fort.76  !tcm v50.66.01  Added for Time Varying Bathymetry
      IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
         ALLOCATE(DPOUT_g(NP_G))
      ENDIF
      NSCOUGB = NSCOUGE !jgf51.51: So we don't share a spool counter.
      BathyDescript % lun                  = 76
      BathyDescript % specifier            = NOUTGE
      BathyDescript % startTimeStep           = NTCYSGE
      BathyDescript % endTimeStep             = NTCYFGE
      BathyDescript % outputTimeStepIncrement = NSPOOLGE
      BathyDescript % spoolCounter            => NSCOUGB
      BathyDescript % filepos                 => IGBP
      BathyDescript % array                => DP
      BathyDescript % array_g              => DPOut_g
      BathyDescript % field_name           = 'Bathy'
      BathyDescript % writerFormats(1:4)   = (/ 1, 3, 4, 5 /)
      if (nddt.eq.0) then
         BathyDescript % writeFlag         = .false.
      endif
      ! jgf51.21.24: Don't use sparse output or XDMF for time varying
      ! bathy, even if they were specified for fulldomain elevation file
      if ((abs(noutge).eq.4).or.(abs(noutge).eq.7)) then
         BathyDescript % specifier         = ASCII
      endif

      IF(USE_TVW)THEN
         IF((MNPROC.GT.1).AND.(MYPROC.EQ.0))THEN
             ALLOCATE(TVW_G(NP_G))
         ENDIF
         weirElevDescript % lun                    = 77
         weirElevDescript % Specifier              = NOUT_TVW
         weirELevDescript % startTimeStep          = NTCYS_TVW
         weirElevDescript % endTimeStep            = NTCYF_TVW
         weirElevDescript % outputTimeStepIncrement = NSPOOL_TVW
         weirElevDescript % spoolCounter           => NSCOU_TVW
         weirElevDescript % filepos                => IG_TVW
         weirElevDescript % num_fd_records         = NP_G
         weirElevDescript % num_records_this       = NP
         weirElevDescript % imap                   => NODES_LG
         weirElevDescript % array                  => TVW
         weirElevDescript % array_g                => TVW_G
         weirElevDescript % ConsiderWetDry         = .FALSE.
         weirElevDescript % field_name             = 'weir_dz'
         weirElevDescript % writerFormats(1:4)     = (/ 1, 3, 4, 5 /)
         weirElevDescript % alternate_value        = 0.0D0
      ENDIF

      ! fort.93
C         write(16,*) 'Allocating IceDescript' !tcm v49.64.01 added for ice
      IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
         ALLOCATE(CICEOUT_g(NP_G))
      ENDIF
      IceDescript % lun                  = 93
      IceDescript % specifier            = NOUTGW
      IceDescript % startTimeStep           = NTCYSGW
      IceDescript % endTimeStep             = NTCYFGW
      IceDescript % outputTimeStepIncrement = NSPOOLGW
      IceDescript % spoolCounter            => NSCOUGI
      IceDescript % filepos                 => IGIP
      IceDescript % writerFormats(1:5)   = (/ 1, 3, 4, 5, 7 /)
      IceDescript % array                => CICEOUT
      IceDescript % array_g              => CICEOUT_g
      IceDescript % field_name           = 'Ice'
      if (ncice.eq.0) then
         IceDescript % writeFlag         = .false.
      endif

C... v49.64.01 tcm -- added for ice stations
      ! fort.91
      IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
         ALLOCATE(RMICE00_g(NSTAM_G))
      ENDIF
      NSCOUI = NSCOUM !jgf51.51: So we don't share a spool counter.
      IceStaDescript % lun                  = 91
      IceStaDescript % specifier            = NOUTM
      IceStaDescript % num_fd_records       = NSTAM_G
      IceStaDescript % num_records_this     = NSTAM
      IceStaDescript % imap                 => IMAP_STAM_LG
      IceStaDescript % array                => RMICE00
      IceStaDescript % array_g              => RMICE00_g
      IceStaDescript % field_name           = 'IceSta'
      IceStaDescript % isStation            = .true.
      IceStaDescript % startTimeStep           = NTCYSM
      IceStaDescript % endTimeStep             = NTCYFM
      IceStaDescript % outputTimeStepIncrement = NSPOOLM
      IceStaDescript % spoolCounter            => NSCOUI
      IceStaDescript % filepos                 => IICESTP
      IceStaDescript % elements                => NNM
      IceStaDescript % interp_fac1             => STAIM1
      IceStaDescript % interp_fac2             => STAIM2
      IceStaDescript % interp_fac3             => STAIM3
      IceStaDescript % interped_array       => CICEOUT
      IF (ICS.eq.1) THEN
         IceStaDescript % x_coord           => XEM ! use orig coord
         IceStaDescript % y_coord           => YEM
      ELSE
         IceStaDescript % x_coord           => SLEM ! radians
         IceStaDescript % y_coord           => SFEM
      ENDIF
      if (ncice.eq.0) then
         IceStaDescript % writeFlag         = .false.
      endif

      ! Casey 090302: Copied these lines for the rads.64 file.
      ! tcm v50.75 removed ifdef cswan to allow for use whenever nrs=3 or nrs=4
      ! rads.64
      ! Radiation Stress
      IF ( (MNPROC.gt.1).and.(MyProc.eq.0).and.(nrs.ne.0) ) THEN
         ALLOCATE(RSNXOUT_g(NP_G))
         ALLOCATE(RSNYOUT_g(NP_G))
      ENDIF
      RSDescript % lun                  = 164
      if (nrs.eq.0) then
         RSDescript % specifier            = OFF
      else
         RSDescript % specifier            = NOUTGW
      endif
      RSDescript % startTimeStep           = NTCYSGW
      RSDescript % endTimeStep             = NTCYFGW
      RSDescript % outputTimeStepIncrement = NSPOOLGW
      RSDescript % spoolCounter            => NSCOURS
      RSDescript % filepos                 => IGRadS
      RSDescript % writerFormats(1:5)   = (/ 1, 3, 4, 5, 7 /)
      RSDescript % num_items_per_record = 2
      RSDescript % array                => RSNXOUT
      RSDescript % array_g              => RSNXOUT_g
      RSDescript % array2               => RSNYOUT
      RSDescript % array2_g             => RSNYOUT_g
      RSDescript % alternate_value      = -99999.0
      RSDescript % considerWetDry       = .true.
      RSDescript % field_name           = 'RadStress'
      RSDescript % file_basename            = 'rads'
      RSDescript % file_extension       = 64
      if (nrs.eq.0) then
         RSDescript % writeFlag         = .false.
      endif
#ifndef CSWAN
      ! jgf52.30.02: Turn off rads.64 if compiled without SWAN but
      ! NRS erroneously set to 3.
      if (nrs.eq.3) then
         RSDescript % writeFlag         = .false.
      endif
#endif

      ! fort.81
      IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
         ALLOCATE(CC00_g(NP_G))
      ENDIF
      ConcStaDescript % lun                  = 81
      ConcStaDescript % specifier            = NOUTC
      ConcStaDescript % num_fd_records       = NSTAC_G
      ConcStaDescript % num_records_this     = NSTAC
      ConcStaDescript % imap                 => IMAP_STAC_LG
      ConcStaDescript % array                => CC00
      ConcStaDescript % array_g              => CC00_g
      ConcStaDescript % field_name           = 'ConcSta'
      ConcStaDescript % isStation            = .true.
      ConcStaDescript % divideByDepth        = .true.
      ConcStaDescript % interped_array       => CH1
      IF (ICS.eq.1) THEN
         ConcStaDescript % x_coord           => XEC ! use orig coord
         ConcStaDescript % y_coord           => YEC
      ELSE
         ConcStaDescript % x_coord           => SLEC ! radians
         ConcStaDescript % y_coord           => SFEC
      ENDIF
      if (im.ne.10) then
         ConcStaDescript % writeFlag         = .false.
      endif

      ! fort.83
      IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
         ALLOCATE(CH1_g(NP_G))
      ENDIF
      ConcDescript % lun                  = 83
      ConcDescript % specifier            = NOUTGC
      ConcDescript % startTimeStep           = NTCYSGC
      ConcDescript % endTimeStep             = NTCYFGC
      ConcDescript % outputTimeStepIncrement = NSPOOLGC
      ConcDescript % spoolCounter            => NSCOUGC
      ConcDescript % filepos                 => IGCP
      ConcDescript % array                => CH1
      ConcDescript % array_g              => CH1_g
      ConcDescript % field_name           = 'Conc'
      ConcDescript % divideByDepth        = .true.
      if (im.ne.10) then
         ConcDescript % writeFlag         = .false.
      endif

C     maxele.63
      IF ( mnproc.gt.1 ) THEN
         ALLOCATE(ETAMAX_g(NP_G))
         ALLOCATE (ETAMAX_Time_G(NP_G))
         ETAMAX_g(:) = -99999.d0
         ETAMAX_Time_G(:) = -99999.d0
      ENDIF
      EtaMaxDescript % lun                 = 311
      EtaMaxDescript % specifier            = NOUTGE
      EtaMaxDescript % num_items_per_record = 2   !overloading this value to account for time stamps along with values
      EtaMaxDescript % array                => ETAMAX
      EtaMaxDescript % array_g              => ETAMAX_G
      EtaMaxDescript % array2               => ETAMAX_Time
      EtaMaxDescript % array2_g             => ETAMAX_Time_G
      EtaMaxDescript % field_name           = 'EtaMax'
      EtaMaxDescript % file_basename        = 'maxele'
      EtaMaxDescript % file_extension       = 63
      EtaMaxDescript % writeFlag            = .false.
      EtaMaxDescript % minmax_timestamp     = .true.
      EtaMaxDescript % readMaxMin           = .true.
      if ( EtaMaxDescript % specifier .eq. XDMF ) then
         EtaMaxDescript % specifier = ASCII
      endif
      ETAMAX(:) = -99999.d0
      ETAMAX_Time(:) = -99999.d0

C     maxvel.63
      IF (mnproc.gt.1)  THEN
         ALLOCATE(UMAX_g(NP_G))
         ALLOCATE (UMAX_Time_G(NP_G))
         UMAX_g(:) = 0.d0
         UMAX_Time_G(:) = -99999.d0
      ENDIF
      UMaxDescript % lun                    = 312
      UMaxDescript % num_items_per_record   = 2
      UMaxDescript % specifier              = NOUTGV
      UMaxDescript % array                  => UMAX
      UMaxDescript % array_g                => UMAX_G
      UMaxDescript % array2                 => UMAX_Time
      UMaxDescript % array2_g               => UMAX_Time_G
      UMaxDescript % field_name             = 'UMax'
      UMaxDescript % file_basename          = 'maxvel'
      UMaxDescript % file_extension         = 63
      UMaxDescript % writeFlag              = .false.
      UMaxDescript % minmax_timestamp       = .true.
      UMaxDescript % readMaxMin             = .true.
      if ( UMaxDescript % specifier .eq. XDMF ) then
         UMaxDescript % specifier = ASCII
      endif
      UMAX(:) = 0.d0
      UMAX_Time(:) = -99999.d0

C     prmin.63
      IF (mnproc.gt.1) THEN
         ALLOCATE(PRMIN_g(NP_G))
         ALLOCATE (PRMIN_Time_G(NP_G))
         PRMIN_g(:) = 99999.d0
         PRMIN_Time_G(:) = -99999.d0
      ENDIF
      PrMinDescript % lun                   = 313
      PrMinDescript % num_items_per_record  = 2
      PrMinDescript % specifier             = NOUTGW
      PrMinDescript % array                 => PRMIN
      PrMinDescript % array_g               => PRMIN_G
      PrMinDescript % array2                => PRMIN_Time
      PrMinDescript % array2_g              => PRMIN_Time_G
      PrMinDescript % field_name            = 'PrMin'
      PrMinDescript % file_basename         = 'minpr'
      PrMinDescript % file_extension        = 63
      PrMinDescript % writeFlag             = .false.
      PrMinDescript % minmax_timestamp      = .true.
      if (nws.ne.0) then
         PrMinDescript % readMaxMin         = .true.
      endif
      if ( PrMinDescript % specifier .eq. XDMF ) then
         PrMinDescript % specifier = ASCII
      endif
      PRMIN(:) = 99999.d0
      PRMIN_Time(:) = -99999.d0

C     maxwvel.63
      IF (mnproc.gt.1) THEN
         ALLOCATE(WVNOUTMAX_g(NP_G))
         ALLOCATE (WVNOUTMAX_Time_G(NP_G))
         WVNOUTMAX_g(:) = 0.d0
         WVNOUTMAX_Time_G(:) = -99999.d0
      ENDIF
      WVMaxDescript % lun                   = 314
      WVMaxDescript % specifier             = NOUTGW
      WVMaxDescript % num_items_per_record  = 2
      WVMaxDescript % array                 => WVNOUTMAX
      WVMaxDescript % array_g               => WVNOUTMAX_G
      WVMaxDescript % array2                => WVNOUTMAX_Time
      WVMaxDescript % array2_g              => WVNOUTMAX_Time_G
      WVMaxDescript % field_name            = 'WVMax'
      WVMaxDescript % file_basename         = 'maxwvel'
      WVMaxDescript % file_extension        = 63
      WVMaxDescript % writeFlag             = .false.
      WVMaxDescript % minmax_timestamp      = .true.
      if (nws.ne.0) then
         WVMaxDescript % readMaxMin         = .true.
      endif
      if ( WVMaxDescript % specifier .eq. XDMF ) then
         WVMaxDescript % specifier = ASCII
      endif
      WVNOUTMAX(:) = 0.d0
      WVNOUTMAX_Time(:) = -99999.d0
C
C     maxrs.63
      IF (mnproc.gt.1) THEN
         ALLOCATE(RSNMAX_g(NP_G))
         ALLOCATE (RSNMAX_Time_G(NP_G))
         RSNMAX_g(:) = 0.d0
         RSNMAX_Time_G(:) = -99999.d0
      ENDIF
      RSMaxDescript % lun                   = 315
      RSMaxDescript % specifier             = NOUTGW
      RSMaxDescript % num_items_per_record  = 2
      RSMaxDescript % array2                => RSNMAX_Time
      RSMaxDescript % array2_g              => RSNMAX_Time_G
      RSMaxDescript % array                 => RSNMAX
      RSMaxDescript % array_g               => RSNMAX_G
      RSMaxDescript % field_name            = 'RSMax'
      RSMaxDescript % file_basename         = 'maxrs'
      RSMaxDescript % file_extension        = 63
      RSMaxDescript % writeFlag             = .false.
      RSMaxDescript % minmax_timestamp      = .true.
      if (nrs.ne.0) then
         RSMaxDescript % readMaxMin         = .true.
      endif
      if ( RSMaxDescript % specifier .eq. XDMF ) then
         RSMaxDescript % specifier = ASCII
      endif
      RSNMAX(:) = 0.d0
      RSNMAX_Time(:) = -99999.d0

!     ESLNodes.63
      IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
         ALLOCATE(ESLONOFF_g(NP_G))
         ESLONOFF_g(:) = 0.d0
      ENDIF
      ESLDescript % lun                  = 323
      ESLDescript % specifier            = ASCII
      ESLDescript % array                => ESLONOFF
      ESLDescript % array_g              => ESLONOFF_G
      ESLDescript % field_name           = 'ESLActive'
      ESLDescript % file_basename        = 'ESLNodes'
      ESLDescript % file_extension       = 63
      ESLDescript % writeFlag            = .false.
      ESLONOFF(:) = 0.d0
      !
      !  D E T A I L E D    I N U N D A T I O N    O U T P U T
      !
      ! inundationtime.63 (works like a min/max file)
      if (mnproc.gt.1) then
         allocate(inundationTime_g(NP_G))
         allocate(inundationTime_onset_g(NP_G))
         inundationTime_g(:) = 0.d0
         inundationTime_onset_g(:) = -99999.d0
      endif
      InundationTimeDescript % lun                  = 400
      InundationTimeDescript % specifier            = NOUTGE
      InundationTimeDescript % num_items_per_record = 2
      InundationTimeDescript % array                => inundationTime
      InundationTimeDescript % array_g              => inundationTime_g
      InundationTimeDescript % array2               => inundationTime_onset
      InundationTimeDescript % array2_g             => inundationTime_onset_g
      InundationTimeDescript % field_name           = 'inundationTime'
      InundationTimeDescript % file_basename        = 'inundationtime'
      InundationTimeDescript % file_extension       = 63
      InundationTimeDescript % writeFlag            = .false.
      InundationTimeDescript % minmax_timestamp     = .true.
      InundationTimeDescript % readMaxMin           = inundationOutput
      if ( InundationTimeDescript % specifier .eq. XDMF ) then
         InundationTimeDescript % specifier = ASCII
      endif
      !
      ! maxinundepth.63 (works like a min/max file)
      if (mnproc.gt.1) then
         allocate(MaxInunDepth_g(NP_G))
         allocate(MaxInunDepth_time_g(NP_G))
         MaxInunDepth_g(:) = -99999.d0
         MaxInunDepth_time_g(:) = -99999.d0
      endif
      MaxInunDepthDescript % lun                  = 401
      MaxInunDepthDescript % specifier            = NOUTGE
      MaxInunDepthDescript % num_items_per_record = 2
      MaxInunDepthDescript % array                => maxInunDepth
      MaxInunDepthDescript % array_g              => maxInunDepth_g
      MaxInunDepthDescript % array2               => maxInunDepth_time
      MaxInunDepthDescript % array2_g             => maxInunDepth_time_g
      MaxInunDepthDescript % field_name           = 'maxInunDepth'
      MaxInunDepthDescript % file_basename        = 'maxinundepth'
      MaxInunDepthDescript % file_extension       = 63
      MaxInunDepthDescript % writeFlag            = .false.
      MaxInunDepthDescript % minmax_timestamp     = .true.
      MaxInunDepthDescript % readMaxMin           = inundationOutput
      if ( MaxInunDepthDescript % specifier .eq. XDMF ) then
         MaxInunDepthDescript % specifier = ASCII
      endif
      !
      ! endrisinginun.63 (1 if water is rising at end of simulation)
      IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
         ALLOCATE(endRisingInun_g(NP_G))
         endRisingInun_g(:) = 0
      ENDIF
      EndRisingInunDescript % lun                  = 403
      EndRisingInunDescript % specifier            = NOUTGE
      EndRisingInunDescript % filepos              => fauxFilePosition
      EndRisingInunDescript % iarray               => endRisingInun
      EndRisingInunDescript % iarray_g             => endRisingInun_g
      EndRisingInunDescript % field_name           = 'endRisingInun'
      EndRisingInunDescript % file_basename        = 'endrisinginun'
      EndRisingInunDescript % file_extension       = 63
      EndRisingInunDescript % isInteger            = .true.
      EndRisingInunDescript % writeFlag            = inundationOutput
      ! write_flag set to .false. immediately after file initialization
      if ( EndRisingInunDescript % specifier .eq. XDMF ) then
         EndRisingInunDescript % specifier = ASCII
      endif
      ! jgf52.08.03: We write this as a regular file, rather than
      ! a min/max file, and we never want to append, we always want
      ! to overwrite.
      if ( EndRisingInunDescript % specifier .gt.0 ) then
         EndRisingInunDescript % specifier = - EndRisingInunDescript % specifier
      endif
      !
      ! everdried.63 (works like a min/max file)
      if (mnproc.gt.1) then
         allocate(EverDried_g(NP_G))
         allocate(EverDried_time_g(NP_G))
         EverDried_g(:) = 1.d0
         EverDried_time_g(:) = -99999.d0
      endif
      EverDriedDescript % lun                  = 404
      EverDriedDescript % specifier            = NOUTGE
      EverDriedDescript % num_items_per_record = 2
      EverDriedDescript % initial_value        = 0.d0
      EverDriedDescript % array                => EverDried
      EverDriedDescript % array_g              => EverDried_g
      EverDriedDescript % array2               => EverDried_time
      EverDriedDescript % array2_g             => EverDried_time_g
      EverDriedDescript % field_name           = 'everDried'
      EverDriedDescript % file_basename        = 'everdried'
      EverDriedDescript % file_extension       = 63
      EverDriedDescript % writeFlag            = .false.
      EverDriedDescript % minmax_timestamp     = .true.
      EverDriedDescript % readMaxMin           = inundationOutput
      if ( EverDriedDescript % specifier .eq. XDMF ) then
         EverDriedDescript % specifier = ASCII
      endif

#ifdef CSWAN
      CALL initSwanOutput() !...couple2swan.F --> performs same process as above
#endif
      !
      ! create file names
      do i=1, numOutputDescript2D
         ! set the file extension to be the same as the logical unit number (lun)
         ! by default if it doesn't have a special value already set
         ! in the code above
         if ( ptr(i) % descript2D % file_extension .eq. -99999 ) then
            ptr(i) % descript2D % file_extension = ptr(i) % descript2D % lun
         endif
         ! write the integer file extension to a character string
         write(extString,'(i0)') ptr(i) % descript2D % file_extension
         !
         ! construct the file name based on the type of run that is occurring
         ! (serial or parallel, write to local files, etc)
         !
         ! jgf: netCDF and XDMF files are always fulldomain...
         !
         if (((mnproc.gt.1).and.(WRITE_LOCAL_FILES.eqv..false.)).or.
     &        (ptr(i) % descript2D % specifier.eq.NETCDF3).or.
     &        (ptr(i) % descript2D % specifier.eq.NETCDF4).or.
     &        (ptr(i) % descript2D % specifier.eq.XDMF) ) then
            ptr(i) % descript2D % file_name =
     &         trim(globaldir) // '/' //
     &         trim(ptr(i) % descript2D % file_basename) //
     &         '.' // trim(extString)
         else
            ptr(i) % descript2D % file_name =
     &         trim(localdir) // '/' //
     &         trim(ptr(i) % descript2D % file_basename) //
     &         '.' // trim(extString)
         endif

      end do
      !
      ! determine whether writer processors will be used for each type
      ! of data
      if (mnwproc.gt.0) then
         do i=1, numOutputDescript2D
            ! if the chosen output format is supported by
            ! writer processors, then use writer processors for
            ! this output type
            if ( ptr(i) % descript2D % writeFlag .eqv. .true.) then
               if ( any(ptr(i) % descript2D % writerFormats(:) .eq.
     &               abs(ptr(i) % descript2D % specifier) )  ) then
                  ptr(i) % descript2D % useWriter = .true.
                  write(scratchMessage,'(a,a,a)') 'Output file ',
     &               trim(ptr(i) % descript2D % file_name),
     &             ' will be written via dedicated writer processor(s).'
                  call allMessage(INFO,scratchMessage)
               endif
            endif
         end do
      endif
      !
      ! jgf49.43.11: if netcdf is not required for a particular output file,
      ! the associated subroutine call will do nothing and just return
      if (myProc.eq.0) then
         do i=1, numOutputDescript2D
            nerr = .false.
            if ( ptr(i) % descript2D % writeFlag .eqv. .true.) then
               select case(abs(ptr(i) % descript2D % specifier))
               case(NETCDF3, NETCDF4)
#ifdef ADCNETCDF
                  call initNetCDFOutputFile(ptr(i) % descript2D, nerr)
                  if (nerr.eqv..true.) then
                     call allMessage(ERROR,
     &             'There was an issue with NetCDF file initialization.')
                     call terminate()
                  endif
                  ptr(i) % descript2D % initialized = .true.
#endif
               ! jgf51.21.24: Initialize XDMF output files
               case(XDMF)
#ifdef ADCXDMF
                  if (ptr(i) % descript2D % useWriter.eqv..false.) then
                     call readControlFile(controlFileName,.false.)
                     call initOutputXDMF(ptr(i) % descript2D)
                     call writeControlXDMF(ptr(i) % descript2D % xdmfMD
     &                                     % xdmfFortranObj)
                     ptr(i) % descript2D % initialized = .true.
#ifdef CMPI
                  else
                     ! make writer read mesh; the writing of the
                     ! mesh to the XDMF output file will be performed
                     ! when actual output data are written to the file
                     call sendInitWriterXDMF(ptr(i) % descript2D)
#endif
                  endif
#endif
               case default
                  ! special initialization not required for ASCII or SPARSE_ASCII
               end select
            endif
         end do
      endif
      !
      ! jgf52.08.03: If inundationOutput is .true. and endrisinginun.63
      ! has been initialized, turn it off so it doesn't get written
      ! during time stepping.
      EndRisingInunDescript % writeFlag = .false.
      !
      !
      !  R E A D   M I N   /   M A X   F I L E S   A T   H O T S T A R T
      !
      ! jgf48.4635 If the maxele.63 (etc) files are present, we need to
      ! load them up, to preserve high water marks across hotstarts.
      ! If they are not present, thats okay too.
      !
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !  Begin TCM v51.20.01 Mods for MaxMin Files to local Variables
      !      Including reading time stamp for max/min when available.
      ! Note:  Global initalizes the TimeStamp Values to 0.0d0
      !   If timestamp values are not found, then they will return
      !   with an initial value of -99999.d0 to begin the hstart runs.
      !   This includes the case of the max/min files not being present.
      !
      ! jgf52.08.03: Moved the reading of min/max files from hotstart()
      ! subroutine to here.
      !
      ! jgf52.08.03: Look for ascii-formatted files and load those if
      ! found; otherwise try to load netcdf-formatted files using new
      ! netcdf reading capability.
      !
      ! jgf52.08.20: Check the timestamp on the min/max data; only load
      ! min/max data at hotstart if the timestamp is earlier than the
      ! current adcirc time.
      if (ihot.ne.0) then
         do i=1, numOutputDescript2D
            fileFound = .false.
            if ( ptr(i) % descript2D % readMaxMin.eqv..true.) then
               call logMessage(INFO,'Attempting to read ' //
     &            trim(ptr(i) % descript2D % file_name) //
     &            ' to continue the min/max record.')
               call readAndMapToSubdomainMaxMin(ptr(i) % descript2D,
     &            timeloc, fileFound)
#ifdef ADCNETCDF
               if ((fileFound.eqv..false.).and.
     &              (netcdf_avail.eqv..true.) ) then
                  ! check to see if the file exists before trying to
                  ! initialize it
                  inquire(file=trim(ptr(i) % descript2D % file_name)//
     &             '.nc',exist=fileFound)
                  if (fileFound.eqv..true.) then
                      nerr = .false.
                     if (myproc.eq.0) then
C                      ! initialize it if it wasn't initialized above
                       if (ptr(i) % descript2D % initialized.eqv..false.) then
                           call initNetCDFOutputFile(ptr(i) % descript2D, nerr)
                           if (nerr.eqv..true.) then
                              call allMessage(ERROR,
     &                        "Failed to initialize "//
     &                        trim(ptr(i)%descript2D%file_name)//".")
                              call terminate()
                           endif
                        endif
                     endif
                     ptr(i) % descript2D % initialized = .true.
                     ! if there wasn't an initialization error, read in
                     ! the min/max values for use in this simulation
                     call readAndMapToSubdomainMaxMinNetCDF(
     &                  ptr(i) % descript2D, timeloc)
                  else
                     call logMessage(INFO,'The file "'//
     &                   trim(ptr(i) % descript2D % file_name)//
     &                   '" was not found.')
                  endif
               endif
#endif
            endif
         enddo
         IF (LoadEleSlopeLim) THEN
            CALL readAndMapToSubdomainMaxMin(ESLDescript,timeloc,fileFound)
            ! must set the elemental slope limiter active values as well
            !
            ! jgf51.46: Really should only turn on the slope limiter if
            ! the slope limiter was turned on in the run that produced
            ! these values (as opposed to just recording the places
            ! where the slope was exceeded). However, this isn't really
            ! practical because it would require knowledge of the
            ! elemental_slope_limiter nodal attribute from the
            ! previous simulation. We don't have access to that.
            !
            ! Also, we don't want to turn the slope limiter on in places
            ! where the analyst does not want it to be turned on, and only
            ! wants to record the places where the slope is exceeded. So
            ! I am adding a check for positive values of the nodal
            ! attribute, and only turning on the slope limiter there.
            DO I = 1,NP
               if ( (ESLONOFF(I).eq.1.d0).and.
     &         (elemental_slope_limiter_grad_max(i).ge.0.d0) ) then
                  elemental_slope_limiter_active(I) = .true.
               else
                  elemental_slope_limiter_active(I) = .false.
               endif
            ENDDO
         ENDIF
      endif ! ihot.ne.0
      !
      ! End  TCM v51.20.01 Mods for MaxMin Files to Local Variables
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C----------------------------------------------------------------------
      end subroutine initOutput2D
C----------------------------------------------------------------------

C----------------------------------------------------------------------
C       S U B R O U T I N E   W R I T E   O U T P U T   2 D
C----------------------------------------------------------------------
C
C     R.L. 8/22/05 Subroutine to write primary 2D model output not
C     including hotstart and harmonic analysis.
C
C----------------------------------------------------------------------
      SUBROUTINE writeOutput2D(IT,TimeLoc)
      USE SIZES, ONLY : SZ, INPUTDIR, NBYTE, MNWPROC, MYPROC, MNPROC,
     &                  GLOBALDIR, OFF
      USE GLOBAL
      USE MESH, ONLY : NP, NE, DP, NM, ICS
      USE GLOBAL_IO, ONLY: storeOne, storeTwo,
     &    packOne, unpackOne, packTwo, unpackTwo
      USE NodalAttributes, ONLY : LoadEleSlopeLim
#ifdef CMPI
      USE WRITER, ONLY : sendDataToWriter, flush_writers, writer_init
#endif
      IMPLICIT NONE
      INTEGER, intent(in) :: IT
      REAL(8), intent(in) :: TimeLoc
      INTEGER :: i,j

      call setMessageSource("writeOutput2D")
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      do i=1,numOutputDescript2D
         !
         ! if these data should not be written at all, due to circumstances
         ! particular to each data type, then skip this
         ! output data type and go to the next one
         if ( ptr(i) % descript2D % writeFlag.eqv..false. ) then
            cycle
         endif
         !
         ! if this output was not selected in the input file, go to the next one
         if (ptr(i) % descript2D % specifier.eq.OFF) then
            cycle
         endif
         !
         ! if we are in the time step window where output should be produced,
         ! increment the spool counter for this output type
         if ( (it.gt.ptr(i) % descript2D % startTimeStep).and.
     &        (it.le.ptr(i) % descript2D % endTimeStep) ) then
            ptr(i) % descript2D % spoolCounter =
     &          ptr(i) % descript2D % spoolCounter + 1
         endif
         !
         ! if the spool counter for this output type does not match the
         ! output time step increment for this output type, then go to
         ! the next data type
         if (ptr(i) % descript2D % spoolCounter.ne.
     &       ptr(i) % descript2D % outputTimeStepIncrement) then
            cycle
         endif
         !
         ! if this is station data, compute the solution at the station
         ! location via spatial interpolation using the three surrrounding
         ! nodes
         if ( ptr(i) % descript2D % isStation .eqv. .true. ) then
            call stationDataInterp(ptr(i) % descript2D)
         endif
         !
         ! write out the array based on the specified file format and the
         ! type of data, using writer processors if appropriate
         !
         !  W R I T E   U S I N G   D E D I C A T E D   W R I T E R
         !
         if ( ptr(i) % descript2D % useWriter .eqv. .true. ) then
            if ( ptr(i) % descript2D % num_items_per_record .eq. 1 ) then
#ifdef CMPI
               call sendDataToWriter(ptr(i) % descript2D, timeLoc, it,
     &               packOne)
            else
               call sendDataToWriter(ptr(i) % descript2D, timeLoc, it,
     &               packTwo)
#endif
            endif
         else
            !
            !  W R I T E   U S I N G   P R O C E S S O R   0
            !
            if ( ptr(i) % descript2D % num_items_per_record .eq. 1 ) then
                  call writeOutArray(TimeLoc, IT, ptr(i) % descript2D,
     &                               packOne, unpackOne)
            else
               call writeOutArray(TimeLoc, IT, ptr(i) % descript2D,
     &                            packTwo, unpackTwo)
            endif
         endif
         !
         ! set the spool counter back to zero, now that we've written these data
         ptr(i) % descript2D % spoolCounter = 0
         !
         ! special cases:
         !
         ! sponge only outputs once
         if ( ptr(i) % descript2D % lun .eq. 92 ) then
            ptr(i) % descript2D % writeFlag = .false.
         endif

      end do

#ifdef CMPI
Csb 02/09/2007
C...  Communicate writer processors and let them write out what
C...  they have.
      !write(6,'(a)') 'DEBUG: Calling flush_writers()'
      CALL FLUSH_WRITERS()
#endif

Csb 11/11/2006
C...   OUTPUT MAXIMUM/MINIMUM GLOBAL DATA AT THE LAST TIME STEP.
C...
C     jgf48.4636 Updated to newer data structure, so that we can
C     compare with values from hot start, so that mins and maxes
C     can be preserved across hot starts.
C     tcmv51.20.01 Changed from pack/unpack one to pack/unpack two
C       as max/min are now initialized from file during hotstart
C       for all processors, not just myproc=0.
      if (it.eq.nt) then
         CALL writeOutArrayMinMax(311, TimeLoc, IT, 2, EtaMaxDescript,
     &        packTwo, unpackTwo)
         CALL writeOutArrayMinMax(312, TimeLoc, IT, 2, UMaxDescript,
     &        packTwo, unpackTwo)
         IF (NWS.NE.0) THEN
            CALL writeOutArrayMinMax(313, TimeLoc, IT, 2, PrMinDescript,
     &           packTwo, unpackTwo)
            CALL writeOutArrayMinMax(314, TimeLoc, IT, 2, WVMaxDescript,
     &           packTwo, unpackTwo)
            IF (NRS.NE.0) THEN
               CALL writeOutArrayMinMax(315, TimeLoc, IT, 2, RSMaxDescript,
     &              packTwo, unpackTwo)
            ENDIF
         ENDIF
         !
         ! OUTPUT THE NODES WHERE THE ELEMENTAL SLOPE LIMITER HAS BEEN
         ! TURNED ON DURING THE RUN
         IF (LoadEleSlopeLim) THEN
            CALL writeOutArrayMinMax(323, TimeLoc, IT, 1, ESLDescript,
     &               packOne, unpackOne)  !no time stamps so only one pack/unpack
         ENDIF

C        OUTPUT THE SWAN MAXIMUM FILES
#ifdef CSWAN
         IF(Swan_OutputHS)
     &      CALL writeOutArrayMinMax(316,TimeLoc,IT,1,SwanHSMaxDescript,packOne,unpackOne)
         IF(Swan_OutputDIR)
     &      CALL writeOutArrayMinMax(317,TimeLoc,IT,1,SwanDIRMaxDescript,packOne,unpackOne)
         IF(Swan_OutputTM01)
     &      CALL writeOutArrayMinMax(318,TimeLoc,IT,1,SwanTM01MaxDescript,packOne,unpackOne)
         IF(Swan_OutputTM02)
     &      CALL writeOutArrayMinMax(321,TimeLoc,IT,1,SwanTM02MaxDescript,packOne,unpackOne)
         IF(Swan_OutputTMM10)
     &      CALL writeOutArrayMinMax(322,TimeLoc,IT,1,SwanTMM10MaxDescript,packOne,unpackOne)
         IF(Swan_OutputTPS)
     &      CALL writeOutArrayMinMax(319,TimeLoc,IT,1,SwanTPSMaxDescript,packOne,unpackOne)
         IF(Swan_OutputWind)
     &      CALL writeOutArrayMinMax(320,TimeLoc,IT,1,SwanWindMaxDescript,packTwo,unpackTwo)
#endif
         ! Write the inundationtime.63, maxinundepth.63, everDried.63, and
         ! endrisinginun.63 files at the end of the simulation.
         if (inundationOutput.eqv..true.) then
            call writeOutArrayMinMax(400, TimeLoc, IT, 2,
     &         InundationTimeDescript, packTwo, unpackTwo)
            call writeOutArrayMinMax(401, TimeLoc, IT, 2,
     &         MaxInunDepthDescript, packTwo, unpackTwo)
            call writeOutArray(timeLoc, it,
     &          EndRisingInunDescript, packOne, unpackOne)
            call writeOutArrayMinMax(404, TimeLoc, IT, 2,
     &         everDriedDescript, packTwo, unpackTwo)
         endif

      ENDIF

#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN
 2120 FORMAT(2X,1pE20.10E3,5X,I10)
 2453 FORMAT(2x, i8, 2x, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3)
 2454 FORMAT(2X,I8,2(2X,1pE20.10E3))
C-----------------------------------------------------------------------
      END SUBROUTINE writeOutput2D
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C      S U B R O U T I N E   S T A T I O N   D A T A   I N T E R P
C-----------------------------------------------------------------------
C     jgf51.21.24: Spatially interpolate the solution at a station using
C     the solution at the three surrounding nodes; performs the
C     the interpolation on both components for output data that contains
C     two components.
C-----------------------------------------------------------------------
      subroutine stationDataInterp(descript)
      use global
      implicit none
      type(OutputDataDescript_t), intent(inout) :: descript

      call setMessageSource("stationDataInterp")
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      call stationArrayInterp(descript, descript % interped_array,
     &       descript % array)
      if ( descript % num_items_per_record .eq. 2 ) then
         call stationArrayInterp(descript, descript % interped_array2,
     &         descript % array2)
      endif

#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

C-----------------------------------------------------------------------
      end subroutine stationDataInterp
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C      S U B R O U T I N E   S T A T I O N   A R R A Y   I N T E R P
C-----------------------------------------------------------------------
C     jgf51.21.24: Performs station interpolations on a single array.
C-----------------------------------------------------------------------
      subroutine stationArrayInterp(descript, fromArray, toArray)
      use sizes, only : sz
      use global, only : nodecode, OutputDataDescript_t, IFNLFA, ETA2, NOFF,
     &                   allMessage, DEBUG, setMessageSource, unsetMessageSource
      use mesh, only : dp, nm
      implicit none
      type(OutputDataDescript_t), intent(inout) :: descript
      real(sz), intent(in) :: fromArray(:)
      real(sz), intent(out) :: toArray(:)
      integer :: e
      real(sz) :: d1, d2, d3
      real(sz) :: H2N1, H2N2, H2N3
      integer :: ncele
      integer :: i, j

      call setMessageSource("stationArrayInterp")
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      do j=1, descript % num_records_this
         e = descript % elements(j)
         d1 = fromArray(nm(e,1))
         d2 = fromArray(nm(e,2))
         d3 = fromArray(nm(e,3))
         !
         ! if the value should be divided by the total depth
         if ( descript % divideByDepth .eqv. .true.) then
            H2N1=DP(NM(e,1))+IFNLFA*ETA2(nm(e,1))
            H2N2=DP(NM(e,2))+IFNLFA*ETA2(NM(e,2))
            H2N3=DP(NM(e,3))+IFNLFA*ETA2(NM(e,3))
            d1=d1/H2N1
            d2=d2/H2N2
            d3=d3/H2N3
         endif
         !
         ! perform spatial interpolation
         toArray(j) = d1 * descript % interp_fac1(j) +
     &                d2 * descript % interp_fac2(j) +
     &                d3 * descript % interp_fac3(j)
         !
         ! if the station output should have a special value in a
         ! dry element, apply the special value if necessary
         ! jgf52.04.03: Removed NOFF from consideration when
         ! determining wet/dry state for station interpolation purposes.
         if (descript % considerWetDry .eqv..true. ) then
            ncele = nodecode(nm(e,1)) * nodecode(nm(e,2)) *
     &             nodecode(nm(e,3))
            if (ncele.eq.0) then
               toArray(j) = descript % alternate_value
            endif
         endif

      end do

#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

C----------------------------------------------------------------------
      end subroutine stationArrayInterp
C----------------------------------------------------------------------



C----------------------------------------------------------------------
C       S U B R O U T I N E
C          W R I T E  H A R M O N I C   A N A L Y S I S   O U T P U T
C----------------------------------------------------------------------
C     jgf49.44: Subroutine to write harmonic analysis output files
C     in serial or in parallel using globalio.
C----------------------------------------------------------------------
      SUBROUTINE writeHarmonicAnalysisOutput(ITIME)
      USE SIZES, ONLY : SZ, MNHARF, LOCALDIR, GLOBALDIR,
     &                  WRITE_LOCAL_FILES, MNPROC, MYPROC, MNP,
     &                  WRITE_LOCAL_HARM_FILES
C RJW added MNPROC and MYPROC
C removed from global  myProc, MNPROC,
      USE GLOBAL, ONLY : OutputDataDescript_t, NSTAE_G,
     &                   NOUTE, NSTAE, IMAP_STAE_LG, SLEL,
     &                   XEL, SFEL, YEL, NSTAV_G, NSTAV,
     &                   IMAP_STAV_LG, SLEV, XEV, SFEV,
     &                   YEV, screenMessage, DEBUG,
     &                   setMessageSource, unsetMessageSource,
     &                   NE_G, NP_G, NODES_LG
      USE MESH, ONLY : NP, NE, ICS
      USE GLOBAL_IO, ONLY:
     &    open_gbl_file, collectFullDomainArray,
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
      LOGICAL, SAVE :: FirstCall = .true.
C
      call setMessageSource("writeHarmonicAnalysisOutput")
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")

#endif
      ! JUST RETURN if harmonic analysis was not specified
      ! tcm v49.64.02 -added itime le ithas to be
      ! consistent with what is in solveHarmonicAnalysis
      IF ((IHARIND.eq.0).or.(ITIME.LE.ITHAS)) THEN
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
         RETURN
      ENDIF
C
C     This subroutine should only be called once, at the end of a run, so
C     a check of FirstCall is not necessary. However, we'll still do it, in
C     case the subroutine is ever called more than once in some future
C     scenario.
      IF (FirstCall) THEN
         !
         ! Chris S. 05/08/14 Add explicit filename to each magnitude
         ! (only do for u magnitude for velocity) descriptor
         ! instead of building fn within writeOutHarmonicArrays() subroutine
         IF (NHASE.ne.0) THEN
            ! fort.51 magnitude
            HAElevStaMagDescript % specifier            = NHASE
            HAElevStaMagDescript % initial_value        = 0.0
            HAElevStaMagDescript % num_items_per_record = MNHARF
            HAElevStaMagDescript % num_fd_records       = NSTAE_G
            HAElevStaMagDescript % num_records_this     = NSTAE
            HAElevStaMagDescript % imap                 => IMAP_STAE_LG
            HAElevStaMagDescript % array2D              => EMAG
            HAElevStaMagDescript % ConsiderWetDry       = .FALSE.
            HAElevStaMagDescript % file_name            = 'fort.51'
            HAElevStaMagDescript % field_name           = 'HAElevStaMag'
            HAElevStaMagDescript % array2D_g            => null()
            IF ((MNPROC.gt.1).AND.(WRITE_LOCAL_HARM_FILES.eqv..false.)) THEN
               IF (MyProc.eq.0) THEN
                  ALLOCATE(EMAG_g(MNHARF,NSTAE_G))
               ENDIF
               HAElevStaMagDescript % array2D_g         => EMAG_g
            ENDIF
            IF (ICS.eq.2) THEN
               HAElevStaMagDescript % x_coord           => SLEL
               HAElevStaMagDescript % y_coord           => SFEL
            ELSE
               HAElevStaMagDescript % x_coord           => XEL
               HAElevStaMagDescript % y_coord           => YEL
            ENDIF
            ! fort.51 phase
            HAElevStaPhaseDescript % specifier            = NHASE
            HAElevStaPhaseDescript % initial_value        = 0.0
            HAElevStaPhaseDescript % num_items_per_record = MNHARF
            HAElevStaPhaseDescript % num_fd_records       = NSTAE_G
            HAElevStaPhaseDescript % num_records_this     = NSTAE
            HAElevStaPhaseDescript % imap                => IMAP_STAE_LG
            HAElevStaPhaseDescript % array2D              => PHASEDE
            HAElevStaPhaseDescript % array2D_g            => null()
            HAElevStaPhaseDescript % ConsiderWetDry       = .FALSE.
            HAElevStaPhaseDescript % field_name     = 'HAElevStaPhase'
            IF ((MNPROC.gt.1).AND.(WRITE_LOCAL_HARM_FILES.eqv..false.)) THEN
               IF (MyProc.eq.0) THEN
                  ALLOCATE(PHASEDE_g(MNHARF,NSTAE_G))
               ENDIF
               HAElevStaPhaseDescript % array2D_g         => PHASEDE_g
            ENDIF
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
            HAVelStaUMagDescript % specifier            = NHASV
            HAVelStaUMagDescript % initial_value        = 0.0
            HAVelStaUMagDescript % num_items_per_record = MNHARF
            HAVelStaUMagDescript % num_fd_records       = NSTAV_G
            HAVelStaUMagDescript % num_records_this     = NSTAV
            HAVelStaUMagDescript % imap                 => IMAP_STAV_LG
            HAVelStaUMagDescript % array2D              => UMAG
            HAVelStaUMagDescript % array2D_g            => null()
            HAVelStaUMagDescript % ConsiderWetDry       = .FALSE.
            HAVelStaUMagDescript % alternate_value      = 0.0
            HAVelStaUMagDescript % field_name           = 'HAVelStaUMag'
            HAVelStaUMagDescript % file_name            = 'fort.52'
            IF ((MNPROC.gt.1).AND.(WRITE_LOCAL_HARM_FILES.eqv..false.)) THEN
               IF (MyProc.eq.0) THEN
                  ALLOCATE(UMAG_g(MNHARF,NSTAV_G))
               ENDIF
               HAVelStaUMagDescript % array2D_g            => UMAG_g
            ENDIF
            IF (ICS.eq.2) THEN
               HAVelStaUMagDescript % x_coord           => SLEV
               HAVelStaUMagDescript % y_coord           => SFEV
            ELSE
               HAVelStaUMagDescript % x_coord           => XEV
               HAVelStaUMagDescript % y_coord           => YEV
            ENDIF
            ! fort.52 v velocity magnitude
            HAVelStaVMagDescript % specifier            = NHASV
            HAVelStaVMagDescript % initial_value        = 0.0
            HAVelStaVMagDescript % num_items_per_record = MNHARF
            HAVelStaVMagDescript % num_fd_records       = NSTAV_G
            HAVelStaVMagDescript % num_records_this     = NSTAV
            HAVelStaVMagDescript % imap                 => IMAP_STAV_LG
            HAVelStaVMagDescript % array2D              => VMAG
            HAVelStaVMagDescript % array2D_g            => null()
            HAVelStaVMagDescript % ConsiderWetDry       = .FALSE.
            HAVelStaVMagDescript % alternate_value      = 0.0
            HAVelStaVMagDescript % field_name           = 'HAVelStaVMag'
            IF ((MNPROC.gt.1).AND.(WRITE_LOCAL_HARM_FILES.eqv..false.)) THEN
               IF (MyProc.eq.0) THEN
                  ALLOCATE(VMAG_g(MNHARF,NSTAV_G))
               ENDIF
               HAVelStaVMagDescript % array2D_g            => VMAG_g
            ENDIF
            IF (ICS.eq.2) THEN
               HAVelStaVMagDescript % x_coord           => SLEV
               HAVelStaVMagDescript % y_coord           => SFEV
            ELSE
               HAVelStaVMagDescript % x_coord           => XEV
               HAVelStaVMagDescript % y_coord           => YEV
            ENDIF
            ! fort.52 u velocity phase
            HAVelStaUPhaseDescript % specifier            = NHASV
            HAVelStaUPhaseDescript % initial_value        = 0.0
            HAVelStaUPhaseDescript % num_items_per_record = MNHARF
            HAVelStaUPhaseDescript % num_fd_records       = NSTAV_G
            HAVelStaUPhaseDescript % num_records_this     = NSTAV
            HAVelStaUPhaseDescript % imap         => IMAP_STAV_LG
            HAVelStaUPhaseDescript % array2D              => PHASEDU
            HAVelStaUPhaseDescript % array2D_g            => null()
            HAVelStaUPhaseDescript % ConsiderWetDry       = .FALSE.
            HAVelStaUPhaseDescript % alternate_value      = 0.0
            HAVelStaUPhaseDescript % field_name   = 'HAVelStaUPhase'
            IF ((MNPROC.gt.1).AND.(WRITE_LOCAL_HARM_FILES.eqv..false.)) THEN
               IF (MyProc.eq.0) THEN
                  ALLOCATE(PHASEDU_g(MNHARF,NSTAV_G))
               ENDIF
               HAVelStaUPhaseDescript % array2D_g            => PHASEDU_g
            ENDIF
            IF (ICS.eq.2) THEN
               HAVelStaUPhaseDescript % x_coord           => SLEV
               HAVelStaUPhaseDescript % y_coord           => SFEV
            ELSE
               HAVelStaUPhaseDescript % x_coord           => XEV
               HAVelStaUPhaseDescript % y_coord           => YEV
            ENDIF
            ! fort.52 v velocity phase
            HAVelStaVPhaseDescript % specifier            = NHASV
            HAVelStaVPhaseDescript % initial_value        = 0.0
            HAVelStaVPhaseDescript % num_items_per_record = MNHARF
            HAVelStaVPhaseDescript % num_fd_records       = NSTAV_G
            HAVelStaVPhaseDescript % num_records_this     = NSTAV
            HAVelStaVPhaseDescript % imap            => IMAP_STAV_LG
            HAVelStaVPhaseDescript % array2D              => PHASEDV
            HAVelStaVPhaseDescript % array2D_g            => null()
            HAVelStaVPhaseDescript % ConsiderWetDry       = .FALSE.
            HAVelStaVPhaseDescript % alternate_value      = 0.0
            HAVelStaVPhaseDescript % field_name      = 'HAVelStaVPhase'
            IF ((MNPROC.gt.1).AND.(WRITE_LOCAL_HARM_FILES.eqv..false.)) THEN
               IF (MyProc.eq.0) THEN
                  ALLOCATE(PHASEDV_g(MNHARF,NSTAV_G))
               ENDIF
               HAVelStaVPhaseDescript % array2D_g            => PHASEDV_g
            ENDIF
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
            HAElevMagDescript % specifier            = NHAGE
            HAElevMagDescript % initial_value        = 0.0
            HAElevMagDescript % num_items_per_record = MNHARF
            HAElevMagDescript % num_fd_records       = NP_G
            HAElevMagDescript % num_records_this     = NP
            HAElevMagDescript % imap                 => NODES_LG
            HAElevMagDescript % array2D              => EMAGT
            HAElevMagDescript % array2D_g            => null()
            HAElevMagDescript % ConsiderWetDry       = .FALSE.
            HAElevMagDescript % alternate_value      = -99999.0
            HAElevMagDescript % field_name           = 'HAElevMag'
            HAElevMagDescript % file_name            = 'fort.53'
            IF ((MNPROC.gt.1).AND.(WRITE_LOCAL_HARM_FILES.eqv..false.)) THEN
               IF (MyProc.eq.0) THEN
                  ALLOCATE(EMAGT_g(MNHARF,NP_G))
               ENDIF
               HAElevMagDescript % array2D_g            => EMAGT_g
            ENDIF
            ! fort.53 Full Domain Elevation Phase
            HAElevPhaseDescript % specifier            = NHAGE
            HAElevPhaseDescript % initial_value        = 0.0
            HAElevPhaseDescript % num_items_per_record = MNHARF
            HAElevPhaseDescript % num_fd_records       = NP_G
            HAElevPhaseDescript % num_records_this     = NP
            HAElevPhaseDescript % imap                 => NODES_LG
            HAElevPhaseDescript % array2D              => PHASEDEN
            HAElevPhaseDescript % array2D_g            => null()
            HAElevPhaseDescript % ConsiderWetDry       = .FALSE.
            HAElevPhaseDescript % alternate_value      = -99999.0
            HAElevPhaseDescript % field_name           = 'HAElevPhase'
            IF ((MNPROC.gt.1).AND.(WRITE_LOCAL_HARM_FILES.eqv..false.)) THEN
               IF (MyProc.eq.0) THEN
                  ALLOCATE(PHASEDEN_g(MNHARF,NP_G))
               ENDIF
               HAElevPhaseDescript % array2D_g            => PHASEDEN_g
            ENDIF

         ENDIF
C
C
         IF (NHAGV.ne.0) THEN
            ! fort.54 Full Domain U Velocity Magnitude
            HAVelUMagDescript % specifier            = NHAGV
            HAVelUMagDescript % initial_value        = 0.0
            HAVelUMagDescript % num_items_per_record = MNHARF
            HAVelUMagDescript % num_fd_records       = NP_G
            HAVelUMagDescript % num_records_this     = NP
            HAVelUMagDescript % imap                 => NODES_LG
            HAVelUMagDescript % array2D              => UMAGT
            HAVelUMagDescript % array2D_g            => null()
            HAVelUMagDescript % ConsiderWetDry       = .FALSE.
            HAVelUMagDescript % alternate_value      = -99999.0
            HAVelUMagDescript % field_name           = 'HAVelUMag'
            HAVelUMagDescript % file_name            = 'fort.54'
            IF ((MNPROC.gt.1).AND.(WRITE_LOCAL_HARM_FILES.eqv..false.)) THEN
               IF (MyProc.eq.0) THEN
                  ALLOCATE(UMAGT_g(MNHARF,NP_G))
               ENDIF
               HAVelUMagDescript % array2D_g            => UMAGT_g
            ENDIF
            ! fort.54 Full Domain V Velocity Magnitude
            HAVelVMagDescript % specifier            = NHAGV
            HAVelVMagDescript % initial_value        = 0.0
            HAVelVMagDescript % num_items_per_record = MNHARF
            HAVelVMagDescript % num_fd_records       = NP_G
            HAVelVMagDescript % num_records_this     = NP
            HAVelVMagDescript % imap                 => NODES_LG
            HAVelVMagDescript % array2D              => VMAGT
            HAVelVMagDescript % array2D_g            => null()
            HAVelVMagDescript % ConsiderWetDry       = .FALSE.
            HAVelVMagDescript % alternate_value      = -99999.0
            HAVelVMagDescript % field_name           = 'HAVelVMag'
            IF ((MNPROC.gt.1).AND.(WRITE_LOCAL_HARM_FILES.eqv..false.)) THEN
               IF (MyProc.eq.0) THEN
                  ALLOCATE(VMAGT_g(MNHARF,NP_G))
               ENDIF
               HAVelVMagDescript % array2D_g            => VMAGT_g
            ENDIF
            ! fort.54 Full Domain U Velocity Phase
            HAVelUPhaseDescript % specifier            = NHAGV
            HAVelUPhaseDescript % initial_value        = 0.0
            HAVelUPhaseDescript % num_items_per_record = MNHARF
            HAVelUPhaseDescript % num_fd_records       = NP_G
            HAVelUPhaseDescript % num_records_this     = NP
            HAVelUPhaseDescript % imap                 => NODES_LG
            HAVelUPhaseDescript % array2D              => PHASEDUT
            HAVelUPhaseDescript % array2D_g            => null()
            HAVelUPhaseDescript % ConsiderWetDry       = .FALSE.
            HAVelUPhaseDescript % alternate_value      = -99999.0
            HAVelUPhaseDescript % field_name           = 'HAVelUPhase'
            IF ((MNPROC.gt.1).AND.(WRITE_LOCAL_HARM_FILES.eqv..false.)) THEN
               IF (MyProc.eq.0) THEN
                  ALLOCATE(PHASEDUT_g(MNHARF,NP_G))
               ENDIF
               HAVelUPhaseDescript % array2D_g            => PHASEDUT_g
            ENDIF
            ! fort.54 Full Domain V Velocity Phase
            HAVelVPhaseDescript % specifier            = NHAGV
            HAVelVPhaseDescript % initial_value        = 0.0
            HAVelVPhaseDescript % num_items_per_record = MNHARF
            HAVelVPhaseDescript % num_fd_records       = NP_G
            HAVelVPhaseDescript % num_records_this     = NP
            HAVelVPhaseDescript % imap                 => NODES_LG
            HAVelVPhaseDescript % array2D              => PHASEDVT
            HAVelVPhaseDescript % array2D_g            => null()
            HAVelVPhaseDescript % ConsiderWetDry       = .FALSE.
            HAVelVPhaseDescript % alternate_value      = -99999.0
            HAVelVPhaseDescript % field_name           = 'HAVelVPhase'
            IF ((MNPROC.gt.1).AND.(WRITE_LOCAL_HARM_FILES.eqv..false.)) THEN
               IF (MyProc.eq.0) THEN
                  ALLOCATE(PHASEDVT_g(MNHARF,NP_G))
               ENDIF
               HAVelVPhaseDescript % array2D_g            => PHASEDVT_g
            ENDIF
         ENDIF
C
C        ! fort.55 means and variance calculations

C        ! fort.55 will use the same output format (ascii, netcdf, etc)
C        ! as the full domain elevation harmonic analysis file (fort.53)
         IF (CHARMV.eqv..true.) THEN

            ! MVDescript is used over and over for different arrays
            MVDescript % specifier            = NHAGE
            MVDescript % initial_value        = 0.0
            MVDescript % num_items_per_record = 1
            MVDescript % num_fd_records       = NP_G
            MVDescript % num_records_this     = NP
            MVDescript % imap                 => NODES_LG
            MVDescript % array                => ELAV
            MVDescript % array_g              => null()
            MVDescript % ConsiderWetDry       = .false.
            MVDescript % alternate_value      = -99999.0
            MVDescript % field_name           = 'HAMeansVariances'
            IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) .and.
     &         WRITE_LOCAL_HARM_FILES.eqv..false. ) THEN
               ALLOCATE(eav_g(NP_g),eavdif_g(NP_G),
     &                  esq_g(NP_G),evadif_g(NP_G))
               ALLOCATE(uav_g(NP_g),uavdif_g(NP_G),
     &                  usq_g(NP_G),uvadif_g(NP_G))
               ALLOCATE(vav_g(NP_g),vavdif_g(NP_G),
     &                  vsq_g(NP_G),vvadif_g(NP_G))
               MVDescript % array_g              => ELAV_g
            ENDIF
         ENDIF
      ENDIF  ! firstcall
C
C     Collect fulldomain data if running in parallel, and we are not
C     writing subdomain harmonic analysis output files.
#ifdef CMPI
      if (write_local_harm_files.eqv..false.) then
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
      endif
#endif
c
c     Now write data to file; format according to the file type specifier.
C
c     Harmonic analysis at elevation stations (fort.51).
      call writeOutHarmonicArrays(51, 1, HAElevStaMagDescript,
     &        HAElevStaPhaseDescript)
C     Harmonic analysis at velocity stations (fort.52).
      call writeOutHarmonicArrays(52, 2, HAVelStaUMagDescript,
     &        HAVelStaUPhaseDescript, HAVelStaVMagDescript,
     &        HAVelStaVPhaseDescript)
C     Harmonic analysis at all nodes for elevation (fort.53).
      call writeOutHarmonicArrays(53, 1, HAElevMagDescript,
     &        HAElevPhaseDescript)
C     Harmonic analysis at all nodes for velocity (fort.54).
      call writeOutHarmonicArrays(54, 2,HAVelUMagDescript,
     &        HAVelUPhaseDescript, HAVelVMagDescript,
     &        HAVelVPhaseDescript)
C
C     Means and variance calculations, if specified (fort.55).
      if ((CHARMV.eqv..true.).and.(FMV.gt.1.0d-3)) then
         IF ((MNPROC.gt.1).and.(MyProc.eq.0)
     &        .and.(WRITE_LOCAL_HARM_FILES.eqv..false.)) THEN
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
         IF ((MNPROC.eq.1).or.(WRITE_LOCAL_HARM_FILES.eqv..true.)) THEN
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
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
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
     & MYPROC, HARMDIR, WRITE_LOCAL_HARM_FILES
      USE HARM, ONLY : nf, nfreq, hafreq, haff, haface, namefr
      USE GLOBAL, ONLY : OutputDataDescript_t, screenMessage, ERROR,
     &                   DEBUG, setMessageSource, unsetMessageSource,
     &                   allMessage, scratchMessage
      IMPLICIT NONE
      INTEGER, intent(in) :: lun  ! adcirc logical unit number for i/o
      INTEGER, intent(in) :: sets ! how many sets of mags/phases to write
      type(OutputDataDescript_t), intent(in) :: mag1
      type(OutputDataDescript_t), intent(in) :: ph1
      ! jgf51.44: The following two arguments will be present if nsets=2
      type(OutputDataDescript_t), optional, intent(in) :: mag2
      type(OutputDataDescript_t), optional, intent(in) :: ph2
C
      INTEGER :: i,j,n           ! loop counters
C
      call setMessageSource("writeOutHarmonicArrays")

#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
C     write data according to format specifier from fort.15 (e.g., NHASE)
      SELECT CASE (mag1%specifier)
      CASE(0) ! no output; do nothing
      CASE(1) ! ascii text
         IF ( (MNPROC.gt.1).and.(MyProc.eq.0)
     &        .and.(WRITE_LOCAL_HARM_FILES.eqv..false.)) THEN
            OPEN(lun,FILE=TRIM(GLOBALDIR)//'/'//trim(mag1%file_name),
     &         STATUS='REPLACE',ACCESS='SEQUENTIAL',ACTION='WRITE')
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
         IF ((MNPROC.eq.1).or.(WRITE_LOCAL_HARM_FILES.eqv..true.)) THEN
            OPEN(lun,FILE=TRIM(LOCALDIR)//'/'//trim(mag1%file_name),
     &         STATUS='REPLACE',ACCESS='SEQUENTIAL',ACTION='WRITE')
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
         call allMessage(ERROR,
     & 'ERROR: Harmonic analysis output is not yet available in NetCDF.')
      CASE DEFAULT
         write(scratchMessage,'(a,i0,a)')
     &  'Invalid harmonic analysis format specifier: ',mag1%specifier,'.'
         call allMessage(ERROR,scratchMessage)
      END SELECT

 2    FORMAT(I2)
 6635 format(2x,1pE16.8E3,1x,0pf11.4)
 6636 format(2x,1pE16.8E3,1x,0pf11.4,2x,1pE16.8E3,1x,0pf11.4)
 3679 format(1x,e20.10,1x,f10.7,1x,f12.8,1x,a10)

#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
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
      USE GLOBAL, ONLY : RUNDES, RUNID, RDES4, RID4,
     &   RDES8, RID8, DTDP, IHOT, DEBUG, scratchMessage,
     &   logMessage, ERROR, allMessage
      USE GLOBAL_3DVS, ONLY : NFEN
      USE MESH, ONLY : AGRID, AID4, AID8

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
      call setMessageSource("initOutput3D")
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
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

#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
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
     &    unsetMessageSource, allMessage, OutputDataDescript_t, NP_G,
     &    NODES_LG
      USE MESH, ONLY : NP, ICS, NM
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
C
      call setMessageSource("writeOutput3D")
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
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
         SigTStaDescript % lun                  =  41
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
         SalStaDescript % lun                  =  41
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
         TempStaDescript % lun                  =  41
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
         RealQStaDescript % lun                  =  42
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
         ImaginaryQStaDescript % lun                  =  42
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
         WZStaDescript % lun                  =  42
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
         Q20StaDescript % lun                  =  43
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

         LStaDescript % specifier            =  I3DST
         LStaDescript % lun                  =  43
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
         EVStaDescript % lun                  =  43
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
         QSurfKp1Descript % lun                  =  47
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
         SigTDescript % lun                  =  44
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
         SalDescript % lun                  =  44
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
         TempDescript % lun                  =  44
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
         RealQDescript % lun                  =  45
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
         ImaginaryQDescript % lun                  =  45
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
         WZDescript % lun                  =  45
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
         Q20Descript % lun                  =  46
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
         LDescript % lun                  =  46
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
         EVDescript % lun                  =  46
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
            CALL initNetCDFOutputFile(SigTStaDescript,
     &           netcdf_error(1), SalStaDescript, TempStaDescript)
            CALL initNetCDFOutputFile(RealQStaDescript,
     &           netcdf_error(2), ImaginaryQStaDescript, WZStaDescript)
            CALL initNetCDFOutputFile(Q20StaDescript,
     &           netcdf_error(3), LStaDescript, EVStaDescript)
            CALL initNetCDFOutputFile(SigTDescript,
     &           netcdf_error(4), SalDescript, TempDescript)
            CALL initNetCDFOutputFile(RealQDescript,
     &           netcdf_error(5), ImaginaryQDescript, WZDescript)
            CALL initNetCDFOutputFile(Q20Descript,
     &           netcdf_error(6), LDescript, EVDescript)
            CALL initNetCDFOutputFile(QSurfKp1Descript,
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

#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C-----------------------------------------------------------------------
      end subroutine writeOutput3D
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   W R I T E  O U T  A R R A Y
C-----------------------------------------------------------------------
C     jgf48.03 This subroutine was created to write out a column
C     vector (i.e., nodal data such as water surface elevation or
C     pressure) to a file.
C
C     jgf51.21.24: Updated arguments to take into account the more
C     comprehensive information found in the OutputDataDescript_t.
C-----------------------------------------------------------------------
      SUBROUTINE writeOutArray(TimeLoc, it, descript, pack_cmd,
     &                         unpack_cmd)
      USE SIZES
      USE GLOBAL
      USE GLOBAL_IO, ONLY : collectFullDomainArray, storeOne, storeTwo,
     &                      writeSparse
      USE MESH, ONLY: labels
#ifdef ADCNETCDF
      USE NETCDFIO, ONLY : writeOutArrayNetCDF
#endif
#ifdef ADCXDMF
      USE XDMFIO, ONLY : writeOutArrayXDMF
#endif
      IMPLICIT NONE
C     args
      REAL(8), intent(in) :: TimeLoc ! seconds since cold start
      INTEGER, intent(in) :: it   ! number of time steps since cold start
      type(OutputDataDescript_t), intent(inout) :: descript !describes output data
      EXTERNAL :: pack_cmd   ! subroutine used to pack data on subdomain
      EXTERNAL :: unpack_cmd ! subroutine used to unpack data on proc 0
      INTEGER :: I           ! loop counter

      call setMessageSource("writeOutArray")
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
      call allMessage(DEBUG,TRIM(descript%file_name))
#endif

C     collect up the data from subdomains if running in parallel
      IF ((MNPROC.gt.1).and.(WRITE_LOCAL_FILES.eqv..false.)) THEN
         CALL collectFullDomainArray(descript, pack_cmd, unpack_cmd)
      ENDIF

C     write data according to format specifier from fort.15 (e.g., NOUTE)
      SELECT CASE (ABS(descript % specifier))

      CASE(ASCII)

         IF ( (MNPROC.gt.1).and.(MyProc.eq.0)
     &        .and.(.not.WRITE_LOCAL_FILES)) THEN

            OPEN(descript % lun, file=trim(descript % file_name),
     &           ACCESS='SEQUENTIAL',POSITION='APPEND')
            WRITE(descript % lun,2120) TimeLoc,IT
            IF (descript % num_items_per_record .eq. 1) THEN
               DO I=1, descript % num_fd_records
                  if ( descript % isInteger .eqv. .true. ) then
                     write(descript % lun,2452) labels(I), descript % iarray_g(I)
                  else
                     write(descript % lun,2453) I, descript % array_g(I)
                  endif
               ENDDO
            ENDIF
            IF (descript % num_items_per_record .eq. 2) THEN
               DO I=1, descript % num_fd_records
                  WRITE(descript % lun,2454) labels(I), descript % array_g(I),
     &                               descript % array2_g(I)
               ENDDO
            ENDIF
            CLOSE(descript % lun)
         ENDIF

         IF ((MNPROC.eq.1).or.(WRITE_LOCAL_FILES)) THEN
            OPEN(descript % lun, file=trim(descript % file_name),
     &           ACCESS='SEQUENTIAL',POSITION='APPEND')
            WRITE(descript % lun,2120) TimeLoc,IT
            IF (descript % num_items_per_record .eq. 1) THEN
!......TCM - v48.4618 -- Fixed dry node output for serial run case
               IF (descript % ConsiderWetDry .EQV. .TRUE.) THEN
                  DO I=1, descript % num_records_this
                     ! jgf52.25: Station data that should take wet dry state
                     ! into account have already done so in subroutine stationArrayInterp()
                     ! above. For stations, it doesn't make sense to check the nodecode.
                     if ( (descript%isStation.eqv..false.).and.(NODECODE(I).EQ.0)) then
                        if ( descript % isInteger .eqv. .true. ) then
                           write(descript % lun,2452) labels(I), int(descript % alternate_value)
                        else
                           write(descript % lun,2453) labels(I), descript % alternate_value
                        endif
                     else
                        if ( descript % isInteger .eqv. .true. ) then
                           write(descript % lun,2452) I, descript % iarray(i)
                        else
                           write(descript % lun,2453) I, descript % array(i)
                        endif
                     endif
                  end do
               else
                  DO I=1, descript % num_records_this
                     if ( descript % isInteger .eqv. .true. ) then
                        write(descript % lun,2452) labels(I), descript % iarray(I)
                     else
                        write(descript % lun,2453) labels(I), descript % array(I)
                     endif
                  END DO
               ENDIF
            ENDIF
            IF (descript % num_items_per_record .eq. 2) THEN
               IF ( descript % isStation .eqv. .false. ) THEN
                  DO I=1, descript % num_records_this
                     WRITE(descript % lun,2454) labels(I), descript % array(I),
     &                               descript % array2(I)
                  END DO
               ELSE
                  DO I=1, descript % num_records_this
                     WRITE(descript % lun,2454) I, descript % array(I),
     &                               descript % array2(I)
                  END DO
               ENDIF
            ENDIF
            CLOSE(descript % lun)
         ENDIF
         descript % filepos = descript % filepos + 1
     &                        + descript % num_records_this

      CASE(BINARY) ! nonportable

         IF ( (MNPROC.gt.1).and.(MyProc.eq.0)
     &         .and.(.not.WRITE_LOCAL_FILES)) THEN
            OPEN(descript%lun,FILE=trim(descript%file_name),
     &           ACCESS='DIRECT',RECL=NBYTE)
            WRITE(descript % lun,REC=descript % filepos+1) TimeLoc
            WRITE(descript % lun,REC=descript % filepos+2) IT
            descript%filepos = descript%filepos + 2
            IF ( descript % num_items_per_record .eq. 1 ) THEN
               DO I=1, descript % num_fd_records
                  WRITE(descript % lun,REC=descript % filepos+I)
     &               descript % array_g(I)
               END DO
            ENDIF
            IF ( descript % num_items_per_record .eq. 2 ) THEN
               DO I=1, descript % num_fd_records
                  WRITE(descript % lun,REC=descript % filepos+2*I-1)
     &                descript % array_g(I)
                  WRITE(descript % lun,REC=descript % filepos+2*I)
     &                descript % array2_g(I)
               END DO
            ENDIF
            CLOSE(descript%lun)
         ENDIF

         IF ((MNPROC.eq.1).or.(WRITE_LOCAL_FILES)) THEN
            OPEN(descript % lun,FILE=TRIM(descript % file_name),
     &           ACCESS='DIRECT',RECL=NBYTE)
            WRITE(descript % lun,REC=descript % filepos+1) TimeLoc
            WRITE(descript % lun,REC=descript % filepos+2) IT
            descript % filepos = descript % filepos + 2
            IF ( descript % num_items_per_record .eq. 1 ) THEN
               IF ((trim(descript % field_name) .eq. 'Elev').and.
     &             (descript % ConsiderWetDry .EQV. .TRUE.)) THEN
                  DO I=1, descript % num_records_this
                     if(NODECODE(I).EQ.1) THEN
                        WRITE(descript % lun,REC=descript % filepos+I) descript % array(I)
                     ELSE
                        WRITE(descript % lun,REC=descript % filepos+I)
     &                    descript % alternate_value !-99999.0 for dry nodes
                     ENDIF
                  END DO
               ELSE
                  DO I=1, descript % num_records_this
                     WRITE(descript % lun,REC=descript % filepos+I) descript % array(I)
                  END DO
               ENDIF
            ENDIF
            IF ( descript % num_items_per_record .eq. 2 ) THEN
               DO I=1, descript % num_records_this
                  !tcmv48.4618 -- changed from array_g to array
                  WRITE(descript % lun,REC=descript % filepos+2*I-1) descript % array(I)
                  !tcmv48.4618 -- changed from array2_g to array2
                  WRITE(descript % lun,REC=descript % filepos+2*I)   descript % array2(I)
               END DO
            ENDIF
            CLOSE(descript % lun)
         ENDIF
         descript % filepos = descript % filepos + descript % num_records_this

      CASE(SPARSE_ASCII)
         if(descript % num_items_per_record.eq.1)then
             CALL writeSparse(descript, timeLoc, it, storeOne)
         elseif(descript % num_items_per_record.eq.2)then
             CALL writeSparse(descript, timeLoc, it, storeTwo)
         endif

      CASE(NETCDF3, NETCDF4) ! (portable)
#ifdef ADCNETCDF
         IF (MYPROC.EQ.0) THEN
            CALL writeOutArrayNetCDF(descript % lun, TimeLoc, it, descript)
         ENDIF
#else
         call allMessage(ERROR, 'NetCDF was specified.')
         call allMessage(ERROR,
     &   'This executable was not compiled with NetCDF support.')
         call terminate()
#endif

      CASE(XDMF)
#ifdef ADCXDMF
         IF (MYPROC.EQ.0) THEN
            CALL writeOutArrayXDMF(TimeLoc, it, descript)
         ENDIF
#endif
      CASE DEFAULT
         write(scratchMessage,'(a,i0)') 'Invalid output specifier: ',
     &       abs(descript % specifier)
         call allMessage(ERROR, scratchMessage)
      END SELECT

#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

 2    FORMAT(I2)
 2120 FORMAT(2X,1pE20.10E3,5X,I10)
 2452 FORMAT(2x, i8, 2x, i8, 2x, i8, 2x, i8)
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
      SUBROUTINE writeOutArrayMinMax(lun, TimeLoc, it, nrecs_this,
     &                                descript, pack_cmd, unpack_cmd)

      USE SIZES
      USE GLOBAL
      USE MESH, ONLY : NP, labels
      USE GLOBAL_IO, ONLY : collectFullDomainArray, HEADER_MAX,
     &                      open_minmax_file
#ifdef ADCNETCDF
      USE NetCDFIO, ONLY :  initNetCDFOutputFile, writeOutArrayNetCDF
#endif
      IMPLICIT NONE

C     args
      INTEGER, intent(in) :: lun ! logical unit number of file to write to
      REAL(8), intent(in) :: TimeLoc ! seconds since cold start
      INTEGER, intent(in) :: it   ! number of time steps since cold start
      INTEGER, intent(in) :: nrecs_this  !number of records to write tcm v51.20.06
      type(OutputDataDescript_t), intent(inout) :: descript !describes output data
      EXTERNAL :: pack_cmd   ! subroutine used to pack data on subdomain
      EXTERNAL :: unpack_cmd ! subroutine used to unpack data on proc 0
      CHARACTER(len=3) :: cmp ! comparison to make vs hotstart val ('min' or 'max')
      LOGICAL :: nerr       ! error code from netcdf initialization
      INTEGER :: I          ! loop counter
C
      call setMessageSource("writeOutArrayMinMax")
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C     collect up the data from subdomains if running in parallel
      IF ((MNPROC.gt.1).and.(.not.WRITE_LOCAL_FILES)) THEN
         CALL collectFullDomainArray(descript, pack_cmd, unpack_cmd)
      ENDIF

C     write data in ascii text only (other output formats are not supported)
      SELECT CASE (ABS(descript % specifier))

      CASE(OFF,ASCII,BINARY,SPARSE_ASCII) ! just write ascii text

         call logMessage(INFO,'Opening min/max file "'//
     &          trim(descript%file_name)//'" for writing.')

         IF ((MNPROC.gt.1).and.(MyProc.eq.0).and.
     &        (.not.WRITE_LOCAL_FILES)) THEN
            CALL OPEN_MINMAX_FILE(lun, descript%file_name,
     &                             NP_G, NP, nrecs_this,HEADER_MAX)
            OPEN(lun,FILE=TRIM(descript%file_name),
     &           ACCESS='SEQUENTIAL',POSITION='APPEND')
            WRITE(lun,2120) TimeLoc,IT
            DO I=1, descript % num_fd_records
               if ( descript % isInteger .eqv. .true. ) then
                  write(descript % lun,2452) labels(I), descript % iarray_g(I)
               else
                  write(lun,2453) labels(I), descript % array_g(I)
               endif
            ENDDO

            !tcm v51.20.01 additions for time stamp output
            if (descript%minmax_timestamp) then
               WRITE(lun,2120) TimeLoc,IT
               DO I=1, descript % num_fd_records
                  WRITE(lun,2453) labels(I), descript % array2_g(I)
               ENDDO
            endif

            CLOSE(lun)
         ENDIF

         IF ((MNPROC.eq.1).or.(WRITE_LOCAL_FILES)) THEN
            CALL OPEN_MINMAX_FILE(lun, descript%file_name, NP_G, NP,
     &           nrecs_this, HEADER_MAX)
            OPEN(lun,FILE=TRIM(descript%file_name),
     &           ACCESS='SEQUENTIAL',POSITION='APPEND')
            WRITE(lun,2120) TimeLoc,IT
            DO I=1, descript % num_records_this
               if ( descript % isInteger .eqv. .true. ) then
                  write(descript % lun,2452) labels(I), descript % iarray(I)
               else
                  write(lun,2453) labels(I), descript % array(I)
               endif
            ENDDO
            !tcm v51.20.01 additions for time stamp output
            if (descript%minmax_timestamp) then
               WRITE(lun,2120) TimeLoc,IT
               DO I=1, descript % num_records_this
                  WRITE(lun,2453) labels(I), descript % array2(I)
               ENDDO
            endif

            close(lun)  !tcm v51.20.01 added close
         ENDIF
Casey 120830: Reversed the next two lines.
      CASE(NETCDF3,NETCDF4) !netcdf (portable)
#ifdef ADCNETCDF
         IF (MYPROC.EQ.0) THEN
         call logMessage(INFO,'Opening min/max file "'//
     &          trim(descript%file_name)//'" for writing.')
            descript % writeFlag = .true.
            ! jgf52.08.03: TODO: Is there ever a case where
            ! we've gotten to this point and any of the min/max files
            ! have not been initialized?
            if ( descript % initialized.eqv..false.) then
               call initNetCDFOutputFile(descript, nerr)
               descript % initialized = .true.
            endif
            CALL writeOutArrayNetCDF(lun, TimeLoc, it, descript)
         ENDIF
#else
         call allMessage(ERROR,'NetCDF is not available.')
#endif
      CASE DEFAULT
         write(scratchMessage,'(a,i0)') 'Invalid output specifier: ',
     &       abs(descript % specifier)
         call allMessage(ERROR, scratchMessage)
      END SELECT

#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
 2    FORMAT(I2)
 2120 FORMAT(2X,1pE20.10E3,5X,I10)
 2452 FORMAT(2x, i8, 2x, i8, 2x, i8, 2x, i8)
 2453 FORMAT(2x, i8, 2x, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3)
 2454 FORMAT(2X,I8,2(2X,1pE20.10E3))
C-----------------------------------------------------------------------
       END SUBROUTINE writeOutArrayMinMax
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C        S U B R O U T I N E    W R I T E   H O T S T A R T
C-----------------------------------------------------------------------
C     jgf48.03 This subroutine was created from code in timestep.F to
C     write hotstart files.
C-----------------------------------------------------------------------
      SUBROUTINE writeHotstart(TimeLoc, IT)

      USE SIZES
      USE GLOBAL
      USE HARM
      USE GLOBAL_IO, ONLY: packOne, unpackOne,
     &    packTwo, unpackTwo, packMbyNP, unpackMbyNP, HEADER_MAX,
     &    collectFullDomainArray
      USE MESH, ONLY : NP, NE
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
      CHARACTER*9 :: itstr  !tcm v51.26 added for time-stamped hot start file name adjustments
C
      LOGICAL ncerror
      INTEGER numHotstartWrites ! number writes to hot start files
      INTEGER nextLun           ! next LUN to write to, after initial write
C
C     !jgf49.44: harmonic analysis vars
      INTEGER N ! number of stations
      INTEGER NSTAEX  ! num of elev stations in either fulldomain or subdomain
      INTEGER NSTAVX  ! num of vel stations in either fulldomain or subdomain
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

         HotstartVelDescript % specifier            = NHSTAR
         HotstartVelDescript % initial_value        = 0.0
         HotstartVelDescript % num_items_per_record = 2
         HotstartVelDescript % num_fd_records       = NP_G
         HotstartVelDescript % num_records_this     = NP
         HotstartVelDescript % imap                 => NODES_LG
         HotstartVelDescript % array                => UU2
         HotstartVelDescript % array2               => VV2
         HotstartVelDescript % array_g              => UU2_g
         HotstartVelDescript % array2_g             => VV2_g

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
         NodeCodeDescript % initial_value        = 0.d0
         NodeCodeDescript % num_items_per_record = 1
         NodeCodeDescript % num_fd_records       = NP_G
         NodeCodeDescript % num_records_this     = NP
         NodeCodeDescript % imap                 => NODES_LG
         NodeCodeDescript % iarray               => NODECODE
         NodeCodeDescript % iarray_g             => NODECODE_g
         NodeCodeDescript % isInteger            = .true.

         NOFFDescript % specifier            = NHSTAR
         NOFFDescript % initial_value        = 0.d0
         NOFFDescript % num_items_per_record = 1
         NOFFDescript % num_fd_records       = NE_G
         NOFFDescript % num_records_this     = NE
         NOFFDescript % imap                 => IMAP_EL_LG
         NOFFDescript % iarray               => NOFF
         NOFFDescript % iarray_g             => NOFF_g
         NOFFDescript % isInteger            = .true.
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
               IF ( ncerror.eqv..true. ) CALL terminate()
               IF (IHARIND.eq.1) THEN
                  CALL initNetCDFHotstartHarmonic(hss%lun,
     &            HarmElevFDLVDescript, HarmElevSLVDescript,
     &            HarmUVelFDLVDescript, HarmVVelFDLVDescript,
     &            HarmUVelSLVDescript, HarmVVelSLVDescript, ncerror)
                  IF ( ncerror.eqv..true. ) CALL terminate()
                  IF (CHARMV.eqv..true.) THEN
                     CALL initNetCDFHotstartHarmonicMeansVariances(
     &                  hss%lun, ELAVDescript, ELVADescript,
     &                  XVELAVDescript, YVELAVDescript, XVELVADescript,
     &                  YVELVADescript,ncerror)
                     IF ( ncerror.eqv..true. ) CALL terminate()
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
                  IF ( ncerror.eqv..true. ) CALL terminate()
                  IF (IHARIND.eq.1) THEN
                     CALL initNetCDFHotstartHarmonic(nextLun,
     &               HarmElevFDLVDescript, HarmElevSLVDescript,
     &               HarmUVelFDLVDescript, HarmVVelFDLVDescript,
     &               HarmUVelSLVDescript, HarmVVelSLVDescript, ncerror)
                     IF ( ncerror.eqv..true. ) CALL terminate()
                     IF (CHARMV.eqv..true.) THEN
                        CALL initNetCDFHotstartHarmonicMeansVariances(
     &                     nextLun, ELAVDescript, ELVADescript,
     &                     XVELAVDescript, YVELAVDescript,
     &                     XVELVADescript, YVELVADescript, ncerror)
                        IF ( ncerror.eqv..true. ) CALL terminate()
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
         CALL collectFullDomainArray(NodeCodeDescript,
     &             packOne, unpackOne)
         CALL collectFullDomainArray(NOFFDescript,
     &             packOne, unpackOne)

         IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
            if (any(NOFF_G.gt.1)) then
               call allMessage(ERROR,
     &            'Some values of NOFF are greater than 1.')
               call terminate()
            endif
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
      !tcm v51.26 mod for time-stamped nhstar=-1
      CASE(-1,1,67,68) ! 1=nonportable binary
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
            !tcm v51.26 mod for time-stamped nhstar=-1
            IF (NHSTAR.eq.-1) THEN
               hss % filename(:) = ' '
               hss % filename = 'fort.68_'
               itstr(:) = ' '
               WRITE(itstr,'(I9.9)') IT
               hss % filename(9:17) = itstr
            ENDIF
            OPEN(hss % lun ,FILE=TRIM(LOCALDIR)//'/'//trim(hss % filename),
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
            IF ((MNPROC.eq.1).or.(WRITE_LOCAL_HOT_START_FILES.eqv..true.)) THEN
               NSTAEX = NSTAE
               NSTAVX = NSTAV
            ELSE
               NSTAEX = NSTAE_G
               NSTAVX = NSTAV_G
            ENDIF
            !
            ! In serial (or in parallel with a fulldomain hotstart file)
            ! only processor 0 will write these data to the hotstart file.
            ! In parallel, with subdomain hotstart files, all processors
            ! will write these data to subdomain hotstart files.
            IF ((MyPROC.eq.0).or.(WRITE_LOCAL_HOT_START_FILES.eqv..true.)) THEN
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
            ENDIF
            !
            ! In serial, or in parallel with subdomain hotstart files,
            ! write load vectors to hotstart file(s).
            !
            if ( (mnproc.eq.1).or.
     &            (write_local_hot_start_files.eqv..true.) ) then
               IF(NHASE.ne.0) THEN
                  do n=1,NSTAE
                     do i=1,mm
                        WRITE(hss % lun,REC=IHOTSTP) STAELV(I,N)
                        IHOTSTP=IHOTSTP+1
                     end do
                  end do
               ENDIF
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
               IF(NHAGE.ne.0) THEN
                  do n=1,np
                     do i=1,mm
                        WRITE(hss % lun,REC=IHOTSTP) GLOELV(I,N)
                        IHOTSTP=IHOTSTP+1
                     end do
                  end do
               ENDIF
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
            endif
            !
            ! In parallel, if fulldomain hotstart files are required,
            ! processor 0 will write the full domain load vectors to the
            ! fulldomain hotstart file.
            IF ( (mnproc.gt.1).and.
     &           (write_local_hot_start_files.eqv..false.).and.
     &           (myproc.eq.0)) THEN
               IF(NHASE.ne.0) THEN
                  do n=1,NSTAE_G
                     do i=1,mm
                        WRITE(hss % lun,REC=IHOTSTP)
     &                        HarmElevSLVDescript % array2D_g(i,n)
                        IHOTSTP=IHOTSTP+1
                     end do
                  end do
               ENDIF
               IF(NHASV.ne.0) THEN
                  do N=1,NSTAV_G
                     do i=1,mm
                        WRITE(hss % lun,REC=IHOTSTP)
     &                           HarmUVelSLVDescript % array2D_g(i,n)
                        IHOTSTP=IHOTSTP+1
                        WRITE(hss % lun,REC=IHOTSTP)
     &                           HarmVVelSLVDescript % array2D_g(i,n)
                        IHOTSTP=IHOTSTP+1
                     end do
                  end do
               ENDIF
               IF(NHAGE.ne.0) THEN
                  do n=1,np_g
                     do i=1,mm
                        WRITE(hss % lun,REC=IHOTSTP)
     &                        HarmElevFDLVDescript % array2D_g(i,n)
                        IHOTSTP=IHOTSTP+1
                     end do
                  end do
               ENDIF
               IF(NHAGV.ne.0) THEN
                  do n=1,np_g
                     do i=1,mm
                        WRITE(hss % lun,REC=IHOTSTP)
     &                       HarmUVelFDLVDescript % array2D_g(i,n)
                        IHOTSTP=IHOTSTP+1
                        WRITE(hss % lun,REC=IHOTSTP)
     &                       HarmVVelFDLVDescript % array2D_g(i,n)
                        IHOTSTP=IHOTSTP+1
                     end do
                  end do
               ENDIF
            ENDIF
         ENDIF  ! IHARIND
C
         if ((CHARMV.eqv..true.).and.(IHARIND.eq.1).and.(IT.gt.ITMV)) then
            if ( (MNPROC.eq.1).or.
     &           (WRITE_LOCAL_HOT_START_FILES.eqv..true.)) then
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
C
            IF ((MNPROC.gt.1).and.
     &           (write_local_hot_start_files.eqv..false.).and.
     &           (MYPROC.eq.0)) then
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
         OPEN(hss%lun,FILE=TRIM(INPUTDIR)//'/'//trim(hss % filename),
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
            OPEN(hss%lun,FILE=TRIM(GBLINPUTDIR)//'/'//trim(hss % filename),
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
            hss % filename(:) = ' '
            hss % filename = 'fort.68'
         ENDIF
         !tcm v51.26 mod for time-stamped nhstar=-1  This section is just a failsafe
         IF (NHSTAR.eq.-1) THEN
            hss % filename(:) = ' '
            hss % filename = 'fort.68_'
            itstr(:) = ' '
            WRITE(itstr,'(I9.9)') IT+NHSINC  !setting the name for the next hot start file
            hss % filename(9:17) = itstr
         ENDIF
      ELSE
         IF ((NHSTAR.EQ.68).OR.(NHSTAR.EQ.268).OR.(NHSTAR.EQ.368).OR.
     &       (NHSTAR.EQ.568)) THEN
            call allMessage(INFO,"NHSTAR.EQ.68; ADCIRC stopping.")
            TEST_HOTSTART = .TRUE.
         ENDIF
         hss % lun=67
         IF((NHSTAR.EQ.1).OR.(NHSTAR.EQ.2)) THEN
            hss % filename(:) = ' '
            hss % filename = 'fort.67'
         ENDIF
         !tcm v51.26 mod for time-stamped nhstar=-1
         IF (NHSTAR.eq.-1) THEN
            hss % filename(:) = ' '
            hss % filename = 'fort.68_'
            itstr(:) = ' '
            WRITE(itstr,'(I9.9)') IT+NHSINC  !setting the name for the next hot start file
            hss % filename(9:17) = itstr
         ENDIF
      ENDIF
      IF (TEST_HOTSTART.eqv..true.) THEN
#ifdef CMPI
          CALL MSG_FINI()
#endif
         STOP
      ENDIF

#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
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
      USE GLOBAL_IO, ONLY: collectFullDomainArray, unpackOne
      USE MESH, ONLY : NP
      USE GLOBAL_IO, ONLY: packOne, unpackOne
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

      call setMessageSource("writeWarnElev")
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
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
         DebugElevDescript % lun                  = 69
         DebugElevDescript % specifier            = ASCII ! TODO: make configurable
         DebugElevDescript % initial_value        = 0.0
         DebugElevDescript % num_items_per_record = 1
         DebugElevDescript % num_fd_records       = NP_G
         DebugElevDescript % num_records_this     = NP
         DebugElevDescript % imap                 => NODES_LG
         DebugElevDescript % array                => ETA2
         DebugElevDescript % array_g              => ETA_DEBUG_g
         DebugElevDescript % ConsiderWetDry       = .TRUE.
         DebugElevDescript % filepos              => IGDP
         DebugElevDescript % alternate_value      = -99999.0
         DebugElevDescript % field_name           = 'DebugElev'
         DebugElevDescript % writeFlag            = .true.
         DebugElevDescript % isStation            = .false.
         DebugElevDescript % divideByDepth        = .false.
         allocate(DebugElevDescript % writerFormats(numFormats))
         DebugElevDescript % writerFormats(:)     = -99999
         DebugElevDescript % useWriter            = .false.
         DebugElevDescript % file_extension       = 69
         DebugElevDescript % file_basename        = 'fort'
         DebugElevDescript % file_name            = 'fort.69'

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
      CALL writeOutArray(timeLoc, IT, DebugElevDescript, packOne, unpackOne)
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
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE WriteWarnElev
C-----------------------------------------------------------------------

      SUBROUTINE terminate(NO_MPI_FINALIZE)
#ifdef CMPI
      USE MESSENGER
#endif
      USE GLOBAL, ONLY : setMessageSource, unsetMessageSource,
     &   allMessage, DEBUG, ECHO, INFO, WARNING, ERROR
      IMPLICIT NONE
      LOGICAL, OPTIONAL :: NO_MPI_FINALIZE
C
      call setMessageSource("terminate")
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      call allMessage(INFO,"ADCIRC Terminating.")

#ifdef CMPI
      subdomainFatalError = .true.
      IF (PRESENT(NO_MPI_FINALIZE)) THEN
        CALL MSG_FINI(NO_MPI_FINALIZE)
      ELSE
        CALL MSG_FINI()
      ENDIF
#endif
      STOP
C
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      call unsetMessageSource()
      END SUBROUTINE terminate

C-----------------------------------------------------------------------
C     S U B R O U T I N E  i n i t S w a n O u t p u t
C-----------------------------------------------------------------------
C THIS ROUTINE WILL PERFORM THE ASSIGNMENT AND ALLOCATIONS FOR THE
C SWAN OUTPUT VARIABLES USED IN WRITE_OUTPUT.F
C-----------------------------------------------------------------------
#ifdef CSWAN
      SUBROUTINE initSwanOutput()
          USE SIZES, ONLY : MNPROC, OFF, ASCII, BINARY, NETCDF3,
     &                       NETCDF4, XDMF, MYPROC
          USE GLOBAL
          USE MESH
          IMPLICIT NONE

      call setMessageSource("initSwanOutput")
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

C.....Allocate SWAN specific local/global arrays if that output is selected
      IF(SWAN_OutputHS)THEN
        IF(MYPROC.EQ.0.AND.MNPROC.GT.1)THEN
            ALLOCATE(SWAN_HSOut_g(1:NP_G))
            ALLOCATE(SWAN_HSMaxOut_g(1:NP_G))
        ENDIF
        ALLOCATE(SWAN_HSOut(NP))
        ALLOCATE(SWAN_HSMaxOut(NP))
        SWAN_HSMaxOut(:) = -99999D0
      ENDIF
      IF(SWAN_OutputDIR)THEN
        IF(MYPROC.EQ.0.AND.MNPROC.GT.1)THEN
            ALLOCATE(SWAN_DIROut_g(NP_G))
            ALLOCATE(SWAN_DIRMaxOut_g(NP_G))
        ENDIF
        ALLOCATE(SWAN_DIROut(NP))
        ALLOCATE(SWAN_DIRMaxOut(NP))
        SWAN_DIRMaxOut(:) = -99999D0
      ENDIF
      IF(SWAN_OutputTM01)THEN
        IF(MYPROC.EQ.0.AND.MNPROC.GT.1)THEN
            ALLOCATE(SWAN_TM01Out_g(NP_G))
            ALLOCATE(SWAN_TM01MaxOut_g(NP_G))
        ENDIF
        ALLOCATE(SWAN_TM01Out(NP))
        ALLOCATE(SWAN_TM01MaxOut(NP))
        SWAN_TM01MaxOut(:) = -99999D0
      ENDIF
      IF(SWAN_OutputTPS)THEN
        IF(MYPROC.EQ.0.AND.MNPROC.GT.1)THEN
            ALLOCATE(SWAN_TPSOut_g(NP_G))
            ALLOCATE(SWAN_TPSMaxOut_g(NP_G))
        ENDIF
        ALLOCATE(SWAN_TPSOut(NP))
        ALLOCATE(SWAN_TPSMaxOut(NP))
        SWAN_TPSMaxOut(:) = -99999D0
      ENDIF
      IF(SWAN_OutputTM02)THEN
        IF(MYPROC.EQ.0.AND.MNPROC.GT.1)THEN
            ALLOCATE(SWAN_TM02Out_g(NP_G))
            ALLOCATE(SWAN_TM02MaxOut_g(NP_G))
        ENDIF
        ALLOCATE(SWAN_TM02Out(NP))
        ALLOCATE(SWAN_TM02MaxOut(NP))
        SWAN_TM02MaxOut(:) = -99999D0
      ENDIF
      IF(SWAN_OutputTMM10)THEN
        IF(MYPROC.EQ.0.AND.MNPROC.GT.1)THEN
            ALLOCATE(SWAN_TMM10Out_g(NP_G))
            ALLOCATE(SWAN_TMM10MaxOut_g(NP_G))
        ENDIF
        ALLOCATE(SWAN_TMM10Out(NP))
        ALLOCATE(SWAN_TMM10MaxOut(NP))
        SWAN_TMM10MaxOut(:) = -99999D0
      ENDIF
      IF(SWAN_OutputWind)THEN
        IF(MYPROC.EQ.0.AND.MNPROC.GT.1)THEN
            ALLOCATE(Swan_WindXOut_g(NP_G),Swan_WindYOut_g(NP_G))
            ALLOCATE(Swan_WindMaxOut_g(NP_G))
        ENDIF
        ALLOCATE(SWAN_WindXOut(NP),SWAN_WindYOut(NP))
        ALLOCATE(SWAN_WindMaxOut(NP))
        SWAN_WindMaxOut(:) = -99999D0
      ENDIF
C.....Significant Wave Height
      IF(SWAN_OutputHS)THEN
          SwanHSDescript % specifier        = NOUTGW
      ELSE
          SwanHSDescript % specifier        = OFF
      ENDIF
      SwanHSDescript % lun                  = 301
      SwanHSDescript % filepos              => SWAN_HS_POS
      SwanHSDescript % initial_value        = 0.0
      SwanHSDescript % num_items_per_record = 1
      SwanHSDescript % num_fd_records       = NP_G
      SwanHSDescript % num_records_this     = NP
      SwanHSDescript % imap                 => NODES_LG
      SwanHSDescript % array                => SWAN_HSOut
      SwanHSDescript % array_g              => SWAN_HSOut_g
      SwanHSDescript % ConsiderWetDry       = .FALSE.
      SwanHSDescript % alternate_value      = -99999.0
      SwanHSDescript % field_name           = "swan_HS"
      SwanHSDescript % file_basename        = "swan_HS"
      SwanHSDescript % file_extension       = 63
      SwanHSDescript % startTimeStep        = NTCYSGW
      SwanHSDescript % endTimeStep          = NTCYFGW
      SwanHSDescript % outputTimeStepIncrement = NSPOOLGW
      SwanHSDescript % spoolCounter         => SWAN_HS_SPOOL
      SwanHSDescript % writerFormats(1:4)   = (/ 1, 3, 4, 5 /)
      SwanHSDescript % writeFlag            = .true.

      IF(Swan_OutputHS)THEN
          SwanHSMaxDescript % specifier        = NOUTGW
      ELSE
          SwanHSMaxDescript % specifier        = OFF
      ENDIF
      SwanHSMaxDescript % lun                  = 316
      SwanHSMaxDescript % initial_value        = 0.0
      SwanHSMaxDescript % num_items_per_record = 1
      SwanHSMaxDescript % num_fd_records       = NP_G
      SwanHSMaxDescript % num_records_this     = NP
      SwanHSMaxDescript % imap                 => NODES_LG
      SwanHSMaxDescript % array                => SWAN_HSMaxOut
      SwanHSMaxDescript % array_g              => SWAN_HSMaxOut_g
      SwanHSMaxDescript % ConsiderWetDry       = .FALSE.
      SwanHSMaxDescript % alternate_value      = -99999.0
      SwanHSMaxDescript % field_name           = "swan_HS_max"
      SwanHSMaxDescript % file_basename        = "swan_HS_max"
      SwanHSMaxDescript % file_extension       = 63
      SwanHSMaxDescript % writeFlag            = .false.

C.....Mean Wave Direction
      IF(SWAN_OutputDIR)THEN
          SwanDIRDescript % specifier        = NOUTGW
      ELSE
          SwanDIRDescript % specifier        = OFF
      ENDIF
      SwanDIRDescript % lun                  = 302
      SwanDIRDescript % filepos              => SWAN_DIR_POS
      SwanDIRDescript % initial_value        = 0.0
      SwanDIRDescript % num_items_per_record = 1
      SwanDIRDescript % num_fd_records       = NP_G
      SwanDIRDescript % num_records_this     = NP
      SwanDIRDescript % imap                 => NODES_LG
      SwanDIRDescript % array                => SWAN_DIROut
      SwanDIRDescript % array_g              => SWAN_DIROut_g
      SwanDIRDescript % ConsiderWetDry       = .FALSE.
      SwanDIRDescript % alternate_value      = -99999.0
      SwanDIRDescript % file_basename        = "swan_DIR"
      SwanDIRDescript % file_extension       = 63
      SwanDIRDescript % field_name           = "swan_DIR"
      SwanDIRDescript % startTimeStep        = NTCYSGW
      SwanDIRDescript % endTimeStep          = NTCYFGW
      SwanDIRDescript % outputTimeStepIncrement = NSPOOLGW
      SwanDIRDescript % spoolCounter         => SWAN_DIR_SPOOL
      SwanDIRDescript % writerFormats(1:4)   = (/ 1, 3, 4, 5 /)
      SwanDIRDescript % writeFlag            = .true.

      IF(SWAN_OutputDir)THEN
          SwanDIRMaxDescript % specifier        = NOUTGW
      ELSE
          SwanDIRMaxDescript % specifier        = OFF
      ENDIF
      SwanDIRMaxDescript % lun                  = 317
      SwanDIRMaxDescript % initial_value        = 0.0
      SwanDIRMaxDescript % num_items_per_record = 1
      SwanDIRMaxDescript % num_fd_records       = NP_G
      SwanDIRMaxDescript % num_records_this     = NP
      SwanDIRMaxDescript % imap                 => NODES_LG
      SwanDIRMaxDescript % array                => SWAN_DIRMaxOut
      SwanDIRMaxDescript % array_g              => SWAN_DIRMaxOut_g
      SwanDIRMaxDescript % ConsiderWetDry       = .FALSE.
      SwanDIRMaxDescript % alternate_value      = -99999.0
      SwanDIRMaxDescript % file_basename        = "swan_DIR_max"
      SwanDIRMaxDescript % file_extension       = 63
      SwanDIRMaxDescript % field_name           = "swan_DIR_max"
      SwanDIRMaxDescript % writeFlag            = .false.

C.....Mean Wave Period (TM01)
      IF(SWAN_OutputTM01)THEN
          SwanTM01Descript % specifier        = NOUTGW
      ELSE
          SwanTM01Descript % specifier        = OFF
      ENDIF
      SwanTM01Descript % lun                  = 303
      SwanTM01Descript % filepos              => SWAN_TM01_POS
      SwanTM01Descript % initial_value        = 0.0
      SwanTM01Descript % num_items_per_record = 1
      SwanTM01Descript % num_fd_records       = NP_G
      SwanTM01Descript % num_records_this     = NP
      SwanTM01Descript % imap                 => NODES_LG
      SwanTM01Descript % array                => SWAN_TM01Out
      SwanTM01Descript % array_g              => SWAN_TM01Out_g
      SwanTM01Descript % ConsiderWetDry       = .FALSE.
      SwanTM01Descript % alternate_value      = -99999.0
      SwanTM01Descript % file_basename        = "swan_TM01"
      SwanTM01Descript % file_extension       = 63
      SwanTM01Descript % field_name           = "swan_TM01"
      SwanTM01Descript % startTimeStep        = NTCYSGW
      SwanTM01Descript % endTimeStep          = NTCYFGW
      SwanTM01Descript % outputTimeStepIncrement = NSPOOLGW
      SwanTM01Descript % spoolCounter         => SWAN_TM01_SPOOL
      SwanTM01Descript % writerFormats(1:4)   = (/ 1, 3, 4, 5 /)
      SwanTM01Descript % writeFlag            = .true.

      IF(SWAN_OutputTM01)THEN
          SwanTM01MaxDescript % specifier        = NOUTGW
      ELSE
          SwanTM01MaxDescript % specifier        = OFF
      ENDIF
      SwanTM01MaxDescript % lun                  = 318
      SwanTM01MaxDescript % initial_value        = 0.0
      SwanTM01MaxDescript % num_items_per_record = 1
      SwanTM01MaxDescript % num_fd_records       = NP_G
      SwanTM01MaxDescript % num_records_this     = NP
      SwanTM01MaxDescript % imap                 => NODES_LG
      SwanTM01MaxDescript % array                => SWAN_TM01MaxOut
      SwanTM01MaxDescript % array_g              => SWAN_TM01MaxOut_g
      SwanTM01MaxDescript % ConsiderWetDry       = .FALSE.
      SwanTM01MaxDescript % alternate_value      = -99999.0
      SwanTM01MaxDescript % file_basename        = "swan_TM01_max"
      SwanTM01MaxDescript % file_extension       = 63
      SwanTM01MaxDescript % field_name           = "swan_TM01_max"
      SwanTM01MaxDescript % writeFlag            = .false.

C.....Peak Wave Period (TPS)
      IF(SWAN_OutputTPS)THEN
          SwanTPSDescript % specifier        = NOUTGW
      ELSE
          SwanTPSDescript % specifier        = OFF
      ENDIF
      SwanTPSDescript % lun                  = 304
      SwanTPSDescript % filepos              => SWAN_TPS_POS
      SwanTPSDescript % initial_value        = 0.0
      SwanTPSDescript % num_items_per_record = 1
      SwanTPSDescript % num_fd_records       = NP_G
      SwanTPSDescript % num_records_this     = NP
      SwanTPSDescript % imap                 => NODES_LG
      SwanTPSDescript % array                => SWAN_TPSOut
      SwanTPSDescript % array_g              => SWAN_TPSOut_g
      SwanTPSDescript % ConsiderWetDry       = .FALSE.
      SwanTPSDescript % alternate_value      = -99999.0
      SwanTPSDescript % file_basename        = "swan_TPS"
      SwanTPSDescript % file_extension       = 63
      SwanTPSDescript % field_name           = "swan_TPS"
      SwanTPSDescript % startTimeStep        = NTCYSGW
      SwanTPSDescript % endTimeStep          = NTCYFGW
      SwanTPSDescript % outputTimeStepIncrement = NSPOOLGW
      SwanTPSDescript % spoolCounter         => SWAN_TPS_SPOOL
      SwanTPSDescript % writerFormats(1:4)   = (/ 1, 3, 4, 5 /)
      SwanTPSDescript % writeFlag            = .true.

      IF(SWAN_OutputTPS)THEN
          SwanTPSMaxDescript % specifier        = NOUTGW
      ELSE
          SwanTPSMaxDescript % specifier        = OFF
      ENDIF
      SwanTPSMaxDescript % lun                  = 319
      SwanTPSMaxDescript % initial_value        = 0.0
      SwanTPSMaxDescript % num_items_per_record = 1
      SwanTPSMaxDescript % num_fd_records       = NP_G
      SwanTPSMaxDescript % num_records_this     = NP
      SwanTPSMaxDescript % imap                 => NODES_LG
      SwanTPSMaxDescript % array                => SWAN_TPSMaxOut
      SwanTPSMaxDescript % array_g              => SWAN_TPSMaxOut_g
      SwanTPSMaxDescript % ConsiderWetDry       = .FALSE.
      SwanTPSMaxDescript % alternate_value      = -99999.0
      SwanTPSMaxDescript % file_basename        = "swan_TPS_max"
      SwanTPSMaxDescript % file_extension       = 63
      SwanTPSMaxDescript % field_name           = "swan_TPS_max"
      SwanTPSMaxDescript % writeFlag            = .false.

C.....SWAN Wind Values (WINDX,WINDY)
      IF(SWAN_OutputWind)THEN
          SwanWindDescript % specifier        = NOUTGW
      ELSE
          SwanWindDescript % specifier        = OFF
      ENDIF
      SwanWindDescript % lun                  = 305
      SwanWindDescript % filepos              => SWAN_WIND_POS
      SwanWindDescript % initial_value        = 0.0
      SwanWindDescript % num_items_per_record = 2
      SwanWindDescript % num_fd_records       = NP_G
      SwanWindDescript % num_records_this     = NP
      SwanWindDescript % imap                 => NODES_LG
      SwanWindDescript % array                => SWAN_WindXOut
      SwanWindDescript % array_g              => SWAN_WindXOut_g
      SwanWindDescript % array2               => SWAN_WindYOut
      SwanWindDescript % array2_g             => SWAN_WindYOut_g
      SwanWindDescript % ConsiderWetDry       = .FALSE.
      SwanWindDescript % alternate_value      = -99999.0
      SwanWindDescript % file_basename        = "swan_WIND"
      SwanWindDescript % file_extension       = 64
      SwanWindDescript % field_name           = "swan_WIND"
      SwanWindDescript % startTimeStep        = NTCYSGW
      SwanWindDescript % endTimeStep          = NTCYFGW
      SwanWindDescript % outputTimeStepIncrement = NSPOOLGW
      SwanWindDescript % spoolCounter         => SWAN_WIND_SPOOL
      SwanWindDescript % writerFormats(1:4)   = (/ 1, 3, 4, 5 /)
      SwanWindDescript % writeFlag            = .true.

      IF(SWAN_OutputWind)THEN
          SwanWindMaxDescript % specifier        = NOUTGW
      ELSE
          SwanWindMaxDescript % specifier        = OFF
      ENDIF
      SwanWindMaxDescript % lun                  = 320
      SwanWindMaxDescript % initial_value        = 0.0
      SwanWindMaxDescript % num_items_per_record = 1
      SwanWindMaxDescript % num_fd_records       = NP_G
      SwanWindMaxDescript % num_records_this     = NP
      SwanWindMaxDescript % imap                 => NODES_LG
      SwanWindMaxDescript % array                => SWAN_WindMaxOut
      SwanWindMaxDescript % array_g             => SWAN_WindMaxOut_g
      SwanWindMaxDescript % ConsiderWetDry       = .FALSE.
      SwanWindMaxDescript % alternate_value      = -99999.0
      SwanWindMaxDescript % file_basename        = "swan_WIND_max"
      SwanWindMaxDescript % file_extension       = 63
      SwanWindMaxDescript % field_name           = "swan_WIND_max"
      SwanWindMaxDescript % writeFlag            = .false.

C.....Mean Wave Period (TM02)
      IF(SWAN_OutputTM02)THEN
          SwanTM02Descript % specifier        = NOUTGW
      ELSE
          SwanTM02Descript % specifier        = OFF
      ENDIF
      SwanTM02Descript % lun                  = 306
      SwanTM02Descript % filepos              => SWAN_TM02_POS
      SwanTM02Descript % initial_value        = 0.0
      SwanTM02Descript % num_items_per_record = 1
      SwanTM02Descript % num_fd_records       = NP_G
      SwanTM02Descript % num_records_this     = NP
      SwanTM02Descript % imap                 => NODES_LG
      SwanTM02Descript % array                => SWAN_TM02Out
      SwanTM02Descript % array_g              => SWAN_TM02Out_g
      SwanTM02Descript % ConsiderWetDry       = .FALSE.
      SwanTM02Descript % alternate_value      = -99999.0
      SwanTM02Descript % file_basename        = "swan_TM02"
      SwanTM02Descript % file_extension       = 63
      SwanTM02Descript % field_name           = "swan_TM02"
      SwanTM02Descript % startTimeStep        = NTCYSGW
      SwanTM02Descript % endTimeStep          = NTCYFGW
      SwanTM02Descript % outputTimeStepIncrement = NSPOOLGW
      SwanTM02Descript % spoolCounter         => SWAN_TM02_SPOOL
      SwanTM02Descript % writerFormats(1:4)   = (/ 1, 3, 4, 5 /)
      SwanTM02Descript % writeFlag            = .true.

      IF(SWAN_OutputTM02)THEN
          SwanTM02MaxDescript % specifier        = NOUTGW
      ELSE
          SwanTM02MaxDescript % specifier        = OFF
      ENDIF
      SwanTM02MaxDescript % lun                  = 321
      SwanTM02MaxDescript % initial_value        = 0.0
      SwanTM02MaxDescript % num_items_per_record = 1
      SwanTM02MaxDescript % num_fd_records       = NP_G
      SwanTM02MaxDescript % num_records_this     = NP
      SwanTM02MaxDescript % imap                 => NODES_LG
      SwanTM02MaxDescript % array                => SWAN_TM02MaxOut
      SwanTM02MaxDescript % array_g             => SWAN_TM02MaxOut_g
      SwanTM02MaxDescript % ConsiderWetDry       = .FALSE.
      SwanTM02MaxDescript % alternate_value      = -99999.0
      SwanTM02MaxDescript % file_basename        = "swan_TM02_max"
      SwanTM02MaxDescript % file_extension       = 63
      SwanTM02MaxDescript % field_name           = "swan_TM02_max"
      SwanTM02MaxDescript % writeFlag            = .false.

C.....Mean Wave Period (TMM10)
      IF(SWAN_OutputTMM10)THEN
          SwanTMM10Descript % specifier        = NOUTGW
      ELSE
          SwanTMM10Descript % specifier        = OFF
      ENDIF
      SwanTMM10Descript % lun                  = 307
      SwanTMM10Descript % filepos              => SWAN_TMM10_POS
      SwanTMM10Descript % initial_value        = 0.0
      SwanTMM10Descript % num_items_per_record = 1
      SwanTMM10Descript % num_fd_records       = NP_G
      SwanTMM10Descript % num_records_this     = NP
      SwanTMM10Descript % imap                 => NODES_LG
      SwanTMM10Descript % array                => SWAN_TMM10Out
      SwanTMM10Descript % array_g              => SWAN_TMM10Out_g
      SwanTMM10Descript % ConsiderWetDry       = .FALSE.
      SwanTMM10Descript % alternate_value      = -99999.0
      SwanTMM10Descript % file_basename        = "swan_TMM10"
      SwanTMM10Descript % file_extension       = 63
      SwanTMM10Descript % field_name           = "swan_TMM10"
      SwanTMM10Descript % startTimeStep        = NTCYSGW
      SwanTMM10Descript % endTimeStep          = NTCYFGW
      SwanTMM10Descript % outputTimeStepIncrement = NSPOOLGW
      SwanTMM10Descript % spoolCounter         => SWAN_TMM10_SPOOL
      SwanTMM10Descript % writerFormats(1:4)   = (/ 1, 3, 4, 5 /)
      SwanTMM10Descript % writeFlag            = .true.

      IF(SWAN_OutputTMM10)THEN
          SwanTMM10MaxDescript % specifier        = NOUTGW
      ELSE
          SwanTMM10MaxDescript % specifier        = OFF
      ENDIF
      SwanTMM10MaxDescript % lun                  = 322
      SwanTMM10MaxDescript % initial_value        = 0.0
      SwanTMM10MaxDescript % num_items_per_record = 1
      SwanTMM10MaxDescript % num_fd_records       = NP_G
      SwanTMM10MaxDescript % num_records_this     = NP
      SwanTMM10MaxDescript % imap                 => NODES_LG
      SwanTMM10MaxDescript % array               => SWAN_TMM10MaxOut
      SwanTMM10MaxDescript % array_g             => SWAN_TMM10MaxOut_g
      SwanTMM10MaxDescript % ConsiderWetDry       = .FALSE.
      SwanTMM10MaxDescript % alternate_value      = -99999.0
      SwanTMM10MaxDescript % file_basename        = "swan_TMM10_max"
      SwanTMM10MaxDescript % file_extension       = 63
      SwanTMM10MaxDescript % field_name           = "swan_TMM10_max"
      SwanTMM10MaxDescript % writeFlag            = .false.


#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      call unsetMessageSource()
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE initSwanOutput
C-----------------------------------------------------------------------
#endif

C-----------------------------------------------------------------------
C      S U B R O U T I N E   I N I T   I N U N D A T I O N   O U T P U T
C-----------------------------------------------------------------------
C     jgf52.08.01: Initialize output variables and arrays for detailed
C     inundation output. These variables must be initialized prior
C     to the call to initOutput2D() because the initiallyDry array
C     depends on the coldstart values in the nnodecode array, which
C     are not available after the call to hstart().
C
C     The initiallydry.63 file must be written during the initialization
C     process (coldstart and hotstart), rather than waiting for the first
C     call to writeOutput2D().
C-----------------------------------------------------------------------
      subroutine initInundationOutput()
      use sizes, only : mnp, globaldir, localdir, write_local_files,
     &    OFF, ASCII, NETCDF3, NETCDF4, XDMF, mnproc, myproc, numFormats
      use global, only : nnodecode, np_g, OutputDataDescript_t, noutge,
     &    fauxFilePosition, nodes_lg
      use global_io, only : writeDomainHeader, packOne, unpackOne
      use mesh, only : np
#ifdef ADCNETCDF
      use netcdfio, only : initNetCDFOutputFile
#endif
      implicit none
      logical :: nerr = .false. !.true. if there was an issue with netcdf init
      call setMessageSource("initInundationOutput")
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      ! jgf52.08.01: Additions for detailed inundation output.
      allocate(initiallyDry(mnp))
      allocate(endRisingInun(mnp))
      allocate(maxInunDepth(mnp))
      allocate(maxInunDepth_time(mnp))
      allocate(penultimateNodeCode(mnp))
      allocate(inundationTime(mnp))
      allocate(inundationTime_onset(mnp))
      allocate(everDried(mnp))
      allocate(everDried_time(mnp))

      initiallyDry(:) = 0

      where (nnodecode.eq.0)
         initiallyDry = 1
      end where
      endRisingInun(:) = 0
      maxInunDepth(:) = -99999.d0
      maxInunDepth_time(:) = 0.d0
      penultimateNodeCode(:) = 0
      inundationTime(:) = 0.d0
      inundationTime_onset(:) = -99999.d0
      everDried(:) = 1.d0
      everDried_time = -99999.d0
      !
      ! initiallydry.63 (nodes that are deemed dry upon cold start)
      IF ( (MNPROC.gt.1) .and. (MyProc.eq.0) ) THEN
         allocate(initiallyDry_g(NP_G))
         initiallyDry_g(:) = 0  ! all wet by default
      ENDIF
      InitiallyDryDescript % lun                  = 402
      InitiallyDryDescript % specifier            = NOUTGE
      InitiallyDryDescript % iarray               => initiallyDry
      InitiallyDryDescript % iarray_g             => initiallyDry_g
      InitiallyDryDescript % field_name           = 'initiallyDry'
      InitiallyDryDescript % file_basename        = 'initiallydry'
      InitiallyDryDescript % file_extension       = 63
      InitiallyDryDescript % filepos              => fauxFilePosition
      InitiallyDryDescript % isInteger            = .true.
      InitiallyDryDescript % writeFlag            = .false.
      InitiallyDryDescript % initial_value        = 0.0
      InitiallyDryDescript % num_items_per_record = 1
      InitiallyDryDescript % num_fd_records       = np_g
      InitiallyDryDescript % num_records_this     = np
      InitiallyDryDescript % imap                 => nodes_lg
      InitiallyDryDescript % considerWetDry       = .false.
      InitiallyDryDescript % alternate_value      = -99999.0
      InitiallyDryDescript % isStation            = .false.
      InitiallyDryDescript % divideByDepth        = .false.
      allocate(InitiallyDryDescript % writerFormats(numFormats))
      InitiallyDryDescript % writerFormats(:)     = -99999
      InitiallyDryDescript % useWriter            = .false.
      InitiallyDryDescript % initialized          = .false.
      InitiallyDryDescript % minmax_timestamp     = .false.
      InitiallyDryDescript % readMaxMin           = .false.

      if ( WRITE_LOCAL_FILES.eqv..true. ) then
         InitiallyDryDescript % file_name =
     &       trim(localdir)//'/'//'initiallydry.63'
      else
         InitiallyDryDescript % file_name =
     &        trim(globaldir)//'/'//'initiallydry.63'
      endif
      if ( InitiallyDryDescript % specifier .eq. XDMF ) then
         InitiallyDryDescript % specifier = ASCII
      endif
      !
      select case(abs(InitiallyDryDescript % specifier))
      case(OFF)
         ! don't need to initialize or write anything
      case(ASCII)
         if ( WRITE_LOCAL_FILES.eqv..true. ) then
            call writeDomainHeader(402,
     &      trim(localdir)//'/'//'initiallydry.63',
     &             NP_G, NP, 'initDry   ')
         else
            if (myproc.eq.0) then
               call writeDomainHeader(402,
     &         trim(globaldir)//'/'//'initiallydry.63',
     &             NP_G, NP, 'initDry   ')
            endif
         endif
#ifdef ADCNETCDF
      case(NETCDF3,NETCDF4)
         ! jgf52.08.05: Only processor 0 needs to initialize the file.
         if (myproc.eq.0) then
            call initNetCDFOutputFile(InitiallyDryDescript, nerr)
            if (nerr.eqv..true.) then
               call allMessage(ERROR,
     &      'There was an issue with NetCDF file initialization for "'//
     &        trim(InitiallyDryDescript%file_name)//'".')
               call terminate()
            endif
         endif
         InitiallyDryDescript % initialized = .true.
#endif
      case default
         call allMessage(ERROR,'The initiallydry.63 file cannot '//
     &      'be written in the specified format.')
      end select
      !
      ! Now write the data rather than waiting for the call to
      ! writeOutput2D() at the end of the timestepping loop.
      call writeOutArray(0.d0, 1, InitiallyDryDescript,
     &     packOne, unpackOne)

#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      end subroutine initInundationOutput
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C                     S U B R O U T I N E
C         C O L L E C T   I N U N D A T I O N   D A T A
C-----------------------------------------------------------------------
C     jgf52.08.01: Collect data used for inundation output.
C-----------------------------------------------------------------------
      subroutine collectInundationData(timeloc, it)
      use sizes, only : mnp
      use global, only : nodecode, dtdp, eta1, eta2, nt, inunThresh
      use mesh, only : dp
      implicit none
      real(8), intent(in) :: timeloc ! time in s since cold start
      integer, intent(in) :: it     ! timesteps since cold start

      call setMessageSource("collectInundationData")
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      !
      ! Accumulate total time inundated in initally dry areas
      ! that are currently wet and have a total water depth greater
      ! than the specified threshold.
      !
      ! jgf52.08.03: The original code I received compared eta2 (water
      ! surface elevation relative to the mesh datum) to the inundation
      ! threshold (e.g., inunThresh=0.6), but I don't think that
      ! makes sense. So I changed this code so it compares the inundation
      ! threshold to the total water depth.
      where ( (initiallyDry.eq.1).and.(nodecode.eq.1).and.
     &            ((eta2+dp).ge.inunThresh) )
         inundationTime = inundationTime + dtdp
      end where
      !
      ! record the time of onset of inundation beyond the threshold
      where ( (initiallyDry.eq.1).and.(nodecode.eq.1).and.
     &         ((eta2+dp).ge.inunThresh).and.
     &         (inundationTime_onset.lt.0.d0) )
            inundationTime_onset = timeloc
      end where
      !
      ! jgf52.08.01: I don't understand why we are accumulating
      ! inundation time in normally wet areas but I am keeping this
      ! here for backward compatibility.
      where ( initiallyDry.eq.0 )
         inundationTime = inundationTime + dtdp
         inundationTime_onset = 0.d0
      end where
      !
      ! Record the maximum depth of inundation in initially dry areas
      ! that are currently wet with a total water depth greater than the
      ! threshold value.
      where ( (initiallyDry.eq.1).and.(nodecode.eq.1).and.
     &             ((eta2+dp).gt.maxInunDepth) )
         maxInunDepth = eta2 + dp
         ! record time of occurrence for the maximum inundation
         maxInunDepth_time = timeloc
      end where
      !
      ! On the next to last time step, record the wet/dry state
      if (it.eq.nt-1) then
         endRisingInun(:) = 0
         penultimateNodeCode(:) = nodecode(:)
      endif
      !
      ! For nodes that were initially dry, but wet on the next-to-last
      ! time step, and are experiencing increasing water level on the
      ! last time step set a flag.
      ! jgf52.08.03: Added a check to make sure the node is also wet
      ! on the very last time step, just in case.
      if (it.eq.nt) then
         where ( (initiallyDry.eq.1).and.(penultimateNodeCode.eq.1)
     &           .and. (nodecode.eq.1) .and. (eta2.gt.eta1) )
            endRisingInun = 1
         end where
      endif
      !
      ! jgf52.08.15: Adding everDried data set for use with keeping
      ! track of which nodes have dried because this affects the
      ! computation of harmonic analysis results
      where ( nodecode.eq.0 )
         everDried = -99999.d0
         ! record the total amount of time this node was dry during the
         ! simulation
         everDried_time = everDried_time + dtdp
      end where

#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      end subroutine collectInundationData
C-----------------------------------------------------------------------



C-----------------------------------------------------------------------
C                     S U B R O U T I N E
C            C O L L E C T   M I N   M A X   D A T A
C-----------------------------------------------------------------------
C     jgf52.08.01: Collect data used for min/max output. Removed
C     from subroutine timestep() and placed here.
C-----------------------------------------------------------------------
      subroutine collectMinMaxData(timeloc)
      use global, only : nodecode, eta2, uu2, vv2, pr2,
     &  rsnx2, rsny2, etamax, etamax_time, umax, umax_time, prmin,
     &  prmin_time, wvnoutmax, wvnoutmax_time, rsnmax, rsnmax_time,
     &  wvnxout, wvnyout
      use mesh, only : np
      implicit none
      real(8), intent(in) :: timeloc ! time in s since cold start
      real(sz), allocatable, save :: minMaxTemp(:) ! temp value for min/max comparisons
      integer :: i

      call setMessageSource("collectMinMaxData")
#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ! jgf52.08.20: Allocate space for temp calculations.
      if (allocated(minMaxTemp).eqv..false.) then
         allocate(minMaxTemp(np))
      endif

C...  jgf50.44: Collect time that min/max occurs as well.
      where ( (nodecode.eq.1).and.(eta2.gt.etamax) )
         etamax = eta2
         etamax_time = timeLoc !tcm v51.20.01
      end where

      minMaxTemp = uu2**2+vv2**2
      where ( minMaxTemp.gt.umax**2 )
         umax = sqrt(uu2**2+vv2**2)
         umax_time = timeLoc   !tcm v51.20.01
      end where

      where (pr2.lt.prmin)
         prmin = pr2
         prmin_time = timeLoc !tcm v51.20.01
      end where

      minMaxTemp = wvnxout**2+wvnyout**2
      where (minMaxTemp.gt.wvnoutmax**2)
         wvnoutmax = sqrt(wvnxout**2+wvnyout**2)
         wvnoutmax_time = timeLoc
      end where

      minMaxTemp = rsnx2**2+rsny2**2
      where (minMaxTemp.gt.rsnmax**2)
         rsnmax = sqrt(rsnx2**2+rsny2**2)
         rsnmax_time = timeLoc
      end where

#if defined(WRITE_OUTPUT_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      end subroutine collectMinMaxData
C-----------------------------------------------------------------------


      !-----------------------------------------------------------------
      end module write_output
      !-----------------------------------------------------------------
