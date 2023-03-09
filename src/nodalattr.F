C******************************************************************************
C     PADCIRC VERSION 46.00 xx/xx/2006
C     last changes in this file VERSION 46.00
C
C     Written for ADCIRC v46.00 by Jason G. Fleming.
C******************************************************************************
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     M O D U L E   N O D A L  A T T R I B U T E S
C-----------------------------------------------------------------------
C
C     jgf46.00 This module manages nodal attribute data, including
C     bottom friction, tau0, startdry, directional wind speed reduction,
C     and etc. Will read the Nodal Attributes File (unit 13) and
C     initialize the nodal attribute arrays.
C
C     Handling data by label rather than an integer encoding should
C     result in increased transparency as well as ease the transition to
C     HDF5/NetCDF i/o. The labels were chosen according to the
C     guidelines of the CF Standard. Creating labels according to CF
C     Standard Guidelines should enhance interoperability with other
C     simulation frameworks.
C
C     To use a nodal attribute contained in the fort.13 file, the
C     corresponding attribute name must appear in the fort.15 file.  A
C     list of nodal attributes is read in from the fort.15 file if the
C     fort.15 parameter NWP > 0. This also signals ADCIRC to look for a
C     fort.13 file.
C
C     Summary of the file format for the Nodal Attributes File:
C
C     AGRID                       ! user's comment line - should be a cross
C                                 !    reference to the grid file
C     NumOfNodes                  ! number of nodes, must match NP
C                                 !    from grid file
C     NAttr                       ! number of attributes contained in this file
C
C     do i=1, NAttr
C        AttrName(i)              ! nodal attribute name (see
C                                 !    valid names below)
C        Units(i)                 ! physical units (ft, m/s, none)
C        ValuesPerNode(i)         ! number of values at each node for
C                                 !   a particular attribute
C        DefaultVal(i)            ! default value(s) for the nodal attribute
C     end do
C
C     do i=1, NAttr
C        AttrName(i)              ! label of the attribute, again
C        NumNodesNotDefaultVal(i) ! number of nodes with non-default values
C        do j=1, NumNodesNotDefault(i)
C           n, (Attr(n,k), k=1,ValuesPerNode(i))
C        end do
C     end do
C
C
C
C     Valid labels are as follows:
C
C     ADCIRC Variable:       CF-Style Label:
C      Tau0                  "primitive_weighting_in_continuity_equation"
C      StartDry              "surface_submergence_state"
C      Fric                  "quadratic_friction_coefficient_at_sea_floor"
C      z0Land                "surface_directional_effective_roughness_length"
C                            (z0Land has ValuesPerNode = 12)
C      VCanopy               "surface_canopy_coefficient"
C      BK,BAlpha,BDelX,POAN  "bridge_pilings_friction_parameters"
C                            (bridge_pilings... has ValuesPerNode=4)
C      ManningsN             "mannings_n_at_sea_floor"
C      Chezy                 "chezy_friction_coefficient_at_sea_floor"
C      Z0b_var               "bottom_roughness_length"
C      GeoidOffset           "sea_surface_height_above_geoid"
C      EVM                   "average_horizontal_eddy_viscosity_in_sea_water_wrt_depth"
C      EVC                   "average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth"
C      Tao0MinMax            "min_and_max_primitive_weighting_in_continuity_equation"
C      River_et_WSE          "initial_river_elevation"
C      IT_Fric               "internal_tide_friction"
C      OverlandReductionFactor  "overland_reduction_factor"
C      subgridBarrier        "subgrid_barrier"
C      CondensedNodes        "condensed_nodes"
C
C-----------------------------------------------------------------------
      MODULE NodalAttributes
      USE SIZES
      USE GLOBAL, ONLY : scratchMessage, allMessage, logMessage,
     &                    DEBUG, ECHO, INFO, WARNING, ERROR

Casey 100210: Allow SWAN to handle wave refraction as a nodal attribute.
C             I've placed these changes outside the #ifdef CSWAN flags
C             because we want to be able to use the same fort.13 files
C             for both ADCIRC and SWAN+ADCIRC runs.  This way, the new
C             nodal attribute will be processed but only applied when
C             ADCIRC is coupled to SWAN.
      LOGICAL              :: LoadSwanWaveRefrac
      LOGICAL              :: FoundSwanWaveRefrac
      CHARACTER(LEN=80)    :: SwanWaveRefracUnits
      INTEGER              :: SwanWaveRefracNoOfVals
      REAL(SZ)             :: SwanWaveRefracDefVal
      REAL(SZ),ALLOCATABLE :: SwanWaveRefrac(:)
C
Corbitt 120321: Allow Advection to be Turned on Locally instead of Globally
      LOGICAL              :: LoadAdvectionState
      LOGICAL              :: FoundAdvectionState
      CHARACTER(LEN=80)    :: AdvectionStateUnits
      INTEGER              :: AdvectionStateNoOfVals
      REAL(SZ)             :: AdvectionStateDefVal
      REAL(SZ),ALLOCATABLE :: AdvectionState(:)
Csb   Specify node groups for 1D condensation.
      LOGICAL              :: LoadCondensedNodes
      LOGICAL              :: FoundCondensedNodes
      CHARACTER(LEN=80)    :: CondensedNodesUnits
      INTEGER              :: CondensedNodesNoOfVals
      INTEGER              :: CondensedNodesDefVal(4)
      INTEGER,ALLOCATABLE  :: CondensedNodes(:,:)
      INTEGER,ALLOCATABLE  :: NCondensedNodes(:)
      INTEGER,ALLOCATABLE  :: ListCondensedNodes(:,:)
      INTEGER              :: NListCondensedNodes
      INTEGER,ALLOCATABLE  :: NNodesListCondensedNodes(:)
C
C     The following flags are .true. if the corresponding data are
C     required for the run, according to the unit 15 control file
      LOGICAL LoadTau0
      LOGICAL LoadStartDry
      LOGICAL LoadDirEffRLen
      LOGICAL LoadCanopyCoef
      LOGICAL LoadQuadraticFric
      LOGICAL LoadBridgePilings
      LOGICAL LoadChezy
      LOGICAL LoadManningsN
      LOGICAL LoadGeoidOffset
      LOGICAL LoadEVM
      LOGICAL LoadEVC
      LOGICAL LoadTau0MinMax
      LOGICAL LoadZ0b_var
      LOGICAL LoadEleSlopeLim ! zc: elemental slope limiter
      LOGICAL LoadRiver_et_WSE  ! tcm 20140502 v51.27
      LOGICAL LoadIT_Fric       ! wjp: Internal tides
      LOGICAL LoadOverlandReductionFactor 
      LOGICAL LoadSubgridBarrier
C
C     The following flags are .true. if there are data with the
C     corresponding label in the unit 13 file.
      LOGICAL FoundTau0
      LOGICAL FoundStartDry
      LOGICAL FoundDirEffRLen
      LOGICAL FoundCanopyCoef
      LOGICAL FoundQuadraticFric
      LOGICAL FoundBridgePilings
      LOGICAL FoundChezy
      LOGICAL FoundManningsN
      LOGICAL FoundGeoidOffset
      LOGICAL FoundEVM
      LOGICAL FoundEVC
      LOGICAL FoundTau0MinMax
      LOGICAL FoundZ0b_var
      LOGICAL FoundEleSlopeLim
      LOGICAL FoundRiver_et_WSE ! tcm 20140502 v51.27      
      LOGICAL FoundIT_Fric  ! wjp internal tides
      LOGICAL FoundOverlandReductionFactor
      LOGICAL FoundSubgridBarrier
C
C     These variables hold the strings which describe the attribute's
C     units. These data are loaded from the file, but not used as of
C     v46.00.
      CHARACTER(len=80) Tau0Units
      CHARACTER(len=80) StartDryUnits
      CHARACTER(len=80) DirEffRLenUnits
      CHARACTER(len=80) CanopyCoefUnits
      CHARACTER(len=80) QuadraticFricUnits
      CHARACTER(len=80) BridgePilingsUnits
      CHARACTER(len=80) ChezyUnits
      CHARACTER(len=80) ManningsNUnits
      CHARACTER(len=80) GeoidOffsetUnits
      CHARACTER(len=80) EVMUnits
      CHARACTER(len=80) EVCUnits
      CHARACTER(len=80) Tau0MinMaxUnits
      CHARACTER(len=80) Z0b_varUnits
      CHARACTER(len=80) EleSlopeLimUnits
      CHARACTER(len=80) River_et_WSEUnits ! tcm 20140502 v51.27
      CHARACTER(len=80) IT_FricUnits  ! wjp internal tides
      CHARACTER(len=80) OverlandReductionFactorUnits
      CHARACTER(len=80) subgridBarrierUnits
C
C     These variables hold the number of values per node for each
C     attribute.
      INTEGER Tau0NoOfVals
      INTEGER StartDryNoOfVals
      INTEGER DirEffRLenNoOfVals
      INTEGER CanopyCoefNoOfVals
      INTEGER QuadraticFricNoOfVals
      INTEGER BridgePilingsNoOfVals
      INTEGER ChezyNoOfVals
      INTEGER ManningsNNoOfVals
      INTEGER GeoidOffsetNoOfVals
      INTEGER EVMNoOfVals
      INTEGER EVCNoOfVals
      INTEGER Tau0MinMaxNoOfVals
      INTEGER Z0b_varNoOfVals
      INTEGER EleSlopeLimNoofVals
      INTEGER River_et_WSENoOfVals  ! tcm 20140502 v51.27      
      INTEGER IT_FricNoOfVals ! wjp internal tides
      INTEGER OverlandReductionFactorNoOfVals
      INTEGER subgridBarrierNoOfVals
C
C     These variables hold the default values for each attribute.
      REAL(SZ) Tau0DefVal
      REAL(SZ) StartDryDefVal
      REAL(SZ) DirEffRLenDefVal(12)
      REAL(SZ) CanopyCoefDefVal
      REAL(SZ) QuadraticFricDefVal
      REAL(SZ) BridgePilingsDefVal(4)
      REAL(SZ) ChezyDefVal
      REAL(SZ) ManningsNDefVal
      REAL(SZ) GeoidOffsetDefVal
      REAL(SZ) EVMDefVal
      REAL(SZ) EVCDefVal
      REAL(SZ) Tau0MinMaxDefVal(2)
      REAL(SZ) Z0b_varDefVal
      REAL(SZ) EleSlopeLimDefVal
      REAL(SZ) River_et_WSEDefVal ! tcm 20140502 v51.27
      REAL(SZ) IT_FricDefVal(3)  ! wjp internal tides
      REAL(SZ) OverlandReductionFactorDefaultValue
      REAL(SZ) SubgridBarrierDefaultValue(2)
C
      INTEGER NumOfNodes    ! number of nodes listed in unit 13 file, cf. NP
      INTEGER NAttr         ! number of nodal attributes in the unit 13 file
C
C     The following variables are inputs from the unit 15 model param. file
      INTEGER NWP     ! number of nodal attributes to read from file
      INTEGER NoLiBF  ! nonlinear bottom friction indicator
      REAL(SZ) Tau0   ! primitive continuity eqn. weight
      REAL(SZ) Tau    ! linear friction coefficient (1/sec)
      REAL(SZ) CF     ! 2DDI bottom fric. coef., effect varies based on NoLiBF
      REAL(SZ) HBreak ! break depth for NOLIBF .eq. 2
      REAL(SZ) FTheta ! dimless param. for NOLIBF .eq. 2
      REAL(SZ) FGamma ! dimless param. for NOLIBF .eq.
      REAL(SZ) ESLM   ! horizontal eddy viscosity (length^2/time)
      REAL(SZ) ESLC   ! horizontal eddy diffusivity (length^2/time)
      INTEGER  IFLINBF! flag to turn on linear bottom friction
      INTEGER  IFNLBF ! flag to turn on nonlinear bottom friction
      INTEGER  IFHYBF ! flag to turn on hybrid bottom friction
      REAL(SZ) BFCdLLimit ! lower limit of quadratic bottom friction
      REAL(SZ) Tau0FullDomainMin ! lower limit of tau0 if time varying
      REAL(SZ) Tau0FullDomainMax ! upper limit of tau0 if time varying
C     jgf48.42 Used to control the back-loaded time averaged tau0
      REAL(SZ), PARAMETER :: AlphaTau0 = 0.25d0
C
C     Nodal attributes.
      REAL(SZ), ALLOCATABLE :: STARTDRY(:) ! 1=nodes below geoid initially dry
      REAL(SZ), ALLOCATABLE :: FRIC(:)     ! bottom friction coefficient
      REAL(SZ), ALLOCATABLE :: FRIC_ORIG(:)
      REAL(SZ), ALLOCATABLE, TARGET :: TAU0VAR(:)! primitive equation weighting
      REAL(SZ), ALLOCATABLE :: Tau0Temp(:) ! used in time varying tau0
C     jjw&sb46.39.sb01: Base (original) primitive equation weighting.
C     Tau0Var may be optimized later based on Tau0Base
      REAL(SZ), ALLOCATABLE :: TAU0BASE(:)
C     jgf47.33 Used for time averaged tau0
      REAL(SZ), ALLOCATABLE :: LastTau0(:)

C     jgf47.06: Added variables to trigger calculations and output of tau0
C     jgf47.30: Added "FullDomain" or "High Res Areas Only" distinction
C     jgf47.31: Added time averaging to time varying tau0
C     jgf48.42: Added backloaded time averaged tau0
      LOGICAL :: HighResTimeVaryingTau0    ! .true. if Tau0 == -3.x in fort.15
      LOGICAL :: FullDomainTimeVaryingTau0 ! .true. if Tau0 == -5.x in fort.15
      LOGICAL :: OutputTau0      ! .true. if Tau0Dig2 == -1 in fort.15
      LOGICAL :: TimeAveragedTau0 ! .true. if Tau0 == -6.x in fort.15
      LOGICAL :: BackLoadedTimeAveragedTau0 ! .true. if Tau0 == -7.x
      LOGICAL,ALLOCATABLE :: elemental_slope_limiter_active(:) ! .true. if elemental_slope_limiter_grad_max
                                                               ! has been exceeded, initialized as .false.
      LOGICAL,ALLOCATABLE :: elemental_slope_limiter_max_exceeded(:) ! .true. if maximum gradient has been
                                                                   ! exceeded, used for gradient warnings

      REAL(SZ), ALLOCATABLE :: z0land(:,:) ! directional wind speed red. fac.
      REAL(SZ), ALLOCATABLE :: vcanopy(:)  ! canopy coefficient
C     The following attribute contains BK(I),BALPHA(I),BDELX(I), and POAN(I)
      REAL(SZ), ALLOCATABLE :: BridgePilings(:,:)
      REAL(SZ), ALLOCATABLE :: Chezy(:)
      REAL(SZ), ALLOCATABLE :: ManningsN(:)
      REAL(SZ), ALLOCATABLE :: GeoidOffset(:)
      REAL(SZ), ALLOCATABLE :: EVM(:)
      REAL(SZ), ALLOCATABLE :: EVC(:)
      REAL(SZ), ALLOCATABLE :: Tau0MinMax(:,:)  ! (node,i); i=1(min), i=2(max)
      REAL(SZ), ALLOCATABLE :: Z0b_var(:)  ! patially varying 3D bottom roughness length
      REAL(SZ), ALLOCATABLE :: elemental_slope_limiter_grad_max(:)
      REAL(SZ), ALLOCATABLE :: River_et_WSE(:)  ! tcm 20140502 v51.27
      REAL(SZ), ALLOCATABLE :: IT_Fric(:,:)  ! WJP: Linear friction coefficient for internal tides
      REAL(SZ), ALLOCATABLE :: OverlandReductionFactor(:)
      REAL(SZ), ALLOCATABLE :: subgridBarrier(:,:)
      REAL(SZ), ALLOCATABLE :: subgridBarrierFriction(:)
      LOGICAL,ALLOCATABLE   :: subgridBarrierOvertopping(:)
C
C------------------ For absorbing layer, DW
      REAL(SZ), pointer:: absorblayer_sigma(:,:), absorblayer_zero(:,:)

      REAL(SZ), pointer:: absorblayer_sigma_eta(:,:)
      REAL(SZ), pointer:: absorblayer_sigma_mnx(:,:)
      REAL(SZ), pointer:: absorblayer_sigma_mny(:,:)


      LOGICAL LoadAbsLayerSigma
      LOGICAL FoundAbsLayerSigma

      INTEGER AbsLayerSigmaNoOfVals  
      REAL(SZ) AbsLayerSigmaDefVal(2)
      
      CHARACTER (LEN=80):: AbsLayerSigmaUnits

      INTEGER:: NumNodesAbsLayer(3)
      INTEGER, pointer:: AbsLayerNodesID(:)

      INTEGER, pointer:: AbsLayer_Eta_NodesID(:)
      INTEGER, pointer:: AbsLayer_MnX_NodesID(:)
      INTEGER, pointer:: AbsLayer_MnY_NodesID(:)

      INTEGER:: AbsLayerType !  0 - Absorbing BCs,  1 - Harmonic forcing (Tides)  !
                             ! -1 - Read from file                                !

! For wave generator layer of type 1
      INTEGER:: AbsLayerNBF

      ! Water elevation
      REAL(SZ), ALLOCATABLE:: AbsLayerEtaAMIG(:), 
     &                AbsLayerEtaFF(:), AbsLayerEtaFACE(:)
      REAL(SZ), ALLOCATABLE:: AbsLayerEtaEMO(:,:), AbsLayerEtaEFA(:,:)

      ! Tentatively for Velocity or Dischage
      REAL(SZ), ALLOCATABLE:: AbsLayerQAMIG(:), 
     &                AbsLayerQFF(:), AbsLayerQFACE(:)
      REAL(SZ), ALLOCATABLE:: AbsLayerQxEMO(:,:), AbsLayerQxEFA(:,:)
      REAL(SZ), ALLOCATABLE:: AbsLayerQyEMO(:,:), AbsLayerQyEFA(:,:)

C
C.....CPB 9/2022: for computing tidal velocities in
C                 apply2dinternalwavedrag
      REAL(SZ), ALLOCATABLE, PRIVATE :: UT(:), VT(:)

C------------------
C     Make i,j,k private DW
      INTEGER,private :: i        ! node loop counter
      INTEGER,private :: j        ! attribute values loop counter
      INTEGER,private :: k        ! attribute loop counter
C END DW
      !
      !  X D M F
      !
      type nodalAttr_t
         character(len=1024) :: attrName ! name of the nodal attr 
         character(len=1024) :: units    ! physical units of the nodal attr
         integer :: numVals  ! number of values at each node for this nodal attr
         integer :: numNodesNotDefault ! number of nodes with values different from the default
         real(8), allocatable :: defaultVals(:) ! default value(s) for real valued attributes 
         real(8), allocatable :: nonDefaultVals(:,:) ! nondefault value(s) for real valued attributes 
         integer, allocatable :: nonDefaultNodes(:) ! node numbers where nondefault vals occur
         real(8), allocatable :: xdmfArray(:)
         real(8), allocatable :: xdmfMatrix(:,:)
      end type nodalAttr_t
      ! variable capable of holding all nodal attributes in the file
      type(nodalAttr_t), allocatable :: na(:)
      character(len=1024) :: nodalAttributesComment ! comment line at the top
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CONTAINS !- - - - - - - - - - - - - - - - - - - - - - - - - - - -
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

C     ----------------------------------------------------------------
C     S U B R O U T I N E     I N I T  N  A  M O D U L E
C     ----------------------------------------------------------------
C
C     jgf46.00 Subroutine to initialize the variables in the nodal
C     attributes module.
C
C     ----------------------------------------------------------------
      SUBROUTINE InitNAModule()
      IMPLICIT NONE
C
Casey 100210: Make changes compact.
      LoadSwanWaveRefrac     = .FALSE.
      FoundSwanWaveRefrac    = .FALSE.
      SwanWaveRefracNoOfVals = 1
      SwanWaveRefracDefVal   = 0.D0
Corbitt 120321:
      LoadAdvectionState     = .False.
      FoundAdvectionState    = .FALSE.
      AdvectionStateNoOfVals = 1
      AdvectionStateDefVal   = 0.D0
Csb
      LoadCondensedNodes     = .False.
      FoundCondensedNodes    = .FALSE.
      CondensedNodesNoOfVals = 3
      CondensedNodesDefVal   = 0
C
      LoadTau0           = .False.
      LoadStartDry       = .False.
      LoadDirEffRLen     = .False.
      LoadManningsN      = .False.
      LoadQuadraticFric  = .False.
      LoadChezy          = .False.
      LoadBridgePilings  = .False.
      LoadCanopyCoef     = .False.
      LoadGeoidOffset    = .False.
      LoadEVM            = .False.
      LoadEVC            = .False.
      LoadTau0MinMax     = .False.
      LoadZ0b_var        = .False.
      LoadEleSlopeLim    = .False.
      LoadRiver_et_WSE   = .False.
      LoadIT_Fric        = .False.  ! WJP
      LoadOverlandReductionFactor = .False.
      LoadSubgridBarrier = .False.
C
      FoundTau0           = .False.
      FoundStartDry       = .False.
      FoundDirEffRLen     = .False.
      FoundManningsN      = .False.
      FoundQuadraticFric  = .False.
      FoundChezy          = .False.
      FoundBridgePilings  = .False.
      FoundCanopyCoef     = .False.
      FoundGeoidOffset    = .False.
      FoundEVM            = .False.
      FoundEVC            = .False.
      FoundTau0MinMax     = .False.
      FoundZ0b_var        = .False.
      FoundEleSlopeLim    = .False.
      FoundRiver_et_WSE   = .False.      
      FoundIT_Fric        = .False. ! WJP
      FoundOverlandReductionFactor = .False.
      FoundSubgridBarrier = .False.
C
C
      Tau0NoOfVals          = 1
      StartDryNoOfVals      = 1
      DirEffRLenNoOfVals    = 12
      QuadraticFricNoOfVals = 1
      ChezyNoOfVals         = 1
      ManningsNNoOfVals     = 1
      BridgePilingsNoOfVals = 4
      CanopyCoefNoOfVals    = 1
      GeoidOffsetNoOfVals   = 1
      EVMNoOfVals           = 1
      EVCNoOfVals           = 1
      Tau0MinMaxNoOfVals    = 2
      Z0b_varNoOfVals       = 1
      EleSlopeLimNoOfVals   = 1
      River_et_WSENoOfVals  = 1
      IT_FricNoOfVals       = 1  ! WJP, could be 3 also, so may change
      OverlandReductionFactorNoOfVals = 1
      SubgridBarrierNoOfVals = 1 ! 1 or 2 
C
      Tau0DefVal             = 0.0
      StartDryDefVal         = 0.0
      DO j=1, DirEffRLenNoOfVals
         DirEffRLenDefVal(j) = 0.0
      END DO
      CanopyCoefDefVal       = 1.0 ! jgf49.1001 default is now full wind stress
      QuadraticFricDefVal    = 0.0
      DO j=1, BridgePilingsNoOfVals
         BridgePilingsDefVal(j) = 0.0
      END DO
      ChezyDefVal            = 0.0
      ManningsNDefVal        = 0.0
      GeoidOffsetDefVal      = 0.0
      EVMDefVal              = 0.0
      EVCDefVal              = 0.0
      DO j=1, Tau0MinMaxNoOfVals
         Tau0DefVal = 0.0
      END DO
      Z0b_varDefVal          = 0.001
      EleSlopeLimDefVal      = 0D0
      River_et_WSEDefVal     = 0.d0      
      IT_FricDefVal          = 0.0 ! WP
      SubgridBarrierDefaultValue = 999999D0
C
      HighResTimeVaryingTau0     = .False.
      FullDomainTimeVaryingTau0  = .False.
      OutputTau0                 = .False.
      TimeAveragedTau0           = .False.
      BackLoadedTimeAveragedTau0 = .False.
C
      HBREAK=1.d0
      FTHETA=1.d0
      FGAMMA=1.d0
C
C   kmd48.33bc this resets the ESLM to 0 and if using constant eddy
C              viscosity it eliminates what the user specified in
C              the input file. The InitNAModule originally come before
C              the read_input call but now it appears after the read_input
C              call.
!      ESLM=0.0
      ESLC=0.0

C---------------- DW ------------------
      LoadAbsLayerSigma = .FALSE. ;
      FoundAbsLayerSigma = .FALSE. ;

      NumNodesAbsLayer = 0 ;
      AbsLayerSigmaNoOfVals = 2 ;
      AbsLayerSigmaDefVal = 0.D0 ;

      AbsLayerType = 0 ;
      AbsLayerNBF = 0 ;

      NULLIFY(absorblayer_sigma) ; 
      NULLIFY(absorblayer_sigma_eta) ; 
      NULLIFY(absorblayer_sigma_mnx) ; 
      NULLIFY(absorblayer_sigma_mny) ;
 
      NULLIFY(AbsLayer_Eta_NodesID) ;
      NULLIFY(AbsLayer_MnX_NodesID) ;
      NULLIFY(AbsLayer_MnY_NodesID) ;
C--------------------------------------
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE InitNAModule
C     ----------------------------------------------------------------

C     ----------------------------------------------------------------
C     S U B R O U T I N E    R E A D   N O D A L   A T T R   X D M F 
C     ----------------------------------------------------------------
C     jgf51.21.24: Read the nodal attributes in XDMF format.
C     ----------------------------------------------------------------
      SUBROUTINE readNodalAttrXDMF()
      USE MESH, ONLY : replaceNullsWithSpaces
      IMPLICIT NONE
#ifndef ADCXDMF
      call allMessage(ERROR,
     & 'An XDMF formatted nodal attributes (fort.13) file was specified.')
      call allMessage(ERROR,
     & 'This ADCIRC executable was not compiled with XDMF support.')
      call na_terminate()
#else
      include 'adcirc_Xdmf.f'
      integer*8 :: xdmfFortranObj ! object that receives the data
      !
      integer :: startIndex
      integer :: arrayStride
      integer :: valueStride
      integer :: attributeIndex
      integer :: attributeType
      integer :: informationIndex
      logical :: isNodalAttribute
      integer :: gridIndex = 0
      integer :: numAttributes
      integer :: numInformation
      integer :: numValues
      integer :: openMaps ! 1 if Maps should be opened within another object
      integer :: openAttributes ! 1 if Attributes should be opened within another object 
      integer :: openInformations ! 1 if Informations should be opened within another object
      integer :: openSets ! 1 if sets should be opened within another object
      integer, parameter :: keyLength = 1024
      integer, parameter :: valueLength = 1024
      character(len=keyLength) :: itemKey
      character(len=valueLength) :: itemValue
      integer, parameter :: nameLength = 256
      character(len=nameLength) :: itemName
      character(len=1024) :: defaultValuesString
      real(8), allocatable :: diff(:)
      logical, allocatable :: areDefaultValues(:)
      logical :: naFound 
      integer :: i, j, nattrCount, nonDefaultCount
      !
      startIndex = 0
      arrayStride = 1
      valueStride = 1
      openMaps = 1
      openAttributes = 1
      openInformations = 1
      openSets = 1
      !
      NAFound = .False.
C
C     Determine if the Nodal Attributes File exists.
      INQUIRE(FILE=TRIM(naFileName),EXIST=NAFound)
C
      IF (.not.NAFound) THEN
         write(scratchMessage,'(a)') 'The XDMF nodal attributes file' //
     &       trim(naFileName)//'was not found.'
         call allMessage(ERROR, scratchMessage)
         call na_terminate()
      ENDIF

      write(6,'(a)') 'INFO: Reading data from the "' 
     &    // trim(adjustl(naFileName)) // '" XDMF file.'
      call xdmfinit(xdmfFortranObj)
      call xdmfRead(xdmfFortranObj, trim(adjustl(naFileName))//char(0))
      call xdmfOpenDomainGrid(xdmfFortranObj, XDMF_GRID_TYPE_UNSTRUCTURED, 
     &   gridIndex, openMaps, openAttributes, openInformations, openSets)
      !
      call xdmfRetrieveNumAttributes(xdmfFortranObj, numAttributes)
      nAttr = numAttributes - 1 ! the depth is included so don't count it
      
      allocate(na(nAttr))
      write(6,'("INFO: Grid ",i0," contains ",i0," nodal attributes.")') gridIndex, nAttr
      !
      ! populate the names of the nodal attributes
      nattrCount = 1
      do attributeIndex=0, numAttributes - 1
         call xdmfRetrieveAttributeName(xdmfFortranObj, attributeIndex,   
     &       itemName, nameLength)
         call replaceNullsWithSpaces(itemName)
         if (trim(itemName).eq.'depth') then
            cycle
         endif
         write(6,'("INFO: Grid ",i0," Attribute ",i0," is named ",a)') gridIndex, attributeIndex, trim(itemName)
         na(nattrCount)%attrName = trim(itemName)
         SELECT CASE (trim(na(nattrCount)%AttrName))
         case("depth")
            cycle
         CASE("primitive_weighting_in_continuity_equation")
            FoundTau0 = .True.
         CASE("surface_submergence_state")
            FoundStartDry = .True.
         CASE("quadratic_friction_coefficient_at_sea_floor")
            FoundQuadraticFric = .True.
         CASE("surface_directional_effective_roughness_length")
            FoundDirEffRLen = .True.
         CASE("surface_canopy_coefficient")
            FoundCanopyCoef = .True.
         CASE("bridge_pilings_friction_parameters")
            FoundBridgePilings = .True.
         CASE("mannings_n_at_sea_floor")
            FoundManningsN = .True.
         CASE("bottom_roughness_length")
            FoundZ0b_var = .True.
         CASE("chezy_friction_coefficient_at_sea_floor")
            FoundChezy = .True.
         CASE("sea_surface_height_above_geoid")
            FoundGeoidOffset = .True.
         CASE
     &    ("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")
            FoundEVM = .True.
         CASE
     &   ("average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth")
            FoundEVC = .true.
         CASE
     &   ("min_and_max_primitive_weighting_in_continuity_equation")
            FoundTau0MinMax = .True.
Casey 100210: Allow SWAN to handle wave refraction as a nodal attribute.
         CASE("wave_refraction_in_swan")
            FoundSwanWaveRefrac = .TRUE.
Corbitt 120321: Allow advection to be turned on locally instead of globally
         CASE("advection_state")
            FoundAdvectionState = .TRUE.
Csb
         CASE("condensed_nodes")
            FoundCondensedNodes = .TRUE.
         CASE("elemental_slope_limiter")
            FoundEleSlopeLim = .TRUE.
         CASE("initial_river_elevation")            
            FoundRiver_et_WSE     = .TRUE.
         CASE("overland_reduction_factor")
            FoundOverlandReductionFactor = .TRUE.
         CASE DEFAULT
            scratchMessage = "Unrecognized nodal attribute : " //
     &         trim(na(nattrCount) % attrName)
            call allMessage(WARNING, scratchMessage)
         END SELECT
         nattrCount = nattrCount + 1  
      END DO
C
C     Determine if there are any attributes required by the fort.15 file
C     that are not in the nodal attributes file.
      call checkForMissingNodalAttributes()
      !
      ! use the names to populate the number of values for each attribute
      call xdmfRetrieveNumInformation(xdmfFortranObj, numInformation)
      write(6,'(a,i0)') 'numInformation=',numInformation
      do informationIndex=0,numInformation-1
         call xdmfRetrieveInformation(xdmfFortranObj, informationIndex, itemKey, keyLength, itemValue, valueLength)
         call replaceNullsWithSpaces(itemKey)
         call replaceNullsWithSpaces(itemValue)
         select case(trim(itemKey))
         case('nodalAttributesComment')
            nodalAttributesComment = trim(itemValue)
         case('numMeshNodes')
            read(itemValue,*) numOfNodes
         case default
            do i=1,nAttr
               if (trim(itemKey).eq. trim(na(i)%attrName) // ' number_of_values') then
                  read(itemValue,*) na(i)%numVals
                  allocate(na(i)%defaultVals(na(i)%numVals))
                  exit
               endif
            end do
         end select
      end do
C
C     Allocate memory for nodal attributes now that we have the
C     numOfNodes value
      call allocateNodalAttributes()
      !
      ! populate the units and the default values
      call xdmfRetrieveNumInformation(xdmfFortranObj, numInformation)
      do informationIndex=0,numInformation-1
         call xdmfRetrieveInformation(xdmfFortranObj, informationIndex, itemKey, keyLength, itemValue, valueLength)
         call replaceNullsWithSpaces(itemKey)
         call replaceNullsWithSpaces(itemValue)
         do i=1,nAttr
            if (trim(itemKey).eq. trim(na(i)%attrName)// ' units') then
               na(i) % units = trim(itemValue)
            endif  
            if (trim(itemKey).eq. trim(na(i)%attrName) // ' default_values') then
               write(6,*) 'default_values'//trim(itemValue)
               read(itemValue,*) (na(i)%defaultVals(j),j=1,na(i)%numVals)
            endif
         end do
      end do
      !
      ! populate the data
      do attributeIndex=0, numAttributes - 1
         call xdmfRetrieveAttributeName(xdmfFortranObj, attributeIndex,  
     &      itemName, nameLength)
         call replaceNullsWithSpaces(itemName)
         if (trim(itemName).eq.'depth') then
            cycle
         endif
         do i=1,nAttr
            if (trim(itemName).eq.trim(na(i)%attrName)) then
               write(6,'(a)') 'loading nodal attribute data for '//trim(itemName)
               if (na(i)%numVals.eq.1) then
                  allocate(na(i)%xdmfArray(numOfNodes))
                  attributeType = XDMF_ATTRIBUTE_TYPE_SCALAR
                  numValues = numOfNodes * na(i) % numVals
                  call xdmfRetrieveAttributeValues(xdmfFortranObj, attributeIndex, 
     &                na(i)%xdmfArray, XDMF_ARRAY_TYPE_FLOAT64, numValues, 
     &                   startIndex, arrayStride, valueStride)
               else
                  attributeType = XDMF_ATTRIBUTE_TYPE_MATRIX
                  numValues = numOfNodes * na(i) % numVals
                  allocate(na(i)%xdmfMatrix(na(i)%numVals,numOfNodes))
                  call xdmfRetrieveAttributeValues(xdmfFortranObj, attributeIndex, 
     &                na(i)%xdmfMatrix, XDMF_ARRAY_TYPE_FLOAT64, numValues, 
     &                   startIndex, arrayStride, valueStride)
               endif
               exit
            endif
         end do
      end do
      !
      ! now place the data into the same data structures that are used
      ! by the ascii nodal attribute reading subroutine
      do i=1,nAttr
         SELECT CASE (trim(na(i)%AttrName))
         CASE("primitive_weighting_in_continuity_equation")
            IF (LoadTau0) THEN
               tau0var = na(i)%xdmfArray
            ENDIF
         CASE("surface_submergence_state")
            IF (LoadStartDry) THEN
               startdry  = na(i)%xdmfArray
            ENDIF
         CASE("quadratic_friction_coefficient_at_sea_floor")
            IF (LoadQuadraticFric) THEN
               fric = na(i)%xdmfArray
            ENDIF
         CASE("surface_directional_effective_roughness_length")
            IF (LoadDirEffRLen) THEN
               z0land  = na(i)%xdmfMatrix
            ENDIF
         CASE("surface_canopy_coefficient")
            IF (LoadCanopyCoef) THEN
               vcanopy = na(i)%xdmfArray
            ENDIF
         CASE("mannings_n_at_sea_floor")
            IF (LoadManningsN) THEN
               manningsn  = na(i)%xdmfArray
            ENDIF
         CASE("bottom_roughness_length")
            IF (LoadZ0b_var) THEN
               z0b_var = na(i)%xdmfArray
            ENDIF
         CASE("chezy_friction_coefficient_at_sea_floor")
            IF (LoadChezy) THEN
               chezy = na(i)%xdmfArray
            ENDIF
         CASE("sea_surface_height_above_geoid")
            IF (LoadGeoidOffset) THEN
               geoidOffset = na(i)%xdmfArray
            ENDIF
         CASE
     &   ("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")
            IF (LoadEVM) THEN
               evm = na(i)%xdmfArray
            ENDIF
         CASE
     &   ("average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth")
            IF (LoadEVC) THEN
               evc  = na(i)%xdmfArray
            ENDIF
         CASE
     &   ("min_and_max_primitive_weighting_in_continuity_equation")
            IF (LoadTau0MinMax) THEN
               tau0minmax  = na(i)%xdmfMatrix
            ENDIF
Casey 100210: Allow SWAN to handle wave refraction as a nodal attribute.
         CASE("wave_refraction_in_swan")
            IF (LoadSwanWaveRefrac) THEN
               swanWaveRefrac = na(i)%xdmfArray
            ENDIF
Corbitt 120321: Allow Advection to be handled locally instead of globally.
         CASE("advection_state")
            IF (LoadAdvectionState) THEN
               advectionState = na(i)%xdmfArray
            ENDIF
         CASE("elemental_slope_limiter")
            IF( LoadEleSlopeLim ) THEN
                elemental_slope_limiter_grad_max = na(i)%xdmfArray
            ENDIF
         CASE("initial_river_elevation")
            IF( LoadRiver_et_WSE ) THEN
               River_et_WSE = na(i)%xdmfArray
            ENDIF
C 02.24.2018 WJP: Add spatially varying linear internal tide friction coefficient
         CASE("internal_tide_friction")
            IF( LoadIT_Fric ) THEN
               IT_Fric = na(i)%xdmfMatrix
            ENDIF
         CASE("overland_reduction_factor")
            IF(LoadOverlandReductionFactor)THEN
                OverlandReductionFactor = na(i)%xdmfArray
            ENDIF
         CASE("subgrid_barrier")
             IF(LoadSubgridBarrier)THEN
                SubgridBarrier = na(i)%xdmfMatrix
             ENDIF
! DW---- for sponge layers ---
         CASE ("sponge_generator_layer") 
            IF (LoadAbsLayerSigma) THEN
               !WJP: need to add
            ENDIF
! DW-------------------------
         CASE DEFAULT
            scratchMessage = "Unrecognized nodal attribute : " //
     &         trim(na(i) % attrName)
            call allMessage(WARNING, scratchMessage)
         END SELECT
      END DO
#endif
      !-----------------------------------------------------------------------
      end subroutine readNodalAttrXDMF
      !-----------------------------------------------------------------------

C     ----------------------------------------------------------------
C                       S U B R O U T I N E 
C    C H E C K   F O R   M I S S I N G   N O D A L   A T T R I B U T E S
C     ----------------------------------------------------------------
C     jgf51.21.24: If a nodal attribute was specified in the control
C     (fort.15) file, but is not present in the nodal attributes (fort.13)
C     file, then we write out an error message and quit.
C     ----------------------------------------------------------------
      subroutine checkForMissingNodalAttributes()
      implicit none
C
C     Determine if there are any attributes required by the fort.15 file
C     that are not in the nodal attributes file.
      IF(((LoadTau0).and.(.not.FoundTau0)).or.
     &   ((LoadStartDry).and.(.not.FoundStartDry)).or.
     &   ((LoadQuadraticFric).and.
     &    (.not.FoundQuadraticFric)).or.
     &   ((LoadDirEffRLen).and.
     &    (.not.FoundDirEffRLen)).or.
     &   ((LoadCanopyCoef).and.
     &    (.not.FoundCanopyCoef)).or.
     &   ((LoadBridgePilings).and.
     &    (.not.FoundBridgePilings)).or.
     &   ((LoadManningsN).and.
     &    (.not.FoundManningsN)).or.
     &   ((LoadZ0b_var).and.
     &    (.not.FoundZ0b_var)).or.
     &   ((LoadGeoidOffset).and.
     &    (.not.FoundGeoidOffset)).or.
     &   ((LoadChezy).and.(.not.FoundChezy)).or.
     &   ((LoadEVM).and.(.not.FoundEVM)).or.
     &   ((LoadEVC).and.(.not.FoundEVC)).or.
Casey 100210: Allow SWAN to handle wave refraction as a nodal attribute.
     &   ((LoadSwanWaveRefrac).and.(.not.FoundSwanWaveRefrac)).or.
Corbitt 120321: Allow advection to be turned on locally instead of globally
     &   ((LoadAdvectionState).and.(.not.FoundAdvectionState)).or.
     &   ((LoadCondensedNodes).and.(.not.FoundCondensedNodes)).or.
     &   ((LoadEleSlopeLim).and.(.not.FoundEleSlopeLim)).or.
     &   ((LoadRiver_et_WSE).and.(.not.FoundRiver_et_WSE)).or.     
     &   ((LoadIT_Fric).and.(.not.FoundIT_Fric)).or. 
     &   ((LoadOverlandReductionFactor).and.
     &    (.not.FoundOverlandReductionFactor)).or.
!-- DW for sponge
     &   ((LoadAbsLayerSigma).and.(.not.FoundAbsLayerSigma)).or. 
     &   ((LoadTau0MinMax).and.(.not.FoundTau0MinMax))) THEN
         WRITE(scratchMessage,1111)
 1111   FORMAT('Nodal Attributes file (unit 13) does '
     &        'not contain all the attributes listed in the '
     &        'model parameter file (unit 15).')
         call allMessage(ERROR, scratchMessage)
         call na_terminate()
      ENDIF
C     ----------------------------------------------------------------
      end subroutine checkForMissingNodalAttributes
C     ----------------------------------------------------------------

C     ----------------------------------------------------------------
C                     S U B R O U T I N E  
C          A L L O C A T E   N O D A L   A T T R I B U T E S
C     ----------------------------------------------------------------
C     Allocate memory to hold our data.
      subroutine allocateNodalAttributes()
      use mesh,only:ne
      USE GLOBAL, ONLY: NTIF
      implicit none
      ALLOCATE(TAU0VAR(NumOfNodes),TAU0BASE(NumOfNodes)) ! jjw&sb46.39sb01
      ALLOCATE(STARTDRY(NumOfNodes))
      ALLOCATE(FRIC(NumOfNodes))
      ALLOCATE(z0land(NumOfNodes,DirEffRLenNoOfVals))
      ALLOCATE(vcanopy(NumOfNodes))
      ALLOCATE(BridgePilings(NumOfNodes,BridgePilingsNoOfVals))
      ALLOCATE(GeoidOffset(NumOfNodes))
      ALLOCATE(Chezy(NumOfNodes))
      ALLOCATE(ManningsN(NumOfNodes))
      ALLOCATE(Z0b_var(NumOfNodes))
      ALLOCATE(EVM(NumOfNodes))
      ALLOCATE(EVC(NumOfNodes))
      ALLOCATE(Tau0MinMax(NumOfNodes,Tau0MinMaxNoOfVals))
Casey 100210: Allow SWAN to handle wave refraction as a nodal attribute.
      ALLOCATE(SwanWaveRefrac(NumOfNodes))
Corbitt 120321: Allow advection to be turned on locally instead of globally
      ALLOCATE(AdvectionState(NumOfNodes))
      ALLOCATE(CondensedNodes(NumOfNodes,CondensedNodesNoOfVals))
      ALLOCATE(River_et_WSE(NumOfNodes))  ! tcm 20140502 v51.27
      ALLOCATE(OverlandReductionFactor(NumOfNodes))
      ALLOCATE(elemental_slope_limiter_grad_max(NumOfNodes))
      ALLOCATE(elemental_slope_limiter_active(NumOfNodes))
      ALLOCATE(elemental_slope_limiter_max_exceeded(NumOfNodes))
      elemental_slope_limiter_active(:) = .FALSE.  !...Alloate and initialize to keep
                                                   !   track of elemental_slope_limiter nodes
      elemental_slope_limiter_max_exceeded(:) = .FALSE.
      IF ( LoadIT_Fric ) THEN
         ALLOCATE(IT_Fric(NumOfNodes,IT_FricNoOfVals))
      END IF   
      OverlandReductionFactor(:) = 1D0/30D0
      IF ( LoadAbsLayerSigma ) THEN
         ALLOCATE( absorblayer_sigma(NumOfNodes,AbsLayerSigmaNoOfVals) )
      ENDIF
      IF( LoadSubgridBarrier )THEN 
          ALLOCATE(subgridBarrier(NumOfNodes,SubgridBarrierNoOfVals))
          ALLOCATE(subgridBarrierOvertopping(NumOfNodes))
      ENDIF
C     ----------------------------------------------------------------
      end subroutine allocateNodalAttributes
C     ----------------------------------------------------------------

C     ----------------------------------------------------------------
C     S U B R O U T I N E     R E A D  N O D A L  A T T R
C     ----------------------------------------------------------------
C
C     jgf46.00 Subroutine to read the nodal attributes file (unit 13).
C
C     ----------------------------------------------------------------
      SUBROUTINE ReadNodalAttr(NScreen, ScreenUnit, MyProc, NAbOut)
      USE MESH, ONLY : NP
      IMPLICIT NONE
      INTEGER, intent(in) :: NScreen ! nonzero for debug info to screen
      INTEGER, intent(in) :: ScreenUnit ! i/o for screen
      INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
      INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16
C
      LOGICAL NAFound  ! .true. if Nodal Attributes File (fort.13) exists
      INTEGER ErrorIO  ! zero if file opened successfully
      CHARACTER(len=80) AttrName ! string where the attribute name is stored
      CHARACTER(len=80) header   ! string where alphanumeric file id is stored
      INTEGER NumNodesNotDefault ! number of individual nodes to specify
      LOGICAL SkipDataSet ! .true. if a data set in unit 13 is not needed
      CHARACTER(len=80) Skipped ! data in unit 13 we do not need
      INTEGER L                 ! line counter
C
      ! temp array; used to load a real from the file,
      ! then convert to integer
      REAL(sz), ALLOCATABLE :: real_loader(:)
C
      NAFound = .False.
      SkipDataSet = .False.
C
C     Check to make sure that NWP is a valid number.
      IF (NWP.LT.0) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NWP =',NWP
            WRITE(ScreenUnit,9728)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NWP =',NWP
         WRITE(16,9728)
         WRITE(16,9973)
 9728    FORMAT(/,1X,'Your selection of NWP (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         CALL EXIT(1)          ! We're toast.
      ENDIF
C
C     Check to see if there are nodal attributes to be read in. If not,
C     simply return.
      IF (NWP.EQ.0) THEN
         WRITE(16,231) NWP
 231     FORMAT(/,5X,'NWP = ',I2,
     &        /,9X,'A Nodal Attributes File (unit 13)',
     &        /,9X,'will not be used.')
         RETURN
      ENDIF
C
C     Read the unit 15 control file to determine what data must be
C     loaded from nodal attributes file.
      WRITE(16,235) NWP
 235  FORMAT(/,9X,'Need to load ',I2,' nodal attribute(s):')
      DO k=1,NWP
         READ(15,*) AttrName
         WRITE(16,'(14X,A80)') AttrName
         SELECT CASE (AttrName)
         CASE("primitive_weighting_in_continuity_equation")
            LoadTau0 = .True.
         CASE("surface_submergence_state")
            LoadStartDry = .True.
         CASE("quadratic_friction_coefficient_at_sea_floor")
            LoadQuadraticFric = .True.
         CASE("surface_directional_effective_roughness_length")
            LoadDirEffRLen = .True.
         CASE("surface_canopy_coefficient")
            LoadCanopyCoef = .True.
         CASE("bridge_pilings_friction_parameters")
            LoadBridgePilings = .True.
         CASE("mannings_n_at_sea_floor")
            LoadManningsN = .True.
         CASE("chezy_friction_coefficient_at_sea_floor")
            LoadChezy = .True.
         CASE("bottom_roughness_length")
            LoadZ0b_var = .True.
         CASE("sea_surface_height_above_geoid")
            LoadGeoidOffset = .True.
         CASE
     &   ("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")
            LoadEVM = .True.
         CASE
     &   ("average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth")
            LoadEVC = .True.
         CASE
     &   ("min_and_max_primitive_weighting_in_continuity_equation")
            LoadTau0MinMax = .True.
Casey 100210: Allow SWAN to handle wave refraction as a nodal attribute.
         CASE("wave_refraction_in_swan")
            LoadSwanWaveRefrac = .TRUE.
Corbitt 120321: Allow advection to be turned on locally instead of globally
         CASE("advection_state")
            LoadAdvectionState = .True.
         CASE("condensed_nodes")
            LoadCondensedNodes = .True.
         CASE("elemental_slope_limiter")
            LoadEleSlopeLim = .TRUE.
         CASE("initial_river_elevation")  ! tcm 20140502 v51.27
            LoadRiver_et_WSE = .TRUE.
C 09232016 WP: Add spatially varying linear internal tide friction coefficient
         CASE("internal_tide_friction")
            LoadIT_Fric = .True.
         CASE("subgrid_barrier")
            LoadSubgridBarrier = .True.
         CASE("overland_reduction_factor")
            LoadOverlandReductionFactor = .True.
! DW---- for sponge layers ---
         CASE ("sponge_generator_layer") 
            LoadAbsLayerSigma = .TRUE.
         CASE DEFAULT
            WRITE(16,1000)          ! unit 15 Model Parameter file
            WRITE(16,1021) AttrName ! contains invalid name
            IF (NScreen.ne.0.and.MyProc.eq.0) THEN
               WRITE(ScreenUnit,1000)
               WRITE(ScreenUnit,1021) AttrName
            ENDIF
         END SELECT
      ENDDO
C
      WRITE(16,232) NWP
 232  FORMAT(/,5X,'NWP = ',I2,
     &     /,9X,'Must read Nodal Attributes File (unit 13).')
      !
      !  R E A D   N O D A L   A T T R I B U T E S   X D M F  
      !
      if (naType.eq.XDMF) then
         call readNodalAttrXDMF()
         return
      endif
      !
      !  R E A D   N O D A L   A T T R I B U T E S   A S C I I
      !
C     Determine if the Nodal Attributes File exists.
      INQUIRE(FILE=TRIM(INPUTDIR)//'/'//'fort.13',EXIST=NAFound)
C
      IF (.not.NAFound) THEN
         WRITE(16,1001)         ! Nodal Attributes file
         WRITE(16,1011)         ! was not found.
         WRITE(16,9973)         ! execution terminated
         IF (NScreen.ne.0.and.MyProc.eq.0) THEN
            WRITE(ScreenUnit,1001)
            WRITE(ScreenUnit,1011)
            WRITE(ScreenUnit,9973)      ! execution terminated
         ENDIF
         CALL EXIT(1)
      ENDIF
C
C     Now open the nodal attributes (unit 13) file.
      WRITE(16,240)
 240  FORMAT(/,9X,'Nodal Attributes File (unit 13) was found.',
     &     ' Opening file.')
      OPEN(UNIT=13, FILE=TRIM(INPUTDIR)//'/'//'fort.13',
     &     IOSTAT=ErrorIO)
      IF ( ErrorIO .GT. 0 ) THEN
         WRITE(16,1001)         ! Nodal attribute file
         WRITE(16,1005)         ! exists but can't be opened
         WRITE(16,9973)         ! execution terminated
         IF (NScreen.ne.0.and.MyProc.eq.0) THEN
            WRITE(ScreenUnit,1001)
            WRITE(ScreenUnit,1005)
            WRITE(ScreenUnit,9973)
         ENDIF
         CALL EXIT(1)                 ! We're toast.
      ENDIF
C
C     Read each attribute name, units, number of values, and default value
      READ(13,'(A80)') header
      WRITE(16,250)
 250  FORMAT(/,9X,'User comment line from unit 13:')
      WRITE(16,'(14X,A80,/)') header
      READ(13,*) NumOfNodes     ! number of nodes according to unit 13     
C      
C     ERROR CHECK: If a nodal attributes file is being used, check to
C     see that the number of nodes in the nodal attribute file is the
C     same as the number of nodes in the grid file.
C     jgf52.18: Make this check immediately so execution can 
C     be stopped immediately if this is not the right nodal attributes
C     file for this mesh.
      if ( (nwp.ne.0).and.(numOfNodes.ne.np) ) then
         write(scratchMessage,9901) np, numOfNodes
         call allMessage(ERROR,scratchMessage)
         call na_terminate()
 9901   format('The number of nodes in the mesh file (unit 14) is ',i0,
     &  'while the number of nodes listed on the 2nd line of the '//
     &  'nodal attributes file (unit 13) is ',i0,'.'//
     &  'These numbers must match. '//
     &  'Is this the right nodal attributes file for this mesh file?')
      endif
C      
      READ(13,*) NAttr          ! number of attributes in the unit 13 file
      DO k=1, NAttr
         READ(13,*) AttrName
         WRITE(16,'(9X,A80)') AttrName
         WRITE(16,260)
 260     FORMAT(14X,'was found!',/)
         SELECT CASE (AttrName)
         CASE("primitive_weighting_in_continuity_equation")
            FoundTau0 = .True.
            READ(13,'(A80)') Tau0Units
            READ(13,*) Tau0NoOfVals
            READ(13,*) Tau0DefVal
         CASE("surface_submergence_state")
            FoundStartDry = .True.
            READ(13,'(A80)') StartDryUnits
            READ(13,*) StartDryNoOfVals
            READ(13,*) StartDryDefVal
         CASE("quadratic_friction_coefficient_at_sea_floor")
            FoundQuadraticFric = .True.
            READ(13,'(A80)') QuadraticFricUnits
            READ(13,*) QuadraticFricNoOfVals
            READ(13,*) QuadraticFricDefVal
         CASE("surface_directional_effective_roughness_length")
            FoundDirEffRLen = .True.
            READ(13,'(A80)') DirEffRLenUnits
            READ(13,*) DirEffRLenNoOfVals
            READ(13,*)
     &      (DirEffRLenDefVal(j),j=1,DirEffRLenNoOfVals)
         CASE("surface_canopy_coefficient")
            FoundCanopyCoef = .True.
            READ(13,'(A80)') CanopyCoefUnits
            READ(13,*) CanopyCoefNoOfVals
            READ(13,*) CanopyCoefDefVal
         CASE("bridge_pilings_friction_parameters")
            FoundBridgePilings = .True.
            READ(13,'(A80)') BridgePilingsUnits
            READ(13,*) BridgePilingsNoOfVals
            READ(13,*)
     &       (BridgePilingsDefVal(j),j=1,BridgePilingsNoOfVals)
         CASE("mannings_n_at_sea_floor")
            FoundManningsN = .True.
            READ(13,'(A80)') ManningsNUnits
            READ(13,*) ManningsNNoOfVals
            READ(13,*) ManningsNDefVal
         CASE("bottom_roughness_length")
            FoundZ0b_var = .True.
            READ(13,'(A80)') Z0b_varUnits
            READ(13,*) Z0b_varNoOfVals
            READ(13,*) Z0b_varDefVal
         CASE("chezy_friction_coefficient_at_sea_floor")
            FoundChezy = .True.
            READ(13,'(A80)') ChezyUnits
            READ(13,*) ChezyNoOfVals
            READ(13,*) ChezyDefVal
         CASE("sea_surface_height_above_geoid")
            FoundGeoidOffset = .True.
            READ(13,'(A80)') GeoidOffsetUnits
            READ(13,*) GeoidOffsetNoOfVals
            READ(13,*) GeoidOffsetDefVal
         CASE
     &    ("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")
            FoundEVM = .True.
            READ(13,'(A80)') EVMUnits
            READ(13,*) EVMNoOfVals
            READ(13,*) EVMDefVal
         CASE
     &   ("average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth")
            READ(13,'(A80)') EVCUnits
            READ(13,*) EVCNoOfVals
            READ(13,*) EVCDefVal
         CASE
     &   ("min_and_max_primitive_weighting_in_continuity_equation")
            FoundTau0MinMax = .True.
            READ(13,'(A80)') Tau0MinMaxUnits
            READ(13,*) Tau0MinMaxNoOfVals
            READ(13,*) (Tau0MinMaxDefVal(j),j=1,Tau0MinMaxNoOfVals)
Casey 100210: Allow SWAN to handle wave refraction as a nodal attribute.
         CASE("wave_refraction_in_swan")
            FoundSwanWaveRefrac = .TRUE.
            READ(13,'(A80)') SwanWaveRefracUnits
            READ(13,*) SwanWaveRefracNoOfVals
            READ(13,*) SwanWaveRefracDefVal
Corbitt 120321: Allow advection to be turned on locally instead of globally
         CASE("advection_state")
            FoundAdvectionState = .TRUE.
            READ(13,'(A80)') AdvectionStateUnits
            READ(13,*) AdvectionStateNoOfVals
            READ(13,*) AdvectionStateDefVal
         CASE("condensed_nodes")
            FoundCondensedNodes = .TRUE.
            READ(13,'(A80)') CondensedNodesUnits
            READ(13,*) CondensedNodesNoOfVals
            READ(13,*) (CondensedNodesDefVal(j),
     &                  j=1,CondensedNodesNoOfVals)
         CASE("elemental_slope_limiter")
            FoundEleSlopeLim = .TRUE.
            READ(13,'(A80)') EleSlopeLimUnits
            READ(13,*) EleSlopeLimNoOfVals
            READ(13,*) EleSlopeLimDefVal
         CASE("initial_river_elevation")  ! tcm 20140502 v51.27
            FoundRiver_et_WSE = .TRUE.
            READ(13,'(A80)') River_et_WSEUnits
            READ(13,*) River_et_WSENoOfVals
            READ(13,*) River_et_WSEDefVal            
C 09232016 WP: Add spatially varying linear internal tide friction
         CASE("internal_tide_friction")
            FoundIT_Fric = .True.              
            READ(13,'(A80)') IT_FricUnits
            READ(13,*) IT_FricNoOfVals
            READ(13,*) (IT_FricDefVal(j),j=1,IT_FricNoOfVals)
         CASE("overland_reduction_factor")
            FoundOverlandReductionFactor = .TRUE.
            READ(13,'(A80)') OverlandReductionFactorUnits
            READ(13,*) OverlandReductionFactorNoOfVals
            READ(13,*) OverlandReductionFactorDefaultValue
         CASE("subgrid_barrier")
            FoundSubgridBarrier = .TRUE.
            READ(13,'(A80)') SubgridBarrierUnits
            READ(13,*) SubgridBarrierNoOfVals
            READ(13,*) (SubgridBarrierDefaultValue(j),j=1,SubgridBarrierNoOfVals)
C---- DW ----- for sponge layer ----------
         CASE ("sponge_generator_layer")  
            FoundAbsLayerSigma = .TRUE.
            READ(13,'(A80)') AbsLayerSigmaUnits 
            READ(13,*) AbsLayerSigmaNoOfVals
            READ(13,*) 
     &           (AbsLayerSigmaDefVal(j),j=1,AbsLayerSigmaNoOfVals)
         CASE DEFAULT
            WRITE(16,1001)          ! Nodal Attributes file
            WRITE(16,1021) AttrName ! contains invalid name
            IF (NScreen.ne.0.and.MyProc.eq.0) THEN
               WRITE(ScreenUnit,1001)
               WRITE(ScreenUnit,1021) AttrName
            ENDIF
            READ(13,'(A80)') Skipped  ! skip the Units for the invalid name
            READ(13,'(A80)') Skipped  ! skip the NoOfVals for invalid name
            READ(13,'(A80)') Skipped  !jgf51.40: skip the default value
         END SELECT
      END DO
C
C     Determine if there are any attributes required by the fort.15 file
C     that are not in the nodal attributes file.
      call checkForMissingNodalAttributes()
C
C     Allocate memory for nodal attributes
      call allocateNodalAttributes()
C
C     Now read each of the attributes required by the model parameter
C     (unit 15) file and skip past the others.
      WRITE(16,270) NWP
 270  FORMAT(/,9X,'Now reading ',I2,' nodal attribute(s).')
      DO k=1, NAttr
         WRITE(16,280) k
 280     FORMAT(/,9X,'Attribute ',I2,':')
         READ(13,*) AttrName
         READ(13,*) NumNodesNotDefault
         WRITE(16,'(14X,A80)') AttrName
         SELECT CASE (AttrName)
         CASE("primitive_weighting_in_continuity_equation")
            IF (LoadTau0) THEN
               CALL LoadAttrVec(TAU0VAR, Tau0DefVal,
     &              NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("surface_submergence_state")
            IF (LoadStartDry) THEN
               CALL LoadAttrVec(STARTDRY, StartDryDefVal,
     &               NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("quadratic_friction_coefficient_at_sea_floor")
            IF (LoadQuadraticFric) THEN
               CALL LoadAttrVec(FRIC, QuadraticFricDefVal,
     &              NumNodesNotDefault, NScreen, MyProc, NAbOut)
               FRIC_ORIG = FRIC
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("surface_directional_effective_roughness_length")
            IF (LoadDirEffRLen) THEN
               CALL LoadAttrMat(z0land, DirEffRLenNoOfVals,
     &              DirEffRLenDefVal, NumNodesNotDefault,
     &              NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("surface_canopy_coefficient")
            IF (LoadCanopyCoef) THEN
               CALL LoadAttrVec(vcanopy, CanopyCoefDefVal,
     &               NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("bridge_pilings_friction_parameters")
            IF (LoadBridgePilings) THEN
               CALL LoadAttrMat(BridgePilings, BridgePilingsNoOfVals,
     &              BridgePilingsDefVal,  NumNodesNotDefault,
     &              NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("mannings_n_at_sea_floor")
            IF (LoadManningsN) THEN
               CALL LoadAttrVec(ManningsN, ManningsNDefVal,
     &              NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("bottom_roughness_length")
            IF (LoadZ0b_var) THEN
               CALL LoadAttrVec(Z0b_var, Z0b_varDefVal,
     &              NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("chezy_friction_coefficient_at_sea_floor")
            IF (LoadChezy) THEN
               CALL LoadAttrVec(Chezy, ChezyDefVal,
     &               NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("sea_surface_height_above_geoid")
            IF (LoadGeoidOffset) THEN
               CALL LoadAttrVec(GeoidOffset, GeoidOffsetDefVal,
     &               NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE
     &   ("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")
            IF (LoadEVM) THEN
               CALL LoadAttrVec(EVM, EVMDefVal,
     &               NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE
     &   ("average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth")
            IF (LoadEVC) THEN
               CALL LoadAttrVec(EVC, EVCDefVal,
     &               NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE
     &   ("min_and_max_primitive_weighting_in_continuity_equation")
            IF (LoadTau0MinMax) THEN
               CALL LoadAttrMat(Tau0MinMax, Tau0MinMaxNoOfVals,
     &              Tau0MinMaxDefVal, NumNodesNotDefault,
     &              NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
Casey 100210: Allow SWAN to handle wave refraction as a nodal attribute.
         CASE("wave_refraction_in_swan")
            IF (LoadSwanWaveRefrac) THEN
               CALL LoadAttrVec(SwanWaveRefrac, SwanWaveRefracDefVal,
     &              NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .TRUE.
            ENDIF
Corbitt 120321: Allow Advection to be handled locally instead of globally.
         CASE("advection_state")
            IF (LoadAdvectionState) THEN
               CALL LoadAttrVec(AdvectionState,AdvectionStateDefVal,
     &              NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .TRUE.
            ENDIF
         CASE("condensed_nodes")
            IF (LoadCondensedNodes) THEN
               CALL LoadAttrMatNodeNo(CondensedNodes,
     &              CondensedNodesNoOfVals,
     &              CondensedNodesDefVal,
     &              NumNodesNotDefault, NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .TRUE.
            ENDIF
         CASE("elemental_slope_limiter")
            IF( LoadEleSlopeLim ) THEN
                CALL LoadAttrVec(elemental_slope_limiter_grad_max,
     &              EleSlopeLimDefVal, NumNodesNotDefault, NScreen,
     &              MyProc, NAbOut)
            ELSE
                SkipDataSet = .TRUE.
            ENDIF
         CASE("initial_river_elevation")   ! tcm 20140502 v51.27
            IF( LoadRiver_et_WSE ) THEN
                CALL LoadAttrVec(River_et_WSE,
     &              River_et_WSEDefVal, NumNodesNotDefault, NScreen,
     &              MyProc, NAbOut)
            ELSE
                SkipDataSet = .TRUE.
            ENDIF
C 09232016 WP : Spatially varying linear internal tide friction
         CASE("internal_tide_friction")
            IF (LoadIT_Fric) THEN
               CALL LoadAttrMat(IT_Fric, IT_FricNoOfVals,
     &              IT_FricDefVal, NumNodesNotDefault,
     &              NScreen, MyProc, NAbOut)
            ELSE
               SkipDataSet = .True.
            ENDIF
         CASE("overland_reduction_factor")
            IF(LoadOverlandReductionFactor)THEN
                CALL LoadAttrVec(OverlandReductionFactor,
     &              OverlandReductionFactorDefaultValue,
     &              NumNodesNotDefault,nScreen,MyProc,
     &              NAbOut)
            ELSE
              SkipDataSet = .True.
            ENDIF
         CASE("subgrid_barrier")
             IF(LoadSubgridBarrier)THEN
                 CALL LoadAttrMat(SubgridBarrier,SubgridBarrierNoOfVals,
     &              subgridBarrierDefaultValue,NumNodesNotDefault,
     &              nScreen,MyProc,NAbOut)
             ELSE
               SkipDataSet = .True.
             ENDIF
C          WP
!------- DW ---- Absorbing layer -----------
         CASE ("sponge_generator_layer") 
            IF (  LoadAbsLayerSigma ) THEN
               !
               NumNodesAbsLayer(:) = NumNodesNotDefault ;
               IF ( .NOT. ASSOCIATED(AbsLayerNodesID) ) THEN
                  ALLOCATE( AbsLayerNodesID(NumNodesAbsLayer(1)) ) ;
               END IF

               CALL LoadAttrMat(absorblayer_sigma, 
     &              AbsLayerSigmaNoOfVals, AbsLayerSigmaDefVal,
     &              NumNodesNotDefault, NScreen, MyProc, NAbOut, 
     &              NodesNotDef = AbsLayerNodesID ) ;
               !
               absorblayer_sigma_eta => absorblayer_sigma ;
               absorblayer_sigma_mnx => absorblayer_sigma ;
               absorblayer_sigma_mny => absorblayer_sigma ;
               
               ! WJP: for only eta sponge shit 
               !ALLOCATE( absorblayer_zero(NumOfNodes,AbsLayerSigmaNoOfVals))
               !absorblayer_zero = absorblayer_sigma; 
               !absorblayer_zero(:,1) = 0d0
               !absorblayer_sigma_mnx => absorblayer_zero;
               !absorblayer_sigma_mny => absorblayer_zero;
               !if (NumNodesNotDefault > 0) then
               !   write(16,*) 'Set absorblayer sigma flux to zero:'
               !   write(16,*) absorblayer_sigma_eta(1,1),
     &         !               absorblayer_sigma_mnx(1,1), 
     &         !               absorblayer_sigma_mny(1,1) 
               !endif
               AbsLayer_Eta_NodesID => AbsLayerNodesID ;
               AbsLayer_MnX_NodesID => AbsLayerNodesID ;
               AbsLayer_MnY_NodesID => AbsLayerNodesID ;
            ELSE
                SkipDataSet = .TRUE. 
            END IF
!--------DW -----------------------------------
         CASE DEFAULT
            SkipDataSet = .True.
            WRITE(16,1001)      ! Nodal Attributes file
            WRITE(16,1021) AttrName ! contains invalid name
            IF (NScreen.ne.0.and.MyProc.eq.0) THEN
               WRITE(ScreenUnit,1001)
               WRITE(ScreenUnit,1021) AttrName
            ENDIF
         END SELECT
         IF (SkipDataSet) THEN
            DO L=1, NumNodesNotDefault
               READ(13,*) Skipped
            END DO
            WRITE(16,'(9X,A8)') 'Skipped.'
            SkipDataSet = .False.
         ELSE
            WRITE(16,'(/,9X,A18,A80)') 'Finished loading ', AttrName
         ENDIF
      END DO
C
 1000 FORMAT('ERROR: The Model Parameter File (unit 15)')
 1001 FORMAT('ERROR: The Nodal Attributes File (unit 13)')
 1002 FORMAT('ERROR: The legacy StartDry File (unit 12)')
 1003 FORMAT('ERROR: Spatially Varying Fric. Coeff. File (unit 21)')
C
 1005 FORMAT('exists but cannot be opened.')
 1011 FORMAT('was not found.')
 1021 FORMAT('contains invalid name: ',A80)
 9972 FORMAT(////,1X,'!!!!!!!!!! INPUT ERROR !!!!!!!!!',/)
 9973 FORMAT(/,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE ReadNodalAttr
C     ----------------------------------------------------------------


C     ----------------------------------------------------------------
C         S U B R O U T I N E     L O A D  A T T R  V E C
C     ----------------------------------------------------------------
C
C     jgf46.00 Subroutine to set a single set of nodal attributes to
C     their user-specified default values, then read the nondefault
C     values from the Nodal Attributes File (unit 13). This subroutine
C     is used for nodal attributes with only one value per node, hence
C     the suffix "vec" in the name.
C
C     ----------------------------------------------------------------
      SUBROUTINE LoadAttrVec(AttributeData, Default, NumNodesNotDef,
     &     NScreen, MyProc, NAbOut)
      USE MESH,ONLY:node_dict,find
      IMPLICIT NONE
      REAL(SZ), intent(out), dimension(NumOfNodes) :: AttributeData
      REAL(SZ), intent(in):: Default ! default value for all nodes
      INTEGER, intent(in) :: NumNodesNotDef ! number of nodes specified in file
      INTEGER, intent(in) :: NScreen ! 1 for debug info to screen (unit 6)
      INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
      INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16
      REAL(SZ)            :: TempData 
      INTEGER             :: NodeNum ! node label listed in the file
      INTEGER             :: DictNode ! node position in array
C
C     Set all values to user-specified default values.
      IF (NABOUT.EQ.0) WRITE(16,1001) Default
      DO i=1, NumOfNodes
         AttributeData(i) = Default
      END DO
C
      IF (NABOUT.EQ.0) WRITE(16,1005)
      DO i=1, NumNodesNotDef
         READ(13,*) NodeNum, TempData
         DictNode = find(node_dict,NodeNum)
         AttributeData(DictNode) = TempData
         IF (NABOUT.EQ.0)
     &       WRITE(16,1010) NodeNum, TempData
      END DO
C
 1001 FORMAT(/,10X,'Set all nodes to the default value of ',E16.8,/)
 1005 FORMAT(/,10X,'Now setting the following nodes to these values:',
     &     /,10X,'NODE',5X,'DATA',5X/)
 1010 FORMAT(7X,I8,6X,E16.8)
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE LoadAttrVec
C     ----------------------------------------------------------------


C     ----------------------------------------------------------------
C         S U B R O U T I N E     L O A D  A T T R  M A T
C     ----------------------------------------------------------------
C
C     jgf46.00 Subroutine to load a single set of nodal attributes from
C     the Nodal Attributes File (unit 13) if there is more than one
C     value per node.
C
C     ----------------------------------------------------------------
      SUBROUTINE LoadAttrMat(AttributeData, NumCol, Default,
     &     NumNodesNotDef, NScreen, MyProc, NAbOut, NodesNotDef )
      USE MESH,ONLY:node_dict,find
      IMPLICIT NONE
      INTEGER, intent(in) :: NumCol  ! number of columns in the matrix
      REAL(SZ), intent(out),
     &     dimension(NumOfNodes,NumCol) :: AttributeData
      REAL(SZ), intent(in), dimension(NumCol) :: Default ! default values
      INTEGER, intent(in) :: NumNodesNotDef  ! number of nodes spec. in file
      INTEGER, intent(in) :: NScreen ! 1 for debug info to screen (unit 6)
      INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
      INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16
      REAL(SZ)            :: TempAttData(NumCol) 
      INTEGER             :: NodeNum ! node label listed in the file
      INTEGER             :: DictNode ! node position in array
C------ DW
      INTEGER, optional::  NodesNotDef(:)
C------
C
C     Set all nodes to user-specified default values.
      IF (NABOUT.EQ.0) WRITE(16,1001)
      DO i=1, NumOfNodes
         DO j=1, NumCol
            AttributeData(i,j)=Default(j)
         END DO
      END DO
C
      IF (NABOUT.EQ.0) WRITE(16,1005)
      IF (NABOUT.EQ.0) WRITE(16,*) NumCol,ubound(AttributeData),NumNodesNotDef
      DO i=1, NumNodesNotDef
         READ(13,*) NodeNum, (TempAttData(j),j=1,NumCol)
         DictNode = find(node_dict,NodeNum)
         IF ( present(NodesNotDef) ) THEN
            NodesNotDef(i) = DictNode 
         END IF
         AttributeData(DictNode,:) = TempAttData(:)
         IF (NABOUT.EQ.0) WRITE(16,1010) NodeNum,
     &        (AttributeData(DictNode,j),j=1,NumCol)
      END DO
C
 1001 FORMAT(/,10X,'Set all nodes to the default values of ',/,
     &     99E16.8,/)
 1005 FORMAT(/,10X,'Now setting the following nodes to these values:',
     &     /,10X,'NODE',5X,'DATA',5X/)
 1010 FORMAT(7X,I8,6X,12(1X,E16.8))
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE LoadAttrMat
C     ----------------------------------------------------------------



C     ----------------------------------------------------------------
C         S U B R O U T I N E     L O A D  A T T R  M A T  N O D E  N O
C     ----------------------------------------------------------------
C
C     Subroutine to load a single set of nodal attributes from
C     the Nodal Attributes File (unit 13) if there is more than one
C     value per node.
C
C     Added for integer node no attributes. Based on LoadAttrMat. SB
C
C     ----------------------------------------------------------------
      SUBROUTINE LoadAttrMatNodeNo(AttributeData, NumCol, Default,
     &     NumNodesNotDef, NScreen, MyProc, NAbOut, NodesNotDef )
      USE MESH,ONLY:node_dict,find
      IMPLICIT NONE
      INTEGER, intent(in) :: NumCol  ! number of columns in the matrix
      INTEGER, intent(out),
     &     dimension(NumOfNodes,NumCol) :: AttributeData
      INTEGER, intent(in), dimension(NumCol) :: Default ! default values
      INTEGER, intent(in) :: NumNodesNotDef  ! number of nodes spec. in file
      INTEGER, intent(in) :: NScreen ! 1 for debug info to screen (unit 6)
      INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
      INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16
      INTEGER             :: TempAttData(NumCol) 
      INTEGER             :: NodeNum ! node label listed in the file
      INTEGER             :: DictNode ! node position in array
C------ DW
      INTEGER, optional::  NodesNotDef(:)
C------
C
C     Set all nodes to user-specified default values.
      IF (NABOUT.EQ.0) WRITE(16,1001)
      DO i=1, NumOfNodes
         DO j=1, NumCol
            AttributeData(i,j)=Default(j)
         END DO
      END DO
C
      IF (NABOUT.EQ.0) WRITE(16,1005)
      IF (NABOUT.EQ.0) WRITE(16,*) NumCol,ubound(AttributeData),NumNodesNotDef
      DO i=1, NumNodesNotDef
         READ(13,*) NodeNum, (TempAttData(j),j=1,NumCol)
         DictNode = find(node_dict,NodeNum)
         DO j=1,NumCol
            TempAttData(j) = find(node_dict,TempAttData(j))
         END DO
         IF ( present(NodesNotDef) ) THEN
            NodesNotDef(i) = DictNode 
         END IF
         AttributeData(DictNode,:) = TempAttData(:)
         IF (NABOUT.EQ.0) WRITE(16,1010) NodeNum,
     &        (AttributeData(DictNode,j),j=1,NumCol)
      END DO
C
 1001 FORMAT(/,10X,'Set all nodes to the default values of ',/,
     &     99I11,/)
 1005 FORMAT(/,10X,'Now setting the following nodes to these values:',
     &     /,10X,'NODE',5X,'DATA',5X/)
 1010 FORMAT(7X,I8,6X,12(1X,I11))
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE LoadAttrMatNODENO
C     ----------------------------------------------------------------



C     ----------------------------------------------------------------
C         S U B R O U T I N E     I N I T  N O D A L  A T T R
C     ----------------------------------------------------------------
C
C     jgf46.00 Subroutine to initialize and error check the nodal
C     attributes read in from the Nodal Attributes File (unit 13).
C
C     ----------------------------------------------------------------
      SUBROUTINE InitNodalAttr(DP, NP, G, NScreen, ScreenUnit,
     &     MyProc, NAbOut, Z0B)
      USE GLOBAL, ONLY : C3D, C2DDI, CBaroclinic
      IMPLICIT NONE
      INTEGER, intent(in) :: NP ! number of nodes in the grid file
      REAL(SZ), intent(in), dimension(NP) :: DP ! array of bathymetric depths
      REAL(SZ), intent(in):: G  ! gravitational acceleration
      INTEGER, intent(in) :: NScreen ! nonzero for debug info to screen
      INTEGER, intent(in) :: ScreenUnit ! i/o for debug info to screen
      INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
      INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16
      REAL(SZ), intent(in) :: Z0B ! 3D bottom friction
      INTEGER Tau0Dig1  ! determines the tau0 scheme
      INTEGER Tau0Dig2  ! determines whether tau0 is being output
      LOGICAL :: invalidCanopyCoefficient ! .true. if any value is neither 0 nor 1

      IF (Tau0.lt.0) THEN
         Tau0Dig1 = INT(Tau0)   ! jgf47.30 truncate the fractional part
         !jgf47.34 round away from zero by subtracting 0.5d0
         Tau0Dig2 = INT( (Tau0 - REAL(Tau0Dig1))*10.0d0 - 0.5d0)
      ELSE
         Tau0Dig1 = 0
         Tau0Dig2 = 0
      ENDIF

C     ERROR CHECK: If a nodal attributes file is being used, check to
C     see that the number of nodes in the nodal attribute file is the
C     same as the number of nodes in the grid file.
      IF (NWP.NE.0.AND.NumOfNodes.NE.NP) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9900)
         WRITE(16,9900)
 9900    FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!',
     &        //,1X,'The number of nodes in the grid file (unit 14) and'
     &        /,1X,'the nodal attributes file (unit 13) must match.',
     &        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
         CALL EXIT(1)                 ! We're toast.
      ENDIF
C
C     ERROR CHECK: If Chezy, Manning's or Quadratic friction was loaded
C     from the nodal attributes file, NOLIBF must be 1.
      IF ((LoadChezy.or.LoadManningsN.or.LoadQuadraticFric).and.
     &     NoLiBF.ne.1) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9800)
         WRITE(16,9800)
 9800    FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!',
     &        //,1X,'Nonlinear bottom friction coefficients were loaded'
     &        /,1X,'from the nodal attributes file (unit 13), so ',
     &        /,1X,'NoLiBF must be set to 1. It is set to ',i2,' in',
     &        /,1X,'the model parameter (unit 15) file.',
     &        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
         CALL EXIT(1)                 ! We're toast.
      ENDIF
C
C     ERROR CHECK: If Tau0=-3.x or -6.x in fort.15, then tau0 MUST be loaded
C     from nodal attributes file.
      IF ( ((Tau0Dig1.eq.-3).or.(Tau0Dig1.eq.-6))
     &          .and.(.not.LoadTau0) ) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9700) Tau0
         WRITE(16,9700)
 9700    FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!',
     &        //,1X,'Spatially and temporally varying tau0 was '
     &        /,1X,'specified in the fort.15 file with Tau0=',E9.2,
     &        /,1X,'but the base value was not specified in the '
     &        /,1X,'nodal attributes file (unit 13). Please ',
     &        /,1X,'load the base value using',
     &        /,1X,'primitive_weighting_in_continuity_equation.',
     &        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
         CALL EXIT(1)                 ! We're toast.
      ENDIF
C
C     ERROR CHECK: The canopy coefficient value should be either 0 or 1;
C     log a warning message if any other values were read.
      if ( loadCanopyCoef.eqv..true. ) then
         invalidCanopyCoefficient = .false.
         do i=1,np
            if ( (vcanopy(i).ne.0.d0).and.(vcanopy(i).ne.1.d0 ) ) then
               invalidCanopyCoefficient = .true.
               exit
            endif
         end do
         if ( invalidCanopyCoefficient.eqv..true. ) then 
            call allMessage(WARNING,'Some surface canopy coefficient '//
     &        'nodal attribute data are neither 0 nor 1.')        
         endif
      endif
C
C     I N I T    S T A R T D R Y
 
      IF (NWP.eq.0) THEN
         ALLOCATE(STARTDRY(NP))
      ENDIF
      IF (LoadStartDry.eqv..false.) THEN
         DO I=1, NP
            STARTDRY(I) = 0.0D0
         ENDDO
      ENDIF
C
C     I N I T     T A U 0
      IF (NWP.eq.0) THEN
         ALLOCATE(TAU0VAR(NP),TAU0BASE(NP))
         ALLOCATE(Tau0MinMax(NP,Tau0MinMaxNoOfVals))
      ENDIF
C
C     jgf46.25 If input tau0 is positive, set all nodes to that value.
      IF (.not.LoadTau0) THEN
         IF (Tau0.ge.0) THEN
            DO I=1,NP
              Tau0Var(I)=Tau0
            END DO
            WRITE(16,7) Tau0
 7          FORMAT(/,5X,
     &           'A SPATIALLY CONSTANT WEIGHTING COEFFICIENT (Tau0)'
     &           ,/,5X,' WILL BE USED IN THE GENERALIZED WAVE',
     &           ' CONTINUITY EQUATION.',
     &           /,5X,'Tau0 = ',E15.8,2X,'1/sec',/)
         ELSE
C           If input tau0 is negative, set value using hardcoded scheme
C           based on depth
            DO I=1,NP
               Tau0Var(I)=Tau0NodalValue(Tau0,DP(I))
            ENDDO
            ! jgf47.30.TODO: This logging needs to be cleaned up.
            IF(Tau0.eq.-2) THEN
               WRITE(16,6) ! spatially vary tau0 according to hard coded scheme
               WRITE(16,62) ! description of scheme
 62            FORMAT(/,5X,'IF DEPTH > 200           Tau0 = 0.005',
     &              /,5X,'IF 200   > DEPTH > 1     Tau0 = 1/DEPTH  ',
     &              /,5X,'IF 1     > DEPTH         Tau0 = 1.0 ')
            ENDIF
            IF (.not.( (Tau0Dig1.eq.-3) .or. (Tau0Dig1.eq.-5) ) ) THEN
               WRITE(16,6) ! spatially vary tau0 according to hard coded scheme
               WRITE(16,61) ! description of scheme
 61            FORMAT(/,5X,' IF DEPTH GE 10           -> TAU0 = 0.005',
     &              /,5X,' IF DEPTH LT 10           -> TAU0 = 0.020',/)
            ENDIF
         ENDIF
      ENDIF
 6    FORMAT(/,5X,'A SPATIALLY VARIABLE WEIGHTING COEFFICIENT (Tau0)'
     &     ,/,5X,' WILL BE USED IN THE GENERALIZED WAVE',
     &     ' CONTINUITY EQUATION.',
     &     /,5x,'THIS VALUE WILL BE DETERMINED AS FOLLOWS:')
C
C     jgf46.27 If we have already loaded the tau0 values directly from
C     the nodal attributes file, check to see if the default value was
C     negative. If so, this indicates that nodal values of tau0 that
C     were not explicitly set in the nodal attributes file should be set
C     according to one of the hard-coded tau0 schemes.
      IF (LoadTau0.and.Tau0DefVal.lt.0) THEN
         DO I=1,NP
            IF (Tau0Var(I).lt.0) THEN
               Tau0Var(I)=Tau0NodalValue(Tau0DefVal,DP(I))
            ENDIF
         ENDDO
      ENDIF
C
C     jgf47.06 Activate time varying tau0 and output tau0 if these options
C     were selected.

C     jgf47.30 Use Joannes' scheme for steady Tau0 in deep water and
C     other coarsely gridded areas, and time varying Tau0 in high
C     resolution areas.
      IF ( (Tau0Dig1.eq.-3).or.(Tau0Dig1.eq.-6)
     &      .or. (Tau0Dig1.eq.-7) ) THEN
         HighResTimeVaryingTau0 = .True.
         DO I=1, NP
            Tau0Base(I) = Tau0Var(I)
         ENDDO
      ENDIF
C
C     jgf47.11 Also allow the min and max tau0 to be set from the fort.15,
C     bypassing the use of the fort.13 file for this purpose.
C     jgf47.30 Changed to emphasize full domain time varying tau0
      IF ( Tau0Dig1.eq.-5 ) THEN
         FullDomainTimeVaryingTau0 = .True.
         IF ( .not.LoadTau0MinMax ) THEN
            DO I=1, NP
               Tau0MinMax(I,1) = Tau0FullDomainMin
               Tau0MinMax(I,2) = Tau0FullDomainMax
            ENDDO
         ENDIF
      ENDIF
C
C     jgf47.30: Output of tau0 is now activated by having a 0.1 fraction
C     for tau0
      IF ( Tau0Dig2.eq.-1 ) THEN
         OutputTau0 = .True.
      ENDIF
C
C     jjw&sb46.38.sb01 If tau0 is loaded from nodal attributes file and
C     Tau0 is -3, time-varing tau0 optimizer will be applied in timestep.F
      IF (HighResTimeVaryingTau0 .or. FullDomainTimeVaryingTau0) THEN
         ALLOCATE(Tau0Temp(NP))
         WRITE(16,8) ! jgf47.30.TODO: This logging should be consolidated.
      ENDIF
 8    FORMAT(/,5X,'A SPATIALLY TEMPORALLY VARIABLE OPTIMIZED '
     &     ,/,5X,' WEIGHTING COEFFICIENT (Tau0) WILL BE USED '
     &     ,/,5X,' IN THE GENERALIZED WAVE CONTINUITY EQUATION.',/)
C
C     jgf47.33 Enable time averaging of tau0 if requested.
      IF (Tau0Dig1.eq.-6) THEN
         TimeAveragedTau0 = .true.
      ENDIF
C
C     jgf48.42 Enable back loaded time averaging of tau0 if requested.
      IF (Tau0Dig1.eq.-7) THEN
         BackLoadedTimeAveragedTau0 = .true.
      ENDIF
C
C     jgf48.46 Allocate array to hold previous value tau0 for use in
C     time averaging, if necessary.
      IF ( TimeAveragedTau0 .or. BackLoadedTimeAveragedTau0 ) THEN
         ALLOCATE(LastTau0(NP))
         DO I=1, NP
            LastTau0(I) = Tau0Base(I)
         ENDDO
      ENDIF
C
C     I N I T   B O T T O M   F R I C T I O N
      IF(NOLIBF.EQ.0) THEN
         IFNLBF=0
         IFLINBF=1
         IFHYBF=0
      ENDIF
      IF(NOLIBF.EQ.1.or.NOLIBF.EQ.3) THEN
         IFNLBF=1
         IFLINBF=0
         IFHYBF=0
      ENDIF
      IF(NOLIBF.EQ.2) THEN
         IFNLBF=0
         IFLINBF=0
         IFHYBF=1
      ENDIF
C
C     Initialize bottom friction if it was not loaded from unit 13.
      IF(C2DDI) THEN
         IF((.not.LoadQuadraticFric).and.(.not.LoadManningsN).and.
     &     (.not.LoadChezy)) THEN
            IF (NoLiBF.eq.0) CF=Tau
C     If a nodal attributes file was read, FRIC was allocated there.
            IF (NWP.eq.0) THEN
               ALLOCATE(FRIC(NP))
               ALLOCATE(FRIC_ORIG(NP))
            ENDIF
            DO I=1,NP
               FRIC(I)=CF
            END DO
            FRIC_ORIG = FRIC
         ENDIF
C
C        jgf47.04 If a depth-dependent friction parameterization is used,
C        the value from the fort.15 file is used as a floor for the
C        minimum equivalent quadratic friction value.
         IF (LoadManningsN) THEN
            BFCdLLimit = CF
         ENDIF
      ENDIF
C
      IF(C3D) THEN
C     Initialize 3D bottom roughness if it was not loaded from unit 13.
        IF((.not.LoadZ0b_var)) THEN
C     If a nodal attributes file was read, Z0b_var was allocated there.
          IF (NWP.eq.0) THEN
             ALLOCATE(Z0b_var(NP))
             ALLOCATE(FRIC(NP))
             ALLOCATE(FRIC_ORIG(NP))
          ENDIF
          DO I=1,NP
             Z0b_var(I)=Z0B
             FRIC(I)=CF
          END DO
          FRIC_ORIG = FRIC
        ENDIF
C     jgf47.04 If a depth-dependent friction parameterization is used,
C     the value from the fort.15 file is used as a floor for the
C     minimum equivalent quadratic friction value.
         IF (LoadZ0b_var.OR.LoadManningsN) THEN
            BFCdLLimit = CF
         ENDIF
      ENDIF
C
C     Initialize bridge pilings.
C
      IF (LoadBridgePilings) THEN
         DO I=1, NP
            IF (BridgePilings(I,1).ne.0) THEN ! only for nodes w/piers
               BridgePilings(I,3) = 4.d0 *
     &              BridgePilings(I,3) / BridgePilings(I,4)
            ENDIF
         END DO
      ENDIF
C
C     I N I T   E D D Y   V I S C O S I T Y  &  D I F F U S I V I T Y
      IF (.not.LoadEVM) THEN
         IF (NWP.eq.0) THEN
            ALLOCATE(EVM(NP))
         ENDIF
         DO I=1,NP
            EVM(I)=ESLM
         END DO
      ENDIF
      IF (.not.LoadEVC.and.ESLC.ne.0) THEN
         IF (NWP.eq.0) THEN
            ALLOCATE(EVC(NP))
         ENDIF
         DO I=1,NP
            EVC(I)=ESLC
         END DO
      ENDIF
C...  DW 
      IF ( NWP .eq. 0 ) THEN
         NumNodesAbsLayer = 0
      END IF
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE InitNodalAttr
C     ----------------------------------------------------------------


C     ----------------------------------------------------------------
C      F U N C T I O N   T A U 0  N O D A L  V A L U E
C     ----------------------------------------------------------------
C
C     jgf46.27 Function to calculate tau0 based on the scheme selection
C     and the depth. This assumes that Scheme is negative.
C
C     ----------------------------------------------------------------
      REAL(SZ) FUNCTION Tau0NodalValue(Scheme, Depth)
      IMPLICIT NONE
      REAL(SZ) Scheme
      REAL(SZ) Depth
C
      IF (Scheme.eq.-2.d0) THEN
C     Smoothly varying tau0 with depth.
         IF(Depth.GE.200.) Tau0NodalValue=0.005
         IF((Depth.LT.200.).AND.(Depth.GE.1.)) THEN
            Tau0NodalValue=1./Depth
         ENDIF
         IF(Depth.LT.1.) Tau0NodalValue=1.0
      ELSE
C     Abrupt variation in tau0 with depth.
         IF(Depth.LE.10.) Tau0NodalValue=0.020d0
         IF(Depth.GT.10.) Tau0NodalValue=0.005d0
      ENDIF
C     ----------------------------------------------------------------
      END FUNCTION Tau0NodalValue
C     ----------------------------------------------------------------


C     ----------------------------------------------------------------
C      S U B R O U T I N E
C         C A L C U L A T E  T I M E  V A R Y I N G  T A U 0
C     ----------------------------------------------------------------
C
C     jgf47.08 Subroutine to calculate a new tau0 value. Called from
C     GWCE_New in timestep.F each time the GWCE matrix is reset (i.e.,
C     upon startup and whenever wetting and/or drying occurs in any
C     subdomain. Based on Casey050711.
C
C     ----------------------------------------------------------------
      SUBROUTINE CalculateTimeVaryingTau0(TK, NNeigh, NeiTab, NP)
      IMPLICIT NONE
      REAL(SZ), intent(in) :: TK(:)       ! bottom friction
      INTEGER, intent(in) :: NNeigh(:)   ! number of neighbor nodes
      INTEGER, intent(in) :: NeiTab(:,:) ! table of neighbor nodes
      INTEGER, intent(in) :: NP           ! number of nodes in the domain
      REAL(SZ) CaseySum  ! sum of tau0temp values around a particular node

Casey 050711 : Made changes for averaged variable Tau0.
cjjw46.39.sb01 :  "high/low LIMITED" variable G.
      !jgf47.30: Distinction between fulldomain and hi res only
      IF ( FullDomainTimeVaryingTau0 ) THEN
         DO i = 1, NP
            Tau0Temp(i)=Tau0MinMax(i,1)+1.5*TK(i)
            IF (Tau0Temp(i).lt.Tau0MinMax(i,1)) THEN
               Tau0Temp(i)=Tau0MinMax(i,1)
            ENDIF
            IF(Tau0Temp(i).gt.Tau0MinMax(i,2)) THEN
               Tau0Temp(i)=Tau0MinMax(i,2)
            ENDIF
         ENDDO
      ENDIF
      IF ( HighResTimeVaryingTau0 ) THEN
         DO i = 1, NP
            IF(Tau0Base(i).lt.0.025) THEN
               Tau0Temp(I)=Tau0Base(i) ! not time varying
            ELSE
               Tau0Temp(i)=Tau0Base(i)+1.5*TK(i) ! time varying
               IF (Tau0Temp(i).gt.0.2) Tau0Temp(i)=0.2 ! ceiling
            ENDIF
         ENDDO
      ENDIF
      ! smoothing
      DO I=1, NP
         CaseySum = 0.0
         DO J=1,NNeigh(I)
            CaseySum = CaseySum + Tau0Temp(NeiTab(I,J))
         ENDDO
         TAU0VAR(I) = CaseySum / NNeigh(I)
      ENDDO
C
C     jgf47.33 Perform time averaging of tau0 if requested.
      IF (TimeAveragedTau0) THEN
         DO I=1, NP
            TAU0VAR(I) = 0.5d0*TAU0VAR(I) + 0.5d0*LastTau0(I)
            LastTau0(I) = TAU0VAR(I)
         ENDDO
      ENDIF
C
C     jgf48.42 Perform backloaded time averaging of tau0 if requested.
      IF (BackLoadedTimeAveragedTau0) THEN
         DO I=1, NP
            TAU0VAR(I) = AlphaTau0*TAU0VAR(I)
     &                 + (1.d0-AlphaTau0)*LastTau0(I)
            LastTau0(I) = TAU0VAR(I)
         ENDDO
      ENDIF
C     ----------------------------------------------------------------
      END SUBROUTINE CalculateTimeVaryingTau0
C     ----------------------------------------------------------------



C     ----------------------------------------------------------------
C     S U B R O U T I N E
C     A P P L Y  2 D  B O T T O M  F R I C T I O N
C     ----------------------------------------------------------------
C
C     jgf46.00 Subroutine to apply 2D bottom friction from turbulent
C     viscous effects as well as bridge pilings. This is used in the
C     time stepping loop.
C
C     ----------------------------------------------------------------
C
C     sb46.28sb02 Lower limit of Cd was added as an argument.
C     jgf47.04 Argument for lower limit of Cd was removed; this value
C     is now specified by the user in fort.15.
C
C     ----------------------------------------------------------------
      SUBROUTINE Apply2DBottomFriction(UU1, VV1, DP, ETA2, G,
     &     IFNLFA, NP, TK)
      USE GLOBAL, ONLY: MLD, DTDP
      USE INTERNALTIDE, ONLY: UBar, VBar
C      USE, INTRINSIC :: IEEE_ARITHMETIC !jgfdebug ieee_is_nan()
      IMPLICIT NONE
      INTEGER, intent(in) :: NP                   ! number of nodes in grid
      REAL(SZ), intent(in), dimension(NP) :: UU1  ! x-dir velocities
      REAL(SZ), intent(in), dimension(NP) :: VV1  ! y-dir velocities
      REAL(SZ), intent(in), dimension(NP) :: DP   ! bathymetric depths
      REAL(SZ), intent(in), dimension(NP) :: ETA2 ! water surf. elevations
      REAL(SZ), intent(in) :: G                   ! gravitational constant
      INTEGER, intent(in) :: IFNLFA               ! nonlin. finite amp. flag
      REAL(SZ), intent(inout), dimension(NP) :: TK! depth avg. fric.
C
      REAL(SZ) UV1   ! velocity magnitude (speed)
      REAL(SZ) H1    ! total depth
      REAL(SZ) Fr
      REAL(SZ) FricBP
      REAL(SZ) BK    ! BK(1) is pier shape factor
      REAL(SZ) BALPHA! BALPHA(2) is constriction fraction
      REAL(SZ) BDELX ! BDELX(3) is effective delx
      REAL(SZ) MLDrat ! MLD ratio to depth
      REAL(SZ) HT, WCH, ManN


      !...Reset the friction to original values
      !   when subgrid barriers might be used
      IF(LoadSubgridBarrier)THEN
        DO I=1,NP
            IF(subgridBarrierNoOfVals.EQ.1)CYCLE
            IF(LoadManningsN.AND.ManningsN(I).GT.0D0)CYCLE
            IF(subgridbarrierOvertopping(I))THEN
                FRIC(I) = subgridBarrier(I,2)
            ELSE
                FRIC(I) = FRIC_orig(I)
            ENDIF
        ENDDO
      ENDIF

C
C     Step 0. Convert Manning's N to Cd, if necessary.
      IF (LoadManningsN) THEN
         DO I=1, NP
            HT = DP(I)+IFNLFA*ETA2(I)
            IF(LoadSubgridBarrier)THEN
                IF(subgridBarrierOvertopping(I).AND.
     &             SubgridBarrierNoOfVals.EQ.2)THEN
                    FRIC(I)=g*subgridBarrier(I,2)**2.d0
     &                   /( ( HT )**(1.d0/3.d0) ) 
                    IF(FRIC(I).LT.BFCdLLimit) THEN
                       FRIC(I) = BFCdLLimit
                    ENDIF
                ELSE
                    IF (ManningsN(I).le.0d0) cycle ! Can use Cd where man <= 0
                    FRIC(I)=g*ManningsN(I)**2.d0
     &                   /( ( HT )**(1.d0/3.d0) ) ! sb46.28sb02
                    IF(FRIC(I).LT.BFCdLLimit) THEN
                       FRIC(I) = BFCdLLimit
                    ENDIF
                ENDIF
            ELSE
                IF (ManningsN(I).le.0d0) cycle ! Can use Cd where man <= 0
                IF (ManningsN(I).gt.1.d0) THEN
                   WCH = floor(ManningsN(I))
                   ManN = ManningsN(I) - WCH
                   FRIC(I)=g*((1.D0+2.D0*HT/WCH)**(2.D0/3.D0)*ManN)**2.d0
     &               /( HT**(1.d0/3.d0) )
                   IF(FRIC(I).LT.BFCdLLimit) THEN
                      FRIC(I) = BFCdLLimit
                   ENDIF
                ELSE
                   FRIC(I)=g*ManningsN(I)**2.d0
     &               /( ( HT )**(1.d0/3.d0) ) ! sb46.28sb02
                   !sb46.28sb02  Lower limit is applied here.
                   IF(FRIC(I).LT.BFCdLLimit) THEN
                      FRIC(I) = BFCdLLimit
                   ENDIF
                ENDIF
            ENDIF
         ENDDO
      ELSEIF (LoadChezy) THEN
         DO I=1,NP
            FRIC(I) = G/(Chezy(I)**2)
         END DO
      ENDIF

C
      NOLIBFRIC: IF (NOLIBF.eq.3) THEN
C     ... Combined Cd for bottom friction and momentum dispersion
         DO I = 1,NP
            UV1 = SQRT(UU1(I)*UU1(I)+VV1(I)*VV1(I))
            H1  = DP(I) + IFNLFA*ETA2(I)
            IF (HBREAK.gt.0d0) then
C              ... Cd as a power function [a(x^b-1)] of x = MLD
               FRIC(I) = CF + HBREAK*(max(4d-3,MLD(I))**FTHETA - 1d0)
            ELSE
C              ... Cd directly from the 3D model (MLD here)
               FRIC(I) = CF + MLD(I)
            ENDIF
            TK(I) = MIN(FRIC(I) * UV1 / H1, 2.0/DTDP)
         ENDDO
C
      ELSE
C
C     .. Step 1. Apply friction arising from turbulent viscous interaction
C        with the sea floor.
         DO I=1, NP
            UV1 = SQRT(UU1(I)*UU1(I)+VV1(I)*VV1(I))
            H1  = DP(I) + IFNLFA*ETA2(I)
            TK(I) = MIN(FRIC(I)*
     &                     (IFLINBF +                                          ! linear
     &                     (UV1/H1) * (IFNLBF +                               ! nonlinear
     &                     IFHYBF*(1+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)) ), ! hybrid 
     &                  2.0/DTDP)
         END DO
C
      ENDIF NOLIBFRIC
C
C     Step 2. Apply friction arising from flow interaction with bridge
C     pilings, if required.
      IF (LoadBridgePilings) THEN 
         DO I=1, NP
            UV1=SQRT(UU1(I)*UU1(I)+VV1(I)*VV1(I))
            H1=DP(I)+IFNLFA*ETA2(I)
            Fr=UV1*UV1/(G*H1)
            BK = BridgePilings(I,1)
            BALPHA = BridgePilings(I,2)
            BDELX = BridgePilings(I,3)
            FricBP=(H1/BDELX)*BK*(BK+5.d0*Fr*Fr-0.6d0)
     &           *(BALPHA+15.d0*BALPHA**4)
            TK(I) = MIN(TK(I)+FricBP*UV1/H1, 2.0/DTDP)
         END DO
      ENDIF
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE Apply2DBottomFriction
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
C     S U B R O U T I N E
C     A P P L Y  3 D  B O T T O M  F R I C T I O N
C     ----------------------------------------------------------------
C
C     jgf46.00 Subroutine to apply 3D bottom friction from turbulent
C     viscous effects as well as bridge pilings. This is used in the
C     time stepping loop.
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
      SUBROUTINE Apply3DBottomFriction(Q, SIGMA, DP, ETA2, G,
     &     IFNLFA, NP, TK, NFEN, Z0B)
      IMPLICIT NONE
      INTEGER, intent(in) :: NP, NFEN             ! number of nodes in grid Horizontal and Vertical
      COMPLEX(SZ), intent(in), dimension(NP,NFEN) :: Q  ! x-dir velocities
      REAL(SZ), intent(in), dimension(NFEN) :: SIGMA  ! x-dir velocities
      REAL(SZ), intent(in), dimension(NP) :: DP   ! bathymetric depths
      REAL(SZ), intent(in), dimension(NP) :: ETA2 ! water surf. elevations
      REAL(SZ), intent(in) :: G                   ! gravitational constant
      INTEGER, intent(in) :: IFNLFA               ! nonlin. finite amp. flag
      REAL(SZ), intent(inout), dimension(NP) :: TK! depth avg. fric.
      REAL(SZ), intent(in) :: Z0B
C
      INTEGER NH
      REAL(SZ) Z0B1  ! velocity magnitude (speed)
      REAL(SZ) UV1   ! velocity magnitude (speed)
      REAL(SZ) H1    ! total depth
      REAL(SZ) Fr
      REAL(SZ) FricBP
      REAL(SZ) BK    ! BK(1) is pier shape factor
      REAL(SZ) BALPHA! BALPHA(2) is constriction fraction
      REAL(SZ) BDELX ! BDELX(3) is effective delx
C
C Determine the bottom roughness length either from fort.15, from Manning's n
C or as read in from nodal attributes
      DO NH=1,NP
         H1=DP(NH)+IFNLFA*ETA2(NH)
         IF (LoadZ0B_var) THEN
            Z0B1 = Z0B_var(NH)
         ELSEIF (LoadManningsN) THEN
            Z0B1 = ( H1 )* exp(-(1.0D0+
     &             ( (0.41D0*( H1 )**(1.0D0/6.0D0) )/
     &                            (ManningsN(NH)*sqrt(g)) ) ))
         ELSE
            Z0B1 = Z0B
         ENDIF

         FRIC(NH)= (1.D0 / ( (1.D0/0.41D0) *
     &   LOG((ABS( ( ( SIGMA(2)-SIGMA(1) )/2.d0 ) *(H1) ) + Z0B1 )/Z0B1)
     &       ) )**2.D0

         TK(NH)= FRIC(NH) * ABS(Q(NH,1))

         IF (LoadBridgePilings) THEN
            Fr=ABS(Q(NH,1))*ABS(Q(NH,1))/(G*H1)
            BK = BridgePilings(I,1)
            BALPHA = BridgePilings(I,2)
            BDELX = BridgePilings(I,3)
            FricBP=(H1/BDELX)*BK*(BK+5.d0*Fr*Fr-0.6d0)
     &           *(BALPHA+15.d0*BALPHA**4)
            TK(I)=TK(I)+FricBP*ABS(Q(NH,1))/H1
         ENDIF
      ENDDO
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE Apply3DBottomFriction
C     ----------------------------------------------------------------

C     ----------------------------------------------------------------
C     S U B R O U T I N E
C     A P P L Y  2 D  I N T E R N A L  W A V E  D R A G
C     ----------------------------------------------------------------
C
C     WJP March 30 2016: Adjusted for TKM input matrix and inserting
C     the scalar version in here too
C     WJP 02.24.2018 put into v53.dev
C     WJP 07.24.2018 updated for baroclinic coupling type 
C     ----------------------------------------------------------------
      SUBROUTINE Apply2DInternalWaveDrag(NP, TK, TKM, U, V, DP,
     & IT)      
      USE GLOBAL, ONLY: CBaroclinic, CORIF, DTDP, NB, NM, IDEN,
     &                  IFSPROTS, RAMP
      USE MESH, ONLY: DRVMAP2DSPVEC, UVECTMP, VVECTMP, RVELF
      USE INTERNALTIDE, ONLY: UNTIDE, UBar, VBar
      IMPLICIT NONE
      INTEGER, intent(in) :: NP                   ! number of nodes in grid
      REAL(SZ), intent(in), dimension(:) :: TK   ! scalar friction
      REAL(SZ), intent(inout), dimension(:,:) :: TKM
      INTEGER,  intent(in), optional :: IT          ! temporal interations
      REAL(SZ), intent(in), dimension(NP),optional :: U    ! x-dir velocities
      REAL(SZ), intent(in), dimension(NP),optional :: V    ! y-dir velocities
      REAL(SZ), intent(in), dimension(NP),optional :: DP   ! depth
      REAL(SZ), allocatable :: Hx(:), Hy(:)
      INTEGER  :: I, J 
      REAL(SZ) :: omega = 1.4052d-4 ![rad/s] M2 tidal frequency
      REAL(SZ) :: hr = 3600.0_SZ ! means one hour in seconds
      REAL(SZ) :: SmallV = 1d-10 ! A small value
      REAL(SZ) :: Fit, alpha2, gamma2, slope2
      REAL(SZ) :: Nmdiff, Nbdiff, Fdiff, adju, adjv, adjuv
      logical, save ::  first_entry = .TRUE.
      ! CPB: for tidal velocities
      REAL(SZ) :: cosAMIGT, sinAMIGT, DissTidal, DissTotal, jjwgamma,
     &            Cxx, Cyy, Cxy
      INTEGER :: NCYC
      REAL(SZ),DIMENSION(:),ALLOCATABLE :: Udmy, Vdmy      
C     Compute TKM by adding wave drag contribution to 
C     current value of the bottom friction, TK
      LoadITFric: IF (LoadIT_Fric) THEN
         IF ( FIRST_ENTRY ) THEN
            FIRST_ENTRY = .FALSE.
            IF ( IFSPROTS.eq.1 ) THEN
C              Rotate the IT_Fric first time step if mesh rotation used
C              (unnecessary for IT_FricNoOfVals == 1)
               IF ( IT_FricNoOfVals.eq.3 ) THEN
                  allocate(Hx(NP),Hy(NP))
                  DO I = 1,NP
                     if (CBaroclinic.and.abs(IDEN).ge.5) then
                        Hx(I) = IT_Fric(I,2)
                        Hy(I) = IT_Fric(I,3)
                     else
                        Hx(I) = sqrt(IT_Fric(I,1))
                        Hy(I) = sign(sqrt(IT_Fric(I,2)),IT_Fric(I,3))
                     endif
                     UVECTMP(I) = Hx(I); 
                     VVECTMP(I) = Hy(I);
                  ENDDO
                  CALL DRVMAP2DSPVEC( Hx, Hy,
     &                UVECTMP, VVECTMP, NP, FWD = .TRUE. ) ;
                  DO I = 1,NP
                     if (CBaroclinic.and.abs(IDEN).ge.5) then
                        IT_Fric(I,2) = Hx(I)
                        IT_Fric(I,3) = Hy(I)
                     else
                        IT_Fric(I,1) = Hx(I)*Hx(I)
                        IT_Fric(I,2) = Hy(I)*Hy(I)
                        IT_Fric(I,3) = Hx(I)*Hy(I)
                     endif
                  ENDDO           
                  deallocate(Hx,Hy)
               ELSEIF ( IT_FricNoOfVals.eq.5 ) THEN
                  ! need do something else
               ENDIF
            ENDIF
         ENDIF
         Baroc: IF (CBaroclinic.and.abs(IDEN).ge.5) THEN
C           2DDI ADCIRC with baroclinic coupling
            if (abs(IDEN).ne.7) then
               CALL UNTIDE(U,V,IT)
            endif
            IF (IT_FricNoOfVals == 5) THEN
C           Use the "Nonlocal" Nycander (2005) type tensor wave drag
C           1st column is coefficient, 2nd and 3rd columns are the Hx and Hy slopes,
C           4th and 5th columns are the Jx and Jy slopes
               DO I = 1,NP
                  IF (omega > NB(I).or.omega < CORIF(I)) THEN
                     ! No freely propagating internal tides
                     TKM(1:2,I) = TK(I)
                     TKM(3,I)   = 0.0d0
                  ELSE 
                     ! Have freely propagating internal tides
                     ! Get coefficient
                     Fit = IT_Fric(I,1)*NB(I)/DP(I)*
     &                     sqrt(1d0 - CORIF(I)**2/omega**2)
                     ! Check criticality of slope
                     alpha2 = (omega**2 - CORIF(I)**2) / 
     &                        (NB(I)**2 - omega**2)
                     gamma2 = (IT_Fric(I,2)**2 + IT_Fric(I,3)**2)/alpha2
                     ! Divide by criticality if > 1 
                     Fit = Fit/max(1d0,gamma2)                   
                     TKM(1,I) = TK(I) + 2*Fit*IT_Fric(I,2)*IT_Fric(I,4)
                     TKM(2,I) = TK(I) + 2*Fit*IT_Fric(I,3)*IT_Fric(I,5)
                     TKM(3,I) = Fit*(IT_Fric(I,2)*IT_Fric(I,5) +  
     &                               IT_Fric(I,3)*IT_Fric(I,4))
                  ENDIF 
               ENDDO
            ELSEIF (IT_FricNoOfVals == 3.or.IT_FricNoOfVals == 1) THEN
C           For the "Local" Lyard et al. (2004) type tensor wave drag
C           1st column is coefficient, 2nd and 3rd columns are the Hx and Hy slopes
C           For the scalar Zaron & Egbert (2006) type
C           1st column is coefficient, 2nd column is the squared hypotenuse slope
               DO I = 1,NP
                  adju = 1d0; adjv = 1d0;
                  NBdiff = (NB(I)**2d0 - omega**2d0)
                  NMdiff = (NM(I)**2d0 - omega**2d0)
                  Fdiff  = (omega**2d0 - CORIF(I)**2d0)
                  IF (NBdiff.le.SmallV.or.Nmdiff.le.SmallV.or.
     &                Fdiff.le.SmallV) THEN
                     ! No freely propagating internal tides
                     TKM(1:2,I) = TK(I)
                     TKM(3,I)   = 0.0d0
                  ELSE 
                     ! Have freely propagating internal tides
                     ! Get the full coefficient
                     IF (NBdiff*NMdiff.LT.0d0) THEN
                        Fit = 0d0
                     ELSE
                        Fit = IT_Fric(I,1)*sqrt(NBdiff*Nmdiff)/omega
                     ENDIF
                     ! Check criticality of slope
                     alpha2 = Fdiff / NBdiff
                     IF (alpha2.LT.0d0) Fit = 0d0
                     if (IT_FricNoOfVals == 3) then
                        slope2 = IT_Fric(I,2)**2d0 + IT_Fric(I,3)**2d0
                     else
                        slope2 = IT_Fric(I,2)
                     endif
                     gamma2 = slope2 / alpha2
                     ! Divide coefficient by criticality if > 1
                     ! (saturates out at high criticalities) 
                     Fit = Fit / max(1d0,gamma2)                   
                     ! Calculate adjustment based on tidal residual
                     if (abs(IDEN).ne.7) then
                        ! Calculate the adjustment factors based on the delayed moving mean
                        IF (ABS(U(I)).GT.1d-3.AND.
     &                      ABS(V(I)).GT.1d-3) THEN
                           adju = 1d0-UBar(I)/U(I)
                           adjv = 1d0-VBar(I)/V(I)
                        ENDIF
                        adjuv = 0.5d0*(adju + adjv)
                     endif
                     ! Implement the friction with possible adjustment
                     IF (IT_FricNoOfVals == 3) then
                        TKM(1,I) = TK(I) + adju*Fit*IT_Fric(I,2)**2d0
                        TKM(2,I) = TK(I) + adjv*Fit*IT_Fric(I,3)**2d0
                        TKM(3,I) = adjuv*Fit*IT_Fric(I,2)*IT_Fric(I,3)
                     else
                        TKM(1:2,I) = TK(I) + adjuv*Fit*IT_Fric(I,2)
                        TKM(3,I)   = 0.0d0
                     endif
                  ENDIF 
               ENDDO
            endif
         ELSE
C           Normal 2DDI ADCIRC without baroclinic coupling
            IF (IT_FricNoOfVals == 3) THEN
C              Use the "Nonlocal" Nycander (2005) or "Local" Lyard et al. (2004) 
C              type tensor wave drag
C              Using three columns of data for the different segments of the matrix
               DO I=1,NP 
                  TKM(1,I) = TK(I) + IT_Fric(I,1)
                  TKM(2,I) = TK(I) + IT_Fric(I,2)
                  TKM(3,I) = IT_Fric(I,3) 
               END DO
            ELSEIF (IT_FricNoOfVals == 1) THEN
C              For the scalar type Zaron & Egbert (2006) type
               DO I=1,NP 
                  TKM(1,I) = TK(I) + IT_Fric(I,1) 
                  TKM(2,I) = TKM(1,I)
                  TKM(3,I) = 0.0d0
               ENDDO
            ENDIF
         ENDIF Baroc
      ELSE 
C       In case we did not load ITFric
        TKM(1,:) = TK
        TKM(2,:) = TK
        TKM(3,:) = 0.0d0
      ENDIF LoadITFric

      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE Apply2DInternalWaveDrag
C     ----------------------------------------------------------------

!     ----------------------------------------------------------------
!>    @brief Subroutine to calculate the land wind reduction factor
!>    based on a table of directional wind drag values. Originally
!>    written into the hstart.F file by jjw in jjw-42.06j. This is used
!>    in hstart.F and timestep.F.
!>
!>    The methology comes from: 
!>    Simiu, E. and R. Scanlan, 1986: Wind Effects on Structures. Wiley Interscience, 604 pp.
!>
!>    Additionally, the z0 coefficient supplied by the user is modified
!>    by subtracting a factor of depth/30 to account for inundation
!>    
!>    Values are linearly interpolated between heading bins to the actual 
!>    wind direction.  
!>
!>    @param[in] NodeNumber index of node under consideration
!>    @param[in] WindDragCo current wind drag coefficient being applied
!>    @param[in] WindMag    current wind speed magnitude at the specified node
!>    @param[in] BathymetricDepth bathymetric depth used to calculate total depth
!>    @param[in] Elevation water surface elevation above datum (eta2)
!>    @param[in,out] WindX wind velocity, x-component
!>    @param[in,out] WindY wind velocity, y-component
!     ----------------------------------------------------------------
      SUBROUTINE ApplyDirectionalWindReduction(NodeNumber, WindDragCo,
     &     WindMag, BathymetricDepth, Elevation, WindX, WindY)
        USE GLOBAL,ONLY:h0
        USE CONSTANTS,ONLY:RAD2DEG,G
        implicit none
        intrinsic               :: atan2
        integer,  intent(in)    :: NodeNumber 
        real(SZ), intent(in)    :: WindDragCo 
        real(SZ), intent(in)    :: WindMag    
        real(SZ), intent(in)    :: BathymetricDepth 
        real(SZ), intent(in)    :: Elevation  
        real(SZ), intent(inout) :: WindX   
        real(SZ), intent(inout) :: WindY  
        
        real(SZ),parameter      :: CutoffDepth = 999999.d0          !< Depth above which reduction is not applied
        real(SZ),parameter,dimension(12) :: z0Angles=(/(i*30,i=0,5),(i*30-180,i=0,5)/)  !< centers of directional bins
        
        integer                 :: idir        !< code for wind direction bins to interpolate between
        integer                 :: idir2       !< code for wind direction bins to interpolate between
        real(SZ)                :: z0m         !< marine roughness coefficient based on Garratt's formula
        real(SZ)                :: angle       !< direction wind is coming from
        real(SZ)                :: z0l         !< drag for a particular node, for particular direction
        real(SZ)                :: TotalDepth  !< bathymetric depth + sea surface elevation
        real(SZ)                :: fr          !< land wind reduction factor

!       if windspeed is zero, exit now
        if((WindX.eq.0d0).and.(WindY.eq.0d0))return

!       compute direction  that the wind is going to
        angle=atan2(WindY,WindX)*RAD2DEG
        idir=0
        idir2=0
        if(angle.lt.-150d0)then
            idir=7
            idir2=8
        elseif(angle.lt.-120d0)then
            idir=8
            idir2=9
        elseif(angle.lt. -90d0)then
            idir=9
            idir2=10
        elseif(angle.lt. -60d0)then
            idir=10
            idir2=11
        elseif(angle.lt. -30d0)then
            idir=11
            idir2=12
        elseif(angle.lt.   0d0)then
            idir=12
            idir2=1
        elseif(angle.lt.  30d0)then
            idir=1
            idir2=2
        elseif(angle.lt.  60d0)then
            idir=2
            idir2=3
        elseif(angle.lt.  90d0)then
            idir=3
            idir2=4
        elseif(angle.lt. 120d0)then
            idir=4
            idir2=5
        elseif(angle.lt. 150d0)then
            idir=5
            idir2=6
        elseif(angle.le. 180d0)then
            idir=6
            idir2=7
        endif


!       compute marine roughness coefficient based on Garratt's formula
        z0m=(0.018d0/G)*WindDragCo*WindMag**2.d0

        !tga 2020-04 updated code to linearly interpolate z0 values, 
        !code previously did nearest neighbor (i.e. binned) z0 values.  
        !Define roughness by linearly (in angle-z0 space) interpolating
        !between binned values.  
        !The 30d0 here comes from assuming the directional bins are 30 
        !degrees apiece. 
        z0l=(z0land(NodeNumber,idir2)-z0land(NodeNumber,idir))/30d0
     &      *(angle-z0Angles(idir))+z0land(NodeNumber,idir)

!       apply overland flooding correction
        TotalDepth = BathymetricDepth + Elevation
        if( (TotalDepth.gt.2D0*h0).and.
     &      (BathymetricDepth.lt.CutOffDepth)) then
           z0l=z0l-TotalDepth*OverlandReductionFactor(NodeNumber)
        endif

!       compute land wind reduction factor
!       Reduction factor is bounded to not exceed 1, i.e. assumes the 
!       land roughness is never less than the water.  
        if(z0l.gt.z0m) then
           fr=min(1.0D0,(z0l/z0m)**0.0706d0 * log(10.d0/z0l) / log(10.d0/z0m))
        else
           fr=1.000d0
        endif

!       adjust time interpolated wind field
        WindX = fr*WindX
        WindY = fr*WindY
        
        RETURN
!     ----------------------------------------------------------------
      END SUBROUTINE ApplyDirectionalWindReduction
!     ----------------------------------------------------------------


C     ----------------------------------------------------------------
C     S U B R O U T I N E
C     A P P L Y   C A N O P Y   C O E F F I C I E N T
C     ----------------------------------------------------------------
C     jgf49.1001 Subroutine to apply the canopy coefficient. Was
C     originally included in the subroutine ApplyDirectionalWindReduction;
C     this combination implicitly assumed that the two attributes would
C     always be used together.
C     ----------------------------------------------------------------
      SUBROUTINE ApplyCanopyCoefficient(NodeNumber, WindX, WindY)
      IMPLICIT NONE
      INTEGER,  intent(in) :: NodeNumber ! index of node under consideration
      REAL(SZ), intent(inout) :: WindX   ! x-dir component of wind velocity
      REAL(SZ), intent(inout) :: WindY   ! x-dir component of wind velocity
C
      WindX = vcanopy(NodeNumber)*WindX
      WindY = vcanopy(NodeNumber)*WindY
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE ApplyCanopyCoefficient
C     ----------------------------------------------------------------

C*************************************************************************
C     Subroutine to determine local advection state for each element
C
C                                   Corbitt 120322
C*************************************************************************
      SUBROUTINE ADVECTLOCAL(IE)
      USE SIZES, ONLY : SZ
      USE GLOBAL, ONLY : IFNLCAT, IFNLCATE, IFNLCT, IFNLCTE
      USE MESH, ONLY : NM, DP
      Implicit None
      Integer  :: NM1, NM2, NM3,IE

      NM1=NM(IE,1)
      NM2=NM(IE,2)
      NM3=NM(IE,3)
      IF ((DP(NM1).GE.AdvectionState(NM1)).AND.
     &    (DP(NM2).GE.AdvectionState(NM2)).AND.
     &    (DP(NM3).GE.AdvectionState(NM3))) THEN
             IFNLCT = IFNLCTE
             IFNLCAT = IFNLCATE
      ELSE
         IFNLCT = 0
         IFNLCAT = 0
      ENDIF
C***********************************************************************
      END SUBROUTINE ADVECTLOCAL
C***********************************************************************


C     ----------------------------------------------------------------
C                      S U B R O U T I N E
C         R E A D  L E G A C Y  S T A R T  D R Y  F I L E
C     ----------------------------------------------------------------
C
C     jgf46.00 Subroutine to load up the legacy startdry file (unit
C     12). This is just a cut-and-paste from the section of the
C     READ_INPUT subroutine that did the same thing. This subroutine is
C     never called. It is vestigial and listed here purely as reference
C     material.
C
C     ----------------------------------------------------------------
      SUBROUTINE ReadLegacyStartDryFile(NP, NScreen, ScreenUnit,
     &     MyProc, NAbOut)
      IMPLICIT NONE
      INTEGER, intent(in) :: NP ! number of nodes in grid file
      INTEGER, intent(in) :: NScreen ! nonzero for debug info to screen
      INTEGER, intent(in) :: ScreenUnit ! i/o for debug info to screen
      INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
      INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16

      INTEGER JKI          ! node number from file
      INTEGER NE2          ! number of elements, according to fort.12 file
      INTEGER NP2          ! number of nodes, according to fort.12 file

      CHARACTER(len=80) AGRID2 ! users comment/description line
      REAL(SZ) DUM1, DUM2  ! data that we want to skip

      OPEN(12,FILE=TRIM(INPUTDIR)//'/'//'fort.12')
C
C...  READ STARTDRY INFORMATION FROM UNIT 12
      READ(12,'(A80)') AGRID2
      WRITE(16,2038) AGRID2
2038  FORMAT(5X,'STARTDRY FILE IDENTIFICATION : ',A80,/)
      READ(12,*) NE2,NP2
C
C...  CHECK THAT NE2 AND NP2 MATCH WITH GRID FILE
C      IF((NE2.NE.NE).OR.(NP2.NE.NP)) THEN
       IF(NP2.NE.NP) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9900)
         WRITE(16,9900)
 9900    FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!',
     &        //,1X,'THE PARAMETER NE2 AND NP2 MUST MATCH NE AND NP ',
     &        /,1X,'USER MUST CHECK FORT.12 INPUT FILE ',
     &        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
         CALL EXIT(1)
      ENDIF
C
C...  READ IN STARTDRY CODE VALUES
      DO I=1,NP
         READ(12,*) JKI,DUM1,DUM2,STARTDRY(JKI)
         IF(JKI.NE.I) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99805)
            WRITE(16,99805)
99805       FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &           'INPUT ERROR  !!!!!!!!!',
     &           //,1X,'YOUR NODE NUMBERING IS NOT SEQUENTIAL ',
     &           'CHECK YOUR UNIT 12 INPUT FILE CAREFULLY',//)
         ENDIF
      END DO
C
C...  CLOSE UNIT 12 FILE
      CLOSE(12)
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE ReadLegacyStartDryFile
C     ----------------------------------------------------------------


C     ----------------------------------------------------------------
C                       S U B R O U T I N E
C       R E A D  L E G A C Y  B O T T O M  F R I C T I O N  F I L E
C     ----------------------------------------------------------------
C
C     jgf46.00 Subroutine to load up the legacy Spatially Varying
C     Friction Coefficient File (unit 21). This is just a cut-and-paste
C     from the section of the READ_INPUT subroutine that did the same
C     thing. This subroutine is never called. It is vestigial and listed
C     here purely as reference material.
C
C     ----------------------------------------------------------------
      SUBROUTINE ReadLegacyBottomFrictionFile(NP, NScreen, ScreenUnit,
     &     MyProc, NAbOut)
      IMPLICIT NONE
      INTEGER, intent(in) :: NP ! number of nodes in grid file
      INTEGER, intent(in) :: NScreen ! nonzero for debug info to screen
      INTEGER, intent(in) :: ScreenUnit ! i/o for debug info to screen
      INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
      INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16

      CHARACTER(len=80) AFRIC  ! user's comment/description line
      INTEGER NHG    ! node number from file

      OPEN(21,FILE=TRIM(INPUTDIR)//'/'//'fort.21')
      READ(21,'(A80)') AFRIC
      DO I=1,NP
         READ(21,*) NHG,FRIC(NHG)
         IF(NHG.NE.I) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99803)
            WRITE(16,99803)
99803       FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ',
     &           'INPUT ERROR  !!!!!!!!!',//,1X,
     &           'YOUR NODAL FRICTION NUMBERING IS NOT SEQUENTIAL ',
     &           /,1X,'CHECK YOUR UNIT 21 INPUT FILE CAREFULLY',//,1X,
     &           '!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
            CALL EXIT(1)
         ENDIF
      END DO
      WRITE(16,3601) AFRIC
 3601 FORMAT(/,5X,'FRICTION FILE IDENTIFICATN : ',A80,/)
      IF(NABOUT.NE.1) THEN
         WRITE(16,2080)
 2080    FORMAT(/,10X,'NODE',5X,'BOTTOM FRICTION FRIC',5X,/)
         DO I=1,NP
            WRITE(16,2087) I,FRIC(I)
 2087       FORMAT(7X,I6,6X,E17.10)
         END DO
      ELSE
         WRITE(16,3504)
 3504    FORMAT(/,5X,'NODAL BOTTOM FRICTION VALUES ARE AVAILABLE',
     &        /,6X,' IN UNIT 21 INPUT FILE')
      ENDIF
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE ReadLegacyBottomFrictionFile
C     ----------------------------------------------------------------


C     ----------------------------------------------------------------
C        S U B R O U T I N E   N A _ T E R M I N A T E
C     ----------------------------------------------------------------
C     Provide clean termination when an error occurs.
C     ----------------------------------------------------------------
      SUBROUTINE na_terminate(NO_MPI_FINALIZE)
#ifdef CMPI
      USE MESSENGER
#endif
      USE GLOBAL, ONLY : setMessageSource, unsetMessageSource, allMessage,
     &   DEBUG, ECHO, INFO, WARNING, ERROR  
      IMPLICIT NONE 
      LOGICAL, OPTIONAL :: NO_MPI_FINALIZE
C
      call setMessageSource("terminate")
#if defined(NODALATTR_TRACE) || defined(ALL_TRACE)
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
      CALL EXIT(1) 
C
#if defined(NODALATTR_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      call unsetMessageSource()
C     ----------------------------------------------------------------
      END SUBROUTINE na_terminate
C     ----------------------------------------------------------------


C     ----------------------------------------------------------------
C        S U B R O U T I N E   P R E P  C O N D E N S E D   N O D E S
C     ----------------------------------------------------------------
C     Parepare list of condensed nodes
C     ----------------------------------------------------------------
      SUBROUTINE PrepCondensedNodes()
C     ----------------------------------------------------------------
      USE SIZES,ONLY : SZ,MNVEL
      USE MESH,ONLY : NP,LBArray_Pointer,X,Y
      USE BOUNDARIES,ONLY : CSIICN,SIIICN
      IMPLICIT NONE
      INTEGER :: I,J,K,II,I1,I2,J1,J2
      INTEGER :: JMAX
      REAL(SZ) :: X1,Y1,X2,Y2,DX,DY,LEN
      LOGICAL :: ON_THE_SAME_BOUNDARY

      IF (.NOT.LoadCondensedNodes) RETURN

      NListCondensedNodes = 0
      JMAX = 0
      DO I=1,NP
         IF(CondensedNodes(I,1).NE.0) THEN
            NListCondensedNodes = NListCondensedNodes + 1
            DO J=1,CondensedNodesNoOfVals
               IF(CondensedNodes(I,J).EQ.0) EXIT
               IF(JMAX.LT.J) JMAX = J
            ENDDO
         ENDIF
      ENDDO

      ALLOCATE(NCondensedNodes(NumOfNodes))
      ALLOCATE(ListCondensedNodes(NListCondensedNodes,JMAX+1))
      ALLOCATE(NNodesListCondensedNodes(NListCondensedNodes))

      NCondensedNodes(:) = 0

      K = 0
      DO I=1,NP
         IF(CondensedNodes(I,1).GT.0) THEN
            K = K + 1
            ListCondensedNodes(K,1) = I
            DO J=1,CondensedNodesNoOfVals
               IF(CondensedNodes(I,J).LE.0) EXIT
               II = CondensedNodes(I,J)
               NCondensedNodes(I) = J+1
               ListCondensedNodes(K,J+1) = II
               NNodesListCondensedNodes(K) = J+1
               CondensedNodes(II,1) = -I ! To tell that node II is one of condensed nodes
            ENDDO
            DO J=1,CondensedNodesNoOfVals
               IF(CondensedNodes(I,J).LE.0) EXIT
               II = CondensedNodes(I,J)
               NCondensedNodes(II) = NCondensedNodes(I)
            ENDDO
         ENDIF
      ENDDO

      ALLOCATE(CSIICN(MNVEL),SIIICN(MNVEL))
      CSIICN = 0.D0
      SIIICN = 0.D0

      DO I=1,NP
         IF((CondensedNodes(I,1) > 0).AND.
     &      (NCondensedNodes(I).EQ.2)) THEN
            I1 = I
            I2 = CondensedNodes(I,1)

            J1 = LBArray_Pointer(I1)
            J2 = LBArray_Pointer(I2)
            IF(J1<1.OR.J2<1) CYCLE
            IF(J1.EQ.(J2+1).OR.
     &         J1.EQ.(J2-1)) THEN ! Whether Nodes J1 and J2 are on the same boundary.
               ON_THE_SAME_BOUNDARY = .TRUE.
            ELSE
               ON_THE_SAME_BOUNDARY = .FALSE.
            END IF
            IF(ON_THE_SAME_BOUNDARY) CYCLE

            X1 = X(I1)
            Y1 = Y(I1)
            X2 = X(I2)
            Y2 = Y(I2)
            DX = X1 - X2
            DY = Y1 - Y2
            LEN = SQRT(DX*DX + DY*DY)
            DX = DX/LEN
            DY = DY/LEN

            CSIICN(J1) = DX
            SIIICN(J1) = DY

            CSIICN(J2) = DX
            SIIICN(J2) = DY
         END IF
      END DO

      END SUBROUTINE PrepCondensedNodes
C     ----------------------------------------------------------------

C-----------------------------------------------------------------------
      END MODULE NodalAttributes
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------


