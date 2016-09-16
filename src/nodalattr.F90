!******************************************************************************
!     PADCIRC VERSION 46.00 xx/xx/2006
!     last changes in this file VERSION 46.00

!     Written for ADCIRC v46.00 by Jason G. Fleming.
!******************************************************************************

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     M O D U L E   N O D A L  A T T R I B U T E S
!-----------------------------------------------------------------------

!     jgf46.00 This module manages nodal attribute data, including
!     bottom friction, tau0, startdry, directional wind speed reduction,
!     and etc. Will read the Nodal Attributes File (unit 13) and
!     initialize the nodal attribute arrays.

!     Handling data by label rather than an integer encoding should
!     result in increased transparency as well as ease the transition to
!     HDF5/NetCDF i/o. The labels were chosen according to the
!     guidelines of the CF Standard. Creating labels according to CF
!     Standard Guidelines should enhance interoperability with other
!     simulation frameworks.

!     To use a nodal attribute contained in the fort.13 file, the
!     corresponding attribute name must appear in the fort.15 file.  A
!     list of nodal attributes is read in from the fort.15 file if the
!     fort.15 parameter NWP > 0. This also signals ADCIRC to look for a
!     fort.13 file.

!     Summary of the file format for the Nodal Attributes File:

!     AGRID                       ! user's comment line - should be a cross
!                                 !    reference to the grid file
!     NumOfNodes                  ! number of nodes, must match NP
!                                 !    from grid file
!     NAttr                       ! number of attributes contained in this file

!     do i=1, NAttr
!        AttrName(i)              ! nodal attribute name (see
!                                 !    valid names below)
!        Units(i)                 ! physical units (ft, m/s, none)
!        ValuesPerNode(i)         ! number of values at each node for
!                                 !   a particular attribute
!        DefaultVal(i)            ! default value(s) for the nodal attribute
!     end do

!     do i=1, NAttr
!        AttrName(i)              ! label of the attribute, again
!        NumNodesNotDefaultVal(i) ! number of nodes with non-default values
!        do j=1, NumNodesNotDefault(i)
!           n, (Attr(n,k), k=1,ValuesPerNode(i))
!        end do
!     end do



!     Valid labels are as follows:

!     ADCIRC Variable:       CF-Style Label:
!      Tau0                  "primitive_weighting_in_continuity_equation"
!      StartDry              "surface_submergence_state"
!      Fric                  "quadratic_friction_coefficient_at_sea_floor"
!      z0Land                "surface_directional_effective_roughness_length"
!                            (z0Land has ValuesPerNode = 12)
!      VCanopy               "surface_canopy_coefficient"
!      BK,BAlpha,BDelX,POAN  "bridge_pilings_friction_parameters"
!                            (bridge_pilings... has ValuesPerNode=4)
!      ManningsN             "mannings_n_at_sea_floor"
!      Chezy                 "chezy_friction_coefficient_at_sea_floor"
!      Z0b_var               "bottom_roughness_length"
!      GeoidOffset           "sea_surface_height_above_geoid"
!      EVM        "average_horizontal_eddy_viscosity_in_sea_water_wrt_depth"
!      EVC        "average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth"
!      Tao0MinMax "min_and_max_primitive_weighting_in_continuity_equation"
!      River_et_WSE "initial_river_elevation"

!-----------------------------------------------------------------------
    MODULE NodalAttributes
    USE SIZES
    USE GLOBAL, ONLY : scratchMessage, allMessage, logMessage, &
    DEBUG, ECHO, INFO, WARNING, ERROR

! Casey 100210: Allow SWAN to handle wave refraction as a nodal attribute.
!             I've placed these changes outside the #ifdef CSWAN flags
!             because we want to be able to use the same fort.13 files
!             for both ADCIRC and SWAN+ADCIRC runs.  This way, the new
!             nodal attribute will be processed but only applied when
!             ADCIRC is coupled to SWAN.
    LOGICAL              :: LoadSwanWaveRefrac
    LOGICAL              :: FoundSwanWaveRefrac
    CHARACTER(LEN=80)    :: SwanWaveRefracUnits
    INTEGER :: SwanWaveRefracNoOfVals
    REAL(SZ)             :: SwanWaveRefracDefVal
    REAL(SZ),ALLOCATABLE :: SwanWaveRefrac(:)

! Corbitt 120321: Allow Advection to be Turned on Locally instead of Globally
    LOGICAL              :: LoadAdvectionState
    LOGICAL              :: FoundAdvectionState
    CHARACTER(LEN=80)    :: AdvectionStateUnits
    INTEGER :: AdvectionStateNoOfVals
    REAL(SZ)             :: AdvectionStateDefVal
    REAL(SZ),ALLOCATABLE :: AdvectionState(:)

!     The following flags are .true. if the corresponding data are
!     required for the run, according to the unit 15 control file
    LOGICAL :: LoadTau0
    LOGICAL :: LoadStartDry
    LOGICAL :: LoadDirEffRLen
    LOGICAL :: LoadCanopyCoef
    LOGICAL :: LoadQuadraticFric
    LOGICAL :: LoadBridgePilings
    LOGICAL :: LoadChezy
    LOGICAL :: LoadManningsN
    LOGICAL :: LoadGeoidOffset
    LOGICAL :: LoadEVM
    LOGICAL :: LoadEVC
    LOGICAL :: LoadTau0MinMax
    LOGICAL :: LoadZ0b_var
    LOGICAL LoadEleSlopeLim ! zc: elemental slope limiter
    LOGICAL :: LoadRiver_et_WSE  ! tcm 20140502 v51.27

!     The following flags are .true. if there are data with the
!     corresponding label in the unit 13 file.
    LOGICAL :: FoundTau0
    LOGICAL :: FoundStartDry
    LOGICAL :: FoundDirEffRLen
    LOGICAL :: FoundCanopyCoef
    LOGICAL :: FoundQuadraticFric
    LOGICAL :: FoundBridgePilings
    LOGICAL :: FoundChezy
    LOGICAL :: FoundManningsN
    LOGICAL :: FoundGeoidOffset
    LOGICAL :: FoundEVM
    LOGICAL :: FoundEVC
    LOGICAL :: FoundTau0MinMax
    LOGICAL :: FoundZ0b_var
    LOGICAL :: FoundEleSlopeLim
    LOGICAL :: FoundRiver_et_WSE ! tcm 20140502 v51.27

!     These variables hold the strings which describe the attribute's
!     units. These data are loaded from the file, but not used as of
!     v46.00.
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

!     These variables hold the number of values per node for each
!     attribute.
    INTEGER :: Tau0NoOfVals
    INTEGER :: StartDryNoOfVals
    INTEGER :: DirEffRLenNoOfVals
    INTEGER :: CanopyCoefNoOfVals
    INTEGER :: QuadraticFricNoOfVals
    INTEGER :: BridgePilingsNoOfVals
    INTEGER :: ChezyNoOfVals
    INTEGER :: ManningsNNoOfVals
    INTEGER :: GeoidOffsetNoOfVals
    INTEGER :: EVMNoOfVals
    INTEGER :: EVCNoOfVals
    INTEGER :: Tau0MinMaxNoOfVals
    INTEGER :: Z0b_varNoOfVals
    INTEGER :: EleSlopeLimNoofVals
    INTEGER :: River_et_WSENoOfVals  ! tcm 20140502 v51.27

!     These variables hold the default values for each attribute.
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

    INTEGER :: NumOfNodes    ! number of nodes listed in unit 13 file, cf. NP
    INTEGER :: NAttr         ! number of nodal attributes in the unit 13 file

!     The following variables are inputs from the unit 15 model param. file
    INTEGER :: NWP     ! number of nodal attributes to read from file
    INTEGER :: NoLiBF  ! nonlinear bottom friction indicator
    REAL(SZ) Tau0   ! primitive continuity eqn. weight
    REAL(SZ) Tau    ! linear friction coefficient (1/sec)
    REAL(SZ) CF     ! 2DDI bottom fric. coef., effect varies based on NoLiBF
    REAL(SZ) HBreak ! break depth for NOLIBF == 2
    REAL(SZ) FTheta ! dimless param. for NOLIBF == 2
    REAL(SZ) FGamma ! dimless param. for NOLIBF == 
    REAL(SZ) ESLM   ! horizontal eddy viscosity (length^2/time)
    REAL(SZ) ESLC   ! horizontal eddy diffusivity (length^2/time)
    INTEGER ::  IFLINBF! flag to turn on linear bottom friction
    INTEGER ::  IFNLBF ! flag to turn on nonlinear bottom friction
    INTEGER ::  IFHYBF ! flag to turn on hybrid bottom friction
    REAL(SZ) BFCdLLimit ! lower limit of quadratic bottom friction
    REAL(SZ) Tau0FullDomainMin ! lower limit of tau0 if time varying
    REAL(SZ) Tau0FullDomainMax ! upper limit of tau0 if time varying
!     jgf48.42 Used to control the back-loaded time averaged tau0
    REAL(SZ), PARAMETER :: AlphaTau0 = 0.25d0

!     Nodal attributes.
    REAL(SZ), ALLOCATABLE :: STARTDRY(:) ! 1=nodes below geoid initially dry
    REAL(SZ), ALLOCATABLE :: FRIC(:)     ! bottom friction coefficient
    REAL(SZ), ALLOCATABLE, TARGET :: TAU0VAR(:)! primitive equation weighting
    REAL(SZ), ALLOCATABLE :: Tau0Temp(:) ! used in time varying tau0
!     jjw&sb46.39.sb01: Base (original) primitive equation weighting.
!     Tau0Var may be optimized later based on Tau0Base
    REAL(SZ), ALLOCATABLE :: TAU0BASE(:)
!     jgf47.33 Used for time averaged tau0
    REAL(SZ), ALLOCATABLE :: LastTau0(:)

!     jgf47.06: Added variables to trigger calculations and output of tau0
!     jgf47.30: Added "FullDomain" or "High Res Areas Only" distinction
!     jgf47.31: Added time averaging to time varying tau0
!     jgf48.42: Added backloaded time averaged tau0
    LOGICAL :: HighResTimeVaryingTau0    ! .TRUE. if Tau0 == -3.x in fort.15
    LOGICAL :: FullDomainTimeVaryingTau0 ! .TRUE. if Tau0 == -5.x in fort.15
    LOGICAL :: OutputTau0      ! .TRUE. if Tau0Dig2 == -1 in fort.15
    LOGICAL :: TimeAveragedTau0 ! .TRUE. if Tau0 == -6.x in fort.15
    LOGICAL :: BackLoadedTimeAveragedTau0 ! .TRUE. if Tau0 == -7.x
    LOGICAL,ALLOCATABLE :: elemental_slope_limiter_active(:) ! .TRUE. if elemental_slope_limiter_grad_max
! has been exceeded, initialized as .false.
    LOGICAL,ALLOCATABLE :: elemental_slope_limiter_max_exceeded(:) ! .TRUE. if maximum gradient has been
! exceeded, used for gradient warnings

    REAL(SZ), ALLOCATABLE :: z0land(:,:) ! directional wind speed red. fac.
    REAL(SZ), ALLOCATABLE :: vcanopy(:)  ! canopy coefficient
!     The following attribute contains BK(I),BALPHA(I),BDELX(I), and POAN(I)
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

    INTEGER :: i        ! node loop counter
    INTEGER :: j        ! attribute values loop counter
    INTEGER :: k        ! attribute loop counter

!  X D M F

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

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    CONTAINS !- - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ----------------------------------------------------------------
!     S U B R O U T I N E     I N I T  N  A  M O D U L E
!     ----------------------------------------------------------------

!     jgf46.00 Subroutine to initialize the variables in the nodal
!     attributes module.

!     ----------------------------------------------------------------
    SUBROUTINE InitNAModule()
    IMPLICIT NONE

! Casey 100210: Make changes compact.
    LoadSwanWaveRefrac     = .FALSE. 
    FoundSwanWaveRefrac    = .FALSE. 
    SwanWaveRefracNoOfVals = 1
    SwanWaveRefracDefVal   = 0.D0
! Corbitt 120321:
    LoadAdvectionState     = .FALSE. 
    FoundAdvectionState    = .FALSE. 
    AdvectionStateNoOfVals = 1
    AdvectionStateDefVal   = 0.D0

    LoadTau0           = .FALSE. 
    LoadStartDry       = .FALSE. 
    LoadDirEffRLen     = .FALSE. 
    LoadManningsN      = .FALSE. 
    LoadQuadraticFric  = .FALSE. 
    LoadChezy          = .FALSE. 
    LoadBridgePilings  = .FALSE. 
    LoadCanopyCoef     = .FALSE. 
    LoadGeoidOffset    = .FALSE. 
    LoadEVM            = .FALSE. 
    LoadEVC            = .FALSE. 
    LoadTau0MinMax     = .FALSE. 
    LoadZ0b_var        = .FALSE. 
    LoadEleSlopeLim    = .FALSE. 
    LoadRiver_et_WSE     = .FALSE. 

    FoundTau0           = .FALSE. 
    FoundStartDry       = .FALSE. 
    FoundDirEffRLen     = .FALSE. 
    FoundManningsN      = .FALSE. 
    FoundQuadraticFric  = .FALSE. 
    FoundChezy          = .FALSE. 
    FoundBridgePilings  = .FALSE. 
    FoundCanopyCoef     = .FALSE. 
    FoundGeoidOffset    = .FALSE. 
    FoundEVM            = .FALSE. 
    FoundEVC            = .FALSE. 
    FoundTau0MinMax     = .FALSE. 
    FoundZ0b_var        = .FALSE. 
    FoundEleSlopeLim    = .FALSE. 
    FoundRiver_et_WSE     = .FALSE. 


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
    River_et_WSENoOfVals    = 1

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
    River_et_WSEDefVal       = 0.d0

    HighResTimeVaryingTau0     = .FALSE. 
    FullDomainTimeVaryingTau0  = .FALSE. 
    OutputTau0                 = .FALSE. 
    TimeAveragedTau0           = .FALSE. 
    BackLoadedTimeAveragedTau0 = .FALSE. 

    HBREAK=1.d0
    FTHETA=1.d0
    FGAMMA=1.d0

!   kmd48.33bc this resets the ESLM to 0 and if using constant eddy
!              viscosity it eliminates what the user specified in
!              the input file. The InitNAModule originally come before
!              the read_input call but now it appears after the read_input
!              call.
!      ESLM=0.0
    ESLC=0.0

    RETURN
!     ----------------------------------------------------------------
    END SUBROUTINE InitNAModule
!     ----------------------------------------------------------------

!     ----------------------------------------------------------------
!     S U B R O U T I N E    R E A D   N O D A L   A T T R   X D M F
!     ----------------------------------------------------------------
!     jgf51.21.24: Read the nodal attributes in XDMF format.
!     ----------------------------------------------------------------
    SUBROUTINE readNodalAttrXDMF()
    USE MESH, ONLY : replaceNullsWithSpaces
    IMPLICIT NONE
#ifndef ADCXDMF
    call allMessage(ERROR, &
    'An XDMF formatted nodal attributes (fort.13) file was specified.')
    call allMessage(ERROR, &
    'This ADCIRC executable was not compiled with XDMF support.')
    call na_terminate()
#else
    include 'adcirc_Xdmf.f90'
    integer*8 :: xdmfFortranObj ! object that receives the data

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

    startIndex = 0
    arrayStride = 1
    valueStride = 1
    openMaps = 1
    openAttributes = 1
    openInformations = 1
    openSets = 1

    NAFound = .FALSE. 

!     Determine if the Nodal Attributes File exists.
    INQUIRE(FILE=TRIM(naFileName),EXIST=NAFound)

    IF ( .NOT. NAFound) THEN
        write(scratchMessage,'(a)') 'The XDMF nodal attributes file' // &
        trim(naFileName)//'was not found.'
        call allMessage(ERROR, scratchMessage)
        call na_terminate()
    ENDIF

    write(6,'(a)') 'INFO: Reading data from the "' &
    // trim(adjustl(naFileName)) // '" XDMF file.'
    call xdmfinit(xdmfFortranObj)
    call xdmfRead(xdmfFortranObj, trim(adjustl(naFileName))//char(0))
    call xdmfOpenDomainGrid(xdmfFortranObj, XDMF_GRID_TYPE_UNSTRUCTURED, &
    gridIndex, openMaps, openAttributes, openInformations, openSets)

    call xdmfRetrieveNumAttributes(xdmfFortranObj, numAttributes)
    nAttr = numAttributes - 1 ! the depth is included so don't count it
          
    allocate(na(nAttr))
    write(6,'("INFO: Grid ",i0," contains ",i0," nodal attributes.")') gridIndex, nAttr

! populate the names of the nodal attributes
    nattrCount = 1
    do attributeIndex=0, numAttributes - 1
        call xdmfRetrieveAttributeName(xdmfFortranObj, attributeIndex, &
        itemName, nameLength)
        call replaceNullsWithSpaces(itemName)
        if (trim(itemName) == 'depth') then
            cycle
        endif
        write(6,'("INFO: Grid ",i0," Attribute ",i0," is named ",a)') gridIndex, attributeIndex, trim(itemName)
        na(nattrCount)%attrName = trim(itemName)
        SELECT CASE (trim(na(nattrCount)%AttrName))
        case("depth")
        cycle
        CASE("primitive_weighting_in_continuity_equation")
        FoundTau0 = .TRUE. 
        CASE("surface_submergence_state")
        FoundStartDry = .TRUE. 
        CASE("quadratic_friction_coefficient_at_sea_floor")
        FoundQuadraticFric = .TRUE. 
        CASE("surface_directional_effective_roughness_length")
        FoundDirEffRLen = .TRUE. 
        CASE("surface_canopy_coefficient")
        FoundCanopyCoef = .TRUE. 
        CASE("bridge_pilings_friction_parameters")
        FoundBridgePilings = .TRUE. 
        CASE("mannings_n_at_sea_floor")
        FoundManningsN = .TRUE. 
        CASE("bottom_roughness_length")
        FoundZ0b_var = .TRUE. 
        CASE("chezy_friction_coefficient_at_sea_floor")
        FoundChezy = .TRUE. 
        CASE("sea_surface_height_above_geoid")
        FoundGeoidOffset = .TRUE. 
        CASE &
        ("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")
        FoundEVM = .TRUE. 
        CASE &
        ("average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth")
        FoundEVC = .TRUE. 
        CASE &
        ("min_and_max_primitive_weighting_in_continuity_equation")
        FoundTau0MinMax = .TRUE. 
    ! Casey 100210: Allow SWAN to handle wave refraction as a nodal attribute.
        CASE("wave_refraction_in_swan")
        FoundSwanWaveRefrac = .TRUE. 
    ! Corbitt 120321: Allow advection to be turned on locally instead of globally
        CASE("advection_state")
        FoundAdvectionState = .TRUE. 
        CASE("elemental_slope_limiter")
        FoundEleSlopeLim = .TRUE. 
        CASE("initial_river_elevation")
        FoundRiver_et_WSE     = .TRUE. 
        CASE DEFAULT
        scratchMessage = "Unrecognized nodal attribute : " // &
        trim(na(nattrCount) % attrName)
        call allMessage(WARNING, scratchMessage)
        END SELECT
        nattrCount = nattrCount + 1
    END DO

!     Determine if there are any attributes required by the fort.15 file
!     that are not in the nodal attributes file.
    call checkForMissingNodalAttributes()

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
            if (trim(itemKey) == trim(na(i)%attrName) // ' number_of_values') then
                read(itemValue,*) na(i)%numVals
                allocate(na(i)%defaultVals(na(i)%numVals))
                exit
            endif
        end do
        end select
    end do

!     Allocate memory for nodal attributes now that we have the
!     numOfNodes value
    call allocateNodalAttributes()

! populate the units and the default values
    call xdmfRetrieveNumInformation(xdmfFortranObj, numInformation)
    do informationIndex=0,numInformation-1
        call xdmfRetrieveInformation(xdmfFortranObj, informationIndex, itemKey, keyLength, itemValue, valueLength)
        call replaceNullsWithSpaces(itemKey)
        call replaceNullsWithSpaces(itemValue)
        do i=1,nAttr
            if (trim(itemKey) == trim(na(i)%attrName)// ' units') then
                na(i) % units = trim(itemValue)
            endif
            if (trim(itemKey) == trim(na(i)%attrName) // ' default_values') then
                write(6,*) 'default_values'//trim(itemValue)
                read(itemValue,*) (na(i)%defaultVals(j),j=1,na(i)%numVals)
            endif
        end do
    end do

! populate the data
    do attributeIndex=0, numAttributes - 1
        call xdmfRetrieveAttributeName(xdmfFortranObj, attributeIndex, &
        itemName, nameLength)
        call replaceNullsWithSpaces(itemName)
        if (trim(itemName) == 'depth') then
            cycle
        endif
        do i=1,nAttr
            if (trim(itemName) == trim(na(i)%attrName)) then
                write(6,'(a)') 'loading nodal attribute data for '//trim(itemName)
                if (na(i)%numVals == 1) then
                    allocate(na(i)%xdmfArray(numOfNodes))
                    attributeType = XDMF_ATTRIBUTE_TYPE_SCALAR
                    numValues = numOfNodes * na(i) % numVals
                    call xdmfRetrieveAttributeValues(xdmfFortranObj, attributeIndex, &
                    na(i)%xdmfArray, XDMF_ARRAY_TYPE_FLOAT64, numValues, &
                    startIndex, arrayStride, valueStride)
                else
                    attributeType = XDMF_ATTRIBUTE_TYPE_MATRIX
                    numValues = numOfNodes * na(i) % numVals
                    allocate(na(i)%xdmfMatrix(na(i)%numVals,numOfNodes))
                    call xdmfRetrieveAttributeValues(xdmfFortranObj, attributeIndex, &
                    na(i)%xdmfMatrix, XDMF_ARRAY_TYPE_FLOAT64, numValues, &
                    startIndex, arrayStride, valueStride)
                endif
                exit
            endif
        end do
    end do

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
        CASE &
        ("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")
        IF (LoadEVM) THEN
            evm = na(i)%xdmfArray
        ENDIF
        CASE &
        ("average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth")
        IF (LoadEVC) THEN
            evc  = na(i)%xdmfArray
        ENDIF
        CASE &
        ("min_and_max_primitive_weighting_in_continuity_equation")
        IF (LoadTau0MinMax) THEN
            tau0minmax  = na(i)%xdmfMatrix
        ENDIF
    ! Casey 100210: Allow SWAN to handle wave refraction as a nodal attribute.
        CASE("wave_refraction_in_swan")
        IF (LoadSwanWaveRefrac) THEN
            swanWaveRefrac = na(i)%xdmfArray
        ENDIF
    ! Corbitt 120321: Allow Advection to be handled locally instead of globally.
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
        CASE DEFAULT
        scratchMessage = "Unrecognized nodal attribute : " // &
        trim(na(i) % attrName)
        call allMessage(WARNING, scratchMessage)
        END SELECT
    END DO
#endif
!-----------------------------------------------------------------------
    end subroutine readNodalAttrXDMF
!-----------------------------------------------------------------------

!     ----------------------------------------------------------------
!                       S U B R O U T I N E
!    C H E C K   F O R   M I S S I N G   N O D A L   A T T R I B U T E S
!     ----------------------------------------------------------------
!     jgf51.21.24: If a nodal attribute was specified in the control
!     (fort.15) file, but is not present in the nodal attributes (fort.13)
!     file, then we write out an error message and quit.
!     ----------------------------------------------------------------
    subroutine checkForMissingNodalAttributes()
    implicit none

!     Determine if there are any attributes required by the fort.15 file
!     that are not in the nodal attributes file.
    IF(((LoadTau0) .AND. ( .NOT. FoundTau0)) .OR. &
    ((LoadStartDry) .AND. ( .NOT. FoundStartDry)) .OR. &
    ((LoadQuadraticFric) .AND. &
    ( .NOT. FoundQuadraticFric)) .OR. &
    ((LoadDirEffRLen) .AND. &
    ( .NOT. FoundDirEffRLen)) .OR. &
    ((LoadCanopyCoef) .AND. &
    ( .NOT. FoundCanopyCoef)) .OR. &
    ((LoadBridgePilings) .AND. &
    ( .NOT. FoundBridgePilings)) .OR. &
    ((LoadManningsN) .AND. &
    ( .NOT. FoundManningsN)) .OR. &
    ((LoadZ0b_var) .AND. &
    ( .NOT. FoundZ0b_var)) .OR. &
    ((LoadGeoidOffset) .AND. &
    ( .NOT. FoundGeoidOffset)) .OR. &
    ((LoadChezy) .AND. ( .NOT. FoundChezy)) .OR. &
    ((LoadEVM) .AND. ( .NOT. FoundEVM)) .OR. &
    ((LoadEVC) .AND. ( .NOT. FoundEVC)) .OR. &
! Casey 100210: Allow SWAN to handle wave refraction as a nodal attribute.
    ((LoadSwanWaveRefrac) .AND. ( .NOT. FoundSwanWaveRefrac)) .OR. &
! Corbitt 120321: Allow advection to be turned on locally instead of globally
    ((LoadAdvectionState) .AND. ( .NOT. FoundAdvectionState)) .OR. &
    ((LoadEleSlopeLim) .AND. ( .NOT. FoundEleSlopeLim)) .OR. &
    ((LoadRiver_et_WSE) .AND. ( .NOT. FoundRiver_et_WSE)) .OR. &
    ((LoadTau0MinMax) .AND. ( .NOT. FoundTau0MinMax))) THEN
        WRITE(scratchMessage,1111)
        1111 FORMAT('Nodal Attributes file (unit 13) does ' &
        'not contain all the attributes listed in the ' &
        'model parameter file (unit 15).')
        call allMessage(ERROR, scratchMessage)
        call na_terminate()
    ENDIF
!     ----------------------------------------------------------------
    end subroutine checkForMissingNodalAttributes
!     ----------------------------------------------------------------

!     ----------------------------------------------------------------
!                     S U B R O U T I N E
!          A L L O C A T E   N O D A L   A T T R I B U T E S
!     ----------------------------------------------------------------
!     Allocate memory to hold our data.
    subroutine allocateNodalAttributes()
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
! Casey 100210: Allow SWAN to handle wave refraction as a nodal attribute.
    ALLOCATE(SwanWaveRefrac(NumOfNodes))
! Corbitt 120321: Allow advection to be turned on locally instead of globally
    ALLOCATE(AdvectionState(NumOfNodes))
    ALLOCATE(River_et_WSE(NumOfNodes))  ! tcm 20140502 v51.27
    ALLOCATE(elemental_slope_limiter_grad_max(NumOfNodes))
    ALLOCATE(elemental_slope_limiter_active(NumOfNodes))
    ALLOCATE(elemental_slope_limiter_max_exceeded(NumOfNodes))
    elemental_slope_limiter_active(:) = .FALSE.  !...Alloate and initialize to keep
!   track of elemental_slope_limiter nodes
    elemental_slope_limiter_max_exceeded(:) = .FALSE. 
!     ----------------------------------------------------------------
    end subroutine allocateNodalAttributes
!     ----------------------------------------------------------------

!     ----------------------------------------------------------------
!     S U B R O U T I N E     R E A D  N O D A L  A T T R
!     ----------------------------------------------------------------

!     jgf46.00 Subroutine to read the nodal attributes file (unit 13).

!     ----------------------------------------------------------------
    SUBROUTINE ReadNodalAttr(NScreen, ScreenUnit, MyProc, NAbOut)
    USE MESH, ONLY : NP
    IMPLICIT NONE
    INTEGER, intent(in) :: NScreen ! nonzero for debug info to screen
    INTEGER, intent(in) :: ScreenUnit ! i/o for screen
    INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
    INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16

    LOGICAL :: NAFound  ! .TRUE. if Nodal Attributes File (fort.13) exists
    INTEGER :: ErrorIO  ! zero if file opened successfully
    CHARACTER(len=80) AttrName ! string where the attribute name is stored
    CHARACTER(len=80) header   ! string where alphanumeric file id is stored
    INTEGER :: NumNodesNotDefault ! number of individual nodes to specify
    LOGICAL :: SkipDataSet ! .TRUE. if a data set in unit 13 is not needed
    CHARACTER(len=80) Skipped ! data in unit 13 we do not need
    INTEGER :: L                 ! line counter

! temp array; used to load a real from the file,
! then convert to integer
    REAL(sz), ALLOCATABLE :: real_loader(:)

    NAFound = .FALSE. 
    SkipDataSet = .FALSE. 

!     Check to make sure that NWP is a valid number.
    IF (NWP < 0) THEN
        IF(NSCREEN /= 0 .AND. MYPROC == 0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NWP =',NWP
            WRITE(ScreenUnit,9728)
            WRITE(ScreenUnit,9973)
        ENDIF
        WRITE(16,9972)
        WRITE(16,*) 'NWP =',NWP
        WRITE(16,9728)
        WRITE(16,9973)
        9728 FORMAT(/,1X,'Your selection of NWP (a UNIT 15 input ', &
        'parameter) is not an allowable value')
        STOP            ! We're toast.
    ENDIF

!     Check to see if there are nodal attributes to be read in. If not,
!     simply return.
    IF (NWP == 0) THEN
        WRITE(16,231) NWP
        231 FORMAT(/,5X,'NWP = ',I2, &
        /,9X,'A Nodal Attributes File (unit 13)', &
        /,9X,'will not be used.')
        RETURN
    ENDIF

!     Read the unit 15 control file to determine what data must be
!     loaded from nodal attributes file.
    WRITE(16,235) NWP
    235 FORMAT(/,9X,'Need to load ',I2,' nodal attribute(s):')
    DO k=1,NWP
        READ(15,*) AttrName
        WRITE(16,'(14X,A80)') AttrName
        SELECT CASE (AttrName)
        CASE("primitive_weighting_in_continuity_equation")
        LoadTau0 = .TRUE. 
        CASE("surface_submergence_state")
        LoadStartDry = .TRUE. 
        CASE("quadratic_friction_coefficient_at_sea_floor")
        LoadQuadraticFric = .TRUE. 
        CASE("surface_directional_effective_roughness_length")
        LoadDirEffRLen = .TRUE. 
        CASE("surface_canopy_coefficient")
        LoadCanopyCoef = .TRUE. 
        CASE("bridge_pilings_friction_parameters")
        LoadBridgePilings = .TRUE. 
        CASE("mannings_n_at_sea_floor")
        LoadManningsN = .TRUE. 
        CASE("chezy_friction_coefficient_at_sea_floor")
        LoadChezy = .TRUE. 
        CASE("bottom_roughness_length")
        LoadZ0b_var = .TRUE. 
        CASE("sea_surface_height_above_geoid")
        LoadGeoidOffset = .TRUE. 
        CASE &
        ("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")
        LoadEVM = .TRUE. 
        CASE &
        ("average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth")
        LoadEVC = .TRUE. 
        CASE &
        ("min_and_max_primitive_weighting_in_continuity_equation")
        LoadTau0MinMax = .TRUE. 
    ! Casey 100210: Allow SWAN to handle wave refraction as a nodal attribute.
        CASE("wave_refraction_in_swan")
        LoadSwanWaveRefrac = .TRUE. 
    ! Corbitt 120321: Allow advection to be turned on locally instead of globally
        CASE("advection_state")
        LoadAdvectionState = .TRUE. 
        CASE("elemental_slope_limiter")
        LoadEleSlopeLim = .TRUE. 
        CASE("initial_river_elevation")  ! tcm 20140502 v51.27
        LoadRiver_et_WSE = .TRUE. 
        CASE DEFAULT
        WRITE(16,1000)          ! unit 15 Model Parameter file
        WRITE(16,1021) AttrName ! contains invalid name
        IF (NScreen /= 0 .AND. MyProc == 0) THEN
            WRITE(ScreenUnit,1000)
            WRITE(ScreenUnit,1021) AttrName
        ENDIF
        END SELECT
    ENDDO

    WRITE(16,232) NWP
    232 FORMAT(/,5X,'NWP = ',I2, &
    /,9X,'Must read Nodal Attributes File (unit 13).')

!  R E A D   N O D A L   A T T R I B U T E S   X D M F

    if (naType == XDMF) then
        call readNodalAttrXDMF()
        return
    endif

!  R E A D   N O D A L   A T T R I B U T E S   A S C I I

!     Determine if the Nodal Attributes File exists.
    INQUIRE(FILE=TRIM(INPUTDIR)//'/'//'fort.13',EXIST=NAFound)

    IF ( .NOT. NAFound) THEN
        WRITE(16,1001)         ! Nodal Attributes file
        WRITE(16,1011)         ! was not found.
        WRITE(16,9973)         ! execution terminated
        IF (NScreen /= 0 .AND. MyProc == 0) THEN
            WRITE(ScreenUnit,1001)
            WRITE(ScreenUnit,1011)
            WRITE(ScreenUnit,9973)      ! execution terminated
        ENDIF
        STOP
    ENDIF

!     Now open the nodal attributes (unit 13) file.
    WRITE(16,240)
    240 FORMAT(/,9X,'Nodal Attributes File (unit 13) was found.', &
    ' Opening file.')
    OPEN(UNIT=13, FILE=TRIM(INPUTDIR)//'/'//'fort.13', &
    IOSTAT=ErrorIO)
    IF ( ErrorIO > 0 ) THEN
        WRITE(16,1001)         ! Nodal attribute file
        WRITE(16,1005)         ! exists but can't be opened
        WRITE(16,9973)         ! execution terminated
        IF (NScreen /= 0 .AND. MyProc == 0) THEN
            WRITE(ScreenUnit,1001)
            WRITE(ScreenUnit,1005)
            WRITE(ScreenUnit,9973)
        ENDIF
        STOP                   ! We're toast.
    ENDIF

!     Read each attribute name, units, number of values, and default value
    READ(13,'(A80)') header
    WRITE(16,250)
    250 FORMAT(/,9X,'User comment line from unit 13:')
    WRITE(16,'(14X,A80,/)') header
    READ(13,*) NumOfNodes     ! number of nodes according to unit 13

!     ERROR CHECK: If a nodal attributes file is being used, check to
!     see that the number of nodes in the nodal attribute file is the
!     same as the number of nodes in the grid file.
!     jgf52.18: Make this check immediately so execution can
!     be stopped immediately if this is not the right nodal attributes
!     file for this mesh.
    if ( (nwp /= 0) .AND. (numOfNodes /= np) ) then
        write(scratchMessage,9901) np, numOfNodes
        call allMessage(ERROR,scratchMessage)
        call na_terminate()
        9901 format('The number of nodes in the mesh file (unit 14) is ',i0, &
        'while the number of nodes listed on the 2nd line of the '// &
        'nodal attributes file (unit 13) is ',i0,'.'// &
        'These numbers must match. '// &
        'Is this the right nodal attributes file for this mesh file?')
    endif

    READ(13,*) NAttr          ! number of attributes in the unit 13 file
    DO k=1, NAttr
        READ(13,*) AttrName
        WRITE(16,'(9X,A80)') AttrName
        WRITE(16,260)
        260 FORMAT(14X,'was found!',/)
        SELECT CASE (AttrName)
        CASE("primitive_weighting_in_continuity_equation")
        FoundTau0 = .TRUE. 
        READ(13,'(A80)') Tau0Units
        READ(13,*) Tau0NoOfVals
        READ(13,*) Tau0DefVal
        CASE("surface_submergence_state")
        FoundStartDry = .TRUE. 
        READ(13,'(A80)') StartDryUnits
        READ(13,*) StartDryNoOfVals
        READ(13,*) StartDryDefVal
        CASE("quadratic_friction_coefficient_at_sea_floor")
        FoundQuadraticFric = .TRUE. 
        READ(13,'(A80)') QuadraticFricUnits
        READ(13,*) QuadraticFricNoOfVals
        READ(13,*) QuadraticFricDefVal
        CASE("surface_directional_effective_roughness_length")
        FoundDirEffRLen = .TRUE. 
        READ(13,'(A80)') DirEffRLenUnits
        READ(13,*) DirEffRLenNoOfVals
        READ(13,*) &
        (DirEffRLenDefVal(j),j=1,DirEffRLenNoOfVals)
        CASE("surface_canopy_coefficient")
        FoundCanopyCoef = .TRUE. 
        READ(13,'(A80)') CanopyCoefUnits
        READ(13,*) CanopyCoefNoOfVals
        READ(13,*) CanopyCoefDefVal
        CASE("bridge_pilings_friction_parameters")
        FoundBridgePilings = .TRUE. 
        READ(13,'(A80)') BridgePilingsUnits
        READ(13,*) BridgePilingsNoOfVals
        READ(13,*) &
        (BridgePilingsDefVal(j),j=1,BridgePilingsNoOfVals)
        CASE("mannings_n_at_sea_floor")
        FoundManningsN = .TRUE. 
        READ(13,'(A80)') ManningsNUnits
        READ(13,*) ManningsNNoOfVals
        READ(13,*) ManningsNDefVal
        CASE("bottom_roughness_length")
        FoundZ0b_var = .TRUE. 
        READ(13,'(A80)') Z0b_varUnits
        READ(13,*) Z0b_varNoOfVals
        READ(13,*) Z0b_varDefVal
        CASE("chezy_friction_coefficient_at_sea_floor")
        FoundChezy = .TRUE. 
        READ(13,'(A80)') ChezyUnits
        READ(13,*) ChezyNoOfVals
        READ(13,*) ChezyDefVal
        CASE("sea_surface_height_above_geoid")
        FoundGeoidOffset = .TRUE. 
        READ(13,'(A80)') GeoidOffsetUnits
        READ(13,*) GeoidOffsetNoOfVals
        READ(13,*) GeoidOffsetDefVal
        CASE &
        ("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")
        FoundEVM = .TRUE. 
        READ(13,'(A80)') EVMUnits
        READ(13,*) EVMNoOfVals
        READ(13,*) EVMDefVal
        CASE &
        ("average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth")
        READ(13,'(A80)') EVCUnits
        READ(13,*) EVCNoOfVals
        READ(13,*) EVCDefVal
        CASE &
        ("min_and_max_primitive_weighting_in_continuity_equation")
        FoundTau0MinMax = .TRUE. 
        READ(13,'(A80)') Tau0MinMaxUnits
        READ(13,*) Tau0MinMaxNoOfVals
        READ(13,*) (Tau0MinMaxDefVal(j),j=1,Tau0MinMaxNoOfVals)
    ! Casey 100210: Allow SWAN to handle wave refraction as a nodal attribute.
        CASE("wave_refraction_in_swan")
        FoundSwanWaveRefrac = .TRUE. 
        READ(13,'(A80)') SwanWaveRefracUnits
        READ(13,*) SwanWaveRefracNoOfVals
        READ(13,*) SwanWaveRefracDefVal
    ! Corbitt 120321: Allow advection to be turned on locally instead of globally
        CASE("advection_state")
        FoundAdvectionState = .TRUE. 
        READ(13,'(A80)') AdvectionStateUnits
        READ(13,*) AdvectionStateNoOfVals
        READ(13,*) AdvectionStateDefVal
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
        CASE DEFAULT
        WRITE(16,1001)          ! Nodal Attributes file
        WRITE(16,1021) AttrName ! contains invalid name
        IF (NScreen /= 0 .AND. MyProc == 0) THEN
            WRITE(ScreenUnit,1001)
            WRITE(ScreenUnit,1021) AttrName
        ENDIF
        READ(13,'(A80)') Skipped  ! skip the Units for the invalid name
        READ(13,'(A80)') Skipped  ! skip the NoOfVals for invalid name
        READ(13,'(A80)') Skipped  !jgf51.40: skip the default value
        END SELECT
    END DO

!     Determine if there are any attributes required by the fort.15 file
!     that are not in the nodal attributes file.
    call checkForMissingNodalAttributes()

!     Allocate memory for nodal attributes
    call allocateNodalAttributes()

!     Now read each of the attributes required by the model parameter
!     (unit 15) file and skip past the others.
    WRITE(16,270) NWP
    270 FORMAT(/,9X,'Now reading ',I2,' nodal attribute(s).')
    DO k=1, NAttr
        WRITE(16,280) k
        280 FORMAT(/,9X,'Attribute ',I2,':')
        READ(13,*) AttrName
        READ(13,*) NumNodesNotDefault
        WRITE(16,'(14X,A80)') AttrName
        SELECT CASE (AttrName)
        CASE("primitive_weighting_in_continuity_equation")
        IF (LoadTau0) THEN
            CALL LoadAttrVec(TAU0VAR, Tau0DefVal, &
            NumNodesNotDefault, NScreen, MyProc, NAbOut)
        ELSE
            SkipDataSet = .TRUE. 
        ENDIF
        CASE("surface_submergence_state")
        IF (LoadStartDry) THEN
            CALL LoadAttrVec(STARTDRY, StartDryDefVal, &
            NumNodesNotDefault, NScreen, MyProc, NAbOut)
        ELSE
            SkipDataSet = .TRUE. 
        ENDIF
        CASE("quadratic_friction_coefficient_at_sea_floor")
        IF (LoadQuadraticFric) THEN
            CALL LoadAttrVec(FRIC, QuadraticFricDefVal, &
            NumNodesNotDefault, NScreen, MyProc, NAbOut)
        ELSE
            SkipDataSet = .TRUE. 
        ENDIF
        CASE("surface_directional_effective_roughness_length")
        IF (LoadDirEffRLen) THEN
            CALL LoadAttrMat(z0land, DirEffRLenNoOfVals, &
            DirEffRLenDefVal, NumNodesNotDefault, &
            NScreen, MyProc, NAbOut)
        ELSE
            SkipDataSet = .TRUE. 
        ENDIF
        CASE("surface_canopy_coefficient")
        IF (LoadCanopyCoef) THEN
            CALL LoadAttrVec(vcanopy, CanopyCoefDefVal, &
            NumNodesNotDefault, NScreen, MyProc, NAbOut)
        ELSE
            SkipDataSet = .TRUE. 
        ENDIF
        CASE("bridge_pilings_friction_parameters")
        IF (LoadBridgePilings) THEN
            CALL LoadAttrMat(BridgePilings, BridgePilingsNoOfVals, &
            BridgePilingsDefVal,  NumNodesNotDefault, &
            NScreen, MyProc, NAbOut)
        ELSE
            SkipDataSet = .TRUE. 
        ENDIF
        CASE("mannings_n_at_sea_floor")
        IF (LoadManningsN) THEN
            CALL LoadAttrVec(ManningsN, ManningsNDefVal, &
            NumNodesNotDefault, NScreen, MyProc, NAbOut)
        ELSE
            SkipDataSet = .TRUE. 
        ENDIF
        CASE("bottom_roughness_length")
        IF (LoadZ0b_var) THEN
            CALL LoadAttrVec(Z0b_var, Z0b_varDefVal, &
            NumNodesNotDefault, NScreen, MyProc, NAbOut)
        ELSE
            SkipDataSet = .TRUE. 
        ENDIF
        CASE("chezy_friction_coefficient_at_sea_floor")
        IF (LoadChezy) THEN
            CALL LoadAttrVec(Chezy, ChezyDefVal, &
            NumNodesNotDefault, NScreen, MyProc, NAbOut)
        ELSE
            SkipDataSet = .TRUE. 
        ENDIF
        CASE("sea_surface_height_above_geoid")
        IF (LoadGeoidOffset) THEN
            CALL LoadAttrVec(GeoidOffset, GeoidOffsetDefVal, &
            NumNodesNotDefault, NScreen, MyProc, NAbOut)
        ELSE
            SkipDataSet = .TRUE. 
        ENDIF
        CASE &
        ("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")
        IF (LoadEVM) THEN
            CALL LoadAttrVec(EVM, EVMDefVal, &
            NumNodesNotDefault, NScreen, MyProc, NAbOut)
        ELSE
            SkipDataSet = .TRUE. 
        ENDIF
        CASE &
        ("average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth")
        IF (LoadEVC) THEN
            CALL LoadAttrVec(EVC, EVCDefVal, &
            NumNodesNotDefault, NScreen, MyProc, NAbOut)
        ELSE
            SkipDataSet = .TRUE. 
        ENDIF
        CASE &
        ("min_and_max_primitive_weighting_in_continuity_equation")
        IF (LoadTau0MinMax) THEN
            CALL LoadAttrMat(Tau0MinMax, Tau0MinMaxNoOfVals, &
            Tau0MinMaxDefVal, NumNodesNotDefault, &
            NScreen, MyProc, NAbOut)
        ELSE
            SkipDataSet = .TRUE. 
        ENDIF
    ! Casey 100210: Allow SWAN to handle wave refraction as a nodal attribute.
        CASE("wave_refraction_in_swan")
        IF (LoadSwanWaveRefrac) THEN
            CALL LoadAttrVec(SwanWaveRefrac, SwanWaveRefracDefVal, &
            NumNodesNotDefault, NScreen, MyProc, NAbOut)
        ELSE
            SkipDataSet = .TRUE. 
        ENDIF
    ! Corbitt 120321: Allow Advection to be handled locally instead of globally.
        CASE("advection_state")
        IF (LoadAdvectionState) THEN
            CALL LoadAttrVec(AdvectionState,AdvectionStateDefVal, &
            NumNodesNotDefault, NScreen, MyProc, NAbOut)
        ELSE
            SkipDataSet = .TRUE. 
        ENDIF
        CASE("elemental_slope_limiter")
        IF( LoadEleSlopeLim ) THEN
            CALL LoadAttrVec(elemental_slope_limiter_grad_max, &
            EleSlopeLimDefVal, NumNodesNotDefault, NScreen, &
            MyProc, NAbOut)
        ELSE
            SkipDataSet = .TRUE. 
        ENDIF
        CASE("initial_river_elevation")   ! tcm 20140502 v51.27
        IF( LoadRiver_et_WSE ) THEN
            CALL LoadAttrVec(River_et_WSE, &
            River_et_WSEDefVal, NumNodesNotDefault, NScreen, &
            MyProc, NAbOut)
        ELSE
            SkipDataSet = .TRUE. 
        ENDIF
        CASE DEFAULT
        SkipDataSet = .TRUE. 
        WRITE(16,1001)      ! Nodal Attributes file
        WRITE(16,1021) AttrName ! contains invalid name
        IF (NScreen /= 0 .AND. MyProc == 0) THEN
            WRITE(ScreenUnit,1001)
            WRITE(ScreenUnit,1021) AttrName
        ENDIF
        END SELECT
        IF (SkipDataSet) THEN
            DO L=1, NumNodesNotDefault
                READ(13,*) Skipped
            END DO
            WRITE(16,'(9X,A8)') 'Skipped.'
            SkipDataSet = .FALSE. 
        ELSE
            WRITE(16,'(/,9X,A18,A80)') 'Finished loading ', AttrName
        ENDIF
    END DO

    1000 FORMAT('ERROR: The Model Parameter File (unit 15)')
    1001 FORMAT('ERROR: The Nodal Attributes File (unit 13)')
    1002 FORMAT('ERROR: The legacy StartDry File (unit 12)')
    1003 FORMAT('ERROR: Spatially Varying Fric. Coeff. File (unit 21)')

    1005 FORMAT('exists but cannot be opened.')
    1011 FORMAT('was not found.')
    1021 FORMAT('contains invalid name: ',A80)
    9972 FORMAT(////,1X,'!!!!!!!!!! INPUT ERROR !!!!!!!!!',/)
    9973 FORMAT(/,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)

    RETURN
!     ----------------------------------------------------------------
    END SUBROUTINE ReadNodalAttr
!     ----------------------------------------------------------------


!     ----------------------------------------------------------------
!         S U B R O U T I N E     L O A D  A T T R  V E C
!     ----------------------------------------------------------------

!     jgf46.00 Subroutine to set a single set of nodal attributes to
!     their user-specified default values, then read the nondefault
!     values from the Nodal Attributes File (unit 13). This subroutine
!     is used for nodal attributes with only one value per node, hence
!     the suffix "vec" in the name.

!     ----------------------------------------------------------------
    SUBROUTINE LoadAttrVec(AttributeData, Default, NumNodesNotDef, &
    NScreen, MyProc, NAbOut)
    IMPLICIT NONE
    REAL(SZ), intent(out), dimension(NumOfNodes) :: AttributeData
    REAL(SZ), intent(in):: Default ! default value for all nodes
    INTEGER, intent(in) :: NumNodesNotDef ! number of nodes specified in file
    INTEGER, intent(in) :: NScreen ! 1 for debug info to screen (unit 6)
    INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
    INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16

    INTEGER :: NodeNum            ! node number listed in the file

!     Set all values to user-specified default values.
    IF (NABOUT == 0) WRITE(16,1001) Default
    DO i=1, NumOfNodes
        AttributeData(i) = Default
    END DO

    IF (NABOUT == 0) WRITE(16,1005)
    DO i=1, NumNodesNotDef
        READ(13,*) NodeNum, AttributeData(NodeNum)
        IF (NABOUT == 0) &
        WRITE(16,1010) NodeNum, AttributeData(NodeNum)
    END DO

    1001 FORMAT(/,10X,'Set all nodes to the default value of ',E16.8,/)
    1005 FORMAT(/,10X,'Now setting the following nodes to these values:', &
    /,10X,'NODE',5X,'DATA',5X/)
    1010 FORMAT(7X,I8,6X,E16.8)

    RETURN
!     ----------------------------------------------------------------
    END SUBROUTINE LoadAttrVec
!     ----------------------------------------------------------------


!     ----------------------------------------------------------------
!         S U B R O U T I N E     L O A D  A T T R  M A T
!     ----------------------------------------------------------------

!     jgf46.00 Subroutine to load a single set of nodal attributes from
!     the Nodal Attributes File (unit 13) if there is more than one
!     value per node.

!     ----------------------------------------------------------------
    SUBROUTINE LoadAttrMat(AttributeData, NumCol, Default, &
    NumNodesNotDef, NScreen, MyProc, NAbOut)
    IMPLICIT NONE
    INTEGER, intent(in) :: NumCol  ! number of columns in the matrix
    REAL(SZ), intent(out), &
    dimension(NumOfNodes,NumCol) :: AttributeData
    REAL(SZ), intent(in), dimension(NumCol) :: Default ! default values
    INTEGER, intent(in) :: NumNodesNotDef  ! number of nodes spec. in file
    INTEGER, intent(in) :: NScreen ! 1 for debug info to screen (unit 6)
    INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
    INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16

    INTEGER :: NodeNum            ! node number listed in the file

!     Set all nodes to user-specified default values.
    IF (NABOUT == 0) WRITE(16,1001)
    DO i=1, NumOfNodes
        DO j=1, NumCol
            AttributeData(i,j)=Default(j)
        END DO
    END DO

    IF (NABOUT == 0) WRITE(16,1005)
    DO i=1, NumNodesNotDef
        READ(13,*) NodeNum, (AttributeData(NodeNum,j),j=1,NumCol)
        IF (NABOUT == 0) WRITE(16,1010) NodeNum, &
        (AttributeData(NodeNum,j),j=1,NumCol)
    END DO

    1001 FORMAT(/,10X,'Set all nodes to the default values of ',/, &
    &      99E16.8,/)
    1005 FORMAT(/,10X,'Now setting the following nodes to these values:', &
    /,10X,'NODE',5X,'DATA',5X/)
    1010 FORMAT(7X,I6,6X,12(1X,E16.8))

    RETURN
!     ----------------------------------------------------------------
    END SUBROUTINE LoadAttrMat
!     ----------------------------------------------------------------



!     ----------------------------------------------------------------
!         S U B R O U T I N E     I N I T  N O D A L  A T T R
!     ----------------------------------------------------------------

!     jgf46.00 Subroutine to initialize and error check the nodal
!     attributes read in from the Nodal Attributes File (unit 13).

!     ----------------------------------------------------------------
    SUBROUTINE InitNodalAttr(DP, NP, G, NScreen, ScreenUnit, &
    MyProc, NAbOut)
    USE GLOBAL, ONLY : C3D, C2DDI
    USE GLOBAL_3DVS, ONLY : Z0B
    IMPLICIT NONE
    INTEGER, intent(in) :: NP ! number of nodes in the grid file
    REAL(SZ), intent(in), dimension(NP) :: DP ! array of bathymetric depths
    REAL(SZ), intent(in):: G  ! gravitational acceleration
    INTEGER, intent(in) :: NScreen ! nonzero for debug info to screen
    INTEGER, intent(in) :: ScreenUnit ! i/o for debug info to screen
    INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
    INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16
    INTEGER :: Tau0Dig1  ! determines the tau0 scheme
    INTEGER :: Tau0Dig2  ! determines whether tau0 is being output
    LOGICAL :: invalidCanopyCoefficient ! .TRUE. if any value is neither 0 nor 1

    IF (Tau0 < 0) THEN
        Tau0Dig1 = INT(Tau0)   ! jgf47.30 truncate the fractional part
    ! jgf47.34 round away from zero by subtracting 0.5d0
        Tau0Dig2 = INT( (Tau0 - REAL(Tau0Dig1))*10.0d0 - 0.5d0)
    ELSE
        Tau0Dig1 = 0
        Tau0Dig2 = 0
    ENDIF

!     ERROR CHECK: If a nodal attributes file is being used, check to
!     see that the number of nodes in the nodal attribute file is the
!     same as the number of nodes in the grid file.
    IF (NWP /= 0 .AND. NumOfNodes /= NP) THEN
        IF(NSCREEN /= 0 .AND. MYPROC == 0) WRITE(ScreenUnit,9900)
        WRITE(16,9900)
        9900 FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!', &
        //,1X,'The number of nodes in the grid file (unit 14) and' &
        /,1X,'the nodal attributes file (unit 13) must match.', &
        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
        STOP                   ! We're toast.
    ENDIF

!     ERROR CHECK: If Chezy, Manning's or Quadratic friction was loaded
!     from the nodal attributes file, NOLIBF must be 1.
    IF ((LoadChezy .OR. LoadManningsN .OR. LoadQuadraticFric) .AND. &
    NoLiBF /= 1) THEN
        IF(NSCREEN /= 0 .AND. MYPROC == 0) WRITE(ScreenUnit,9800)
        WRITE(16,9800)
        9800 FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!', &
        //,1X,'Nonlinear bottom friction coefficients were loaded' &
        /,1X,'from the nodal attributes file (unit 13), so ', &
        /,1X,'NoLiBF must be set to 1. It is set to ',i2,' in', &
        /,1X,'the model parameter (unit 15) file.', &
        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
        STOP                   ! We're toast.
    ENDIF

!     ERROR CHECK: If Tau0=-3.x or -6.x in fort.15, then tau0 MUST be loaded
!     from nodal attributes file.
    IF ( ((Tau0Dig1 == -3) .OR. (Tau0Dig1 == -6)) &
     .AND. ( .NOT. LoadTau0) ) THEN
        IF(NSCREEN /= 0 .AND. MYPROC == 0) WRITE(ScreenUnit,9700) Tau0
        WRITE(16,9700)
        9700 FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!', &
        //,1X,'Spatially and temporally varying tau0 was ' &
        /,1X,'specified in the fort.15 file with Tau0=',E9.2, &
        /,1X,'but the base value was not specified in the ' &
        /,1X,'nodal attributes file (unit 13). Please ', &
        /,1X,'load the base value using', &
        /,1X,'primitive_weighting_in_continuity_equation.', &
        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
        STOP                   ! We're toast.
    ENDIF

!     ERROR CHECK: The canopy coefficient value should be either 0 or 1;
!     log a warning message if any other values were read.
    if ( loadCanopyCoef.eqv. .TRUE. ) then
        invalidCanopyCoefficient = .FALSE. 
        do i=1,np
            if ( (vcanopy(i) /= 0.d0) .AND. (vcanopy(i) /= 1.d0 ) ) then
                invalidCanopyCoefficient = .TRUE. 
                exit
            endif
        end do
        if ( invalidCanopyCoefficient.eqv. .TRUE. ) then
            call allMessage(WARNING,'Some surface canopy coefficient '// &
            'nodal attribute data are neither 0 nor 1.')
        endif
    endif

!     I N I T    S T A R T D R Y
     
    IF (NWP == 0) THEN
        ALLOCATE(STARTDRY(NP))
    ENDIF
    IF (LoadStartDry.eqv. .FALSE. ) THEN
        DO I=1, NP
            STARTDRY(I) = 0.0D0
        ENDDO
    ENDIF

!     I N I T     T A U 0
    IF (NWP == 0) THEN
        ALLOCATE(TAU0VAR(NP),TAU0BASE(NP))
        ALLOCATE(Tau0MinMax(NP,Tau0MinMaxNoOfVals))
    ENDIF

!     jgf46.25 If input tau0 is positive, set all nodes to that value.
    IF ( .NOT. LoadTau0) THEN
        IF (Tau0 >= 0) THEN
            DO I=1,NP
                Tau0Var(I)=Tau0
            END DO
            WRITE(16,7) Tau0
            7 FORMAT(/,5X, &
            'A SPATIALLY CONSTANT WEIGHTING COEFFICIENT (Tau0)' &
            ,/,5X,' WILL BE USED IN THE GENERALIZED WAVE', &
            ' CONTINUITY EQUATION.', &
            /,5X,'Tau0 = ',E15.8,2X,'1/sec',/)
        ELSE
        !           If input tau0 is negative, set value using hardcoded scheme
        !           based on depth
            DO I=1,NP
                Tau0Var(I)=Tau0NodalValue(Tau0,DP(I))
            ENDDO
        ! jgf47.30.TODO: This logging needs to be cleaned up.
            IF(Tau0 == -2) THEN
                WRITE(16,6) ! spatially vary tau0 according to hard coded scheme
                WRITE(16,62) ! description of scheme
                62 FORMAT(/,5X,'IF DEPTH > 200           Tau0 = 0.005', &
                /,5X,'IF 200   > DEPTH > 1     Tau0 = 1/DEPTH  ', &
                /,5X,'IF 1     > DEPTH         Tau0 = 1.0 ')
            ENDIF
            IF ( .NOT. ( (Tau0Dig1 == -3) .OR. (Tau0Dig1 == -5) ) ) THEN
                WRITE(16,6) ! spatially vary tau0 according to hard coded scheme
                WRITE(16,61) ! description of scheme
                61 FORMAT(/,5X,' IF DEPTH GE 10           -> TAU0 = 0.005', &
                /,5X,' IF DEPTH LT 10           -> TAU0 = 0.020',/)
            ENDIF
        ENDIF
    ENDIF
    6 FORMAT(/,5X,'A SPATIALLY VARIABLE WEIGHTING COEFFICIENT (Tau0)' &
    ,/,5X,' WILL BE USED IN THE GENERALIZED WAVE', &
    ' CONTINUITY EQUATION.', &
    /,5x,'THIS VALUE WILL BE DETERMINED AS FOLLOWS:')

!     jgf46.27 If we have already loaded the tau0 values directly from
!     the nodal attributes file, check to see if the default value was
!     negative. If so, this indicates that nodal values of tau0 that
!     were not explicitly set in the nodal attributes file should be set
!     according to one of the hard-coded tau0 schemes.
    IF (LoadTau0 .AND. Tau0DefVal < 0) THEN
        DO I=1,NP
            IF (Tau0Var(I) < 0) THEN
                Tau0Var(I)=Tau0NodalValue(Tau0DefVal,DP(I))
            ENDIF
        ENDDO
    ENDIF

!     jgf47.06 Activate time varying tau0 and output tau0 if these options
!     were selected.

!     jgf47.30 Use Joannes' scheme for steady Tau0 in deep water and
!     other coarsely gridded areas, and time varying Tau0 in high
!     resolution areas.
    IF ( (Tau0Dig1 == -3) .OR. (Tau0Dig1 == -6) &
     .OR. (Tau0Dig1 == -7) ) THEN
        HighResTimeVaryingTau0 = .TRUE. 
        DO I=1, NP
            Tau0Base(I) = Tau0Var(I)
        ENDDO
    ENDIF

!     jgf47.11 Also allow the min and max tau0 to be set from the fort.15,
!     bypassing the use of the fort.13 file for this purpose.
!     jgf47.30 Changed to emphasize full domain time varying tau0
    IF ( Tau0Dig1 == -5 ) THEN
        FullDomainTimeVaryingTau0 = .TRUE. 
        IF ( .NOT. LoadTau0MinMax ) THEN
            DO I=1, NP
                Tau0MinMax(I,1) = Tau0FullDomainMin
                Tau0MinMax(I,2) = Tau0FullDomainMax
            ENDDO
        ENDIF
    ENDIF

!     jgf47.30: Output of tau0 is now activated by having a 0.1 fraction
!     for tau0
    IF ( Tau0Dig2 == -1 ) THEN
        OutputTau0 = .TRUE. 
    ENDIF

!     jjw&sb46.38.sb01 If tau0 is loaded from nodal attributes file and
!     Tau0 is -3, time-varing tau0 optimizer will be applied in timestep.F
    IF (HighResTimeVaryingTau0 .OR. FullDomainTimeVaryingTau0) THEN
        ALLOCATE(Tau0Temp(NP))
        WRITE(16,8) ! jgf47.30.TODO: This logging should be consolidated.
    ENDIF
    8 FORMAT(/,5X,'A SPATIALLY TEMPORALLY VARIABLE OPTIMIZED ' &
    ,/,5X,' WEIGHTING COEFFICIENT (Tau0) WILL BE USED ' &
    ,/,5X,' IN THE GENERALIZED WAVE CONTINUITY EQUATION.',/)

!     jgf47.33 Enable time averaging of tau0 if requested.
    IF (Tau0Dig1 == -6) THEN
        TimeAveragedTau0 = .TRUE. 
    ENDIF

!     jgf48.42 Enable back loaded time averaging of tau0 if requested.
    IF (Tau0Dig1 == -7) THEN
        BackLoadedTimeAveragedTau0 = .TRUE. 
    ENDIF

!     jgf48.46 Allocate array to hold previous value tau0 for use in
!     time averaging, if necessary.
    IF ( TimeAveragedTau0 .OR. BackLoadedTimeAveragedTau0 ) THEN
        ALLOCATE(LastTau0(NP))
        DO I=1, NP
            LastTau0(I) = Tau0Base(I)
        ENDDO
    ENDIF

!     I N I T   B O T T O M   F R I C T I O N
    IF(NOLIBF == 0) THEN
        IFNLBF=0
        IFLINBF=1
        IFHYBF=0
    ENDIF
    IF(NOLIBF == 1) THEN
        IFNLBF=1
        IFLINBF=0
        IFHYBF=0
    ENDIF
    IF(NOLIBF == 2) THEN
        IFNLBF=0
        IFLINBF=0
        IFHYBF=1
    ENDIF

!     Initialize bottom friction if it was not loaded from unit 13.
    IF(C2DDI) THEN
        IF(( .NOT. LoadQuadraticFric) .AND. ( .NOT. LoadManningsN) .AND. &
        ( .NOT. LoadChezy)) THEN
            IF (NoLiBF == 0) CF=Tau
        !     If a nodal attributes file was read, FRIC was allocated there.
            IF (NWP == 0) THEN
                ALLOCATE(FRIC(NP))
            ENDIF
            DO I=1,NP
                FRIC(I)=CF
            END DO
        ENDIF
    
    !        jgf47.04 If a depth-dependent friction parameterization is used,
    !        the value from the fort.15 file is used as a floor for the
    !        minimum equivalent quadratic friction value.
        IF (LoadManningsN) THEN
            BFCdLLimit = CF
        ENDIF
    ENDIF

    IF(C3D) THEN
    !     Initialize 3D bottom roughness if it was not loaded from unit 13.
        IF(( .NOT. LoadZ0b_var)) THEN
        !     If a nodal attributes file was read, Z0b_var was allocated there.
            IF (NWP == 0) THEN
                ALLOCATE(Z0b_var(NP))
                ALLOCATE(FRIC(NP))
            ENDIF
            DO I=1,NP
                Z0b_var(I)=Z0B
                FRIC(I)=CF
            END DO
        ENDIF
    !     jgf47.04 If a depth-dependent friction parameterization is used,
    !     the value from the fort.15 file is used as a floor for the
    !     minimum equivalent quadratic friction value.
        IF (LoadZ0b_var .OR. LoadManningsN) THEN
            BFCdLLimit = CF
        ENDIF
    ENDIF

!     Initialize bridge pilings.

    IF (LoadBridgePilings) THEN
        DO I=1, NP
            IF (BridgePilings(I,1) /= 0) THEN ! only for nodes w/piers
                BridgePilings(I,3) = 4.d0 * &
                BridgePilings(I,3) / BridgePilings(I,4)
            ENDIF
        END DO
    ENDIF

!     I N I T   E D D Y   V I S C O S I T Y  &  D I F F U S I V I T Y
    IF ( .NOT. LoadEVM) THEN
        IF (NWP == 0) THEN
            ALLOCATE(EVM(NP))
        ENDIF
        DO I=1,NP
            EVM(I)=ESLM
        END DO
    ENDIF
    IF ( .NOT. LoadEVC .AND. ESLC /= 0) THEN
        IF (NWP == 0) THEN
            ALLOCATE(EVC(NP))
        ENDIF
        DO I=1,NP
            EVC(I)=ESLC
        END DO
    ENDIF

    RETURN
!     ----------------------------------------------------------------
    END SUBROUTINE InitNodalAttr
!     ----------------------------------------------------------------


!     ----------------------------------------------------------------
!      F U N C T I O N   T A U 0  N O D A L  V A L U E
!     ----------------------------------------------------------------

!     jgf46.27 Function to calculate tau0 based on the scheme selection
!     and the depth. This assumes that Scheme is negative.

!     ----------------------------------------------------------------
    REAL(SZ) FUNCTION Tau0NodalValue(Scheme, Depth)
    IMPLICIT NONE
    REAL(SZ) Scheme
    REAL(SZ) Depth

    IF (Scheme == -2.d0) THEN
    !     Smoothly varying tau0 with depth.
        IF(Depth >= 200.) Tau0NodalValue=0.005
        IF((Depth < 200.) .AND. (Depth >= 1.)) THEN
            Tau0NodalValue=1./Depth
        ENDIF
        IF(Depth < 1.) Tau0NodalValue=1.0
    ELSE
    !     Abrupt variation in tau0 with depth.
        IF(Depth <= 10.) Tau0NodalValue=0.020d0
        IF(Depth > 10.) Tau0NodalValue=0.005d0
    ENDIF
!     ----------------------------------------------------------------
    END FUNCTION Tau0NodalValue
!     ----------------------------------------------------------------


!     ----------------------------------------------------------------
!      S U B R O U T I N E
!         C A L C U L A T E  T I M E  V A R Y I N G  T A U 0
!     ----------------------------------------------------------------

!     jgf47.08 Subroutine to calculate a new tau0 value. Called from
!     GWCE_New in timestep.F each time the GWCE matrix is reset (i.e.,
!     upon startup and whenever wetting and/or drying occurs in any
!     subdomain. Based on Casey050711.

!     ----------------------------------------------------------------
    SUBROUTINE CalculateTimeVaryingTau0(TK, NNeigh, NeiTab, NP)
    IMPLICIT NONE
    REAL(SZ), intent(in) :: TK(:)       ! bottom friction
    INTEGER, intent(in) :: NNeigh(:)   ! number of neighbor nodes
    INTEGER, intent(in) :: NeiTab(:,:) ! table of neighbor nodes
    INTEGER, intent(in) :: NP           ! number of nodes in the domain
    REAL(SZ) CaseySum  ! sum of tau0temp values around a particular node

! Casey 050711 : Made changes for averaged variable Tau0.
! jw46.39.sb01 :  "high/low LIMITED" variable G.
! jgf47.30: Distinction between fulldomain and hi res only
    IF ( FullDomainTimeVaryingTau0 ) THEN
        DO i = 1, NP
            Tau0Temp(i)=Tau0MinMax(i,1)+1.5*TK(i)
            IF (Tau0Temp(i) < Tau0MinMax(i,1)) THEN
                Tau0Temp(i)=Tau0MinMax(i,1)
            ENDIF
            IF(Tau0Temp(i) > Tau0MinMax(i,2)) THEN
                Tau0Temp(i)=Tau0MinMax(i,2)
            ENDIF
        ENDDO
    ENDIF
    IF ( HighResTimeVaryingTau0 ) THEN
        DO i = 1, NP
            IF(Tau0Base(i) < 0.025) THEN
                Tau0Temp(I)=Tau0Base(i) ! not time varying
            ELSE
                Tau0Temp(i)=Tau0Base(i)+1.5*TK(i) ! time varying
                IF (Tau0Temp(i) > 0.2) Tau0Temp(i)=0.2 ! ceiling
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

!     jgf47.33 Perform time averaging of tau0 if requested.
    IF (TimeAveragedTau0) THEN
        DO I=1, NP
            TAU0VAR(I) = 0.5d0*TAU0VAR(I) + 0.5d0*LastTau0(I)
            LastTau0(I) = TAU0VAR(I)
        ENDDO
    ENDIF

!     jgf48.42 Perform backloaded time averaging of tau0 if requested.
    IF (BackLoadedTimeAveragedTau0) THEN
        DO I=1, NP
            TAU0VAR(I) = AlphaTau0*TAU0VAR(I) &
            + (1.d0-AlphaTau0)*LastTau0(I)
            LastTau0(I) = TAU0VAR(I)
        ENDDO
    ENDIF
!     ----------------------------------------------------------------
    END SUBROUTINE CalculateTimeVaryingTau0
!     ----------------------------------------------------------------



!     ----------------------------------------------------------------
!     S U B R O U T I N E
!     A P P L Y  2 D  B O T T O M  F R I C T I O N
!     ----------------------------------------------------------------

!     jgf46.00 Subroutine to apply 2D bottom friction from turbulent
!     viscous effects as well as bridge pilings. This is used in the
!     time stepping loop.

!     ----------------------------------------------------------------

!     sb46.28sb02 Lower limit of Cd was added as an argument.
!     jgf47.04 Argument for lower limit of Cd was removed; this value
!     is now specified by the user in fort.15.

!     ----------------------------------------------------------------
    SUBROUTINE Apply2DBottomFriction(UU1, VV1, DP, ETA2, G, &
    IFNLFA, NP, TK)
    USE SIZES
    IMPLICIT NONE
    INTEGER, intent(in) :: NP                   ! number of nodes in grid
    REAL(SZ), intent(in), dimension(NP) :: UU1  ! x-dir velocities
    REAL(SZ), intent(in), dimension(NP) :: VV1  ! y-dir velocities
    REAL(SZ), intent(in), dimension(NP) :: DP   ! bathymetric depths
    REAL(SZ), intent(in), dimension(NP) :: ETA2 ! water surf. elevations
    REAL(SZ), intent(in) :: G                   ! gravitational constant
    INTEGER, intent(in) :: IFNLFA               ! nonlin. finite amp. flag
    REAL(SZ), intent(inout), dimension(NP) :: TK! depth avg. fric.

    REAL(SZ) UV1   ! velocity magnitude (speed)
    REAL(SZ) H1    ! total depth
    REAL(SZ) Fr
    REAL(SZ) FricBP
    REAL(SZ) BK    ! BK(1) is pier shape factor
    REAL(SZ) BALPHA! BALPHA(2) is constriction fraction
    REAL(SZ) BDELX ! BDELX(3) is effective delx

!     Step 0. Convert Manning's N to Cd, if necessary.
    IF (LoadManningsN) THEN
        DO I=1, NP
            FRIC(I)=g*ManningsN(I)**2.d0 &
            /( ( DP(I)+IFNLFA*ETA2(I) )**(1.d0/3.d0) ) ! sb46.28sb02
        ! b46.28sb02  Lower limit is applied here.
            IF(FRIC(I) < BFCdLLimit) THEN
                FRIC(I) = BFCdLLimit
            ENDIF
        ENDDO
    ENDIF

!     ... Convert Chezy to Cd, if necessary.
    IF (LoadChezy) THEN
        DO I=1,NP
            FRIC(I)=G/(Chezy(I)**2)
        END DO
    ENDIF


!     Step 1. Apply friction arising from turbulent viscous interaction
!     with the sea floor.
    DO I=1, NP
        UV1=SQRT(UU1(I)*UU1(I)+VV1(I)*VV1(I))
        H1=DP(I)+IFNLFA*ETA2(I)
        TK(I)= FRIC(I)* &
        ( IFLINBF +        & ! linear
        (UV1/H1) * (IFNLBF  & ! nonlinear
        + IFHYBF*(1+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA))) ! hybrid
    END DO

!     Step 2. Apply friction arising from flow interaction with bridge
!     pilings, if required.
    IF (LoadBridgePilings) THEN
        DO I=1, NP
            UV1=SQRT(UU1(I)*UU1(I)+VV1(I)*VV1(I))
            H1=DP(I)+IFNLFA*ETA2(I)
            Fr=UV1*UV1/(G*H1)
            BK = BridgePilings(I,1)
            BALPHA = BridgePilings(I,2)
            BDELX = BridgePilings(I,3)
            FricBP=(H1/BDELX)*BK*(BK+5.d0*Fr*Fr-0.6d0) &
            *(BALPHA+15.d0*BALPHA**4)
            TK(I)=TK(I)+FricBP*UV1/H1
        END DO
    ENDIF

    RETURN
!     ----------------------------------------------------------------
    END SUBROUTINE Apply2DBottomFriction
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!     S U B R O U T I N E
!     A P P L Y  3 D  B O T T O M  F R I C T I O N
!     ----------------------------------------------------------------

!     jgf46.00 Subroutine to apply 3D bottom friction from turbulent
!     viscous effects as well as bridge pilings. This is used in the
!     time stepping loop.
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
    SUBROUTINE Apply3DBottomFriction(Q, SIGMA, DP, ETA2, G, &
    IFNLFA, NP, TK, NFEN)
    USE SIZES
    USE GLOBAL_3DVS, ONLY : Z0B
    IMPLICIT NONE
    INTEGER, intent(in) :: NP, NFEN             ! number of nodes in grid Horizontal and Vertical
    COMPLEX(SZ), intent(in), dimension(NP,NFEN) :: Q  ! x-dir velocities
    REAL(SZ), intent(in), dimension(NFEN) :: SIGMA  ! x-dir velocities
    REAL(SZ), intent(in), dimension(NP) :: DP   ! bathymetric depths
    REAL(SZ), intent(in), dimension(NP) :: ETA2 ! water surf. elevations
    REAL(SZ), intent(in) :: G                   ! gravitational constant
    INTEGER, intent(in) :: IFNLFA               ! nonlin. finite amp. flag
    REAL(SZ), intent(inout), dimension(NP) :: TK! depth avg. fric.

    INTEGER :: NH
    REAL(SZ) Z0B1  ! velocity magnitude (speed)
    REAL(SZ) UV1   ! velocity magnitude (speed)
    REAL(SZ) H1    ! total depth
    REAL(SZ) Fr
    REAL(SZ) FricBP
    REAL(SZ) BK    ! BK(1) is pier shape factor
    REAL(SZ) BALPHA! BALPHA(2) is constriction fraction
    REAL(SZ) BDELX ! BDELX(3) is effective delx

! Determine the bottom roughness length either from fort.15, from Manning's n
! or as read in from nodal attributes
    DO NH=1,NP
        H1=DP(NH)+IFNLFA*ETA2(NH)
        IF (LoadZ0B_var) THEN
            Z0B1 = Z0B_var(NH)
        ELSEIF (LoadManningsN) THEN
            Z0B1 = ( H1 )* exp(-(1.0D0+ &
            ( (0.41D0*( H1 )**(1.0D0/6.0D0) )/ &
            (ManningsN(NH)*sqrt(g)) ) ))
        ELSE
            Z0B1 = Z0B
        ENDIF

        FRIC(NH)= (1.D0 / ( (1.D0/0.41D0) * &
        LOG((ABS( ( ( SIGMA(2)-SIGMA(1) )/2.d0 ) *(H1) ) + Z0B1 )/Z0B1) &
        ) )**2.D0

        TK(NH)= FRIC(NH) * ABS(Q(NH,1))

        IF (LoadBridgePilings) THEN
            Fr=ABS(Q(NH,1))*ABS(Q(NH,1))/(G*H1)
            BK = BridgePilings(I,1)
            BALPHA = BridgePilings(I,2)
            BDELX = BridgePilings(I,3)
            FricBP=(H1/BDELX)*BK*(BK+5.d0*Fr*Fr-0.6d0) &
            *(BALPHA+15.d0*BALPHA**4)
            TK(I)=TK(I)+FricBP*ABS(Q(NH,1))/H1
        ENDIF
    ENDDO

    RETURN
!     ----------------------------------------------------------------
    END SUBROUTINE Apply3DBottomFriction
!     ----------------------------------------------------------------


!     ----------------------------------------------------------------
!     S U B R O U T I N E
!     A P P L Y  D I R E C T I O N A L  W I N D  R E D U C T I O N
!     ----------------------------------------------------------------

!     jgf46.00 Subroutine to calculate the land wind reduction factor
!     based on a table of directional wind drag values. Originally
!     written into the hstart.F file by jjw in jjw-42.06j. This is used
!     in hstart.F and timestep.F.

!     jgf49.1001 Extracted the application of the canopy coefficient and
!     placed in the ApplyCanopyCoefficient subroutine.
!     ----------------------------------------------------------------
    SUBROUTINE ApplyDirectionalWindReduction(NodeNumber, WindDragCo, &
    WindMag, BathymetricDepth, Elevation, CutOffDepth, G, &
    WindX, WindY)
    USE SIZES
    IMPLICIT NONE
    INTEGER,  intent(in) :: NodeNumber ! index of node under consideration
    REAL(SZ), intent(in) :: WindDragCo ! wind drag coefficient
    REAL(SZ), intent(in) :: WindMag    ! wind magnitude
    REAL(SZ), intent(in) :: BathymetricDepth ! a.k.a. dp(i),depth below geoid
    REAL(SZ), intent(in) :: Elevation  ! a.k.a. eta2(i)
    REAL(SZ), intent(in) :: CutOffDepth! a.k.a. h0, user-spec. min. depth
    REAL(SZ), intent(in) :: G          ! gravitational constant

    REAL(SZ), intent(inout) :: WindX   ! x-dir component of wind velocity
    REAL(SZ), intent(inout) :: WindY   ! x-dir component of wind velocity

    REAL(SZ) z0m   ! marine roughness coefficient based on Garratt's formula
    REAL(SZ) angle ! direction wind is coming from
    INTEGER :: idir   ! code for wind direction
    REAL(SZ) z0l   ! drag for a particular node, for particular direction
    REAL(SZ) TotalDepth  ! bathymetric depth + sea surface elevation
    REAL(SZ) fr    ! land wind reduction factor

!     compute marine roughness coefficient based on Garratt's formula
    z0m=(0.018d0/G)*WindDragCo*WindMag**2.d0

!     compute direction  that the wind is coming from
    if((WindX == 0) .AND. (WindY == 0))then
        angle=0.d0
    else
        angle=atan2(WindY,WindX)
    endif
    angle=360.*angle/(2*3.141592654d0)
    idir=0
    if((angle > -15.) .AND. (angle <= 15))  idir=1
    if((angle > 15.) .AND. (angle <= 45))   idir=2
    if((angle > 45.) .AND. (angle <= 75))   idir=3
    if((angle > 75.) .AND. (angle <= 105))  idir=4
    if((angle > 105.) .AND. (angle <= 135)) idir=5
    if((angle > 135.) .AND. (angle <= 165)) idir=6
    if((angle > 165.) .AND. (angle <= 180)) idir=7
    if((angle > -45.) .AND. (angle <= -15)) idir=12
    if((angle > -75.) .AND. (angle <= -45)) idir=11
    if((angle > -105.) .AND. (angle <= -75)) idir=10
    if((angle > -135.) .AND. (angle <= -105)) idir=9
    if((angle > -165.) .AND. (angle <= -135)) idir=8
    if((angle >= -180.) .AND. (angle <= -165)) idir=7

!     define land roughness from usace values
    z0l=z0land(NodeNumber,idir)

!     reset z0l depending on situation
    if(z0l <= 0.006) then
    !     coe set their value to a marine value -> reset to correct marine value
        z0l=z0m
    else
    !     coe set their value to a land value -> proceed with checking this value
        TotalDepth = BathymetricDepth + Elevation
        if( (TotalDepth > 2*CutOffDepth) .AND. &
        (BathymetricDepth < 20)) then
        !     compute adjusted z0l to account for overland flooding - do this only
        !     in the case where the water column is greater than twice h0 and
        !     you are not in a river (I assume that rivers are deeper than 20m
        !     and have z0l>0.006)
            z0l=z0l-TotalDepth/30. ! correction for overland flooding
        endif
    endif

!     compute land wind reduction factor
    if(z0l > 0.0001) then
        fr=(z0m/z0l)**0.0706d0
    else
        fr=1.000d0
    endif
    if(fr > 1.0000d0) fr=1.0000d0
!     adjust time interpolated wind field
    WindX = fr*WindX
    WindY = fr*WindY

    RETURN
!     ----------------------------------------------------------------
    END SUBROUTINE ApplyDirectionalWindReduction
!     ----------------------------------------------------------------


!     ----------------------------------------------------------------
!     S U B R O U T I N E
!     A P P L Y   C A N O P Y   C O E F F I C I E N T
!     ----------------------------------------------------------------
!     jgf49.1001 Subroutine to apply the canopy coefficient. Was
!     originally included in the subroutine ApplyDirectionalWindReduction;
!     this combination implicitly assumed that the two attributes would
!     always be used together.
!     ----------------------------------------------------------------
    SUBROUTINE ApplyCanopyCoefficient(NodeNumber, WindX, WindY)
    IMPLICIT NONE
    INTEGER,  intent(in) :: NodeNumber ! index of node under consideration
    REAL(SZ), intent(inout) :: WindX   ! x-dir component of wind velocity
    REAL(SZ), intent(inout) :: WindY   ! x-dir component of wind velocity

    WindX = vcanopy(NodeNumber)*WindX
    WindY = vcanopy(NodeNumber)*WindY

    RETURN
!     ----------------------------------------------------------------
    END SUBROUTINE ApplyCanopyCoefficient
!     ----------------------------------------------------------------


!     ----------------------------------------------------------------
!                      S U B R O U T I N E
!         R E A D  L E G A C Y  S T A R T  D R Y  F I L E
!     ----------------------------------------------------------------

!     jgf46.00 Subroutine to load up the legacy startdry file (unit
!     12). This is just a cut-and-paste from the section of the
!     READ_INPUT subroutine that did the same thing. This subroutine is
!     never called. It is vestigial and listed here purely as reference
!     material.

!     ----------------------------------------------------------------
    SUBROUTINE ReadLegacyStartDryFile(NP, NScreen, ScreenUnit, &
    MyProc, NAbOut)
    IMPLICIT NONE
    INTEGER, intent(in) :: NP ! number of nodes in grid file
    INTEGER, intent(in) :: NScreen ! nonzero for debug info to screen
    INTEGER, intent(in) :: ScreenUnit ! i/o for debug info to screen
    INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
    INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16

    INTEGER :: JKI          ! node number from file
    INTEGER :: NE2          ! number of elements, according to fort.12 file
    INTEGER :: NP2          ! number of nodes, according to fort.12 file

    CHARACTER(len=80) AGRID2 ! users comment/description line
    REAL(SZ) DUM1, DUM2  ! data that we want to skip

    OPEN(12,FILE=TRIM(INPUTDIR)//'/'//'fort.12')

!...  READ STARTDRY INFORMATION FROM UNIT 12
    READ(12,'(A80)') AGRID2
    WRITE(16,2038) AGRID2
    2038 FORMAT(5X,'STARTDRY FILE IDENTIFICATION : ',A80,/)
    READ(12,*) NE2,NP2

!...  CHECK THAT NE2 AND NP2 MATCH WITH GRID FILE
!      IF((NE2.NE.NE).OR.(NP2.NE.NP)) THEN
    IF(NP2 /= NP) THEN
        IF(NSCREEN /= 0 .AND. MYPROC == 0) WRITE(ScreenUnit,9900)
        WRITE(16,9900)
        9900 FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!', &
        //,1X,'THE PARAMETER NE2 AND NP2 MUST MATCH NE AND NP ', &
        /,1X,'USER MUST CHECK FORT.12 INPUT FILE ', &
        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
        STOP
    ENDIF

!...  READ IN STARTDRY CODE VALUES
    DO I=1,NP
        READ(12,*) JKI,DUM1,DUM2,STARTDRY(JKI)
        IF(JKI /= I) THEN
            IF(NSCREEN /= 0 .AND. MYPROC == 0) WRITE(ScreenUnit,99805)
            WRITE(16,99805)
            99805 FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ', &
            'INPUT ERROR  !!!!!!!!!', &
            //,1X,'YOUR NODE NUMBERING IS NOT SEQUENTIAL ', &
            'CHECK YOUR UNIT 12 INPUT FILE CAREFULLY',//)
        ENDIF
    END DO

!...  CLOSE UNIT 12 FILE
    CLOSE(12)

    RETURN
!     ----------------------------------------------------------------
    END SUBROUTINE ReadLegacyStartDryFile
!     ----------------------------------------------------------------


!     ----------------------------------------------------------------
!                       S U B R O U T I N E
!       R E A D  L E G A C Y  B O T T O M  F R I C T I O N  F I L E
!     ----------------------------------------------------------------

!     jgf46.00 Subroutine to load up the legacy Spatially Varying
!     Friction Coefficient File (unit 21). This is just a cut-and-paste
!     from the section of the READ_INPUT subroutine that did the same
!     thing. This subroutine is never called. It is vestigial and listed
!     here purely as reference material.

!     ----------------------------------------------------------------
    SUBROUTINE ReadLegacyBottomFrictionFile(NP, NScreen, ScreenUnit, &
    MyProc, NAbOut)
    IMPLICIT NONE
    INTEGER, intent(in) :: NP ! number of nodes in grid file
    INTEGER, intent(in) :: NScreen ! nonzero for debug info to screen
    INTEGER, intent(in) :: ScreenUnit ! i/o for debug info to screen
    INTEGER, intent(in) :: MyProc  ! in parallel, only MyProc=0 i/o to screen
    INTEGER, intent(in) :: NAbOut  ! 1 to abbrev. output to unit 16

    CHARACTER(len=80) AFRIC  ! user's comment/description line
    INTEGER :: NHG    ! node number from file

    OPEN(21,FILE=TRIM(INPUTDIR)//'/'//'fort.21')
    READ(21,'(A80)') AFRIC
    DO I=1,NP
        READ(21,*) NHG,FRIC(NHG)
        IF(NHG /= I) THEN
            IF(NSCREEN /= 0 .AND. MYPROC == 0) WRITE(ScreenUnit,99803)
            WRITE(16,99803)
            99803 FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ', &
            'INPUT ERROR  !!!!!!!!!',//,1X, &
            'YOUR NODAL FRICTION NUMBERING IS NOT SEQUENTIAL ', &
            /,1X,'CHECK YOUR UNIT 21 INPUT FILE CAREFULLY',//,1X, &
            '!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
            STOP
        ENDIF
    END DO
    WRITE(16,3601) AFRIC
    3601 FORMAT(/,5X,'FRICTION FILE IDENTIFICATN : ',A80,/)
    IF(NABOUT /= 1) THEN
        WRITE(16,2080)
        2080 FORMAT(/,10X,'NODE',5X,'BOTTOM FRICTION FRIC',5X,/)
        DO I=1,NP
            WRITE(16,2087) I,FRIC(I)
            2087 FORMAT(7X,I6,6X,E17.10)
        END DO
    ELSE
        WRITE(16,3504)
        3504 FORMAT(/,5X,'NODAL BOTTOM FRICTION VALUES ARE AVAILABLE', &
        /,6X,' IN UNIT 21 INPUT FILE')
    ENDIF

    RETURN
!     ----------------------------------------------------------------
    END SUBROUTINE ReadLegacyBottomFrictionFile
!     ----------------------------------------------------------------


!     ----------------------------------------------------------------
!        S U B R O U T I N E   N A _ T E R M I N A T E
!     ----------------------------------------------------------------
!     Provide clean termination when an error occurs.
!     ----------------------------------------------------------------
    SUBROUTINE na_terminate(NO_MPI_FINALIZE)
#ifdef CMPI
    USE MESSENGER
#endif
    USE GLOBAL, ONLY : setMessageSource, unsetMessageSource, allMessage, &
    DEBUG, ECHO, INFO, WARNING, ERROR
    IMPLICIT NONE
    LOGICAL, OPTIONAL :: NO_MPI_FINALIZE

    call setMessageSource("terminate")
#if defined(NODALATTR_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif

    call allMessage(INFO,"ADCIRC Terminating.")

#ifdef CMPI
    subdomainFatalError = .TRUE. 
    IF (PRESENT(NO_MPI_FINALIZE)) THEN
        CALL MSG_FINI(NO_MPI_FINALIZE)
    ELSE
        CALL MSG_FINI()
    ENDIF
#endif
    STOP

#if defined(NODALATTR_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
    call unsetMessageSource()
!     ----------------------------------------------------------------
    END SUBROUTINE na_terminate
!     ----------------------------------------------------------------


!-----------------------------------------------------------------------
    END MODULE NodalAttributes
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


