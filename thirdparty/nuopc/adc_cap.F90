!>
!! @mainpage ADCIRC NUOPC Cap
!! @author Saeed Moghimi (moghimis@gmail.com)
!! @date 15/1/17 Original documentation
!---------------------------------------------
!LOG-----------------
!
!
!

!!
!! @tableofcontents
!!
!! @section Overview Overview
!!
!! **This ADCIRC cap is under development with version 53.dev of ADCIRC.**
!!
!! This document describes the ADCIRC "cap", which is a small software layer that is
!! required when the [ADCIRC model] (http://adcirc.org/)
!! is used in [National Unified Operation Prediction Capability]
!!
!! The code and documentation is based on work by Fei Liu (fei.liu@gmail.com)
!! and Rocky Dunlap (rocky.dunlap@noaa.gov) [Link1]
!! (https://esgf.esrl.noaa.gov/projects/couplednems/ADCIRC_cap) and
!! [Link2](http://earthsystemmodeling.org/nuopc/caps/ADCIRC/) .

!! @subsection CapSubroutines Cap Subroutines
!!
!! The following table summarizes the NUOPC-required subroutines that appear in the
!! cap.  The "Phase" column says whether the subroutine is called during the
!! initialization, run, or finalize part of the coupled system run.
!!
!!
!! Phase    | ADCIRC Cap Subroutine                                                |  Description
!! ---------|--------------------------------------------------------------------|-------------------------------------------------------------
!! Init     | [InitializeP0] (@ref ADCIRC_cap_mod::initializep0)                   | Sets the Initialize Phase Definition (IPD) version to use
!! Init     | [InitializeAdvertise] (@ref ADCIRC_cap_mod::initializeadvertise)     | Advertises standard names of import and export fields
!! Init     | [InitializeRealize] (@ref ADCIRC_cap_mod::initializerealize)         | Creates an ESMF_Grid for the ADCIRC grid as well as ESMF_Fields for import and export fields
!! Init     | [SetClock] (@ref ADCIRC_cap_mod::setclock)                           | Before the run, sets up the timestep interval
!! Run      | [ModelAdvance_Slow] (@ref ADCIRC_cap_mod::modeladvance_slow)         | Advances the model by a timestep by calling ADCIRC_Run
!! Final    | [ADCIRC_Model_Finalize] (@ref ADCIRC_cap_mod::ADCIRC_model_finalize)     | Cleans up by calling ADCIRC_Finalize
!!
!! @subsection DomainCreation Domain Creation
!!
!! The ADCIRC cap currently only supports the triangular grid
!! for 2D variables.
!!
!! No changes were made to the native ADCIRC mechanisms for setting up the mesh.
!! This procedure happens as part of the call to `ADCIRC_Init()`.  After ADCIRC sets up
!! its mesh, the mesh structure is translated into an `ESMF_Mesh` object.  Setting
!! up the `ESMF_Mesh` is handled as part of the [RealizeFieldsProvidingGrid]
!! (@ref adc_cap::RealizeFieldsProvidingGrid) subroutine.  All of the details of the triangular
!! grid are be provided to ESMF via reading input file or by including representative variables.
!! SEE:
!! @code
!!    USE MESSENGER, ONLY  : MPI_COMM_ADCIRC
!!    USE ADCIRC_Mod, ONLY : ADCIRC_Init
!!    USE ADCIRC_Mod, ONLY : ADCIRC_Run
!!    USE ADCIRC_Mod, ONLY : ADCIRC_Final
!!    use MESH   , only: np,ne,nm,slam,sfea
!!    use GLOBAL , only: IMAP_EL_LG,NODES_LG
!!    use GLOBAL , only: ETA2, UU2, VV2
!!    use SIZES  , only: ROOTDIR
!! @endcode
!!
!! The ESMF_Mesh is set up in several steps inside [RealizeFieldsProvidingGrid]
!! (@ref adc_cap::RealizeFieldsProvidingGrid). A mesh data type
!! [`meshdata`](@ref adc_cap::meshdata)
!! wrote by [Ali Samii](https://github.com/samiiali) were used to read and create mesh
!! data structure.

!! @subsection Initialization Initialization
!!
!! The ADCIRC cap calls into the native `ADCIRC_Init(adc_comm)` subroutine in the
!! [RealizeFieldsProvidingGrid](@ref adc_cap::RealizeFieldsProvidingGrid)
!! subroutine.
!! The only parameter passed is the MPI communicator.  The global MPI communicator
!! is managed by ESMF and is retrieved from the current `ESMF_VM`.
!!
!! @subsection Run Run
!!
!! The time stepping needs to be setteled <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<.
!! The ADCIRC cap calls into the native `ADCIRC_Run(ZZZZZZZZZZZZ)` subroutine in the
!! [ModelAdvance] (@ref adc_cap::ModelAdvance) subroutine. The internal
!! ADCIRC timestepping loop has been >>>>>>> disabled or not <<<<<< when using
!! the NUOPC cap since the
!! main driver loop is provided by the NUOPC infrastructure at a level above
!! the ADCIRC model.  Therefore, a call into `ADCIRC_Run()` will only advance the
!! ADCIRC model by a single time step or number of time steps for one couling period.
!!
!! @subsection Finalization Finalization
!!
!! To finalize the model, the ADCIRC cap simple calls into the native `ADCIRC_Final()`
!! subroutine.  >>>>>>>>>>>>>>> NOT decided where and how !!!!!!!!!
!!
!!
!! @section ModelFields Model Fields
!!
!! The following tables list the import and export fields currently set up in the ADCIRC cap.
!! *** We assumed all the fields are 1D vectors both for import and export ***
!!
!! @subsection ExportField Export Fields
!!
!! Standard Name                        | Units      | Model Variable  | File         | Description             | Notes
!! -------------------------------------|------------|-----------------|--------------|-------------------------|-----------------
!! surface_eastward_sea_water_velocity  | m s-1      | UU2             | global.F     | current u component     |
!! surface_northward_sea_water_velocity | m s-1      | VV2             | global.F     | current v component     |
!! sea_surface_height_above_sea_level   | m          | ETA2            | global.F     | sea surface elevation   |
!!
!!
!! From CF convention page
!! sea_surface_height_above_sea_level
!!       alias: sea_surface_height
!!      sea_level means mean sea level, which is close to the geoid in sea areas.
!!      "Sea surface height" is a time-varying quantity. The standard name for the
!!       height of the sea surface above the geoid is sea_surface_height_above_geoid.
!!       The standard name for the height of the sea surface above the reference
!!       ellipsoid is sea_surface_height_above_reference_ellipsoid.
!!
!!
!! The following tables list the import and export fields currently set up in the ADCIRC cap.
!! @subsection ImportFields Import Fields
!!
!! Standard Name               |Short Name | Units                | Model Variable  | File         | Description                | Notes
!! ----------------------------|-----------|----------------------|-----------------|--------------|----------------------------|-----------------
!! eastward_radiation_stress   |sxx        | N.m^-2/rho ->  m2s-2 | ADCIRC_SXX      | global.F     |                            |
!! northward_radiation_stress  |syy        | N.m^-2/rho ->  m2s-2 | ADCIRC_SXY      | global.F     |                            |
!! cross_radiation_stress      |sxy        | N.m^-2/rho ->  m2s-2 | ADCIRC_SXY      | global.F     |                            |
!!
!!expFieldName         expFieldStdName
!!sxx                           eastward_radiation_stress
!!syy                           northward_radiation_stress
!!sxy                           cross_radiation_stress
!!ADCIRC accepts wave-driven stresses "in units of velocity squared
!!    (consistent with the units of gravity).  Stress in these units is obtained
!!    by dividing stress in units of force/area by the reference density of water."
!!    SO WE MUST DIVIDE BY RHO!
!!
!! TODO : make sure about the unit!!!!
!!         ADCIRC_SXX(I,2) = ADCIRC_SXX(I,2) / RHO
!!         ADCIRC_SXY(I,2) = ADCIRC_SXY(I,2) / RHO
!!         ADCIRC_SYY(I,2) = ADCIRC_SYY(I,2) / RHO
!-------------------------------------------------------
!!
!! @subsection MemoryManagement Memory Management
!!
!! For coupling fields, the ADCIRC cap has the capability to either reference the internal
!! ADCIRC data arrays directly, or to allow ESMF to allocate separate memory for the
!! coupling fields ????????????????????? Really   ??????????.  Currently, the ADCIRC
!! cap is set up not to ?????? reference internal data
!! arrays directly.  During the [XXXXXXXXXX] (@ref adc_cap::XXXXXXX)
!! phase, fields are copied from the import state into ADCIRC internal data arrays.
!! After the model integration timestep completes, internal data arrays are copied
!! into the export state so they can be transferred for coupling.
!!
!! @subsection IO I/O
!!
!! The ADCIRC cap implements a subroutine called `dumpADCIRCInternal` that writes
!! the fields from XXXXXX_flux.F90 to NetCDF files. <<<<<< NICE to USE >>>>>>>


module adc_cap

  !-----------------------------------------------------------------------------
  ! ADCIRC Component.
  !-----------------------------------------------------------------------------
  use mpi
  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS        => SetServices,    &
    model_label_SetClock    => label_SetClock, &
    model_label_CheckImport => label_CheckImport, &    
    model_label_Advance     => label_Advance,  &
    model_label_Finalize    => label_Finalize

  USE MESSENGER, ONLY  : MPI_COMM_ADCIRC,UPDATER

  USE ADCIRC_Mod, ONLY : ADCIRC_Init
  USE ADCIRC_Mod, ONLY : ADCIRC_Run
  USE ADCIRC_Mod, ONLY : ADCIRC_Final
!++GML added dp(bathymetric depth) 20210512
  use MESH   , only: np,ne,nm,slam,sfea,dp
  use GLOBAL , only: IMAP_EL_LG,NODES_LG

  !TODO::  Need to carefully check units in particular pressure
  use GLOBAL , only: ETA2, UU2, VV2      ! Export water level and velocity fileds to wave model
  USE GLOBAL,  ONLY: RSNX1, RSNY1        ! Import wave 2D forces from wave model
  USE GLOBAL,  ONLY: RSNX2, RSNY2        ! Import wave 2D forces from wave model
 ! USE GLOBAL,  ONLY: WVNX2, WVNY2, PRN2  ! Import wind and pressure variables
 ! USE GLOBAL,  ONLY: WVNX1, WVNY1, PRN1
  USE WIND,  ONLY: WVNX2, WVNY2, PRN2  ! Import wind and pressure variables
  USE WIND,  ONLY: WVNX1, WVNY1, PRN1
 ! USE GLOBAL,  ONLY: WTIMINC             ! wind time interval  may be set in ATM.cap or ........  <<:TODO:
  USE WIND,  ONLY: WTIMINC             ! wind time interval  may be set in ATM.cap or ........  <<:TODO:
  USE GLOBAL,  ONLY: RSTIMINC            ! wave time interval
  use GLOBAL,  ONLY: RhoWat0, NUOPC4MET, NUOPC4WAV, NWS, g
  use GLOBAL,  only: ITHS, NT, DTDP, ITIME
  use GLOBAL,  only: allMessage

  use SIZES  , only: ROOTDIR
  use couple2swan_modif, only: ADCIRC_SXX, ADCIRC_SXY, ADCIRC_SYY
  use couple2swan_modif, only: ComputeWaveDrivenForces, InterpoWeight

  use adc_mod, only: meshdata
  use adc_mod, only: create_parallel_esmf_mesh_from_meshdata
  use adc_mod, only: extract_parallel_data_from_mesh_orig
  !use adc_mod, only: adc_cpl_int,adc_cpl_num,adc_cpl_den  !time info is getting from driver
  !use adc_mod, only: read_config
  
  implicit none

  private
  public SetServices

  type fld_list_type
    character(len=64) :: stdname
    character(len=64) :: shortname
    character(len=64) :: unit
    logical           :: assoc    ! is the farrayPtr associated with internal data
    logical           :: connected
    real(ESMF_KIND_R8), dimension(:), pointer :: farrayPtr
  end type fld_list_type

  integer,parameter :: fldsMax = 100
  integer :: fldsToAdc_num = 0
  type (fld_list_type) :: fldsToAdc(fldsMax)
  integer :: fldsFrAdc_num = 0
  type (fld_list_type) :: fldsFrAdc(fldsMax)

  type(meshdata),save  :: mdataIn, mdataOut
  
  character(len=2048):: info
  integer :: dbrc     ! temporary debug rc value

  !to test feild halo update.
  type (ESMF_RouteHandle), save :: ATM_HaloRouteHandel
  
  !real(ESMF_KIND_R8)      :: WaveCouplingIntervalSec, WindCouplingIntervalSec  !in seconds
  !type(ESMF_TimeInterval) :: WaveCouplingInterval, WindCouplingInterval

  logical, save            :: first_exchange = .true.
  logical, save            :: first_exchange1 = .true.
  integer, save            :: iunit_log = 10000


!  real,parameter :: wave_force_limmit = 0.05
!++ GML  0.01 work for AKUT mesh
  real,parameter :: wave_force_limmit = 0.001
!++
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  !> NUOPC SetService method is the only public entry point.
  !! SetServices registers all of the user-provided subroutines
  !! in the module with the NUOPC layer.
  !!
  !! @param model an ESMF_GridComp object
  !! @param rc return code
  subroutine SetServices(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    character(len=*),parameter   :: subname='(adc_cap:SetServices)'
    ! Local variables
    integer                      :: num,i
    rc = ESMF_SUCCESS
    
    ! read config file
    !call read_config()   !information will be read 
                          !from nems.configure by NEMS driver as time slots  
    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(model, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeP2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    !Assume no need to change clock settings
    ! attach specializing method(s)
!    call NUOPC_CompSpecialize(model, specLabel=model_label_SetClock, &
!      specRoutine=SetClock, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    !comment out for now to avoid over writing NOUPC check import
!    call NUOPC_CompSpecialize(model, specLabel=model_label_CheckImport, &
!      specPhaseLabel="RunPhase1", specRoutine=CheckImport, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_CompSpecialize(model, specLabel=model_label_Finalize, &
      specRoutine=ADCIRC_model_finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ADCIRC_FieldsSetup()
!
    do num = 1,fldsToAdc_num
        !print *,  "fldsFrAdc_num  ", fldsToAdc_num
        !print *,  "fldsToAdc(num)%shortname  ", fldsToAdc(num)%shortname
        !print *,  "fldsToAdc(num)%stdname  ", fldsToAdc(num)%stdname
     write(info,*) subname,"fldsToAdc(num)%stdname  ", fldsToAdc(num)%stdname
   end do
!
   do num = 1,fldsFrAdc_num
     !print *,  "fldsFrAdc_num  ", fldsFrAdc_num
     !print *,  "fldsFrAdc(num)%shortname  ", fldsFrAdc(num)%shortname
     !print *,  "fldsFrAdc(num)%stdname  ", fldsFrAdc(num)%stdname
     write(info,*) subname,"fldsFrAdc(num)%stdname  ", fldsFrAdc(num)%stdname
   end do
!
    
    write(info,*) subname,' --- adc SetServices completed --- '
    !print *,      subname,' --- adc SetServices completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)
  end subroutine
  
  !-----------------------------------------------------------------------------
    !> First initialize subroutine called by NUOPC.  The purpose
    !! is to set which version of the Initialize Phase Definition (IPD)
    !! to use.
    !!
    !! For this ADCIRC cap, we are using IPDv01.
    !!
    !! @param model an ESMF_GridComp object
    !! @param importState an ESMF_State object for import fields
    !! @param exportState an ESMF_State object for export fields
    !! @param clock an ESMF_Clock object
    !! @param rc return code

   subroutine InitializeP1(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Local Variables
    type(ESMF_VM)                :: vm
    integer                      :: esmf_comm,adc_comm,ierr
    character(len=*),parameter   :: subname='(adc_cap:AdvertiseFields)'

    rc = ESMF_SUCCESS

    ! details Get current ESMF VM.
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! details Get MPI_communicator from ESMF VM.
    call ESMF_VMGet(vm, mpiCommunicator=esmf_comm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call MPI_Comm_dup(esmf_comm, adc_comm, ierr)
    ! Duplicate the MPI communicator not to interfere with ESMF communications.
    ! The duplicate MPI communicator can be used in any MPI call in the user
    ! code. Here the MPI_Barrier() routine is called.
    call MPI_Barrier(adc_comm, ierr)
    !Initialize adcirc before setting up fields
    
    !NUOPC4MET = .true.
    !NUOPC4WAV = .true.
        
    CALL ADCIRC_Init(adc_comm)

    !WTIMINC = 3*3600.0
    !RSTIMINC = adc_cpl_int +  adc_cpl_num / adc_cpl_den
    !print *,   'WTIMINC > ',WTIMINC,'  RSTIMINC > ',RSTIMINC
!
    call ADCIRC_AdvertiseFields(importState, fldsToAdc_num, fldsToAdc, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ADCIRC_AdvertiseFields(exportState, fldsFrAdc_num, fldsFrAdc, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

        write(info,*) subname,' --- initialization phase 1 completed --- '
        !print *,      subname,' --- initialization phase 1 completed --- '
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
!
  end subroutine


  !> Advertises a set of fields in an ESMF_State object by calling
  !! NUOPC_Advertise in a loop.
  !!
  !! @param state the ESMF_State object in which to advertise the fields
  !! @param nfields number of fields
  !! @param field_defs an array of fld_list_type listing the fields to advertise
  !! @param rc return code
  subroutine ADCIRC_AdvertiseFields(state, nfields, field_defs, rc)
    type(ESMF_State), intent(inout)             :: state
    integer,intent(in)                          :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    integer, intent(inout)                      :: rc

    integer                                     :: i
    character(len=*),parameter  :: subname='(adc_cap:ADCIRC_AdvertiseFields)'

    rc = ESMF_SUCCESS

    do i = 1, nfields
      !print *, 'Advertise: '//trim(field_defs(i)%stdname)//'---'//trim(field_defs(i)%shortname)
      call ESMF_LogWrite('Advertise: '//trim(field_defs(i)%stdname), ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call NUOPC_Advertise(state, &
        standardName=field_defs(i)%stdname, &
        name=field_defs(i)%shortname, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    enddo
    !print *,      subname,' --- IN   --- '

  end subroutine ADCIRC_AdvertiseFields
  !
  !----------------------------------------------------------------------------------
  subroutine ADCIRC_FieldsSetup
    integer                     :: rc
    character(len=*),parameter  :: subname='(adc_cap:ADCIRC_FieldsSetup)'

    !--------- import fields to Sea Adc -------------
    !TODO: Consider moving these lines to driver to avoid doing it in both CAPS
!    call NUOPC_FieldDictionaryAddEntry("eastward_radiation_stress",  "mx", rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!        line=__LINE__, &
!        file=__FILE__)) &
!        return  ! bail out

!    call NUOPC_FieldDictionaryAddEntry("northward_radiation_stress", "mx", rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!        line=__LINE__, &
!        file=__FILE__)) &
!        return  ! bail out

!    call NUOPC_FieldDictionaryAddEntry("cross_radiation_stress",    "mx", rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!        line=__LINE__, &
!        file=__FILE__)) &
!        return  ! bail out


     !--- kf fixes
     ! sxx
     if (.not.NUOPC_FieldDictionaryHasEntry( &
                  "eastward_wave_radiation_stress")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_wave_radiation_stress", &
          canonicalUnits="N m-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      ! sxy
      if (.not.NUOPC_FieldDictionaryHasEntry( &
                  "eastward_northward_wave_radiation_stress")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_northward_wave_radiation_stress", &
          canonicalUnits="N m-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      ! syy
      if (.not.NUOPC_FieldDictionaryHasEntry( &
                   "northward_wave_radiation_stress")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="northward_wave_radiation_stress", &
          canonicalUnits="N m-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

!   !------ immport field from ww3 to adc
!    call fld_list_add(num=fldsToAdc_num, fldlist=fldsToAdc, stdname="eastward_radiation_stress", shortname= "sxx")
!    call fld_list_add(num=fldsToAdc_num, fldlist=fldsToAdc, stdname="northward_radiation_stress",shortname= "syy")
!    call fld_list_add(num=fldsToAdc_num, fldlist=fldsToAdc, stdname="cross_radiation_stress",    shortname= "sxy")
    !c- kf fixes
    call fld_list_add(num=fldsToAdc_num, fldlist=fldsToAdc, stdname="eastward_wave_radiation_stress", shortname= "sxx")
    call fld_list_add(num=fldsToAdc_num, fldlist=fldsToAdc, stdname="northward_wave_radiation_stress",shortname= "syy")
    call fld_list_add(num=fldsToAdc_num, fldlist=fldsToAdc, stdname="eastward_northward_wave_radiation_stress",shortname= "sxy")

    !--------- import fields from atm to Adc -------------
    call fld_list_add(num=fldsToAdc_num, fldlist=fldsToAdc, stdname= "air_pressure_at_sea_level", shortname= "pmsl" )
    call fld_list_add(num=fldsToAdc_num, fldlist=fldsToAdc, stdname= "inst_merid_wind_height10m", shortname= "imwh10m" )
    call fld_list_add(num=fldsToAdc_num, fldlist=fldsToAdc, stdname="inst_zonal_wind_height10m" , shortname= "izwh10m" )
    !--------- export fields from Sea Adc -------------
    call fld_list_add(num=fldsFrAdc_num, fldlist=fldsFrAdc, stdname="sea_surface_height_above_sea_level",  shortname= "zeta" )
    call fld_list_add(num=fldsFrAdc_num, fldlist=fldsFrAdc, stdname="surface_eastward_sea_water_velocity", shortname= "velx" )
    call fld_list_add(num=fldsFrAdc_num, fldlist=fldsFrAdc, stdname="surface_northward_sea_water_velocity",shortname= "vely" )

!NEMS hycod standard names>>>
!https://esgf.esrl.noaa.gov/projects/couplednems/coupling_fields
!ocn_current_zonal
!ocncz
!m s-1	Ocean current X component.	 	 	 	
!ocn_current_merid
!ocncm
!m s-1	Ocean current Y component.	 	



   !
    write(info,*) subname,' --- Passed--- '
    !print *,      subname,' --- Passed --- '
  end subroutine ADCIRC_FieldsSetup

    !---------------------------------------------------------------------------------

  subroutine fld_list_add(num, fldlist, stdname, data, shortname, unit)
    ! ----------------------------------------------
    ! Set up a list of field information
    ! ----------------------------------------------
    integer,             intent(inout)  :: num
    type(fld_list_type), intent(inout)  :: fldlist(:)
    character(len=*),    intent(in)     :: stdname
    real(ESMF_KIND_R8), dimension(:), optional, target :: data
    character(len=*),    intent(in),optional :: shortname
    character(len=*),    intent(in),optional :: unit

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(adc_cap:fld_list_add)'

    ! fill in the new entry

    num = num + 1
    if (num > fldsMax) then
      call ESMF_LogWrite(trim(subname)//": ERROR num gt fldsMax "//trim(stdname), &
        ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=rc)
      return
    endif

    fldlist(num)%stdname        = trim(stdname)
    if (present(shortname)) then
       fldlist(num)%shortname   = trim(shortname)
    else
       fldlist(num)%shortname   = trim(stdname)
    endif

    if (present(data)) then
      fldlist(num)%assoc        = .true.
      fldlist(num)%farrayPtr    => data
    else
      fldlist(num)%assoc        = .false.
    endif

    if (present(unit)) then
       fldlist(num)%unit        = unit
    endif


    write(info,*) subname,' --- Passed--- '
    !print *,      subname,' --- Passed --- '
  end subroutine fld_list_add

  !-----------------------------------------------------------------------------
    !> Called by NUOPC to realize import and export fields.

    !! The fields to import and export are stored in the fldsToAdc and fldsFrAdc
    !! arrays, respectively.  Each field entry includes the standard name,
    !! information about whether the field's grid will be provided by the cap,
    !! and optionally a pointer to the field's data array.  Currently, all fields
    !! are defined on the same mesh defined by the cap.
    !! The fields are created by calling adc_cap::adcirc_XXXXXXXXXXXXXXXXXXX.
    !!
    !! @param model an ESMF_GridComp object
    !! @param importState an ESMF_State object for import fields
    !! @param exportState an ESMF_State object for export fields
    !! @param clock an ESMF_Clock object
    !! @param rc return code
!
  subroutine InitializeP2(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables    
    type(ESMF_Field)        :: field
    !Saeed added
    type(meshdata)               :: mdata
    type(ESMF_Mesh)              :: ModelMesh,meshIn,meshOut
    type(ESMF_VM)                :: vm
    integer                      :: localPet, petCount
    character(len=*),parameter   :: subname='(adc_cap:RealizeFieldsProvidingGrid)'

    rc = ESMF_SUCCESS

    !> \details Get current ESMF VM.
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Get query local pet information for handeling global node information
    call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    ! call ESMF_VMPrint(vm, rc=rc)

    !! Assign VM to mesh data type.
    mdata%vm = vm
    !print *,localPet,"< LOCAL pet, ADC ..1.............................................. >> "
    ! create a Mesh object for Fields
    !call extract_parallel_data_from_mesh(ROOTDIR, mdata, localPet)
    call extract_parallel_data_from_mesh_orig(ROOTDIR, mdata, localPet)
    !    print *,"ADC ..2.............................................. >> "
    call create_parallel_esmf_mesh_from_meshdata(mdata,ModelMesh )
    !    print *,"ADC ..3.............................................. >> "
    !
    
    if (.false.) then
        call ESMF_MeshWrite(ModelMesh, filename="adc_mesh.nc", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
    end if

    !print *,"ADC >> "
    !print *,"NumNd", mdata%NumNd
    !print *,"NumOwnedNd", mdata%NumOwnedNd
    !print *,"NumEl", mdata%NumEl
    !print *,"NumND_per_El", mdata%NumND_per_El


    meshIn  = ModelMesh ! for now out same as in
    meshOut = meshIn

    mdataIn  = mdata
    mdataOut = mdata

    !print *,"ADC ..4.............................................. >> "
    call ADCIRC_RealizeFields(importState, meshIn , mdata, fldsToAdc_num, fldsToAdc, "Adcirc import", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    !
    call ADCIRC_RealizeFields(exportState, meshOut, mdata, fldsFrAdc_num, fldsFrAdc, "Adcirc export", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(info,*) subname,' --- initialization phase 2 completed --- '
    !print *,      subname,' --- initialization phase 2 completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)
  end subroutine
  

 !> Adds a set of fields to an ESMF_State object.  Each field is wrapped
  !! in an ESMF_Field object.  Memory is either allocated by ESMF or
  !! an existing ADCIRC pointer is referenced.
  !!
  !! @param state the ESMF_State object to add fields to
  !! @param grid the ESMF_Grid object on which to define the fields
  !! @param nfields number of fields
  !! @param field_defs array of fld_list_type indicating the fields to add
  !! @param tag used to output to the log
  !! @param rc return code
  subroutine ADCIRC_RealizeFields(state, mesh, mdata, nfields, field_defs, tag, rc)

    type(ESMF_State), intent(inout)             :: state
    type(ESMF_Mesh), intent(in)                 :: mesh
    type(meshdata)                              :: mdata
    integer, intent(in)                         :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    character(len=*), intent(in)                :: tag
    integer, intent(inout)                      :: rc


    type(ESMF_Field)                            :: field
    integer                                     :: i
    character(len=*),parameter  :: subname='(adc_cap:ADCIRC_RealizeFields)'

    rc = ESMF_SUCCESS

    do i = 1, nfields
        field = ESMF_FieldCreate(name=field_defs(i)%shortname, mesh=mesh, &
          typekind=ESMF_TYPEKIND_R8, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        if (NUOPC_IsConnected(state, fieldName=field_defs(i)%shortname)) then

            call NUOPC_Realize(state, field=field, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail out

            call ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is connected.", &
              ESMF_LOGMSG_INFO, &
              line=__LINE__, &
              file=__FILE__, &
              rc=dbrc)

            !print *,      subname,' --- Connected --- '
            field_defs(i)%connected = .true.
        else
            call ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is not connected.", &
              ESMF_LOGMSG_INFO, &
              line=__LINE__, &
              file=__FILE__, &
              rc=dbrc)
            ! TODO: Initialize the value in the pointer to 0 after proper restart is setup
            !if(associated(field_defs(i)%farrayPtr) ) field_defs(i)%farrayPtr = 0.0
            ! remove a not connected Field from State
            call ESMF_StateRemove(state, (/field_defs(i)%shortname/), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail out
            !print *,      subname," Field ", field_defs(i)%stdname ,' --- Not-Connected --- '
            field_defs(i)%connected = .false.

        endif
    enddo

        write(info,*) subname,' --- OUT--- '
        !print *,      subname,' --- OUT --- '
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=rc)
  end subroutine ADCIRC_RealizeFields
  !-----------------------------------------------------------------------------

!  subroutine SetClock_mine_not_active(model, rc)
!    type(ESMF_GridComp)  :: model
!    integer, intent(out) :: rc

!    ! local variables
!    type(ESMF_Clock)              :: clock
!    type(ESMF_TimeInterval)       :: ADCTimeStep

!    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
!    call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    ! initialize internal clock
    ! - on entry, the component clock is a copy of the parent clock
    ! - the parent clock is on the slow timescale atm timesteps
    ! - reset the component clock to have a timeStep that is for adc-wav of the parent
    !   -> timesteps
    
    !call ESMF_TimeIntervalSet(ADCTimeStep, s=     adc_cpl_int, sN=adc_cpl_num, sD=adc_cpl_den, rc=rc) ! 5 minute steps
    !TODO: use nint !!?
!    call ESMF_TimeIntervalSet(ADCTimeStep, s= adc_cpl_int , rc=rc) ! 5 minute steps
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

!    call NUOPC_CompSetClock(model, clock, ADCTimeStep, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
      
!    print *, "ADC Timeinterval1 = "
!    call ESMF_TimeIntervalPrint(ADCTimeStep, options="string", rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!        line=__LINE__, &
!        file=__FILE__)) &
!        return  ! bail out    

!  end subroutine
  !-----------------------------------------------------------------------------
  ! From ADCIRC model uses same clock as parent gridComp
!  subroutine SetClock_not_active(model, rc)
!    type(ESMF_GridComp)  :: model
!    integer, intent(out) :: rc
!    
!    ! local variables
!    type(ESMF_Clock)              :: clock
!    type(ESMF_TimeInterval)       :: ADCTimeStep, timestep
!    character(len=*),parameter  :: subname='(adc_cap:SetClock)'

!    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
!    call ESMF_GridCompGet(model, clock=clock, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    !call ESMF_TimeIntervalSet(ADCTimeStep, s=     adc_cpl_int, sN=adc_cpl_num, sD=adc_cpl_den, rc=rc) ! 5 minute steps
    ! tcraig: dt is the ADCIRC thermodynamic timestep in seconds
!   call ESMF_TimeIntervalSet(timestep, s=adc_cpl_int, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

!    call ESMF_ClockSet(clock, timestep=timestep, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
      
    ! initialize internal clock
    ! here: parent Clock and stability timeStep determine actual model timeStep
!    call ESMF_TimeIntervalSet(ADCTimeStep, s=adc_cpl_int, rc=rc) 
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
!    call NUOPC_CompSetClock(model, clock, ADCTimeStep, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
    
!  end subroutine


  !-----------------------------------------------------------------------------

  !> Called by NUOPC to advance the ADCIRC model a single timestep >>>>>
  !! <<<<<<<<<  TODO: check! this is not what we want!!!.
  !!
  !! This subroutine copies field data out of the cap import state and into the
  !! model internal arrays.  Then it calls ADCIRC_Run to make a NN timesteps.
  !! Finally, it copies the updated arrays into the cap export state.
  !!
  !! @param model an ESMF_GridComp object
  !! @param rc return code
  !-----------------------------------------------------------------------------
  subroutine ModelAdvance(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)           :: clock
    type(ESMF_State)           :: importState, exportState
    type(ESMF_Time)            :: currTime
    type(ESMF_TimeInterval)    :: timeStep
    character(len=*),parameter :: subname='(adc_cap:ModelAdvance)'
    !tmp vector
    real(ESMF_KIND_R8), pointer:: tmp(:)

    ! exports
    real(ESMF_KIND_R8), pointer:: dataPtr_zeta(:)
    real(ESMF_KIND_R8), pointer:: dataPtr_velx(:)
    real(ESMF_KIND_R8), pointer:: dataPtr_vely(:)

    !imports
    real(ESMF_KIND_R8), pointer:: dataPtr_sxx(:)
    real(ESMF_KIND_R8), pointer:: dataPtr_syy(:)
    real(ESMF_KIND_R8), pointer:: dataPtr_sxy(:)

    real(ESMF_KIND_R8), pointer:: dataPtr_pmsl(:)
    real(ESMF_KIND_R8), pointer:: dataPtr_imwh10m(:)
    real(ESMF_KIND_R8), pointer:: dataPtr_izwh10m(:)


    type(ESMF_StateItem_Flag)  :: itemType
    type(ESMF_Mesh)            :: mesh
    type(ESMF_Field)           :: lfield
    character(len=128)         :: fldname,timeStr
    integer                    :: i1,num
    integer                    :: ITIME_BGN_ADC, ITIME_END_ADC
    integer                    :: nCplADC
    real(ESMF_KIND_R8)         :: timeStepAbs
    ! local variables for Get methods
    integer :: YY, MM, DD, H, M, S
    integer :: ss,ssN,ssD
    logical :: wave_forcing, meteo_forcing, surge_forcing

    type(ESMF_Time) :: BeforeCaribbeanTime,AfterCaribbeanTime

    ! DW
    INTEGER, save:: ienter = 0 ;
!++ GML
    double precision :: xiplus, xineg, xi

    xi = HUGE(xi)
    xiplus  = 2 * xi     ! yields +Infinity
    xineg   =-2 * xi     ! yields -Infinity
!++
    rc = ESMF_SUCCESS
    dbrc = ESMF_SUCCESS
    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    ! Because of the way that the internal Clock was set in SetClock(),
    ! its timeStep is likely smaller than the parent timeStep. As a consequence
    ! the time interval covered by a single parent timeStep will result in
    ! multiple calls to the ModelAdvance() routine. Every time the currTime
    ! will come in by one internal timeStep advanced. This goes until the
    ! stopTime of the internal Clock has been reached.

    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing ADC from: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimePrint(currTime + timeStep, &
      preString="-------------------------------> to: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimeGet(currTime, yy=YY, mm=MM, dd=DD, h=H, m=M, s=S, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    write(info, *)  "ADC currTime = ", YY, "/", MM, "/", DD," ", H, ":", M, ":", S
!   call allMessage(1,info)  ! Dec 2020

    call ESMF_TimeGet(currTime, timeStringISOFrac=timeStr , rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    !global values for adcirc time step number
    ITIME_BGN_ADC = ITHS + 1   !if hot start then ITHS>0
    ITIME_END_ADC = NT         !NT is set in read_input.F

    call ESMF_TimeIntervalGet(timeStep, s=ss,sN=ssN,sD=ssD,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    timeStepAbs = real(ss) + real(ssN/ssD)
    if (mod(timeStepAbs,DTDP) .EQ. 0) then
      nCplADC = nint(timeStepAbs/DTDP)
    else
      nCplADC = nint(timeStepAbs/DTDP)+1
    endif


     !print *, '  nCplADC = ', nCplADC
     !print *, '  ADCCouplingTimeInterval = ', timeStepAbs
    !-----------------------------------------
    !   IMPORT
    !-----------------------------------------
    !Get and fill imported fields
    ! <<<<< RECEIVE and UN-PACK SXX
   wave_forcing= .true.
   
   do num = 1,fldsToAdc_num
      if (fldsToAdc(num)%shortname == 'sxx') wave_forcing = wave_forcing .and. fldsToAdc(num)%connected
      if (fldsToAdc(num)%shortname == 'syy') wave_forcing = wave_forcing .and. fldsToAdc(num)%connected
      if (fldsToAdc(num)%shortname == 'sxy') wave_forcing = wave_forcing .and. fldsToAdc(num)%connected
   end do

    if (wave_forcing) then

        ! Wave time step
        RSTIMINC = nint(timeStepAbs)  !TODO: Get it from coupler based on the time slots.
                                      !TODO: This implemetation wotks for one time slot for wind and wave right now.
                                      !TODO: 
        !RSTIMINC = adc_cpl_int +  adc_cpl_num / adc_cpl_den
        !print *, ' in cap   ....> RSTIMINC > ', RSTIMINC
        call State_getFldPtr_(ST=importState,fldname='sxx',fldptr=dataPtr_sxx, &
          rc=rc,dump=.true.,timeStr=timeStr)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        !-----------------------------------------
        ! <<<<< RECEIVE and UN-PACK SYY
        call State_getFldPtr(ST=importState,fldname='syy',fldptr=dataPtr_syy,rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        !-----------------------------------------
        ! <<<<< RECEIVE and UN-PACK SXY
        call State_getFldPtr(ST=importState,fldname='sxy',fldptr=dataPtr_sxy,rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        !print *, 'size dataPtr_sxy > ', size(dataPtr_sxy)
        !print *, 'in cap maxval(RSNX2)', maxval(RSNX2)
        !print *, 'in cap maxval(RSNY2)', maxval(RSNY2)
        ! Allocate arrays for radiation stresses.
        IF(.NOT.ALLOCATED(ADCIRC_SXX)) ALLOCATE(ADCIRC_SXX(1:NP,1:2))
        IF(.NOT.ALLOCATED(ADCIRC_SXY)) ALLOCATE(ADCIRC_SXY(1:NP,1:2))
        IF(.NOT.ALLOCATED(ADCIRC_SYY)) ALLOCATE(ADCIRC_SYY(1:NP,1:2))

        ! Allocate arrays for radiation stresses.
        IF(.NOT.ALLOCATED(RSNX1)) ALLOCATE(RSNX1(1:NP))
        IF(.NOT.ALLOCATED(RSNY1)) ALLOCATE(RSNY1(1:NP))

        IF(.NOT.ALLOCATED(RSNX2)) ALLOCATE(RSNX2(1:NP))
        IF(.NOT.ALLOCATED(RSNY2)) ALLOCATE(RSNY2(1:NP))

        ! updaye last rad str
        ADCIRC_SXX(:,1) = ADCIRC_SXX(:,2)
        ADCIRC_SYY(:,1) = ADCIRC_SYY(:,2)
        ADCIRC_SXY(:,1) = ADCIRC_SXY(:,2)

        ! Fill owned nodes from imported data to model variable
        ! devide by water density to convert from N.m-2 to m2s-2
        do i1 = 1, mdataOut%NumOwnedNd, 1
            ADCIRC_SXX(mdataOut%owned_to_present_nodes(i1),1) = dataPtr_sxx(i1)
        end do
        ADCIRC_SXX(:,2) = ADCIRC_SXX(:,1)

        do i1 = 1, mdataOut%NumOwnedNd, 1
            ADCIRC_SYY(mdataOut%owned_to_present_nodes(i1),1) = dataPtr_syy(i1)
        end do
        ADCIRC_SYY(:,2) = ADCIRC_SYY(:,1)

        do i1 = 1, mdataOut%NumOwnedNd, 1
            ADCIRC_SXY(mdataOut%owned_to_present_nodes(i1),1) = dataPtr_sxy(i1)
        end do
        ADCIRC_SXY(:,2) = ADCIRC_SXY(:,1)

        ! Ghost nodes update before calculating the wave forces
        call UPDATER( ADCIRC_SXX(:,1), ADCIRC_SYY(:,1), ADCIRC_SXY(:,1),3)
        call UPDATER( ADCIRC_SXX(:,2), ADCIRC_SYY(:,2), ADCIRC_SXY(:,2),3)

        !print *, 'Hard Coded >>>>>  SXX >>>>>> '
        !mask
        where(abs(ADCIRC_SXX).gt. 1e6)  ADCIRC_SXX =  1e-10
        where(abs(ADCIRC_SYY).gt. 1e6)  ADCIRC_SYY =  1e-10
        where(abs(ADCIRC_SXY).gt. 1e6)  ADCIRC_SXY =  1e-10
        !max values
        !where((ADCIRC_SXX).gt. 1e3)  ADCIRC_SXX =  1e3
        !where((ADCIRC_SYY).gt. 1e3)  ADCIRC_SYY =  1e3
        !where((ADCIRC_SXY).gt. 1e3)  ADCIRC_SXY =  1e3
        !where((ADCIRC_SXX).lt.-1e3)  ADCIRC_SXX = -1e3
        !where((ADCIRC_SYY).lt.-1e3)  ADCIRC_SYY = -1e3
        !where((ADCIRC_SXY).lt.-1e3)  ADCIRC_SXY = -1e3




!    NOTE: ADCIRC accepts wave-driven stresses "in units of velocity squared
!    (consistent with the units of gravity).  Stress in these units is obtained
!    by dividing stress in units of force/area by the reference density of water."


!
!
!From WW3 source : w3iogomd.ftn  line: 1756
!constants.ftn:34:!      DWAT      Real  Global   Density of water               (kg/m3)
!constants.ftn:63:      REAL, PARAMETER         :: DWAT   = 1000.
!
!w3iogomd.ftn:1753:      SXX    = SXX * DWAT * GRAV
!w3iogomd.ftn:1754:      SYY    = SYY * DWAT * GRAV
!w3iogomd.ftn:1755:      SXY    = SXY * DWAT * GRAV
!
!    Rad.Str info from netcdf header
!		 sxx:long_name = "Radiation stress component Sxx" ;
!		 sxx:standard_name = "radiation_stress_component_sxx" ;
!		 sxx:globwave_name = "significant_wave_height" ;
!		 sxx:units = "N m-1" ;
!		 sxx:_FillValue = 9.96921e+36f ;
!		 sxx:scale_factor = 1.f ;
!		 sxx:add_offset = 0.f ;
!		 sxx:valid_min = -3000 ;
!	 	 sxx:valid_max = 3000 ;
 
!    Therefore we need to divide sxx/rho to change its unit to m3s-2
!    in force calculation we do d(sxx)/dx therefore the final force 
!    unit will be m2s-2 which is the correct one.   

        ADCIRC_SXX = ADCIRC_SXX / RhoWat0
        ADCIRC_SYY = ADCIRC_SYY / RhoWat0
        ADCIRC_SXY = ADCIRC_SXY / RhoWat0
        
        !print *, 'in cap maxval(ADCIRC_SXX)', maxval(ADCIRC_SXX)

        InterpoWeight = 0.0  !avoid time interpolation

        RSNX1 = RSNX2
        RSNY1 = RSNY2

        ! Calculate wave forces
        call ComputeWaveDrivenForces

        !iunit_log = iunit_log + 1
        !open(unit = iunit_log, ACTION = "write", STATUS ="replace" )
        !write(iunit_log, *) RSNX2, RSNY2
        !close(iunit_log)
!++ GML 20210308 exclude infinity doesn't work yet
!        if ( first_exchange1 ) then
!          RSNX2 = 0.0d0
!          RSNY2 = 0.0d0
!          first_exchange1 = .false.
!        end if  
!        where(RSNX2 .eq. xiplus) RSNX2 = 0.0d0
!        where(RSNX2 .eq. xineg) RSNX2 = 0.0d0
!        where(RSNY2 .eq. xiplus) RSNY2 = 0.0d0
!        where(RSNY2 .eq. xineg) RSNY2 = 0.0d0
!!
!        do i1 = 1,  NP
!           if(RSNX2(i1) == xineg)then
!           write(*,*) 'check -infinity'
!           write(*,*) i1, RSNX2(i1), xineg
!           RSNX2(i1) = 0.0d0
!           endif
!           if(RSNY2(i1) == xineg)then
!           RSNY2(i1) = 0.0d0
!           endif
!           if(RSNX2(i1) == xiplus)then
!           write(*,*) 'check infinity'
!           write(*,*) i1, RSNX2(i1), xiplus
!           RSNX2(i1) = 0.0d0
!           endif
!           if(RSNY2(i1) == xiplus)then
!           RSNY2(i1) = 0.0d0
!           endif
!        enddo
!++
        write(info,*) subname,' --- wave data exchange OK / wave feilds are all connected --- / Model advances '
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
        !print *, info
!++ GML 20210308
        !RSNX2 = 0.0001
        !RSNY2 = 0.0001
!        if(RSNX2 - 1 == RSNX2)then
!        write(*,*) 'check RSNX2=', RSNX2
!        RSNX2 = 0.0
!        endif
!        if(RSNY2 - 1 == RSNY2)then
!        write(*,*) 'check RSNY2=', RSNY2
!        RSNY2 = 0.0
!        endif
!        where(abs(RSNX2).gt. wave_force_limmit) RSNX2 = 0.0d0
!        where(abs(RSNY2).gt. wave_force_limmit) RSNY2 = 0.0d0
!        RSNX2 = 0.0
!        RSNY2 = 0.0
!            where(RSNX2.gt. wave_force_limmit) RSNX2 =  wave_force_limmit
!            where(RSNY2.gt. wave_force_limmit) RSNY2 =  wave_force_limmit
!
!            where(RSNX2.le. (-1.0 * wave_force_limmit)) RSNX2 =  wave_force_limmit
!            where(RSNY2.le. (-1.0 * wave_force_limmit)) RSNY2 =  wave_force_limmit
!++
        ! initialize time reach to Caribean Islands
        !call ESMF_TimeSet(BeforeCaribbeanTime, yy=2008, mm=9, dd=6 , h=12, m=0, s=0, rc=rc)
        !call ESMF_TimeSet(AfterCaribbeanTime , yy=2008, mm=9, dd=12, h=12, m=0, s=0, rc=rc)

!!!!#ifdef NO_COMPILE00000
        !-----------------------
   !     if ((currTime > BeforeCaribbeanTime) .and. (currTime < AfterCaribbeanTime)) then
            write(info,*) subname, 'in cap after maxval(RSNX2)', maxval(RSNX2)
            call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

            write(info,*) subname, 'in cap after maxval(RSNY2)', maxval(RSNY2)
            call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
! GML
            !print *, 'Hard Coded >>>>>>>>>>>  where(abs(RSNX2).gt. wave_force_limmit) RSNX2 =  wave_force_limmit'
            !write(info,*) subname,'Hard Coded >>>>>>>>>>>  where(abs(RSNX2).gt. wave_force_limmit) RSNX2 =  wave_force_limmit'
            !call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

            !where(RSNX2.gt. wave_force_limmit) RSNX2 =  wave_force_limmit
            !where(RSNY2.gt. wave_force_limmit) RSNY2 =  wave_force_limmit

            !where(RSNX2.le. (-1.0 * wave_force_limmit)) RSNX2 =  -1.0 * wave_force_limmit
            !where(RSNY2.le. (-1.0 * wave_force_limmit)) RSNY2 =  -1.0 * wave_force_limmit

!            where(RSNX2.gt. wave_force_limmit) RSNX2 =  0.0d0
!            where(RSNY2.gt. wave_force_limmit) RSNY2 =  0.0d0
!
!            where(RSNX2.le. (-1.0 * wave_force_limmit)) RSNX2 =  0.0d0
!            where(RSNY2.le. (-1.0 * wave_force_limmit)) RSNY2 =  0.0d0
!stop
!++
!!!!#endif

    !        endif

    else
        NUOPC4WAV = .false.
        write(info,*) subname,' --- no wave forcing exchange / waves are not all connected --- '
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
        !print *, info
        !stop
    endif        
    !-----------------------------------------
    !   IMPORT from ATM
    !-----------------------------------------
   meteo_forcing= .true.
   do num = 1,fldsToAdc_num
      if (fldsToAdc(num)%shortname == 'pmsl')    meteo_forcing = meteo_forcing .and. fldsToAdc(num)%connected
      if (fldsToAdc(num)%shortname == 'imwh10m') meteo_forcing = meteo_forcing .and. fldsToAdc(num)%connected
      if (fldsToAdc(num)%shortname == 'izwh10m') meteo_forcing = meteo_forcing .and. fldsToAdc(num)%connected
   end do
    
    if ( meteo_forcing) then
        !NWS = 39   ! over write NWS option to be sure we incldue wind forcing
        
        !WTIMINC = wrf_int 
        !WTIMINC = wrf_int +  wrf_num / wrf_den
        WTIMINC =  nint(timeStepAbs)  !TODO: Get it from coupler based on the time slots.
                                      !TODO: This implemetation wotks for one time slot for wind and wave right now.
                                      !TODO: 
        !Get and fill imported fields
        ! <<<<< RECEIVE and UN-PACK pmsl
        call State_getFldPtr_(ST=importState,fldname='pmsl',fldptr=dataPtr_pmsl,rc=rc,dump=.false.,timeStr=timeStr)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        !print *, 'size > 5 dataPtr_pmsl>', size(dataPtr_pmsl) !, dataPtr_pmsl(5:)
        
        !-----------------------------------------
        ! <<<<< RECEIVE and UN-PACK imwh10m    V-Y wind comp
        call State_getFldPtr_(ST=importState,fldname='imwh10m',fldptr=dataPtr_imwh10m,rc=rc,dump=.false.,timeStr=timeStr)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        !-----------------------------------------
        ! <<<<< RECEIVE and UN-PACK izwh10m    U-X  wind comp
        call State_getFldPtr_(ST=importState,fldname='izwh10m',fldptr=dataPtr_izwh10m,rc=rc,dump=.false.,timeStr=timeStr)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        
        write(info,*) subname,' --- meteo forcing exchange OK / atm feilds are all connected --- / Model advances '
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
        !print *, info

        ! Allocate arrays for radiation stresses.
        IF(.NOT.ALLOCATED(WVNX1)) ALLOCATE(WVNX1(1:NP))
        IF(.NOT.ALLOCATED(WVNY1)) ALLOCATE(WVNY1(1:NP))
        IF(.NOT.ALLOCATED(PRN1) ) ALLOCATE(PRN1 (1:NP))

        IF(.NOT.ALLOCATED(WVNX2)) ALLOCATE(WVNX2(1:NP))
        IF(.NOT.ALLOCATED(WVNY2)) ALLOCATE(WVNY2(1:NP))
        IF(.NOT.ALLOCATED(PRN2) ) ALLOCATE(PRN2 (1:NP))

        !print *, 'maxval(WVNX2)', maxval(WVNX2)
!        PRINT*, "Enter ModelAdvance() : ", ienter ; 
!        ienter = ienter + 1 ;
!        PRINT*, " " ; 

        WVNX1 = WVNX2   
        WVNY1 = WVNY2  
        PRN1  = PRN2


!        WRITE(*,'(A,4E)') "  In ModelAdvance() 1:" , MAXVAL(WVNX1), MAXVAL(WVNX2), & 
!            MAXVAL(WVNY1), MAXVAL(WVNY2) ; 

        !call UPDATER( dataPtr_izwh10m(:), dataPtr_imwh10m(:), dataPtr_pmsl(:),3)
       
        ! Fill owned nodes from imported data to model variable
        !TODO: unit check
        do i1 = 1, mdataOut%NumOwnedNd, 1
            WVNX2(mdataOut%owned_to_present_nodes(i1)) =  dataPtr_izwh10m(i1) !* 0.0  !zonal is u-comp or x-comp
        end do

        do i1 = 1, mdataOut%NumOwnedNd, 1
            WVNY2(mdataOut%owned_to_present_nodes(i1)) =  dataPtr_imwh10m(i1) !* 0.0  !Meridionalis v-comp or y-comp
        end do
        
        do i1 = 1, mdataOut%NumOwnedNd, 1
!            PRN2(mdataOut%owned_to_present_nodes(i1) ) = dataPtr_pmsl(i1) / (1025 * 9.81)    !convert Pascal to mH2O
!++ GML  in ADCIRC global.F RhoWat0=1000.D0; g = 9.80665
            PRN2(mdataOut%owned_to_present_nodes(i1) ) = dataPtr_pmsl(i1) / (1000 * 9.80665)    !convert Pascal to mH2O
!++        
          !if ( abs(dataPtr_pmsl(i1) ).gt. 1e11)  then
          !  STOP '  dataPtr_pmsl > mask '     
          !end if
        end do
         
        call UPDATER( WVNX2(:), WVNY2(:), PRN2(:),3)

 !       WRITE(*,'(A,4E)') "  In ModelAdvance() 2:" , MAXVAL(WVNX1), MAXVAL(WVNX2), & 
 !           MAXVAL(WVNY1), MAXVAL(WVNY2) ; 

 !       if (first_exchange .and. sum(PRN1) .le. 1.0) then
 !      ! DW:  band-aid fix by zeroing out wind velocity and atm pressure.
 !      !      Todo: look at atmesh code to see if there could be a better fix.  
        if ( first_exchange .or. sum(PRN1)  .le. 1.0) then
          WVNX2 = 1e-10
          WVNY2 = 1e-10
          WVNX1 = 1e-10
          WVNY1 = 1e-10

          PRN2 = 10.0  !hard coded to handel the 1st exchange zeros problem :TODO! Need to resolve this!
          PRN1 = PRN2
!++ GML 20210308
          RSNX2 = 0.0d0
          RSNY2 = 0.0d0
!++
          first_exchange = .false.
        end if  


    
 !       WRITE(*,'(A,4E)') "  In ModelAdvance() 3:" , MAXVAL(WVNX1), MAXVAL(WVNX2), & 
 !           MAXVAL(WVNY1), MAXVAL(WVNY2) ; 


        !if (sum(PRN1) .eq. 0.0 ) then
        !  PRN1 = 10000.0
        !end if  
    
        !if (sum(WVNX2) .eq. 0.0 ) then
        !  WVNX2 = 8.0
        !end if  
    
        !if (sum(WVNX1) .eq. 0.0 ) then
        !  WVNX1 = 8.0
        !end if  
        

        !if (sum(WVNY2) .eq. 0.0 ) then
        !  WVNY2 = -8.0
        !end if  
    
        !if (sum(WVNY1) .eq. 0.0 ) then
        !  WVNY1 = -8.0
        !end if  


        !where(abs(PRN1).gt. 1e11)  PRN1 =  1e4
        !where(abs(PRN2).gt. 1e11)  PRN2 =  1e4

        
        ! where(abs(WVNX1).gt. 1e6)  WVNX1 =  8.0
        !where(abs(WVNX2).gt. 1e6)  WVNX2 =  8.0

        !where(abs(WVNY1).gt. 1e6)  WVNY1 =  -8.0
        !where(abs(WVNY2).gt. 1e6)  WVNY2 =  -8.0

            
        !PRN2 = 10000.0
        !PRN1 = 10000.0
        !WVNX2 =  8.0
        !WVNX1 =  8.0
        !WVNY2 = -8.0
        !WVNY1 = -8.0       
           
        !where(dataPtr_pmsl .gt. 1e20)  dataPtr_pmsl =  10e4
        !where(dataPtr_pmsl .lt. 8e4 )  dataPtr_pmsl =  10e4        
        
        !print *, 'size(dataPtr_pmsl) > ',  size(dataPtr_pmsl)
        !print *, 'size(PRN2)         > ',  size(PRN2)


        !
        !where((WVNX2).gt. 20)  WVNX2 =  20
        !where((WVNY2).gt. 20)  WVNY2 =  20
        !where((PRN2) .gt. 1e20)   PRN2 =  1e4                

        !where((WVNX1).gt. 20)  WVNX1 =  20
        !where((WVNY1).gt. 20)  WVNY1 =  20
        !where((PRN1).gt.  1e20)   PRN1 =  1e4  
        
        !PRN1 =  1e4  
        !PRN2 =  1e4  
        
        ! Ghost nodes update for meteo infos
        !call UPDATER( WVNX2(:), WVNY2(:), PRN2(:),3)

        !print *, 'in cap before maxval(WVNX2)', maxval(WVNX2)
        !print *, 'in cap before maxval(WVNY2)', maxval(WVNY2)
         !WVNX2 = 10.0
         !WVNY2 = 10.0
         !PRN2  = 10000.0

         !WVNX2 = 0.0
         !WVNY2 = 0.0
         !PRN2  = 10000.0

        !print *, 'in cap after maxval(WVNX2)', maxval(WVNX2)
        !print *, 'in cap after maxval(WVNY2)', maxval(WVNY2)    


    else
        NUOPC4MET = .false.
        write(info,*) subname,' --- no meteo forcing exchange / atm feilds are not all connected --- '
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
        !print *, info
        !stop
    endif
! GML 
!   surge_forcing= .false.
!++
   surge_forcing= .true.
   do num = 1,fldsFrAdc_num
      if (fldsFrAdc(num)%shortname == 'zeta') surge_forcing = surge_forcing .and. fldsFrAdc(num)%connected
      if (fldsFrAdc(num)%shortname == 'velx') surge_forcing = surge_forcing .and. fldsFrAdc(num)%connected
      if (fldsFrAdc(num)%shortname == 'vely') surge_forcing = surge_forcing .and. fldsFrAdc(num)%connected
   end do

    !------------------------------------------
    !---------------  RUN  --------------------
    !------------------------------------------
    write(info,*) subname,' --- run phase 2 called --- '
    !print *,      subname,' --- run phase 2 called --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

    CALL ADCIRC_Run(nCplADC)

    write(info,*) subname,' --- run phase 3 called ---  nCplADC = ',nCplADC
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
    !-------------------------------------------


    if (surge_forcing) then
      !-----------------------------------------
      !   EXPORT
      !-----------------------------------------
      !pack and send exported fields
      allocate (tmp(mdataOut%NumOwnedNd))

      ! >>>>> PACK and send ZETA
      call State_getFldPtr_(ST=exportState,fldname='zeta',fldptr=dataPtr_zeta, &
        rc=rc,dump=.true.,timeStr=timeStr)
      !call State_getFldPtr(ST=exportState,fldname='zeta',fldptr=dataPtr_zeta, &
      !  rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      !fill only owned nodes for tmp vector
      do i1 = 1, mdataOut%NumOwnedNd, 1
          tmp(i1) = ETA2(mdataOut%owned_to_present_nodes(i1))
!++ GML
!          tmp(i1) = 0.0d0
      end do
!++ GML
!      where(abs(tmp).gt. 10.d0)  tmp =  0.0d0
!++
      !assign to field
      dataPtr_zeta = tmp
      !----------------------------------------
      ! >>>>> PACK and send VELX
      call State_getFldPtr(ST=exportState,fldname='velx',fldptr=dataPtr_velx,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      !fill only owned nodes for tmp vector
      do i1 = 1, mdataOut%NumOwnedNd, 1
          tmp(i1) = UU2(mdataOut%owned_to_present_nodes(i1))
!++ GML
!          if(dp(i1)<0.5)then
!          tmp(i1) = 0.0d0
!          endif
      end do
!++ GML
!      where(abs(tmp).gt. 10.d0)  tmp =  0.0d0
!++
      !assign to field
      dataPtr_velx = tmp
      !----------------------------------------
      ! >>>>> PACK and send VELY
      call State_getFldPtr(ST=exportState,fldname='vely',fldptr=dataPtr_vely,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      !fill only owned nodes for tmp vector
      do i1 = 1, mdataOut%NumOwnedNd, 1
          tmp(i1) = VV2(mdataOut%owned_to_present_nodes(i1))
!++ GML
!          if(dp(i1)<0.5)then
!          tmp(i1) = 0.0d0
!          endif
      end do
!++ GML
!      where(abs(tmp).gt. 10.d0)  tmp =  0.0d0
!++
      !assign to field
      dataPtr_vely = tmp
    else
      write(info,*) subname,' --- no surge forcing for wave. 1way coupled WW3 -> ADC  ---'
      call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
      !print *, info
    endif

    !----------------------------------------
!
!
  end subroutine

!-----------------------------------------------------------
  !> Retrieve a pointer to a field's data array from inside an ESMF_State object.
  !!
  !! @param ST the ESMF_State object
  !! @param fldname name of the fields
  !! @param fldptr pointer to 1D array
  !! @param rc return code
  subroutine State_GetFldPtr_(ST, fldname, fldptr, rc, dump,timeStr)
    type(ESMF_State), intent(in) :: ST
    character(len=*), intent(in) :: fldname
    real(ESMF_KIND_R8), pointer, intent(in) :: fldptr(:)
    integer, intent(out), optional :: rc
    logical, intent(in), optional  :: dump
    character(len=128),intent(inout), optional :: timeStr
    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
    character(len=*),parameter :: subname='(adc_cap:State_GetFldPtr)'

    call ESMF_StateGet(ST, itemName=trim(fldname), field=lfield, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (present(rc)) rc = lrc

    !!Thu 01 Jun 2017 01:04:19 PM UTC  did not show any promise
    !! added to test the halo update
    !!TODO: this should not be here. It should initilize once
    !call ESMF_FieldHaloStore ( lfield, routehandle = ATM_HaloRouteHandel, rc=lrc)
    !if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    !! 
    !!Halo update
    !call ESMF_FieldHalo      ( lfield, ATM_HaloRouteHandel, rc=lrc)
    !if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    !!
    !!TODO: this should not be here. It should finalize once
    !call ESMF_FieldHaloRelease (routehandle = ATM_HaloRouteHandel, rc=lrc)
    !if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return 
    !if (present(rc)) rc = lrc
    !write(info,*) ' --- ATM  halo routehandel in work >>>>>  ---'
    !call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

    if (dump) then
       if (.not. present(timeStr)) timeStr="_"
        call ESMF_FieldWrite(lfield, &
        fileName='field_ocn_'//trim(fldname)//trim(timeStr)//'.nc', &
        rc=rc,overwrite=.true.)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
    endif
  end subroutine State_GetFldPtr_


  !> Retrieve a pointer to a field's data array from inside an ESMF_State object.
  !!
  !! @param ST the ESMF_State object
  !! @param fldname name of the fields
  !! @param fldptr pointer to 1D array
  !! @param rc return code
  subroutine State_GetFldPtr(ST, fldname, fldptr, rc)
    type(ESMF_State), intent(in) :: ST
    character(len=*), intent(in) :: fldname
    real(ESMF_KIND_R8), pointer, intent(in) :: fldptr(:)
    integer, intent(out), optional :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
    character(len=*),parameter :: subname='(WAV:State_GetFldPtr)'

    call ESMF_StateGet(ST, itemName=trim(fldname), field=lfield, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (present(rc)) rc = lrc
  end subroutine State_GetFldPtr

  subroutine CheckImport_not_comp(model, rc)
    type(ESMF_GridComp)   :: model
    integer, intent(out)  :: rc
    
    ! This is the routine that enforces correct time stamps on import Fields
    
    ! local variables
    type(ESMF_Clock)        :: driverClock
    type(ESMF_Time)         :: startTime, currTime
    type(ESMF_State)        :: importState
    type(ESMF_Field)        :: field
    logical                 :: atCorrectTime

    rc = ESMF_SUCCESS
    return
    
!    ! query Component for the driverClock
!    call NUOPC_ModelGet(model, driverClock=driverClock, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
    
!    ! get the start time and current time out of the clock
!    call ESMF_ClockGet(driverClock, startTime=startTime, &
!      currTime=currTime, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
    
   
  end subroutine

  !-----------------------------------------------------------------------------
  


  !-----------------------------------------------------------------------------
  !> Called by NUOPC at the end of the run to clean up.  The cap does
  !! this simply by calling ADCIRC_Final.
  !!
  !! @param model the ESMF_GridComp object
  !! @param rc return code
    subroutine ADCIRC_model_finalize(model, rc)

        ! input arguments
        type(ESMF_GridComp)  :: model
        integer, intent(out) :: rc

        ! local variables
        type(ESMF_Clock)     :: clock
        type(ESMF_Time)                        :: currTime
        character(len=*),parameter  :: subname='(adc_cap:adc_model_finalize)'

        rc = ESMF_SUCCESS

        write(info,*) subname,' --- finalize called --- '
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

        call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

        call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

        CALL ADCIRC_Final()

        write(info,*) subname,' --- finalize completed --- '
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

    end subroutine ADCIRC_model_finalize





end module
