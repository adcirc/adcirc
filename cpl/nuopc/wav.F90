!>
!! @mainpage ADCIRC NUOPC Cap
!! @author Saeed Moghimi (moghimis@gmail.com)
!! @date 15/1/17 Original documentation
!------------------------------------------------------
!LOG-----------------
! Will keep the import part in case we want to change this for a new wave model.
!
!
module WAV
  !-----------------------------------------------------------------------------
  ! WAV Component.
  !-----------------------------------------------------------------------------
  use mpi
  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS        => SetServices,    &
    model_label_SetClock    => label_SetClock, &
    model_label_Advance     => label_Advance,  &
    model_label_CheckImport => label_CheckImport, &    
    model_label_Finalize    => label_Finalize

  use wav_mod, only: meshdata
  use wav_mod, only: create_parallel_esmf_mesh_from_meshdata
  !use wav_mod, only: extract_parallel_data_from_mesh_orig
  !use wav_mod, only: ww3_cpl_int,ww3_cpl_num,ww3_cpl_den
  use wav_mod, only: SXX, SYY, SXY
  use wav_mod, only: read_config

  !read from netcdf file
  use wav_mod, only: init_ww3_nc, read_ww3_nc , ww3_from_file
  use wav_mod, only: construct_meshdata_from_netcdf
  
  implicit none
  private
  
  public SetServices
  
   type fld_list_type
        character(len=64) :: stdname
        character(len=64) :: shortname
        character(len=64) :: unit
        logical           :: assoc    ! is the farrayPtr associated with internal data
        real(ESMF_KIND_R8), dimension(:), pointer :: farrayPtr
    end type fld_list_type

    integer,parameter :: fldsMax = 100
    integer :: fldsToWav_num = 0
    type (fld_list_type) :: fldsToWav(fldsMax)
    integer :: fldsFrWav_num = 0
    type (fld_list_type) :: fldsFrWav(fldsMax)

  type(meshdata),save  :: mdataInw, mdataOutw

  character(len=2048):: info
  integer :: dbrc     ! temporary debug rc value


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

    ! Local Variables
    type(ESMF_VM)                :: vm
    character(len=*),parameter   :: subname='(WAV:SetServices)'

    rc = ESMF_SUCCESS
    
    ! readd config file
    call read_config()
    
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
    
    ! attach specializing method(s)
!    call NUOPC_CompSpecialize(model, specLabel=model_label_SetClock, &
!      specRoutine=SetClock, rc=rc)
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

    call NUOPC_CompSpecialize(model, specLabel=model_label_CheckImport, &
      specPhaseLabel="RunPhase1", specRoutine=CheckImport, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_CompSpecialize(model, specLabel=model_label_Finalize, &
      specRoutine=WAV_model_finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out


    if (ww3_from_file) then
      call init_ww3_nc()
      write(info,*) subname,' --- Read wave info from file --- '
      !print *,      info
      call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)
    endif

    write(info,*) subname,' --- adc SetServices completed --- '
    !print *,      subname,' --- adc SetServices completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)
  end subroutine
  
  !-----------------------------------------------------------------------------
    !> First initialize subroutine called by NUOPC.  The purpose
    !! is to set which version of the Initialize Phase Definition (IPD)
    !! to use.
    !!
    !! For this WAV cap, we are using IPDv01.
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
    integer              :: num,i
    character(len=*),parameter  :: subname='(WAV:AdvertiseFields)'

    rc = ESMF_SUCCESS


    call WAV_FieldsSetup()
!

      do num = 1,fldsToWav_num
          !print *,  "fldsToWav_num  ", fldsToWav_num
          !print *,  "fldsToWav(num)%shortname  ", fldsToWav(num)%shortname
          !print *,  "fldsToWav(num)%stdname  ", fldsToWav(num)%stdname

          write(info,*) subname,  "fldsToWav(num)%shortname  ", fldsToWav(num)%shortname
          call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
     end do

      call WAV_AdvertiseFields(importState, fldsToWav_num, fldsToWav, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

!----------------------------------------------------------------
    do num = 1,fldsFrWav_num
        !print *,  "fldsFrWav_num  ", fldsFrWav_num
        !print *,  "fldsFrWav(num)%shortname  ", fldsFrWav(num)%shortname
        !print *,  "fldsFrWav(num)%stdname  ", fldsFrWav(num)%stdname
        write(info,*) subname,"fldsFrWav(num)%stdname  ", fldsFrWav(num)%stdname
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

    end do
!
    call WAV_AdvertiseFields(exportState, fldsFrWav_num, fldsFrWav, rc)
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
    subroutine WAV_AdvertiseFields(state, nfields, field_defs, rc)
        type(ESMF_State), intent(inout)             :: state
        integer,intent(in)                          :: nfields
        type(fld_list_type), intent(inout)          :: field_defs(:)
        integer, intent(inout)                      :: rc

        integer                                     :: i
        character(len=*),parameter  :: subname='(WAV:WAV_AdvertiseFields)'

        rc = ESMF_SUCCESS

        do i = 1, nfields
          !print *, 'Advertise: '//trim(field_defs(i)%stdname)//'---'//trim(field_defs(i)%shortname)
          write(info,*) subname,'Advertise: '//trim(field_defs(i)%stdname)//'---'//trim(field_defs(i)%shortname)
          call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)


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
        write(info,*) subname,'---- Done WAV_AdvertiseFields '
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
    end subroutine WAV_AdvertiseFields




   subroutine WAV_FieldsSetup
        integer                     :: rc
        character(len=*),parameter  :: subname='(WAV:WAV_FieldsSetup)'



        !check if read from file (1way coupled) or from wave model (2way coupled)
        if (.not. ww3_from_file) then
          !--------- import fields to wav -------------
          !                     num,        fldlist,       stdname,                       transferOffer, data, shortname
          call fld_list_add(num=fldsToWav_num, fldlist=fldsToWav, stdname="sea_surface_height_above_sea_level",  shortname= "zeta" )
          call fld_list_add(num=fldsToWav_num, fldlist=fldsToWav, stdname="surface_eastward_sea_water_velocity", shortname= "velx" )
          call fld_list_add(num=fldsToWav_num, fldlist=fldsToWav, stdname="surface_northward_sea_water_velocity",shortname= "vely" )
       else
          write(info,*) subname,  "No fields for WW3 importState advertised. Only WW3 -> ADC from file."
          call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
       endif


!       !--------- export fields from Sea Adc -------------
        ! First add Non-standard fields to dictionary
        ! exportable field: wav_x a Non-standard field

        call NUOPC_FieldDictionaryAddEntry("eastward_radiation_stress",  "mx", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

        call NUOPC_FieldDictionaryAddEntry("northward_radiation_stress", "mx", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

        call NUOPC_FieldDictionaryAddEntry("cross_radiation_stress",    "mx", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

         call fld_list_add(num=fldsFrWav_num, fldlist=fldsFrWav, stdname="eastward_radiation_stress", shortname= "sxx")
         call fld_list_add(num=fldsFrWav_num, fldlist=fldsFrWav, stdname="northward_radiation_stress",shortname= "syy")
         call fld_list_add(num=fldsFrWav_num, fldlist=fldsFrWav, stdname="cross_radiation_stress",    shortname= "sxy")

        write(info,*) subname,' --- Passed--- '
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
        !print *,      subname,' --- Passed --- '
    end subroutine WAV_FieldsSetup


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
        character(len=*), parameter :: subname='(WAV:fld_list_add)'

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
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
        !print *,      subname,' --- Passed --- '
    end subroutine fld_list_add



  
  !-----------------------------------------------------------------------------
    !> Called by NUOPC to realize import and export fields.

    !! The fields to import and export are stored in the fldsToWav and fldsFrWav
    !! arrays, respectively.  Each field entry includes the standard name,
    !! information about whether the field's grid will be provided by the cap,
    !! and optionally a pointer to the field's data array.  Currently, all fields
    !! are defined on the same mesh defined by the cap.
    !! The fields are created by calling WAV::adcirc_XXXXXXXXXXXXXXXXXXX.
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
    type(ESMF_TimeInterval) :: WAVTimeStep
    type(ESMF_Field)        :: field
    !Saeed added
    type(meshdata)               :: mdataw
    type(ESMF_Mesh)              :: ModelMesh,meshIn,meshOut
    type(ESMF_VM)                :: vm
    integer                      :: localPet, petCount
    character(len=*),parameter   :: subname='(WAV:RealizeFieldsProvidingGrid)'

    rc = ESMF_SUCCESS

    !print *,"WAV ..1.............................................. >> "
    !> \details Get current ESMF VM.
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !print *,"WAV ..2.............................................. >> "
    ! Get query local pet information for handeling global node information
    call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    ! call ESMF_VMPrint(vm, rc=rc)

    !print *,localPet,"< LOCAL pet, WAV ..3.............................................. >> "
    !! Assign VM to mesh data type.
    mdataw%vm = vm

    
    if (.not. ww3_from_file) then
      ! create a Mesh object for Fields
      print *, ' not the porpuse of this cap. STOP! >>'
      write(info,*) subname,' not the porpuse of this cap. STOP! >>'
      call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
      stop
      !call extract_parallel_data_from_mesh_orig('.', mdataw, localPet)
    else
      call construct_meshdata_from_netcdf(mdataw)
    endif
    
    call create_parallel_esmf_mesh_from_meshdata(mdataw,ModelMesh )
    !

    call ESMF_MeshWrite(ModelMesh, filename="wav_mesh.nc", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
      



    meshIn  = ModelMesh ! for now out same as in
    meshOut = meshIn

    mdataInw  = mdataw
    mdataOutw = mdataw

    !print *,"..................................................... >> "
    !print *,"NumNd", mdataw%NumNd
    !print *,"NumOwnedNd", mdataw%NumOwnedNd
    !print *,"NumEl", mdataw%NumEl
    !print *,"NumND_per_El", mdataw%NumND_per_El
    !print *,"NumOwnedNd mdataOutw", mdataOutw%NumOwnedNd

        call WAV_RealizeFields(importState, meshIn , mdataw, fldsToWav_num, fldsToWav, "WAV import", rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
!
        call WAV_RealizeFields(exportState, meshOut, mdataw, fldsFrWav_num, fldsFrWav, "WAV export", rc)
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
  !! an existing WAV pointer is referenced.
  !!
  !! @param state the ESMF_State object to add fields to
  !! @param grid the ESMF_Grid object on which to define the fields
  !! @param nfields number of fields
  !! @param field_defs array of fld_list_type indicating the fields to add
  !! @param tag used to output to the log
  !! @param rc return code
  subroutine WAV_RealizeFields(state, mesh, mdata, nfields, field_defs, tag, rc)

    type(ESMF_State), intent(inout)             :: state
    type(ESMF_Mesh), intent(in)                 :: mesh
    type(meshdata)                              :: mdata
    integer, intent(in)                         :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    character(len=*), intent(in)                :: tag
    integer, intent(inout)                      :: rc


    type(ESMF_Field)                            :: field
    integer                                     :: i
    character(len=*),parameter  :: subname='(WAV:WAV_RealizeFields)'

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
        endif
    enddo

        write(info,*) subname,' --- passed--- '
        !print *,      subname,' --- OUT --- '
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=rc)
  end subroutine WAV_RealizeFields
  !-----------------------------------------------------------------------------

!  subroutine SetClock_mine(model, rc)
!    type(ESMF_GridComp)  :: model
!    integer, intent(out) :: rc

!    ! local variables
!    type(ESMF_Clock)              :: clock
!    type(ESMF_TimeInterval)       :: WAVTimeStep

!    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
!    call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out


    ! initialize internal clock
    ! - on entry, the component clock is a copy of the parent clock
    ! - the parent clock is on the slow timescale hwrf timesteps
    ! - reset the component clock to have a timeStep that is for adc-wav of the parent
    !   -> timesteps
    
    !call ESMF_TimeIntervalSet(WAVTimeStep, s=ww3_cpl_int, sN=ww3_cpl_num, sD=ww3_cpl_den, rc=rc) ! 5 minute steps
    !TODO: use nint !!?

!    call ESMF_TimeIntervalSet(WAVTimeStep, s=ww3_cpl_int, rc=rc) ! 5 minute steps
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!       line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
!    call NUOPC_CompSetClock(model, clock, WAVTimeStep, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out


!    print *, "WAV Timeinterval1 = "
!    call ESMF_TimeIntervalPrint(WAVTimeStep, options="string", rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!        line=__LINE__, &
!        file=__FILE__)) &
!        return  ! bail out    


    ! initialize internal clock
    ! - on entry, the component clock is a copy of the parent clock
    ! - the parent clock is on the slow timescale hwrf timesteps
    ! - reset the component clock to have a timeStep that is for adc-wav of the parent
    !   -> timesteps
    
    !call ESMF_TimeIntervalSet(WAVTimeStep, s=     adc_cpl_int, sN=adc_cpl_num, sD=adc_cpl_den, rc=rc) ! 5 minute steps
    !TODO: use nint !!?


!  end subroutine

  !-----------------------------------------------------------------------------
  ! From CICE model uses same clock as parent gridComp
!  subroutine SetClock(model, rc)
!    type(ESMF_GridComp)  :: model
!    integer, intent(out) :: rc
    
    ! local variables
!    type(ESMF_Clock)              :: clock
!    type(ESMF_TimeInterval)       :: WAVTimeStep, timestep
!    character(len=*),parameter  :: subname='(wav_cap:SetClock)'

!    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
!    call ESMF_GridCompGet(model, clock=clock, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
    !call ESMF_TimeIntervalSet(WAVTimeStep, s=ww3_cpl_int, sN=ww3_cpl_num, sD=ww3_cpl_den, rc=rc) ! 5 minute steps
    ! tcraig: dt is the cice thermodynamic timestep in seconds
!    call ESMF_TimeIntervalSet(timestep, s=ww3_cpl_int, rc=rc)
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
!    call ESMF_TimeIntervalSet(WAVTimeStep, s=ww3_cpl_int, rc=rc) 
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
!    call NUOPC_CompSetClock(model, clock, WAVTimeStep, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
    
!  end subroutine
    !-----------------------------------------------------------------------------

      !> Called by NUOPC to advance the WAV model a single timestep >>>>>
      !! <<<<<<<<<  TODO: check! this is not what we want!!!.
      !!
      !! This subroutine copies field data out of the cap import state and into the
      !! model internal arrays.  Then it calls WAV_Run to make a NN timesteps.
      !! Finally, it copies the updated arrays into the cap export state.
      !!
      !! @param model an ESMF_GridComp object
      !! @param rc return code
      !-----------------------------------------------------------------------------
  subroutine ModelAdvance(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_State)              :: importState, exportState
    type(ESMF_Time)               :: currTime
    type(ESMF_TimeInterval)       :: timeStep
    character(len=*),parameter    :: subname='(WAV:ModelAdvance)'
    !tmp vector
    real(ESMF_KIND_R8), pointer   :: tmp(:)

    !imports
    real(ESMF_KIND_R8), pointer   :: dataPtr_zeta(:)
    real(ESMF_KIND_R8), pointer   :: dataPtr_velx(:)
    real(ESMF_KIND_R8), pointer   :: dataPtr_vely(:)

    ! exports
    real(ESMF_KIND_R8), pointer   :: dataPtr_sxx(:)
    real(ESMF_KIND_R8), pointer   :: dataPtr_syy(:)
    real(ESMF_KIND_R8), pointer   :: dataPtr_sxy(:)

    type(ESMF_StateItem_Flag)     :: itemType
    type(ESMF_Mesh)               :: mesh
    type(ESMF_Field)              :: lfield
    character(len=128)            :: fldname,timeStr
    integer                       :: i1
    ! local variables for Get methods
    integer :: YY, MM, DD, H, M, S

    rc = ESMF_SUCCESS
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
      preString="------>Advancing WAV from: ", rc=rc)
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
      preString="------------------WAV-------------> to: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimeGet(currTime, yy=YY, mm=MM, dd=DD, h=H, m=M, s=S, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    !print *      , "currTime = ", YY, "/", MM, "/", DD," ", H, ":", M, ":", S

    call ESMF_TimeGet(currTime, timeStringISOFrac=timeStr , rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    if (.not. ww3_from_file) then
      !-----------------------------------------
      !   IMPORT
      !-----------------------------------------
      !Get and fill imported fields
      ! <<<<< RECEIVE and UN-PACK ZETA
      !call State_getFldPtr(ST=importState,fldname='zeta',fldptr=dataPtr_zeta,rc=rc)
      call State_getFldPtr_(ST=importState,fldname='zeta',fldptr=dataPtr_zeta, &
        rc=rc,dump=.true.,timeStr=timeStr)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      !fill owned nodes from imported data to model variable
      !do i1 = 1, mdataOutw%NumOwnedNd, 1
      !    WLEV(mdataOut%owned_to_present_nodes(i1))=dataPtr_zeta(i1)
      !end do
      !-----------------------------------------
      ! <<<<< RECEIVE and UN-PACK VELX
      call State_getFldPtr(ST=importState,fldname='velx',fldptr=dataPtr_velx,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      !fill owned nodes from imported data to model variable
      !do i1 = 1, mdataOutw%NumOwnedNd, 1
      !    UU(mdataOutw%owned_to_present_nodes(i1))=dataPtr_velx(i1)
      !end do
      !-----------------------------------------
      ! <<<<< RECEIVE and UN-PACK VELY
      call State_getFldPtr(ST=importState,fldname='vely',fldptr=dataPtr_vely,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      !fill owned nodes from imported data to model variable
      !do i1 = 1, mdataOutw%NumOwnedNd, 1
      !    VV(mdataOutw%owned_to_present_nodes(i1))=dataPtr_vely(i1)
      !end do
     endif

    !-----------------------------------------
    !   EXPORT
    !-----------------------------------------
    if (ww3_from_file) then
      !update sxx, syy, sxy from nearset time in ww3 netcdf file

      !TODO: update file name!!!!

      call read_ww3_nc(currTime)

    end if

    !pack and send exported fields
    allocate (tmp(mdataOutw%NumOwnedNd))

    ! >>>>> PACK and send SXX
    !call State_getFldPtr(ST=exportState,fldname='sxx',fldptr=dataPtr_sxx,rc=rc)
    call State_getFldPtr_(ST=exportState,fldname='sxx',fldptr=dataPtr_sxx, &
      rc=rc,dump=.true.,timeStr=timeStr)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

      !print *, 'mdataOutw%NumOwnedNd > ',mdataOutw%NumOwnedNd, 'SXX > ', SXX(1:10,1)


    !fill only owned nodes for tmp vector
    do i1 = 1, mdataOutw%NumOwnedNd, 1
        tmp(i1) = SXX(mdataOutw%owned_to_present_nodes(i1),1)
        !tmp(i1) = 0.0 * real(i1) / 5000.0
    end do
    !assign to field
    dataPtr_sxx = tmp
    !----------------------------------------
    ! >>>>> PACK and send SYY
    call State_getFldPtr(ST=exportState,fldname='syy',fldptr=dataPtr_syy,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !fill only owned nodes for tmp vector
    do i1 = 1, mdataOutw%NumOwnedNd, 1
        tmp(i1) = SYY(mdataOutw%owned_to_present_nodes(i1),1)
    end do
    !assign to field
    dataPtr_syy = tmp
    !----------------------------------------
    ! >>>>> PACK and send SXY
    call State_getFldPtr(ST=exportState,fldname='sxy',fldptr=dataPtr_sxy,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !fill only owned nodes for tmp vector
    do i1 = 1, mdataOutw%NumOwnedNd, 1
        tmp(i1) = SXY(mdataOutw%owned_to_present_nodes(i1),1)
    end do
    !assign to field
    dataPtr_sxy = tmp
    !----------------------------------------


    !! TODO:  not a right thing to do. we need to fix the grid object mask <<<<<<
    !where(dataPtr_sxx.gt.3e4) dataPtr_sxx = 0.0
    !where(dataPtr_syy.gt.3e4) dataPtr_syy = 0.0
    !where(dataPtr_sxy.gt.3e4) dataPtr_sxy = 0.0
!
  end subroutine
!
  !-----------------------------------------------------------------------
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
    character(len=*),parameter :: subname='(wav_cap:State_GetFldPtr)'

    call ESMF_StateGet(ST, itemName=trim(fldname), field=lfield, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (present(rc)) rc = lrc

    if (dump) then
       if (.not. present(timeStr)) timeStr="_"
        call ESMF_FieldWrite(lfield, &
        fileName='field_wave_'//trim(fldname)//trim(timeStr)//'.nc', &
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



  subroutine CheckImport(model, rc)
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
!    
!    ! get the start time and current time out of the clock
!    call ESMF_ClockGet(driverClock, startTime=startTime, &
!      currTime=currTime, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
    
   
  end subroutine


  !-----------------------------------------------------------------------
  !> Called by NUOPC at the end of the run to clean up.  The cap does
  !! this simply by calling WAV_Final.
  !!
  !! @param gcomp the ESMF_GridComp object
  !! @param rc return code
    subroutine WAV_model_finalize(gcomp, rc)

        ! input arguments
        type(ESMF_GridComp)  :: gcomp
        integer, intent(out) :: rc

        ! local variables
        type(ESMF_Clock)     :: clock
        type(ESMF_Time)                        :: currTime
        character(len=*),parameter  :: subname='(WAV:wav_model_finalize)'

        rc = ESMF_SUCCESS

        write(info,*) subname,' --- finalize called --- '
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

        call NUOPC_ModelGet(gcomp, modelClock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

        call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

        write(info,*) subname,' --- finalize completed --- '
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

    end subroutine WAV_model_finalize

end module
