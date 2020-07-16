!>
!! @mainpage ADCIRC NUOPC Cap
!! @author Saeed Moghimi (moghimis@gmail.com)
!! @date 15/1/17 Original documentation
!------------------------------------------------------
!LOG-----------------
!
!
!

!############# NOTES #############
! 1- HWRF staggerd grid on output? Is there interpolation
!    on the cell centers or it is out puted on staggered system?
! 2- Here we assume all the variables on the grid center
! 3- No time interpolation on the fields is applided here.
! 4- This may force us to use mediator. to read two wind field and
! update models on every coupling time step?
! 5- Or pass two wind filed and let the models take care of the time interpolation
!    using their own internal routines (makes more sense to me)
! 6- might make sense to prepare wind for ADCIRC and let it to read
!    and pass it to wave at every coupling time step (I do not like this because
!    it may limit our ability to couple HWRF to our system down the road
!    however may make the life easier).
! 7- >>>>>>>>>
!    A> We read HWRF and initialize both models approperiatly (check if they
!       one or two time step meteo forcing.
!    B> Read HWRF and pass to models every meteo time step and they will do
!       the time interpolation


module hwrf_cap

    use ESMF
    use NUOPC
  use NUOPC_Model, &
    model_routine_SS      => SetServices,    &
    model_label_SetClock  => label_SetClock, &
    model_label_Advance   => label_Advance,  &
    model_label_Finalize  => label_Finalize

    use hwrf_mod, only     : read_hwrf_nc
    use hwrf_mod, only     : LONS, LATS, TIMES
    use hwrf_mod, only     : UGRD10, VGRD10, PRMSL
    use hwrf_mod, only     : nlat, nlon, ntime
    use hwrf_mod, only     : FILE_NAME
    use hwrf_mod, only     : wrf_int,wrf_num,wrf_den

    ! hwrf data directory
    use hwrf_mod, only: wrf_dir,wrf_nam,update_hwrf_filename
    use hwrf_mod, only: read_config
    
    
    implicit none

    public SetServices

  type fld_list_type
    character(len=64) :: stdname
    character(len=64) :: shortname
    character(len=64) :: unit
    logical           :: assoc    ! is the farrayPtr associated with internal data
    real(ESMF_KIND_R8), dimension(:), pointer :: farrayPtr
  end type fld_list_type

  integer,parameter :: fldsMax = 100
  integer :: fldsToWrf_num = 0
  type (fld_list_type) :: fldsToWrf(fldsMax)
  integer :: fldsFrWrf_num = 0
  type (fld_list_type) :: fldsFrWrf(fldsMax)

  integer :: dbrc     ! temporary debug rc value
  character(len=2048):: info


contains

    subroutine SetServices(model, rc)
        type(ESMF_GridComp)  :: model
        integer, intent(out) :: rc
        character(len=*),parameter   :: subname='(hwrf_cap:SetServices)'
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
    
    !following MOM5 cap: we assume no need to change the clock settings
    ! attach specializing method(s)
    !call NUOPC_CompSpecialize(model, specLabel=model_label_SetClock, &
    !  specRoutine=SetClock, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_CompSpecialize(model, specLabel=model_label_Finalize, &
      specRoutine=HWRF_model_finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(info,*) subname,' --- hwrf SetServices completed --- '
    print *,      subname,' --- hwrf SetServices completed --- '
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
    integer                                     :: num
    character(len=*),parameter   :: subname='(hwrf_cap:InitializeP0)'        

    rc = ESMF_SUCCESS

    FILE_NAME =  TRIM(wrf_dir)//'/hwrf_grid.nc'

    print *, FILE_NAME
    call read_hwrf_nc()

    call HWRF_FieldsSetup()

    !no fields to import for now
    !
    !call HWRF_AdvertiseFields(importState, fldsToWrf_num, fldsToWrf, rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out

    do num = 1,fldsFrWrf_num
      print *,  "fldsFrWrf_num  ", fldsFrWrf_num
      print *,  "fldsFrWrf(num)%shortname  ", fldsFrWrf(num)%shortname
      print *,  "fldsFrWrf(num)%stdname  ", fldsFrWrf(num)%stdname

      write(info,*) subname,  "fldsFrWrf(num)%shortname  ", fldsFrWrf(num)%shortname
      call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
    end do


    call HWRF_AdvertiseFields(exportState, fldsFrWrf_num, fldsFrWrf, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(info,*) subname,' --- hwrf advertize completed --- '
    print *,      subname,' --- hwrf advertize completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)
  end subroutine



  !> Advertises a set of fields in an ESMF_State object by calling
  !! NUOPC_Advertise in a loop.
  !!
  !! @param state the ESMF_State object in which to advertise the fields
  !! @param nfields number of fields
  !! @param field_defs an array of fld_list_type listing the fields to advertise
  !! @param rc return code
  subroutine HWRF_AdvertiseFields(state, nfields, field_defs, rc)
    type(ESMF_State), intent(inout)             :: state
    integer,intent(in)                          :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    integer, intent(inout)                      :: rc

    integer                                     :: i
    character(len=*),parameter  :: subname='(hwrf_cap:HWRF_AdvertiseFields)'

    rc = ESMF_SUCCESS

    do i = 1, nfields
      print *, 'Advertise: '//trim(field_defs(i)%stdname)//'---'//trim(field_defs(i)%shortname)
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
    print *,      subname,' --- IN   --- '

  end subroutine HWRF_AdvertiseFields

 !----------------------------------------------------------------------------------
  subroutine HWRF_FieldsSetup
    integer                     :: rc
    character(len=*),parameter  :: subname='(hwrf_cap:HWRF_FieldsSetup)'

    !--------- import fields to hwrf  -------------
    
    !--------- export fields from hwrf -------------
       call fld_list_add(num=fldsFrWrf_num, fldlist=fldsFrWrf, stdname="air_pressure_at_sea_level" , shortname= "pmsl" )
       call fld_list_add(num=fldsFrWrf_num, fldlist=fldsFrWrf, stdname="inst_merid_wind_height10m" , shortname= "imwh10m" )
       call fld_list_add(num=fldsFrWrf_num, fldlist=fldsFrWrf, stdname="inst_zonal_wind_height10m" , shortname="izwh10m" )

  !
    write(info,*) subname,' --- Passed--- '
    print *,      subname,' --- Passed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)     
  end subroutine HWRF_FieldsSetup


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
    character(len=*), parameter :: subname='(hwrf_cap:fld_list_add)'

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
    print *,      subname,' --- Passed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)     
  end subroutine fld_list_add




 !-----------------------------------------------------------------------------
    !> Called by NUOPC to realize import and export fields.

    !! The fields to import and export are stored in the fldsToAdc and fldsFrAdc
    !! arrays, respectively.  Each field entry includes the standard name,
    !! information about whether the field's grid will be provided by the cap,
    !! and optionally a pointer to the field's data array.  Currently, all fields
    !! are defined on the same grid defined by the cap.
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
        character(len=*),parameter   :: subname='(hwrf_cap:RealizeFieldsProvidingGrid)' 
       
        integer(ESMF_KIND_I4), pointer :: gmask(:,:) , index1(:)
        integer                        :: i,j

        type(ESMF_Grid) :: ModelGrid


        rc = ESMF_SUCCESS

        print *, '<<<<< before  hwrf grid create'
        ModelGrid = CreateGrid_ModelGrid(rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

        
        
      
        !print *, ' min prmsl >>>', minval(PRMSL(1,:,:))
        
        
        !gmask = nint(PRMSL(1,:,:))
        !gmask = where(PRMSL(1,:,:).gt.9.8e20,0,1)
        
        !index1 = minloc(PRMSL, MASK=(PRMSL > 1e15)) 
        
        !print *, 'index1  >>  ',index1(1)
         
         
         
         
         
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! add mask
        !call ESMF_GridGetItem(ModelGrid, &!localDE=0, & 
        !    staggerloc=ESMF_STAGGERLOC_CENTER, &
        !    itemflag=ESMF_GRIDITEM_MASK, farrayPtr=gmask, rc=rc)




        !!TODO: I do not know howto maskout nan values. I want to add a mask array to grid
        !!TODO: however it is not working.
        ALLOCATE(gmask(nlon,nlat))
        do i=1,nlon
            do j=1,nlat
                gmask(i,j) = 0         !not masked
                if (PRMSL(1,i,j) .gt. 1e10) then
                    gmask(i,j)   = 1   !masked
                endif
            enddo
        enddo
        
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        print *, '>>>>> pass  hwrf grid create'

!        height = ESMF_FieldCreate(name="pmsl", &
!            grid=ModelGrid, &
!            farray=PRMSL, &
!            indexflag=ESMF_INDEX_GLOBAL, &
!            rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!            line=__LINE__, &
!            file=__FILE__)) &
!            return  ! bail out
!
!        call NUOPC_Realize(exportState, field=height, rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!            line=__LINE__, &
!            file=__FILE__)) &
!            return  ! bail out


    !call HWRF_RealizeFields(importState, GridIn , mdata, fldsToAdc_num, fldsToAdc, "Adcirc import", rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    !
    call HWRF_RealizeFields(exportState, ModelGrid, fldsFrWrf_num, fldsFrWrf, "hwrf export", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out


    write(info,*) subname,' --- initialization phase 2 completed --- '
    print *,      subname,' --- initialization phase 2 completed --- '
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
  subroutine HWRF_RealizeFields(state, grid, nfields, field_defs, tag, rc)

    type(ESMF_State), intent(inout)             :: state
    type(ESMF_Grid), intent(in)                 :: grid
    integer, intent(in)                         :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    character(len=*), intent(in)                :: tag
    integer, intent(inout)                      :: rc


    type(ESMF_Field)                            :: field
    integer                                     :: i
    character(len=*),parameter  :: subname='(hwrf_cap:HWRF_RealizeFields)'

    rc = ESMF_SUCCESS

    do i = 1, nfields
        field = ESMF_FieldCreate(name=field_defs(i)%shortname, grid=grid, &
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

            print *,      subname, field_defs(i)%stdname ,' --- Connected --- '

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
            print *,      subname," Field ", field_defs(i)%stdname ,' --- Not-Connected --- '

        endif
    enddo

        write(info,*) subname,' --- OUT--- '
        print *,      subname,' --- OUT --- '
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=rc)
  end subroutine HWRF_RealizeFields
  !--------------------------------------


  subroutine SetClock(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_TimeInterval)       :: HWRFTimeStep
    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! initialize internal clock
    ! here: parent Clock and stability timeStep determine actual model timeStep
    !TODO: stabilityTimeStep should be read in from configuation
    !TODO: or computed from internal Grid information
!    call ESMF_TimeIntervalSet(HWRFTimeStep, s=wrf_int, sN=wrf_num, sD=wrf_den, rc=rc) ! 5 minute steps
    call ESMF_TimeIntervalSet(HWRFTimeStep, s=wrf_int, rc=rc) ! 5 minute steps

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetClock(model, clock, HWRFTimeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
  end subroutine



!    subroutine RealizeFieldsProvidingGrid(model, importState, exportState, clock, rc)
!        type(ESMF_GridComp)  :: model
!        type(ESMF_State)     :: importState, exportState
!        type(ESMF_Clock)     :: clock
!        integer, intent(out) :: rc
!!        character(len=*),parameter   :: subname='(hwrf_cap:RealizeFieldsProvidingGrid)'                
!
 !       type(ESMF_Grid) :: ModelGrid
 !       type(ESMF_Field) :: height
!
!        rc = ESMF_SUCCESS
!
!        print *, '<<<<< before  hwrf grid create'
!        ModelGrid = CreateGrid_ModelGrid(rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!!            line=__LINE__, &
! !           file=__FILE__)) &
!            return  ! bail out
!
 !       print *, '>>>>> pass  hwrf grid create'

!        height = ESMF_FieldCreate(name="pmsl", &
!            grid=ModelGrid, &
!            farray=PRMSL, &
!            indexflag=ESMF_INDEX_GLOBAL, &
!            rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!            line=__LINE__, &
!            file=__FILE__)) &
!            return  ! bail out
!
!        call NUOPC_Realize(exportState, field=height, rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!            line=__LINE__, &
!            file=__FILE__)) &
!            return  ! bail out
!!    write(info,*) subname,' --- hwrf SetServices completed --- '
 !   print *,      subname,' --- hwrf SetServices completed --- '
 !   call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)
!    end subroutine


    function CreateGrid_ModelGrid(rc)
        type(ESMF_Grid) :: CreateGrid_ModelGrid
        integer, intent(out), optional :: rc
        character(len=*),parameter   :: subname='(hwrf_cap:CreateGrid_ModelGrid)'        

        real(ESMF_KIND_R8), pointer :: coordX(:,:), coordY(:,:)
        integer :: i, j

        print *, '<0>',nlon,nlat
        rc = ESMF_SUCCESS
        CreateGrid_ModelGrid = ESMF_GridCreateNoPeriDim(name="ModelGrid", &
            minIndex=(/1, 1/), &
            maxIndex=(/nlon, nlat/), &
            indexflag=ESMF_INDEX_GLOBAL, &
            !minCornerCoord=(/1.0_ESMF_KIND_R8, 1.0_ESMF_KIND_R8/), &
            !maxCornerCoord=(/100.0_ESMF_KIND_R8, 100.0_ESMF_KIND_R8/), &
            rc=rc)

        print *, '<1> HWRF LONS ', minval(LONS),maxval(LONS)

        print *, '<1> HWRF LATS ', minval(LATS),maxval(LATS)
        ! add coordinates
        call ESMF_GridAddCoord(CreateGrid_ModelGrid, &
            staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        print *, '<2>'

        call ESMF_GridGetCoord(CreateGrid_ModelGrid, coordDim=1, &
            staggerloc=ESMF_STAGGERLOC_CENTER, &
            farrayPtr=coordX, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        print *, '<3>'

        call ESMF_GridGetCoord(CreateGrid_ModelGrid, coordDim=2, &
            staggerloc=ESMF_STAGGERLOC_CENTER, &
            farrayPtr=coordY, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        print *, '<4>'
        ! set coordinates
        do i=1,nlon
            do j=1,nlat
                coordX(i,j) = LONS(i)
                coordY(i,j) = LATS(j)
            enddo
        enddo

    
    write(info,*) subname,' --- hwrf SetServices completed --- '
    print *,      subname,' --- hwrf SetServices completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)
       
    end function


 !    function CreateGrid_Uniform(rc)
!        type(ESMF_Grid) :: CreateGrid_Uniform
!        integer, intent(out), optional :: rc
!
!        rc = ESMF_SUCCESS
!        CreateGrid_Uniform = ESMF_GridCreateNoPeriDimUfrm(name="ModelGrid", &
!            minIndex=(/1, 1/), &
!            maxIndex=(/sw_nx, sw_ny/), &
!            !indexflag=ESMF_INDEX_GLOBAL, &
!            minCornerCoord=(/0.0_ESMF_KIND_R8, 0.0_ESMF_KIND_R8/), &
!            maxCornerCoord=(/(sw_nx-1)*sw_dx, (sw_ny-1)*sw_dy/), &
!            staggerloclist=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
!            coordSys=ESMF_COORDSYS_CART, &
!            rc=rc)
!
!    end function

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
    type(ESMF_Time)            :: currTime,startTime
    type(ESMF_TimeInterval)    :: timeStep,deltaTime
    character(len=*),parameter :: subname='(hwrf_cap:ModelAdvance)'
    !tmp vector
    !real(ESMF_KIND_R8), pointer:: tmp(:,:)

    ! exports
    real(ESMF_KIND_R8), pointer:: dataPtr_pmsl(:,:)
    real(ESMF_KIND_R8), pointer:: dataPtr_imwh10m(:,:)
    real(ESMF_KIND_R8), pointer:: dataPtr_izwh10m(:,:)

    !imports
    !real(ESMF_KIND_R8), pointer:: dataPtr_x(:)
 

    type(ESMF_StateItem_Flag)  :: itemType
    type(ESMF_Grid)            :: grid
    type(ESMF_Field)           :: lfield
    character(len=128)         :: fldname,timeStr
    integer                    :: i1
    integer                    :: ITIME_BGN_ADC, ITIME_END_ADC
    integer                    :: ITHS, NT, DTDP, ITIME
    integer                    :: nCplWRF
    real(ESMF_KIND_R8)         :: timeStepAbs
    ! local variables for Get methods
    integer :: YY, MM, DD, H, M, S, delta_d, delta_h, delta_m, delta_h_tot
    integer :: ss,ssN,ssD


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
      preString="------>Advancing HWRF from: ", rc=rc)
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
      preString="------------------HWRF-------------> to: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimeGet(currTime, yy=YY, mm=MM, dd=DD, h=H, m=M, s=S, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    print *      , "currTime = ", YY, "/", MM, "/", DD," ", H, ":", M, ":", S


    print *, "HWRF Timeinterval1 = "
    call ESMF_TimeIntervalPrint(timeStep, options="string", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out      

    call ESMF_TimeGet(currTime, timeStringISOFrac=timeStr , rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

        ! Generate file name.
        !Start time
        call ESMF_ClockGet(clock, startTime=startTime, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        call ESMF_TimeGet(startTime, yy=YY, mm=MM, dd=DD, h=H, m=M, s=S, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out        
        
        !to construct HWRF file name we need elapsed time from start date in hours
        deltaTime = currTime - startTime
        call ESMF_TimeIntervalGet (deltaTime, d=delta_d, h=delta_h, m=delta_m) 
        

        !TODO:  consider change this to read new files every hwrf output time steps
        !FILE_NAME="hwrf_data/andrew04l.1992082100.hwrfprs.d123.0p06.f000.grb2.nc"
        !TODO: data assumed to be hourly! check if data provided in min
        delta_h_tot = nint (delta_d * 24.0 + delta_h + delta_m / 60.0)   
        !call update_hwrf_filename (YY, MM, DD, delta_h_tot)
        
        
        
        !!!!  TODO: need to remove
        FILE_NAME =  TRIM(wrf_dir)//'/hwrf_grid.nc'
        print *, FILE_NAME
        call read_hwrf_nc()       
       
        print *, 'min max pmsl', minval(PRMSL(1,:,:)), maxval(PRMSL(1,:,:))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! time managment stuff    
!    ITHS = 0
!    NT   = 100
!    DTDP = 300
!    ITIME= 1
 !   
!    !global values for adcirc time step number
 !   ITIME_BGN_ADC = ITHS + 1   !if hot start then ITHS>0
 !   ITIME_END_ADC = NT         !NT is set in read_input.F
!
 !   call ESMF_TimeIntervalGet(timeStep, s=ss,sN=ssN,sD=ssD,rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
 !       line=__LINE__, &
 !       file=__FILE__)) &
!        return  ! bail out
!    timeStepAbs = real(ss) + real(ssN/ssD)
!    if (mod(timeStepAbs,DTDP) .EQ. 0) then
!      nCplWRF = int(timeStepAbs/DTDP)
!    else
!      nCplWRF = int(timeStepAbs/DTDP)+1
!    endif

!    write(info,*) subname,' --- run phase 2 called --- '
 !   print *,      subname,' --- run phase 2 called --- '
 !   call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

    !CALL ADCIRC_Run(nCplWRF)

 !   write(info,*) subname,' --- run phase 3 called ---  nCplWRF = ',nCplWRF
!    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    !-----------------------------------------
    !   IMPORT
    !-----------------------------------------
    !Get and fill imported fields


    !-----------------------------------------
    !   EXPORT
    !-----------------------------------------
    !pack and send exported fields
    !allocate (tmp(nlon, nlat))

    ! >>>>> PACK and send ZETA
    call State_getFldPtr(ST=exportState,fldname='pmsl',fldptr=dataPtr_pmsl, &
      rc=rc,dump=.true.,timeStr=timeStr)


    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !assign to field
    dataPtr_pmsl(:,:) = PRMSL(1,:,:) + 33333.0
    
    
    print *, 'size(PRMSL(1,:,:))  >>',  size(PRMSL (:,:,1))
    print *, 'size(UGRD10(1,:,:))  >>', size(UGRD10(:,:,1))
    !----------------------------------------
    ! >>>>> PACK and send VELX
    !call State_getFldPtr(ST=exportState,fldname='izwh10m',fldptr=dataPtr_izwh10m,rc=rc)
    call State_getFldPtr(ST=exportState,fldname='izwh10m',fldptr=dataPtr_izwh10m, &
      rc=rc,dump=.true.,timeStr=timeStr)  
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    !assign to field
    dataPtr_izwh10m(:,:) =  UGRD10(:,:,1) * 0.0 + 10.0 
    !----------------------------------------
    ! >>>>> PACK and send VELY
    !call State_getFldPtr(ST=exportState,fldname='imwh10m',fldptr=dataPtr_imwh10m,rc=rc)
    call State_getFldPtr(ST=exportState,fldname='imwh10m',fldptr=dataPtr_imwh10m, &
      rc=rc,dump=.true.,timeStr=timeStr)    
    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !assign to field
    dataPtr_imwh10m(:,:) = VGRD10(:,:,1) * 0.0 + 5.0
    !----------------------------------------
!
!
    !! TODO:  not a right thing to do. we need to fix the grid object mask <<<<<<
    where(dataPtr_pmsl.gt.9.8e10) dataPtr_pmsl = 0.0
    where(dataPtr_imwh10m.gt.9.8e10) dataPtr_imwh10m = 0.0      !v-y-comp
    where(dataPtr_izwh10m.gt.9.8e10) dataPtr_izwh10m = 0.0      !u-x-comp

  end subroutine

!-----------------------------------------------------------
  !> Retrieve a pointer to a field's data array from inside an ESMF_State object.
  !!
  !! @param ST the ESMF_State object
  !! @param fldname name of the fields
  !! @param fldptr pointer to 1D array
  !! @param rc return code
  subroutine State_GetFldPtr(ST, fldname, fldptr, rc, dump,timeStr)
    type(ESMF_State), intent(in) :: ST
    character(len=*), intent(in) :: fldname
    real(ESMF_KIND_R8), pointer, intent(in) :: fldptr(:,:)
    integer, intent(out), optional :: rc
    logical, intent(in), optional  :: dump
    character(len=128),intent(inout), optional :: timeStr
    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
    character(len=*),parameter :: subname='(hwrf_cap:State_GetFldPtr)'

    call ESMF_StateGet(ST, itemName=trim(fldname), field=lfield, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (present(rc)) rc = lrc

    if (dump) then
       if (.not. present(timeStr)) timeStr="_"
        call ESMF_FieldWrite(lfield, &
        fileName='field_hwrf_'//trim(fldname)//trim(timeStr)//'.nc', &
        rc=rc,overwrite=.true.)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
    endif
  end subroutine State_GetFldPtr


  !-----------------------------------------------------------------------------
  !> Called by NUOPC at the end of the run to clean up.  The cap does
  !! this simply by calling ADCIRC_Final.
  !!
  !! @param gcomp the ESMF_GridComp object
  !! @param rc return code
    subroutine HWRF_model_finalize(gcomp, rc)

        ! input arguments
        type(ESMF_GridComp)  :: gcomp
        integer, intent(out) :: rc

        ! local variables
        type(ESMF_Clock)     :: clock
        type(ESMF_Time)                        :: currTime
        character(len=*),parameter  :: subname='(adc_cap:adc_model_finalize)'

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

        !CALL ADCIRC_Final()

        write(info,*) subname,' --- hwrf finalize completed --- '
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

    end subroutine HWRF_model_finalize





end module


