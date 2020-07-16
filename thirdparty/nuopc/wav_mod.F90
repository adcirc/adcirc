!>
!! @mainpage ADCIRC NUOPC Cap
!! @author Saeed Moghimi (moghimis@gmail.com)
!! @date 15/1/17 Original documentation
!------------------------------------------------------
!LOG-----------------
!
!
!

module wav_mod

  !-----------------------------------------------------------------------------
  ! ADCIRC mesh utility
  !-----------------------------------------------------------------------------
  use mpi
  use ESMF
  use NUOPC
  USE netcdf

  implicit none
    ! reading data time management info WW3 <-----> ADC exchange
    integer                  :: ww3_cpl_int,ww3_cpl_num,ww3_cpl_den
    logical                  :: ww3_from_file
    character (len = 280)    :: wav_dir, wav_nam !, ww3_grd
    character (len = 280)    :: FILE_NAME
    character (len =2048)    :: info

    ! info for reading ww3 netcdf file
    integer               :: nnode,nelem , ntime, noel
    real(ESMF_KIND_R8), allocatable     :: LONS(:), LATS(:),TIMES(:)
    integer           , allocatable     :: TRI(:,:)
    real(ESMF_KIND_R8), allocatable     :: SXX (:,:), SYY(:,:), SXY (:,:)

    ! netcdf vars
    integer :: ncid, NOD_dimid, rec_dimid, ELM_dimid, NOE_dimid
    integer :: LON_varid, LAT_varid, rec_varid, tri_varid
    integer :: SXX_varid, SYY_varid, SXY_varid

    !> \author Ali Samii - 2016
    !! See: https://github.com/samiiali
    !! \brief This object stores the data required for construction of a parallel or serial
    !! ESMF_Mesh from <tt>fort.14, fort.18, partmesh.txt</tt> files.
    !!
   type meshdata
        !> \details vm is an ESMF_VM object.  ESMF_VM is just an ESMF virtual machine class,
        !! which we will use to get the data about the local PE and PE count.
        type(ESMF_VM)                      :: vm
        !> \details This array contains the node coordinates of the mesh. For
        !! example, in a 2D mesh, the \c jth coordinate of the \c nth node
        !! is stored in location <tt> 2*(n-1)+j</tt> of this array.
        real(ESMF_KIND_R8), allocatable    :: NdCoords(:)
        !> \details This array contains the elevation of different nodes of the mesh
        real(ESMF_KIND_R8), allocatable    :: bathymetry(:)
        !> \details Number of nodes present in the current PE. This is different from the
        !! number of nodes owned by this PE (cf. NumOwnedNd)
        integer(ESMF_KIND_I4)              :: NumNd
        !> \details Number of nodes owned by this PE. This is different from the number of
        !! nodes present in the current PE (cf. NumNd)
        integer(ESMF_KIND_I4)              :: NumOwnedNd
        !> \details Number of elements in the current PE. This includes ghost elements and
        !! owned elements. However, we do not bother to distinguish between owned
        !! element and present element (as we did for the nodes).
        integer(ESMF_KIND_I4)              :: NumEl
        !> \details Number of nodes of each element, which is simply three in 2D ADCIRC.
        integer(ESMF_KIND_I4)              :: NumND_per_El
        !> \details Global node numbers of the nodes which are present in the current PE.
        integer(ESMF_KIND_I4), allocatable :: NdIDs(:)
        !> \details Global element numbers which are present in the current PE.
        integer(ESMF_KIND_I4), allocatable :: ElIDs(:)
        !> \details The element connectivity array, for the present elements in the current PE.
        !! The node numbers are the local numbers of the present nodes. All the element
        !! connectivities are arranged in this one-dimensional array.
        integer(ESMF_KIND_I4), allocatable :: ElConnect(:)
        !> \details The number of the PE's which own each of the nodes present this PE.
        !! This number is zero-based.
        integer(ESMF_KIND_I4), allocatable :: NdOwners(:)
        !> \details An array containing the element types, which are all triangles in our
        !! application.
        integer(ESMF_KIND_I4), allocatable :: ElTypes(:)
        !> \details This is an array, which maps the indices of the owned nodes to the indices of the present
        !! nodes. For example, assume we are on <tt>PE = 1</tt>, and we have four nodes present, and the
        !! first and third nodes belong to <tt>PE = 0</tt>. So we have:
        !! \code
        !! NumNd = 4
        !! NumOwnedNd = 2
        !! NdOwners = (/0, 1, 0, 1/)
        !! NdIDs = (/2, 3, 5, 6/)
        !! owned_to_present = (/2, 4/)    <-- Because the first node owned by this PE is actually
        !!                                    the second node present on this PE, and so on.
        !! \endcode
        integer(ESMF_KIND_I4), allocatable :: owned_to_present_nodes(:)
    end type meshdata

!-----------------------------------------------------------------------------
  contains

!-----------------------------------------------------------------------
!- Sub !!!????
!-----------------------------------------------------------------------
    SUBROUTINE init_ww3_nc()
      IMPLICIT NONE
      character (len = *), parameter :: NOD_NAME    = "node"
      character (len = *), parameter :: NOE_NAME    = "noel"
      character (len = *), parameter :: ELM_NAME    = "element"
      character (len = *), parameter :: LAT_NAME    = "latitude"
      character (len = *), parameter :: LON_NAME    = "longitude"
      character (len = *), parameter :: REC_NAME    = "time"
      character (len = *), parameter :: SXX_NAME    = "sxx"
      character (len = *), parameter :: SYY_NAME    = "syy"
      character (len = *), parameter :: SXY_NAME    = "sxy"
      character (len = *), parameter :: TRI_NAME    = "tri"
      character (len = 140)          :: units
      character (len = *), parameter :: subname='(wav_mod:init_ww3_nc)'

      logical :: THERE
      integer :: lat, lon,i, iret, rc

      FILE_NAME =  TRIM(wav_dir)//'/'//TRIM(wav_nam)
      print *, ' FILE_NAME  > ', FILE_NAME
      INQUIRE( FILE= FILE_NAME, EXIST=THERE )
      if ( .not. THERE)  stop 'WW3 netcdf file does not exist!'

      ncid = 0
      ! Open the file.
      call check(  nf90_open(trim(FILE_NAME), NF90_NOWRITE, ncid))

      ! Get ID of unlimited dimension
      call check( nf90_inquire(ncid, unlimitedDimId = rec_dimid) )

      ! Get ID of limited dimension
      call check( nf90_inq_dimid(ncid, NOD_NAME, NOD_dimid) )
      call check( nf90_inq_dimid(ncid, ELM_NAME, ELM_dimid) )
      call check( nf90_inq_dimid(ncid, NOE_NAME, NOE_dimid) )

      ! How many values of "nodes" are there?
      call check(nf90_inquire_dimension(ncid, NOD_dimid, len = nnode) )
      call check(nf90_inquire_dimension(ncid, ELM_dimid, len = nelem) )
      call check(nf90_inquire_dimension(ncid, NOE_dimid, len = noel) )
      ! What is the name of the unlimited dimension, how many records are there?
      call check(nf90_inquire_dimension(ncid, rec_dimid, len = ntime))

      !print *,  ' nelem  > ',nelem , ' noel  > ' ,noel,  ' ntime > ',ntime

      ! Get the varids of the pressure and temperature netCDF variables.
      call check( nf90_inq_varid(ncid, LAT_NAME,    LAT_varid) )
      call check( nf90_inq_varid(ncid, LON_NAME,    LON_varid) )
      call check( nf90_inq_varid(ncid, REC_NAME,    rec_varid) )
      call check( nf90_inq_varid(ncid, SXX_NAME,    SXX_varid) )
      call check( nf90_inq_varid(ncid, SYY_NAME,    SYY_varid) )
      call check( nf90_inq_varid(ncid, SXY_NAME,    SXY_varid) )
      call check( nf90_inq_varid(ncid, TRI_NAME,    TRI_varid) )

      !allocate vars
      if(.not. allocated(LATS))  allocate (LATS  (1:nnode))
      if(.not. allocated(LONS))  allocate (LONS  (1:nnode))
      if(.not. allocated(TIMES)) allocate (TIMES (1:ntime))
      if(.not. allocated(TRI))   allocate (TRI (1:noel ,1:nelem))
      ! read vars
      call check(nf90_get_var(ncid, LAT_varid, LATS ))
      call check(nf90_get_var(ncid, LON_varid, LONS ))
      call check(nf90_get_var(ncid, rec_varid, TIMES))
      !call check(nf90_get_var(ncid, SXX_varid, SXX  ))
      !TODO: Why the order is other way???? Might change the whole forcing fields!!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< IMPORTANT <<<<<
      ! plot input and out put to be sure we are not scrambling the data. the same for HWRF netcdf file
      call check(nf90_get_var(ncid, TRI_varid, TRI, start = (/1,1/),count = (/noel,nelem/)  ))

      write(info,*) subname,' --- init ww3 netcdf file  --- '
      !print *, info
      call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)

    END SUBROUTINE

!-----------------------------------------------------------------------
!- Sub !!!????
!-----------------------------------------------------------------------
    SUBROUTINE read_ww3_nc(currTime)
      IMPLICIT NONE
      type(ESMF_Time),intent(in)     :: currTime
      type(ESMF_Time)                :: refTime
      type(ESMF_TimeInterval)        :: dTime


      character (len = 140)          :: units

      integer, parameter :: NDIMS = 2
      integer    :: start(NDIMS),count(NDIMS)
      logical    :: THERE
      real       :: delta_d_all (ntime) , delta_d_ref
      !integer   :: dimids(NDIMS)
      character (len = *), parameter :: subname='(wav_mod:read_ww3_nc)'
      character  :: c1,c2,c3,c4,c5,c6,c7
      integer    :: yy,mm,dd,hh,min,ss
      integer    :: d_d, d_h, d_m, d_s
      integer    :: lat, lon,it, iret, rc

      rc = ESMF_SUCCESS

      !units = "days since 1990-01-01 00:00:00"
      call check(nf90_get_att(ncid,rec_varid,'units',units))
      READ(units,'(a4,a7,i4,a1,i2,a1,i2,a1,i2,a1,i2,a1,i2)',iostat=iret)  &
                   c1,c2,yy,c3,mm,c4,dd,c5,hh,c6,min,c7,ss

      if (iret .ne. 0) then
        print *, 'Fatal error: A non valid time units string was provided'
        stop 'wav_mod: read_ww3_nc'
      end if

      call ESMF_TimeSet(refTime, yy=yy, mm=mm, dd=dd, h=hh, m=min, s=ss, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      dTime = currTime - refTime
      call ESMF_TimeIntervalGet (dTime, d=d_d, h=d_h, m=d_m, s=d_s, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      delta_d_ref = d_d + d_h /24.0 + d_m / (24.0 * 60.0) +  d_s / (24.0 * 3600.0)

      do it = 1, ntime
        delta_d_all(it) =  delta_d_ref - TIMES (it)
      end do

      it = minloc(abs(delta_d_all),dim=1)

      !if (abs(delta_d_all (it)) .gt. 7200.) then
      !   write(info,*) subname,' --- STOP WAVcap: Time dif is GT 2 hours ---  '
      !   call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)
      !   stop                  ' --- STOP WAVcap: Time dif is GT 2 hours ---  '
      !endif

      !it = 1
      !print *, 'wave file it index > ',it
      write(info,*) subname,'WAV file it index > ',it
      !print *, info
      call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)

      
      call ESMF_TimePrint(refTime, preString="WAV refTime=  ", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out


      !alocate vars
      if(.not. allocated(SXX))   allocate (SXX (1:nnode,1))
      if(.not. allocated(SYY))   allocate (SYY (1:nnode,1))
      if(.not. allocated(SXY))   allocate (SXY (1:nnode,1))

      start = (/ 1   , it/)
      count = (/nnode, 1 /)  !for some reason the order here is otherway around?!

      !print *, start+count
      !print *,size(SXX(ntime,:))
      call check( nf90_get_var(ncid,SXX_varid, SXX, start, count) )
      call check( nf90_get_var(ncid,SYY_varid, SYY, start, count) )
      call check( nf90_get_var(ncid,SXY_varid, SXY, start, count) )

      !print *,FILE_NAME , '   HARD CODED for NOWWWW>>>>>     Time index from wave file is > ', it, SXX(1:10,1)

      write(info,*) subname,' --- read ww3 netcdf file  --- '
      !print *, info
      call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)

    !
    END SUBROUTINE


    subroutine construct_meshdata_from_netcdf(the_data)
        implicit none
        integer                               :: i1
        integer, parameter                    :: dim1=2, spacedim=2, NumND_per_El=3
        type(meshdata), intent(inout)         :: the_data
        the_data%NumEl = nelem
        the_data%NumNd = nnode
        allocate(the_data%NdIDs(the_data%NumNd))
        allocate(the_data%ElIDs(the_data%NumEl))
        allocate(the_data%NdCoords(dim1*the_data%NumNd))
        allocate(the_data%bathymetry(the_data%NumNd))
        allocate(the_data%ElConnect(NumND_per_El*the_data%NumEl))
        allocate(the_data%NdOwners(the_data%NumNd))
        allocate(the_data%ElTypes(the_data%NumEl))
        allocate(the_data%owned_to_present_nodes(the_data%NumNd))

        do i1 = 1, the_data%NumNd, 1
                the_data%NdIDs(i1)                 = i1
                the_data%NdCoords((i1-1)*dim1 + 1) = LONS(i1)
                the_data%NdCoords((i1-1)*dim1 + 2) = LATS(i1)
        end do
        do i1 = 1, the_data%NumEl, 1
                the_data%ElIDs(i1)                        =   i1
                the_data%ElConnect((i1-1)*NumND_per_El+1) = TRI(1,i1)
                the_data%ElConnect((i1-1)*NumND_per_El+2) = TRI(2,i1)
                the_data%ElConnect((i1-1)*NumND_per_El+3) = TRI(3,i1)
        end do
        !We have only one node therefore:
        the_data%NdOwners = 0                  !process 0 owns all the nodes
        the_data%NumOwnedND = the_data%NumNd   !number of nodes = number of owned nodes
        the_data%owned_to_present_nodes = the_data%NdIDs

        the_data%ElTypes = ESMF_MESHELEMTYPE_TRI

        close(14)
    end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop
    end if
  end subroutine check


!-------------------------------------------------------------------------
    FUNCTION Replace_Text (s,text,rep)  RESULT(outs)
    CHARACTER(*)        :: s,text,rep
    CHARACTER(LEN(s)+300) :: outs     ! provide outs with extra char len
    INTEGER             :: i, nt, nr

    outs = s ; nt = LEN_TRIM(text) ; nr = LEN_TRIM(rep)
        DO
           i = INDEX(outs,text(:nt)) ; IF (i == 0) EXIT
           outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
        END DO
    END FUNCTION Replace_Text


    subroutine update_ww3_filename (YY, MM, DD, H)
        integer             :: YY, MM, DD, H
        CHARACTER(len=280)      :: inps     ! provide outs with extra 100 char len
        CHARACTER(len=4)        :: year
        CHARACTER(len=2)        :: mon,day
        CHARACTER(len=3)        :: hours

        ! example:  wav_nam: ww3.Constant.YYYYMMDD_sxy.nc
        inps = trim(wav_nam)

        write(year,"(I4.4)") YY
        inps =  Replace_Text (inps,'YYYY',year)

        write(mon,"(I2.2)")  MM
        inps =  Replace_Text (inps,'MM',mon)

        write(day,"(I2.2)")  DD
        inps =  Replace_Text (inps,'DD',day)

        !past hours from start date
        !write(hours,"(I3.3)") H
        !inps =  Replace_Text (inps,'HHH',hours)

        FILE_NAME =  TRIM(wav_dir)//'/'//TRIM(inps)

    END subroutine update_ww3_filename



    !> \author Ali Samii - 2016
    !! See: https://github.com/samiiali
    !> @details Using the data available in <tt> fort.14, fort.18, partmesh.txt</tt> files
    !! this function extracts the scalars and arrays required for construction of a
    !! meshdata object.
    !! After calling this fucntion, one can call create_parallel_esmf_mesh_from_meshdata()
    !! or create_masked_esmf_mesh_from_data() to create an ESMF_Mesh.
    !! @param vm This is an ESMF_VM object, which will be used to obtain the \c localPE
    !! and \c peCount of the \c MPI_Communicator.
    !! @param global_fort14_dir This is the directory path (relative to the executable
    !! or an absolute path) which contains the global \c fort.14 file (not the fort.14
    !! after decomposition).
    !! @param the_data This is the output meshdata object.
    !!

    !> \details As the name of this function suggests, this funciton creates a parallel
    !! ESMF_Mesh from meshdata object. This function should be called collectively by
    !! all PEs for the parallel mesh to be created. The function, extract_parallel_data_from_mesh()
    !! should be called prior to calling this function.
    !! \param the_data This the input meshdata object.
    !! \param out_esmf_mesh This is the ouput ESMF_Mesh object.
    subroutine create_parallel_esmf_mesh_from_meshdata(the_data, out_esmf_mesh)
        implicit none
        type(ESMF_Mesh), intent(out)                  :: out_esmf_mesh
        type(meshdata), intent(in)                    :: the_data
        integer, parameter                            :: dim1=2, spacedim=2, NumND_per_El=3
        integer                                       :: rc
        out_esmf_mesh=ESMF_MeshCreate(parametricDim=dim1, spatialDim=spacedim, &
            nodeIDs=the_data%NdIDs, nodeCoords=the_data%NdCoords, &
            nodeOwners=the_data%NdOwners, elementIDs=the_data%ElIDs, &
            elementTypes=the_data%ElTypes, elementConn=the_data%ElConnect, &
            rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

    end subroutine


    subroutine read_config()
    character(ESMF_MAXPATHLEN)    :: fname ! config file name
    type(ESMF_Config)             :: cf     ! the Config itself
    integer                       :: rc

    rc = ESMF_SUCCESS

   !Initiate reading resource file
    cf = ESMF_ConfigCreate(rc=rc)  ! Create the empty Config
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    fname = "config.rc" ! Name the Resource File
    call ESMF_ConfigLoadFile(cf, fname, rc=rc) ! Load the Resource File
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

   ! This subroutine is not used with NEMS system for time info reading. 
   ! Because the time interval information is passed via nems.configure file 
   ! with time slot definitation.
   
   ! read time coupling interval info
   ! call ESMF_ConfigGetAttribute(cf, ww3_cpl_int, label="cpl_int:",default=300, rc=rc)
   ! call ESMF_ConfigGetAttribute(cf, ww3_cpl_num, label="cpl_num:",default=0  , rc=rc)
   ! call ESMF_ConfigGetAttribute(cf, ww3_cpl_den, label="cpl_den:",default=1  , rc=rc)
   ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
   !   line=__LINE__, &
   !   file=__FILE__)) &
   !   return  ! bail out

   ! From file or couple to model
    call ESMF_ConfigGetAttribute(cf, ww3_from_file, label="ww3_from_file:",default=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ConfigGetAttribute(cf, wav_dir, label="wav_dir:",default='ww3_inp/'  , rc=rc)
    call ESMF_ConfigGetAttribute(cf, wav_nam, label="wav_nam:", &
         default='ww3.Constant.YYYYMMDD_sxy.nc'  , rc=rc)
    !call ESMF_ConfigGetAttribute(cf, ww3_grd, label="ww3_grd:",default='ww3_grid.nc/'  , rc=rc)


    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
        
    call ESMF_ConfigDestroy(cf, rc=rc) ! Destroy the Config
        
    end subroutine read_config







!
!
!
!    !
!    !> \author Ali Samii - 2016
!    !! See: https://github.com/samiiali
!    !> @details Using the data available in <tt> fort.14, fort.18, partmesh.txt</tt> files
!    !! this function extracts the scalars and arrays required for construction of a
!    !! meshdata object.
!    !! After calling this fucntion, one can call create_parallel_esmf_mesh_from_meshdata()
!    !! or create_masked_esmf_mesh_from_data() to create an ESMF_Mesh.
!    !! @param vm This is an ESMF_VM object, which will be used to obtain the \c localPE
!    !! and \c peCount of the \c MPI_Communicator.
!    !! @param global_fort14_dir This is the directory path (relative to the executable
!    !! or an absolute path) which contains the global \c fort.14 file (not the fort.14
!    !! after decomposition).
!    !! @param the_data This is the output meshdata object.
!    !!
!    subroutine extract_parallel_data_from_mesh(global_fort14_dir, the_data,localPet)
!        implicit none
!        type(meshdata), intent(inout)         :: the_data
!        character(len=*), intent(in)          :: global_fort14_dir
!        character(len=200)                    :: partmesh_filename
!        integer, intent(in)                   :: localPet
!        integer                               :: i1, j1, i_num, num_global_nodes,io,garbage2
!        integer, allocatable                  :: local_node_numbers(:), local_elem_numbers(:), node_owner(:)
!        integer, parameter                    :: dim1=2, NumND_per_El=3
!
!    print *,"WAVm ..1.............................................. >> "
!        the_data%NumNd = np
!        the_data%NumEl = ne
!        allocate(the_data%NdIDs(the_data%NumNd))
!        allocate(local_node_numbers(the_data%NumNd))
!        allocate(the_data%ElIDs(the_data%NumEl))
!        allocate(local_elem_numbers(the_data%NumEl))
!        allocate(the_data%NdCoords(dim1*the_data%NumNd))
!        allocate(the_data%bathymetry(the_data%NumNd))
!        allocate(the_data%ElConnect(NumND_per_El*the_data%NumEl))
!        allocate(the_data%NdOwners(the_data%NumNd))
!        allocate(the_data%ElTypes(the_data%NumEl))
!    print *,"WAVm ..2.............................................. >> "
!        local_elem_numbers = IMAP_EL_LG
!        the_data%ElIDs = abs(local_elem_numbers)
!        local_node_numbers = NODES_LG
!        the_data%NumOwnedND = 0
!    print *,"WAVm ..3.............................................. >> "
!        do i1 = 1, the_data%NumNd, 1
!            if (local_node_numbers(i1) > 0) then
!                the_data%NumOwnedND = the_data%NumOwnedND + 1
!            end if
!        end do
!        the_data%NdIDs = abs(local_node_numbers)
!        allocate(the_data%owned_to_present_nodes(the_data%NumOwnedND))
!    print *,"WAVm ..4............................................. >> "
!        !> @details Read partmesh file to get global node information
!        !print *, 'size local_elem_numbers', size(local_elem_numbers),size(IMAP_EL_LG)
!        !print *, 'size local_node_numbers', size(local_node_numbers),size(NODES_LG)
!
!        !partmesh_filename = trim(global_fort14_dir//"/partmesh.txt")
!        open(unit=10099, file = TRIM(global_fort14_dir)//'/'//'partmeshw.txt', &
!            form='FORMATTED', status='OLD', action='READ')
!
!        !! Very ugly way of finding global element numebr
!        !! TODO: Saeed: need to find the reprsentive value for this
!        num_global_nodes = 0
!        do
!            read(unit=10099,fmt=*,iostat=io) garbage2
!            if (io/=0) exit
!            num_global_nodes = num_global_nodes + 1
!        end do
!        rewind(unit=10099)
!        !
!        !print *, 'size num_global_nodes', num_global_nodes,localPet
!
!        allocate(node_owner(num_global_nodes))
!        read(unit=10099, fmt=*) node_owner
!        close(10099)
!
!        do i1 = 1, the_data%NumNd, 1
!                the_data%NdCoords((i1-1)*dim1 + 1) = slam(i1)
!                the_data%NdCoords((i1-1)*dim1 + 2) = sfea(i1)
!        end do
!        do i1 = 1, the_data%NumEl, 1
!                the_data%ElConnect((i1-1)*NumND_per_El+1) = nm (i1,1)
!                the_data%ElConnect((i1-1)*NumND_per_El+2) = nm (i1,2)
!                the_data%ElConnect((i1-1)*NumND_per_El+3) = nm (i1,3)
!        end do
!
!        do i1= 1, the_data%NumNd, 1
!            the_data%NdOwners(i1) = node_owner(the_data%NdIDs(i1)) - 1
!        end do
!
!        j1 = 0
!        do i1 = 1, the_data%NumNd, 1
!            if (the_data%NdOwners(i1) == localPet) then
!                j1 = j1 + 1
!                the_data%owned_to_present_nodes(j1) = i1
!            end if
!        end do
!        the_data%ElTypes = ESMF_MESHELEMTYPE_TRI
!
!        !TODO: Saeed: Check if I need to dealocate arrays here!
!
!
!    end subroutine extract_parallel_data_from_mesh
!!
!
!    subroutine extract_parallel_data_from_mesh_orig(global_fort14_dir, the_data,localPet)
!        implicit none
!
!        type(meshdata), intent(inout)         :: the_data
!        character(len=*), intent(in)          :: global_fort14_dir
!        integer, intent(in)                   :: localPet
!        character(len=6)                      :: PE_ID, garbage1
!        
!        character(len=200)                    :: fort14_filename, fort18_filename, partmesh_filename
!        integer                               :: i1, j1, i_num, petCount, num_global_nodes, garbage2, garbage3
!        integer, allocatable                  :: local_node_numbers(:), local_elem_numbers(:), node_owner(:)
!        integer, parameter                    :: dim1=2, NumND_per_El=3
!
!        write(PE_ID, "(A,I4.4)") "WE", localPet
!        fort14_filename = TRIM(global_fort14_dir)//'/'//PE_ID//"/fort.14"
!        fort18_filename = TRIM(global_fort14_dir)//'/'//PE_ID//"/fort.18"
!        partmesh_filename = TRIM(global_fort14_dir)//'/'//'partmeshw.txt'
!
!
!        
!
!        open(unit=23414, file=fort14_filename, form='FORMATTED', status='OLD', action='READ')
!        open(unit=23418, file=fort18_filename, form='FORMATTED', status='OLD', action='READ')
!        open(unit=234100, file=partmesh_filename, form='FORMATTED', status='OLD', action='READ')
!
!        read(unit=23414, fmt=*)
!        read(unit=23414, fmt=*) the_data%NumEl, the_data%NumNd
!        allocate(the_data%NdIDs(the_data%NumNd))
!        allocate(local_node_numbers(the_data%NumNd))
!        allocate(the_data%ElIDs(the_data%NumEl))
!        allocate(local_elem_numbers(the_data%NumEl))
!        allocate(the_data%NdCoords(dim1*the_data%NumNd))
!        allocate(the_data%bathymetry(the_data%NumNd))
!        allocate(the_data%ElConnect(NumND_per_El*the_data%NumEl))
!        allocate(the_data%NdOwners(the_data%NumNd))
!        allocate(the_data%ElTypes(the_data%NumEl))
!
!        read(unit=23418, fmt=*)
!        read(unit=23418, fmt=*)
!        read(unit=23418, fmt=*) local_elem_numbers
!        the_data%ElIDs = abs(local_elem_numbers)
!        read(unit=23418, fmt=*) garbage1, num_global_nodes, garbage2, garbage3
!        read(unit=23418, fmt=*) local_node_numbers
!        the_data%NumOwnedND = 0
!        do i1 = 1, the_data%NumNd, 1
!            if (local_node_numbers(i1) > 0) then
!                the_data%NumOwnedND = the_data%NumOwnedND + 1
!            end if
!        end do
!        the_data%NdIDs = abs(local_node_numbers)
!        allocate(node_owner(num_global_nodes))
!        allocate(the_data%owned_to_present_nodes(the_data%NumOwnedND))
!        read(unit=234100, fmt=*) node_owner
!
!        do i1 = 1, the_data%NumNd, 1
!            read(unit=23414, fmt=*) local_node_numbers(i1), &
!                the_data%NdCoords((i1-1)*dim1 + 1), &
!                the_data%NdCoords((i1-1)*dim1 + 2), &
!                the_data%bathymetry(i1)
!        end do
!        do i1 = 1, the_data%NumEl, 1
!            read(unit=23414, fmt=*) local_elem_numbers(i1), i_num, &
!                the_data%ElConnect((i1-1)*NumND_per_El+1), &
!                the_data%ElConnect((i1-1)*NumND_per_El+2), &
!                the_data%ElConnect((i1-1)*NumND_per_El+3)
!        end do
!
!        do i1= 1, the_data%NumNd, 1
!            the_data%NdOwners(i1) = node_owner(the_data%NdIDs(i1)) - 1
!        end do
!
!        j1 = 0
!        do i1 = 1, the_data%NumNd, 1
!            if (the_data%NdOwners(i1) == localPet) then
!                j1 = j1 + 1
!                the_data%owned_to_present_nodes(j1) = i1
!            end if
!        end do
!        the_data%ElTypes = ESMF_MESHELEMTYPE_TRI
!
!        close(23414)
!        close(23418)
!        close(234100)
!    end subroutine extract_parallel_data_from_mesh_orig
!


end module
