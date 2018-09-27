!>
!! @mainpage ADCIRC NUOPC Cap
!! @author Saeed Moghimi (moghimis@gmail.com)
!! @date 15/1/17 Original documentation
!------------------------------------------------------
!LOG-----------------
!
!
!

module adc_mod

  !-----------------------------------------------------------------------------
  ! ADCIRC mesh utility
  !-----------------------------------------------------------------------------
  use mpi
  use ESMF
  use NUOPC

  use MESH    , only: np,ne,nm,slam,sfea
  use GLOBAL  , only: IMAP_EL_LG,NODES_LG
  use GLOBAL  , only: ETA2, UU2, VV2  ! Export water level and velocity fileds to wave model
  USE GLOBAL  , only: RSNX2, RSNY2    ! Import wave 2D forces from wave model
  use SIZES   , only: ROOTDIR


  implicit none
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

    ! reading data time management info WW3 <-----> ADC exchange
    integer                       :: adc_cpl_int,adc_cpl_num,adc_cpl_den

  !-----------------------------------------------------------------------------
  contains
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
    !
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
    subroutine extract_parallel_data_from_mesh(global_fort14_dir, the_data,localPet)
        implicit none
        type(meshdata), intent(inout)         :: the_data
        character(len=*), intent(in)          :: global_fort14_dir
        character(len=200)                    :: partmesh_filename
        integer, intent(in)                   :: localPet
        integer                               :: i1, j1, i_num, num_global_nodes,io,garbage2
        integer, allocatable                  :: local_node_numbers(:), local_elem_numbers(:), node_owner(:)
        integer, parameter                    :: dim1=2, NumND_per_El=3

        the_data%NumNd = np
        the_data%NumEl = ne
        allocate(the_data%NdIDs(the_data%NumNd))
        allocate(local_node_numbers(the_data%NumNd))
        allocate(the_data%ElIDs(the_data%NumEl))
        allocate(local_elem_numbers(the_data%NumEl))
        allocate(the_data%NdCoords(dim1*the_data%NumNd))
        allocate(the_data%bathymetry(the_data%NumNd))
        allocate(the_data%ElConnect(NumND_per_El*the_data%NumEl))
        allocate(the_data%NdOwners(the_data%NumNd))
        allocate(the_data%ElTypes(the_data%NumEl))

        local_elem_numbers = IMAP_EL_LG
        the_data%ElIDs = abs(local_elem_numbers)
        local_node_numbers = NODES_LG
        the_data%NumOwnedND = 0

        do i1 = 1, the_data%NumNd, 1
            if (local_node_numbers(i1) > 0) then
                the_data%NumOwnedND = the_data%NumOwnedND + 1
            end if
        end do
        the_data%NdIDs = abs(local_node_numbers)
        allocate(the_data%owned_to_present_nodes(the_data%NumOwnedND))

        !> @details Read partmesh file to get global node information
        !print *, 'size local_elem_numbers', size(local_elem_numbers),size(IMAP_EL_LG)
        !print *, 'size local_node_numbers', size(local_node_numbers),size(NODES_LG)

        !partmesh_filename = trim(global_fort14_dir//"/partmesh.txt")
        open(unit=10099, file = TRIM(global_fort14_dir)//'/'//'partmesh.txt', &
            form='FORMATTED', status='OLD', action='READ')

        !! Very ugly way of finding global element numebr
        !! TODO: Saeed: need to find the reprsentive value for this
        num_global_nodes = 0
        do
            read(unit=10099,fmt=*,iostat=io) garbage2
            if (io/=0) exit
            num_global_nodes = num_global_nodes + 1
        end do
        rewind(unit=10099)
        !
        !print *, 'size num_global_nodes', num_global_nodes,localPet

        allocate(node_owner(num_global_nodes))
        read(unit=10099, fmt=*) node_owner
        close(10099)


        !print *, 'ADC SLAM > ',size(slam),'>>>>>>',slam
        !print *, 'ADC SFEA > ',size(sfea),'>>>>>>',sfea
        
        do i1 = 1, the_data%NumNd, 1
                the_data%NdCoords((i1-1)*dim1 + 1) = slam(i1)
                the_data%NdCoords((i1-1)*dim1 + 2) = sfea(i1)
        end do
        do i1 = 1, the_data%NumEl, 1
                the_data%ElConnect((i1-1)*NumND_per_El+1) = nm (i1,1)
                the_data%ElConnect((i1-1)*NumND_per_El+2) = nm (i1,2)
                the_data%ElConnect((i1-1)*NumND_per_El+3) = nm (i1,3)
        end do

        do i1= 1, the_data%NumNd, 1
            the_data%NdOwners(i1) = node_owner(the_data%NdIDs(i1)) - 1
        end do

        j1 = 0
        do i1 = 1, the_data%NumNd, 1
            if (the_data%NdOwners(i1) == localPet) then
                j1 = j1 + 1
                the_data%owned_to_present_nodes(j1) = i1
            end if
        end do
        the_data%ElTypes = ESMF_MESHELEMTYPE_TRI

        !TODO: Saeed: Check if I need to dealocate arrays here!


    end subroutine extract_parallel_data_from_mesh
!
    subroutine extract_parallel_data_from_mesh_orig(global_fort14_dir, the_data,localPet)
        implicit none

        type(meshdata), intent(inout)         :: the_data
        character(len=*), intent(in)          :: global_fort14_dir
        integer, intent(in)                   :: localPet
        character(len=6)                      :: PE_ID, garbage1
        
        character(len=200)                    :: fort14_filename, fort18_filename, partmesh_filename
        integer                               :: i1, j1, i_num, petCount, num_global_nodes, garbage2, garbage3
        integer, allocatable                  :: local_node_numbers(:), local_elem_numbers(:), node_owner(:)
        integer, parameter                    :: dim1=2, NumND_per_El=3

        write(PE_ID, "(A,I4.4)") "PE", localPet
        fort14_filename = TRIM(global_fort14_dir)//'/'//PE_ID//"/fort.14"
        fort18_filename = TRIM(global_fort14_dir)//'/'//PE_ID//"/fort.18"
        partmesh_filename = TRIM(global_fort14_dir)//'/'//'partmesh.txt'


        

        open(unit=23514, file=fort14_filename, form='FORMATTED', status='OLD', action='READ')
        open(unit=23518, file=fort18_filename, form='FORMATTED', status='OLD', action='READ')
        open(unit=235100, file=partmesh_filename, form='FORMATTED', status='OLD', action='READ')

        read(unit=23514, fmt=*)
        read(unit=23514, fmt=*) the_data%NumEl, the_data%NumNd
        allocate(the_data%NdIDs(the_data%NumNd))
        allocate(local_node_numbers(the_data%NumNd))
        allocate(the_data%ElIDs(the_data%NumEl))
        allocate(local_elem_numbers(the_data%NumEl))
        allocate(the_data%NdCoords(dim1*the_data%NumNd))
        allocate(the_data%bathymetry(the_data%NumNd))
        allocate(the_data%ElConnect(NumND_per_El*the_data%NumEl))
        allocate(the_data%NdOwners(the_data%NumNd))
        allocate(the_data%ElTypes(the_data%NumEl))

        read(unit=23518, fmt=*)
        read(unit=23518, fmt=*)
        read(unit=23518, fmt=*) local_elem_numbers
        the_data%ElIDs = abs(local_elem_numbers)
        read(unit=23518, fmt=*) garbage1, num_global_nodes, garbage2, garbage3
        read(unit=23518, fmt=*) local_node_numbers
        the_data%NumOwnedND = 0
        do i1 = 1, the_data%NumNd, 1
            if (local_node_numbers(i1) > 0) then
                the_data%NumOwnedND = the_data%NumOwnedND + 1
            end if
        end do
        the_data%NdIDs = abs(local_node_numbers)
        allocate(node_owner(num_global_nodes))
        allocate(the_data%owned_to_present_nodes(the_data%NumOwnedND))
        read(unit=235100, fmt=*) node_owner

        !print *, 'ADC SLAM > ',size(slam),'>>>>>>',minval(slam),maxval(slam)
        !print *, 'ADC SFEA > ',size(sfea),'>>>>>>',minval(sfea),maxval(sfea)

        do i1 = 1, the_data%NumNd, 1
            read(unit=23514, fmt=*) local_node_numbers(i1), &
                the_data%NdCoords((i1-1)*dim1 + 1), &
                the_data%NdCoords((i1-1)*dim1 + 2), &
                the_data%bathymetry(i1)
        end do
        !i1=100
        !print *,'ADC > lon,lat >>', the_data%NdCoords((i1-1)*dim1 + 1),  the_data%NdCoords((i1-1)*dim1 + 2)
        !i1=400
        !print *,'ADC > lon,lat >>', the_data%NdCoords((i1-1)*dim1 + 1),  the_data%NdCoords((i1-1)*dim1 + 2)
        
        do i1 = 1, the_data%NumEl, 1
            read(unit=23514, fmt=*) local_elem_numbers(i1), i_num, &
                the_data%ElConnect((i1-1)*NumND_per_El+1), &
                the_data%ElConnect((i1-1)*NumND_per_El+2), &
                the_data%ElConnect((i1-1)*NumND_per_El+3)
        end do

        do i1= 1, the_data%NumNd, 1
            the_data%NdOwners(i1) = node_owner(the_data%NdIDs(i1)) - 1
        end do

        j1 = 0
        do i1 = 1, the_data%NumNd, 1
            if (the_data%NdOwners(i1) == localPet) then
                j1 = j1 + 1
                the_data%owned_to_present_nodes(j1) = i1
            end if
        end do
        the_data%ElTypes = ESMF_MESHELEMTYPE_TRI

        close(23514)
        close(23518)
        close(235100)
    end subroutine extract_parallel_data_from_mesh_orig


    subroutine read_config()
    ! This subroutine is not used with NEMS system. Because the 
    ! time interval information is passed via nems.configure file 
    ! with time slot definitation.

    
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

   ! read time coupling interval info
   
    call ESMF_ConfigGetAttribute(cf, adc_cpl_int, label="cpl_int:",default=300, rc=rc)
    call ESMF_ConfigGetAttribute(cf, adc_cpl_num, label="cpl_num:",default=0  , rc=rc)
    call ESMF_ConfigGetAttribute(cf, adc_cpl_den, label="cpl_den:",default=1  , rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
        
    call ESMF_ConfigDestroy(cf, rc=rc) ! Destroy the Config
        
    end subroutine read_config



end module
