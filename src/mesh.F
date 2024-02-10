!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!                      M O D U L E   M E S H
!     ------------------------------------------------------------------
!     jgf51.21.08: Created this mesh module in order to modularize
!     mesh related data. Modularity gives us greater flexibility in
!     reading meshes in different file formats (such as NetCDF or XDMF)
!     or even to read meshes that were originally developed and formatted
!     for other unstructured mesh models (such as DG ADCIRC, RiCOM, FVCOM,
!     SUNTANS, or unstructured SWAN).
!
!     The variables and subroutines in this module were refactored
!     out of the other parts of the code, particularly from the global
!     module.
!     ------------------------------------------------------------------
      module mesh
!     ------------------------------------------------------------------
      use hashtable
      use sizes, only: meshFileName, mne, mnp, mnope, mneta, myproc,
     &   mnbou, mnvel, nprec, meshtype, ascii, xdmf, mnei, GBLINPUTDIR
      use global, only: DEBUG, INFO, ECHO, WARNING, ERROR,
     &   scratchMessage, logMessage, screenMessage, allMessage,
     &   setMessageSource, openFileForRead, unsetMessageSource,
     &   scratchFormat
      use constants, only: deg2rad, rad2deg, pi, omega, rearth

      implicit none
       
      character(len=80) :: agrid
      character(4)              :: aid4(20)
      character(8)              :: aid8(10)
      integer :: np   ! number of nodes in the mesh
      integer :: ne   ! number of elements in the mesh
      integer :: ics  ! mesh coordinate system (1=cartesian, 2=geographic,
                      ! 20-24 geographic w/ corrected equations,
                      ! 20 = equal-area cylindrical,
                      ! 21 = CPP (equi-distant) cylindrical,
                      ! 22 = Mercator (conformal) cylindrical,
                      ! 23 = Miller cylindrical,
                      ! 24 = Gall-Stereographic cylindrical)
      real(8),allocatable, target ::  dp(:)  ! bathymetric depth
      real(8),allocatable         :: dp_orig(:) ! original bathymetric depth
      integer,allocatable :: nm(:,:)          ! element table size(ne,3)
      real(8),allocatable :: x(:)   ! x cartesian node locations x(np)
      real(8),allocatable :: y(:)   ! y cartesian node locations y(np)
      real(8),allocatable :: slam(:) ! longitude node locations slam(np)
      real(8),allocatable :: sfea(:) ! latitude node locations sfea(np)
      type(ipair), allocatable, target :: node_dict(:) ! map node labels to numbers
      integer, allocatable, target :: labels(:) ! map node numbers to labels
      integer :: n1, n2, n3 ! nodes 1, 2, 3 for each element

      integer, allocatable :: lbarray_pointer(:)
      integer, allocatable :: nneigh(:) !(MNP) number of neighboring nodes for each node
      integer, allocatable :: nneighele(:)
      integer, allocatable :: mju(:)
      integer, allocatable :: nodele(:)
      integer, allocatable :: neitab(:,:) !(MNP,MNEI) table of neighbor nodes for each node
      integer, allocatable :: neitabele(:,:) !(MNP,MNei)table of neighboring elements for each node
      real(8), allocatable :: areas(:) !(MNE) = 2*Element Area
      real(8), allocatable :: totalarea(:)
      real(8), allocatable :: bcxy(:,:)  ! element barycenters for kdtree
      real(8), allocatable :: rmax(:)    ! element radii for kdtree
      real(8), allocatable :: csi(:),sii(:)
      real(8) :: slam0,sfea0 !center point of CPP spherical projection
      integer  :: neimin,neimax
 
      real(8) :: theta ! rotation around y-axis
      real(8) :: beta  ! rotation around z-axis
      !c  FDXE(1:3,1:NE) -- b1, b2, b3
      !c  FDYE(1:3,1:NE) -- a1, a2, a3
      REAL(8), ALLOCATABLE:: FDXE(:,:), FDYE(:,:) 
      INTEGER:: NELNDTAB
      INTEGER, ALLOCATABLE:: PERBC_ELNDTAB(:,:) 
      
      REAL(8), POINTER, DIMENSION(:,:):: EXV, EYV ! Nov 22, 2019
      REAL(8), POINTER, DIMENSION(:,:):: EXVO, EYVO
      REAL(8), POINTER, DIMENSION(:,:):: EXVR, EYVR
      REAL(8), ALLOCATABLE, DIMENSION(:):: SFEAR, SLAMR
      !c Matrix mapping two Cartesian frames
      REAL(8):: ROTMAT(3,3) = reshape( (/ 1.0, 0.0, 0.0, 
     &                             0.0, 1.0, 0.0,
     &                             0.0, 0.0, 1.0 /), (/3, 3/)) ;  
      !c 
      !c RVELF: V(\varphi,\theta)  -> V'(\varphi',\theta')
      !c RVELI: V'(\varphi,\theta) -> V(\varphi,\theta)
      REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: RVELF, RVELI   
      
      ! Temporary vectors
      REAL(8), ALLOCATABLE, DIMENSION(:):: UVECTMP, VVECTMP
 
      ! The correction factors
      REAL(8), allocatable, DIMENSION(:) :: SFac, SFMX, SFMY, SFCT,
     &    SFCX, SFCY, YCSFAC, TANPHI ! spherical coordinate tan(lat)/R ratio
!...  Element averages of these...
      REAL(8), allocatable, DIMENSION(:) :: SfacEle, SFCYEle, SFCXEle,
     &    SFMXEle, SFMYEle, YCSFACEle, TANPHIEle 

      ! Overload a functions 
      INTERFACE SPCOORSROTS
        MODULE PROCEDURE SPCOORSROTS1, SPCOORSROTS0
      END INTERFACE

      INTERFACE DRVSPCOORSROTS
        MODULE PROCEDURE DRVSPCOORSROTS0, DRVSPCOORSROTS1
      END INTERFACE

      INTERFACE GET_ROTMAT 
         MODULE PROCEDURE DRV_ROTMAT_ZNVEC, DRV_ROTMAT_ZXZROT
      END INTERFACE

!     ----------
      contains
!     ----------

!     ------------------------------------------------------------------
!             S U B R O U T I N E    R E A D   M E S H
!     ------------------------------------------------------------------
!     jgf51.21.11: Read the ADCIRC mesh in various formats
!     ------------------------------------------------------------------
      subroutine readMesh()
      implicit none

      call setMessageSource("readMesh")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      select case(meshType)
      case(ASCII)
         call read14()
      case(XDMF)
         call readMeshXDMF()
      end select

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

!     ------------------------------------------------------------------
      end subroutine readMesh
!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
!        S U B R O U T I N E    R E A D  1 4   F I N D   D I M S
!     ------------------------------------------------------------------
!     jgf51.21.11: Read the ADCIRC mesh in various formats
!     ------------------------------------------------------------------
      subroutine read14_findDims ()
      use global, only : nabout
      use boundaries
      implicit none
      integer :: ios ! i/o status
      integer :: lineNum ! line number currently being read
      integer :: i, j, k
      integer :: neta_count, nvel_count
      integer :: nvell_max, nvdll_max
      integer, parameter :: iunit = 14

      call setMessageSource("read14FindDims")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      !
      ! initializations
      lineNum = 1
      numSimpleFluxBoundaries = 0
      numExternalFluxBoundaries = 0
      numInternalFluxBoundaries = 0
      numInternalFluxBoundariesWithPipes = 0
      !
      ! reading the file
      read(iunit,'(A80)',err=10,end=20,iostat=ios) agrid
      lineNum = lineNum + 1
      call logMessage(INFO,"Mesh file comment line: "//trim(agrid))
      call logMessage(INFO,"Reading mesh file dimensions.")
      read(unit=iunit,fmt=*,err=10,end=20,iostat=ios) ne, np
      mne = ne
      mnp = np
      do k = 1, np
         read(unit=iunit,fmt=*,err=10,end=20,iostat=ios)
         lineNum = lineNum + 1
      enddo
      do k = 1, ne
         read(unit=iunit,fmt=*,err=10,end=20,iostat=ios)
         lineNum = lineNum + 1
      enddo
      read(unit=iunit,fmt=*,err=10,end=20,iostat=ios) nope  ! total number of elevation boundaries
      mnope = nope
      lineNum = lineNum + 1
      read(unit=iunit,fmt=*,err=10,end=20,iostat=ios) neta  ! total number of nodes on elevation boundaries
      mneta = neta
      lineNum = lineNum + 1
      call allocateElevationBoundaryLengths()
      call allocateAdcircElevationBoundaryArrays()
      neta_count = 0
      nvdll_max = 0
      do k = 1, nope
         read(unit=iunit,fmt=*,err=10,end=20,iostat=ios) nvdll(k) ! number of nodes on the kth elevation boundary segment
         lineNum = lineNum + 1
         nvdll_max = max(nvdll_max,nvdll(k))
         do j = 1, nvdll(k)
            read(unit=iunit,fmt=*,err=10,end=20,iostat=ios)
            lineNum = lineNum + 1
            neta_count = neta_count + 1
         enddo
      enddo
      if ( neta_count.ne.neta ) then
         write(scratchMessage,1020) neta, neta_count
1020    format('Number of open boundary nodes was set to ',
     &           i0,' but ',i0,' were found.')
         call allMessage(WARNING,scratchMessage)
      endif
      read(unit=iunit,fmt=*,err=10,end=20,iostat=ios) nbou ! total number of flux boundaries
      lineNum = lineNum + 1
      read(unit=iunit,fmt=*,err=10,end=20,iostat=ios) nvel ! total number of nodes on flux boundaries
      lineNum = lineNum + 1
      mnbou = nbou
      mnvel = nvel*2
      call allocateFluxBoundaryLengths()
      call allocateAdcircFluxBoundaryArrays()
      nvel_count = 0
      nvell_max = 0
      do k = 1, nbou
         read(unit=iunit,fmt=*,err=10,end=20,iostat=ios)
     &       nvell(k), ibtype_orig(k)  ! number of nodes and type of kth flux boundary
         lineNum = lineNum + 1
         nvell_max = max(nvell_max,nvell(k))
         do j = 1, nvell(k)
            read(unit=iunit,fmt=*,err=10,end=20,iostat=ios)
            lineNum = lineNum + 1
            nvel_count = nvel_count + 1
         enddo
         ! count the total number of each type of boundary for later
         ! use in memory allocation
         select case(ibtype_orig(k))
         case(0,1,2,10,11,12,20,21,22,30,32,52,94,102,112,122,152)
             numSimpleFluxBoundaries = numSimpleFluxBoundaries + 1
         case(3,13,23)
             numExternalFluxBoundaries = numExternalFluxBoundaries + 1
         case(4,24,64)
             numInternalFluxBoundaries = numInternalFluxBoundaries + 1
             ! jgf51.52.24: added nvell(k) for boundaries with backside
             ! nodes, since these are to be included in nvel.
             nvel_count = nvel_count + nvell(k)
         case(5,25)
             numInternalFluxBoundariesWithPipes = numInternalFluxBoundariesWithPipes + 1
             nvel_count = nvel_count + nvell(k)
         case default
             write(scratchMessage, 1030) ibtype_orig(k)
1030        format('The boundary type ',i0,
     &           ' was found in the file but is not valid.')
             call allMessage(ERROR, scratchMessage)
             call terminate()
         end select
      enddo
      if ( nvel_count.ne.nvel) then
         write(scratchMessage, 1040) nvel, nvel_count
1040     format('Number of land boundary nodes was set to ',i0,' but ',
     &      i0,' were found.')
         call allMessage(WARNING, scratchMessage)
         call logMessage(DEBUG,
     &      'Here is the summary of land boundary node information:')
         write(scratchMessage,1050) nvel, nvel_count
1050     format('NVEL (specified number of land boundary nodes) = ',
     &       i0,'. Counted number of land boundary nodes = ',i0,'.')
         call logMessage(DEBUG,scratchMessage)
         if (nabout.eq.DEBUG) then
            do k=1,nbou
               write(scratchMessage,1060) k, ibtype_orig(k), k,
     &            nvell(k), sum(nvell(1:k))

1060          format('ibtype(',i0,')=',i0,', nvell(',i0,')=',i0,
     &                ', total=',i0,'.')
               call logMessage(DEBUG,scratchMessage)
            end do
         endif
      endif
      rewind(iunit)
      call logMessage(INFO,'Finished reading mesh file dimensions.')

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
         !
         ! jump to here on error and end conditions during read
 10   write(scratchMessage,1070) lineNum, ios
1070  format('Reading line ',i0,' gave the following error code: ',
     &   i0,'.')
      call allMessage(ERROR, scratchMessage)
      close(iunit)
      call terminate()
 20   write(scratchMessage,1080) lineNum
1080  format('Reached premature end of file on line ',i0,'.')
      call allMessage(ERROR,scratchMessage)
      close(iunit)
      call terminate()
      !-----+---------+---------+---------+---------+---------+---------+
      end subroutine read14_findDims
      !-----+---------+---------+---------+---------+---------+---------+

      !-----+---------+---------+---------+---------+---------+---------+
      ! READ14
      !-----+---------+---------+---------+---------+---------+---------+
      subroutine read14 ()
      use boundaries
      use global, only : openFileForRead, nabout, nperseg, nnperbc,
     &                   iperconn

      implicit none
      integer :: i, j, k, jn, je, nhy
      integer, parameter :: iunit = 14
      integer :: ios     ! i/o status
      integer :: lineNum ! line number currently being read

      call setMessageSource("read14")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      !
      ! initialization
      nfluxf = 0
      !
      call openFileForRead(iunit,trim(meshFileName),ios)
      call read14_findDims ()
      !
      ! Write a summary of the different boundary types if the logging
      ! level is set to DEBUG.
      if (nabout.eq.DEBUG) then
         !
         ! ocean boundaries
         scratchFormat = '("Number of elevation specified boundaries '
     &    // '(nope): ",i0,".")'
         write(scratchMessage,scratchFormat) nope
         call logMessage(DEBUG,trim(scratchMessage))
         !
         ! land boundaries
         scratchFormat = '("Number of simple flux specified boundaries '
     &     // '(0,1,2,etc): ",i0,".")'
         write(scratchMessage,scratchFormat) numSimpleFluxBoundaries
         call logMessage(DEBUG,trim(scratchMessage))
         !
         ! external over flow (land) boundaries
         scratchFormat = '("Number of external flux boundaries (3,etc):'
     &      // ' ",i0,".")'
         write(scratchMessage,scratchFormat) numExternalFluxBoundaries
         call logMessage(DEBUG,trim(scratchMessage))
         !
         ! levee boundaries
         scratchFormat = '("Number of internal flux boundaries (4,etc):'
     &      // ' ",i0,".")'
         write(scratchMessage,scratchFormat) numInternalFluxBoundaries
         call logMessage(DEBUG,trim(scratchMessage))
         !
         ! levee boundaries with culverts
         scratchFormat = '("Number of internal flux boundaries with '
     &      // 'pipes (5,etc): ",i0,".")'
         write(scratchMessage,scratchFormat) numInternalFluxBoundariesWithPipes
         call logMessage(DEBUG,trim(scratchMessage))
      endif
      !
      call allocateNodalAndElementalArrays()
      ! allocate boundary parameter arrays (barrier height, backface nodes, etc)
      call allocateBoundaryArrays()
      call allocateFluxBoundaryArrayTemporaries()
      !
      !  N O D E   T A B L E
      !

      call logMessage(INFO,'Reading mesh file coordinates, '
     &   // ' connectivity, and boundary data.')
      lineNum = 1
      read(unit=iunit,fmt='(a80)',err=10,end=20,iostat=ios) agrid
      lineNum = lineNum + 1
      read(unit=iunit,fmt=*,err=10,end=20,iostat=ios) ne, np
      lineNum = lineNum + 1
      do k = 1, np
         read(unit=iunit,fmt=*,err=10,end=20,iostat=ios) labels(k), slam(k),
     &     sfea(k), dp(k)
         lineNum = lineNum + 1
         call add_ipair(node_dict, labels(k), k)
      enddo
      !
      !  E L E M E N T   T A B L E
      !
      do k = 1, ne
         read(unit=iunit,fmt=*,err=10,end=20,iostat=ios)
     &      je, nhy, n1, n2, n3
         nm(k, 1) = find(node_dict, n1)
         nm(k, 2) = find(node_dict, n2)
         nm(k, 3) = find(node_dict, n3)
         lineNum = lineNum + 1
      enddo

      !
      !  B O U N D A R I E S
      !
      read(unit=iunit,fmt=*,err=10,end=20,iostat=ios) nope
      lineNum = lineNum + 1
      read(unit=iunit,fmt=*,err=10,end=20,iostat=ios) neta
      lineNum = lineNum + 1
      do k = 1, nope
         read(unit=iunit,fmt=*,err=10,end=20,iostat=ios) nvdll(k)
         lineNum = lineNum + 1
         elevationBoundaries(k)%indexNum = k
         do j = 1, nvdll(k)
            read(unit=iunit,fmt=*,err=10,end=20,iostat=ios) n1
            elevationBoundaries(k)%nodes(j) = find(node_dict, n1)
            elevationBoundaries(k)%xdmf_nodes(j) = elevationBoundaries(k)%nodes(j) - 1
            lineNum = lineNum + 1
         enddo
      enddo
      read(unit=iunit,fmt=*,err=10,end=20,iostat=ios) nbou
      lineNum = lineNum + 1
      read(unit=iunit,fmt=*,err=10,end=20,iostat=ios) nvel
      lineNum = lineNum + 1
      sfCount = 1
      efCount = 1
      ifCount = 1
      ifwpCount = 1
      NPERSEG = 0  ;
      NNPERBC = 0  ;

      do k = 1, nbou
         read(unit=iunit,fmt=*,err=10,end=20,iostat=ios) nvell(k), ibtype_orig(k)
         lineNum = lineNum + 1
         select case(ibtype_orig(k))
         ! flux boundaries (mainland, island, river)
         ! For NON BARRIER type boundaries, read the node numbers for the Kth
         ! boundary segment
         case(0,1,2,10,11,12,20,21,22,30,32,52,102,112,122,152)
            simpleFluxBoundaries(sfCount)%indexNum = k
            do j = 1, nvell(k)
               read(unit=iunit,fmt=*,err=10,end=20,iostat=ios) n1
               simpleFluxBoundaries(sfCount)%nodes(j) = find(node_dict, n1)
               simpleFluxBoundaries(sfCount)%xdmf_nodes(j) =
     &                  simpleFluxBoundaries(sfCount)%nodes(j) - 1
               lineNum = lineNum + 1
            end do
            sfCount = sfCount + 1
         ! For EXTERNAL BARRIER type boundaries, read the node number, the
         ! barrier elevation above the Geoid and the coefficient of free
         ! surface super-critical flow for the nodes in the Kth boundary
         ! segment.
          case(3,13,23)
            externalFluxBoundaries(efCount)%indexNum = k
            do j = 1, nvell(k)
               read(unit=iunit,fmt=*,err=10,end=20,iostat=ios)
     &              n1,
     &              externalFluxBoundaries(efCount)%barlanht(j),
     &              externalFluxBoundaries(efCount)%barlancfsp(j)
               externalFluxBoundaries(efCount)%nodes(j) =
     &              find(node_dict, n1)
                 externalFluxBoundaries(efCount)%xdmf_nodes(j) =
     &              externalFluxBoundaries(efCount)%nodes(j) - 1
               lineNum = lineNum + 1
            end do
            efCount = efCount + 1
         !
         ! For INTERNAL BARRIER type boundaries, read the node number, the
         ! paired node number, the barrier elevation above the Geoid and the
         ! coefficients of supercritical and subcritical flow for the nodes
         ! in the Kth boundary segment.
         case(4,24,64)
            internalFluxBoundaries(ifCount)%indexNum = k
            do j = 1, nvell(k)
               read(unit=iunit,fmt=*,err=10,end=20,iostat=ios)
     &                      n1,
     &                      n2,
     &                      internalFluxBoundaries(ifCount)%barinht(j),
     &                      internalFluxBoundaries(ifCount)%barincfsb(j),
     &                      internalFluxBoundaries(ifCount)%barincfsp(j)
               internalFluxBoundaries(ifCount)%nodes(j) =
     &             find(node_dict, n1)
               internalFluxBoundaries(ifCount)%ibconn(j) = 
     &             find(node_dict, n2)
               internalFluxBoundaries(ifCount)%xdmf_nodes(j) =
     &             internalFluxBoundaries(ifCount)%nodes(j) - 1
               internalFluxBoundaries(ifCount)%xdmf_ibconn(j) =
     &            internalFluxBoundaries(ifCount)%ibconn(j) - 1
               lineNum = lineNum + 1
            end do
            ifCount = ifCount + 1
         !
         ! For INTERNAL BARRIER WITH CROSS BARRIER PIPES type boundaries,
         ! read the node number, the paired node number, the barrier
         ! elevation above the Geoid, the coefficients of supercritical and
         ! subcritical flow, the cross barrier pipe height, the cross barrier
         ! pipe coefficient and the cross barrier pipe diameter for the nodes
         ! in the Kth boundary segment.
         case(5,25)
            internalFluxBoundariesWithPipes(ifwpCount)%indexNum = k
            do j = 1, nvell(k)
               read(unit=iunit,fmt=*,err=10,end=20,iostat=ios)
     &                        n1,
     &                        n2,
     &                        internalFluxBoundariesWithPipes(ifwpCount)%barinht(j),
     &                        internalFluxBoundariesWithPipes(ifwpCount)%barincfsb(j),
     &                        internalFluxBoundariesWithPipes(ifwpCount)%barincfsp(j),
     &                        internalFluxBoundariesWithPipes(ifwpCount)%pipeht(j),
     &                        internalFluxBoundariesWithPipes(ifwpCount)%pipecoef(j),
     &                        internalFluxBoundariesWithPipes(ifwpCount)%pipediam(j)
               internalFluxBoundariesWithPipes(ifwpCount)%nodes(j) =
     &                        find(node_dict, n1)
               internalFluxBoundariesWithPipes(ifwpCount)%ibconn(j) = 
     &                        find(node_dict, n2)          
               internalFluxBoundariesWithPipes(ifwpCount)%xdmf_nodes(j) =
     &                        internalFluxBoundariesWithPipes(ifwpCount)%nodes(j) - 1
               internalFluxBoundariesWithPipes(ifwpCount)%xdmf_ibconn(j) =
     &                        internalFluxBoundariesWithPipes(ifwpCount)%ibconn(j) - 1
               lineNum = lineNum + 1
            end do
            ifwpCount = ifwpCount + 1
         case(94)
!...        Periodic boundary, read in a pair of nodes 
            !c Allocate memory, if there is a periodic bc segment c!
            IF ( .NOT. ALLOCATED(IPERCONN) ) THEN
               ALLOCATE( IPERCONN(nvell(k),2) ) ;
               IPERCONN = 0 ;
            END IF
            !
            NPERSEG = NPERSEG + 1 ;
            simpleFluxBoundaries(sfCount)%indexNum = k
            do j = 1, nvell(k)
               NNPERBC =  NNPERBC + 1 ;
               read(unit=iunit,fmt=*,err=10,end=20,iostat=ios) n1, n2
               IPERCONN(NNPERBC,1) = find(node_dict, n1) 
               IPERCONN(NNPERBC,2) = find(node_dict, n2) 
               simpleFluxBoundaries(sfCount)%nodes(j) = 
     &                  IPERCONN(NNPERBC,1)
               simpleFluxBoundaries(sfCount)%xdmf_nodes(j) =
     &                  simpleFluxBoundaries(sfCount)%nodes(j) - 1
               lineNum = lineNum + 1
            end do
            sfCount = sfCount + 1
         case default
            scratchFormat = '("IBTYPE ",i0," is not allowed.")'
            write(scratchMessage, scratchFormat) ibtype_orig(k)
            call allMessage(ERROR,scratchFormat)
            call terminate()
         end select
      end do
      close(14)
      !
      ! populate the adcirc arrays that are used during execution
      call populateADCIRCNativeArrays()
      !
      call logMessage(INFO,'Finished reading mesh file coordinates, '
     & // 'connectivity, and boundary data.')

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
       !
      ! jump to here on error condition during read
10    write(scratchMessage,140) lineNum, ios
140   format('Reading line ',i0,' gave the following error code: ',i0,'.')
      call allMessage(ERROR,scratchMessage)
      close(iunit)
      call terminate()
      !
      ! jump to here on end condition during read
20    write(scratchMessage,150) lineNum
150   format('Reached premature end of file on line ',i0,'.')
      call allMessage(ERROR,scratchMessage)
      close(iunit)
      call close_dict( node_dict )
      call terminate()

      !----------------------------------------------------------------
      end subroutine read14
      !----------------------------------------------------------------



      !-----+---------+---------+---------+---------+---------+---------+
      ! READ14MESH() : Read just mesh
      !  input:
      !      meshfileName = file name
      !  output:
      !      nn = (/ np, ne /)
      !      vx(2,np)   - coordinates of FE nodes
      !      bxy(np)    - bathymetry
      !      etov(3,ne) - Element conectivity table
      !      (optional) nodeLael - node label
      !-----+---------+---------+---------+---------+---------+---------+
      subroutine read14femesh( meshfileName, nn, vx, etov, bxy, nodeLabel )
      use global, only : openFileForRead, nabout

      implicit none

      ! dummy
      CHARACTER (LEN=*):: meshfileName
      INTEGER:: nn(:)

      REAL (8), allocatable:: vx(:,:), bxy(:) ! node, bathymetry
      INTEGER, allocatable:: etov(:,:) ! element connectity table
      INTEGER, allocatable, optional:: nodeLabel(:)

      ! local
      integer :: i, j, k, jn, je, nhy
      integer, parameter :: iunit = 714
      integer :: ios     ! i/o status
      integer :: lineNum ! line number currently being read
      CHARACTER (LEN=80):: agridtmp

      INTEGER:: ne, np, nfluxftmp
      INTEGER, allocatable:: labelstmp(:)

      ! hash table
      type(ipair), allocatable, target :: node_dict_tmp(:) ! map node labels to numbers

      call setMessageSource( trim(meshfileName) )
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      !
      ! initialization
      nfluxftmp = 0
      !
      call openFileForRead(iunit,trim(meshFileName),ios)

      ! Get dimension
      !  - reading the file
      read(iunit,'(A80)',err=110,end=120,iostat=ios) agridtmp ! header
      lineNum = lineNum + 1

      read(unit=iunit,fmt=*,err=110,end=120,iostat=ios) ne, np ! no. elements, nodes
      lineNum = lineNum + 1

      nn(1) = np
      nn(2) = ne

      ! Allocate memory for the
      IF ( allocated( vx ) ) THEN
         DEALLOCATE( vx )
      END IF
      IF ( allocated( etov ) ) THEN
         DEALLOCATE( etov )
      END IF
      IF ( allocated( bxy ) ) THEN
         DEALLOCATE( bxy )
      END IF
      ALLOCATE( vx( 2, np ), bxy(np), etov(3, ne) )

      IF ( present( nodeLabel ) ) THEN
         IF ( allocated( nodeLabel ) ) THEN
            DEALLOCATE( nodeLabel )
         END IF
         ALLOCATE( nodeLabel( np ) )
      END IF
      ! Intialize a hash table
      call dict( node_dict_tmp, np )

      !  N O D E   T A B L E
      !
      !
      ALLOCATE( labelstmp(np) )
      do k = 1, np
         read(unit=iunit,fmt=*,err=110,end=120,iostat=ios) labelstmp(k), vx(1,k),
     &     vx(2,k), bxy(k)
         lineNum = lineNum + 1
         call add_ipair( node_dict_tmp, labelstmp(k), k)
      enddo
      IF ( present( nodeLabel ) ) THEN
         DO k = 1, np
            nodeLabel(k) = labelstmp(k)
         END DO
      END IF

      !
      !  E L E M E N T   T A B L E
      !
      do k = 1, ne
         read(unit=iunit,fmt=*,err=110,end=120,iostat=ios)
     &      je, nhy, n1, n2, n3
         etov(1,k) = find(node_dict_tmp, n1)
         etov(2,k) = find(node_dict_tmp, n2)
         etov(3,k) = find(node_dict_tmp, n3)
         lineNum = lineNum + 1
      enddo

      ! close file
      close(iunit)
      call close_dict( node_dict_tmp ) !

      ! free memory
      deallocate( labelstmp )

!      !
!      ! populate the adcirc arrays that are used during execution
!     call populateADCIRCNativeArrays()
!      !
!      call logMessage(INFO,'Finished reading mesh file coordinates, '
!     & // 'connectivity, and boundary data.')

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage( DEBUG, "Return.")
#endif
      call unsetMessageSource()

      return
      !
      ! jump to here on error condition during read
110    write(scratchMessage,140) lineNum, ios
140   format('Reading line ',i0,' gave the following error code: ', i0, '.')
      call allMessage(ERROR,scratchMessage)
      close(iunit)
      call terminate()
      !
      ! jump to here on end condition during read
120    write(scratchMessage,150) lineNum
150   format('Reached premature end of file on line ',i0,'.')
      call allMessage(ERROR,scratchMessage)
      close(iunit)

      ! free memory
      deallocate( labelstmp )

      call close_dict( node_dict_tmp )
      call terminate()

      !----------------------------------------------------------------
      end subroutine read14femesh
      !----------------------------------------------------------------


!----------------------------------------------------------------------
!         S U B R O U T I N E   R E A D   M E S H   X D M F
!----------------------------------------------------------------------
!     Reads the mesh node table, element table, and boundaries from an
!     XDMF file, allocating memory along the way and populating adcirc-style
!     data structures.
!----------------------------------------------------------------------
      subroutine readMeshXDMF()
      use global, only : nabout
      use boundaries
      implicit none

#ifndef ADCXDMF
      call allMessage(ERROR,'An XMDF mesh file was specified.')
      call allMessage(ERROR,'This ADCIRC executable was not compiled '
     & // 'with XDMF support.')
      call terminate()
#else
      include 'adcirc_Xdmf.f'
      integer*8 :: xdmfFortranObj
      integer, allocatable :: xdmf_nm(:,:)   ! 0-offset connectivity array
      integer, allocatable :: setSize(:)
      integer, parameter :: keyLength = 256
      integer, parameter :: valueLength = 256
      integer, parameter :: tagLength = 256
      integer, parameter :: nameLength = 256
      character(len=keyLength) :: itemKey
      character(len=valueLength) :: itemValue
      character(len=tagLength) :: itemTag
      character(len=nameLength) :: itemName
      integer, parameter :: gridCollectionIndex = 0
      integer :: gridIndex = 0
      integer :: typeHolder
      integer :: topologyPropertyIndex
      integer :: geometryIndex
      integer :: startIndex
      integer :: arrayStride
      integer :: valueStride
      integer :: setIndex
      integer :: setPropertyIndex
      integer :: mySetSize
      integer :: attributeIndex
      integer :: attributeDataType
      integer :: attributePropertyIndex
      integer :: numAttributeProperties
      integer :: numAttributes
      integer :: numTopologyProperties
      integer :: numSetProperties
      integer :: numSets
      integer :: topologyDataType
      integer :: setDataType
      integer :: numElementValues
      integer :: setType
      integer :: geometryType
      integer :: geometryDataType
      integer :: topologyType
      integer :: fluxCount
      integer :: elevCount
      integer :: infoIndex
      integer :: propertyIndex
      integer, allocatable :: boundaryTypes(:)
      logical, allocatable :: elevationBoundary(:) ! .true. if the boundary is elevation-specified
      integer, allocatable :: setData(:)
      integer, allocatable :: firstSetAttributeIndex(:) ! the index of the first attribute that corresponds to each set
      character(len=256) :: setDataTypeString
      character(len=256) :: geometryTypeString
      character(len=256) :: topologyDataTypeString
      character(len=256) :: topologyTypeString
      character(len=256) :: geometryDataTypeString
      character(len=256) :: setTypeString
      character(len=256) :: gridName
      integer :: numContained
      integer :: numInformations
      integer :: openAttributes
      integer :: openInformations
      integer :: openMaps
      integer :: openSets
      integer :: attStart
      integer :: numGridCollections
      real(8) :: timeSec
      logical :: fileFound
      integer i, j, k
      real(8), allocatable :: tempCoord(:,:)
      real(8),allocatable :: dmyrealarray(:)
      !

      call setMessageSource("readMeshXDMF")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      startIndex = 0
      arrayStride = 1
      valueStride = 1
      openMaps = 1
      openAttributes = 1
      openInformations = 1
      openSets = 1
      fileFound = .false.
      !
      inquire(file=trim(adjustl(meshFileName)),exist=fileFound)
      if (fileFound.eqv..false.) then
         call allMessage(ERROR,"The file " // trim(adjustl(meshFileName)) // " was not found.")
         call terminate()
      endif
      call logMessage(INFO,'Reading mesh data from the ' //
     &  trim(adjustl(meshFileName)) // ' XDMF file.')
      call xdmfInit(xdmfFortranObj)
      call xdmfRead(xdmfFortranObj, trim(adjustl(meshFileName))//char(0))
      !
      !  A G R I D
      !
      call xdmfOpenDomainGrid(xdmfFortranObj, XDMF_GRID_TYPE_UNSTRUCTURED,
     &   gridIndex, openMaps, openAttributes, openInformations, openSets)
      call xdmfRetrieveDomainGridName(xdmfFortranObj, XDMF_GRID_TYPE_UNSTRUCTURED,
     &    gridIndex, gridName, nameLength)
      call replaceNullsWithSpaces(gridName)
      agrid(1:80) = gridName(1:80)
      !
      !  S I Z E S
      !
      call xdmfRetrieveGeometryNumPoints(xdmfFortranObj, np)
      scratchFormat = '("The Geometry contains ",i0," points (nodes).")'
      write(scratchMessage,scratchFormat) np
      call logMessage(INFO,trim(scratchMessage))
      call xdmfRetrieveTopologyNumElements(xdmfFortranObj, ne)
      scratchFormat = '("The Topology contains ",i0," elements.")'
      write(scratchMessage,scratchFormat) ne
      call logMessage(INFO,trim(scratchMessage))
      !
      mnp = np
      mne = ne
      call allocateNodalAndElementalArrays()
      !
      !  N O D E   T A B L E
      !
      call xdmfRetrieveGeometryType(xdmfFortranObj, geometryType)
      select case(geometryType)
      case(301)
         geometryTypeString = 'XDMF_GEOMETRY_TYPE_XYZ'
      case(302)
         geometryTypeString = 'XDMF_GEOMETRY_TYPE_XY'
      case default
         scratchFormat = '("Unrecognized geometry type ",i0,".")'
         write(scratchMessage,scratchFormat) geometryType
         call allMessage(WARNING,trim(scratchMessage))
      end select
      !
      call xdmfRetrieveGeometryValueType(xdmfFortranObj, geometryDataType)
      call createDataTypeString(geometryDataType, geometryDataTypeString)
      !
      call xdmfRetrieveGeometrySize(xdmfFortranObj, numContained)
      !
      allocate(tempCoord(2,np))
      tempCoord = -99999.d0
      call logMessage(INFO,'Reading nodal table.')
      call xdmfRetrieveGeometryValues(xdmfFortranObj, tempCoord,
     &    geometryDataType, 2*np, startIndex, arrayStride, valueStride)
      !
      slam = tempCoord(1,:)
      sfea = tempCoord(2,:)
      deallocate(tempCoord)
      !
      ! depth
      call xdmfRetrieveNumAttributes(xdmfFortranObj, numAttributes)
      scratchFormat = '("Grid ",i0," contains ",i0," attributes.")'
      write(scratchMessage,scratchFormat) gridIndex, numAttributes
      call logMessage(INFO,trim(scratchMessage))
      do attributeIndex=0, numAttributes - 1
         call xdmfRetrieveAttributeName(xdmfFortranObj, attributeIndex,
     &       itemName, nameLength)
         call replaceNullsWithSpaces(itemName)
         scratchFormat = '("Grid ",i0," Attribute ",i0," is named ",a)'
         write(scratchMessage,scratchFormat) gridIndex, attributeIndex, trim(itemName)
         call logMessage(INFO,trim(scratchMessage))
         select case(trim(itemName))
         case("depth")
            call xdmfRetrieveAttributeValues(xdmfFortranObj, attributeIndex,
     &         dp, XDMF_ARRAY_TYPE_FLOAT64, np, startIndex, arrayStride, valueStride)
         case default
            ! this is not the attribute we're looking for, at least not yet
         end select
      end do
      !
      if (nabout.eq.DEBUG) then
         scratchFormat =
     &    '("node=",i0," x=",f15.7," y=",f15.7," depth=",f8.3)'
         do i=1,np
            write(scratchMessage,scratchFormat) i, slam(i), sfea(i), dp(i)
            call logMessage(ECHO,trim(scratchMessage))
         end do
      end if
      !
      !  E L E M E N T   T A B L E
      !
      call logMessage(INFO,'Reading element table.')
      call xdmfRetrieveTopologyType(xdmfFortranObj, topologyType)
      ! XDMF supports other topology types, but these are unlikely to
      ! be found in ADCIRC output files
      call createTopologyTypeString(topologyType, topologyTypeString)
      scratchFormat = '("The topology type is ",a,".")'
      write(scratchMessage,scratchFormat) trim(topologyTypeString)
      call logMessage(DEBUG,scratchMessage)
      !
      call xdmfRetrieveTopologyValueType(xdmfFortranObj, topologyDataType)
      call createDataTypeString(topologyDataType, topologyDataTypeString)
      scratchFormat = '("The data type of the topology is ",a,".")'
      write(scratchMessage,scratchFormat) trim(topologyDataTypeString)
      call logMessage(DEBUG,scratchMessage)
      !
      call xdmfRetrieveTopologySize(xdmfFortranObj, numElementValues)
      allocate(xdmf_nm(3,ne))
      call xdmfRetrieveTopologyValues(xdmfFortranObj, xdmf_nm,
     &    XDMF_ARRAY_TYPE_INT32, numElementValues, startIndex,
     &    arrayStride, valueStride)
      !
      ! need to add 1 since XDMF stores arrays as 0 offset but ADCIRC reads
      ! them as 1 offset
      !
      do j=1,ne
         do k=1,3
            nm(j,k) = xdmf_nm(k,j) + 1
         end do
      end do
      deallocate(xdmf_nm)
      !
      if (nabout.eq.DEBUG) then
         scratchFormat = '("Element ",i0," nodes ",i0," ",i0," ",i0)'
         do i=1,ne
            write(scratchMessage, scratchFormat) i, (nm(i,j), j=1,3)
            call logMessage(ECHO,scratchMessage)
         end do
      endif
      !
      ! B O U N D A R I E S
      !
      ! boundaries are saved as sets and attributes
      call logMessage(INFO,'Reading boundaries.')
      CALL xdmfRetrieveNumSets(xdmfFortranObj, numSets)
      scratchFormat = '("The Grid contains ",i0," boundaries.")'
      write(scratchMessage,scratchFormat) numSets
      call logMessage(INFO,scratchMessage)
      !
      ! allocate xdmf-specific arrays to hold boundary parameters
      allocate(elevationBoundary(0:numSets-1))
      allocate(setSize(0:numSets-1))
      allocate(firstSetAttributeIndex(0:numSets-1))
      allocate(boundaryTypes(0:numSets-1))
      !
      ! count the different types of boundaries for use in memory allocation
      !
      call logMessage(DEBUG,'Counting the various boundary types.')
      numSimpleFluxBoundaries = 0
      numExternalFluxBoundaries = 0
      numInternalFluxBoundaries = 0
      numInternalFluxBoundariesWithPipes = 0
      nope = 0
      neta = 0
      nbou = 0
      nvel = 0
      elevationBoundary(:) = .false.
      do setIndex=0, numSets-1
         call xdmfRetrieveNumAttributes(xdmfFortranObj, numAttributes)
         firstSetAttributeIndex(setIndex) = numAttributes
         scratchFormat = '("Opening set ",i0,".")'
         write(scratchMessage,scratchFormat) setIndex
         call logMessage(DEBUG,scratchMessage)
         call xdmfOpenSet(xdmfFortranObj, setIndex, openAttributes,
     &         openInformations)
         !
         ! get the boundary type
         call xdmfRetrieveNumInformation(xdmfFortranObj, numInformations)
         call xdmfRetrieveInformation(xdmfFortranObj, numInformations-1,
     &       itemKey, keyLength, itemValue, valueLength)
         call replaceNullsWithSpaces(itemValue)
         read(itemValue,*) boundaryTypes(setIndex) ! value of either ibtype_orig or ibtypee
         !
         ! determine the number of nodes on this boundary
         call xdmfRetrieveSetSize(xdmfFortranObj, setSize(setIndex), setIndex)
         scratchFormat = '("Set ",i0," contains ",i0," values.")'
         write(scratchMessage,scratchFormat) setIndex, setSize(setIndex)
         call logMessage(DEBUG,scratchMessage)
         !
         call xdmfRetrieveSetNumProperties(xdmfFortranObj, setIndex, numSetProperties)
         do setPropertyIndex=0, numSetProperties-1
            call xdmfRetrieveSetProperty(xdmfFortranObj, setIndex,
     &           setPropertyIndex, itemKey, keyLength, itemValue, valueLength)
            call replaceNullsWithSpaces(itemValue)
            select case(trim(itemValue))
            case("elevation_specified_boundary ")
               scratchFormat = '("The set property is ",a,a,a,".")'
               write(scratchMessage,scratchFormat)'"',trim(itemValue),'"'
               call logMessage(ECHO,scratchMessage)
               call logMessage(DEBUG,
     &              'Found one elevation specified boundary.')
               nope = nope + 1
               neta = neta + setSize(setIndex)
               elevationBoundary(setIndex) = .true.
            case("flux_specified_boundary")
               scratchFormat = '("The set property is ",a,a,a,".")'
               write(scratchMessage,scratchFormat)'"',trim(itemValue),'"'
               call logMessage(ECHO,scratchMessage)
               call logMessage(DEBUG,
     &              'Found one flux specified boundary.')
               nbou = nbou + 1
               nvel = nvel + setSize(setIndex)
               select case(boundaryTypes(setIndex))
               case(0,1,2,10,11,12,20,21,22,30,32,52)
                  numSimpleFluxBoundaries = numSimpleFluxBoundaries + 1
               case(3,13,23)
                  numExternalFluxBoundaries = numExternalFluxBoundaries + 1
               case(4,24,64)
                  numInternalFluxBoundaries = numInternalFluxBoundaries + 1
               case(5,25)
                  numInternalFluxBoundariesWithPipes = numInternalFluxBoundariesWithPipes + 1
               case default
                  scratchFormat = 'Mesh file contains IBTYPE=",i0," '
     &               // 'which is not a valid flux boundary type.")'
                  write(scratchMessage,scratchFormat) boundaryTypes(setIndex)
                  call allMessage(ERROR,scratchMessage)
                  call terminate()
               end select
            case("Node")
               ! do nothing, this property simply indicates that the boundaries
               ! are defined by lists of nodes
            case default
               scratchFormat = '("Unrecognized set property ",a,".")'
               write(scratchMessage,scratchFormat) trim(itemValue)
               call allMessage(WARNING,scratchMessage)
            end select
         end do
      end do
      !
      ! Now that we know how many of each boundary type we have, we can
      ! allocate memory to hold the data/parameters for each boundary of
      ! each type
      scratchFormat = '("Number of elevation boundaries : ",i0,".")'
      write(scratchMessage,scratchFormat) nope
      call logMessage(INFO,scratchMessage)
      scratchFormat = '("Number of elevation boundary nodes : ",i0,".")'
      write(scratchMessage,scratchFormat) neta
      call logMessage(INFO,scratchMessage)
      scratchFormat = '("Total number of flux boundaries : ",i0,".")'
      write(scratchMessage,scratchFormat) nbou
      call logMessage(INFO,scratchMessage)
      scratchFormat = '("Total number of flux boundary nodes : ",i0,".")'
      write(scratchMessage,scratchFormat) nvel
      call logMessage(INFO,scratchMessage)

      scratchFormat = '("Number of simple flux boundaries : ",i0,".")'
      write(scratchMessage,scratchFormat) numSimpleFluxBoundaries
      call logMessage(INFO,scratchMessage)
      scratchFormat = '("Number of external flux boundaries : ",i0,".")'
      write(scratchMessage,scratchFormat) numExternalFluxBoundaries
      call logMessage(INFO,scratchMessage)
      scratchFormat = '("Number of internal flux boundaries : ",i0,".")'
      write(scratchMessage,scratchFormat) numInternalFluxBoundaries
      call logMessage(INFO,scratchMessage)
      scratchFormat = '("Num. int. flux bnd. with cross barrier pipes '
     &                 // ' : ",i0,".")'
      write(scratchMessage,scratchFormat) numInternalFluxBoundariesWithPipes
      call logMessage(INFO,scratchMessage)
      !
      ! populate nvdll (adcirc array representing number of nodes on each
      ! elevation boundary segment) and nvell (adcirc array representing
      ! number of nodes on each flux boundary segment) as these variables
      ! are used in the allocateBoundaryArrays() subroutine
      mnope = nope
      mneta = neta
      mnbou = nbou
      mnvel = nvel*2
      call allocateElevationBoundaryLengths()
      call allocateFluxBoundaryLengths()
      elevCount = 1
      fluxCount = 1
      do setIndex=0, numSets-1
         if (elevationBoundary(setIndex).eqv..true.) then
            nvdll(elevCount) = setSize(setIndex)
            ibtypee(elevCount) = boundaryTypes(setIndex)
            elevCount = elevCount + 1
         else
            nvell(fluxCount) = setSize(setIndex)
            ibtype_orig(fluxCount) = boundaryTypes(setIndex) ! preserve the original
            ibtype(fluxCount) = ibtype_orig(fluxCount)       ! may be changed in adcirc (rivers)
            fluxCount = fluxCount + 1
         endif
      end do
      !
      ! allocate boundary parameter arrays (barrier height, backface nodes, etc)
      call allocateBoundaryArrays()
      call allocateFluxBoundaryArrayTemporaries()
      !
      ! allocate boundary-related variables that are used by adcirc internally
      call allocateAdcircElevationBoundaryArrays()
      call allocateAdcircFluxBoundaryArrays()
      !
      ! iterate over all boundaries, read data relevant to each boundary,
      ! and populate data structures
      elevCount = 1
      fluxCount = 1
      sfCount = 1
      efCount = 1
      ifCount = 1
      ifwpCount = 1
      do setIndex=0,numSets-1
         attStart = firstSetAttributeIndex(setIndex)
         call xdmfRetrieveSetValueType(xdmfFortranObj, setIndex, setDataType)
         call createDataTypeString(setDataType, setDataTypeString)
         !
         ! elevation boundary
         if (elevationBoundary(setIndex).eqv..true.) then
            ! get the node numbers on the boundary
            elevationBoundaries(elevCount)%indexNum = elevCount
            call xdmfRetrieveSetValues(xdmfFortranObj, setIndex,
     &          elevationBoundaries(elevCount)%xdmf_nodes, setDataType,
     &          setSize(setIndex), startIndex, arrayStride, valueStride)
            ! convert node numbers from 0 starting index (xdmf-style) to
            ! 1 starting index (fortran-style)
            elevationBoundaries(elevCount)%nodes =
     &          elevationBoundaries(elevCount)%xdmf_nodes + 1
            elevCount = elevCount + 1
         else
            !
            ! flux boundary type
            select case(ibtype_orig(fluxCount))
            case(0,1,2,10,11,12,20,21,22,30,32,52,102,112,122,152)
               simpleFluxBoundaries(sfCount)%indexNum = fluxCount
               ! get node numbers on the boundary
               call xdmfRetrieveSetValues(xdmfFortranObj, setIndex,
     &             simpleFluxBoundaries(sfCount)%xdmf_nodes, setDataType,
     &             setSize(setIndex), startIndex, arrayStride, valueStride)
               ! convert node numbers from 0 starting index (xdmf-style) to
               ! 1 starting index (fortran-style)
               simpleFluxBoundaries(sfCount)%nodes =
     &              simpleFluxBoundaries(sfCount)%xdmf_nodes + 1
               sfCount = sfCount + 1
            case(3,13,23)
               externalFluxBoundaries(efCount)%indexNum = fluxCount
               ! get the node numbers on the boundary
               call xdmfRetrieveSetValues(xdmfFortranObj, setIndex,
     &             externalFluxBoundaries(efCount)%xdmf_nodes, setDataType,
     &             setSize(setIndex),startIndex, arrayStride, valueStride)
               ! convert node numbers from 0 starting index (xdmf-style) to
               ! 1 starting index (fortran-style)
               externalFluxBoundaries(efCount)%nodes =
     &             externalFluxBoundaries(efCount)%xdmf_nodes + 1
               do i=attStart,attStart+1
                  call xdmfRetrieveAttributeName(xdmfFortranObj, i,
     &                 itemName, nameLength)
                  call replaceNullsWithSpaces(itemName)
                  select case(trim(itemName))
                  case("BARLANHT")   ! barrier height at each node
                     call xdmfRetrieveAttributeValues(xdmfFortranObj, i,
     &                   externalFluxBoundaries(efCount)%barlanht, XDMF_ARRAY_TYPE_FLOAT64,
     &                   nvell(fluxCount), startIndex, arrayStride, valueStride)
                  case("BARLANCFSP") ! coefficient of free surface super critical flow
                     call xdmfRetrieveAttributeValues(xdmfFortranObj, i,
     &                   externalFluxBoundaries(efCount)%barlancfsp, XDMF_ARRAY_TYPE_FLOAT64,
     &                   nvell(fluxCount), startIndex, arrayStride, valueStride)
                  case default
                     scratchFormat =
     &                 '("Unrecognized boundary attribute : ",a,".")'
                     write(scratchMessage,scratchFormat) trim(itemName)
                     call allMessage(ERROR,scratchMessage)
                  end select
               end do
               efCount = efCount + 1
            case(4,24,64)  ! internal barrier boundary (e.g., subgrid scale levee)
               internalFluxBoundaries(ifCount)%indexNum = fluxCount
               ! get the node numbers on the boundary
               call xdmfRetrieveSetValues(xdmfFortranObj, setIndex,
     &             internalFluxBoundaries(ifCount)%xdmf_nodes, setDataType, setSize(setIndex),
     &             startIndex, arrayStride, valueStride)
               ! convert node numbers from 0 starting index (xdmf-style) to
               ! 1 starting index (fortran-style)
               internalFluxBoundaries(ifCount)%nodes =
     &              internalFluxBoundaries(ifCount)%xdmf_nodes + 1
               do i=attStart,attStart+3
                  call xdmfRetrieveAttributeName(xdmfFortranObj, i, itemName, nameLength)
                  call replaceNullsWithSpaces(itemName)
                  select case(trim(itemName))
                  case("IBCONN")      ! paired (back face) nodes
                     allocate(dmyrealarray(1:size(internalFluxBoundaries(ifCount)%xdmf_ibconn)))
                     call xdmfRetrieveAttributeValues(xdmfFortranObj, i,
     &                   dmyrealarray, XDMF_ARRAY_TYPE_INT32,
     &                   nvell(fluxCount), startIndex, arrayStride, valueStride)
                     internalFluxBoundaries(ifCount)%xdmf_ibconn = int(dmyrealarray)
                     deallocate(dmyrealarray)
                  case("BARINHT")     ! barrier height at each node and its paired node
                     call xdmfRetrieveAttributeValues(xdmfFortranObj, i,
     &                  internalFluxBoundaries(ifCount)%barinht, XDMF_ARRAY_TYPE_FLOAT64,
     &                   nvell(fluxCount), startIndex, arrayStride, valueStride)
                  case("BARINCFSB")   ! coefficient of free surface sub critical flow
                     call xdmfRetrieveAttributeValues(xdmfFortranObj, i,
     &                   internalFluxBoundaries(ifCount)%barincfsb, XDMF_ARRAY_TYPE_FLOAT64,
     &                   nvell(fluxCount), startIndex, arrayStride, valueStride)
                  case("BARINCFSP")   ! coefficient of free surface super critical flow
                     call xdmfRetrieveAttributeValues(xdmfFortranObj, i,
     &                   internalFluxBoundaries(ifCount)%barincfsp, XDMF_ARRAY_TYPE_FLOAT64,
     &                   nvell(fluxCount), startIndex, arrayStride, valueStride)
                  case default
                     scratchFormat =
     &                '("Unrecognized boundary attribute : ",a,".")'
                     write(scratchMessage,scratchFormat) trim(itemName)
                     call allMessage(ERROR,scratchMessage)
                  end select
               end do
               internalFluxBoundaries(ifCount)%ibconn =
     &             internalFluxBoundaries(ifCount)%xdmf_ibconn + 1
               ifCount = ifCount + 1
            case(5,25)  ! internal barrier boundary with cross barrier pipes
               internalFluxBoundariesWithPipes(ifwpCount)%indexNum = fluxCount
               ! get the node numbers on the boundary
               call xdmfRetrieveSetValues(xdmfFortranObj, setIndex,
     &             internalFluxBoundariesWithPipes(ifwpCount)%xdmf_nodes,
     &             setDataType, setSize(setIndex),
     &             startIndex, arrayStride, valueStride)
               ! convert node numbers from 0 starting index (xdmf-style) to
               ! 1 starting index (fortran-style)
               internalFluxBoundariesWithPipes(ifwpCount)%nodes =
     &             internalFluxBoundaries(ifwpCount)%xdmf_nodes + 1
               do i=attStart,attStart+6
                  call xdmfRetrieveAttributeName(xdmfFortranObj, i, itemName, nameLength)
                  call replaceNullsWithSpaces(itemName)
                  select case(trim(itemName))
                  case("IBCONN")  ! paired (i.e., back face) nodes
                     allocate(dmyrealarray(1:size(internalFluxBoundaries(ifwpCount)%xdmf_ibconn)))
                     call xdmfRetrieveAttributeValues(xdmfFortranObj, i,
     &                  dmyrealarray,
     &                  XDMF_ARRAY_TYPE_INT32,
     &                   nvell(fluxCount), startIndex, arrayStride, valueStride)
                     internalFluxBoundaries(ifwpCount)%xdmf_ibconn = int(dmyrealarray)
                     deallocate(dmyrealarray)
                  case("BARINHT") ! barrier height at each node and its paired node
                     call xdmfRetrieveAttributeValues(xdmfFortranObj, i,
     &                   internalFluxBoundariesWithPipes(ifwpCount)%barinht,
     &                  XDMF_ARRAY_TYPE_FLOAT64,
     &                   nvell(fluxCount), startIndex, arrayStride, valueStride)
                  case("BARINCFSB") ! coefficient of free surface sub critical flow
                     call xdmfRetrieveAttributeValues(xdmfFortranObj, i,
     &                   internalFluxBoundariesWithPipes(ifwpCount)%barincfsb,
     &                   XDMF_ARRAY_TYPE_FLOAT64,
     &                   nvell(fluxCount), startIndex, arrayStride, valueStride)
                  case("BARINCFSP") ! coefficient of free surface super critical flow
                     call xdmfRetrieveAttributeValues(xdmfFortranObj, i,
     &                   internalFluxBoundariesWithPipes(ifwpCount)%barincfsp,
     &                   XDMF_ARRAY_TYPE_FLOAT64,
     &                   nvell(fluxCount), startIndex, arrayStride, valueStride)
                  case("PIPEHT")    ! barrier height at each node and its paired node
                     call xdmfRetrieveAttributeValues(xdmfFortranObj, i,
     &                   internalFluxBoundariesWithPipes(ifwpCount)%pipeht,
     &                    XDMF_ARRAY_TYPE_FLOAT64,
     &                    nvell(fluxCount), startIndex, arrayStride, valueStride)
                  case("PIPECOEF")   ! coefficient of free surface sub critical flow
                     call xdmfRetrieveAttributeValues(xdmfFortranObj, i,
     &                   internalFluxBoundariesWithPipes(ifwpCount)%pipecoef,
     &                   XDMF_ARRAY_TYPE_FLOAT64,
     &                   nvell(fluxCount), startIndex, arrayStride, valueStride)
                  case("PIPEDIAM")   ! coefficient of free surface super critical flow
                     call xdmfRetrieveAttributeValues(xdmfFortranObj, i,
     &                   internalFluxBoundariesWithPipes(ifwpCount)%pipediam,
     &                   XDMF_ARRAY_TYPE_FLOAT64,
     &                   nvell(fluxCount), startIndex, arrayStride, valueStride)
                  case default
                     scratchFormat =
     &                '("Unrecognized boundary attribute : ",a,".")'
                     write(scratchMessage,scratchFormat) trim(itemName)
                     call allMessage(ERROR,scratchMessage)
                  end select
               end do
               internalFluxBoundariesWithPipes(ifwpCount)%ibconn =
     &             internalFluxBoundariesWithPipes(ifwpCount)%xdmf_ibconn + 1
               ifwpCount = ifwpCount + 1
            case default
                scratchFormat = '("The boundary type ",i0," '
     &                 // 'was found in the files but is not valid.")'
                write(scratchMessage,scratchFormat) ibtype_orig(fluxCount)
                call allMessage(ERROR,scratchMessage)
                call terminate()
            end select
            fluxCount = fluxCount + 1
         endif
      end do ! loop over boundaries
      !
      ! populate the adcirc arrays that are used during execution
      call populateADCIRCNativeArrays()
#endif
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      !----------------------------------------------------------------------
      end subroutine readMeshXDMF
      !----------------------------------------------------------------------

      !----------------------------------------------------------------------
      !                   S U B R O U T I N E
      !    P O P U L A T E   A D C I R C   N A T I V E   A R R A Y S
      !----------------------------------------------------------------------
      ! Use the mesh data that was read in from an external file (from
      ! whatever source and format) and populate the arrays that ADCIRC
      ! will actually use when it is running
      !----------------------------------------------------------------------
      subroutine populateADCIRCNativeArrays()
      use boundaries
      implicit none
      integer :: elevCount
      integer :: fluxCount
      integer :: i, j, k, riverval

      call setMessageSource("populateADCIRCNativeArrays")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      ! populate adcirc-native arrays
      do elevCount=1,nope
         nbdv(elevCount,1:nvdll(elevCount)) = elevationBoundaries(elevCount)%nodes
      end do
      sfCount = 1
      efCount = 1
      ifCount = 1
      ifwpCount = 1
      ! populate adcirc-native arrays
      do fluxCount=1,nbou
         select case(ibtype_orig(fluxCount))  
         case(0,1,2,10,11,12,20,21,22,30,32,52,94,102,112,122,152)
            nbvv(fluxCount,1:nvell(fluxCount)) = simpleFluxBoundaries(sfCount)%nodes
            sfCount = sfCount + 1
         case(3,13,23)
            ! populate adcirc-native arrays
            nbvv(fluxCount,1:nvell(fluxCount)) = externalFluxBoundaries(efCount)%nodes
            barlanhtr(fluxCount,1:nvell(fluxCount)) =
     &               externalFluxBoundaries(efCount)%barlanht
            barlancfspr(fluxCount,1:nvell(fluxCount)) =
     &               externalFluxBoundaries(efCount)%barlancfsp
            efCount = efCount + 1
         case(4,24,64)  ! internal barrier boundary (e.g., subgrid scale levee)
            ! populate adcirc-native arrays
            nbvv(fluxCount,1:nvell(fluxCount)) =
     &             internalFluxBoundaries(ifCount)%nodes
            ibconnr(fluxCount,1:nvell(fluxCount)) =
     &             internalFluxBoundaries(ifCount)%ibconn
            barinhtr(fluxCount,1:nvell(fluxCount)) =
     &             internalFluxBoundaries(ifCount)%barinht
            barincfsbr(fluxCount,1:nvell(fluxCount)) =
     &             internalFluxBoundaries(ifCount)%barincfsb
            barincfspr(fluxCount,1:nvell(fluxCount)) =
     &             internalFluxBoundaries(ifCount)%barincfsp
            ifCount = ifCount + 1
         case(5,25)  ! internal barrier boundary (e.g., subgrid scale levee)
            ! populate adcirc-native arrays
            nbvv(fluxCount,1:nvell(fluxCount))
     &             = internalFluxBoundariesWithPipes(ifwpCount)%nodes
            ibconnr(fluxCount,1:nvell(fluxCount))
     &             = internalFluxBoundariesWithPipes(ifwpCount)%ibconn
            barinhtr(fluxCount,1:nvell(fluxCount))
     &            = internalFluxBoundariesWithPipes(ifwpCount)%barinht
            barincfsbr(fluxCount,1:nvell(fluxCount))
     &             = internalFluxBoundariesWithPipes(ifwpCount)%barincfsb
            barincfspr(fluxCount,1:nvell(fluxCount))
     &             = internalFluxBoundariesWithPipes(ifwpCount)%barincfsp
            pipehtr(fluxCount,1:nvell(fluxCount))
     &             = internalFluxBoundariesWithPipes(ifwpCount)%pipeht
            pipecoefr(fluxCount,1:nvell(fluxCount))
     &             = internalFluxBoundariesWithPipes(ifwpCount)%pipecoef
            pipediamr(fluxCount,1:nvell(fluxCount))
     &             = internalFluxBoundariesWithPipes(ifwpCount)%pipediam
            ifwpCount = ifwpCount + 1
         case default
            write(6,*) "Unrecognized boundary type ",ibtype_orig(fluxCount)
         end select
      end do
      !
      ! initialize ibtype array
      ibtype = ibtype_orig
      ! kmd - add variables for river boundary conditions
      bndbcriver = .false.
      bcriversec=0
      totalbcrivernodes=0
      ! read in the flux boundaries
      do k = 1, nbou
         select case(ibtype_orig(k))
         case(102,112,122,152)
            bndBCRiver=.true.
            bcriversec=bcriversec+1
            bcrnvell(bcriversec)=nvell(k)
            riverval=0
            do j = 1, nvell(k)
               totalbcrivernodes=totalbcrivernodes+1
               riverval=riverval+1
               bcrnbvv(bcriversec,riverval)=nbvv(k,j)
            end do
            !
            ! in all practical respects, these river boundaries are treated
            ! the same as the corresponding boundary types that don't
            ! have a 1 in the hundreds place, so now we modify the ibtype
            ! accordingly for use in ADCIRC
            !
            ! this is why we had to save the orig value in ibtype_orig above
            ibtype(k) = ibtype_orig(k) - 100
         case default
            ! ignore other, non-baroclinic river boundaries
         end select
      end do
      !
      ! jgf51.21.11: if there were any specified flux boundaries,
      ! set an integer flag that causes ADCIRC to look for corresponding
      ! boundary conditions either in the fort.15 or fort.20
      nfluxf = 0
      do k=1,nbou
         select case(ibtype(k))
         case(2,12,22,32,52)
            nfluxf = 1 ! at least one specified flux boundary is present
            exit
         case default
            ! not a specified flux boundary, check the next one
         end select
      end do

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      !-----------------------------------------------------------------
      end subroutine populateADCIRCNativeArrays
      !-----------------------------------------------------------------


!     ------------------------------------------------------------------
!        S U B R O U T I N E    I N I T I A L I Z E   M E S H
!     ------------------------------------------------------------------
!     jgf51.21.11: Perform all the startup calculations required for
!     mesh data, regardless of the file format of the mesh file.
!     These initializations are adjustments based on parameters parsed
!     from the fort.15 control file. This must be called after the
!     ICS, NOLIFA, H0, SLAM0 and SFEA0 parameters have been read from
!     the fort.15 control file and before the kdtree computations that
!     are used to locate the elements that contain the stations.
!     ------------------------------------------------------------------
      subroutine initializeMesh()
      use global, only : h0, nolifa, nbfr, ncor, 
     &                   ntip, ifsprots
      use boundaries

      implicit none

      real(8) x1,x2,x3,x4,y1,y2,y3,y4
      real(8) avgxy,dif1r,dif2r,dif3r
!     v49.48.02 tcm added these local variable for computing
!     multiple values like element radius and areas, etc...
      real(8) x2mx1,x3mx2,x1mx3,y2my1,y3my2,y1my3
      real(8) Lon21,Lon32,Lon13
      real(8) lened1,lened2,lened3
      !real(8) :: slamr, sfear
      real(8) :: slamr1, sfear1, slamr2, sfear2, slamr3, sfear3
      integer :: i, j
      integer :: SFLAG, NMM(3)
      REAL (8), dimension(3):: LONVE, LATVE

      call setMessageSource("initializeMesh")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      ! jgf51.21.11: Initialize the comment line data for use in
      ! non portable binary output files.
      do i=1,20
         j=(i-1)*4+1
         aid4(i)=agrid(j:j+3)
      end do
      do i=1,10
         j=(i-1)*8+1
         aid8(i)=agrid(j:j+7)
      end do
!
!      F I N A L I Z E   N O D E   T A B L E
!
      SLAM0=SLAM0*DEG2RAD
      SFEA0=SFEA0*DEG2RAD
    
!...  IF ICS=1 INPUT NODAL COORDINATES AND BATHYMETRY FROM UNIT 14
!.... IF EITHER NTIP=1 OR NCOR=1, COMPUTE THE INVERSE CPP PROJECTION

      SELECT CASE( ICS )
      ! if the original nodal coordinates were cartesian
      CASE(1)
         ! we already read the original coordinates into slam(:) and sfea(:)
         ! so copy these data to x() and y()
         x = slam
         y = sfea
         ! run the inverse cpp projection on the original x and y
         ! cartesian coordinates to get the lon lat slam and sfea
         DO I=1,NP
            IF((NTIP.GE.1).OR.(NCOR.EQ.1)) THEN
               CALL INVCP(X(I),Y(I),SLAM(I),SFEA(I),SLAM0,SFEA0)
            ENDIF
         END DO
      
      CASE(2,20:24)
      ! if mesh coordinates were in lon lat degrees then compute the projection
      ! and the correction factors
         slam = deg2rad * slam
         sfea = deg2rad * sfea
         do i = 1,np
            call CYLINDERMAP(x(i),y(i),slam(i),sfea(i),slam0,sfea0,ics)
         enddo
          
      END SELECT

      IF ( IFSPROTS .EQ. 1 ) THEN
!        Allocate arrays to be used in coordinate rotations
         CALL ALLOC_SPROTCOORS_MEM( ) ;

!        Read the rotation matrix
         CALL READFORTROTS( ) ; 
        
!        Get the variables assoicated with rotated coordinates
!        The X and Y vectors are edited in this subroutine (rotated)
         CALL COMP_ROTSPCOORD_VARS( ICS ) ; 

      END IF
!...
!...  if wetting and drying will not be used, make sure all bathymetric
!...  depths are > or = to h0.
      if ( (nolifa.eq.0).or.(nolifa.eq.1) ) then
         do i=1,np
            if(dp(i).lt.h0) dp(i)=h0
         end do
      endif
 
!     check that sufficient accuracy is provided by the code to handle
!     the resolution of the input mesh (i.e., very fine meshes may have
!     nodes whose coordinates only differ in the last couple decimal
!     places ... these differences may be lost due to round off error
!     if there is not enough machine precision
      PROCELM: DO I=1,NE

        INNER: DO
          SELECT CASE( ICS )
          CASE (1)
!     c Cartesian coordinate
             EXV(1,I) = X(NM(I,1))
             EXV(2,I) = X(NM(I,2))
             EXV(3,I) = X(NM(I,3))
             EYV(1,I) = Y(NM(I,1))
             EYV(2,I) = Y(NM(I,2))
             EYV(3,I) = Y(NM(I,3)) 

          CASE (2,20:24)
!     c Spherical coordinates c!
           
!     lon (in rad)
             LONVE(1) = SLAM( NM(I,1) ) ;
             LONVE(2) = SLAM( NM(I,2) ) ;
             LONVE(3) = SLAM( NM(I,3) ) ;
           
!     lat (in rad)
              LATVE(1) = SFEA( NM(I,1) ) ;
              LATVE(2) = SFEA( NM(I,2) ) ;
              LATVE(3) = SFEA( NM(I,3) ) ;           
           
              EXV => EXVO ; 
              EYV => EYVO ;           
              CALL CAL_ELXV_SPCOOR( EXV(:,I), 
     &             EYV(:,I), LONVE, LATVE, ICS, SFLAG ) ;
              
              IF ( IFSPROTS .EQ. 1 ) THEN
!     lon in rotated coordinates (in rad)
                 LONVE(1) = SLAMR( NM(I,1) ) ;
                 LONVE(2) = SLAMR( NM(I,2) ) ;
                 LONVE(3) = SLAMR( NM(I,3) ) ;
                 
!     lat in rotated coordinates (in rad)
                 LATVE(1) = SFEAR( NM(I,1) ) ;
                 LATVE(2) = SFEAR( NM(I,2) ) ;
                 LATVE(3) = SFEAR( NM(I,3) ) ;
                 
                 EXV => EXVR ; 
                 EYV => EYVR ;              
                 CALL CAL_ELXV_SPCOOR( EXV(:,I), 
     &                EYV(:,I), LONVE, LATVE, ICS, SFLAG ) ;             
              END IF 

           END SELECT
        
           IF ( IFSPROTS .EQ. 1 ) THEN        
              EXV => EXVR ; 
              EYV => EYVR ;              
           END IF
        
           X1 = EXV(1,I) ; 
           X2 = EXV(2,I) ;
           X3 = EXV(3,I) ;
           Y1 = EYV(1,I) ;
           Y2 = EYV(2,I) ;
           Y3 = EYV(3,I) ;        

!...v49.48.02 tcm added these to reuse the computed values
           X2mX1 = (X2 - X1);
           X3mX2 = (X3 - X2);
           X1mX3 = (X1 - X3);
           Y2mY1 = (Y2 - Y1);
           Y3mY2 = (Y3 - Y2);
           Y1mY3 = (Y1 - Y3);

!.....     For later use for FDX1, FDY1 in gwce, momentum 
           FDXE(1,I) = -Y3mY2 ; !c (Y(NM2)-Y(NM3))*SFacAvg !b1
           FDXE(2,I) = -Y1mY3 ; !c (Y(NM3)-Y(NM1))*SFacAvg !b2
           FDXE(3,I) = -Y2mY1 ; !c (Y(NM1)-Y(NM2))*SFacAvg !b3
           FDYE(1,I) =  X3mX2 ; !c  X(NM3)-X(NM2) !a1
           FDYE(2,I) =  X1mX3 ; !c  X(NM1)-X(NM3) !a2
           FDYE(3,I) =  X2mX1 ; !c  X(NM2)-X(NM1) !a3

           AVGXY=(ABS(X1)+ABS(X2)+ABS(X3)+ABS(Y1)+ABS(Y2)+ABS(Y3))/6.d0
! v49.48.02 tcm -- These are the lengths of each element edge
           LENED1=((X2mX1)**2+(Y2mY1)**2)**0.5d0
           LENED2=((X3mX2)**2+(Y3mY2)**2)**0.5d0
           LENED3=((X1mX3)**2+(Y1mY3)**2)**0.5d0
! v49.48.02 rewrote to use the edge lengths which get used below also
           DIF1R=AVGXY/LENED1
           DIF2R=AVGXY/LENED2
           DIF3R=AVGXY/LENED3
           DIF1R=LOG10(DIF1R)
           DIF2R=LOG10(DIF2R)
           DIF3R=LOG10(DIF3R)

           IF((DIF1R.GT.NPREC).OR.(DIF2R.GT.NPREC).OR.(DIF3R.GT.NPREC)) THEN
              write(scratchMessage,9898) I
              call allMessage(WARNING, scratchMessage)
              WRITE(16,9898) I
 9898         FORMAT('!!!!!!!!!!  WARNING  !!!!!!!!!' ,
     &                'IF THE GRID COORDINATES HAVE 32 BITS ',
     &                '(APPROX 7 DIGITS) OF PRECISION ',
     &                'A ROBUST MODEL SOLUTION CAN NOT BE ',
     &                'GUARANTEED AT ELEMENT NO. ',I0,
     &                'MORE PRECISION MUST BE USED IN THE GRID.')
           ENDIF

!          compute and store 2 x element areas
!          v49.48.02 tcm rewrote 2 x area to reuse values already being used
           AREAS(I)=(X1mX3)*(-Y3mY2)+(X3mX2)*(Y1mY3)
  
!           IF ( IFSPROTS.EQ.1 .and. AREAS(I).lt.0d0 ) THEN
           IF ( AREAS(I).lt.0d0 ) THEN
!             reorder and start loop again
              NMM = NM(I,:)
              NM(I,2) = NMM(3); NM(I,3) = NMM(2);
           ELSE 
              EXIT INNER
           ENDIF
         ENDDO INNER
    
!...     v49.48.02 tcm -- compute element barycenters and element radius
         x1 = X( NM(I,1) )  ; 
         x2 = X( NM(I,2) )  ;
         x3 = X( NM(I,3) )  ;
         y1 = Y( NM(I,1) )  ;
         y2 = Y( NM(I,2) )  ;
         y3 = Y( NM(I,3) )  ;

         bcxy(1,i) = (x1+x2+x3)/3.d0
         bcxy(2,i) = (y1+y2+y3)/3.d0

!...     compute the radius of the circle that circumscribes
!...     the element then scale it by 50% larger to allow for
!...     a buffer later on
         RMAX(I) = 1.5D0*(LENED1*LENED2*LENED3)/(2.D0*AREAS(I))
      END DO PROCELM

!...  IF ICS=1 SET THE SFAC VECTOR EQUAL TO UNITY
!...  IF ICS=2 COMPUTE THE SFAC VECTOR IN ORDER TO ADJUST EQUATIONS TO CPP
!...  COORDINATES
!...  
      IF ( IFSPROTS .EQ. 1 ) THEN
         CALL COMPUTE_CYLINPROJ_SFAC( SLAMR, SFEAR ) ; 
      ELSE
         CALL COMPUTE_CYLINPROJ_SFAC( SLAM, SFEA ) ; 
      END IF

! move to this location by DW
!
!      F I N A L I Z E   E L E M E N T   T A B L E
!
!...
!...  Compute neighbor tables.
!...
!     estimate the number of neighbor nodes around any node = number of
!     elements containing that node.  This is correct for non-boundary
!     nodes and one too small for boundary nodes
      dp_orig = dp
      NNeigh(:) = 0
      DO I=1,NE
        NNeigh(NM(I,1))=NNeigh(NM(I,1))+1
        NNeigh(NM(I,2))=NNeigh(NM(I,2))+1
        NNeigh(NM(I,3))=NNeigh(NM(I,3))+1
      ENDDO
!     determine the maximum NNeigh
      MNei=maxval(nneigh)
 
!     estimate the maximum array space needed for the neighbor table by
!     increasing this number by 2, to provide array space for the node itself
!     and in case the maximum number of nodes occurs at a boundary node
      MNei = MNei+2
!     allocate space for neighbor tables
      call allocateNeighborArrays()
!     compute the neighbor table and redo NNeigh array
      call allMessage(INFO,'THE NEIGHBOR TABLE IS BEING COMPUTED.')
      call neighb()
      write(scratchMessage,1195) NEIMIN,NEIMAX,NEIMAX
      call allMessage(INFO, scratchMessage)
 1195 FORMAT('THE NEIGHBOR TABLE IS COMPLETED. ',
     &   'THE MINIMUM NUMBER OF NEIGHBORS FOR ANY NODE = ',i0,
     &   '. 1+THE MAXIMUM NUMBER OF NEIGHBORS FOR ANY NODE = ',i0,
     &   '. THE PARAMETER MNEI CAN BE SET AS SMALL AS ',i0,'.')

!     check to ensure that the nodal ordering in the connectivity table is
!     counter clockwise around the elements
      do i=1,ne
         if (areas(i).lt.0.0d0) then
            write(scratchMessage,9899) I
            call allMessage(ERROR,scratchMessage)
 9899      format('!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &      'IN THE CONNECTIVITY TABLE, THE NODES AROUND',
     &      ' ELEMENT ',I0,' MUST BE SPECIFIED IN COUNTERCLOCKWISE ORDER',
     &      ' - CHECK INPUT ',
     &      '!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!')
            call terminate()
         endif
      end do
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
!     ------------------------------------------------------------------
      end subroutine initializeMesh
!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
!       S U B R O U T I N E    I N I T I A L I Z E   B O U N D A R I E S
!     ------------------------------------------------------------------
!     jgf51.21.11: Perform all the finishing calculations required for
!     boundary data, regardless of the file format of the mesh file.
!     These finalizations are based on parameters parsed from the
!     fort.15 control file, including ANGINN.
!     ------------------------------------------------------------------
      subroutine initializeBoundaries()
      use global, only : nbfr, c2ddi, ifsprots
      use constants, only: pi
      use boundaries
#ifdef CMPI
      USE MESSENGER, ONLY : msg_imax  !sb msg_imax added 10/13/2022
#endif
      implicit none
      integer  nprbi, nbvk
      integer i, j, k, n, ick, iprbi !local loop counters
      character(len=80) :: exeTerm9973 ! execution terminated
      real(8) xli,xlj,xlk,yli,ylk,ylj
      real(8) xgi,xgj,xgk,ygi,ygk,ygj
      real(8) x1,x2,x3,x4,y1,y2,y3,y4
      real(8) xl, lonJI, lonKJ, lonIK 
      real(8) costheta,costheta1,cross,cross1
      real(8) theta,theta1
      real(8) dotvec
      real(8) vecnorm,vl1x,vl1y,vl2x,vl2y
      integer  jnmm,kmin,n1,n2,n3
      real(8) aemin,ae,aa,a1,a2,a3
      logical :: flg_iblen2   ! tcm v51.15 20130906 added for special case of
                               ! interior boundaries of length 2 in parallel jobs
      logical :: calctheta    ! .true. if theta is to be calculated at a flux boundary node
      real(8) :: ZNGFicNodeDist,ZNGFicNodeDistTemp
      logical :: ibtype52 = .false.
      real(8) :: delx,dely
      real(8) :: slamrI, slamrK, slamrJ, sfearI, sfearJ, sfearK

      call setMessageSource("initializeBoundaries")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      !
      ! initialization
      write(exeTerm9973,9973)
      !
      write(scratchMessage,1852) nope
1852  format('Total number of elevation boundary forcing segments :',
     &          i0,'.')
      call logMessage(ECHO,scratchMessage)

      write(scratchMessage,1854) neta
1854  format('Total number of elevation specified boundary nodes :',
     &        i0,'.')
      call logMessage(ECHO,scratchMessage)
!...
!...  NODE NUMBERS ON EACH ELEVATION BOUNDARY FORCING SEGMENT
!...  NOTE: In the future, if a boundary type is included for each
!...  elevation segment it would be possible to distinguish between
!...  periodic and aperiodic B.C.s This would require another NBD type
!...  array.
      JNMM=0
      DO K=1,NOPE
         write(scratchMessage,281) k,nvdll(k)
 281     format('Total number of nodes on elevation specified',
     &              ' boundary segment ',2x,i0,2x,'=',1x,i0)
         call logMessage(ECHO,scratchMessage)
         DO I=1,NVDLL(K)
            WRITE(16,1855) NBDV(K,I)
 1855      FORMAT(7X,I7)
            NBD(JNMM+I)=NBDV(K,I)
         END DO
         JNMM=JNMM+NVDLL(K)
      END DO !k=1,nope
!...
!...  SET UP TO READ IN TIME SERIES ELEVATION SPECIFIED BOUNDARY
!...  CONDITIONS IF APPROPRIATE
!...
      if((nbfr.eq.0).and.(nope.gt.0)) then
         write(scratchmessage,1871)
         call logMessage(INFO,scratchMessage)
 1871    format('Time series elevation specified values ',
     &          'will be read from unit 19.',
     &          ' Interpolation in time is done to sync the ',
     &          'elevation data with the model time step.')
      endif
!...
!...  FLOW BOUNDARY OUTPUT TO UNIT 16
!...
!...  INTERIOR NODES, LBCODE=-1, COS=0, SIN=1
!...  BOUNDARY NODES, LBCODE=LBCODEI=IBTYPE,
!...  COS & SIN DETERMINED FROM NORMAL DIRECTION IN ALL CASES, ALTHOUGH THIS
!...  INFORMATION IS ONLY USED WHEN NORMAL FLOW IS AN ESSENTIAL B.C. AND
!...  FREE TANGENTIAL SLIP IS ALLOWED.

!...  TOTAL NUMBER OF FLOW BOUNDARY SEGMENTS
      WRITE(16,1112)
      WRITE(16,1878)
1878  FORMAT(//,1X,'FLOW BOUNDARY INFORMATION ',/)
      WRITE(16,1879) NBOU
 1879 FORMAT(//,5X,'THE TOTAL NUMBER OF FLOW BOUNDARY SEGMENTS = ',i0)

!.....INPUT THE TOTAL NUMBER OF FLOW BOUNDARY NODES
      WRITE(16,1881) NVEL
1881  FORMAT(/,5X,'THE TOTAL NUMBER OF FLOW BOUNDARY NODES = ',i0)

      JGW=0
      JME=0
      NFLUXF=0
      NFLUXB=0
      NFLUXIB=0
      NFLUXIBP=0
      NFLUXIB64=0
      NFLUXRBC=0
      NFLUXGBC=0
      NVELEXT=0
      !
      ! loop over the flux boundary segments, logging information about
      ! each segment and setting up boundary arrays
      NBOULOOP: DO K=1,NBOU
!        write out flow boundary information to unit 16
         IF((IBTYPE(K).EQ.4).OR.(IBTYPE(K).EQ.24).OR.
     &      (IBTYPE(K).EQ.64)) THEN
            WRITE(16,28) K,NVELL(K),K,2*NVELL(K)
 28         FORMAT(///,5X,'TOTAL NUMBER OF PAIRS FOR FLOW BOUNDARY',
     &           ' SEGMENT',2X,i0,2X,'=',2X,i0,
     &           /,5X,'TOTAL NUMBER OF NODES FOR FLOW BOUNDARY',
     &           ' SEGMENT',2X,i0,2X,'=',2X,i0)
         ELSE
            WRITE(16,128) K,NVELL(K)
 128        FORMAT(///,5X,'TOTAL NUMBER OF NODES FOR FLOW BOUNDARY',
     &           ' SEGMENT',2X,i0,2X,'=',2X,i0)
         ENDIF
!        continue processing flow boundary information
         SELECT CASE(IBTYPE(K))
         CASE(0)
            WRITE(16,2340)
 2340       FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &           7X,'NO NORMAL FLOW AS AN ESSENTIAL B.C.',/,
     &           7X,'AND FREE TANGENTIAL SLIP',/)
         CASE(1)
            WRITE(16,2341)
 2341       FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BOUNDARY WITH:',/,
     &           7X,'NO NORMAL FLOW AS AN ESSENTIAL B.C.',/,
     &           7X,'AND FREE TANGENTIAL SLIP',/)
         CASE(2)
            NFLUXF=1
            WRITE(16,2342) 
 2342       FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &           7X,'SPECIFIED NORMAL FLOW AS AN ESSENTIAL B.C.',/,
     &           7X,'AND FREE TANGENTIAL SLIP',/)
         CASE(3)
           NFLUXB=1
           WRITE(16,2344)
 2344      FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &            7X,'A BARRIER WHICH ALLOWS FREE SURFACE',
     &               ' SUPERCRITICAL OUTFLOW',/,
     &            7X,'FROM THE DOMAIN ONCE THE BARRIER HAS BEEN',
     &               ' OVERTOPPED',/,
     &            7X,'AND FREE TANGENTIAL SLIP',/)
         CASE(4)
           NFLUXIB=1
           WRITE(16,2345)
 2345      FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BARRIER BOUNDARY:',/,
     &            7X,'WITH CROSS BARRIER FLOW TREATED AS AN ESSENTIAL ',
     &               ' NORMAL FLOW BOUNDARY CONDITION',/,
     &            7X,'WHICH LEAVES/ENTERS THE DOMAIN ON ONE SIDE OF ',
     &               ' THE BARRIER AND ENTERS/LEAVES THE DOMAIN ON THE',
     &          /,7X,'CORRESPONDING OPPOSITE SIDE OF THE BARRIER ',/,
     &            7X,'FLOW RATE AND DIRECTION ARE BASED ON BARRIER ',
     &               ' HEIGHT, SURFACE WATER ELEVATION',/,
     &            7X,'ON BOTH SIDES OF THE BARRIER, BARRIER COEFFICIEN',
     &               'T AND THE APPROPRIATE BARRIER FLOW FORMULA',/,
     &            7X,'FREE TANGENTIAL SLIP IS ALLOWED',/)
         CASE(5)
           NFLUXIB=1
           NFLUXIBP=1
           WRITE(16,2347)
 2347      FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BARRIER BOUNDARY:',/,
     &            7X,'WITH ADDITIONAL CROSS BARRIER PIPES ',
     &               'LOCATED UNDER THE CROWN ',/,
     &            7X,'CROSS BARRIER FLOW IS TREATED AS AN ESSENTIAL',
     &               ' NORMAL FLOW BOUNDARY CONDITION',/,
     &            7X,'WHICH LEAVES/ENTERS THE DOMAIN ON ONE SIDE OF ',
     &               ' THE BARRIER AND ENTERS/LEAVES THE DOMAIN ON THE',
     &          /,7X,'CORRESPONDING OPPOSITE SIDE OF THE BARRIER ',/,
     &            7X,'FLOW RATE AND DIRECTION ARE BASED ON BARRIER ',
     &               ' HEIGHT, SURFACE WATER ELEVATION',/,
     &            7X,'ON BOTH SIDES OF THE BARRIER, BARRIER COEFFICIEN',
     &               'T AND THE APPROPRIATE BARRIER FLOW FORMULA',/,
     &            7X,'IN ADDITION CROSS BARRIER PIPE FLOW RATE AND ',
     &               ' DIRECTION ARE BASED ON PIPE CROWN HEIGHT, ',/,
     &            7X,'SURFACE WATER ELEVATION ON BOTH SIDES OF THE B',
     &               'ARRIER, PIPE FRICTION COEFFICIENT, PIPE DIAMETER',
     &          /,7X,' AND THE APPROPRIATE PIPE FLOW FORMULA',/,
     &            7X,'FREE TANGENTIAL SLIP IS ALLOWED',/)
         CASE(10)
           WRITE(16,2350)
 2350      FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &            7X,'NO NORMAL FLOW AS AN ESSENTIAL B.C.',/,
     &            7X,'AND NO TANGENTIAL SLIP',/)
         CASE(11)
           WRITE(16,2351)
 2351      FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BOUNDARY WITH:',/,
     &            7X,'NO NORMAL FLOW AS AN ESSENTIAL B.C.',/,
     &            7X,'AND NO TANGENTIAL SLIP',/)
         CASE(12)
           NFLUXF=1
           WRITE(16,2352)
 2352      FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &            7X,'SPECIFIED NORMAL FLOW AS AN ESSENTIAL B.C.',/,
     &            7X,'AND NO TANGENTIAL SLIP',/)
         CASE(13)
           NFLUXB=1
           WRITE(16,2353)
 2353      FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &            7X,'A BARRIER WHICH ALLOWS FREE SURFACE',
     &               ' SUPERCRITICAL OUTFLOW',/,
     &            7X,'FROM THE DOMAIN ONCE THE BARRIER HAS BEEN',
     &               ' OVERTOPPED',/,
     &            7X,'AND NO TANGENTIAL SLIP',/)
         CASE(20)
           WRITE(16,2354)
 2354      FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &            7X,'NO NORMAL FLOW AS A NATURAL B.C.',/,
     &            7X,'AND FREE TANGENTIAL SLIP',/)
         CASE(21)
           WRITE(16,2355)
 2355      FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BOUNDARY WITH:',/,
     &            7X,'NO NORMAL FLOW AS A NATURAL B.C.',/,
     &            7X,'AND FREE TANGENTIAL SLIP',/)
         CASE(22)
           NFLUXF=1
           WRITE(16,2356)
 2356      FORMAT(5X,'THIS SEGMENT IS A EXTERNAL BOUNDARY WITH:',/,
     &            7X,'SPECIFIED NORMAL FLOW AS A NATURAL B.C.',/,
     &            7X,'AND FREE TANGENTIAL SLIP',/)
         CASE(23)
           NFLUXB=1
           WRITE(16,2357)
 2357      FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &            7X,'A BARRIER WHICH ALLOWS FREE SURFACE',
     &               ' SUPERCRITICAL OUTFLOW',/,
     &            7X,'FROM THE DOMAIN ONCE THE BARRIER HAS BEEN',
     &               ' OVERTOPPED',/,
     &            7X,' IMPLEMENTED AS A NATURAL BOUNDARY CONDITION',
     &            7X,'FREE TANGENTIAL SLIP IS ALSO ALLOWED',/)
         CASE(24)
           NFLUXIB=1
           WRITE(16,2358)
 2358      FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BARRIER BOUNDARY:',/,
     &            7X,'WITH CROSS BARRIER FLOW TREATED AS A NATURAL ',
     &               ' NORMAL FLOW BOUNDARY CONDITION',/,
     &            7X,'WHICH LEAVES/ENTERS THE DOMAIN ON ONE SIDE OF ',
     &               ' THE BARRIER AND ENTERS/LEAVES THE DOMAIN ON THE',
     &          /,7X,'CORRESPONDING OPPOSITE SIDE OF THE BARRIER ',/,
     &            7X,'FLOW RATE AND DIRECTION ARE BASED ON BARRIER ',
     &               ' HEIGHT, SURFACE WATER ELEVATION',/,
     &            7X,'ON BOTH SIDES OF THE BARRIER, BARRIER COEFFICIEN',
     &               'T AND THE APPROPRIATE BARRIER FLOW FORMULA',/,
     &            7X,'FREE TANGENTIAL SLIP IS ALLOWED',/)
         CASE(25)
           NFLUXIB=1
           NFLUXIBP=1
           WRITE(16,2359)
 2359      FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BARRIER BOUNDARY:',/,
     &            7X,'WITH ADDITIONAL CROSS BARRIER PIPES ',
     &               'LOCATED UNDER THE CROWN ',/,
     &            7X,'CROSS BARRIER FLOW IS TREATED AS A NATURAL',
     &               ' NORMAL FLOW BOUNDARY CONDITION',/,
     &            7X,'WHICH LEAVES/ENTERS THE DOMAIN ON ONE SIDE OF ',
     &               ' THE BARRIER AND ENTERS/LEAVES THE DOMAIN ON THE',
     &          /,7X,'CORRESPONDING OPPOSITE SIDE OF THE BARRIER ',/,
     &            7X,'FLOW RATE AND DIRECTION ARE BASED ON BARRIER ',
     &               ' HEIGHT, SURFACE WATER ELEVATION',/,
     &            7X,'ON BOTH SIDES OF THE BARRIER, BARRIER COEFFICIEN',
     &               'T AND THE APPROPRIATE BARRIER FLOW FORMULA',/,
     &            7X,'IN ADDITION CROSS BARRIER PIPE FLOW RATE AND ',
     &               ' DIRECTION ARE BASED ON PIPE CROWN HEIGHT, ',/,
     &            7X,'SURFACE WATER ELEVATION ON BOTH SIDES OF THE B',
     &               'ARRIER, PIPE FRICTION COEFFICIENT, PIPE DIAMETER',
     &          /,7X,'AND THE APPROPRIATE PIPE FLOW FORMULA',/,
     &            7X,'FREE TANGENTIAL SLIP IS ALLOWED',/)
         CASE(30)
           NFluxRBC=1
           WRITE(16,2360)
 2360      FORMAT(5X,'This segment is an outward radiating boundary:',/,
     &            7X,'the GWCE is forced with normal flux.  A wave ',
     &               'radiation condition',
     &            7X,'is used to relate the time derivative of the ',
     &               'normal flux to',/,
     &            7X,'the time derivative of the elevation.  The ',
     &               'momentum equations',/,
     &            7X,'are used to compute the velocity field the same ',
     &               'as for a nonboundary node.',/)
         CASE(32)
           NFluxRBC=1
           NFluxF=1
           WRITE(16,2356)
 2362      FORMAT(5X,'This segment is a combined specified normal ',
     &               'flux and outward radiating boundary.',/,
     &            7X,'The GWCE is forced with the total normal flux ',
     &               'computed by adding the specified',/,
     &            7X,'normal flux and the flux associated with the ',
     &               'outward radiating wave.',/,
     &            7X,'The latter is determine from a Sommerfeld type ',
     &               'condition, flux=celerity*wave elevation',/,
     &            7X,'The momentum equations are used to compute the',/,
     &            7X,'velocity field the same as for a nonboundary ',
     &               'node',/)
         CASE(40)
           NFluxGBC=1
           WRITE(16,2370)
 2370      FORMAT(5X,'This segment is a zero normal velocity gradient ',
     &               'boundary:',/,
     &            7X,'the GWCE is forced with normal flux',/,
     &            7X,'the momentum eqs are sacrificed in favor of ',
     &               'setting the velocity at a',/,
     &            7X,'boundary node equal to the value at a fictitious',
     &               ' point inside the domain.',/,
     &            7X,'The fictitious point lies on the inward directed',
     &               ' normal to the boundary',/,
     &            7X,'a distance equal to the distance from the ',
     &               'boundary node to its farthest',/,
     &            7X,'neighbor. This should ensure that the fictitious',
     &               ' point does not',/,
     &            7X,'fall into an element that contains the boundary ',
     &               'node.',/,
     &            7X,'The velocity at the fictitious point is ',
     &               'determined by interpolation.',/)
         CASE(41)
           NFluxGBC=1
           WRITE(16,2371)
 2371      FORMAT(5X,'This segment is a zero normal velocity gradient ',
     &               'boundary:',/,
     &            7X,'the GWCE is forced with normal flux',/,
     &            7X,'the momentum eqs are sacrificed in favor of eqs ',
     &               'that set the velocity',/,
     &            7X,'gradient normal to the boundary equal to zero ',
     &               'in the Galerkin sense.',/)
         CASE(52)
            NFLUXF=1
            IBTYPE52=.TRUE.
            WRITE(16,6532)
 6532       FORMAT(5X,'This segment uses a specified periodic normal ',
     &          'flux as a natural b.c. combined with outward ',
     &          7X,'radiating boundary. The GWCE is forced with the ,'
     &          'total normal flux computed by adding the',/,
     &          7X,'specified normal flux and the flux associated ',
     &          'with the outward radiating wave.',/,
     &          7X,'The latter is determined from ',
     &          'flux=celerity*wave elevation',/,
     &          7X,'The momentum equations are used to compute the',/,
     &          7X,'velocity field the same as for a nonboundary ',
     &          'node.',/)
         CASE(64)
           NFLUXIB=1
           NFLUXIB64=1
           WRITE(16,6533)
 6533      FORMAT(5X,'THIS SEGMENT IS AN VERTICAL ELEMENT WALL BOUNDARY:',/,
     &            7X,'EACH NODE PAIR IS TREATED AS IF THEY ARE A SINGL',
     &               'E NODE IF THE WATER HEIGHTS ON BOTH SIDES',/,
     &            7X,'EXCEED THE BARRIER HEIGHTS.',/,
     &            7X,'OTHERWISE, WITH CROSS BARRIER FLOW TREATED ',
     &               ' AS A NATURAL NORMAL FLOW BOUNDARY CONDITION',/,
     &            7X,'WHICH LEAVES/ENTERS THE DOMAIN ON ONE SIDE OF ',
     &               ' THE BARRIER AND ENTERS/LEAVES THE DOMAIN ON THE',
     &          /,7X,'CORRESPONDING OPPOSITE SIDE OF THE BARRIER ',/,
     &            7X,'FLOW RATE AND DIRECTION ARE BASED ON BARRIER ',
     &               ' HEIGHT, SURFACE WATER ELEVATION',/,
     &            7X,'ON BOTH SIDES OF THE BARRIER, BARRIER COEFFICIEN',
     &               'T AND THE APPROPRIATE BARRIER FLOW FORMULA',/,
     &            7X,'FREE TANGENTIAL SLIP IS ALLOWED.',/)
         END SELECT
         !
         ! set the constants used in setting up the boundary arrays
         select case(ibtype(k))
         case(3,13,23)
            nprbi = 1
            npipe = 0
         case(4,24,64)
            nprbi = 2
            npipe = 0
         case(5,25)
            nprbi = 2
            npipe = 1
         case default
            nprbi = 1
            npipe = 0
         end select

!...  PROCESS INFORMATION FOR VARIOUS TYPES OF FLOW BOUNDARY SEGMENTS

         LOOPPRBI: DO IPRBI=1,NPRBI

!...  LOAD PAIRED NODES INTO PRIMARY PROCESSING VECTORS AND RESET
!...  CONNECTING NODES FOR BACK FACE THUS BACK/CONNECTING NODES ARE
!...  BEING LOADED AS PRIMARY NODES AND FRONT NODES ARE RELOADED AS
!...  CONNECTING NODES NOTE THAT THE CLOCKWISE ORIENTATION OF ISLAND
!...  TYPE BOUNDARIES IS BEING MAINTAINED WHEN BACK NODES ARE RELOADED
!...  AS PRIMARY NODES ADDITIONAL INTERNAL BARRIER BOUNDARY INFORMATION
!...  IS ALSO RESET

            IF(IPRBI.EQ.2) THEN
               DO I=1,NVELL(K)
                  NTRAN1(I)=NBVV(K,I)
                  NTRAN2(I)=IBCONNR(K,I)
                  BTRAN3(I)=BARINHTR(K,I)
                  BTRAN4(I)=BARINCFSBR(K,I)
                  BTRAN5(I)=BARINCFSPR(K,I)
                  IF(NPIPE.EQ.1) THEN
                     BTRAN6(I)=PIPEHTR(K,I)
                     BTRAN7(I)=PIPECOEFR(K,I)
                     BTRAN8(I)=PIPEDIAMR(K,I)
                  ENDIF
               END DO  !I=1,NVELL(K)

               DO I=1,NVELL(K)
                  NBVV(K,I)=NTRAN2(NVELL(K)+1-I)
                  IBCONNR(K,I)=NTRAN1(NVELL(K)+1-I)
                  BARINHTR(K,I)=BTRAN3(NVELL(K)+1-I)
                  BARINCFSBR(K,I)=BTRAN4(NVELL(K)+1-I)
                  BARINCFSPR(K,I)=BTRAN5(NVELL(K)+1-I)
                  IF(NPIPE.EQ.1) THEN
                     PIPEHTR(K,I)=BTRAN6(NVELL(K)+1-I)
                     PIPECOEFR(K,I)=BTRAN7(NVELL(K)+1-I)
                     PIPEDIAMR(K,I)=BTRAN8(NVELL(K)+1-I)
                  ENDIF
               END DO  !I=1,NVELL(K)
            ENDIF

!...        WRITE OUT ADDITIONAL HEADER FOR INTERNAL BARRIER BOUNDARIES
            IF((IBTYPE(K).EQ.4).OR.(IBTYPE(K).EQ.24).OR.
     &         (IBTYPE(K).EQ.64)) THEN
               IF(IPRBI.EQ.1) THEN
                  WRITE(16,1842)
 1842            FORMAT(/,5X,'FRONT FACE OF INTERNAL BARRIER BOUNDARY',/)
               ELSE
                  WRITE(16,1843)
 1843            FORMAT(/,5X,'BACK FACE OF INTERNAL BARRIER BOUNDARY',/)
               ENDIF
            ENDIF

!...        WRITE OUT ADDITIONAL HEADER FOR INTERNAL BARRIER BOUNDARIES WITH
!...        CROSS BARRIER PIPES
            IF((IBTYPE(K).EQ.5).OR.(IBTYPE(K).EQ.25)) THEN
               IF(IPRBI.EQ.1) THEN
                  WRITE(16,1844)
 1844             FORMAT(/,5X,'FRONT FACE OF INTERNAL BARRIER BOUNDARY',
     &                 ' WITH CROSS BARRIER PIPES',/)
               ELSE
                  WRITE(16,1845)
 1845             FORMAT(/,5X,'BACK FACE OF INTERNAL BARRIER BOUNDARY',
     &                 ' WITH CROSS BARRIER PIPES',/)
               ENDIF
            ENDIF

!...        WRITE OUT GENERAL HEADER FOR BOUNDARY INFORMATION
            WRITE(16,1841)
 1841      FORMAT('    JGW    JME    ME2GW   NODE #  BNDRY CODE   INNER',
     &      ' ANGLE',7X,'COS',13X,'SIN',9X,'0.667*BNDRY LEN',/)

!...        COMPLETE THE BOUNDARY ARRAY FOR THE Kth FLOW BOUNDARY SEGMENT

            FLG_IBLEN2 = .FALSE.   !ADDED TCM V51.15 20130906
            NBVV(K,0)=NBVV(K,1) !UNCLOSED EXTERNAL
            IF((IBTYPE(K).EQ.1).OR.(IBTYPE(K).EQ.11).OR.
     &            (IBTYPE(K).EQ.21)) THEN
               IF(NBVV(K,NVELL(K)).NE.NBVV(K,1)) THEN !CLOSE AN UNCLOSED INTERNAL

                  IF (NVELL(K).EQ.2) FLG_IBLEN2 = .TRUE.   !ADDED TCM V51.15 20130906

                  NVELL(K)=NVELL(K)+1
                  NBVV(K,NVELL(K))=NBVV(K,1)
                  !tcm 20130906 -- Closing can be problematic for parallel jobs
                  ! where it is possible that only part of the internal boundary
                  ! is in the domain
               ENDIF
            ENDIF
            IF(NBVV(K,NVELL(K)).EQ.NBVV(K,1)) THEN !CLOSED EXTERNAL OR INTERNAL
               NBVV(K,0)=NBVV(K,NVELL(K)-1)
            ENDIF
            NBVV(K,NVELL(K)+1)=NBVV(K,NVELL(K))

!...        PUT BOUNDARY INFORMATION INTO 2 TYPES OF ARRAYS, ONE FOR THE GWCE
!...        B.C. AND ONE FOR THE MOMENTUM EQUATION B.C.
!
!...        THE GWCE ARRAYS INCLUDE EVERY NODE IN THE UNIT 14 FILE, I.E.,
!...        NODES ARE REPEATED WHERE SPECIFIED NORMAL FLOW AND NO NORMAL FLOW
!...        BOUNDARIES MEET AND AT THE BEGINNING AND END OF CLOSED EXTERNAL
!...        BOUNDARIES AND ISLANDS.
!
!...        THE MOMENTUM EQUATION ARRAYS ARE KEYED TO THE GWCE ARRAYS VIA THE
!...        ARRAY ME2GW WHICH INDICATES THE LOCATION IN THE GWCE ARRAYS THAT
!...        THE APPROPRIATE M.E. VALUE LIES.
!
!...        THE M.E. ARRAYS DO NOT REPEAT NODES THAT ARE DUPLICATED IN THE
!...        UNIT 14 FILE, I.E., WHEN SPECIFIED NORMAL FLOW AND NO NORMAL FLOW
!...        BOUNDARIES MEET, THE SPECIFIED NORMAL FLOW BOUNDARY CONDITION
!...        TAKES PRECEDENT.  ALSO THE BEGINNING AND ENDING NODES OF CLOSED
!...        EXTERNAL AND ISLAND BOUNDARIES ARE NOT REPEATED.

            NVELL_LOOP: DO I=1,NVELL(K)
!...           SET UP THE GWCE BOUNDARY ARRAYS WHICH CONSIST OF
!...           BOUNDARY NODE NUMBERS
!...           BOUNDARY CODES
!...           0.66667*LENGTH OF EACH BOUNDARY SEGMENT.  NOTE, THE LENGTH OF THE
!...           LAST BOUNDARY SEGMENT ON EACH BOUNDARY SHOULD BE ZERO
!              jgf46.21 Added IBTYPE=52.
               JGW=JGW+1
               select case(ibtype(k))
               case(0,10,20,2,12,22,52,3,13,23,30,32,40,41)
                  NVELEXT=NVELEXT+1
               end select
               NBVI = NBVV(K,I)
               NBVJ = NBVV(K,I+1)
               NBVK = NBVV(K,I-1)
               NBV(JGW) = NBVI
               IF (NBVJ.gt.ubound(X,1)) THEN
                  write(16,*) 'NBVJ', NBVJ, '> bounds of X',I
               ENDIF
               IF (NBVI.gt.ubound(X,1)) THEN
                  write(16,*) 'NBVI', NBVI, '> bounds of X',I
               ENDIF
               XGJ = X(NBVJ)
               XGI = X(NBVI)
               XGK = X(NBVK)
               YGJ = Y(NBVJ)
               YGI = Y(NBVI)
               YGK = Y(NBVK)
               DELX = XGJ - XGI
               DELY = YGJ - YGI
               BNDLEN2O3(JGW) = 2.D0*(SQRT(DELX*DELX+DELY*DELY))/3.D0

!...           COMPUTE THE INCLUDED ANGLE AND TEST TO DETERMINE WHETHER TO ZERO
!...           TANGENTIAL VELOCITIES

!...           NOTE:.IMPLEMENTATION FOR ICS=2 REQUIRES COMPUTING ALL COORDINATES
!...           IN A LOCALIZED SYSTEM (I.E. THE TRANSFORMATION IS CENTERED AT
!...           X0,Y0)

               IF(ICS.EQ.1) THEN
                  XLI = XGI
                  XLJ = XGJ
                  XLK = XGK
                  YLI = YGI
                  YLJ = YGJ
                  YLK = YGK
               ELSE
                  IF ( IFSPROTS .EQ. 1 ) THEN
                     CALL CYLINDERMAP( XLI, YLI, SLAMR(NBVI),
     &                    SFEAR(NBVI), SLAMR(NBVI),
     &                    SFEAR(NBVI), ICS ) ;
                     CALL CYLINDERMAP(XLJ, YLJ, SLAMR(NBVJ),
     &                    SFEAR(NBVJ), SLAMR(NBVI),
     &                    SFEAR(NBVI), ICS ) ;
                     CALL CYLINDERMAP(XLK, YLK, SLAMR(NBVK),
     &                    SFEAR(NBVK), SLAMR(NBVK),
     &                    SFEAR(NBVI), ICS)
                  ELSE
                     CALL CYLINDERMAP( XLI, YLI, SLAM(NBVI),
     &                    SFEA(NBVI), SLAM(NBVI),
     &                    SFEA(NBVI), ICS ) ;
                     CALL CYLINDERMAP(XLJ, YLJ, SLAM(NBVJ),
     &                    SFEA(NBVJ), SLAM(NBVI),
     &                    SFEA(NBVI), ICS ) ;
                     CALL CYLINDERMAP(XLK, YLK, SLAM(NBVK),
     &                    SFEA(NBVK), SLAM(NBVI),
     &                    SFEA(NBVI), ICS)
                 ENDIF
               ENDIF

!...           NOTE: INTERIOR ANGLE AT ENDS OF BOUNDARIES MUST BE EQUAL, EITHER:
!...           A FICTITIOUSLY LARGE VALUE IF THE BOUNDARY IS NOT CLOSED OR A TRUE
!...           VALUE IF THE BOUNDARY IS CLOSED

               calcTheta = .true.
               ! jgf51.15: special treatment at the start and
               ! end of the boundary
               IF (I.eq.1) THEN ! start of boundary
                  ! jgf51.15: if ADCIRC previously had to close an
                  ! unclosed external boundary
                  IF (NBVV(K,1).EQ.NBVV(K,0)) THEN
                     THETA=-9999999.d0    ! set to a fictitiously large val
                     CROSS=0.d0
                     COSTHETA=COSTSET     !COSTSET = COS(ANGINN*DEG2RAD)
                     calcTheta = .false.
                  ENDIF
               ENDIF
               IF (I.eq.NVELL(K)) THEN ! end of boundary
                  ! set properties equal to those from start of boundary
                  THETA=THETA1
                  CROSS=CROSS1
                  COSTHETA=COSTHETA1
                  calcTheta = .false.
               ENDIF
               IF (calcTheta.eqv..true.) THEN
                  ! jgf51.15: calculate the properties of the interior
                  ! angle of the boundary
                  VL1X = XLK-XLI
                  VL1Y = YLK-YLI
                  VL2X = XLJ-XLI
                  VL2Y = YLJ-YLI
                  DOTVEC = VL1X*VL2X+VL1Y*VL2Y
                  VECNORM = (SQRT(VL1X**2 + VL1Y**2))*
     &                      (SQRT(VL2X**2 + VL2Y**2))
                  COSTHETA=DOTVEC/VECNORM
                  IF(COSTHETA.GT.1.0d0) COSTHETA=1.0d0
                  IF(COSTHETA.LT.-1.0d0) COSTHETA=-1.0d0
                  THETA = RAD2DEG*ACOS(COSTHETA)
                  CROSS = -VL1X*VL2Y+VL2X*VL1Y
                  IF(CROSS.LT.0) THETA = 360.d0 - THETA
               ENDIF
               IF (I.eq.1) THEN
                  ! jgf51.15: record the properties at the
                  ! start of the boundary, for re-use at the end
                  THETA1=THETA
                  CROSS1=CROSS
                  COSTHETA1=COSTHETA
               ENDIF

!...           CHECK WHETHER ANGLE IS LESS THAN MINIMUM ANGLE, IF SO CHANGE THE
!...           BOUNDARY CODE TO ZERO TANGENTIAL VELOCITIES 

               ! jgf53.dev: jgfdebug: FIXME: The 3D code does not appear
               ! to support boundary types with no tangential slip, so
               ! as an interim measure, I am making this boundary type
               ! modification a 2DDI-only feature.
               if (C2DDI.eqv..true.) then
                  LBCODEI(JGW)=IBTYPE(K)
                  IF((COSTHETA.GT.COSTSET).AND.(CROSS.GT.0.0)) THEN
                     IF(IBTYPE(K).EQ.0) LBCODEI(JGW)=10
                     IF(IBTYPE(K).EQ.1) LBCODEI(JGW)=11
                     IF(IBTYPE(K).EQ.2) LBCODEI(JGW)=12
                     IF(IBTYPE(K).EQ.3) LBCODEI(JGW)=13
                     IF((IBTYPE(K).GE.0).AND.(IBTYPE(K).LE.3)) THEN
                        WRITE(16,1856) NBVI,THETA
 1856                   FORMAT(2X,i0,4X,'THE INNER ANGLE = ',F8.2,1X,
     &                    'TANGENTIAL SLIP WILL BE ZEROED')
                     ENDIF
                  ENDIF
               endif

!...           COMPUTE COS AND SIN OF OUTWARD NORMAL REGARDLESS OF BOUNDARY TYPE
               ! TCM V51.15 20130906 ADDED TEST TO HANDLE SPECIAL CASE WHEN THE ORIGINAL
               ! INTERIOR BOUNDARY HAD LENGTH 2 AND IBTYPE=1,11, OR 21 AND WAS FORCED CLOSED
               IF (FLG_IBLEN2 .EQV. .TRUE.) THEN
                  XGJ = XGI
                  YGJ = YGI
               ENDIF
               XL = SQRT((XGK-XGJ)**2+(YGK-YGJ)**2)
               CSII(JGW) = SFMX(NBVI)*(YGJ-YGK)/XL
               SIII(JGW) = (XGK-XGJ)/XL

!...           Compute the location of and the element containing the fictitious
!...           node used for a zero normal velocity gradient boundary condition
!...           (type 40)

               IF(IBType(K).EQ.40) Then
                  ZNGFicNodeDist=0.d0

                  DO N=2,NNeigh(NBVI)
                     ZNGFicNodeDistTemp=SQRT(
     &                     (X(NeiTab(NBVI,N))-X(NBVI))**2
     &                    +(Y(NeiTab(NBVI,N))-Y(NBVI))**2 )

                     IF(ZNGFicNodeDist.LT.ZNGFicNodeDistTemp)
     &                    ZNGFicNodeDist=ZNGFicNodeDistTemp

                  ENDDO  !N=2,NNeigh(NBVI)

                  ZNGFicNodeDist=1.001d0*ZNGFicNodeDist
                  X4=X(NBVI)-ZNGFicNodeDist*CSII(JGW)
                  Y4=Y(NBVI)-ZNGFicNodeDist*SIII(JGW)
                  AEMIN=1.0E+25

                  DO N=1,NE
                     N1=NM(N,1)
                     N2=NM(N,2)
                     N3=NM(N,3)
                     X1=X(N1)
                     X2=X(N2)
                     X3=X(N3)
                     Y1=Y(N1)
                     Y2=Y(N2)
                     Y3=Y(N3)
                     A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
                     A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
                     A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
                     AA=ABS(A1)+ABS(A2)+ABS(A3)
                     AE=ABS(AA-AREAS(N))/AREAS(N)
                     IF(AE.LT.AEMIN) THEN
                        AEMIN=AE
                        NEleZNG(JGW)=N
                     ENDIF
                  ENDDO  !N=1,NE

                  IF(AEMIN.GT.1.0E-5) THEN
                     write(scratchMessage,9770) nbvi
                     call allMessage(WARNING,scratchMessage)
 9770                FORMAT('!!!!!!!!!!  ERROR !!!!!!!!! ',
     &                'Zero Normal Velocity Gradient Boundary Node ',
     &                'Number ',i0,
     &                'does not appear to have a fictitious node ',
     &                'located within the domain. ',
     &                'This should be checked.')
                  ENDIF

!...              Pre-compute the information required to interpolate at zero normal
!...              gradient fictitious nodes

                  N1=NM(NEleZNG(JGW),1)
                  N2=NM(NEleZNG(JGW),2)
                  N3=NM(NEleZNG(JGW),3)
                  X1=X(N1)
                  X2=X(N2)
                  X3=X(N3)
                  Y1=Y(N1)
                  Y2=Y(N2)
                  Y3=Y(N3)
                  ZNGIF1(JGW)=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))
     &                 /AREAS(NEleZNG(JGW))
                  ZNGIF2(JGW)=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))
     &                 /AREAS(NEleZNG(JGW))
                  ZNGIF3(JGW)=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))
     &                 /AREAS(NEleZNG(JGW))
!             WRITE(2,*) 'For Boundary Node # ',NBVV(K,I)
!             WRITE(2,*) 'The fictitious point is ',X4,Y4
!             WRITE(2,*) 'which is located in element # ',NEleZNG(JGW)
!             WRITE(2,*) 'The interpolating factors are ',ZNGIF1(JGW),
!    &                                           ZNGIF2(JGW),ZNGIF3(JGW)
!             WRITE(2,*) ' '
               ENDIF

!...
!...  SET UP THE MOMENTUM EQUATION BOUNDARY ARRAY WHICH CONSISTS OF A
!...  KEY TO THE GWCE BOUNDARY CONDITION ARRAY
!...
               IF(I.EQ.1) THEN  !DEAL WITH FIRST NODE IN L.B. SEG
                  IF(JGW.EQ.1) THEN !VERY FIRST L.B. SEG
                     JME=JME+1  !M.E. USES IT
                     ME2GW(JME)=JGW
                  ENDIF
                  IF(JGW.NE.1) THEN
                     IF(NBV(JGW).NE.NBV(JGW-1)) THEN !L.B. SEGS DON'T OVERLAP
                        JME=JME+1 !M.E. USES IT
                        ME2GW(JME)=JGW
                     ENDIF
                     IF(NBV(JGW).EQ.NBV(JGW-1)) THEN !L.B. SEGS OVERLAP
                        IF((LBCODEI(JGW).EQ.2) .OR. ! M.E. USES IT ONLY
     &                       (LBCODEI(JGW).EQ.12).OR. ! IF IT IS
     &                       (LBCODEI(JGW).EQ.22).OR. ! SPECIFIED FLOW,
     &                       (LBCODEI(JGW).EQ.3) .OR. ! AN OVERFLOW BARRIER
     &                       (LBCODEI(JGW).EQ.13).OR. ! OR A RADIATION
     &                       (LBCODEI(JGW).EQ.23).OR. ! BOUNDARY
     &                       (LBCODEI(JGW).EQ.30).OR.
     &                       (LBCODEI(JGW).EQ.32).OR.
     &                       (LBCODEI(JGW).EQ.40).OR.
     &                       (LBCODEI(JGW).EQ.41).OR.
     &                       (LBCODEI(JGW).EQ.52).OR.
     &                       (LBCODEI(JGW).EQ.64)) ME2GW(JME)=JGW
                     ENDIF
                  ENDIF
               ENDIF
               IF((I.GT.1).AND.(I.LT.NVELL(K))) THEN !IF NOT FIRST OR
                  JME=JME+1     !LAST NODE
                  ME2GW(JME)=JGW !M.E. USES IT
               ENDIF
               IF(I.EQ.NVELL(K)) THEN !DEAL WITH LAST NODE ON BOUNDARY
                  IF((NBV(JGW).NE.NBVV(K,1)).AND. !IF UNCLOSED BOUNDARY
     &                 (NBV(JGW).NE.NBV(1))) THEN !M.E. USES IT
                     JME=JME+1
                     ME2GW(JME)=JGW
                  ENDIF
                  IF(NBVV(K,I).EQ.NBV(1)) THEN !IF OVERLAPS WITH VERY FIRST
                     IF((LBCODEI(JGW).EQ.2) .OR. ! L.B. NODE
     &                    (LBCODEI(JGW).EQ.12).OR. ! M.E. USES IT ONLY IF IT IS
     &                    (LBCODEI(JGW).EQ.22).OR. ! SPECIFIED FLOW,
     &                    (LBCODEI(JGW).EQ.3) .OR. ! AN OVERFLOW BARRIER OR
     &                    (LBCODEI(JGW).EQ.13).OR. ! A RADIATION
     &                    (LBCODEI(JGW).EQ.23).OR. ! BOUNDARY
     &                    (LBCODEI(JGW).EQ.30).OR.
     &                    (LBCODEI(JGW).EQ.32).OR.
     &                    (LBCODEI(JGW).EQ.40).OR.
     &                    (LBCODEI(JGW).EQ.41).OR.
     &                    (LBCODEI(JGW).EQ.52).OR.
     &                    (LBCODEI(JGW).EQ.64)) ME2GW(1)=JGW
                  ENDIF
               ENDIF

!...........   LOAD EXTERNAL BARRIER BOUNDARY INFORMATION INTO THE CORRECT VECTORS
               IF((IBTYPE(K).EQ.3).OR.(IBTYPE(K).EQ.13)
     &          .OR.(IBTYPE(K).EQ.23)) THEN
                  BARLANHT(JGW)=BARLANHTR(K,I)
                  BARLANCFSP(JGW)=BARLANCFSPR(K,I)
               ENDIF

!...........   LOAD INTERNAL BARRIER BOUNDARY INFORMATION INTO THE CORRECT VECTORS
               IF((IBTYPE(K).EQ.4).OR.(IBTYPE(K).EQ.24).OR.
     &            (IBTYPE(K).EQ.64)) THEN
                  IBCONN(JGW)=IBCONNR(K,I)
                  BARINHT(JGW)=BARINHTR(K,I)
                  BARINCFSB(JGW)=BARINCFSBR(K,I)
                  BARINCFSP(JGW)=BARINCFSPR(K,I)
               ENDIF

!...........   LOAD INTERNAL BARRIER WITH PIPES BOUNDARY INFORMATION INTO
!...........   THE CORRECT VECTORS
               IF((IBTYPE(K).EQ.5).OR.(IBTYPE(K).EQ.25)) THEN
                  IBCONN(JGW)=IBCONNR(K,I)
                  BARINHT(JGW)=BARINHTR(K,I)
                  BARINCFSB(JGW)=BARINCFSBR(K,I)
                  BARINCFSP(JGW)=BARINCFSPR(K,I)
                  PIPEHT(JGW)=PIPEHTR(K,I)
                  PIPECOEF(JGW)=PIPECOEFR(K,I)
                  PIPEDIAM(JGW)=PIPEDIAMR(K,I)
               ENDIF

!...........   WRITE OUT BOUNDARY CONDITION ARRAY INFORMATION

               WRITE(16,1857) JGW,JME,ME2GW(JME),NBV(JGW),LBCODEI(JGW),
     &              THETA,CSII(JGW),SIII(JGW),BNDLEN2O3(JGW)
 1857          FORMAT(1X,I6,1X,I6,1X,I6,3X,I8,3X,I4,8X,E12.4,2X,E16.8,
     &              1X,E16.8,2X,E16.8)

!...........   CHECK EXTERNAL BARRIER HEIGHTS AGAINST DEPTHS
               IF((IBTYPE(K).EQ.3).OR.(IBTYPE(K).EQ.13)
     &          .OR.(IBTYPE(K).EQ.23)) THEN
                  IF(BARLANHT(JGW).LT.-DP(NBV(JGW))) THEN
                     write(scratchMessage,8367)
     &                    JGW,NBV(JGW),BARLANHT(JGW),DP(NBV(JGW))
                     call allMessage(ERROR,scratchMessage)
 8367               FORMAT(
     &           '!!!!!!!!!!  FATAL INPUT ERROR   !!!!!!!!!!!! ',
     &           'AT BOUNDARY NODE NO.',i0,' (GLOBAL NODE NO.',
     &           i0, ' AND OF EXTERNAL BARRIER TYPE) ',
     &           'THE EXTERNAL BARRIER HEIGHT = ',E12.5,
     &           'IS EXCEEDED BY THE DEPTH SPECIFIED AT ',
     &           'THE ASSOCIATED GLOBAL NODE = ',E12.5,
     &           'USER MUST SPECIFY CONSISTENT BARRIER HEIGHTS',
     &           ' AND DEPTHS')
                     call allMessage(ERROR, trim(exeTerm9973))
                     CALL terminate()
                  ENDIF
               ENDIF
!...........   CHECK INTERNAL BARRIER HEIGHTS AGAINST DEPTHS
               IF((IBTYPE(K).EQ.4).OR.(IBTYPE(K).EQ.24).OR.
     &            (IBTYPE(K).EQ.64)) THEN
                  IF(BARINHT(JGW).LT.-DP(NBV(JGW))) THEN
                     write(scratchMessage,8368)
     &                    JGW,NBV(JGW),BARINHT(JGW),DP(NBV(JGW))
                     call allMessage(ERROR, scratchMessage)
 8368               FORMAT(
     &           '!!!!!!!!!!  FATAL INPUT ERROR   !!!!!!!!!!!!',
     &           'AT BOUNDARY NODE NO.',i0,' (GLOBAL NODE NO. ',
     &           i0,' AND OF INTERNAL BARRIER TYPE) ',
     &           'THE INTERNAL BARRIER HEIGHT = ',E12.5,
     &           'IS EXCEEDED BY THE DEPTH SPECIFIED AT ',
     &           'THE ASSOCIATED GLOBAL NODE = ',E12.5,
     &           'USER MUST SPECIFY CONSISTENT BARRIER HEIGHTS',
     &           ' AND DEPTHS')
                     call allMessage(ERROR, trim(exeTerm9973))
                     CALL terminate()
                  ENDIF
               ENDIF

!...........   CHECK INTERNAL BARRIER WITH PIPES BARRIER HEIGHTS AGAINST DEPTHS
               IF((IBTYPE(K).EQ.5).OR.(IBTYPE(K).EQ.25)) THEN
                  IF(BARINHT(JGW).LT.-DP(NBV(JGW))) THEN
                     write(scratchMessage,8370)
     &                    JGW,NBV(JGW),BARINHT(JGW),DP(NBV(JGW))
                     call allMessage(ERROR, scratchMessage)
 8370            FORMAT(
     &           '!!!!!!!!!!  FATAL INPUT ERROR   !!!!!!!!!!!!',
     &           'AT BOUNDARY NODE NO.',i0,' (GLOBAL NODE NO. ',
     &           i0,' AND OF INTERNAL BARRIER TYPE) ',
     &           'THE INTERNAL BARRIER HEIGHT = ',E12.5,
     &           'IS EXCEEDED BY THE DEPTH SPECIFIED AT ',
     &           'THE ASSOCIATED GLOBAL NODE = ',E12.5,
     &           'USER MUST SPECIFY CONSISTENT BARRIER HEIGHTS',
     &           ' AND DEPTHS')
                     call allMessage(ERROR, trim(exeTerm9973))
                     CALL terminate()
                  ENDIF
               ENDIF
!...           CHECK INTERNAL BARRIER WITH PIPES PIPE HEIGHTS AGAINST DEPTHS
               IF((IBTYPE(K).EQ.5).OR.(IBTYPE(K).EQ.25)) THEN
                  IF(PIPEHT(JGW).LT.-DP(NBV(JGW))) THEN
                     write(scratchMessage,8372)
     &                    JGW,NBV(JGW),BARINHT(JGW),DP(NBV(JGW))
                     call allMessage(ERROR, scratchMessage)
 8372            FORMAT(
     &           '!!!!!!!!!!  FATAL INPUT ERROR   !!!!!!!!!!!!',
     &           'AT BOUNDARY NODE NO.',i0,' (GLOBAL NODE NO. ',i0,
     &           ' AND OF INTERNAL BARRIER TYPE) ',
     &           'THE BARRIER PIPE HEIGHT = ',E12.5,
     &           'IS EXCEEDED BY THE DEPTH SPECIFIED AT ',
     &           'THE ASSOCIATED GLOBAL NODE = ',E12.5,
     &           'USER MUST SPECIFY CONSISTENT PIPE HEIGHTS',
     &           ' AND DEPTHS')
                     call allMessage(ERROR, trim(exeTerm9973))
                     CALL terminate()
                  ENDIF
               ENDIF

!...  CHECK FOR OVERLAPPING OF AN INTERNAL BARRIER BOUNDARY WITH
!...  ANY EXTERNAL BARRIER BOUNDARY. IF THIS DOES OCCUR, TAKE
!...  APPROPRIATE ACTION
               IF((IBTYPE(K).EQ.4).OR.(IBTYPE(K).EQ.24).OR.
     &              (IBTYPE(K).EQ.64).OR.(IBTYPE(K).EQ.5)
     &              .OR.(IBTYPE(K).EQ.25)) THEN
                  DO ICK=1,NVELEXT
!...  CHECK IF OVERLAP EXISTS
                     IF(NBV(ICK).EQ.NBV(JGW)) THEN
!...  CHECK FOR ILLEGAL OVERLAPS
                        IF((LBCODEI(ICK).EQ.2).OR.
     &                   (LBCODEI(ICK).EQ.3).OR.(LBCODEI(ICK).EQ.12).OR.
     &                   (LBCODEI(ICK).EQ.13)) THEN
                           write(scratchMessage,8567)
     &                             JGW,NBV(JGW),ICK,NBV(ICK)
                           call allMessage(ERROR, scratchMessage)
 8567                      FORMAT(
     &               '!!!!!!!!!!  FATAL INPUT ERROR   !!!!!!!!!!!!',
     &               'BOUNDARY NODE NO. ',i0,' (GLOBAL NODE NO. ',
     &               i0, 'AND OF INTERNAL BARRIER TYPE) ',
     &               'OVERLAPS BOUNDARY NODE NO.',i0,' (GLOBAL NODE',
     &               ' NO.',i0,' )',
     &               'THIS IS AN ILLEGAL TYPE OVERLAP !! - INTERNAL ',
     &               'BARRIER BOUNDARIES CAN ONLY OVERLAP WITH ',
     &               'NO NORMAL FLOW EXTERNAL BOUNDARIES',
     &               '(I.E. IBTYPE=0,10,20)')
                           call allMessage(ERROR, trim(exeTerm9973))
                           call terminate()
                        ENDIF
!...  CHECK FOR OVERLAPS WHICH REQUIRE ADJUSTMENTS OF BOUNDARY
!...  CODE ON THE EXTERNAL BOUNDARY
                        IF(((IBTYPE(K).EQ.4).AND.(LBCODEI(ICK).EQ.0))
     &                       .OR.((IBTYPE(K).EQ.5)
     &                       .AND.(LBCODEI(ICK).EQ.0))) THEN
                           WRITE(16,8568) JGW,ICK,ICK
 8568                      FORMAT(1X,'DUE TO LEGAL OVERLAPPING OF ',
     &               'BOUNDARY NODE',i0,' (WHICH IS AN ESSENTIAL INTER'
     &               ,'NAL BARRIER BOUNDARY NODE)', /,2X,
     &               'AND BOUNDARY NODE',i0,' (WHICH IS AN ESSENTIAL ',
     &               'EXTERNAL NO NORMAL FLOW WITH SLIP BOUNDARY',
     &               ' NODE),',/,2X,
     &               'THE BOUNDARY TYPE FOR BOUNDARY NODE ',i0,
     &               ' IS BEING RESET TO IBTYPE=20',/,2X,
     &               '(NATURAL NO NORMAL FLOW WITH SLIP BOUNDARY) ')
                           LBCODEI(ICK)=20
                        ENDIF
                        IF(((IBTYPE(K).EQ.4).AND.(LBCODEI(ICK).EQ.10))
     &                       .OR.((IBTYPE(K).EQ.5)
     &                       .AND.(LBCODEI(ICK).EQ.10))) THEN
                           WRITE(16,8569) JGW,ICK,ICK
 8569                      FORMAT(1X,'DUE TO LEGAL OVERLAPPING OF ',
     &               'BOUNDARY NODE ',i0,' (WHICH IS AN ESSENTIAL INTER'
     &               ,'NAL BARRIER BOUNDARY NODE)', /,2X,
     &               'AND BOUNDARY NODE',i0,' (WHICH IS AN ESSENTIAL ',
     &               'EXTERNAL NO NORMAL FLOW WITH NO SLIP BOUNDARY',
     &               ' NODE),',/,2X,
     &               'THE BOUNDARY TYPE FOR BOUNDARY NODE ',i0,
     &               ' IS BEING RESET TO IBTYPE=20',/,2X,
     &               '(NATURAL NO NORMAL FLOW WITH SLIP BOUNDARY) ')
                           LBCODEI(ICK)=20
                        ENDIF
                        IF((((IBTYPE(K).EQ.24).OR.(IBTYPE(K).EQ.64))
     &                       .AND.(LBCODEI(ICK).EQ.10))
     &                       .OR.((IBTYPE(K).EQ.25)
     &                       .AND.(LBCODEI(ICK).EQ.10))) THEN
                           WRITE(16,8570) JGW,ICK,ICK
 8570                      FORMAT(1X,'DUE TO LEGAL OVERLAPPING OF ',
     &               'BOUNDARY NODE',i0,' (WHICH IS A NATURAL INTERNAL'
     &               ,' BARRIER BOUNDARY NODE)', /,2X,
     &              'AND BOUNDARY NODE',i0,' (WHICH IS AN ESSENTIAL ',
     &              'EXTERNAL NO NORMAL FLOW WITH NO SLIP BOUNDARY',
     &              ' NODE),',/,2X,
     &              'THE BOUNDARY TYPE FOR BOUNDARY NODE',i0,
     &                   ' IS BEING RESET TO IBTYPE=0',/,2X,
     &               '(ESSENTIAL NO NORMAL FLOW WITH SLIP BOUNDARY) ')
                           LBCODEI(ICK)=0
                        ENDIF
                     ENDIF
                  END DO  !ICK=1,NVELEXT

               ENDIF

            END DO NVELL_LOOP ! I = 1,NVELL(K)
         END DO  LOOPPRBI !IPRBI=1,NPRBI
      END DO NBOULOOP  ! K=1,NBOU

Csb...10/13/2022 moved from an earlier section
      IF (ALLOCATED(ISSUBMERGED64)) THEN
         ISSUBMERGED64(:) = 0
      ENDIF
Csb...

!...  ONCE ALL FLOW BOUNDARY NODES HAVE BEEN PROCESSED, CHECK TO MAKE
!...  SURE THAT JGW LE MNVEL.  NOTE, JME MUST BE < JGW.

      IF(MNVEL.LT.JGW) THEN
         write(scratchMessage,9947)
         call allMessage(ERROR, scratchMessage)
 9947    FORMAT('!!!!!!!!!!  FATAL INPUT ERROR   !!!!!!!!!!!!',
     &       'THE DIMENSION PARAMETER MNVEL IS LESS THAN ',
     &       'THE TOTAL NUMBER OF FLOW BOUNDARY NODES',
     &       'FROM ALL THE SPECIFIED FLOW SEGMENTS COMBINED')
         call allMessage(ERROR, exeTerm9973)
         CALL terminate()
      ENDIF

      NVEL=JGW
      NVELME=JME

!.....TRANSFER FLOW BOUNDARY INFORMATION INTO NODAL ARRAYS
      DO I=1,NP
         LBArray_Pointer(I)=-1
         CSI(I)=0.
         SII(I)=1.
      END DO

      !do i=1,np !jgfdebug
      !   write(6,*) 'i, lbarray_pointer(i)=',i,lbarray_pointer(i)
      !end do
      !write(6,*) 'i=1,nvelme  nvelme=',nvelme
      DO I=1,NVELME
         J=ME2GW(I)
            !write(6,*) 'i, j=me2gw(i)',i,j !jgfdebug
         LBArray_Pointer(NBV(J))=J
            !write(6,*) 'nbv(j)=',nbv(j) !jgfdebug
         CSI(NBV(J))=CSII(J)
         SII(NBV(J))=SIII(J)
      END DO
      !do i=1,np !jgfdebug
      !   write(6,*) 'i, lbarray_pointer(i)=',i,lbarray_pointer(i)
      !end do
      !do i=1,np !jgfdebug
      !   if (lbarray_pointer(i).gt.0) then
      !      write(6,*) 'i, lbcodei(lbarray_pointer(i))=',i,lbcodei(lbarray_pointer(i)) !jgfdebug
      !   endif
      !end do    
      !stop !jgfdebug

      CALL PERBCSETUP( NM(:,:) ) ;

!...  IF ANY EXTERNAL BARRIER BOUNDARIES WERE SPECIFIED, (NFLUXB=1)
!.....WRITE OUT EXTERNAL BARRIER BOUNDARY INFORMATION TO UNIT 16 FILE
!.....NOTE THAT THIS INFORMATION WAS READ IN FROM THE UNIT 14 FILE
      IF(NFLUXB.EQ.1) THEN
!.....  WRITE OUT INFO ON SPECIFIED EXTERNAL BARRIER BOUNDARIES
        WRITE(16,1112)
        WRITE(16,2220)
 2220   FORMAT(//,1X,'EXTERNAL BARRIER BOUNDARY INFORMATION ',/)

!.......OUTPUT ELEVATION OF EXTERNAL BARRIER NODES ABOVE THE GEOID AND
!........THE COEFFICIENT OF FREE SURFACE SUPERCRITICAL FLOW AT
!........DESIGNATED EXTERNAL BARRIER BOUNDARY NODES TO UNIT 16
         WRITE(16,2224)
 2224   FORMAT(//,9X,'NODE',10X,'BARRIER HEIGHT',
     &                     6X,'SUPER-CRIT. EXTERNAL BAR. COEF.',/)
         DO J=1,NVEL
            IF((LBCODEI(J).EQ.3).OR.(LBCODEI(J).EQ.13)
     &         .OR.(LBCODEI(J).EQ.23)) THEN
               WRITE(16,2225) NBV(J),BARLANHT(J),BARLANCFSP(J)
 2225         FORMAT(5X,i0,6X,F14.5,15X,F12.3)
            ENDIF
         END DO
      ENDIF

!...  IF ANY INTERNAL BARRIER BOUNDARIES WERE SPECIFIED, (NFLUXIB=1)
!.....WRITE INTERNAL BARRIER BOUNDARY INFORMATION TO UNIT 16 FILE
      IF(NFLUXIB.EQ.1) THEN

!.....   WRITE OUT INFO ON SPECIFIED INTERNAL BARRIER BOUNDARIES
         WRITE(16,1112)
         WRITE(16,2320)
 2320   FORMAT(//,1X,'INTERNAL BARRIER BOUNDARY INFORMATION ',/)

!....... WRITE CONNECTION NODE NUMBER AND ELEVATION OF THE INTERNAL BARRIER
!........NODES ABOVE THE GEOID AND THE COEFFICIENTS OF FREE SURFACE SUPERCRITICAL
!........AND SUBCRITICAL FLOW AT DESIGNATED INTERNAL BARRIER BOUNDARY NODES
!........TO UNIT 16 (NOTE THAT THIS INFORMATION WAS INPUT FROM THE UNIT 14
!........FILE WITH BOUNDARY NODE INFORMATION)
         WRITE(16,2324)
 2324   FORMAT(//,9X,'NODE',6X,'CONNECTED NODE',6X,'BARRIER HEIGHT',
     &                4X,'SUB-CRIT. INT. BAR. COEF.',
     &                4X,'SUPER-CRIT. INT. BAR. COEF.',/)
         DO J=1,NVEL
            IF((LBCODEI(J).EQ.4).OR.(LBCODEI(J).EQ.24).OR.
     &         (LBCODEI(J).EQ.64)) THEN
               WRITE(16,2325) NBV(J),IBCONN(J),BARINHT(J),
     &                       BARINCFSB(J),BARINCFSP(J)
 2325         FORMAT(5X,i0,7X,i0,6X,F14.5,12X,F12.3,17X,F12.3)
            ENDIF
         END DO
      ENDIF

!jjwm001 - begin add
!...  IF ANY INTERNAL BARRIER BOUNDARIES WITH CROSS BARRIER PIPES
!.....WERE SPECIFIED, (NFLUXIBP=1)
!.....WRITE INTERNAL BARRIER BOUNDARY INFORMATION WITH CROSS
!.....BARRIER PIPE INFORMATION TO UNIT 16 FILE
      IF(NFLUXIBP.EQ.1) THEN
!.....   WRITE OUT INFO ON SPECIFIED INTERNAL BARRIER BOUNDARIES
         WRITE(16,1112)
         WRITE(16,2326)
 2326   FORMAT(//,1X,'INTERNAL BARRIER BOUNDARY WITH CROSS BARRIER',
     &   ' PIPE INFORMATION ',/)

!....... WRITE CONNECTION NODE NUMBER AND ELEVATION OF THE INTERNAL BARRIER
!........NODES ABOVE THE GEOID AND THE COEFFICIENTS OF FREE SURFACE SUPERCRITICAL
!........AND SUBCRITICAL FLOW AT DESIGNATED INTERNAL BARRIER BOUNDARY NODES
!........IN ADDITION TO CROSS BARRIER PIPE CROWN HEIGHT, CROSS BARRIER PIPE
!........COEFFICIENT AND CROSS BARRIER PIPE DIAMETER TO UNIT 16
!........(NOTE THAT THIS INFORMATION WAS INPUT FROM THE UNIT 14 FILE WITH
!........BOUNDARY NODE INFORMATION)
         WRITE(16,2327)
 2327    FORMAT(//,7X,'NODE',4X,'CONNECTED NODE',4X,'BARRIER HEIGHT',
     &                4X,'SUB-CRIT INT BAR COEF',
     &                4X,'SUPER-CRIT INT BAR COEF',
     &                4X,'PIPEHT  ',
     &                4X,'PIPECOEF',
     &                4X,'PIPEDIAM',/)
         DO J=1,NVEL
            IF((LBCODEI(J).EQ.5).OR.(LBCODEI(J).EQ.25)) THEN
               WRITE(16,2328) NBV(J),IBCONN(J),BARINHT(J),
     &                     BARINCFSB(J),BARINCFSP(J),
     &                     PIPEHT(J),PIPECOEF(J),PIPEDIAM(J)
 2328         FORMAT(3X,i0,5X,i0,4X,F14.5,8X,F12.3,12X,F12.3,
     &              2X,F10.5,2X,F10.5,2X,F10.5)
            ENDIF
         END DO
      ENDIF

!sb...
#ifdef CMPI
      nfluxib64_gbl = msg_imax(nfluxib64)
#else
      nfluxib64_gbl = nfluxib64
#endif
!sb...

!jjwm001 - end add
      ! jgf51.11.18: Free memory associated with the reading of the
      ! levee and culvert parameters since it is no longer needed
      call freeFluxBoundaryArrayTemporaries()
!
!...GENERAL PURPOSE FORMAT STATEMENTS for subtly expressed error messages
!...
 1112 FORMAT(/,1X,79('_'))
 9973 FORMAT('!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!')

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
!     ------------------------------------------------------------------
      end subroutine initializeBoundaries
!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
!                          S U B R O U T I N E
!     A L L O C A T E   N O D A L   A N D   E L E M E N T A L   A R R A Y S
!     ------------------------------------------------------------------
!     jgf51.21.09: Mesh related memory allocation for any array that is
!     dimensioned by the number of nodes in the mesh or the number of
!     elements in the mesh.
!     ------------------------------------------------------------------
      subroutine allocateNodalAndElementalArrays()
      implicit none

      call setMessageSource("allocateNodalAndElementalArrays")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      allocate ( slam(mnp),sfea(mnp),x(mnp),y(mnp),labels(mnp))
      call dict(node_dict, mnp)
      allocate ( dp(mnp) )
      allocate ( dp_orig(mnp) )
      allocate ( nm(mne,3))
      allocate ( areas(mne), totalarea(mnp))
      allocate ( nneigh(mnp),mju(mnp),nodele(mnp))
      allocate ( lbarray_pointer(mnp))

      allocate ( csi(mnp),sii(mnp))
!...  v49.48.02 tcm -- Variables related to kdtree searches
!...  which are deallocated at the end of read_input.F
      allocate( bcxy(2,mne), rmax(mne) )

      ALLOCATE ( FDXE(3,MNE), FDYE(3,MNE) ) ;
      FDXE = 0.0D0 ; 
      FDYE = 0.0D0 ;
      
      ALLOCATE( EXVO(3,MNE), EYVO(3,MNE) ) ; 
      EXVO = 0.0D0 ;
      EYVO = 0.0D0 ; 
      EXV => EXVO ; 
      EYV => EYVO ; 

      NULLIFY( EXVR ) ;
      NULLIFY( EYVR ) ; 
        
      ALLOCATE( SFAC(MNP), SFMX(MNP), SFMY(MNP), 
     &          SFCX(MNP), SFCY(MNP), YCSFAC(MNP),
     &          SFCT(MNP), TANPHI(MNP) )

      SFCT = 1.0D0 ; 
      SFCX = 1.0D0 ; 
      SFCY = 1.0D0 ;
        
      SFMX = 1.0D0 ;
      SFMY = 1.0D0 ; 

      TANPHI = 0.0D0 ;
      YCSFAC = 1.0D0 ;
        
      ALLOCATE( SfacEle(MNE), SFCYEle(MNE), SFCXEle(MNE),
     &          SFMXEle(MNE), SFMYEle(MNE), YCSFacEle(MNE),
     &          TANPHIEle(MNE) ) 
      SFacEle = 1.0D0 ; 
      SFCYEle = 1.0D0 ;
      SFCXEle = 1.0D0 ; 
      SFMYEle = 1.0D0 ;
      SFMXEle = 1.0D0 ;
      YCSFacEle = 1.0D0 ;
      TANPHIEle = 0.0D0 ;
      !
      ! initialize to something troublesome to make it easy to spot issues
      slam = -99999.d0
      sfea = -99999.d0
      x = -99999.d0
      y = -99999.d0
      dp = -99999.d0
      dp_orig = -99999.d0
      nm = -99999
      areas = -99999.d0
      totalarea = -99999.d0
      nneigh = -99999
      mju = -99999
      nodele = -99999
      lbarray_pointer = -99999
      csi = -99999.d0
      sii = -99999.d0
      bcxy = -99999.d0
      rmax = -99999.d0

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
!     ------------------------------------------------------------------
      end subroutine allocateNodalAndElementalArrays
!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
!     S U B R O U T I N E   A L L O C A T E   N E I G H B O R   A R R A Y S
!     ------------------------------------------------------------------
!     Allocate space for Arrays needed to determine neighbor tables
!     ------------------------------------------------------------------
      subroutine allocateNeighborArrays()
      implicit none

      call setMessageSource("allocateNeighborArrays")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      allocate ( neitab(mnp,mnei))
      allocate ( neitabele(mnp,mnei))

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
!     ------------------------------------------------------------------
      end subroutine allocateNeighborArrays
!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
!     S U B R O U T I N E   F R E E  M E S H
!     ------------------------------------------------------------------
!     De allocate memory that is no longer needed after the input
!     has been read.
!     ------------------------------------------------------------------
      subroutine freeMesh()
      implicit none

      call setMessageSource("freeMesh")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      deallocate(rmax)
      deallocate(bcxy)

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
!     ------------------------------------------------------------------
      end subroutine freeMesh
!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
!               S U B R O U T I N E    N E I G H B
!     ------------------------------------------------------------------
!      Subroutine to generate neighbor tables from a connectivity table.
!
!      NOTES
!      a node neighbor table is generated with the node itself is listed as
!         neighbor #1 and all other neighbors are sorted and placed in cw
!         order from east
!      a neighbor element table is generated with:
!         entry 1 = element # defined by neighbors 1,2,3
!         entry 2 = element # defined by neighbors 1,3,4
!         entry 3 = element # defined by neighbors 1,4,5
!          .......
!         entry last = element # defined by neighbors 1,nneigh,2
!         a zero area means that the defined triangle lies outside the domain
!
!
!    v1.0   R.L.   6/29/99  used in 3D code
!    v2.0   R.L.   5/23/02  adapted to provide neighbor el table
!
!     -  PARAMETERS WHICH MUST BE SET TO CONTROL THE DIMENSIONING OF ARRAYS
!           ARE AS FOLLOWS:
!
!          MNP = MAXIMUM NUMBER OF NODAL POINTS
!          MNE = MAXIMUM NUMBER OF ELEMENTS
!          MNEI= 1+MAXIMUM NUMBER OF NODES CONNECTED TO ANY ONE NODE IN THE
!                   FINITE ELEMENT GRID
!
!    VARIABLE DEFINITIONS:
!       NE - NUMBER OF ELEMENTS
!       NP - NUMBER OF NODES
!       NM(MNE,3) - NODE NUMBERS ASSOCIATED WITH EACH ELEMENT
!       NNeigh(MNP) NUMBER OF NEIGHBORS FOR EACH NODE
!       NeiTab(MNP,NEIMAX) 2D ARRAY OF NEIGHBORS FOR EACH NODE
!       NeiTabEle(MNP,NEIMAX) 2D ARRAY OF NEIGHBOR ELEMENTS FOR EACH NODE
!       NEIMIN - 1+MINIMUM NUMBER OF NEIGHBORS FOR ANY NODE
!       NEIMAX - 1+MAXIMUM NUMBER OF NEIGHBORS FOR ANY NODE
!     ------------------------------------------------------------------
      SUBROUTINE NEIGHB()
      IMPLICIT NONE
      INTEGER :: N,NN,I,J,JJ,K,JLOW
      INTEGER :: NN1,NN2,NN3,NE1,NE2,NE3
      REAL(8) :: DELX,DELY,DIST
      REAL(8) :: ANGLELOW,ANGLEMORE
      REAL(8), ALLOCATABLE :: ANGLE(:)
      INTEGER, ALLOCATABLE :: NEITEM(:)
      LOGICAL :: FOUND

      INTEGER:: EJ, IJ, II
      INTEGER, ALLOCATABLE::  NeiTabETMP(:,:)

      call setMessageSource("neighb")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      ALLOCATE ( ANGLE(MNEI) )
      ALLOCATE ( NEITEM(MNP) )
      IF ( .NOT. ALLOCATED(NNeighEle) ) THEN
          ALLOCATE ( NNeighEle(MNP) )
      END IF

      ALLOCATE ( NeiTabETMP(MNP,MNEI) ) ;

      DO N=1,NP
         NNeigh(N)=0
         NNeighEle(N)=0
         DO NN=1,MNEI
            NeiTab(N,NN)=0
            NeiTabEle(N,NN)=0
         END DO
      END DO

      DO N=1,NE
         NN1 = NM(N,1)
         NN2 = NM(N,2)
         NN3 = NM(N,3)

         NNeighEle(NN1)=NNeighEle(NN1)+1
         NNeighEle(NN2)=NNeighEle(NN2)+1
         NNeighEle(NN3)=NNeighEle(NN3)+1
         NeiTabEle(NN1,NNeighEle(NN1))=N
         NeiTabEle(NN2,NNeighEle(NN2))=N
         NeiTabEle(NN3,NNeighEle(NN3))=N

         FOUND = .FALSE.
         DO J=1,NNeigh(NN1)
            IF(NN2.EQ.NeiTab(NN1,J)) THEN
               FOUND = .TRUE.
               EXIT
            ENDIF
         END DO
         IF (.NOT.FOUND) THEN
            NNeigh(NN1)=NNeigh(NN1)+1
            NNeigh(NN2)=NNeigh(NN2)+1
            IF((NNeigh(NN1).GT.MNEI-1).OR.(NNeigh(NN2).GT.MNEI-1))THEN
                call neighb_terminateDimensioning()
            ENDIF
            NeiTab(NN1,NNeigh(NN1))=NN2
            NeiTab(NN2,NNeigh(NN2))=NN1

            NeiTabETMP(NN1,NNeigh(NN1))=N ;
            NeiTabETMP(NN2,NNeigh(NN2))=N ;
         ENDIF

         FOUND = .FALSE.
         DO J=1,NNeigh(NN1)
            IF(NN3.EQ.NeiTab(NN1,J)) THEN
               FOUND = .TRUE.
               EXIT
            ENDIF
         END DO
         IF (.NOT.FOUND) THEN
            NNeigh(NN1)=NNeigh(NN1)+1
            NNeigh(NN3)=NNeigh(NN3)+1
            IF((NNeigh(NN1).GT.MNEI-1).OR.(NNeigh(NN3).GT.MNEI-1))THEN
                call neighb_terminateDimensioning()
            ENDIF
            NeiTab(NN1,NNeigh(NN1))=NN3
            NeiTab(NN3,NNeigh(NN3))=NN1

            NeiTabETMP(NN1,NNeigh(NN1))=N ;
            NeiTabETMP(NN3,NNeigh(NN3))=N ;
         ENDIF

         FOUND = .FALSE.
         DO J=1,NNeigh(NN2)
            IF(NN3.EQ.NeiTab(NN2,J)) THEN
               FOUND = .TRUE.
               EXIT
            ENDIF
         END DO
         IF (.NOT.FOUND) THEN
            NNeigh(NN2)=NNeigh(NN2)+1
            NNeigh(NN3)=NNeigh(NN3)+1
            IF((NNeigh(NN2).GT.MNEI-1).OR.(NNeigh(NN3).GT.MNEI-1))THEN
                call neighb_terminateDimensioning()
            ENDIF
            NeiTab(NN2,NNeigh(NN2))=NN3
            NeiTab(NN3,NNeigh(NN3))=NN2

            NeiTabETMP(NN2,NNeigh(NN2))=N ;
            NeiTabETMP(NN3,NNeigh(NN3))=N ;
         ENDIF
      END DO

!
!     INSERT NODE ITSELF IN PLACE #1 and SORT other NEIGHBORS by
!     increasing cw angle from East
!
      DO I=1,NP
         DO J=1,NNeigh(I)
            NEITEM(J)=NeiTab(I,J)
            DELX=X(NEITEM(J))-X(I)
            DELY=Y(NEITEM(J))-Y(I)

            EJ = NeiTabETMP(I,J) ;
            IJ = sum(MERGE( (/ 1, 2, 3 /), (/ 0, 0, 0 /),
     &             NM(EJ,:) - NEITEM(J) == 0 )) ;

            II = sum(MERGE( (/ 1, 2, 3 /), (/ 0, 0, 0 /),
     &             NM(EJ,:) - I == 0 )) ;

            DELX=EXV(IJ,EJ) - EXV(II,EJ) 
            DELY=EYV(IJ,EJ) - EYV(II,EJ) 

            DIST=SQRT(DELX*DELX+DELY*DELY)


            IF(DIST.EQ.0.0d0) THEN
                call neighb_terminateDuplicateNode(I,NEITEM(J))
            ENDIF
            IF(DELY.NE.0.0d0) THEN
               ANGLE(J)=RAD2DEG*ACOS(DELX/DIST)
               IF(DELY.GT.0.0) ANGLE(J)=360.0d0-ANGLE(J)
            ENDIF
            IF(DELY.EQ.0.0d0) THEN
               IF(DELX.GT.0.0d0) ANGLE(J)=0.0d0
               IF(DELX.LT.0.d0) ANGLE(J)=180.0d0
            ENDIF
         END DO
         ANGLEMORE=-1.d0
         DO JJ=1,NNeigh(I)
            ANGLELOW=400.d0
            DO J=1,NNeigh(I)
               IF((ANGLE(J).LT.ANGLELOW).AND.
     &              (ANGLE(J).GT.ANGLEMORE)) THEN
                  ANGLELOW=ANGLE(J)
                  JLOW=J
               ENDIF
            END DO
            NeiTab(I,JJ+1)=NEITEM(JLOW)
            ANGLEMORE=ANGLELOW
         END DO
         NeiTab(I,1)=I
         NNeigh(I)=NNeigh(I)+1
      ENDDO
!
!     MATCH EACH SET OF 3 NODES WITH CORRESPONDING ELEMENT AND REORDER
!     ELEMENTS ACCORDINGLY
!
      DO I=1,NP
         DO K=1,NNeighEle(I)
            NEITEM(K)=NeiTabEle(I,K)
            NeiTabEle(I,K)=0
         END DO
         NN = 0; 
         DO J=2,NNeigh(I)
            NN1=NeiTab(I,1)
            NN3=NeiTab(I,J)
            IF(J.NE.NNeigh(I)) NN2=NeiTab(I,J+1)
            IF(J.EQ.NNeigh(I)) NN2=NeiTab(I,2)
            KLP1: DO K=1,NNeighEle(I)
               IF(NEITEM(K).NE.0) THEN
                  IF(NM(NEITEM(K),1).EQ.NN1) THEN
                     NE1=NM(NEITEM(K),1)
                     NE2=NM(NEITEM(K),2)
                     NE3=NM(NEITEM(K),3)
                  ENDIF
                  IF(NM(NEITEM(K),2).EQ.NN1) THEN
                     NE1=NM(NEITEM(K),2)
                     NE2=NM(NEITEM(K),3)
                     NE3=NM(NEITEM(K),1)
                  ENDIF
                  IF(NM(NEITEM(K),3).EQ.NN1) THEN
                     NE1=NM(NEITEM(K),3)
                     NE2=NM(NEITEM(K),1)
                     NE3=NM(NEITEM(K),2)
                  ENDIF
                  IF((NE2.EQ.NN2).AND.(NE3.EQ.NN3)) THEN
                     NN = NN + 1;
                     NeiTabEle(I,NN)=NEITEM(K)
                     NEITEM(K)=0
                     CYCLE KLP1;
                  ENDIF
               ENDIF
            END DO KLP1
         END DO
      END DO

!
!  DETERMINE THE MAXIMUM AND MINIMUM NUMBER OF NEIGHBORS
!
      NEIMAX = 0
      NEIMIN = 1000
      DO N=1,NP
         IF(NNeigh(N).LT.NEIMIN) NEIMIN=NNeigh(N)
         IF(NNeigh(N).GT.NEIMAX) NEIMAX=NNeigh(N)
      END DO
!
!  WRITE OUT DIAGNOSTIC OUTPUT
!
!     OPEN(333,file='fort.333')
!     DO N=1,NP
!       WRITE(333,331) (NEIGH(N,J),J=1,NNEIGH(N))
!       WRITE(333,331) N,(NEIGHELE(N,J),J=1,NNEIGH(N)-1)
!       WRITE(333,*) ' '
!331    FORMAT(15(1X,I7))
!       END DO
!     CLOSE (333)

!     Deallocate local work arrays
      DEALLOCATE ( ANGLE )
      DEALLOCATE ( NEITEM )
!      DEALLOCATE ( NNEIGHELE )
      DEALLOCATE ( NeiTabETMP ) 

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN

!     ------------------------------------------------------------------
      end subroutine neighb
!     ------------------------------------------------------------------

      subroutine neighb_terminateDimensioning()
          implicit none
        ! TERMINATE PROGRAM IF MAXIMUM NUMBER OF NEIGHBORS SET TOO SMALL
          write(scratchMessage,99311)
          call allMessage(ERROR,scratchMessage)
99311     FORMAT('!!!!!!!!!!  FATAL ERROR !!!!!!!!! ',
     &        'THE DIMENSIONING PARAMETER MNEI IS TOO SMALL. ',
     &        'THERE IS A PROBLEM WITH THE DYNAMIC MEMORY ALLOCATION. ',
     &        '!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!')
        call terminate()
      end subroutine neighb_terminateDimensioning
 
      subroutine neighb_terminateDuplicateNode(i,j)
        implicit none
        integer,intent(in) :: i, j
        write(scratchMessage,99312) I,J
        call allMessage(ERROR,scratchMessage)
99312   FORMAT('!!!!!!!!!!  FATAL ERROR !!!!!!!!! ',
     &       'NODES ',i0,' AND ',i0,' HAVE THE SAME COORDINATES.'
     &       '!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!')
        call terminate()
      end subroutine neighb_terminateDuplicateNode


!     ------------------------------------------------------------------
!           S U B R O U T I N E   T O T A L   A R E A   C A L C
!     ------------------------------------------------------------------
!     jgf51.21.10: Refactored and modularized out of the ADCIRC_Init
!     subroutine in the ADCIRC module.
!
!     Determine the number of active elements (MJU), the total number of
!     elements (NODELE) and the total area of elements (TotalArea)
!     attached to each node.
!
!     The calculation of the number of active elements must wait until
!     the model has initialized to the point that either cstart() or
!     hstart() has already been called and nnodecode is already
!     initialized.
!     ------------------------------------------------------------------
      subroutine totalAreaCalc()
      use global, only : nodecode, nnodecode, noff
      implicit none
      real(8) :: areaEle
      integer nm1, nm2, nm3  ! node numbers around an element
      integer ncele          ! wet/dry indicator (zero if dry)
      integer :: i           ! node loop counter
      integer :: ie          ! element loop counter

      call setMessageSource("totalAreaCalc")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      !
      ! initialize values
      mju(:)=0
      nodele(:)=0
      totalarea(:)=0.d0
      nodecode(:)=nnodecode(:)
      !
      ! assemble by elements
      do ie=1,ne
         nm1=nm(ie,1)
         nm2=nm(ie,2)
         nm3=nm(ie,3)
         ncele=nodecode(nm1)*nodecode(nm2)*nodecode(nm3)*noff(ie)
         mju(nm1)=mju(nm1)+ncele
         mju(nm2)=mju(nm2)+ncele
         mju(nm3)=mju(nm3)+ncele
         nodele(nm1)=nodele(nm1)+1
         nodele(nm2)=nodele(nm2)+1
         nodele(nm3)=nodele(nm3)+1
         areaele=ncele*areas(ie)/2.d0
         totalarea(nm1)=totalarea(nm1)+areaele
         totalarea(nm2)=totalarea(nm2)+areaele
         totalarea(nm3)=totalarea(nm3)+areaele
      end do

      do i=1,np
         if(mju(i).eq.0) mju(i)=1
      end do

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
!     ------------------------------------------------------------------
      end subroutine totalAreaCalc
!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
!                  S U B R O U T I N E   I N V C P
!     ------------------------------------------------------------------
!     Transform from CPP coordinates to lon,lat (lamda,phi) coordinates
!     Lon,Lat is in radians.
!     ------------------------------------------------------------------
      SUBROUTINE INVCP(XXCP,YYCP,RLAMBDA,PHI,RLAMBDA0,PHI0)
      USE CONSTANTS, ONLY : Rearth
      IMPLICIT NONE
      REAL*8 XXCP,YYCP,RLAMBDA,PHI,RLAMBDA0,PHI0

      call setMessageSource("invcp")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      RLAMBDA=RLAMBDA0+XXCP/(Rearth*COS(PHI0))
      PHI=YYCP/Rearth

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
!     ------------------------------------------------------------------
      END SUBROUTINE INVCP
!     ------------------------------------------------------------------

      !----------------------------------------------------------------------
      !                   S U B R O U T I N E
      ! G E T   A T T R I B U T E   C H A R A C T E R I S T I C S   X D M F
      !----------------------------------------------------------------------
      ! Get info about an attribute.
      !----------------------------------------------------------------------
      subroutine getAttributeCharacteristicsXDMF(xdmfFortranObj, attributeIndex)
#ifdef ADCXDMF
      implicit none
      include 'adcirc_Xdmf.f'
#endif
      integer*8, intent(in) :: xdmfFortranObj
      integer, intent(in) :: attributeIndex
#ifdef ADCXDMF
      !
      integer :: typeHolder
      character(len=256) :: logString
      !
      call setMessageSource("getAttributeCharacteristicsXDMF")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      call xdmfRetrieveAttributeType(xdmfFortranObj, attributeIndex, typeHolder)
      call createAttributeTypeString(typeHolder, logString)
      scratchFormat = '("The Attribute type is ",a,".")'
      write(scratchMessage,scratchFormat) trim(logString)
      !
      call xdmfRetrieveAttributeCenter(xdmfFortranObj, attributeIndex, typeHolder)
      call createAttributeCenterString(typeHolder, logString)
      scratchFormat = '("The Attribute Center is ",a,".")'
      write(scratchMessage,scratchFormat) trim(logString)
      call logMessage(DEBUG,scratchMessage)
      !
      call xdmfRetrieveAttributeValueType(xdmfFortranObj, attributeIndex, typeHolder)
      call createDataTypeString(typeHolder, logString)
      scratchFormat = '("The Attribute data type is ",a,".")'
      write(scratchMessage,scratchFormat) trim(logString)
      call logMessage(DEBUG,scratchMessage)
      !
      call xdmfRetrieveAttributeSize(xdmfFortranObj, attributeIndex, typeHolder)
      scratchFormat = '("The Attribute consists of ",i0," values.")'
      write(scratchMessage,scratchFormat) typeHolder
      call logMessage(DEBUG,scratchMessage)
#endif
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      !----------------------------------------------------------------------
      end subroutine getAttributeCharacteristicsXDMF
      !----------------------------------------------------------------------


      !-----------------------------------------------------------------
      !    S U B R O U T I N E   G E T   M E S H   F O R   P R E P
      !-----------------------------------------------------------------
      ! jgf51.21.28: Called by adcprep, which uses this module to read
      ! mesh files in XDMF format. Because of the duplication between
      ! the pre_global module and the variables used here, this subroutine
      ! is necessary to extract parameters for adcprep. Someday adcprep
      ! will be integrated with the rest of the code and this sub will
      ! no longer be necessary.
      !-----------------------------------------------------------------
      subroutine getMeshForPrep(prepAgrid, prepX, prepY, prepDP)
      implicit none
      character(len=80), intent(out) :: prepAgrid
      real(8), intent(out) :: prepX(:)
      real(8), intent(out) :: prepY(:)
      real(8), intent(out) :: prepDP(:)

      call setMessageSource("getMeshForPrep")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      prepAgrid = agrid
      prepX = X
      prepY = Y
      prepDP = dp


#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      !-----------------------------------------------------------------
      end subroutine getMeshForPrep
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      !                   S U B R O U T I N E
      !        G E T   M E S H  S I Z E S  F O R   P R E P
      !-----------------------------------------------------------------
      ! jgf51.21.28: Called by adcprep, which uses this module to get
      ! mesh sizes in XDMF format.
      !-----------------------------------------------------------------
      subroutine getMeshSizesForPrep(prepNP, prepNE)
      implicit none
      integer, intent(out) :: prepNP
      integer, intent(out) :: prepNE

      call setMessageSource("getMeshSizesForPrep")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      prepNP = np
      prepNE = ne

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      !-----------------------------------------------------------------
      end subroutine getMeshSizesForPrep
      !-----------------------------------------------------------------


      !----------------------------------------------------------------------
      !  S U B R O U T I N E   R E P L A C E   N U L L S   W I T H   S P A C E S
      !----------------------------------------------------------------------
      ! The XDMF library is written in C++, which conventionally uses null
      ! characters to terminate strings. As a result, when strings come back
      ! through the XDMF library to Fortran, the strings are padded out to
      ! their full length with null characters.
      !
      ! However, Fortran generally expects the unused portion of the string to
      ! contain spaces, which allows functions like  trim(), len_trim(), and
      ! adjustl() to work properly. As a result, this subroutine is provided
      ! to convert the null characters in a string (from XDMF) to spaces for
      ! conventional use in Fortran.
      !----------------------------------------------------------------------
      subroutine replaceNullsWithSpaces(myString)
      implicit none
      integer :: nullCharLocation
      character(len=*), intent(inout) :: myString

      call setMessageSource("replaceNullsWithSpaces")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      do
         nullCharLocation = index(myString,char(0))
         if (nullCharLocation.ne.0) then
            myString(nullCharLocation:nullCharLocation) = ' '
         else
            exit ! there are no more null characters
         endif
      end do

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      !----------------------------------------------------------------------
      end subroutine replaceNullsWithSpaces
      !----------------------------------------------------------------------


      !----------------------------------------------------------------------
      !                      S U B R O U T I N E
      !     C R E A T E   A T T R I B U T E   T Y P E   S T R I N G
      !----------------------------------------------------------------------
      ! Sets the string that corresponds to the attribute type parameter from Xdmf.f
      !----------------------------------------------------------------------
      subroutine createAttributeTypeString(typeHolder, typeString)
      implicit none
      integer, intent(in) :: typeHolder
      character(len=256), intent(out) :: typeString

      call setMessageSource("createAttributeTypeString")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      !
      select case(typeHolder)
      case(200)
         typeString = 'XDMF_ATTRIBUTE_TYPE_SCALAR'
      case(201)
         typeString = 'XDMF_ATTRIBUTE_TYPE_VECTOR'
      case(202)
         typeString = 'XDMF_ATTRIBUTE_TYPE_TENSOR'
      case(203)
         typeString = 'XDMF_ATTRIBUTE_TYPE_MATRIX'
      case(204)
         typeString = 'XDMF_ATTRIBUTE_TYPE_TENSOR6'
      case(205)
         typeString = 'XDMF_ATTRIBUTE_TYPE_GLOBALID'
      case(206)
         typeString = 'XDMF_ATTRIBUTE_TYPE_NOTYPE'
      case default
         scratchFormat = '("Unrecognized attribute type ",i0,".")'
         write(scratchMessage,scratchFormat) trim(typeString)
         call allMessage(WARNING,scratchMessage)
      end select

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      !----------------------------------------------------------------------
      end subroutine createAttributeTypeString
      !----------------------------------------------------------------------


      !----------------------------------------------------------------------
      !                      S U B R O U T I N E
      !          C R E A T E   S E T   T Y P E   S T R I N G
      !----------------------------------------------------------------------
      ! Sets the string that corresponds to the set type parameter from Xdmf.f
      !----------------------------------------------------------------------
      subroutine createSetTypeString(typeHolder, typeString)
      implicit none
      integer, intent(in) :: typeHolder
      character(len=256), intent(out) :: typeString
      !
      call setMessageSource("createSetTypeString")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      select case(typeHolder)
      case(601)
         typeString = 'XDMF_SET_TYPE_NODE'
      case(602)
         typeString = 'XDMF_SET_TYPE_CELL'
      case(603)
         typeString = 'XDMF_SET_TYPE_FACE'
      case(604)
         typeString = 'XDMF_SET_TYPE_EDGE'
      case default
         scratchFormat = '("Unrecognized set type ",i0,".")'
         write(scratchMessage,scratchFormat) typeHolder
         call allMessage(WARNING,scratchMessage)
      end select

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      !----------------------------------------------------------------------
      end subroutine createSetTypeString
      !----------------------------------------------------------------------


      !----------------------------------------------------------------------
      !                      S U B R O U T I N E
      !     C R E A T E   A T T R I B U T E   C E N T E R  S T R I N G
      !----------------------------------------------------------------------
      ! Sets the string that corresponds to the attribute center parameter from Xdmf.f
      !----------------------------------------------------------------------
      subroutine createAttributeCenterString(typeHolder, typeString)
      implicit none
      integer, intent(in) :: typeHolder
      character(len=256), intent(out) :: typeString

      call setMessageSource("createAttributeCenterString")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      !
      select case(typeHolder)
      case(100)
         typeString = 'XDMF_ATTRIBUTE_CENTER_GRID'
      case(101)
         typeString = 'XDMF_ATTRIBUTE_CENTER_CELL'
      case(102)
         typeString = 'XDMF_ATTRIBUTE_CENTER_FACE'
      case(103)
         typeString = 'XDMF_ATTRIBUTE_CENTER_EDGE'
      case(104)
         typeString = 'XDMF_ATTRIBUTE_CENTER_NODE'
      case default
         scratchFormat = '("Unrecognized attribute center ",i0,".")'
         write(scratchMessage,scratchFormat) trim(typeString)
         call allMessage(WARNING,scratchMessage)
      end select

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      !----------------------------------------------------------------------
      end subroutine createAttributeCenterString
      !----------------------------------------------------------------------

      !----------------------------------------------------------------------
      !  S U B R O U T I N E   C R E A T E   D A T A   T Y P E   S T R I N G
      !----------------------------------------------------------------------
      ! Sets the string that corresponds to the data type parameter from Xdmf.f.
      !----------------------------------------------------------------------
      subroutine createDataTypeString(typeHolder, typeString)
      implicit none
      integer, intent(in) :: typeHolder
      character(len=256), intent(out) :: typeString

      call setMessageSource("createDataTypeString")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      !
      select case(typeHolder)
      case(0)
         typeString = 'XDMF_ARRAY_TYPE_INT8'
      case(1)
         typeString = 'XDMF_ARRAY_TYPE_INT16'
      case(2)
         typeString = 'XDMF_ARRAY_TYPE_INT32'
      case(3)
         typeString = 'XDMF_ARRAY_TYPE_INT64'
      case(4)
         typeString = 'XDMF_ARRAY_TYPE_UINT8'
      case(5)
         typeString = 'XDMF_ARRAY_TYPE_UINT16'
      case(6)
         typeString = 'XDMF_ARRAY_TYPE_UINT32'
      case(7)
         typeString = 'XDMF_ARRAY_TYPE_FLOAT32'
      case(8)
         typeString = 'XDMF_ARRAY_TYPE_FLOAT64'
      case default
         scratchFormat = '("Unrecognized data type ",i0,".")'
         write(scratchMessage,scratchFormat) trim(typeString)
         call allMessage(WARNING,scratchMessage)
      end select

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      !----------------------------------------------------------------------
      end subroutine createDataTypeString
      !----------------------------------------------------------------------


      !----------------------------------------------------------------------
      !                   S U B R O U T I N E
      !      C R E A T E   T O P O L O G Y   T Y P E   S T R I N G
      !----------------------------------------------------------------------
      ! Sets the string that corresponds to the data type parameter from Xdmf.f
      !----------------------------------------------------------------------
      subroutine createTopologyTypeString(typeHolder, typeString)
      implicit none
      integer, intent(in) :: typeHolder
      character(len=256), intent(out) :: typeString

      call setMessageSource("createTopologyTypeString")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      !
      select case(typeHolder)
      case(500)
         typeString = 'XDMF_TOPOLOGY_TYPE_POLYVERTEX'
      case(501)
         typeString = 'XDMF_TOPOLOGY_TYPE_POLYLINE'
      case(502)
         typeString = 'XDMF_TOPOLOGY_TYPE_POLYGON'
      case(503)
         typeString = 'XDMF_TOPOLOGY_TYPE_TRIANGLE'
      case(504)
         typeString = 'XDMF_TOPOLOGY_TYPE_QUADRILATERAL'
      case(507)
         typeString = 'XDMF_TOPOLOGY_TYPE_WEDGE'
      case(509)
         typeString = 'XDMF_TOPOLOGY_TYPE_EDGE_3'
      case(510)
         typeString = 'XDMF_TOPOLOGY_TYPE_TRIANGLE_6'
      case(511)
         typeString = 'XDMF_TOPOLOGY_TYPE_QUADRILATERAL_8'
      case(512)
         typeString = 'XDMF_TOPOLOGY_TYPE_QUADRILATERAL_9'
      case(515)
         typeString = 'XDMF_TOPOLOGY_TYPE_WEDGE_15'
      case(516)
         typeString = 'XDMF_TOPOLOGY_TYPE_WEDGE_18'
      case default
         scratchFormat = '("Unrecognized topology type ",i0,".")'
         write(scratchMessage,scratchFormat) typeHolder
         call logMessage(WARNING,scratchMessage)
      end select
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      !----------------------------------------------------------------------
      end subroutine createTopologyTypeString
      !----------------------------------------------------------------------

!...  Periodic boundary 
!...  Setting
!        (1) ID_PERBC_MAP(i) =  i, for i = interior nodes or primary nodes of a node pairs  
!                            =  p, i = slave nodes of a node pair
!        (2) PERBC_ELNDTAB(:,1:4) - table for switching vertex labels of elements 
!                                   containing slave nodes
!
!        - Switch from the original node numbering to primary nodes of the node pairs
!            ETAB( PERBC_ELNATAB(:,2), PERBC_ELNTAB(:,1) ) =   PERBC_ELNATAB(:,3) ;                      
!
!        - Switch (back) to original node numbering 
!            ETAB( PERBC_ELNATAB(:,2), PERBC_ELNTAB(:,1) ) =   PERBC_ELNATAB(:,4) ;                      
!        (3)   
!            
      SUBROUTINE PERBCSETUP( ETAB )

      use global, only: nperseg, perbc_idn_map, iperconn, nnperbc
      IMPLICIT NONE
      
!      dummy variables
      INTEGER:: ETAB(:,:)
      
!      local variables c!
      INTEGER:: I, J, K, ILOC, IE
      INTEGER:: IC, IA1, IA2, N1, N2, N3

      
      INTEGER:: ID1, ID2, ID3
      INTEGER, ALLOCATABLE:: ptab(:), stab(:), mneitab(:), mneitabele(:)
            
      IF ( NPERSEG > 0 ) THEN
         IF ( .NOT. ALLOCATED(PERBC_IDN_MAP)  ) THEN
            ALLOCATE( PERBC_IDN_MAP(MNP) ) ;
            DO I = 1,MNP
               PERBC_IDN_MAP(I) = I;
            ENDDO
         END IF
!     CALL NEIGHB(NE,NP,NM,NNeigh,NeiTab,NeiTabEle,NeiMin,NeiMax,
!     &     X,Y,NScreen,ScreenUnit)
!        - For maping original indices to primary periodic indices 
         PERBC_IDN_MAP(IPERCONN(1:NNPERBC,2)) = IPERCONN(1:NNPERBC,1) ; 

!        Create a table for mapping a original element table to periodic bc 
         NELNDTAB = SUM( NNeighEle(IPERCONN(1:NNPERBC,2)) )  ;        
         IF ( .NOT. ALLOCATED(PERBC_ELNDTAB) ) THEN
            ALLOCATE( PERBC_ELNDTAB(NELNDTAB,4) ) ;
         END IF  
         
         
         IC = 0 ; 
         DO I = 1, NNPERBC
            ID1 = IPERCONN(I,2) ;
            DO J = 1, NNeighEle(ID1)
               IE = NeiTabEle(ID1,J) ; 
               
               IC = IC + 1 ;
               PERBC_ELNDTAB(IC,1) = IE ; ! IE, VLOC, I_P, I_S 
               PERBC_ELNDTAB(IC,3:4) = IPERCONN(I,1:2) ; 
               
               KLOOP: DO K  = 1, 3
                 IF ( ETAB(IE,K) == ID1 ) THEN
                    ILOC = K ; 
                    EXIT KLOOP ; 
                 END IF
               ENDDO KLOOP

               PERBC_ELNDTAB(IC,2) = ILOC ;               
            END DO
         END DO         

         
!        Modify NeiTab, NeiTabEle 
         ALLOCATE( ptab(Mnei), stab(Mnei), 
     &             mneitab(Mnei*2), mneitabele(2*Mnei) ) ; 
         DO I = 1, NNPERBC
            ptab = 0 ; 
            stab = 0 ; 
            mneitab = 0 ;
            
            ID1 = IPERCONN(I,1) ; ! primary node !
            ID2 = IPERCONN(I,2) ; ! slave node   !
            
!     Modify the table of nodes sharing the primary nodes
            N1 = NNeigh(ID1) ; 
            N2 = NNeigh(ID2) ;
            
            ptab(1:N1) = NeiTab(ID1,1:N1) ;
            stab(1:N2) = NeiTab(ID2,1:N2) ;
             
!           1. Loop through entries in the neighboring-node table associated 
!              with nodes having the slave node as their neigbor. Replace 
!              the slave node number by the primary node number.   
            DO J = 2, N2
               ID3 = NeiTab(ID2,J) ;               
               IF ( COUNT( (IPERCONN(1:NNPERBC,2) 
     &              - ID3) == 0) == 0) THEN ! interior nodes
                  N3 = NNeigh(ID3) ; 
                  KLP2: DO K = 2, N3
                    IF ( NeiTab(ID3,K) == ID2 ) THEN
                       NeiTab(ID3,K) = ID1 ; 
                       EXIT KLP2 ; 
                    END IF
                  END DO KLP2  
               END IF               
            END DO            

!           2. Augment the neighbor tables of the primary nodes   
!           Switch to a primary node number 
            stab(1:N2) = PERBC_IDN_MAP(stab(1:N2)) ;
            
            ia1 = N1 ;
            mneitab(1:mnei) = ptab;
            DO J = 1, N2
               IF ( COUNT(mneitab(1:ia1) - stab(J) == 0) == 0 ) THEN
                  ia1 = ia1 + 1 ;
                  IF (ia1 > ubound(mneitab,1)) then
                      write(16,*) 'ia1 large ',ia1,mnei,
     &                ubound(mneitab,1),ubound(ptab,1) 
                  endif
                  mneitab(ia1) = stab(J) ;
               END IF  
            END DO
            
            IF ( ia1 <= Mnei .AND. ia1 >= N1 ) THEN
               NNeigh(ID1) = ia1 ; 
               NeiTab(ID1,1:ia1) = mneitab(1:ia1) ;                
               IF ( ia1 > NeiMax ) THEN                 
                  WRITE(16,*) "Warning in PREBCSETUP(): No. neighbors =",
     &                 IA1, " > ", " NeiMax = ", NeiMax ;    
                  NeiMax = ia1 ; 
               END IF

               IF ( ia1 == N1 ) THEN
                  WRITE(16,*) "Fatal error in PREBCSETUP(): the primary 
     &            and slave node label are identical for the node ", ID1
                  CALL EXIT(1)
               END IF
            ELSE                              
               WRITE(16,*) "Error in PREBCSETUP(): No. neighbors = ", 
     &              IA1, " > ", " Mnei = ", Mnei ;                   
               WRITE(16,*) " -- Fix: expand the size of NeiTab(:,:), 
     &                                 and NeiTabEle(:,:) tables"   ;  
               
               CALL EXPAND_NeiTAB( ia1 ) ;

               NeiMax = ia1 ; 
               NNeigh(ID1) = ia1 ; 
               NeiTab(ID1,1:ia1) = mneitab(1:ia1) ;                    
            END IF
            !
            
!     Modify the table of elements sharing the primary nodes            
            N1 = NNeighEle(ID1) ;
            N2 = NNeighEle(ID2) ; 
            
            mneitabele(1:N1) = NeiTabEle(ID1,1:N1)
            mneitabele(N1+1:N1+N2) = NeiTabEle(ID2,1:N2)    
            
            NNeighEle(ID1) = N1 + N2 ; 
            NeiTabEle(ID1,1:N1+N2) = mneitabele(1:N1+N2) ;
         END DO
         !
         DEALLOCATE( ptab, stab, mneitab, mneitabele ) ; 
      END IF
            

      RETURN ;
      CONTAINS
      
        SUBROUTINE EXPAND_NeiTAB( ia1 ) 
        IMPLICIT NONE
        
        INTEGER:: ia1 !c dummy c!
        INTEGER, ALLOCATABLE:: NeiTab_TMP(:,:), NeiTabEle_TMP(:,:) ! local variables 

        !
        allocate(NeiTab_TMP(MNP,ia1+4), NeiTabEle_TMP(MNP,ia1+4)) ;     
                       
        NeiTab_TMP(:,1:Mnei) = NeiTab(:,1:Mnei) ; 
        NeiTabEle_TMP(:,1:Mnei) = NeiTabEle(:,1:Mnei) ; 
        
        deallocate( NeiTab, NeiTabEle ) ; 

        Mnei = ia1 + 4 ;
        allocate(NeiTab(MNP,Mnei), NeiTabEle(MNP,Mnei)) ;        

        NeiTab = NeiTab_TMP ; 
        NeiTabEle = NeiTabEle_TMP ; 
        

        deallocate( NeiTab_TMP, NeiTabEle_TMP) ; 
        
        RETURN ;
        END SUBROUTINE EXPAND_NeiTAB

      END SUBROUTINE PERBCSETUP

      SUBROUTINE SWITCH_ELTAB_PERBC( IDIREC )
       use global, only: NPERSEG
       IMPLICIT NONE
       
       INTEGER:: IDIREC
       INTEGER:: I, I1, I2

       IF ( NPERSEG > 0 ) THEN
          SELECT CASE(IDIREC)
          CASE (1)
!     Switch to the element table to the primary node
             DO I = 1, NELNDTAB
                I1 = PERBC_ELNDTAB(I,1) ;
                I2 = PERBC_ELNDTAB(I,2) ; 

                NM(I1,I2) = PERBC_ELNDTAB(I,3) ;
             END DO 

          CASE (2)
!     Switch back to the original element table 
             DO I = 1, NELNDTAB
                I1 = PERBC_ELNDTAB(I,1) ;
                I2 = PERBC_ELNDTAB(I,2) ; 
                NM(I1,I2) = PERBC_ELNDTAB(I,4) ;
             END DO 
          END SELECT
       END IF


       RETURN 
      END SUBROUTINE SWITCH_ELTAB_PERBC
!-----------------------------------------------------------------------
      FUNCTION haversine(deglon1,deglon2,deglat1,deglat2) result (dist)
          ! great circle distance -- adapted from Matlab 
          real(8),intent(in) :: deglat1,deglon1,deglat2,deglon2
          real(8) :: a,c,dist,dlat,dlon,lat1,lat2
 
          dlat = deg2rad*(deglat2-deglat1)
          dlon = deg2rad*(deglon2-deglon1)
          lat1 = deg2rad*(deglat1)
          lat2 = deg2rad*(deglat2)
          a = ( sin(0.5d0*dlat) )**2 + 
     &        cos(lat1) * cos(lat2) * ( sin(0.5d0*dlon) )**2
          c = 2d0*asin( sqrt(a) )
          dist = Rearth*c
      END FUNCTION haversine
!-----------------------------------------------------------------------
      FUNCTION binarysearch(length, array, value, delta) 
          ! Given an array and a value, returns the index of the element that
          ! is closest to, but less than, the given value.
          ! Uses a binary search algorithm.
          ! "delta" is the tolerance used to determine if two values are equal
          ! if ( abs(x1 - x2) <= delta) then
          !    assume x1 = x2
          ! endif
          implicit none
          integer, intent(in) :: length
          real(8), dimension(length), intent(in) :: array
          !f2py depend(length) array
          real(8), intent(in) :: value
          real(8), intent(in), optional :: delta

          integer :: binarysearch

          integer :: left, middle, right, orientation
          real(8) :: d

          if (present(delta) .eqv. .true.) then
             d = delta
          else
             d = 1d-9
          endif

          orientation = 1
          if (array(2) < array(1)) orientation = -1
       
          left = 1
          right = length
          do
             if (left > right) exit
             middle = nint((left+right) / 2.0d0)
             if ( abs(array(middle) - value) <= d) then
                binarySearch = middle
                return
             endif
             select case(orientation)
             case(1) 
                if (array(middle) > value) then
                   right = middle - 1
                else
                   left = middle + 1
                end if
             case(-1)
                if (array(middle) < value) then
                   right = middle - 1
                else
                   left = middle + 1
                end if
             end select
          end do
          binarysearch = right

      END FUNCTION binarysearch
!-----------------------------------------------------------------------
      SUBROUTINE bl_interp(xp, x_array, yp, y_array, x, y, ii, w)
          ! This function uses bilinear interpolation to get the
          ! interpolation weights, w and the indices for interpolation
          ! Assumed to be sampled on a regular grid, with the grid x values
          ! specified by x_array and the grid y values specified by y_array
          implicit none
          integer, intent(in) :: xp, yp           
          real(8), dimension(xp), intent(in) :: x_array
          real(8), dimension(yp), intent(in) :: y_array
          real(8), intent(in) :: x,y
          real(8), intent(out) :: w(4)
          integer, intent(out) :: ii(4)
      
          real(8) :: denom, x1, x2, y1, y2
          real(8) :: x2x1, y2y1, x2x, y2y, xx1, yy1
          integer  :: i, j, ir, jr

          i = binarysearch(xp, x_array, x)
          j = binarysearch(yp, y_array, y)
          ! Make sure we start from one
          if (i.eq.0) then
             i = 1; ir = 1
          else
             ir = i + 1
             if (ir > xp) then
                if (3d0*x_array(1) - 2d0*x_array(2) + 
     &              360d0 < x_array(xp)) then
                   ! Wrap the longitude around
                   ir = 1
                else
                   ir = xp
                endif
             endif
          endif
          x1 = x_array(i)
          x2 = x_array(ir)
          
          if (j.eq.0) then
             j = 1; jr = 1
          else
             jr = j + 1
             if (jr > yp) jr = yp
          endif
          y1 = y_array(j)
          y2 = y_array(jr)

          ! Lon, lat to meters 
          if (ir.eq.i) then
             x2x1 = 1.0d0
             x2x  = 0.0d0
             xx1  = 1.0d0
          else
             x2x1 = haversine(x1,x2,y,y)
             x2x  = haversine(x,x2,y,y)
             xx1  = haversine(x1,x,y,y)
          endif
          if (jr.eq.j) then
             y2y1 = 1.0d0
             y2y  = 0.0d0
             yy1  = 1.0d0
          else
             y2y1 = haversine(x,x,y1,y2)
             y2y  = haversine(x,x,y,y2)
             yy1  = haversine(x,x,y1,y)
          endif 
          denom = x2x1*y2y1
         
          if (denom < 1e-12) then
             write(16,*) 'Denominator is small',denom
             write(16,*) x,y,x1,y1,x2,y2,i,j,ir,jr
          endif
            
          ! Output the indices
          ii   = [i, j, ir, jr]

          ! Output the weights
          w(1) = x2x*y2y/denom
          w(2) = xx1*y2y/denom
          w(3) = x2x*yy1/denom
          w(4) = xx1*yy1/denom

      END SUBROUTINE bl_interp
!-----------------------------------------------------------------------
      SUBROUTINE bl_interp2(xp, yp, x_mat, y_mat, x, y, ii, w, tree)
          ! This function uses bilinear interpolation to get the
          ! interpolation weights, w and the indices for interpolation
          ! Assumed to be sampled on a structured unregular grid x values
          ! specified by x_mat and the grid y values specified by y_mat
          use kdtree2_module, only: kdtree2, kdtree2_result, 
     &                              kdtree2_n_nearest
          implicit none
          integer, intent(in) :: xp, yp           
          real, dimension(xp,yp), intent(in) :: x_mat, y_mat
          real(8), intent(in) :: x,y
          real(8), intent(out) :: w(4)
          integer, intent(out) :: ii(4)
          real(8), dimension(xp) :: x_array
          real(8), dimension(yp) :: y_array
          real(8) :: x11, x21, x12, x22, x_vec(yp)
          real(8) :: y11, y21, y12, y22, y_vec(yp)
          real(8) :: d11, d12, d21, d22, denom
          real(8) :: p(2), a(2), b(2), c(2), d(2)
          integer  :: i, j, jj, ir, jr, ia, xi, yi, ind
          type(kdtree2), pointer :: tree
          type(kdtree2_result) :: kdresult(1)

          ! find the closest point to this node and convert to subscript
          call kdtree2_n_nearest(tp=tree,qv=[x,y],
     &                           nn=1,results=kdresult)
          ind = kdresult(1)%idx  ! location of closest point
          xi  = ind/yp           ! ind2sub
          yi  = mod(ind,yp)      ! ind2sub
          if (yi == 0) then 
             yi = yp
          else
             xi = xi + 1 
          endif
 
          ii   = [-1, 0, 0, 0]; ! Default indices
          ia = 0
          ! loop over points close to the nearest one  
          do i = max(1,xi-1),min(xi+1,xp) !1,xp
             ir = i + 1
             if (ir > xp) then
                ! If wrapping around
                if (3d0*x_mat(1,1) - 2d0*x_mat(2,1) +
     &              360d0 < x_mat(xp,1)) then
                   ir = 1; ia = 1
                else
                   exit
                endif
             endif
             do j = max(1,yi-1),min(yi+1,yp-1) !1,yp-1
                jr = j + 1
                x11 = x_mat(i,j);  x21 = x_mat(ir,j) + ia*360d0
                x12 = x_mat(i,jr); x22 = x_mat(ir,jr) + ia*360d0
                y11 = y_mat(i,j);  y21 = y_mat(ir,j)
                y12 = y_mat(i,jr); y22 = y_mat(ir,jr)

                ! Cycle if outside bbox
                if (x.lt.min(x11,x12).or.x.gt.max(x21,x22).or.
     &              y.lt.min(y11,y21).or.y.gt.max(y12,y22)) cycle

                ! Confirm if in the quadrilateral or not
                p = [x,y];
                a = [x11, y11]; b = [x21, y21];
                c = [x22, y22]; d = [x12, y12];
                if (PtInTri(p,a,b,c).or.PtInTri(p,c,d,a)) then
                   ! Lon, lat to meters 
                   d11 = 1d0/haversine(x11,x,y11,y)
                   d12 = 1d0/haversine(x12,x,y12,y)
                   d21 = 1d0/haversine(x21,x,y21,y)
                   d22 = 1d0/haversine(x22,x,y22,y)
                   denom  = d11 + d12 + d21 + d22 

                   ! Output the indices
                   ii   = [i, j, ir, jr]

                   ! Output the weights
                   w(1) = d11/denom
                   w(2) = d21/denom
                   w(3) = d12/denom
                   w(4) = d22/denom

                   return; 
                endif  
             enddo
          enddo

          return;

      END SUBROUTINE bl_interp2

!-----------------------------------------------------------------------
      SUBROUTINE conserve_interp(xp, x_array, yp, y_array, 
     &           x1, y1, x2, y2, x3, y3, ii, ww)
          ! This function works out the indices of x_array and y_array 
          ! inside the triangle specified by x1, y1, x2, y2, x3, y3
          implicit none
          integer, intent(in) :: xp, yp           
          real(8), dimension(xp), intent(in) :: x_array
          real(8), dimension(yp), intent(in) :: y_array
          real(8), intent(in) :: x1 ,y1, x2, y2, x3, y3
          integer, intent(out) :: ii(4)
          real(8), intent(out) :: ww(4)
      
          integer  :: is, ie, js, je
          real(8) :: x1xm, x2xm, y1ym, y2ym

          ! Get starting and ending indices for loop
          is = binarysearch(xp, x_array, x2)
          js = binarysearch(yp, y_array, y2)
          ie = binarysearch(xp, x_array, x3) 
          je = binarysearch(yp, y_array, y3)
          is = is + 1
          js = js + 1
          if (is > ie) then
             ! Find nearest
             is   = binarysearch(xp, x_array, x1)
             x1xm = haversine(x_array(is),x1,y1,y1)
             x2xm = haversine(x_array(min(xp,is+1)),x1,y1,y1)
             if (x1xm > x2xm) then
                is = is+1
             endif
             ie = is
          endif
          if (js > je) then
             ! Find nearest
             js   = binarysearch(yp, y_array, y1)
             x1xm = haversine(x1,x1,y_array(js),y1)
             x2xm = haversine(x1,x1,y_array(min(yp,js+1)),y1)
             if (x1xm > x2xm) then
                js = js+1
             endif
             je = js
          endif
          !if ((ie-is)+(je-js).le.4) then
          !   ! Revert to linear interp search size similar to grid size
          !call bl_interp(xp, x_array, yp, y_array, x1, y1, ii, ww)
          !else
          !   ! Output the indices
          ii   = [is, js, ie, je]
          ww   = 1d4
          !endif
      END SUBROUTINE conserve_interp
!-----------------------------------------------------------------------
      function PtinTri(p, a, b, c)
        real(8), intent(in)  :: p(2), a(2), b(2), c(2)
        logical :: PtinTri

        if (SameSide(p, a, b, c).and.Sameside(p, b, a, c).and.
     &      SameSide(p, c, a, b)) then
            PtinTri = .true.
        else
            PtinTri = .false.
        endif

      end function PtinTri
!-----------------------------------------------------------------------
      function SameSide(p1, p2, a, b)
        real(8), intent(in)  :: p1(2), p2(2), a(2), b(2)
        real(8) :: cp1, cp2
        logical :: SameSide
  
        ! cross (b-a) with (p1-a) 
        cp1 = (b(1)-a(1))*(p1(2)-a(2)) - (b(2)-a(2))*(p1(1)-a(1))
        ! cross (b-a) with (p2-a) 
        cp2 = (b(1)-a(1))*(p2(2)-a(2)) - (b(2)-a(2))*(p2(1)-a(1))
        if (cp1*cp2 >= 0d0) then
            SameSide = .true.
        else
            SameSide = .false.
        endif

      end function SameSide
!-----------------------------------------------------------------------
      function calcarea(x1,y1,x2,y2,x3,y3)
        real(8), intent(in) :: x1 ,y1, x2, y2, x3, y3
        real(8) :: calcarea

        calcarea = 0.5d0*abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)) 

      end function calcarea
!-----------------------------------------------------------------------
      SUBROUTINE terminate(NO_MPI_FINALIZE)
#ifdef CMPI
      USE MESSENGER
#endif
      IMPLICIT NONE
C
      LOGICAL, OPTIONAL :: NO_MPI_FINALIZE
C
      call setMessageSource("terminate")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
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
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      call unsetMessageSource()
      END SUBROUTINE terminate

!*********************************************************************************C
!                                                                                 C
!  Spherical coordinates higher level functions & subroutines.                    C 
!                                                                                 C
!   Use data declared in the module                                               C
!                                                                                 C
!*********************************************************************************C
      SUBROUTINE READFORTROTS( ) 
        IMPLICIT NONE

        !
        CHARACTER (LEN=1024):: FNAME           
        INTEGER:: FUNIT, errorIO
        CHARACTER (len=100) :: errorVar 
        CHARACTER (len=100) :: inputtype

        INTEGER:: II
        REAL(8):: PHIN, THETAN 
        CHARACTER (LEN=300):: IGNOREMSG, MSG

        ! Euler angle
        REAL(8):: AL, BE, GAM

        call setMessageSource("ReadFortRots")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif
        
        FUNIT = 1002 ;    
        FNAME = TRIM(GBLINPUTDIR)//'/'//'fort.rotm' ;
        errorVar = ""
        call openFileForRead( FUNIT, FNAME, errorIO) ;

        inputtype = "" ; 
        read( FUNIT, *) inputtype ;

        scratchMessage = 
     &       "**** Solve SWE in the rotated coordinate frame ****" ; 
        call allMessage(INFO,scratchMessage)

        errorVar = "inputype" ;
        SELECT CASE(adjustl(trim(inputtype)))
        CASE ("znorth_in_spherical_coors")
           READ(FUNIT, *) PHIN, THETAN, IGNOREMSG ; 

           WRITE(scratchMessage,'(A)') "A vector describing the"//
     &            " north-south axis in spherical coordinates "//
     &            "(\varphi,\theta) is given."
           CALL allMessage(INFO,scratchMessage) 
           WRITE(scratchMessage,'(A,2F20.8)') " (varphin,thetan) = ",
     &                                        PHIN, THETAN
           CALL allMessage(INFO,scratchMessage)

           CALL DRV_ROTMAT_ZNVEC( PHIN, THETAN ) ;           
        CASE ("z-x-z")
           READ(FUNIT, *) AL, BE, GAM, IGNOREMSG ; 
           scratchMessage = "Intrinsic rotation z-x-z (alpha,beta,"//
     &             "gamma) is given." ;
           CALL allMessage(INFO,scratchMessage)
           WRITE(scratchMessage,'(A,3F20.8)') " (alpha,beta,gamma) = ",
     &           AL, BE, GAM
           CALL allMessage(INFO,scratchMessage)

           CALL DRV_ROTMAT_ZXZROT( AL, BE, GAM ) ;          
        CASE ("rotation_matrix")           
           DO II = 1, 3
              READ(FUNIT,*) ROTMAT(II,3) ;
           END DO
           
           scratchMessage = "ROTATION MATRIX: "
           call allMessage(INFO,scratchMessage)
           scratchMessage =  
     &       "****************************************" ;  
           call allMessage(INFO,scratchMessage)
           scratchMessage = "A rotation matrix R: x -> x' is given"
           call allMessage(INFO,scratchMessage)
           DO II = 1, 3 
              WRITE(scratchMessage,'(A,3F12.4,A)') '| ', 
     &            ROTMAT(II,:), ' |' ; 
              call allMessage(INFO,scratchMessage)
           END DO
           scratchMessage =  
     &       "****************************************" ;  
           call allMessage(INFO,scratchMessage)
           
        CASE DEFAULT
           WRITE(scratchMessage,'(A)') 
     &     "Error in readfortrots(): defining"// 
     &     " rotation matrix through "//adjustl(trim(inputtype))// 
     &     " is either not valid or has not been implemented" 
           CALL allMessage(ERROR,scratchMessage)
           scratchMessage = "ADCIRC will stop"
           CALL allMessage(ERROR,scratchMessage)
#if defined(MESH_TRACE) || defined(ALL_TRACE)
           call allMessage(DEBUG,"Return.") 
#endif
           call unsetMessageSource()
           CALL EXIT(1)          
        END SELECT

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") 
#endif
      call unsetMessageSource()

        CLOSE(FUNIT) 

        RETURN       
      END SUBROUTINE READFORTROTS
      
      SUBROUTINE ALLOC_SPROTCOORS_MEM( )  
        IMPLICIT NONE
      
        ALLOCATE( SLAMR(MNP), SFEAR(MNP) ) ;
        ALLOCATE( RVELF(2,2,MNP), RVELI(2,2,MNP) ) ; 
        ALLOCATE( UVECTMP(MNP), VVECTMP(MNP) ) ; 
      
        SLAMR = 0.0D0 ; 
        SFEAR = 0.0D0 ;
        
        RVELF = 0.0D0 ; 
        RVELI = 0.0D0 ;

        UVECTMP = 0.0D0 ; 
        VVECTMP = 0.0D0 ; 
        
! To be used in neighb()
        ALLOCATE( EXVR(3,MNE), EYVR(3,MNE) ) ; 
        EXVR = 0.0D0 ; 
        EYVR = 0.0D0 ;
        
        RETURN ;    
      END SUBROUTINE ALLOC_SPROTCOORS_MEM
      
      SUBROUTINE COMP_ROTSPCOORD_VARS( ICS )
      IMPLICIT NONE
      
        INTEGER:: ICS
      
        !c local c!
        INTEGER:: II, IERR

        call setMessageSource("COMP_ROTSPCOORD_VARS")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif

        CALL CHECK_RTOTMAT( ROTMAT, IERR ) ; 
        
        CALL SPCOORSROTS( ROTMAT, SLAMR, SFEAR, SLAM, SFEA, NP ) ;
        
        CALL SPVECROTSMAT( RVELF, RVELI, ROTMAT, 
     &                     SLAMR, SFEAR, SLAM, SFEA, NP ) ;
          
        DO II = 1, NP
           CALL CYLINDERMAP( X(II), Y(II), 
     &          SLAMR(II), SFEAR(II), SLAM0, SFEA0, ICS ) ;
        END DO

#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.") 
#endif
        call unsetMessageSource()
        
        RETURN 
      END SUBROUTINE COMP_ROTSPCOORD_VARS

      SUBROUTINE GET_ROTSPCOORD_CORIFVAL( CORIF ) 
        IMPLICIT NONE
      
        REAL(8) :: CORIF(:)

        call setMessageSource("GET_ROTSPCOORD_CORIFVAL")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif
        
        CALL CORIFVAL_ROTSPCOORD(CORIF, ROTMAT, SLAMR, SFEAR, NP) ; 

#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.") 
#endif
        call unsetMessageSource()

        RETURN ;
      END SUBROUTINE GET_ROTSPCOORD_CORIFVAL

      ! 
      SUBROUTINE CHECK_RTOTMAT( RTOTS, IERR )
        IMPLICIT NONE

        INTEGER:: IERR
        REAL(8), intent(in):: RTOTS(3,3)

        !c local c!
        REAL(8) :: RTT(3,3), TR, SOFFD

        call setMessageSource("CHECK_RTOTMAT")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif

        RTT = TRANSPOSE( RTOTS ) ; 

        RTT = MATMUL(RTT, RTOTS ) ; 

        ! trace of a matrix RTT 
        TR = RTT(1,1) + RTT(2,2) + RTT(3,3) ; 
       
        RTT(1,1) = 0.0 ;
        RTT(2,2) = 0.0 ;
        RTT(3,3) = 0.0 ;
        
        SOFFD = SUM( MATMUL(ABS(RTT), (/ 1.0D0, 1.0D0, 1.0D0 /)) )  ;

        IERR = 0 ;
        IF ( SOFFD > 1.0e-10 ) THEN
           IERR = 1 ; 
           
           WRITE(6,*) "Error: a given rotation matrix is not"// 
     &          " orthgonal and thus invalid" ;
        END IF

#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.") 
#endif
        call unsetMessageSource()

        RETURN ;
      END SUBROUTINE CHECK_RTOTMAT

      !c Define a transformation matrix from a given vector in   
      !c sherical coordinates (1,\varphi,\theta) (\varphi,\theta 
      !c are in degree) 
      SUBROUTINE DRV_ROTMAT_ZNVEC( VARPHI, THETA )
        IMPLICIT NONE

        REAL(8) :: VARPHI, THETA ! in deg

        !c local c!
        REAL(8) :: RADPHI, RADTHETA, ROTMATOUT(3)

        call setMessageSource("DRV_ROTMAT_ZNVEC")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif

        !c deg to rad
        RADPHI = VARPHI*DEG2RAD ;
        RADTHETA = THETA*DEG2RAD ;   

        CALL GET_ROTMAT_ZNVEC( ROTMAT, RADPHI, RADTHETA )   ;

        WRITE(scratchMessage,'(A)') "ROTATION MATRIX: " ;
        CALL allMessage(INFO,scratchMessage)
        scratchMessage =  
     &       "****************************************" ;  
        call allMessage(INFO,scratchMessage)
        ROTMATOUT = ROTMAT(1,:)
        WRITE(scratchMessage,'(A,3F12.4,A)') '| ', ROTMATOUT, ' |' 
        CALL allMessage(INFO,scratchMessage)
        ROTMATOUT = ROTMAT(2,:)
        WRITE(scratchMessage,'(A,3F12.4,A)') '| ', ROTMATOUT, ' |' 
        CALL allMessage(INFO,scratchMessage)
        ROTMATOUT = ROTMAT(3,:)
        WRITE(scratchMessage,'(A,3F12.4,A)') '| ', ROTMATOUT, ' |' 
        CALL allMessage(INFO,scratchMessage)
        scratchMessage =  
     &       "****************************************" ;  
        call allMessage(INFO,scratchMessage)

#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.") 
#endif
        call unsetMessageSource()

        RETURN 
      END SUBROUTINE DRV_ROTMAT_ZNVEC

      SUBROUTINE DRV_ROTMAT_ZXZROT( ALPHA, BETA, GAMMA )
         IMPLICIT NONE

         REAL(8) :: ALPHA, BETA, GAMMA

         ! local !
         REAL(8) :: RAL, RBE, RGAM

         call setMessageSource("DRV_ROTMAT_ZXZROT")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Enter.")
#endif

         RAL = ALPHA*DEG2RAD ;
         RBE = BETA*DEG2RAD ;
         RGAM = GAMMA*DEG2RAD ;

         CALL GET_ROTMAT_ZXZROT( ROTMAT, RAL, RBE, RGAM )   ;
         
         WRITE(scratchMessage,'(A)') "ROTATION MATRIX: " ;
         CALL allMessage(INFO,scratchMessage)
         scratchMessage =  
     &       "****************************************" ;  
         call allMessage(INFO,scratchMessage)
         WRITE(scratchMessage,'(A,3F12.4,A)') '| ', ROTMAT(1,:), ' |' 
         CALL allMessage(INFO,scratchMessage)
         WRITE(scratchMessage,'(A,3F12.4,A)') '| ', ROTMAT(2,:), ' |' 
         CALL allMessage(INFO,scratchMessage)
         WRITE(scratchMessage,'(A,3F12.4,A)') '| ', ROTMAT(3,:), ' |' 
         CALL allMessage(INFO,scratchMessage)
         scratchMessage =  
     &       "****************************************" ;  
         call allMessage(INFO,scratchMessage)

#if defined(MESH_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.") 
#endif
         call unsetMessageSource()
         
         RETURN ;
      END SUBROUTINE DRV_ROTMAT_ZXZROT

      !c
      !c Get Lon-Lat in the rotated frame
      SUBROUTINE DRVSPCOORSROTS1( LONR, LATR, LONO, LATO, NN )
        IMPLICIT NONE

        INTEGER:: NN
        REAL(8), dimension(:):: LONR, LATR, LONO, LATO

        call setMessageSource("DRVSPCOORSROTS1")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif
        
        CALL SPCOORSROTS1( ROTMAT, LONR, LATR, LONO, LATO, NN )

#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.") 
#endif
        call unsetMessageSource()

        RETURN ;
      END SUBROUTINE DRVSPCOORSROTS1

      SUBROUTINE DRVSPCOORSROTS0( LONR, LATR, LONO, LATO )
        IMPLICIT NONE

        REAL(8) :: LONR, LATR, LONO, LATO
        
        call setMessageSource("DRVSPCOORSROTS0")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif

        CALL SPCOORSROTS0( ROTMAT, LONR, LATR, LONO, LATO )

#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.") 
#endif
        call unsetMessageSource()

        RETURN ;
      END SUBROUTINE DRVSPCOORSROTS0


!*********************************************************************************C
!                                                                                 C
!  Lower level functions & subroutines.  Do not use data declared in the module   C
!                                                                                 C
!*********************************************************************************C
      ! 
      !  Corilois values in the (x',y',z') at (\varphi',\theta')
      ! 
      SUBROUTINE CORIFVAL_ROTSPCOORD( CORIF, ROTM, VARPHIR, THETAR, NN )
         IMPLICIT NONE
      
         !     dummy !
        INTEGER:: NN
        REAL(8) :: ROTM(3,3)
        REAL(8), DIMENSION(:):: CORIF, VARPHIR, THETAR
        
        ! local c!
        INTEGER:: II
        REAL(8) :: EOMEGA2, ROT3(3)

        call setMessageSource("CORIFVAL_ROTSPCOORD")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif
        
        ROT3 = ROTM(:,3) ; 

        !c Angular speed of Earth  
        EOMEGA2 = 2.0D0*omega;
        DO II = 1, NN
           CORIF(II) = ROT3(1)*cos( THETAR(II) )*cos( VARPHIR(II) ) + 
     &          ROT3(2)*cos( THETAR(II) )*sin( VARPHIR(II) ) +
     &          ROT3(3)*sin( THETAR(II) ) ; 
             
           CORIF(II) = EOMEGA2*CORIF(II) ;              
        END DO

#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.") 
#endif
        call unsetMessageSource()
        
        RETURN ;
      END SUBROUTINE CORIFVAL_ROTSPCOORD

      !
      ! FWD = .TRUE.
      !   RVELF*(UUS,VVS) -> (UUD,VVD)
      !
      ! FWD = .FALSE. 
      !   RVELI*(UUS,VVS) -> (UUD,VVD)
      SUBROUTINE DRVMAP2DSPVEC( UUD, VVD, UUS, VVS, NN, FWD )
        IMPLICIT NONE

        INTEGER:: NN
        LOGICAL:: FWD
        REAL(8), DIMENSION(:):: UUD, VVD, UUS, VVS

        call setMessageSource("DRVMAP2DSPVEC")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif

        IF ( FWD ) THEN
           ! forward transform
           CALL MAP2DSPVECS( UUD, VVD, UUS, VVS, RVELF, NN )   
        ELSE
           ! inverse transform
           CALL MAP2DSPVECS( UUD, VVD, UUS, VVS, RVELI, NN )   
        END IF

#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.") 
#endif
        call unsetMessageSource()

        RETURN ;
      END SUBROUTINE DRVMAP2DSPVEC
      !
      ! Map two vectors in spherical coordinates 
      ! 
      !    VR(:,II) = SPRV(:,:,II)*VO(:,II)
      SUBROUTINE MAP2DSPVECS( VPHIR, VTHETAR, VPHIO, VTHETAO, SPRV, NN )
        IMPLICIT NONE

        ! dummy !
        INTEGER:: NN
        REAL(8) :: SPRV(:,:,:)
        REAL(8), DIMENSION(:):: VPHIR, VTHETAR
        REAL(8), DIMENSION(:):: VPHIO, VTHETAO

        ! local !
        INTEGER:: II
        REAL(8) :: VO(2), VR(2)

        call setMessageSource("MAP2DSPVECS")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif
 
        DO  II = 1, NN
           VO = (/ VPHIO(II), VTHETAO(II) /)  ;
           VR = MATMUL( SPRV(1:2,1:2,II), VO ) ; 
           
           VPHIR(II) = VR(1) ;
           VTHETAR(II) = VR(2) ;
        END DO
        
#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.") 
#endif

        call unsetMessageSource()
        RETURN ;
      END SUBROUTINE MAP2DSPVECS
      ! Given a radius vector (\varphi,\theta) (in rad)
      ! of z', the new z after rotation, return
      ! a transformation matrix:
      !        RTOT: (x,y,z) -> (x',y',z')  
      SUBROUTINE GET_ROTMAT_ZNVEC( RTOT, VARPHIN, THETAN )   
        IMPLICIT NONE
      
        REAL(8) :: RTOT(:,:) 
        REAL(8) :: VARPHIN, THETAN

        call setMessageSource("GET_ROTMAT_ZNVEC")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif
      
        RTOT = 0.0D0 ;
      
        RTOT(1,1) = sin( VARPHIN ) ; 
        RTOT(2,1) = sin( THETAN )*cos( VARPHIN ) ;
        RTOT(3,1) = cos( THETAN )*cos( VARPHIN ) ; 
      
        RTOT(1,2) = -cos( VARPHIN ) ; 
        RTOT(2,2) = sin( THETAN )*sin( VARPHIN ) ;
        RTOT(3,2) = cos( THETAN )*sin( VARPHIN ) ; 
        
        RTOT(1,3) = 0.0D0 ;
        RTOT(2,3) = -cos( THETAN ) ; 
        RTOT(3,3) =  sin( THETAN ) ; 

#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.") 
#endif
        call unsetMessageSource()
        
        RETURN ;
      END SUBROUTINE GET_ROTMAT_ZNVEC


      SUBROUTINE GET_ROTMAT_ZXZROT( RTOT, AL, BE, GAM )
        IMPLICIT NONE

        ! dummy !
        REAL(8) :: RTOT(:,:) 
        REAL(8) :: AL, BE, GAM

        ! local !
        REAL(8) :: RTMP(3,3), RM2D(2,2)
        
        call setMessageSource("GET_ROTMAT_ZXZROT")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif

        ! A rotation around z-axis
        RM2D = RTOT2DMAT( AL ) ;

        RTMP = 0.0D0 ;
        RTMP(3,3) = 1.0D0 ;
        RTMP(1:2,1:2) = RM2D ; 

        RTOT = RTMP ; 

        ! A rotation around x'-axis    
        RM2D = RTOT2DMAT( BE ) ;

        RTMP = 0.0D0 ; 
        RTMP(1,1) = 1.0D0 ;
        RTMP(2:3,2:3) = RM2D ; 

        RTOT = MATMUL( RTOT, RTMP ) ; 

        ! A rotation around z''-axis
        RM2D = RTOT2DMAT( GAM ) ; 

        RTMP = 0.0D0 ;
        RTMP(3,3) = 1.0D0 ;
        RTMP(1:2,1:2) = RM2D ; 

        RTOT = MATMUL( RTOT, RTMP ) ; 

#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.") 
#endif
        call unsetMessageSource()

        RETURN ;        
      END SUBROUTINE GET_ROTMAT_ZXZROT

      FUNCTION RTOT2DMAT( ALPHA ) RESULT( RM )
        IMPLICIT NONE

        REAL(8) :: ALPHA
        REAL(8) :: RM(2,2) 

        call setMessageSource("RTOT2DMAT")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif
           
        RM(1,1) = cos(ALPHA) ;
        RM(2,1) = -sin(ALPHA) ;
        RM(1,2) = sin(ALPHA) ;
        RM(2,2) = cos(ALPHA) ;

#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.") 
#endif
        call unsetMessageSource()
        
        RETURN ; 
      END FUNCTION RTOT2DMAT


      ! (e_{r},e_{\varphi},e_{\theta})^{T} = T (i,j,k)^{T}
      !  x' = R x
      ! (e_{r'},e_{\varphi'},e_{\theta'})^{T} = T' (i',j',k')^{T}
      !
      ! V'(r',\varphi',\theta') = T' R T^{T} V(r,\varphi,\theta)
      !
      ! RVELF = T' R T^{T}
      ! RVELI = (T')^{R} R^{T} T
      !
      ! Evaluated the forward and inverse transformation matrices
      ! between two spherical coordinate
      !     RVELF: VEC(LONO,LATO) -> VEC(LONR,LATR)
      !     RVELI: VEC(LONR,LATR) -> VEC(LONO,LATO)
      !
      !  lonr,latr in rad
      !
      SUBROUTINE SPVECROTSMAT( RVELF, RVELI, RTOTS, 
     &                         LONR, LATR, LONO, LATO, NN )
         IMPLICIT NONE

         ! Dummys !
         REAL(8), DIMENSION(:,:,:):: RVELF, RVELI 

         INTEGER :: NN
         REAL(8) :: RTOTS(3,3)
         REAL(8), DIMENSION(:):: LONR, LATR, LONO, LATO          

         ! local !
         INTEGER:: IP
         REAL(8), DIMENSION(3,3):: TP, TM, TRT 

         call setMessageSource("SPVECROTSMAT")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Enter.")
#endif

         DO IP = 1, NN
           TRT = 0.0D0 ; 
        
           ! T' = T(varphi',theta')
           TP = CT2SP_VECTRANSMAT( LONR(IP), LATR(IP) ) ;
 
           ! T = T(varphi,theta)
           TM = CT2SP_VECTRANSMAT( LONO(IP), LATO(IP) ) ; 
           TM = TRANSPOSE( TM ) ; ! T^{T}

           ! T' R T^{T}
           TRT = MATMUL( RTOTS, TM ) ; 
           TRT = MATMUL( TP, TRT ) ; 

           RVELF(1:2,1:2,IP) = TRT(2:3,2:3) ; ! forward transform
           RVELI(:,:,IP) = TRANSPOSE( RVELF(:,:,IP) ) ; ! inverse transform
         END DO

#if defined(MESH_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.") 
#endif
         call unsetMessageSource()

         RETURN ; 
      END SUBROUTINE SPVECROTSMAT

      !
      ! Transformation matrix mapping a vector from the Cartesian to spherical coordinates 
      !
      ! [e_{r},e_{varphi},e_{theta}]^{T} = TM [ i, j, k]^{T}
      FUNCTION CT2SP_VECTRANSMAT( VARPHI, THETA ) RESULT(TM)
        IMPLICIT NONE

        REAL(8) :: TM(3,3)
        REAL(8), intent(in):: VARPHI, THETA

        call setMessageSource("CT2SP_VECTRANSMAT")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif
        
        TM(1,1) = cos( THETA )*cos( VARPHI ) ;
        TM(2,1) = -sin( VARPHI ) ; 
        TM(3,1) = -sin( THETA )*cos( VARPHI ) ; 
        
        TM(1,2) = cos( THETA )*sin( VARPHI ) ; 
        TM(2,2) = cos( VARPHI ) ; 
        TM(3,2) = -sin( THETA )*sin( VARPHI ) ; 
        
        TM(1,3) = sin( THETA ) ; 
        TM(2,3) = 0.0D0 ; 
        TM(3,3) = cos( THETA ) ;    

#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.") 
#endif
        call unsetMessageSource()
        
        RETURN ;   
      END FUNCTION CT2SP_VECTRANSMAT 
      
      SUBROUTINE SPCOORSROTS1( RTOTS, LONR, LATR, LONO, LATO, NN )
        IMPLICIT NONE

        !c dummy c!
        INTEGER :: NN
        REAL(8) :: RTOTS(3,3)
        REAL(8), dimension(:):: LONR, LATR, LONO, LATO

        !c local c!
        INTEGER :: IP
        REAL(8) :: XP(3), XPR(3), LLO, LTO

        call setMessageSource("SPCOORSROTS1")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif

        DO IP = 1, NN
           LLO = LONO( IP ) ; 
           LTO = LATO( IP ) ; 
 
           XP(1) = cos( LTO )*cos( LLO ) ; 
           XP(2) = cos( LTO )*sin( LLO ) ; 
           XP(3) = sin( LTO ) ; 

           XPR = MATMUL( RTOTS, XP) ;

           LATR(IP) = atan2( XPR(3), SQRT(XPR(1)*XPR(1) 
     &                                  + XPR(2)*XPR(2)) )
           LONR(IP) = atan2( XPR(2), XPR(1) ) ;   
        END DO

#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.") 
#endif
        call unsetMessageSource()
        
        RETURN ; 
      END SUBROUTINE SPCOORSROTS1

      SUBROUTINE SPCOORSROTS0( RTOTS, LONR, LATR, LONO, LATO )
        IMPLICIT NONE

        !c dummy c!
        REAL(8) :: RTOTS(3,3)
        REAL(8) :: LONR, LATR, LONO, LATO

        !c local c!
        INTEGER:: IP
        REAL(8) :: XP(3), XPR(3), LLO, LTO

        call setMessageSource("SPCOORSROTS0")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif

        LLO = LONO  ; 
        LTO = LATO  ; 
 
        XP(1) = cos( LTO )*cos( LLO ) ; 
        XP(2) = cos( LTO )*sin( LLO ) ; 
        XP(3) = sin( LTO ) ; 
        
        XPR = MATMUL( RTOTS, XP) ;
        
        LATR = atan2( XPR(3), SQRT(XPR(1)*XPR(1) + XPR(2)*XPR(2)) ) ;
        LONR = atan2( XPR(2), XPR(1) ) ;   

#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.") 
#endif
        call unsetMessageSource()

        RETURN ; 
      END SUBROUTINE SPCOORSROTS0
!******************************************************************************
!                                                                             *
!    Transform from lon,lat (lamda,phi) coordinates into the cylindrical      *
!    mapping coordinates.                                                     *
!    Lon,Lat must be in radians.                                              *
!                                                                             *
!******************************************************************************
      SUBROUTINE CYLINDERMAP(X,Y,RLAMBDA,PHI,RLAMBDA0,PHI0,MTYPE)
      IMPLICIT NONE
      
      INTEGER, INTENT(IN):: MTYPE
      REAL*8 X,Y,RLAMBDA,PHI,RLAMBDA0,PHI0

      call setMessageSource("CYLINDERMAP")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      
      SELECT CASE(MTYPE)
      CASE (2,21) !c CPP
         X= Rearth*(RLAMBDA-RLAMBDA0)*COS(PHI0)
         Y= PHI*Rearth
      CASE (20) !c Equal area map
         X= Rearth*(RLAMBDA-RLAMBDA0)*COS(PHI0)
         Y= Rearth*SIN(PHI)/COS(PHI0) ; 
      CASE (22) ! Mercartor
         X= Rearth*(RLAMBDA-RLAMBDA0)*COS(PHI0)  ;
         Y= Rearth*log(tan(PHI) + 1.0d0/cos(PHI))*cos(PHI0) ;          
      CASE (23) ! Miller, Compromised Mercator'
         X= Rearth*(RLAMBDA-RLAMBDA0)*COS(PHI0)  ;
         Y= 1.25d0*Rearth*log(tan(0.25D0*PI + 0.4d0*PHI))*cos(PHI0)
      CASE (24) ! Gall stereographic
         X = 0.5D0*SQRT(2.D0)*RLAMBDA ;
         Y = Rearth*(1.0 + 0.5D0*SQRT(2.D0))*tan( 0.5D0*PHI )
      CASE DEFAULT
         WRITE(6,'(A)') "Error in cylindermap(): "
         WRITE(6,'(A,I1,A)') "   Mapping type =  ", MTYPE, 
     &                      " has not been implemented" ;
         WRITE(6,'(A)') "  Valid MTYPE = 2, 20-24" ;
      END SELECT
      
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") 
#endif
      call unsetMessageSource()

      RETURN ;
      END SUBROUTINE CYLINDERMAP 

      SUBROUTINE SFAC_ELEAVG( SFAELES, SFNODES, NM, NE ) 

      IMPLICIT NONE

      INTEGER :: NM(:,:), NE
      REAL(8), dimension(:) :: SFAELES, SFNODES
      
      ! local
      INTEGER :: IE, NMI1, NMI2, NMI3
      REAL(8) :: ONETRD 

      call setMessageSource("SFAC_ELEAVG")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      
      ONETRD = 1.0D0/3.0D0 ;
      DO IE = 1, NE
         NMI1 = NM(IE,1) ; 
         NMI2 = NM(IE,2) ; 
         NMI3 = NM(IE,3) ;  
         
         SFAELES(IE) = (SFNODES(NMI1) + 
     &        SFNODES(NMI2) + SFNODES(NMI3))*ONETRD ; 
      END DO

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") 
#endif
      call unsetMessageSource()
      
      RETURN ;
      END SUBROUTINE SFAC_ELEAVG 
    
      SUBROUTINE CAL_ELXV_SPCOOR( XVE, YVE, LONVE, LATVE, ICS, SFLAG )

      USE CONSTANTS, ONLY: RAD2DEG

      IMPLICIT NONE

      !c Dummy 
      !c   XVE, YVE     -- output
      !c   LONVE, LATVE -- input
      INTEGER, optional:: SFLAG
      INTEGER:: ICS
      REAL(8), DIMENSION(:):: XVE, YVE, LONVE, LATVE
      
      !c Local 
      INTEGER :: II, IDX, SPF
      REAL(8) :: XC, YC
      REAL(8) :: Jac1, Jac2
      REAL(8), dimension(3):: LONM, LATM, LONTMP
      REAL(8), dimension(3):: DLX, DL
      REAL(8):: DLED1, DLED2

      call setMessageSource("CAL_ELXV_SPCOOR")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
  
      !c Rad -> deg
      LATM(:) = LATVE(1:3)*RAD2DEG ;

      !c Adjust so that 0 <= lon <= 360 
      LONM(:) = MODULO( LONVE(1:3)*RAD2DEG, 360.0D0 ) ; 
      LONTMP = LONM ; 

      CALL CAL_JAC( Jac1, LONM, LATM )
      CALL CAL_EDGELENGTH( DLED1, LONM, LATM ) ; 
   
      SPF = 0 ; 
      IF (  Jac1 < 0.0D0 .OR. 
     &     (Jac1 > 0.0D0 .AND. DLED1 > 360.0D0 ) ) THEN
         !c An element has a circular index. 
         !c Put a wrapped node on the other side

         !c A wrapped node is on the left of 180
         IF ( COUNT( LONM > 180.0D0 ) == 1 ) THEN
            IDX = sum( merge( (/ 1, 2, 3 /), 
     &           (/ 0, 0, 0 /), LONM > 180.0D0 ) ) ;

            LONM(IDX) = LONM(IDX) - 360.0D0 ; 
         END IF
         
         !c A wrapped node is on the right of 180
         IF ( COUNT( LONM < 180.0D0 ) == 1 ) THEN
            IDX = sum( merge( (/ 1, 2, 3 /), 
     &           (/ 0, 0, 0 /), LONM < 180.0D0 ) ) ;
            
            LONM(IDX) = LONM(IDX) + 360.0D0 ; 
         END IF        

         CALL CAL_JAC( Jac2, LONM, LATM )
         CALL CAL_EDGELENGTH( DLED2, LONM, LATM ) ; 

         SPF = 1 ;
         IF ( Jac2 < 0.0D0 .OR.
     &        (Jac1 > 0.0D0 .AND. DLED2 > DLED1) ) THEN
            LONM = LONTMP ;  
            SPF = 0 ; 
         END IF
      END IF
            
      !c deg -> rad
      LONM = LONM*DEG2RAD ;
      LATM = LATM*DEG2RAD ;
      DO II = 1, 3
         CALL CYLINDERMAP(XC, YC, LONM(II), LATM(II), SLAM0, SFEA0, ICS)

         XVE(II) = XC ;
         YVE(II) = YC ; 
      END DO
      
      IF ( PRESENT( SFLAG ) ) THEN
         SFLAG = SPF ; 
      END IF

#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") 
#endif
      call unsetMessageSource()
      
      RETURN ;

      CONTAINS

         SUBROUTINE CAL_EDGELENGTH( DEDGE, LON, LAT ) 
           IMPLICIT NONE

           REAL(8) :: DEDGE, LON(3), LAT(3)

           !c local c!
           REAL(8), dimension(3):: DLX, DLY

           DLX = (/ LON(2) - LON(1), 
     &          LON(3) - LON(2), LON(1) - LON(3) /) ; 
           DLY = (/ LAT(2) - LAT(1), 
     &          LAT(3) - LAT(2), LAT(1) - LAT(3) /) ; 
           
           DEDGE = SQRT( SUM(DLX*DLX + DLY*DLY) ) ; 

           RETURN ;
         END SUBROUTINE CAL_EDGELENGTH

         SUBROUTINE CAL_JAC( Jac, LON, LAT )
           IMPLICIT NONE

           REAL(8) :: Jac, LON(3), LAT(3)

           REAL(8) :: XR(2), XS(2)

           XR = 0.5D0*(/ LON(2) - LON(1), LAT(2) - LAT(1) /) ; 
           XS = 0.5D0*(/ LON(3) - LON(1), LAT(3) - LAT(1) /) ; 

           Jac = (XR(1)*XS(2) - XR(2)*XS(1)) ; ! Find Jacobian

           RETURN ;
         END SUBROUTINE CAL_JAC

      END SUBROUTINE CAL_ELXV_SPCOOR

      SUBROUTINE COMPUTE_CYLINPROJ_SFAC( SLAMV, SFEAV )
        IMPLICIT NONE

! dummy !
        REAL(8) :: SLAMV(:), SFEAV(:)
       
        ! local !
        INTEGER:: I, MEXP
        REAL(8) :: RFAC1, RFAC2, RFAC3

        call setMessageSource("COMPUTE_CYLINPROJ_SFAC")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif

!.... Set:
!.... SFCT, SFCX, SFCY, SFMX, SFMY, YCSFAC
!.... accordingly in order to determine the forms of equations to be solved
!     
!.....DEFAULT:
!     Cartesian coordinate 

        SELECT CASE( ICS )
        CASE (1)                  !c Cartesian
           SFAC = 1.0D0 ; 
        CASE (2)                !c Original SWE form in spherical coordinates
           DO I=1,NP
              SFAC(I)=COS(SFEA0)/COS(SFEAV(I)) ;
           END DO
         
           SFCX = SFAC ;
           SFMX = SFAC ;
        CASE (20,21,22)         ! Alternative SWE form in spherical coordinates
           MEXP = MOD(ICS,20) ;
!......ICS = 20 : Equal area projection
!......    = 21 : CPP
!.....     = 22 : Mercartor
!......from WP
           DO I = 1, NP
              SFCT(I) = cos(SFEAV(I))**MEXP ;           
              SFCX(I) = cos(SFEA0)*( cos(SFEAV(I))**(MEXP - 1) ) ;
              SFCY(I) = cos(SFEA0)**(MEXP - 1) ; 
            
              SFMX(I) = cos(SFEA0)/COS(SFEAV(I)) ;
              SFMY(I) = SFMX(I)**(MEXP - 1) ;

              TANPHI(I) = (TAN(SFEAV(I))/Rearth) ;
              YCSFAC(I) = cos( SFEAV(I) ) ; 
           END DO
        CASE (23)
!.......ICS = 23: Miller - compromised Mercator
           DO I = 1, NP
              SFCT(I) = cos(SFEAV(I))*COS(4.D0*SFEAV(I)/5.D0) ;       
              SFCX(I) = cos(SFEA0)*COS(4.D0*SFEAV(I)/5.D0) ;           
              SFCY(I) = cos(SFEA0) ;
              
              SFMX(I) = cos(SFEA0)/COS(SFEAV(I)) ;
              SFMY(I) = cos(SFEA0)/COS(4.D0*SFEAV(I)/5.D0) ;
              
              TANPHI(I) = (TAN(SFEAV(I))/Rearth) ;
              YCSFAC(I) = cos( SFEAV(I) ) ; 
           END DO
        CASE (24)
!......ICS = 24: Gall Stereographic projection           
           RFAC1 = 0.5D0*SQRT(2.D0) ;
           RFAC2 = 0.5D0*(1.D0 + 0.5D0*SQRT(2.D0)) ;

           DO I = 1, NP
              RFAC3 = cos(0.5D0*SFEAV(I))**2 ;
  
              SFCT(I) = cos(SFEAV(I))*RFAC3         
              SFCX(I) = RFAC1*RFAC3 ;            
              SFCY(I) = RFAC2 ;
              
              SFMX(I) = RFAC1/cos(SFEAV(I)) ;
              SFMY(I) = RFAC2/RFAC3 ;
              
              TANPHI(I) = (TAN(SFEAV(I))/Rearth) ;
              YCSFAC(I) = cos( SFEAV(I) ) ; 
           END DO
        END SELECT
        
!     COMPUTE ELEMENT AVERAGE FROM NODAL VECTORS ADJUSTING EQUATIONS
!     TO CYLINDER COORDINATES
        CALL SFAC_ELEAVG( SFacEle, SFAC, NM, NE ) ; 
        CALL SFAC_ELEAVG( SFCXEle, SFCX, NM, NE ) ; 
        CALL SFAC_ELEAVG( SFCYEle, SFCY, NM, NE ) ; 
        CALL SFAC_ELEAVG( SFMXEle, SFMX, NM, NE ) ; 
        CALL SFAC_ELEAVG( SFMYEle, SFMY, NM, NE ) ; 
        CALL SFAC_ELEAVG( YCSFacEle, YCSFAC, NM, NE ) ;          
        CALL SFAC_ELEAVG( TANPHIEle, TANPHI, NM, NE ) ;          
        
!      Only need the  averages of these... 
        DEALLOCATE( SFMY, SFCX, SFCY, YCSFAC )

#if defined(MESH_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.") 
#endif
        call unsetMessageSource()
        
        RETURN ; 
      END SUBROUTINE COMPUTE_CYLINPROJ_SFAC

!-----------------------------------------------------------------------
      function findIBTYPEAtNode(I) RESULT(IBTYPE)
!-----------------------------------------------------------------------
      use boundaries,only : LBCODEI
      implicit none
      integer,intent(in) :: i
      integer :: j,ibtype
      ibtype = -1
      j = lbarray_pointer(i)
      if (j > 0) then
         ibtype = lbcodei(j)
      end if
      end function findIBTYPEAtNode
      
!     ------------------------------------------------------------------
      end module mesh
!     ------------------------------------------------------------------
