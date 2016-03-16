!-------------------------------------------------------------------------
! adcirc2xdmf.f90
!-------------------------------------------------------------------------
! Author: Jason Fleming (jason.fleming@seahorsecoastal.com) 
!
! Convert ascii adcirc data file to XDMF format. 
!---------------------------------------------------------------------
!---------------------------------------------------------------------
program adcirc2xdmf
!---------------------------------------------------------------------
!---------------------------------------------------------------------
use netcdf
use adcmesh
use asgsio
use nodalattr
use control
implicit none
include 'Xdmf.f' 
! 
character(1024) :: dataFileName ! full path name of the ascii file to be converted
character(1024) :: controlFileName  ! full path name of the adcirc fort.15 for metadata  
character(1024) :: convertedFileName  ! full path name of the converted files  
integer*8 :: xdmfFortranObj  ! represents pointer to the XdmfFortran object
integer :: informationID     ! information object reference
integer :: topologyID        ! topology reference 
integer :: geometryID        ! geometry reference
integer :: depthID
integer :: attributeID       ! attribute reference
integer :: attributeType     ! attribute type
integer, allocatable :: xdmf_nm(:,:)   ! 0-offset connectivity array
logical :: convertOutputData ! .true. if output data conversion is required
logical :: release           ! .true. to release memory after writing data to HDF5
logical :: writeToHDF5       ! .true. if XdmfAddGrid should immediately write mesh
integer :: unitNumber        ! fortran i/o unit number for ascii data file 
logical :: formatKnown       ! .true. if the ascii format has been discovered already
logical :: sparseAsciiFile   ! .true. if the ascii file has sparse format
real(8) :: defaultValue      ! fill value for sparse ascii format
integer :: lightDataLimit    ! max values to save as light data
integer :: tsinterval          ! spooling interval in time steps
real(8) :: tinterval          ! spooling interval in seconds
integer :: numSnaps          ! unreliable number of datasets in ascii file
integer :: numValues         ! total number of values in a dataset
real(8) :: timeSec            ! time of a dataset, in seconds
integer :: timeStep          ! time step number for a dataset
integer :: ss                ! dataset counter
integer :: nCol              ! number of columns of data in ascii file
integer :: numNodes          ! number of nodes in mesh according to ascii file
integer :: numNodesNonDefault ! in sparse format, the number of nodes in a dataset that do not have the default value
integer :: timeStepID
character(len=256) :: timeStepString
integer :: startingDataset  ! first dataset to convert
integer :: endingDataset    ! last dataset to convert
logical :: meshonly ! .true. if only the mesh is to be converted
!
! levees and boundaries for visualization
integer :: leveeDimensions(3)
integer :: leveeDimensionsID
integer :: leveeGeometryID
integer :: boundaryID
integer, allocatable :: ibtypeGeom(:)
integer :: b
integer :: numCoord
integer :: ind
integer :: nodeNum
real(8) :: leveeHeight
character(len=256) :: boundaryName 
!
character(len=256) :: projection
!
type(xdmfMetaData_t) :: md  ! holds the metadata for whatever data we are writing
!
real(8), allocatable :: data_array1(:)
real(8), allocatable :: data_array3(:,:)
!
integer :: lastSlashPosition ! used for trimming full path from a filename
integer :: lastDotPosition ! to determine file extension
character(2048) :: dataFileExtension ! something like 13, 14, 15, 63, 222 etc

character(80) :: line ! a line of data from the ascii file
character(80) :: topcomment ! comment line at the top of ascii file
integer :: i, j, n 
!
! initializations
startingDataset = 0
endingDataset = 999999
datafileName = 'none'
meshFileName = 'fort.14'
controlFileName = 'fort.15'
projection = 'unknown'
convertOutputData = .true.
verbose = .false.
release = .true.
writeToHDF5 = .true.
meshonly = .false.
lightDataLimit = 10
unitNumber = 90
ss=1  ! initialize the dataset counter
sparseAsciiFile = .false.
formatKnown = .false.
i=0
!
! process command line options
argcount = iargc() ! count up command line options
write(6,'("INFO: There are ",i0," command line arguments.")') argcount
do while (i.lt.argcount)
   i = i + 1
   call getarg(i, cmdlineopt)
   select case(trim(cmdlineopt))
   case("--verbose")
      write(6,'(a)') "INFO: Processing " // trim(cmdlineopt) // "."
      verbose = .true.
   case("--nodalattributes")
      write(6,'(a)') "INFO: Processing " // trim(cmdlineopt) // "."
      convertOutputData = .false.
   case("--meshfile")
      i = i + 1
      call getarg(i, cmdlinearg)
      write(6,'(a)') "INFO: Processing " // trim(cmdlineopt) // & 
         " " // trim(cmdlinearg) // "."
      meshFileName = trim(cmdlinearg)
   case("--datafile")
      i = i + 1
      call getarg(i, dataFileName)     
      write(6,'(a)') "INFO: Processing " // trim(cmdlineopt) // " " // &
         trim(dataFileName) // "."
   case("--controlfile")
      i = i + 1
      call getarg(i, controlFileName)
      write(6,'(a)') "INFO: Processing " // trim(cmdlineopt)  // " " // &
         trim(controlFileName) // "."
   case("--lightdatalimit")
      i = i + 1
      call getarg(i, cmdlinearg)
      write(6,'(a)') "INFO: Processing " // trim(cmdlineopt) // " " // & 
         trim(cmdlinearg) // "."
      read(cmdlinearg,*) lightdatalimit
   case("--starting-dataset")
      i = i + 1
      call getarg(i, cmdlinearg)
      write(6,'(a)') "INFO: Processing " // trim(cmdlineopt) // " " // & 
         trim(cmdlinearg) // "."
      read(cmdlinearg,*) startingDataset
   case("--ending-dataset")
      i = i + 1
      call getarg(i, cmdlinearg)
      write(6,'(a)') "INFO: Processing " // trim(cmdlineopt) // " " // & 
         trim(cmdlinearg) // "."
      read(cmdlinearg,*) endingDataset
   case default
      write(6,'(a)') "WARNING: Command line option '" // &
         TRIM(cmdlineopt) // "' was not recognized."
   end select
end do
! read ADCIRC mesh from ascii file
write(6,'(a,a)') 'INFO: Preparing to read mesh file "'//trim(meshFileName),'".'
call read14()
call constructFluxBoundaryTypesArray() ! create LBCODEI array
!
! read ADCIRC control data from ascii file unless otherwise specified
if (trim(adjustl(controlFileName)).ne."none") then
   write(6,'(a)') 'INFO: Reading control file.'
   call readControlFile(controlFileName,verbose)
   if (verbose.eqv..true.) then
      write(6,'(a)') 'INFO: Echoing control file.'
      open(25,file='echo.15',status='replace')
      call echoControlFile(25)
      close(25)
      write(6,'(a)') 'INFO: Finished echoing control file.'!
   endif
   ! modify projection string if indicated by the fort.15
   select case(ics)
   case(1)
      projection = 'cartesian' ! most likely CPP, but no way of knowing for sure
   case(2)
      projection = 'geographic'
   case default
      projection = 'unknown'
   end select
endif
!
! create a 0-offset element table, i.e., one that assumes the nodes
! are numbered starting from zero ... XDMF2 is zero offset while the 
! ADCIRC ascii format assumes nodes are numbered starting from 1
allocate(xdmf_nm(3,ne))
do i=1,ne
   do j=1,3
      xdmf_nm(j,i) = nm(i,j) - 1 ! store in row major format
   end do
end do
! 
! create xdmf file
!
! call the XdmfFortran constructor; return a reference to the XdmfFortran object
! (the pointer is packaged as a long int in Fortran)
write(6,'(A)') 'INFO: Initializing XDMF object.'
call XdmfInit(xdmfFortranObj)
! 
! call the initHDF5 method; arguments include the ref to the XdmfFortran
! object, the name of the HDF5 file, and whether to release memory after write 
write(6,'(A)') 'INFO: Initializing HDF5 file.'
!
! read ADCIRC control data from ascii file unless otherwise specified
if (trim(adjustl(dataFileName)).eq."none") then
   meshonly = .true.
   convertOutputData = .false.
endif
!
! If the full path to the file was given, the path to the file must 
! trimmed off, leaving just the file name. This way the converted files
! will be written to the directory where the command was executed rather
! than the directory where the file to be converted is located (if these
! locations are different.
lastSlashPosition = index(trim(dataFileName),"/",.true.)
if (meshonly.eqv..true.) then
   lastSlashPosition = index(trim(meshFileName),"/",.true.)
   convertedFileName = trim(meshFileName(lastSlashPosition+1:))
   convertOutputData = .false.
   dataFileExtension = 'none'
else
   lastSlashPosition = index(trim(dataFileName),"/",.true.)
   convertedFileName = trim(dataFileName(lastSlashPosition+1:))
   lastDotPosition = index(trim(convertedFileName),'.',.true.)
   dataFileExtension = trim(convertedFileName(lastDotPosition+1:))
   !write(6,'("DEBUG: Data file extension is ",a,".")') trim(dataFileExtension) !jgfdebug
   !
   ! Check to see if the data file is a nodal attributes file (with a
   ! .13 extension)
   if (trim(dataFileExtension).eq.'13') then
      write(6,'(a)') 'INFO: The data file to be converted is a nodal attributes input file.'
      convertOutputData = .false.
   endif
endif
!
!
write(6,'(a,a,a,a,a)') 'DEBUG: The names of the converted files will be ',trim(convertedFileName),'.h5 and ',trim(convertedFileName),'.xmf.'
call XdmfInitHDF5(xdmfFortranObj, trim(convertedFileName)//'.h5'//char(0), release)
!
! name the temporal grid collection for time varying output data
if (convertOutputData.eqv..true.) then
   informationID = XdmfAddInformation(xdmfFortranObj, 'TimeVaryingOutputData'//CHAR(0), trim(adjustl(agrid))//CHAR(0))
endif
!
! write fort.15 data to the xml file unless otherwise specified
if (trim(adjustl(controlFileName)).ne."none") then
   call writeControlXDMF(xdmfFortranObj)
endif
!
! create a grid collection for time varying data 
if (convertOutputData.eqv..true.) then
   call xdmfAddGridCollection(xdmfFortranObj, "Temporal"//CHAR(0), &
       XDMF_GRID_COLLECTION_TYPE_TEMPORAL)
endif
!
! call the setTopology method; arguments include the type of the topology, 
! the number of values in the connectivity array, the numeric type of
! the connectivity array, and the connectivity values with a 0 offset;
! the call returns the topology ID in case it will be reused. The XDMF
! library interprets the connectivity array as a series of contiguous 
! values (i.e., as a 1D array). 
topologyID = xdmfSetTopology(xdmfFortranObj, XDMF_TOPOLOGY_TYPE_TRIANGLE, ne*3, XDMF_ARRAY_TYPE_INT32, xdmf_nm)
!
!
! call the setGeometry method; arguments include the geometry type, the 
! number of values in the coordinates array, the numeric type of the coordinates,
! and the coordinates array (again, interpreted by XDMF as a series of 
! contiguous values). 
geometryID = xdmfSetGeometry(xdmfFortranObj, XDMF_GEOMETRY_TYPE_XY, np*2, XDMF_ARRAY_TYPE_FLOAT64, xyd(1:2,:))
!
! add the mesh boundaries to the grid
call addBoundaries(xdmfFortranObj)
!
! add bathymetric depth
write(6,'(a)') 'INFO: Writing bathy/topo to XDMF file.'
! write the bathymetric depth
md%variable_name = 'depth'
md%long_name = 'distance from geoid'
md%standard_name = 'depth_below_geoid'
md%coordinates = 'y x'
md%units = 'm'
md%positive = 'downward'
call writeMetaData(xdmfFortranObj, md)
! write projection info 
md%variable_name_id = XdmfAddInformation(xdmfFortranObj, 'projection'//CHAR(0), &
   trim(projection)//CHAR(0))
depthID = XdmfAddAttribute(xdmfFortranObj, trim(md%variable_name)//CHAR(0), &
   XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, np, &
   XDMF_ARRAY_TYPE_FLOAT64, xyd(3,:))
write(6,'(a)') 'INFO: Finished writing bathy/topo to XDMF file.'
!
! read the data from the adcirc data file and write datasets to xdmf
md%createdIDs = .false.
md%positive = 'null'
md%ndset = 0 ! unknown number of data sets
select case(trim(dataFileExtension))
case('63','69')
   if (trim(dataFileName).eq.'maxele.63') then
      md%variable_name = 'zeta_max'
      md%long_name = 'max water surface elevation above geoid'
      md%standard_name = 'max_water_surface_elevation'
      md%coordinates = 'y x'
      md%units = 'm'   
      md%ndset = 1 ! TODO: fix to also handle min/max files with time of occurrence
   else 
      md%variable_name = 'zeta'
      md%long_name = 'water surface elevation above geoid'
      md%standard_name = 'water_surface_elevation'
      md%coordinates = 'time y x'
      md%units = 'm'
   endif
case('64')
   md%variable_name = 'vel'
   md%long_name = 'water column vertically averaged velocity'
   md%standard_name = 'water_velocity'
   md%coordinates = 'time y x'
   md%units = 'm s-1'
   md%positive = 'north/east'
case('73')
   md%variable_name = 'pressure'
   md%long_name = 'air pressure at sea level'
   md%standard_name = 'air_pressure_at_sea_level'
   md%coordinates = 'time y x'
   md%units = 'meters of water'
case('74')
   md%variable_name = 'wvel'
   md%long_name = 'wind_velocity'
   md%standard_name = 'wind'
   md%coordinates = 'time x y'
   md%units = 'm s-1'
   md%positive = 'north/east'
case('93')
   md%variable_name = 'ice'
   md%long_name = 'ice coverage at at sea surface'
   md%standard_name = 'ice_pressure_at_sea_level'
   md%coordinates = 'time y x'
   md%units = 'percent'
case('none')
   ! do nothing, we are just converting a mesh file
case('13')
   call readNodalAttributesFile(dataFileName)
   call writeNodalAttributesXDMF(xdmfFortranObj)
case default  ! includes fort.13 and none
   write(6,'(A)') 'ERROR: adcirc2xdmf cannot convert ' // trim(dataFileName) // ' files.'
   stop
end select
!
! if we are writing plain old ADCIRC data
if (convertOutputData.eqv..true.) then
   call openFileForRead(UnitNumber, trim(dataFileName))
   read(unitNumber,'(A80)') topcomment ! comment line at the top
   ! jgf: Can't rely on the NumSnaps value; in general, it will not
   ! actually reflect the number of datasets in the file.
   read(unitNumber,*) numSnaps, numNodes, tinterval, tsinterval, nCol
   if (np.ne.numNodes) then
      write(6,*) 'ERROR: The output file contains ',NumNodes,        &
         ' nodes, but the mesh file contains ',np,' nodes.'
      write(6,*) 'ERROR: The output file does not correspond to the mesh file.'
      close(unitNumber)
      stop
   endif
   !
   select case(nCol)
   case(1) ! scalar data
      allocate(data_array1(1:numNodes))
      attributeType = XDMF_ATTRIBUTE_TYPE_SCALAR
      numValues = numNodes
   case(2) ! vector data
      ! XDMF doesn't understand vectors with two components, so we
      ! have to allocate for three components, and set the z component
      ! of the vector to zero
      allocate(data_array3(3,1:numNodes))
      attributeType = XDMF_ATTRIBUTE_TYPE_VECTOR
      numValues = 3*numNodes
   case default
      write(6,*) 'ERROR: The ADCIRC output file contains ',nCol,        &
         ' columns, but adcirc2xdmf only supports 1 or 2 column data.'
      write(6,*) 'ERROR: adcirc2xdmf cannot continue.'
      close(unitNumber)
      stop   
   end select
   !
   write(6,'(A)') 'INFO: Adding unstuctured mesh dataset.'
   do   ! jgf: loop until we run out of data
      read(unitNumber,'(A)',END=2,ERR=2) line
      ! if this is the first dataset, we will need to discover the format
      ! of the file: sparse ascii, or full ascii
      if (formatKnown.eqv..false.) then
         ! try to read extra info that are only present in sparse ascii files
         read(line,*,END=1, ERR=1) timeSec, timeStep, numNodesNonDefault, defaultValue
         sparseAsciiFile = .true. ! we only get to this line if the extra values were present
   1     if ( sparseAsciiFile .eqv..false. ) then
            write(6,'(A)') 'INFO: The ascii file is not in sparse format.'
         endif
         formatKnown = .true. ! skip the format discovery for remaining datasets
      endif
      !
      ! we know the format of the file, now use the appropriate read statement
      if (sparseAsciiFile.eqv..false.) then
         read(line,*) timeSec, timeStep ! dataset time in seconds, time step number      
         numNodesNonDefault = numNodes
         defaultValue = -99999.0d0
      else
         read(line,*) timeSec, timeStep, numNodesNonDefault, defaultValue    
      endif
      ! 
      ! now read the data from the ascii ADCIRC file
      select case(nCol)
      case(1) ! scalar data
         data_array1 = defaultValue
         do n=1,numNodesNonDefault
            read(unitNumber,*) j, data_array1(j)
         end do
      case(2) ! 2D vector data
         data_array3 = defaultValue
         do n=1,numNodesNonDefault
            read(unitNumber,*) j, data_array3(1,j), data_array3(2,j)
         end do
         data_array3(3,:) = 0.d0
      case default
         ! already accounted for prior to reading the data
      end select
      !
      ! skip to the next dataset if this dataset was not specified by the user
      if (ss.lt.startingDataset) then
         write(6,fmt='(I4)',advance='no') ss
         ss = ss + 1
         cycle
      endif
      if (ss.gt.endingDataset) then
         exit
      endif
      !
      ! add boundary references if necessary
      if (ss.gt.1) then
         call addBoundaryReferences(xdmfFortranObj)
         call xdmfAddPreviousAttribute(xdmfFortranObj, depthID)
      endif
      !
      ! set the time of this dataset in the XDMF file
      call XdmfSetTime(xdmfFortranObj, timeSec)
      !
      ! set the metadata for this dataset
      call writeMetaData(xdmfFortranObj, md)      
      ! 
      ! now read the data from the ascii ADCIRC file and write it to 
      ! the xdmf file according to the data type (scalar or vector)
      select case(nCol)
      case(1) ! scalar data
         attributeID = XdmfAddAttribute(xdmfFortranObj, trim(md%variable_name)//CHAR(0), &
            XDMF_ATTRIBUTE_CENTER_NODE, attributeType, numValues, &
            XDMF_ARRAY_TYPE_FLOAT64, data_array1)
      case(2) ! 2D vector data
         attributeID = XdmfAddAttribute(xdmfFortranObj, trim(md%variable_name)//CHAR(0), &
            XDMF_ATTRIBUTE_CENTER_NODE, attributeType, numValues, &
            XDMF_ARRAY_TYPE_FLOAT64, data_array3)
      case default
         ! already accounted for prior to reading the data
      end select
      !
      ! set the time step of this dataset in the XDMF file
      write(timeStepString,'(i0)') timeStep
      timeStepID = XdmfAddInformation(xdmfFortranObj, 'IT'//CHAR(0), &
      trim(timeStepString)//CHAR(0))
      ! 
      ! call the addGrid method; creates an unstructured mesh object with the
      ! specified name (2nd arg), then associates the geometry and topology 
      ! created above with this new unstructured mesh, also associates any
      ! informations or attributes with the new mesh, immediately writing
      ! it to the hdf5 file if the last argument is set to .true. 
      write(6,fmt='(I4)',advance='no') ss
      call XdmfAddGrid(xdmfFortranObj,trim(agrid)//char(0), writeToHDF5)
      ss = ss + 1
      !
      !
      ! quit the loop if we are only supposed to write a certain number
      ! of datasets
      if ( (md%ndset.ne.0).and.(ss.gt.md%ndset)) then
         exit
      endif
   end do
   2  close(unitNumber) ! jump to here when all data sets have been read.
   !
   ! close this grid collection, writing it to the heavy data (HDF5) file
   ! if the value of writeToHDF5 is .true.; further grids cannot 
   ! be added to this collection
   call XdmfCloseGridCollection(xdmfFortranObj, writeToHDF5)
else
   !
   ! if we are not writing time-dependent output data, simply add this
   ! grid to the file
   call XdmfAddGrid(xdmfFortranObj,trim(agrid)//char(0), writeToHDF5)
endif
!
! call the write method; arguments include the name of the xml file, 
! the max number of light data values, and whether to release memory
! after writing. This actually writes both the light and heavy data.
write(6,'(/,A)') 'INFO: Writing data to file and releasing memory.'
call XdmfWrite(xdmfFortranObj, trim(convertedFileName)//'.xmf'//char(0), lightDataLimit, release)
!
! call the close method (deletes the XdmfFortran object)
write(6,'(A)') 'INFO: Cleaning up and deleting XDMF object.'
call XdmfClose(xdmfFortranObj)
!
!
! L E V E E S   A N D   B O U N D A R I E S
!
! for visualization
write(6,'(A)') 'INFO: Initializing XDMF object.'
call XdmfInit(xdmfFortranObj)
! 
! call the initHDF5 method; arguments include the ref to the XdmfFortran
! object, the name of the HDF5 file, and whether to release memory after write 
write(6,'(A)') 'INFO: Initializing HDF5 file.'
if (trim(adjustl(dataFileName)).eq.'none') then
   convertedFileName = trim(adjustl(meshFileName))//'_viz_boundaries'
else
   convertedFileName = trim(adjustl(dataFileName))//'_viz_boundaries'
endif
!
call XdmfInitHDF5(xdmfFortranObj, trim(convertedFileName)//'.h5'//char(0), release)
!
! create a grid collection for visualizing spatial data     
call xdmfAddGridCollection(xdmfFortranObj, "Levees and Boundaries"//CHAR(0), &
    XDMF_GRID_COLLECTION_TYPE_SPATIAL)   


!
! simple flux boundaries (mainland, island, river)
do b = 1, numSimpleFluxBoundaries
   n = size(simpleFluxBoundaries(b)%nodes)
   leveeDimensions(1) = 2 ! front base, top front face (no back face)
   leveeDimensions(2) = n ! length of boundary
   leveeDimensions(3) = 1 ! has no thickness
   numCoord = 3*product(leveeDimensions)
   allocate(simpleFluxBoundaries(b)%bGeom(numCoord))
   simpleFluxBoundaries(b)%bGeom = 0.0d0
   ind = 1
   do i=1,n
      nodeNum = simpleFluxBoundaries(b)%nodes(i) + 1
      simpleFluxBoundaries(b)%bGeom(ind)   = xyd(1,nodeNum) ! x 
      simpleFluxBoundaries(b)%bGeom(ind+1) = xyd(2,nodeNum) ! y
      simpleFluxBoundaries(b)%bGeom(ind+2) = - xyd(3,nodeNum) ! depth
      ind = ind + 3
   end do
   do i=1,n
      ! height
      nodeNum = simpleFluxBoundaries(b)%nodes(i) + 1
      simpleFluxBoundaries(b)%bGeom(ind)   = xyd(1,nodeNum) ! x 
      simpleFluxBoundaries(b)%bGeom(ind+1) = xyd(2,nodeNum) ! y
      simpleFluxBoundaries(b)%bGeom(ind+2) = 1.d0 ! hard code to 1m arbitrarily 
      ind = ind + 3
   end do
   leveeDimensionsID = xdmfSetDimensions(xdmfFortranObj, 3, XDMF_ARRAY_TYPE_INT32, leveeDimensions)
   leveeGeometryID = xdmfSetGeometry(xdmfFortranObj, XDMF_GEOMETRY_TYPE_XYZ, numCoord, &
      XDMF_ARRAY_TYPE_FLOAT64, simpleFluxBoundaries(b)%bGeom)

   !allocate(ibtypeGeom(2*product(leveeDimensions)))
   !ibtypeGeom(:) = ibtype_orig(simpleFluxBoundaries(b)%indexNum)  
   !boundaryID = xdmfAddAttribute(xdmfFortranObj, 'ibtype'//CHAR(0), &
   !   XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, 2*product(leveeDimensions), &
   !   XDMF_ARRAY_TYPE_INT32, ibtypeGeom ) 
   !deallocate(ibtypeGeom)
          
   write(boundaryName,'("fluxBoundary",i0,"_ibtype",i0)') &
      simpleFluxBoundaries(b)%indexNum, ibtype_orig(simpleFluxBoundaries(b)%indexNum)  
      
   call xdmfAddGridCurvilinear(xdmfFortranObj,trim(boundaryName)//char(0), writeToHDF5)
end do


!
! external flux boundaries (overflow out of the domain)
do b = 1, numExternalFluxBoundaries
   n = size(externalFluxBoundaries(b)%nodes)
   leveeDimensions(1) = 2 ! front base, top front face (no back face)
   leveeDimensions(2) = n ! length of boundary
   leveeDimensions(3) = 1 ! has no thickness
   numCoord = 3*product(leveeDimensions)
   allocate(externalFluxBoundaries(b)%leveeGeom(numCoord))
   externalFluxBoundaries(b)%leveeGeom = 0.0d0
   ind = 1
   do i=1,n
      nodeNum = externalFluxBoundaries(b)%nodes(i) + 1
      externalFluxBoundaries(b)%leveeGeom(ind)   = xyd(1,nodeNum) ! x 
      externalFluxBoundaries(b)%leveeGeom(ind+1) = xyd(2,nodeNum) ! y
      externalFluxBoundaries(b)%leveeGeom(ind+2) = - xyd(3,nodeNum) ! depth
      ind = ind + 3
   end do
   do i=1,n
      ! height
      leveeHeight = externalFluxBoundaries(b)%barlanht(i)
      nodeNum = externalFluxBoundaries(b)%nodes(i) + 1
      externalFluxBoundaries(b)%leveeGeom(ind)   = xyd(1,nodeNum) ! x 
      externalFluxBoundaries(b)%leveeGeom(ind+1) = xyd(2,nodeNum) ! y
      externalFluxBoundaries(b)%leveeGeom(ind+2) = leveeHeight 
      ind = ind + 3
   end do
   leveeDimensionsID = xdmfSetDimensions(xdmfFortranObj, 3, XDMF_ARRAY_TYPE_INT32, leveeDimensions)
   leveeGeometryID = xdmfSetGeometry(xdmfFortranObj, XDMF_GEOMETRY_TYPE_XYZ, numCoord, &
      XDMF_ARRAY_TYPE_FLOAT64, externalFluxBoundaries(b)%leveeGeom)    
   write(boundaryName,'("fluxBoundary",i0,"_ibtype",i0)') &
      externalFluxBoundaries(b)%indexNum, ibtype_orig(externalFluxBoundaries(b)%indexNum)  
   call xdmfAddGridCurvilinear(xdmfFortranObj,trim(boundaryName)//char(0), writeToHDF5)
end do
!
! levees (internal flux boundaries)
do b = 1, numInternalFluxBoundaries
   n = size(internalFluxBoundaries(b)%nodes)
   leveeDimensions(1) = 4 ! front base, top front face, top back face, base back face
   leveeDimensions(2) = n ! length of boundary
   leveeDimensions(3) = 1 ! has no thickness
   numCoord = 3*product(leveeDimensions)
   allocate(internalFluxBoundaries(b)%ibTypeAttribute(3*product(leveeDimensions)))
   internalFluxBoundaries(b)%ibTypeAttribute = ibtype_orig(internalFluxBoundaries(b)%indexNum)
   allocate(internalFluxBoundaries(b)%leveeGeom(numCoord))
   internalFluxBoundaries(b)%leveeGeom = 0.0d0
   ind = 1
   ! front side
   do i=1,n
      nodeNum = internalFluxBoundaries(b)%nodes(i) + 1
      internalFluxBoundaries(b)%leveeGeom(ind)   = xyd(1,nodeNum) ! x 
      internalFluxBoundaries(b)%leveeGeom(ind+1) = xyd(2,nodeNum) ! y
      internalFluxBoundaries(b)%leveeGeom(ind+2) = - xyd(3,nodeNum) ! depth
      ind = ind + 3
   end do
   do i=1,n
      ! height
      nodeNum = internalFluxBoundaries(b)%nodes(i) + 1
      internalFluxBoundaries(b)%leveeGeom(ind)   = xyd(1,nodeNum) ! x 
      internalFluxBoundaries(b)%leveeGeom(ind+1) = xyd(2,nodeNum) ! y
      internalFluxBoundaries(b)%leveeGeom(ind+2) = internalFluxBoundaries(b)%barinht(i) 
      ind = ind + 3
   end do
   ! back face
   do i=1,n
      nodeNum = internalFluxBoundaries(b)%ibconn(i) + 1
      internalFluxBoundaries(b)%leveeGeom(ind)   = xyd(1,nodeNum) ! x 
      internalFluxBoundaries(b)%leveeGeom(ind+1) = xyd(2,nodeNum) ! y
      internalFluxBoundaries(b)%leveeGeom(ind+2) = internalFluxBoundaries(b)%barinht(i) 
      ind = ind + 3
   end do
   do i=1,n
      nodeNum = internalFluxBoundaries(b)%ibconn(i) + 1
      internalFluxBoundaries(b)%leveeGeom(ind)   = xyd(1,nodeNum) ! x 
      internalFluxBoundaries(b)%leveeGeom(ind+1) = xyd(2,nodeNum) ! y
      internalFluxBoundaries(b)%leveeGeom(ind+2) = - xyd(3,nodeNum) ! depth
      ind = ind + 3
   end do

   leveeDimensionsID = xdmfSetDimensions(xdmfFortranObj, 3, XDMF_ARRAY_TYPE_INT32, leveeDimensions)
   leveeGeometryID = xdmfSetGeometry(xdmfFortranObj, XDMF_GEOMETRY_TYPE_XYZ, numCoord, &
      XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundaries(b)%leveeGeom)    
   !attributeID = XdmfAddAttribute(xdmfFortranObj, 'ibtype'//CHAR(0), &
   !      XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, &
   !      numCoord, XDMF_ARRAY_TYPE_INT32, internalFluxBoundaries(b)%ibtypeAttribute)   
   write(boundaryName,'("fluxBoundary",i0,"_ibtype",i0)') &
      internalFluxBoundaries(b)%indexNum, ibtype_orig(internalFluxBoundaries(b)%indexNum)  
   call xdmfAddGridCurvilinear(xdmfFortranObj,trim(boundaryName)//char(0), writeToHDF5)
end do
!
! levees with culverts (internal barrier boundaries with cross barrier pipes)
do b = 1, numInternalFluxBoundariesWithPipes
   n = size(internalFluxBoundariesWithPipes(b)%nodes)
   leveeDimensions(1) = 4 ! front base, top front face, top back face, base back face
   leveeDimensions(2) = n ! length of boundary
   leveeDimensions(3) = 1 ! has no thickness
   numCoord = 3*product(leveeDimensions)
   allocate(internalFluxBoundariesWithPipes(b)%leveeGeom(numCoord))
   internalFluxBoundariesWithPipes(b)%leveeGeom = 0.0d0
   ind = 1
   ! front side
   do i=1,n
      nodeNum = internalFluxBoundariesWithPipes(b)%nodes(i) + 1
      internalFluxBoundariesWithPipes(b)%leveeGeom(ind)   = xyd(1,nodeNum) ! x 
      internalFluxBoundariesWithPipes(b)%leveeGeom(ind+1) = xyd(2,nodeNum) ! y
      internalFluxBoundariesWithPipes(b)%leveeGeom(ind+2) = - xyd(3,nodeNum) ! depth
      ind = ind + 3
   end do
   do i=1,n
      ! height
      nodeNum = internalFluxBoundariesWithPipes(b)%nodes(i) + 1
      internalFluxBoundariesWithPipes(b)%leveeGeom(ind)   = xyd(1,nodeNum) ! x 
      internalFluxBoundariesWithPipes(b)%leveeGeom(ind+1) = xyd(2,nodeNum) ! y
      internalFluxBoundariesWithPipes(b)%leveeGeom(ind+2) = internalFluxBoundariesWithPipes(b)%barinht(i) 
      ind = ind + 3
   end do
   ! back face
   do i=1,n
      nodeNum = internalFluxBoundariesWithPipes(b)%ibconn(i) + 1
      internalFluxBoundariesWithPipes(b)%leveeGeom(ind)   = xyd(1,nodeNum) ! x 
      internalFluxBoundariesWithPipes(b)%leveeGeom(ind+1) = xyd(2,nodeNum) ! y
      internalFluxBoundariesWithPipes(b)%leveeGeom(ind+2) = internalFluxBoundariesWithPipes(b)%barinht(i)  
      ind = ind + 3
   end do
   do i=1,n
      nodeNum = internalFluxBoundariesWithPipes(b)%ibconn(i) + 1
      internalFluxBoundariesWithPipes(b)%leveeGeom(ind)   = xyd(1,nodeNum) ! x 
      internalFluxBoundariesWithPipes(b)%leveeGeom(ind+1) = xyd(2,nodeNum) ! y
      internalFluxBoundariesWithPipes(b)%leveeGeom(ind+2) = - xyd(3,nodeNum) ! depth
      ind = ind + 3
   end do

   leveeDimensionsID = xdmfSetDimensions(xdmfFortranObj, 3, XDMF_ARRAY_TYPE_INT32, leveeDimensions)
   leveeGeometryID = xdmfSetGeometry(xdmfFortranObj, XDMF_GEOMETRY_TYPE_XYZ, numCoord, &
      XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundariesWithPipes(b)%leveeGeom)    
   write(boundaryName,'("fluxBoundary",i0,"_ibtype",i0)') &
      internalFluxBoundariesWithPipes(b)%indexNum, ibtype_orig(internalFluxBoundariesWithPipes(b)%indexNum)  
   call xdmfAddGridCurvilinear(xdmfFortranObj,trim(boundaryName)//char(0), writeToHDF5)
end do
!
call xdmfCloseGridCollection(xdmfFortranObj, writeToHDF5)
!
write(6,'(/,A)') 'INFO: Writing data to file and releasing memory.'
call XdmfWrite(xdmfFortranObj, trim(convertedFileName)//'.xmf'//char(0), lightDataLimit, release)
!
write(6,'(A)') 'INFO: Cleaning up and deleting XDMF object.'
call XdmfClose(xdmfFortranObj)

!----------------------------------------------------------------------
!----------------------------------------------------------------------
END PROGRAM adcirc2xdmf
!----------------------------------------------------------------------
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!                        S U B R O U T I N E    
!       W R I T E   N O D A L   A T T R I B U T E S   X D M F
!-----------------------------------------------------------------------
! jgf: Writes all nodal attributes to an XDMF file.
!-----------------------------------------------------------------------
subroutine writeNodalAttributesXDMF(xdmfFortranObj)
use nodalattr
implicit none
include 'Xdmf.f'
integer*8, intent(in) :: xdmfFortranObj ! object that receives the data

integer :: attributeID
integer :: attributeType
integer :: nameID
integer :: informationID
integer :: numValsID
integer :: defaultValsID
character(len=1024) :: defaultValsString
integer :: numValues
character(len=256) :: numValsString
character(len=256) :: numMeshNodesString
integer :: unitsID
integer :: i, j
integer :: infoCount
!
do i=1,numNodalAttributes
   ! allocate space
   numValues = na(i)%numVals*numMeshNodes
   select case(na(i)%numVals)
   case(1) ! scalar data
      allocate(na(i)%xdmfArray(1:numMeshNodes))
      attributeType = XDMF_ATTRIBUTE_TYPE_SCALAR
   case(2:) ! matrix data
      allocate(na(i)%xdmfMatrix(na(i)%numVals,1:numMeshNodes))
      attributeType = XDMF_ATTRIBUTE_TYPE_MATRIX
   case default
      write(6,*) 'ERROR: The nodal attribute "',trim(adjustl(na(i)%attrName)), &
         '" has ',na(i)%numVals,' values at each node, but adcirc2xdmf ', &
         'does not recognize this nodal attribute.'
      stop   
   end select
   ! 
   ! now read the data from the ascii ADCIRC file and write it to 
   ! the xdmf file according to the data type (scalar or vector)
   select case(na(i)%numVals)
   case(1) ! scalar data
      ! populate the nodal attribute array at every node
      na(i)%xdmfArray(:) = na(i)%defaultVals(1)
      do j=1,na(i)%numNodesNotDefault
         na(i)%xdmfArray(na(i)%nonDefaultNodes(j)) = na(i)%nonDefaultVals(1,j)
      end do
      attributeID = XdmfAddAttribute(xdmfFortranObj, trim(adjustl(na(i)%attrName))//CHAR(0), &
         XDMF_ATTRIBUTE_CENTER_NODE, attributeType, numValues, &
         XDMF_ARRAY_TYPE_FLOAT64, na(i)%xdmfArray)
   case(2:) ! matrix data
     ! populate the nodal attribute matrix at every node
      do j=1,numMeshNodes
         na(i)%xdmfMatrix(:,j) = na(i)%defaultVals(:)
      end do
      do j=1,na(i)%numNodesNotDefault
         na(i)%xdmfMatrix(:,na(i)%nonDefaultNodes(j)) = na(i)%nonDefaultVals(:,j)
      end do  
      attributeID = XdmfAddAttribute(xdmfFortranObj, trim(adjustl(na(i)%attrName))//CHAR(0), &
         XDMF_ATTRIBUTE_CENTER_NODE, attributeType, numValues, &
         XDMF_ARRAY_TYPE_FLOAT64, na(i)%xdmfMatrix)
   case default
      ! already accounted for prior to reading the data
   end select
end do

infoCount = 0
informationID = XdmfAddInformation(xdmfFortranObj, &
      'nodalAttributesComment'//char(0), trim(adjustl(nodalAttributesComment))//char(0))
write(numMeshNodesString,'(i0)') numMeshNodes
informationID = XdmfAddInformation(xdmfFortranObj, &
      'numMeshNodes'//char(0), trim(adjustl(numMeshNodesString))//char(0))      

do i=1,numNodalAttributes
   write(6,'(A)') 'INFO: Adding nodal attribute to XDMF.'
   ! set the metadata for this nodal attribute
   unitsID = XdmfAddInformation(xdmfFortranObj, &
      trim(adjustl(na(i)%attrName)) // ' units' // CHAR(0), &
      trim(adjustl(na(i)%units))//CHAR(0)) 
   write(numValsString,*) na(i)%numVals 
   numValsID = XdmfAddInformation(xdmfFortranObj, & 
      trim(adjustl(na(i)%attrName)) // ' number_of_values'//CHAR(0), &
      trim(adjustl(numValsString))//CHAR(0))
   write(defaultValsString,*) (na(i)%defaultVals(j), j=1,na(i)%numVals)
   defaultValsID = XdmfAddInformation(xdmfFortranObj,  &
      trim(adjustl(na(i)%attrName)) // ' default_values'//CHAR(0), &
      trim(adjustl(defaultValsString))//CHAR(0))
   !do j=1,3
   !   call xdmfInsertInformationIntoInformation(xdmfFortranObj, infoCount, infoCount+1, .true.)
   !end do
   infoCount = infoCount + 1
end do

write(6,'(a)') 'INFO: Finished writing nodal attributes data to XDMF.'
!-----------------------------------------------------------------------
end subroutine writeNodalAttributesXDMF
!-----------------------------------------------------------------------


!----------------------------------------------------------------------
!       S U B R O U T I N E    W R I T E   M E T A   D A T A 
!----------------------------------------------------------------------
! Write the metadata for the variable of interest.
!----------------------------------------------------------------------
subroutine writeMetaData(xdmfFortranObj, md)
use adcmesh
implicit none
include 'Xdmf.f' 
!
integer*8, intent(in) :: xdmfFortranObj
type(xdmfMetaData_t), intent(inout) :: md
!
if (md%createdIDs.eqv..true.) then
   call XdmfAddPreviousInformation(xdmfFortranObj, md%variable_name_id)
   call XdmfAddPreviousInformation(xdmfFortranObj, md%long_name_id)
   call XdmfAddPreviousInformation(xdmfFortranObj, md%standard_name_id)
   call XdmfAddPreviousInformation(xdmfFortranObj, md%coordinates_id)
  call XdmfAddPreviousInformation(xdmfFortranObj, md%units_id)
   if (trim(md%positive).ne.'null') then
      call XdmfAddPreviousInformation(xdmfFortranObj, md%positive_id)
   endif
else
   md%variable_name_id = XdmfAddInformation(xdmfFortranObj, 'variable_name'//CHAR(0), &
      trim(md%variable_name)//CHAR(0))
   md%long_name_id = XdmfAddInformation(xdmfFortranObj, 'long_name'//CHAR(0), &
     trim(md%long_name)//CHAR(0))
  md%standard_name_id = XdmfAddInformation(xdmfFortranObj, 'standard_name'//CHAR(0), &
     trim(md%standard_name)//CHAR(0))
   md%coordinates_id = XdmfAddInformation(xdmfFortranObj, 'coordinates'//CHAR(0), &
     trim(md%coordinates)//CHAR(0))
   md%units_id = XdmfAddInformation(xdmfFortranObj, 'units'//CHAR(0), &
     trim(md%units)//CHAR(0))
  if (trim(md%positive).ne.'null') then
     md%positive_id = XdmfAddInformation(xdmfFortranObj, 'positive'//CHAR(0), &
        trim(md%positive)//CHAR(0))
   endif
   md%createdIDs = .true.
endif
!----------------------------------------------------------------------
end subroutine writeMetaData
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!       S U B R O U T I N E   A D D   B O U N D A R I E S
!----------------------------------------------------------------------
! Add the boundaries as XdmfSet objects; if any attribute or 
! information objects have been created, they are added to the set and
! then cleared. We can use XdmfAttribute and XdmfInformation objects
! to set boundary types, etc. 
!----------------------------------------------------------------------
subroutine addBoundaries(xdmfFortranObj)
use adcmesh
implicit none
include 'Xdmf.f' 
!
integer*8, intent(in) :: xdmfFortranObj
character(1024) :: fluxBoundaryType ! character represenation of IBTYPE
!
!integer, allocatable :: xdmfElevB(:) ! 0-offset elevation boundary array
!integer, allocatable :: xdmfFluxB(:) ! 0-offset flux boundary array 
integer :: i, j  ! loop counter
!
write(6,'(A)') 'INFO: Adding boundary data.'
! create a zero offset boundary array long enough to hold any elevation boundary
!allocate(xdmfElevB(nvdll_max))
! create a zero offset boundary array long enough to hold any flux boundary
!allocate(xdmfFluxB(nvell_max))
do i=1, nope 
   elevationBoundaries(i)%nodes(:) = elevationBoundaries(i)%nodes(:) - 1
   !xdmfElevB(1:nvdll(i)) = elevationBoundaries(i)%nodes(:) - 1
   !do j=1,nvdll(i)
   !   write(6,'("DEBUG: elevBN: ",i0," xdmfElevBN: ",i0)') elevationBoundaries(i)%nodes(j), xdmfElevB(j)
   !end do
   ! jgftodo: this is hardcoded to 0 since this is the only type that ADCIRC
   ! supports, and most or all mesh files don't even have this value, although
   ! the documentation calls for it 
   elevationBoundaries(i)%informationID = XdmfAddInformation(xdmfFortranObj, &
      'IBTYPEE'//CHAR(0), '0'//CHAR(0))
   elevationBoundaries(i)%setID = XdmfAddSet(xdmfFortranObj, &
      'elevation_specified_boundary'//char(0), XDMF_SET_TYPE_NODE, &
      elevationBoundaries(i)%nodes, nvdll(i), XDMF_ARRAY_TYPE_INT32)
end do
sfCount = 1
efCount = 1
ifCount = 1
ifwpCount = 1 
do i=1, nbou 
   write(fluxBoundaryType,'(i0)') ibtype_orig(i) 
   select case(ibtype_orig(i))
   case(0,1,2,10,11,12,20,21,22,30,52)
      !xdmfFluxB(1:nvell(i)) = simpleFluxBoundaries(sfCount)%nodes(:) - 1
      simpleFluxBoundaries(sfCount)%informationID = &
         XdmfAddInformation(xdmfFortranObj, 'IBTYPE'//CHAR(0), &
         trim(fluxBoundaryType)//CHAR(0))
      simpleFluxBoundaries(sfCount)%nodes(:) = simpleFluxBoundaries(sfCount)%nodes(:) - 1
      simpleFluxBoundaries(sfCount)%setID = XdmfAddSet(xdmfFortranObj, &
         'flux_specified_boundary'//char(0), XDMF_SET_TYPE_NODE, &
         simpleFluxBoundaries(sfCount)%nodes(:), nvell(i), XDMF_ARRAY_TYPE_INT32)
      sfCount = sfCount + 1
   case(3,13,23)
      !xdmfFluxB(1:nvell(i)) = externalFluxBoundaries(efCount)%nodes(:) - 1
      externalFluxBoundaries(efCount)%attributeIDs(1) = &
         XdmfAddAttribute(xdmfFortranObj, 'BARLANHT'//CHAR(0), &
         XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), &
         XDMF_ARRAY_TYPE_FLOAT64, externalFluxBoundaries(efCount)%barlanht)
      externalFluxBoundaries(efCount)%attributeIDs(2) = &
         XdmfAddAttribute(xdmfFortranObj, 'BARLANCFSP'//CHAR(0), &
         XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), &
         XDMF_ARRAY_TYPE_FLOAT64, externalFluxBoundaries(efCount)%barlancfsp)
      externalFluxBoundaries(efCount)%informationID = &
         XdmfAddInformation(xdmfFortranObj, 'IBTYPE'//CHAR(0), &
         trim(fluxBoundaryType)//CHAR(0))
      externalFluxBoundaries(efCount)%nodes(:) = externalFluxBoundaries(efCount)%nodes(:) - 1
      externalFluxBoundaries(efCount)%setID = XdmfAddSet(xdmfFortranObj, &
         'flux_specified_boundary'//char(0), XDMF_SET_TYPE_NODE, &
         externalFluxBoundaries(efCount)%nodes(:), nvell(i), XDMF_ARRAY_TYPE_INT32)
      efCount = efCount + 1
   case(4,24)
      !xdmfFluxB(1:nvell(i)) = internalFluxBoundaries(ifCount)%ibconn(:) - 1
      !xdmfFluxB(1:nvell(i)) = internalFluxBoundaries(ifCount)%nodes(:) - 1
      internalFluxBoundaries(ifCount)%ibconn(:) = internalFluxBoundaries(ifCount)%ibconn(:) - 1
      internalFluxBoundaries(ifCount)%attributeIDs(1) = &
         XdmfAddAttribute(xdmfFortranObj, 'IBCONN'//CHAR(0), &
         XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), &
         XDMF_ARRAY_TYPE_INT32, internalFluxBoundaries(ifCount)%ibconn) ! zero offset
      internalFluxBoundaries(ifCount)%attributeIDs(2) = &
         XdmfAddAttribute(xdmfFortranObj, 'BARINHT'//CHAR(0), &
         XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), &
         XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundaries(ifCount)%barinht)
      internalFluxBoundaries(ifCount)%attributeIDs(3) = &
         XdmfAddAttribute(xdmfFortranObj, 'BARINCFSB'//CHAR(0), &
         XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), &
         XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundaries(ifCount)%barincfsb)
      internalFluxBoundaries(ifCount)%attributeIDs(4) = &
         XdmfAddAttribute(xdmfFortranObj, 'BARINCFSP'//CHAR(0), &
         XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), &
         XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundaries(ifCount)%barincfsp)
      internalFluxBoundaries(ifCount)%informationID = XdmfAddInformation(xdmfFortranObj, &
         'IBTYPE'//CHAR(0), trim(fluxBoundaryType)//CHAR(0))
      internalFluxBoundaries(ifCount)%nodes(:) = internalFluxBoundaries(ifCount)%nodes(:) - 1
      internalFluxBoundaries(ifCount)%setID = XdmfAddSet(xdmfFortranObj, &
         'flux_specified_boundary'//char(0), XDMF_SET_TYPE_NODE, &
         internalFluxBoundaries(ifCount)%nodes, nvell(i), XDMF_ARRAY_TYPE_INT32)
      ifCount = ifCount + 1
   case(5,25)
      !xdmfFluxB(1:nvell(i)) = internalFluxBoundariesWithPipes(ifwpCount)%ibconn(:) - 1   
      internalFluxBoundariesWithPipes(ifwpCount)%ibconn(:) = &
         internalFluxBoundariesWithPipes(ifwpCount)%ibconn(:) - 1   
      internalFluxBoundariesWithPipes(ifwpCount)%attributeIDs(1) = &
         XdmfAddAttribute(xdmfFortranObj, 'IBCONN'//CHAR(0), &
         XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), &
         XDMF_ARRAY_TYPE_INT32, internalFluxBoundariesWithPipes(ifwpCount)%ibconn(:)) ! zero offset
      internalFluxBoundariesWithPipes(ifwpCount)%attributeIDs(2) = &
         XdmfAddAttribute(xdmfFortranObj, 'BARINHT'//CHAR(0), &
         XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), &
         XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundariesWithPipes(ifwpCount)%barinht)
      internalFluxBoundariesWithPipes(ifwpCount)%attributeIDs(3) = &
         XdmfAddAttribute(xdmfFortranObj, 'BARINCFSB'//CHAR(0), &
         XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), &
         XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundariesWithPipes(ifwpCount)%barincfsb)
      internalFluxBoundariesWithPipes(ifwpCount)%attributeIDs(4) = &
         XdmfAddAttribute(xdmfFortranObj, 'BARINCFSP'//CHAR(0), &
         XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), &
         XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundariesWithPipes(ifwpCount)%barincfsp)   
      internalFluxBoundariesWithPipes(ifwpCount)%attributeIDs(5) = &
         XdmfAddAttribute(xdmfFortranObj, 'PIPEHT'//CHAR(0), &
         XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), &
         XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundariesWithPipes(ifwpCount)%pipeht)
      internalFluxBoundariesWithPipes(ifwpCount)%attributeIDs(6) = &
         XdmfAddAttribute(xdmfFortranObj, 'PIPECOEF'//CHAR(0), &
         XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), &
         XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundariesWithPipes(ifwpCount)%pipecoef)
      internalFluxBoundariesWithPipes(ifwpCount)%attributeIDs(7) = &
         XdmfAddAttribute(xdmfFortranObj, 'PIPEDIAM'//CHAR(0), &
         XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), &
         XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundariesWithPipes(ifwpCount)%pipediam)
      internalFluxBoundariesWithPipes(ifwpCount)%informationID = &
         XdmfAddInformation(xdmfFortranObj, 'IBTYPE'//CHAR(0), &
         trim(fluxBoundaryType)//CHAR(0))
      !xdmfFluxB(1:nvell(i)) = internalFluxBoundariesWithPipes(ifwpCount)%nodes(:) - 1
      internalFluxBoundariesWithPipes(ifwpCount)%nodes(:) = &
         internalFluxBoundariesWithPipes(ifwpCount)%nodes(:) - 1
      internalFluxBoundariesWithPipes(ifwpCount)%setID = &
         XdmfAddSet(xdmfFortranObj,'flux_specified_boundary'//char(0), & 
         XDMF_SET_TYPE_NODE, internalFluxBoundariesWithPipes(ifwpCount)%nodes, &
         nvell(i), XDMF_ARRAY_TYPE_INT32)
      ifwpCount = ifwpCount + 1
   case default
      write(6,'("ERROR: File contains IBTYPE=",i0," which is not a valid flux boundary type.")'), ibtype_orig(i)
   end select
end do
!----------------------------------------------------------------------
end subroutine addBoundaries
!----------------------------------------------------------------------


!----------------------------------------------------------------------
! S U B R O U T I N E   A D D   B O U N D A R Y  R E F E R E N C E S
!----------------------------------------------------------------------
! Add boundary references to each grid in the collection.
!----------------------------------------------------------------------
subroutine addBoundaryReferences(xdmfFortranObj)
use adcmesh
implicit none
include 'Xdmf.f' 
!
integer*8, intent(in) :: xdmfFortranObj
character(len=256) :: fluxBoundaryType
integer :: i, j  ! loop counter
!
do i=1, nope 
   elevationBoundaries(i)%informationID = XdmfAddInformation(xdmfFortranObj, &
      'IBTYPEE'//CHAR(0), '0'//CHAR(0))
   elevationBoundaries(i)%setID = XdmfAddSet(xdmfFortranObj, &
      'elevation_specified_boundary'//char(0), XDMF_SET_TYPE_NODE, &
      elevationBoundaries(i)%nodes, nvdll(i), XDMF_ARRAY_TYPE_INT32)
end do
sfCount = 1
efCount = 1
ifCount = 1
ifwpCount = 1 
do i=1, nbou 
   write(fluxBoundaryType,'(i0)') ibtype_orig(i)
   select case(ibtype_orig(i))
   case(0,1,2,10,11,12,20,21,22,30,52)  
      simpleFluxBoundaries(sfCount)%informationID = &
         XdmfAddInformation(xdmfFortranObj, 'IBTYPE'//CHAR(0), &
         trim(fluxBoundaryType)//CHAR(0))
      simpleFluxBoundaries(sfCount)%setID = XdmfAddSet(xdmfFortranObj, &
         'flux_specified_boundary'//char(0), XDMF_SET_TYPE_NODE, &
         simpleFluxBoundaries(sfCount)%nodes(:), nvell(i), XDMF_ARRAY_TYPE_INT32)
      sfCount = sfCount + 1
   case(3,13,23)
      do j=1, externalFluxBoundaries(efCount)%numAttributes
         call xdmfAddPreviousAttribute(xdmfFortranObj, &
               externalFluxBoundaries(efCount)%attributeIDs(j))
      end do
      externalFluxBoundaries(efCount)%informationID = &
         XdmfAddInformation(xdmfFortranObj, 'IBTYPE'//CHAR(0), &
         trim(fluxBoundaryType)//CHAR(0))
      externalFluxBoundaries(efCount)%setID = XdmfAddSet(xdmfFortranObj, &
         'flux_specified_boundary'//char(0), XDMF_SET_TYPE_NODE, &
         externalFluxBoundaries(efCount)%nodes(:), nvell(i), XDMF_ARRAY_TYPE_INT32)
      efCount = efCount + 1
   case(4,24)
      do j=1, internalFluxBoundaries(ifCount)%numAttributes
         call xdmfAddPreviousAttribute(xdmfFortranObj, &
               internalFluxBoundaries(ifCount)%attributeIDs(j))
      end do
      internalFluxBoundaries(ifCount)%informationID = XdmfAddInformation(xdmfFortranObj, &
         'IBTYPE'//CHAR(0), trim(fluxBoundaryType)//CHAR(0))
      internalFluxBoundaries(ifCount)%setID = XdmfAddSet(xdmfFortranObj, &
         'flux_specified_boundary'//char(0), XDMF_SET_TYPE_NODE, &
         internalFluxBoundaries(ifCount)%nodes, nvell(i), XDMF_ARRAY_TYPE_INT32)
      ifCount = ifCount + 1
   case(5,25)
      do j=1, internalFluxBoundariesWithPipes(ifwpCount)%numAttributes
         call xdmfAddPreviousAttribute(xdmfFortranObj, &
               internalFluxBoundariesWithPipes(ifwpCount)%attributeIDs(j))
      end do
      internalFluxBoundariesWithPipes(ifwpCount)%informationID = &
         XdmfAddInformation(xdmfFortranObj, 'IBTYPE'//CHAR(0), &
         trim(fluxBoundaryType)//CHAR(0))
      internalFluxBoundariesWithPipes(ifwpCount)%setID = &
         XdmfAddSet(xdmfFortranObj,'flux_specified_boundary'//char(0), & 
         XDMF_SET_TYPE_NODE, internalFluxBoundariesWithPipes(ifwpCount)%nodes, &
         nvell(i), XDMF_ARRAY_TYPE_INT32)
      ifwpCount = ifwpCount + 1
   case default
      write(6,'("ERROR: File contains IBTYPE=",i0," which is not a valid flux boundary type.")'), ibtype_orig(i)
   end select
end do
!----------------------------------------------------------------------
end subroutine addBoundaryReferences
!----------------------------------------------------------------------

