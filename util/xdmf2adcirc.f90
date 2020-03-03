!-------------------------------------------------------------------------
! xdmf2adcirc.f90
!-------------------------------------------------------------------------
! Author: Jason Fleming (jason.fleming@seahorsecoastal.com) 
!
! Convert XDMF data to ascii adcirc data. 
!---------------------------------------------------------------------
!---------------------------------------------------------------------
program xdmf2adcirc
!---------------------------------------------------------------------
!---------------------------------------------------------------------
use netcdf
use adcmesh
use nodalattr
use control
implicit none
include 'Xdmf.f' 
! 
character(1024) :: datafilename ! full path name of the xmf file to be converted
character(1024) :: controlFileName  ! full path name of the adcirc fort.15 for metadata  
integer*8 :: xdmfFortranObj  ! represents pointer to the XdmfFortran object
!$omp threadprivate(xdmfFortranObj)
integer :: attributeType     ! attribute type
!$omp threadprivate(attributeType)
integer :: startIndex
!$omp threadprivate(startIndex)
integer :: arrayStride
!$omp threadprivate(arrayStride)
integer :: valueStride
!$omp threadprivate(valueStride)
integer :: attributeIndex
!$omp threadprivate(attributeIndex)
integer :: attributeDataType
!$omp threadprivate(attributeDataType)
integer :: attributePropertyIndex
!$omp threadprivate(attributePropertyIndex)
integer :: numAttributes
!$omp threadprivate(numAttributes)
integer :: numAttributeProperties
!$omp threadprivate(numAttributeProperties)
integer :: unitNumber        ! fortran i/o unit number for ascii data file 
!$omp threadprivate(unitNumber)
integer :: tsinterval          ! spooling interval in time steps
!$omp threadprivate(tsinterval)
real(8) :: tinterval          ! spooling interval in seconds
!$omp threadprivate(tinterval)
integer :: numSnaps          ! unreliable number of datasets in ascii file
!$omp threadprivate(numSnaps)
integer :: numValues         ! total number of values in a dataset
!$omp threadprivate(numValues)
real(8) :: timeSec            ! time of a dataset, in seconds
!$omp threadprivate(timeSec)
integer :: timeStep          ! time step number for a dataset
!$omp threadprivate(timeStep)
integer :: ss                ! dataset counter
!$omp threadprivate(ss)
integer :: nCol              ! number of columns of data in ascii file
!$omp threadprivate(nCol)
integer :: numNodes          ! number of nodes in mesh according to ascii file
!$omp threadprivate(numNodes)
real(8), allocatable :: data_array1(:)
!$omp threadprivate(data_array1)
real(8), allocatable :: data_array2(:,:)
!$omp threadprivate(data_array2)
integer :: gridCollectionIndex ! the grid collections are numbered starting from zero
!$omp threadprivate(gridCollectionIndex)
integer :: informationPropertyIndex !
!$omp threadprivate(informationPropertyIndex)
!
! the following four integers act like logicals (1=true, 2=false), and 
! they control whether the associated data should be opened when a grid 
! collection is opened  
integer :: openMaps ! 1 if Maps should be opened within another object
!$omp threadprivate(openMaps)
integer :: openAttributes ! 1 if Attributes should be opened within another object 
!$omp threadprivate(openAttributes)
integer :: openInformations ! 1 if Informations should be opened within another object
!$omp threadprivate(openInformations)
integer :: openSets ! 1 if sets should be opened within another object
!$omp threadprivate(openSets)
!
integer :: informationIndex              ! index of the information   
!$omp threadprivate(informationIndex)
integer, parameter :: keyLength = 256
integer, parameter :: valueLength = 256
integer, parameter :: tagLength = 256
integer, parameter :: nameLength = 256
character(len=keyLength) :: itemKey
character(len=valueLength) :: itemValue
character(len=tagLength) :: itemTag
character(len=nameLength) :: itemName
integer :: gridCollectionType
!$omp threadprivate(gridCollectionType)
integer :: numGridCollections
!$omp threadprivate(numGridCollections)
character(len=256) :: gridCollectionTypeString
character(len=256) :: dataTypeString
character(len=256) :: gridName
integer :: numGrids
!$omp threadprivate(numGrids)
integer :: infoIndex
!$omp threadprivate(infoIndex)
logical :: meshonly
!$omp threadprivate(meshonly)
integer :: typeHolder
!$omp threadprivate(typeHolder)
integer :: arrayIndex
!$omp threadprivate(arrayIndex)
integer :: valueIndex 
!$omp threadprivate(valueIndex)
integer :: startingIndex
!$omp threadprivate(startingIndex)
integer :: gridIndex
!$omp threadprivate(gridIndex)
integer :: propertyIndex
!$omp threadprivate(propertyIndex)
integer :: numContained
!$omp threadprivate(numContained)
integer :: numInformation
!$omp threadprivate(numInformation)
integer :: num_components
!$omp threadprivate(num_components)
real(8), allocatable :: adcirc_data(:)
!$omp threadprivate(adcirc_data)
character(len=1025) :: ascii_datafile_name
character(len=256) :: attributeName
character(len=256) :: myAttributeName
!
character(1024) :: cmdlinearg
character(1024) :: cmdlineopt
character(80) :: line ! a line of data from the ascii file
character(80) :: topcomment ! comment line at the top of ascii file
integer :: argcount
!$omp threadprivate(argcount)
integer :: i, j, k, n, p
!$omp threadprivate(i,j,k,n,p)
integer :: numInformations
!$omp threadprivate(numInformations)
integer :: numProperties
!$omp threadprivate(numProperties)
integer :: numGridCollectionGrids
!$omp threadprivate(numGridCollectionGrids)
integer :: domainNumGrids
!$omp threadprivate(domainNumGrids)
!
! initializations
dataFileName = 'null'
meshFileName = 'fort.14'
controlFileName = 'fort.15'
SS=1  ! initialize the dataset counter
meshonly = .false.
i=0
startIndex = 0
arrayStride = 1
valueStride = 1
openMaps = 1
openAttributes = 1
openInformations = 1
openSets = 1
!
! process command line options
argcount = iargc() ! count up command line options
write(6,'("INFO: There are ",i0," command line arguments.")') argcount
do while (i.lt.argcount)
   i = i + 1
   call getarg(i, cmdlineopt)
   select case(trim(cmdlineopt))
   case("--verbose")
      write(6+CK_LUN,'(a)') "INFO: Processing " // trim(cmdlineopt) // "."
      verbose = .true.
   case("--meshonly")
      write(6+CK_LUN,'(a)') "INFO: Processing " // trim(cmdlineopt) // "."
      meshonly = .true.
   case("--datafile")
      i = i + 1
      call getarg(i, datafilename)     
      write(6+CK_LUN,'(a)') "INFO: Processing " // trim(cmdlineopt) // &
         " " // trim(datafilename) // "."
   case default
      write(6+CK_LUN,'(a)') "WARNING: Command line option '",TRIM(cmdlineopt),"' was not recognized."
   end select
end do
!
write(6,'(a)') 'INFO: Reading data from the "' // trim(adjustl(datafilename)) &
   // '" XDMF file.'
call xdmfinit(xdmfFortranObj)
call xdmfRead(xdmfFortranObj, trim(adjustl(datafilename))//char(0))
!
!   N O D A L   A T T R I B U T E   D A T A
! 
if (trim(datafilename).eq.'fort.13.xmf') then
   call readNodalAttributesXDMF(xdmfFortranObj)
   datafilename = 'xdmf_fort.13'
   call writeNodalAttributesFile(datafilename)
   call xdmfClose(xdmfFortranObj)
   stop
endif
!
!  G R I D   C O L L E C T I O N
!
call xdmfRetrieveNumDomainGridCollections(xdmfFortranObj, numGridCollections)
write(6,'("INFO: Number of GridCollections in this file : ",i0,".")') numGridCollections
if (numGridCollections.gt.0) then
   write(6+CK_LUN,'(a)') 'INFO: Opening GridCollection.'
   gridCollectionIndex = 0 ! the first grid collection (they're numbered from zero)
   openMaps = 1
   openAttributes = 1
   openInformations = 1
   openSets = 1
   call xdmfOpenDomainGridCollection(xdmfFortranObj, gridCollectionIndex, & 
      openMaps, openAttributes, openInformations, openSets)
   !
   ! determine the type of grid collection we are working with
   call xdmfRetrieveGridCollectionType(xdmfFortranObj, gridCollectionType)
   select case(gridCollectionType)
   case(400)
       gridCollectionTypeString = 'XDMF_GRID_COLLECTION_TYPE_SPATIAL'
   case(401)
       gridCollectionTypeString = 'XDMF_GRID_COLLECTION_TYPE_TEMPORAL'
   case default
       write(6+CK_LUN,'("WARNING: The grid collection type code ",i0," is not supported by xdmf2adcirc.")') gridCollectionType
   end select
   write(6+CK_LUN,'(a)') 'INFO: The grid collection type is ' // &
      trim(gridCollectionTypeString) // '.' 
   call xdmfRetrieveGridCollectionNumGrids(xdmfFortranObj, XDMF_GRID_TYPE_UNSTRUCTURED,  numGrids)
   write(6+CK_LUN,'("INFO: Number of Grids contained in the GridCollection: ",i0)'), numGrids
endif
!
!    C O N T R O L   ( F O R T . 1 5 )   D A T A
!  
!call xdmfRetrieveNumInformation(xdmfFortranObj, numInformations)
!write(6,'("DEBUG: Number of Information items for this GridCollection : ",i0,".")') numInformations
!
! Most of these properties are related to the fort.15 control parameters
!do i=0,numInformations-1
!   write(6,'("INFO: Information ",i0,".")') i
!   call xdmfRetrieveInformationNumProperties(xdmfFortranObj, &
!     i, numProperties)
!   write(6,'("INFO: Information ",i0," contains ",i0," Properties.")') &
!      i, numProperties
!   do p=0,numProperties-1
!      call xdmfRetrieveInformationProperty(xdmfFortranObj, i,  &
!         p, itemKey, keyLength, itemValue, valueLength)
!      write(6,'("INFO: Information ",i0," Property ",i0," Key: ",a,".")') &
!         i, p, trim(itemKey)
!      write(6,'("INFO: Information ",i0," Property ",i0," Value: ",a,".")') &
!         i, p, trim(itemValue)   
!   end do
!end do 

!
!   M E S H   D A T A 
!
call readMeshXDMF(xdmfFortranObj)
meshFileName = 'xdmf_fort.14'
call writeMesh()
if (meshonly.eqv..true.) then
   call xdmfClose(xdmfFortranObj)
   stop
endif 

!
!  A D C I R C   O U T P U T   D A T A 
!
call xdmfClearAttributes(xdmfFortranObj)
gridIndex = 0
write(6,'(a)') 'INFO: Opening the first Grid in the GridCollection.'
call xdmfOpenGridCollectionGrid(xdmfFortranObj, XDMF_GRID_TYPE_UNSTRUCTURED,  &
   gridIndex, openMaps, openAttributes, openInformations, openSets)
write(6,'(a)') 'DEBUG: Reading the Grid name.'
call xdmfRetrieveGridCollectionGridName(xdmfFortranObj, XDMF_GRID_TYPE_UNSTRUCTURED, &
   gridIndex,  gridName, nameLength)
call replaceNullsWithSpaces(gridName)
agrid(1:80) = gridName(1:80)
write(6,'("DEBUG: Grid ",i0," of GridCollection ",i0," is named ",a,".")') &
   gridIndex, gridCollectionIndex, '"' // trim(gridName) // '"'
!
call xdmfRetrieveNumAttributes(xdmfFortranObj, numAttributes)
write(6,'("INFO: Grid ",i0," contains ",i0," attributes.")') gridIndex, numAttributes
do attributeIndex=0, numAttributes - 1
   call xdmfRetrieveAttributeName(xdmfFortranObj, attributeIndex,  & 
      myAttributeName, nameLength)
   call replaceNullsWithSpaces(myAttributeName)
   write(6+CK_LUN,'("INFO: Grid ",i0," Attribute ",i0," is named ",a)') gridIndex, attributeIndex, trim(myAttributeName)
   select case(trim(myAttributeName))
   case("zeta")
      write(6+CK_LUN,'(a)') 'INFO: Preparing to write an ADCIRC water surface elevation file.'
      ascii_datafile_name = "fort.63"
      num_components = 1
      exit
   case("zeta_max")
      write(6+CK_LUN,'(a)') 'INFO: Preparing to write an ADCIRC maximum water surface elevation file.'
      ascii_datafile_name = "maxele.63"
      num_components = 1
      exit
   case("vel")
      write(6+CK_LUN,'(a)') 'INFO: Preparing to write an ADCIRC water velocity file.'
      ascii_datafile_name = "fort.64"
      num_components = 2
      exit
   case("pressure")
      write(6+CK_LUN,'(a)') 'INFO: Preparing to write an ADCIRC barometric pressure file.'
      ascii_datafile_name = "fort.73"
      num_components = 1
      exit
   case("wvel")
      write(6+CK_LUN,'(a)') 'INFO: Preparing to write an ADCIRC wind velocity file.'
      ascii_datafile_name = "fort.74"
      num_components = 2
      exit
   case("ice")
      write(6+CK_LUN,'(a)') 'INFO: Preparing to write an ADCIRC ice coverage file.'
      ascii_datafile_name = "fort.93"
      num_components = 1
      exit
   case("wind_max")
      write(6+CK_LUN,'(a)') 'INFO: Preparing to write an ADCIRC maximum wind speed file.'
      ascii_datafile_name = "maxwvel.63"
      num_components =  1     
      exit
   case("dir")
      write(6+CK_LUN,'(a)') 'INFO: Preparing to write a mean wave direction file.'
      ascii_datafile_name = "swan_DIR.63"
      num_components = 1
      exit
   case("hs")
      write(6+CK_LUN,'(a)') 'INFO: Preparing to write a significant wave height file.'
      ascii_datafile_name = "swan_HS.63"
      num_components = 1
      exit
   case("tmm10")
      write(6+CK_LUN,'(a)') 'INFO: Preparing to write a mean absolute wave period file.'
      ascii_datafile_name = "swan_TMM10.63"
      num_components = 1
      exit
   case("tps")
      write(6+CK_LUN,'(a)') 'INFO: Preparing to write a relative peak period file.'
      ascii_datafile_name = "swan_TPS.63"
      num_components = 1
      exit
   case("swan_HS_max")
      write(6+CK_LUN,'(a)') 'INFO: Preparing to write a maximum significant wave height file.'
      ascii_datafile_name = "swan_HS_max.63"
      num_components = 1
      exit
   case("swan_TPS_max")
      write(6+CK_LUN,'(a)') 'INFO: Preparing to write an maximum relative peak wave period file.'
      ascii_datafile_name = "swan_TPS_max.63"
      num_components = 1                              
      exit
   case("depth")
      ! do nothing, this will be written with the mesh
   case default
      write(6+CK_LUN,'("ERROR: Unrecognized variable name: ",a,".")') trim(itemName)
      stop
   end select
end do
write(6,*) "INFO: " // trim(ascii_datafile_name)
if (num_components.eq.1) then
   attributeType = XDMF_ATTRIBUTE_TYPE_SCALAR
   numValues = np
endif
if (num_components.eq.2) then
   attributeType = XDMF_ATTRIBUTE_TYPE_VECTOR
   numValues = 3*np
endif
allocate(adcirc_data(numValues))
!
! open the ascii adcirc file that will hold the data
open(11,file=trim(ascii_datafile_name),status='replace',action='write')
! write header info
write(11,'(A)') trim(agrid)
write(11,1010) numGrids, np, -99999.0, -99999, num_components
timeStep = -99999
do gridIndex=0,numGrids-1
   ! clear out the info from the previous mesh
   call xdmfClearAttributes(xdmfFortranObj)
   call xdmfClearInformations(xdmfFortranObj)
   ! open the current mesh
   call xdmfOpenGridCollectionGrid(xdmfFortranObj, XDMF_GRID_TYPE_UNSTRUCTURED,  &
      gridIndex, openMaps, openAttributes, openInformations, openSets)
   ! time in seconds
   call xdmfRetrieveTime(xdmfFortranObj, timeSec)
   ! time step number
   call xdmfRetrieveNumInformation(xdmfFortranObj,numInformation)
   do informationIndex=0,numInformation-1
      call xdmfRetrieveInformation(xdmfFortranObj, informationIndex, itemKey, keyLength, itemValue, valueLength)
      call replaceNullsWithSpaces(itemKey)
      if (trim(itemKey).eq.'IT') then
         call replaceNullsWithSpaces(itemValue)
         read(itemValue+CK_LUN,*) timeStep
      endif
   end do
   ! the actual adcirc data
   call xdmfRetrieveNumAttributes(xdmfFortranObj, numAttributes)
   ! looking for the attribute of interest
   do attributeIndex=0, numAttributes - 1
      call xdmfRetrieveAttributeName(xdmfFortranObj, attributeIndex,  & 
         attributeName, nameLength)
      call replaceNullsWithSpaces(attributeName)
      ! is this the attribute we are looking for?
      if (trim(attributeName).eq.trim(myAttributeName)) then
         ! read the xdmf dataset
         call xdmfRetrieveAttributeValues(xdmfFortranObj, attributeIndex, & 
            adcirc_data, XDMF_ARRAY_TYPE_FLOAT64, numValues, &
            startIndex, arrayStride, valueStride) 
         ! write the dataset to ascii adcirc format
         write(11+CK_LUN,2120) timeSec, timeStep
         if (num_components.eq.1) then ! scalar
            do k=1,numValues
               write(11+CK_LUN,2453) k, adcirc_data(k)
            end do
         endif
         if (num_components.eq.2) then ! faux 3-component vector
            do k=1,numValues,3
               write(11+CK_LUN,2453) k, adcirc_data(k), adcirc_data(k+1) 
            end do
         endif
      endif
   end do
   write(6+CK_LUN,advance='no',fmt='(I4)') gridIndex+1
enddo
write(6,'(/,A)') "INFO: ... finished writing file."
write(6,'("INFO: Wrote ",i0," data sets.")') numGrids
close(11)
!
 1010 FORMAT(1X,I0,1X,I0,1X,E15.7E3,1X,I0,1X,I0,1X,'FileFmtVersion: ',I0)
 1011 FORMAT(1X,I0,1X,I0,1X,E15.7E3,1X,I0,1X,I0,1X,I0,1X,'FileFmtVersion: ',I0)
 2120 FORMAT(2X,1pE20.10E3,5X,I0)
 2453 FORMAT(2x, i0, 2x, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3)

call xdmfclose(xdmfFortranObj)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
end program xdmf2adcirc
!----------------------------------------------------------------------
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!                     S U B R O U T I N E    
!        R E A D   N O D A L   A T T R I B U T E S   X D M F
!-----------------------------------------------------------------------
! jgf: Reads all nodal attributes from an XDMF file.
!-----------------------------------------------------------------------
subroutine readNodalAttributesXDMF(xdmfFortranObj)
      use CkLunMod, only : CK_LUN
use adcmesh, only : np
use nodalattr
implicit none
include 'Xdmf.f'
integer*8, intent(in) :: xdmfFortranObj ! object that receives the data
!
integer :: insideInformationIndex
integer :: numInsideInformation
integer :: startIndex
integer :: arrayStride
integer :: valueStride
integer :: attributeIndex
integer :: attributeType
integer :: informationIndex
logical :: isNodalAttribute
integer :: gridIndex = 0
!$omp threadprivate(gridIndex)
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
integer :: i, j, k, nattrCount, nonDefaultCount
!
startIndex = 0
arrayStride = 1
valueStride = 1
openMaps = 1
openAttributes = 1
openInformations = 1
openSets = 1
!
call xdmfOpenDomainGrid(xdmfFortranObj, XDMF_GRID_TYPE_UNSTRUCTURED, &
   gridIndex, openMaps, openAttributes, openInformations, openSets)
!
call xdmfRetrieveNumAttributes(xdmfFortranObj, numAttributes)
numNodalAttributes = numAttributes - 1 ! the depth is included so don't count it
allocate(na(numNodalAttributes))
write(6,'("INFO: Grid ",i0," contains ",i0," nodal attributes.")') gridIndex, numNodalAttributes
!
! populate the names of the nodal attributes
nattrCount = 1
do attributeIndex=0, numAttributes - 1
   call xdmfRetrieveAttributeName(xdmfFortranObj, attributeIndex,  & 
      itemName, nameLength)
   call replaceNullsWithSpaces(itemName)
   if (trim(itemName).eq.'depth') then
      cycle
   endif
   write(6+CK_LUN,'("INFO: Grid ",i0," Attribute ",i0," is named ",a)') gridIndex, attributeIndex, trim(itemName)
   na(nattrCount)%attrName = trim(itemName)
   nattrCount = nattrCount + 1
end do
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
      read(itemValue+CK_LUN,*) numMeshNodes
   case default
      do i=1,numNodalAttributes
         if (trim(itemKey).eq. trim(na(i)%attrName) // ' number_of_values') then
            read(itemValue+CK_LUN,*) na(i)%numVals
            allocate(na(i)%defaultVals(na(i)%numVals))
            exit
         endif
      end do
   end select
end do
allocate(diff(numMeshNodes))
allocate(areDefaultValues(numMeshNodes))
!
! populate the units and the default values
call xdmfRetrieveNumInformation(xdmfFortranObj, numInformation)
do informationIndex=0,numInformation-1
   call xdmfRetrieveInformation(xdmfFortranObj, informationIndex, itemKey, keyLength, itemValue, valueLength)
   call replaceNullsWithSpaces(itemKey)
   call replaceNullsWithSpaces(itemValue)
   do i=1,numNodalAttributes
      if (trim(itemKey).eq. trim(na(i)%attrName)// ' units') then
         na(i) % units = trim(itemValue)
      endif  
      if (trim(itemKey).eq. trim(na(i)%attrName) // ' default_values') then
         write(6+CK_LUN,*) 'default_values'//trim(itemValue)
         read(itemValue+CK_LUN,*) (na(i)%defaultVals(j),j=1,na(i)%numVals)
      endif
   end do
end do
!
! populate the data
do attributeIndex=0, numAttributes - 1
   call xdmfRetrieveAttributeName(xdmfFortranObj, attributeIndex,  & 
      itemName, nameLength)
   call replaceNullsWithSpaces(itemName)
   if (trim(itemName).eq.'depth') then
      cycle
   endif
   do i=1,numNodalAttributes
      if (trim(itemName).eq.trim(na(i)%attrName)) then
         write(6+CK_LUN,'(a)') 'loading nodal attribute data for '//trim(itemName)
         if (na(i)%numVals.eq.1) then
            allocate(na(i)%xdmfArray(numMeshNodes))
            attributeType = XDMF_ATTRIBUTE_TYPE_SCALAR
            numValues = numMeshNodes * na(i) % numVals
            call xdmfRetrieveAttributeValues(xdmfFortranObj, attributeIndex, & 
               na(i)%xdmfArray, XDMF_ARRAY_TYPE_FLOAT64, numValues, &
                  startIndex, arrayStride, valueStride)
            ! determine the number of nondefault values
            !
            ! machine precision prevents us from simply checking whether the 
            ! value .ne. the default value
            diff = abs(na(i)%xdmfArray - na(i)%defaultVals(1))
            na(i)%numNodesNotDefault = count(diff.gt.1.e-6)        
            ! now allocate space for the non default values and populate them
            allocate(na(i)%nonDefaultVals(1,na(i)%numNodesNotDefault))
            allocate(na(i)%nonDefaultNodes(na(i)%numNodesNotDefault))
            ! now record the node number and value where the values are not
            ! the default
            nonDefaultCount = 1
            do j=1,numMeshNodes
               if (diff(j).gt.1.e-6) then
                  na(i)%nonDefaultNodes(nonDefaultCount) = j
                  na(i)%nonDefaultVals(1,nonDefaultCount) = na(i)%xdmfArray(j)
                  nonDefaultCount = nonDefaultCount + 1
               endif
            end do
         else
            attributeType = XDMF_ATTRIBUTE_TYPE_MATRIX
            numValues = numMeshNodes * na(i) % numVals
            allocate(na(i)%xdmfMatrix(na(i)%numVals,numMeshNodes))
            call xdmfRetrieveAttributeValues(xdmfFortranObj, attributeIndex, & 
               na(i)%xdmfMatrix, XDMF_ARRAY_TYPE_FLOAT64, numValues, &
                  startIndex, arrayStride, valueStride)
            ! determine the number of nondefault values
            areDefaultValues = .true.
            do j=1,numMeshNodes
               do k=1,na(i)%numVals
                  if (abs(na(i)%xdmfMatrix(k,j)-na(i)%defaultVals(k)).gt.1.e-6) then
                     areDefaultValues(j) = .false.
                  endif
               enddo
            enddo
            ! now allocate space for the non default values and populate them
            na(i)%numNodesNotDefault = count(areDefaultValues.eqv..false.) 
            allocate(na(i)%nonDefaultVals(na(i)%numVals,na(i)%numNodesNotDefault))
            allocate(na(i)%nonDefaultNodes(na(i)%numNodesNotDefault))
            nonDefaultCount = 1
            do j=1,numMeshNodes
               if (areDefaultValues(j).eqv..false.) then
                  na(i)%nonDefaultNodes(nonDefaultCount) = j
                  do k=1,na(i)%numVals
                     na(i)%nonDefaultVals(k,nonDefaultCount) = &
                        na(i)%xdmfMatrix(k,j)
                  end do
                  nonDefaultCount = nonDefaultCount + 1
               endif
            end do 
         endif
         exit
      endif
   end do
end do
!-----------------------------------------------------------------------
end subroutine readNodalAttributesXDMF
!-----------------------------------------------------------------------

!----------------------------------------------------------------------
!     S U B R O U T I N E   R E A D   M E S H   X D M F 
!----------------------------------------------------------------------
! Reads the mesh node table, element table, and boundaries from an 
! XDMF file, allocating memory along the way and populating adcirc-style
! data structures.
!----------------------------------------------------------------------
subroutine readMeshXDMF(xdmfFortranObj)
      use CkLunMod, only : CK_LUN
use adcmesh
implicit none
include 'Xdmf.f'
integer*8, intent(in) :: xdmfFortranObj 
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
!$omp threadprivate(gridIndex)
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
integer i, j, k
real(8), allocatable :: tempCoord(:,:)
!
startIndex = 0
arrayStride = 1
valueStride = 1
openMaps = 1
openAttributes = 1
openInformations = 1
openSets = 1
!
!  A G R I D
!
call xdmfRetrieveNumDomainGridCollections(xdmfFortranObj, numGridCollections)
if (numGridCollections.gt.0) then
   write(6+CK_LUN,'(a)') 'INFO: Opening the first Grid in the GridCollection.' 
   call xdmfOpenGridCollectionGrid(xdmfFortranObj, XDMF_GRID_TYPE_UNSTRUCTURED,  &
      gridIndex, openMaps, openAttributes, openInformations, openSets)
   write(6+CK_LUN,'(a)') 'INFO: Reading the Grid name.'
   call xdmfRetrieveGridCollectionGridName(xdmfFortranObj, XDMF_GRID_TYPE_UNSTRUCTURED, &
      gridIndex, gridName, nameLength)
else
   ! there are no grid collections, there must just be a single grid ... open it
   call xdmfOpenDomainGrid(xdmfFortranObj, XDMF_GRID_TYPE_UNSTRUCTURED, &
      gridIndex, openMaps, openAttributes, openInformations, openSets)
   call xdmfRetrieveDomainGridName(xdmfFortranObj, XDMF_GRID_TYPE_UNSTRUCTURED, &
      gridIndex, gridName, nameLength)    
endif
call replaceNullsWithSpaces(gridName)
agrid(1:80) = gridName(1:80)

!
!  S I Z E S 
!
call xdmfRetrieveGeometryNumPoints(xdmfFortranObj, np)
write(6,'("INFO: The Geometry contains ",i0," points (nodes).")') np
call xdmfRetrieveTopologyNumElements(xdmfFortranObj, ne)
write(6,'("INFO: The Topology contains ",i0," elements.")') ne
!
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
   write(6+CK_LUN,'("WARNING: Unrecognized geometry type ",i0,".")') geometryType
end select
!write(6,'("INFO: The geometry type is ",a,".")') trim(geometryTypeString)
!
call xdmfRetrieveGeometryValueType(xdmfFortranObj, geometryDataType)
call createDataTypeString(geometryDataType, geometryDataTypeString)
!write(6,'(a)') 'INFO: The data type of the coordinates in this geometry is ' // trim(geometryDataTypeString) // '.'
!
call xdmfRetrieveGeometrySize(xdmfFortranObj, numContained)
!
allocate(tempCoord(2,np))
tempCoord = -99999.d0
write(6,'(a)') 'INFO: Reading nodal table.'
call xdmfRetrieveGeometryValues(xdmfFortranObj, tempCoord, & 
   geometryDataType, 2*np, startIndex, arrayStride, valueStride)
!
xyd(1:2,:) = tempCoord(1:2,:)
deallocate(tempCoord)
!
! depth
call xdmfRetrieveNumAttributes(xdmfFortranObj, numAttributes)
write(6,'("INFO: Grid ",i0," contains ",i0," attributes.")') gridIndex, numAttributes
do attributeIndex=0, numAttributes - 1
   call xdmfRetrieveAttributeName(xdmfFortranObj, attributeIndex,  & 
      itemName, nameLength)
   call replaceNullsWithSpaces(itemName)
   write(6+CK_LUN,'("INFO: Grid ",i0," Attribute ",i0," is named ",a)') gridIndex, attributeIndex, trim(itemName)
   select case(trim(itemName))
   case("depth")
      call xdmfRetrieveAttributeValues(xdmfFortranObj, attributeIndex, & 
      xyd(3,:), XDMF_ARRAY_TYPE_FLOAT64, np, startIndex, arrayStride, valueStride)
   case default
      ! this is not the attribute we're looking for, at least not yet
   end select
end do
!
if (verbose.eqv..true.) then
   do i=1,np
      write(6+CK_LUN,'("ECHO: node=",i0," x=",f15.7," y=",f15.7," depth=",f8.3)') &
         i, (xyd(j,i), j=1,3)
   end do
end if
! 
!  E L E M E N T   T A B L E
!
!call xdmfRetrieveTopologyNumProperties(xdmfFortranObj, numContained)
!print *, "Number of  Topology Properties: ", numTopologyProperties
!do topologyPropertyIndex=0, numTopologyProperties-1 
!   call xdmfRetrieveTopologyProperty(xdmfFortranObj, topologyPropertyIndex, itemKey,   keyLength, itemValue, valueLength)
!    print *, "Key: ", itemKey
!    print *, "Value: ", itemValue
! end do
!
write(6,'(a)') 'INFO: Reading element table.'
call xdmfRetrieveTopologyType(xdmfFortranObj, topologyType)
! XDMF supports other topology types, but these are unlikely to
! be found in ADCIRC output files
call createTopologyTypeString(topologyType, topologyTypeString)
!write(6,'("INFO: The topology type is ",a,".")') trim(topologyTypeString)
!
call xdmfRetrieveTopologyValueType(xdmfFortranObj, topologyDataType)
call createDataTypeString(topologyDataType, topologyDataTypeString)
!write(6,'(a)') 'INFO: The data type of the topology is ' // trim(topologyDataTypeString) // '.'
!
call xdmfRetrieveTopologySize(xdmfFortranObj, numElementValues)
allocate(xdmf_nm(3,ne))
call xdmfRetrieveTopologyValues(xdmfFortranObj, xdmf_nm, &
   XDMF_ARRAY_TYPE_INT32, numElementValues, startIndex, arrayStride, valueStride) 
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
if (verbose.eqv..true.) then
   do i=1,ne
      write(6+CK_LUN,'("ECHO: element ",i0," nodes ",i0," ",i0," ",i0)') i, (nm(i,j), j=1,3)
   end do
endif
!
! B O U N D A R I E S
!
! boundaries are saved as sets and attributes
write(6,'(a)') 'INFO: Reading boundaries.'
CALL xdmfRetrieveNumSets(xdmfFortranObj, numSets)
write(6,'("INFO: The Grid contains ",i0," boundaries.")') numSets
!
! allocate xdmf-specific arrays to hold boundary parameters
allocate(elevationBoundary(0:numSets-1))
allocate(setSize(0:numSets-1))
allocate(firstSetAttributeIndex(0:numSets-1))
allocate(boundaryTypes(0:numSets-1))
!
! count the different types of boundaries for use in memory allocation
!
!write(6,'("DEBUG: Counting the various boundary types.")') 
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
   !write(6,'("INFO: Opening set ",i0,".")') setIndex
   call xdmfOpenSet(xdmfFortranObj, setIndex, openAttributes, openInformations) 
   !
   ! get the boundary type
   call xdmfRetrieveNumInformation(xdmfFortranObj, numInformations)
   call xdmfRetrieveInformation(xdmfFortranObj, numInformations-1, itemKey, & 
      keyLength, itemValue, valueLength)
   call replaceNullsWithSpaces(itemValue)
   read(itemValue+CK_LUN,*) boundaryTypes(setIndex) ! value of either ibtype_orig or ibtypee
   !
   ! determine the number of nodes on this boundary
   call xdmfRetrieveSetSize(xdmfFortranObj, setSize(setIndex), setIndex)
   !write(6,'("INFO: Set ",i0," contains ",i0," values.")') setIndex, setSize(setIndex)
   !            
   call xdmfRetrieveSetNumProperties(xdmfFortranObj, setIndex, numSetProperties)
   do setPropertyIndex=0, numSetProperties-1
      call xdmfRetrieveSetProperty(xdmfFortranObj, setIndex, setPropertyIndex, itemKey, keyLength, itemValue, valueLength)
      call replaceNullsWithSpaces(itemValue)
      select case(trim(itemValue))
      case("elevation_specified_boundary ")
         !write(6,'("ECHO: The set property is ",a,".")') '"' // trim(itemValue) // '"'
         !write(6,'("DEBUG: Found one elevation specified boundary.")')
         nope = nope + 1
         neta = neta + setSize(setIndex)
         elevationBoundary(setIndex) = .true.
      case("flux_specified_boundary")
         !write(6,'("ECHO: The set property is ",a,".")') '"' // trim(itemValue) // '"'
         !write(6,'("DEBUG: Found one flux specified boundary.")')
         nbou = nbou + 1
         nvel = nvel + setSize(setIndex)
         select case(boundaryTypes(setIndex))
         case(0,1,2,10,11,12,20,21,22,30,52)
            numSimpleFluxBoundaries = numSimpleFluxBoundaries + 1
         case(3,13,23)
            numExternalFluxBoundaries = numExternalFluxBoundaries + 1
         case(4,24)
            numInternalFluxBoundaries = numInternalFluxBoundaries + 1
         case(5,25)
            numInternalFluxBoundariesWithPipes = numInternalFluxBoundariesWithPipes + 1
         case default
            write(6+CK_LUN,'("ERROR: File contains IBTYPE=",i0," which is not a valid flux boundary type.")'), boundaryTypes(setIndex)
         end select         
      case("Node")
         ! do nothing, this property simply indicates that the boundaries
         ! are defined by lists of nodes
      case default
         write(6+CK_LUN,'("WARNING: Unrecognized set property ",a,".")') trim(itemValue)
      end select
   end do
end do     
! 
! Now that we know how many of each boundary type we have, we can 
! allocate memory to hold the data/parameters for each boundary of 
! each type
write(6,'("INFO: Number of elevation boundaries : ",i0)') nope
write(6,'("INFO: Number of elevation boundary nodes : ",i0)') neta
write(6,'("INFO: Total number of flux boundaries : ",i0)') nbou
write(6,'("INFO: Total number of flux boundary nodes : ",i0)') nvel

write(6,'("INFO: Number of simple flux boundaries : ",i0)') numSimpleFluxBoundaries   
write(6,'("INFO: Number of external flux boundaries : ",i0)') numExternalFluxBoundaries   
write(6,'("INFO: Number of internal flux boundaries : ",i0)') numInternalFluxBoundaries
write(6,'("INFO: Number of internal flux boundaries with cross barrier pipes : ",i0)') numInternalFluxBoundariesWithPipes  
!
! populate nvdll (adcirc array representing number of nodes on each 
! elevation boundary segment) and nvell (adcirc array representing
! number of nodes on each flux boundary segment) as these variables 
! are used in the allocateBoundaryArrays() subroutine
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
   !write(6,'("INFO: The data type of set ",i0," is ",a,".")') setIndex, trim(setDataTypeString)
   !
   ! elevation boundary
   if (elevationBoundary(setIndex).eqv..true.) then
      ! get the node numbers on the boundary
      elevationBoundaries(elevCount)%indexNum = elevCount
      call xdmfRetrieveSetValues(xdmfFortranObj, setIndex, & 
         elevationBoundaries(elevCount)%nodes, setDataType, setSize(setIndex), &
         startIndex, arrayStride, valueStride)
      ! convert node numbers from 0 starting index (xdmf-style) to 
      ! 1 starting index (fortran-style)
      elevationBoundaries(elevCount)%nodes = elevationBoundaries(elevCount)%nodes + 1
      ! populate adcirc-native array
      nbdv(elevCount,1:nvdll(elevCount)) = elevationBoundaries(elevCount)%nodes
      elevCount = elevCount + 1
   else
      !
      ! flux boundary type
      !write(6,'("ECHO: The boundary type of flux boundary ",i0," is ",i0,".")') fluxCount, ibtype_orig(fluxCount)
      select case(ibtype_orig(fluxCount))
      case(0,1,2,10,11,12,20,21,22,30,52)
         !write(6,'("DEBUG: No attributes for this type of boundary.")')
         simpleFluxBoundaries(sfCount)%indexNum = fluxCount
         ! get node numbers on the boundary
         call xdmfRetrieveSetValues(xdmfFortranObj, setIndex, & 
            simpleFluxBoundaries(sfCount)%nodes, setDataType, setSize(setIndex), &
            startIndex, arrayStride, valueStride)
         ! convert node numbers from 0 starting index (xdmf-style) to 
         ! 1 starting index (fortran-style)
         simpleFluxBoundaries(sfCount)%nodes = simpleFluxBoundaries(sfCount)%nodes + 1
         ! populate adcirc-native array
         nbvv(fluxCount,1:nvell(fluxCount)) = simpleFluxBoundaries(sfCount)%nodes  
         sfCount = sfCount + 1           
      case(3,13,23)
         externalFluxBoundaries(efCount)%indexNum = fluxCount
         ! get the node numbers on the boundary
         call xdmfRetrieveSetValues(xdmfFortranObj, setIndex, & 
            externalFluxBoundaries(efCount)%nodes, setDataType, setSize(setIndex), &
            startIndex, arrayStride, valueStride)
         ! convert node numbers from 0 starting index (xdmf-style) to 
         ! 1 starting index (fortran-style)
         externalFluxBoundaries(efCount)%nodes = externalFluxBoundaries(efCount)%nodes + 1
         !write(6,'("DEBUG: Attributes of this boundary:")')
         do i=attStart,attStart+1 
            call xdmfRetrieveAttributeName(xdmfFortranObj, i, itemName, nameLength)
            call replaceNullsWithSpaces(itemName)
            select case(trim(itemName))
            case("BARLANHT")                ! barrier height at each node
               call xdmfRetrieveAttributeValues(xdmfFortranObj, i, & 
                  externalFluxBoundaries(efCount)%barlanht, XDMF_ARRAY_TYPE_FLOAT64, &
                  nvell(fluxCount), startIndex, arrayStride, valueStride)
            case("BARLANCFSP")              ! coefficient of free surface super critical flow
               call xdmfRetrieveAttributeValues(xdmfFortranObj, i, & 
                  externalFluxBoundaries(efCount)%barlancfsp, XDMF_ARRAY_TYPE_FLOAT64, &
                  nvell(fluxCount), startIndex, arrayStride, valueStride)
            case default
               write(6+CK_LUN,'("ERROR: Unrecognized boundary attribute : ",a,".")') trim(itemName)
            end select
         end do
         ! populate adcirc-native arrays
         nbvv(fluxCount,1:nvell(fluxCount)) = externalFluxBoundaries(efCount)%nodes
         barlanht(fluxCount,1:nvell(fluxCount)) = externalFluxBoundaries(efCount)%barlanht
         barlancfsp(fluxCount,1:nvell(fluxCount)) = externalFluxBoundaries(efCount)%barlancfsp
         efCount = efCount + 1
      case(4,24)  ! internal barrier boundary (e.g., subgrid scale levee)
         internalFluxBoundaries(ifCount)%indexNum = fluxCount
         ! get the node numbers on the boundary
         call xdmfRetrieveSetValues(xdmfFortranObj, setIndex, & 
            internalFluxBoundaries(ifCount)%nodes, setDataType, setSize(setIndex), &
            startIndex, arrayStride, valueStride)
         ! convert node numbers from 0 starting index (xdmf-style) to 
         ! 1 starting index (fortran-style)
         internalFluxBoundaries(ifCount)%nodes = internalFluxBoundaries(ifCount)%nodes + 1
         !write(6,'("DEBUG: Attributes of this boundary:")')
         do i=attStart,attStart+3
            call xdmfRetrieveAttributeName(xdmfFortranObj, i, itemName, nameLength)
            call replaceNullsWithSpaces(itemName)
            select case(trim(itemName))
            case("IBCONN")      ! paired (back face) nodes
               call xdmfRetrieveAttributeValues(xdmfFortranObj, i, & 
                  internalFluxBoundaries(ifCount)%ibconn, XDMF_ARRAY_TYPE_INT32, &
                  nvell(fluxCount), startIndex, arrayStride, valueStride)
            case("BARINHT")     ! barrier height at each node and its paired node
               call xdmfRetrieveAttributeValues(xdmfFortranObj, i, & 
                  internalFluxBoundaries(ifCount)%barinht, XDMF_ARRAY_TYPE_FLOAT64, &
                  nvell(fluxCount), startIndex, arrayStride, valueStride)
            case("BARINCFSB")   ! coefficient of free surface sub critical flow
               call xdmfRetrieveAttributeValues(xdmfFortranObj, i, & 
                  internalFluxBoundaries(ifCount)%barincfsb, XDMF_ARRAY_TYPE_FLOAT64, &
                  nvell(fluxCount), startIndex, arrayStride, valueStride)
            case("BARINCFSP")   ! coefficient of free surface super critical flow
               call xdmfRetrieveAttributeValues(xdmfFortranObj, i, & 
                  internalFluxBoundaries(ifCount)%barincfsp, XDMF_ARRAY_TYPE_FLOAT64, &
                  nvell(fluxCount), startIndex, arrayStride, valueStride)
            case default
               write(6+CK_LUN,'("ERROR: Unrecognized boundary attribute : ",a,".")') trim(itemName)
            end select
         end do
         internalFluxBoundaries(ifCount)%ibconn = internalFluxBoundaries(ifCount)%ibconn + 1
         ! populate adcirc-native arrays
         nbvv(fluxCount,1:nvell(fluxCount)) = internalFluxBoundaries(ifCount)%nodes  
         ibconn(fluxCount,1:nvell(fluxCount)) = internalFluxBoundaries(ifCount)%ibconn         
         barinht(fluxCount,1:nvell(fluxCount)) = internalFluxBoundaries(ifCount)%barinht
         barincfsb(fluxCount,1:nvell(fluxCount)) = internalFluxBoundaries(ifCount)%barincfsb
         barincfsp(fluxCount,1:nvell(fluxCount)) = internalFluxBoundaries(ifCount)%barincfsp
         ifCount = ifCount + 1
      case(5,25)  ! internal barrier boundary with cross barrier pipes
         internalFluxBoundariesWithPipes(ifwpCount)%indexNum = fluxCount
         ! get the node numbers on the boundary
         call xdmfRetrieveSetValues(xdmfFortranObj, setIndex, & 
            internalFluxBoundariesWithPipes(ifwpCount)%nodes, setDataType, setSize(setIndex), &
            startIndex, arrayStride, valueStride)
         ! convert node numbers from 0 starting index (xdmf-style) to 
         ! 1 starting index (fortran-style)
         internalFluxBoundariesWithPipes(ifwpCount)%nodes = internalFluxBoundaries(ifwpCount)%nodes + 1
         !write(6,'("DEBUG: Attributes of this boundary:")')
         do i=attStart,attStart+6
            call xdmfRetrieveAttributeName(xdmfFortranObj, i, itemName, nameLength)
            call replaceNullsWithSpaces(itemName)
            select case(trim(itemName))
            case("IBCONN")  ! paired (i.e., back face) nodes
               call xdmfRetrieveAttributeValues(xdmfFortranObj, i, & 
                  internalFluxBoundariesWithPipes(ifwpCount)%ibconn, XDMF_ARRAY_TYPE_INT32, &
                  nvell(fluxCount), startIndex, arrayStride, valueStride)
            case("BARINHT") ! barrier height at each node and its paired node
               call xdmfRetrieveAttributeValues(xdmfFortranObj, i, & 
                  internalFluxBoundariesWithPipes(ifwpCount)%barinht, XDMF_ARRAY_TYPE_FLOAT64, &
                  nvell(fluxCount), startIndex, arrayStride, valueStride)
            case("BARINCFSB") ! coefficient of free surface sub critical flow
               call xdmfRetrieveAttributeValues(xdmfFortranObj, i, & 
                  internalFluxBoundariesWithPipes(ifwpCount)%barincfsb, &
                  XDMF_ARRAY_TYPE_FLOAT64, & 
                  nvell(fluxCount), startIndex, arrayStride, valueStride)
            case("BARINCFSP") ! coefficient of free surface super critical flow
               call xdmfRetrieveAttributeValues(xdmfFortranObj, i, & 
                  internalFluxBoundariesWithPipes(ifwpCount)%barincfsp, &
                  XDMF_ARRAY_TYPE_FLOAT64, &
                  nvell(fluxCount), startIndex, arrayStride, valueStride)
            case("PIPEHT")    ! barrier height at each node and its paired node
               call xdmfRetrieveAttributeValues(xdmfFortranObj, i, & 
                  internalFluxBoundariesWithPipes(ifwpCount)%pipeht, XDMF_ARRAY_TYPE_FLOAT64, &
                  nvell(fluxCount), startIndex, arrayStride, valueStride)
            case("PIPECOEF")   ! coefficient of free surface sub critical flow
               call xdmfRetrieveAttributeValues(xdmfFortranObj, i, & 
                  internalFluxBoundariesWithPipes(ifwpCount)%pipecoef, &
                  XDMF_ARRAY_TYPE_FLOAT64, & 
                  nvell(fluxCount), startIndex, arrayStride, valueStride)
            case("PIPEDIAM")   ! coefficient of free surface super critical flow
               call xdmfRetrieveAttributeValues(xdmfFortranObj, i, & 
                  internalFluxBoundariesWithPipes(ifwpCount)%pipediam, &
                  XDMF_ARRAY_TYPE_FLOAT64, &
                  nvell(fluxCount), startIndex, arrayStride, valueStride)
            case default
               write(6+CK_LUN,'("ERROR: Unrecognized boundary attribute : ",a,".")') trim(itemName)
            end select
         end do
         internalFluxBoundariesWithPipes(ifwpCount)%ibconn = internalFluxBoundariesWithPipes(ifwpCount)%ibconn + 1
         ! populate adcirc-native arrays       
         nbvv(fluxCount,1:nvell(fluxCount)) &
            = internalFluxBoundariesWithPipes(ifwpCount)%nodes
         ibconn(fluxCount,1:nvell(fluxCount)) &
            = internalFluxBoundariesWithPipes(ifwpCount)%ibconn
         barinht(fluxCount,1:nvell(fluxCount)) &
           = internalFluxBoundariesWithPipes(ifwpCount)%barinht
         barincfsb(fluxCount,1:nvell(fluxCount)) &
            = internalFluxBoundariesWithPipes(ifwpCount)%barincfsb
         barincfsp(fluxCount,1:nvell(fluxCount)) &
            = internalFluxBoundariesWithPipes(ifwpCount)%barincfsp
         pipeht(fluxCount,1:nvell(fluxCount)) &
            = internalFluxBoundariesWithPipes(ifwpCount)%pipeht
         pipecoef(fluxCount,1:nvell(fluxCount)) &
            = internalFluxBoundariesWithPipes(ifwpCount)%pipecoef
         pipediam(fluxCount,1:nvell(fluxCount)) &
            = internalFluxBoundariesWithPipes(ifwpCount)%pipediam
         ifwpCount = ifwpCount + 1            
      case default
          write(6+CK_LUN,'("ERROR: The boundary type ",I3," was found in the files but is not valid.")')
          stop
      end select
      fluxCount = fluxCount + 1
   endif
end do ! loop over boundaries
!----------------------------------------------------------------------
end subroutine readMeshXDMF
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!                   S U B R O U T I N E   
! G E T   A T T R I B U T E   C H A R A C T E R I S T I C S   X D M F
!----------------------------------------------------------------------
! Get info about an attribute.
!----------------------------------------------------------------------
subroutine getAttributeCharacteristicsXDMF(xdmfFortranObj, attributeIndex)
implicit none
include 'Xdmf.f'
integer*8, intent(in) :: xdmfFortranObj
integer, intent(in) :: attributeIndex
!
integer :: typeHolder
character(len=256) :: logString
!
call xdmfRetrieveAttributeType(xdmfFortranObj, attributeIndex, typeHolder)
call createAttributeTypeString(typeHolder, logString) 
write(6,'("DEBUG: The Attribute type is ",a,".")') trim(logString)
!
call xdmfRetrieveAttributeCenter(xdmfFortranObj, attributeIndex, typeHolder)
call createAttributeCenterString(typeHolder, logString)
write(6,'("DEBUG: The Attribute Center is ",a,".")') trim(logString)
!
call xdmfRetrieveAttributeValueType(xdmfFortranObj, attributeIndex, typeHolder)
call createDataTypeString(typeHolder, logString)
write(6,'("DEBUG: The Attribute data type is ",a,".")') trim(logString)
!
call xdmfRetrieveAttributeSize(xdmfFortranObj, attributeIndex, typeHolder)
write(6,'("DEBUG: The Attribute consists of ",i0," values.")') typeHolder
!----------------------------------------------------------------------
end subroutine getAttributeCharacteristicsXDMF
!----------------------------------------------------------------------


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
do 
   nullCharLocation = index(myString,char(0)) 
   if (nullCharLocation.ne.0) then
      myString(nullCharLocation:nullCharLocation) = ' ' 
   else
      exit ! there are no more null characters
   endif
end do
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
   write(6+CK_LUN,'("WARNING: Unrecognized attribute type ",i0,".")') trim(typeString)
end select
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
   write(6+CK_LUN,'("WARNING: Unrecognized set type ",i0,".")') typeHolder
end select
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
   write(6+CK_LUN,'("WARNING: Unrecognized attribute center ",i0,".")') trim(typeString)
end select
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
   write(6+CK_LUN,'("WARNING: Unrecognized data type ",i0,".")') trim(typeString)
end select
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
   write(6+CK_LUN,'("WARNING: Unrecognized topology type ",i0,".")') typeHolder
end select
!----------------------------------------------------------------------
end subroutine createTopologyTypeString
!----------------------------------------------------------------------
