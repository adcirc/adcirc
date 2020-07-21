#ifdef ADCXDMF
C----------------------------------------------------------------------
C                     M O D U L E     X D M F I O
C----------------------------------------------------------------------
C     jgf51.21.28: Created to hold XDMF-related output routines 
C     where they can be called by the writer module as well as the
C     write output module.
C----------------------------------------------------------------------
      module xdmfio
C----------------------------------------------------------------------
      use sizes, only : sz
      logical :: release = .true.  ! .true. to release memory after writing data to HDF5
      logical :: writeToHDF5 = .true. ! .true. if XdmfAddGrid should immediately write mesh
      integer :: lightDataLimit = 10 ! max values to save as light data (xml)
      logical :: meshInitialized = .false. ! needed by dedicated writer processors

      contains
C-----------------------------------------------------------------------
C     S U B R O U T I N E   I N I T  O U T P U T   X D M F
C-----------------------------------------------------------------------
C     Writes the ADCIRC mesh to the XDMF output file.
C-----------------------------------------------------------------------
      subroutine initOutputXDMF(descript)
      use global
      use mesh
      use sizes,only:XDMF
      implicit none
      include 'adcirc_Xdmf.f'
      type(OutputDataDescript_t) :: descript
      integer :: informationID
      integer :: i, j
      !
      ! return early if initialization of XDMF file is not required 
      ! for this output file type
      if (abs(descript%specifier).ne.XDMF) then
         return
      endif
#ifndef ADCXDMF
      call allMessage(ERROR,'XDMF output was specified.')
      call allMessage(ERROR,'This executable was not compiled with XDMF support.')
      call terminate()
#else
      !jgfdebug
      write(6,'(a)') 'fn='//trim(descript % file_name)
      write(6,'(a,i0)') 'lun=',descript % lun
      !
      ! read ADCIRC control data from ascii file unless otherwise specified
      !call readControlFile(controlFileName)
      !if (nabout.eq.-1) then
      !      write(6,'(a)') 'INFO: Echoing control file.'
      !      open(25,file='echo.15',status='replace')
      !      call echoControlFile(25)
      !      close(25)
      !      write(6,'(a)') 'INFO: Finished echoing control file.'!
      !endif
      ! 
      ! create xdmf file
      !
      ! call the XdmfFortran constructor; return a reference to the XdmfFortran object
      ! (the pointer is packaged as a long int in Fortran)
      write(6,'(A)') 'INFO: Initializing XDMF object.'
      call xdmfInit(descript % xdmfMD % xdmfFortranObj)
      !
      call xdmfInitHDF5(descript % xdmfMD % xdmfFortranObj,
     &    trim(descript % file_name) //'.h5'//char(0), release)
      !
      ! name the temporal grid collection for time varying output data
      informationID = XdmfAddInformation(descript % xdmfMD % xdmfFortranObj, 
     &      'TimeVaryingOutputData'//CHAR(0), trim(adjustl(agrid))//CHAR(0))
      !
      ! write fort.15 data to the xml file 
      !call writeControlXDMF(xdmfFortranObj)
      !
      ! create a grid collection for time varying data 
      call xdmfAddGridCollection(descript % xdmfMD % xdmfFortranObj, 
     &        "Temporal"//CHAR(0), XDMF_GRID_COLLECTION_TYPE_TEMPORAL)


#endif
      !----------------------------------------------------------------
      end subroutine initOutputXDMF
      !----------------------------------------------------------------

      !----------------------------------------------------------------
      !    S U B R O U T I N E   W R I T E   O U T   A R R A Y   X D M F
      !----------------------------------------------------------------
      ! Writes an array to an existing XDMF file. 
      !----------------------------------------------------------------
      subroutine writeOutArrayXDMF(timeLoc, it, descript)
      use global, only: OutputDataDescript_t, scratchmessage, error, nt,
     &            allMessage
      use mesh, only: NM, NE, agrid, dp, ics, np, x, y, terminate
      use sizes,only: mnproc
      implicit none
#ifdef ADCXDMF
      include 'adcirc_Xdmf.f'
#endif
      type(OutputDataDescript_t) :: descript
      real(8), intent(in) :: timeLoc
      integer, intent(in) :: it
      integer :: attributeID
      integer :: attributeType
      integer :: numValues
      integer :: timeStepID
      integer :: gridCollectionIndex ! the grid collections are numbered starting from zero
      integer :: numGrids
      !
      ! the following four integers act like logicals (1=true, 0=false), and 
      ! they control whether the associated data should be opened when a grid 
      ! collection is opened  
      integer :: openMaps ! 1 if Maps should be opened within another object
      integer :: openAttributes ! 1 if Attributes should be opened within another object 
      integer :: openInformations ! 1 if Informations should be opened within another object
      integer :: openSets ! 1 if sets should be opened within another object     
      character(len=80) :: timeStepString
      logical :: lastOutput
      integer :: timeStepsRemaining
      integer :: i, j
#ifndef ADCXDMF
      call allMessage(ERROR,'XDMF output was specified.')
      call allMessage(ERROR,
     &    'This executable was not compiled with XDMF support.')
      call terminate()
#else
      !jgfdebug
      write(6,'(a)') 'call to write '//trim(descript % file_name)
      ! 
      ! if this is the first Grid to be added to the collection, then
      ! write the mesh and boundary data to the file ... if this is not
      ! the first Grid, then just add references to the mesh and boundaries
      call xdmfRetrieveGridCollectionNumGrids(
     &   descript % xdmfMD % xdmfFortranObj, 
     &   XDMF_GRID_TYPE_UNSTRUCTURED,  numGrids)
      if (numGrids.eq.0) then    
         !
         ! create a 0-offset element table, i.e., one that assumes the nodes
         ! are numbered starting from zero ... XDMF2 is zero offset while the 
         ! ADCIRC ascii format assumes nodes are numbered starting from 1
         allocate(descript % xdmfMD % xdmf_nm(3,ne))
         do i=1,ne
            do j=1,3
               descript % xdmfMD % xdmf_nm(j,i) = nm(i,j) - 1 ! store in row major format
            end do
         end do
         !
         ! call the setTopology method; arguments include the type of the topology, 
         ! the number of values in the connectivity array, the numeric type of
         ! the connectivity array, and the connectivity values with a 0 offset;
         ! the call returns the topology ID in case it will be reused. The XDMF
         ! library interprets the connectivity array as a series of contiguous 
         ! values (i.e., as a 1D array). 
         descript % xdmfMD % topologyID = 
     &       xdmfSetTopology(descript % xdmfMD % xdmfFortranObj, 
     &       XDMF_TOPOLOGY_TYPE_TRIANGLE, ne*3, XDMF_ARRAY_TYPE_INT32, 
     &       descript % xdmfMD % xdmf_nm)
         !
         ! call the setGeometry method; arguments include the geometry type, the 
         ! number of values in the coordinates array, the numeric type of the coordinates,
         ! and the coordinates array (again, interpreted by XDMF as a series of 
         ! contiguous values). 
         allocate(descript % xdmfMD % tempCoord(2,np))
         descript % xdmfMD % tempCoord(1,:) = x
         descript % xdmfMD % tempCoord(2,:) = y
         descript % xdmfMD % geometryID = xdmfSetGeometry(
     &      descript % xdmfMD % xdmfFortranObj, XDMF_GEOMETRY_TYPE_XY,
     &      np*2, XDMF_ARRAY_TYPE_FLOAT64, descript % xdmfMD % tempCoord)
         !
         ! add the mesh boundaries to the grid
         call addBoundaries(descript)
         !
         ! add bathymetric depth
         write(6,'(a)') 'INFO: Writing bathy/topo to XDMF file.'
         ! write the bathymetric depth
         !descript % xdmfmd % createdIDs = .false.
         !descript % xdmfmd % variable_name = 'depth'
         !descript % xdmfmd % long_name = 'distance from geoid'
         !descript % xdmfmd % standard_name = 'depth_below_geoid'
         !descript % xdmfmd % coordinates = 'y x'
         !descript % xdmfmd % units = 'm'
         !descript % xdmfmd % positive = 'downward'
         !call writeMetaData(xdmfFortranObj, md)
         ! write projection info 
         !md%variable_name_id = XdmfAddInformation(xdmfFortranObj, 
         ! &   'projection'//CHAR(0), trim(projection)//CHAR(0))
         descript % xdmfMD % depthID = XdmfAddAttribute(
     &      descript % xdmfMD % xdmfFortranObj, 
     &      'depth'//CHAR(0), 
     &      XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, np, 
     &      XDMF_ARRAY_TYPE_FLOAT64, dp)
         write(6,'(a)') 'INFO: Finished writing bathy/topo to XDMF file.'
      else
         ! we've already added Grids to the collection ... don't write
         ! the mesh and boundary data again, just create references to it 
         !     
         !call xdmfSetPreviousTopology(xdmfFortranObj, topologyID)
         !call xdmfSetPreviousGeometry(xdmfFortranObj, geometryID)
         call addBoundaryReferences(descript)
         call xdmfAddPreviousAttribute(descript % xdmfMD % xdmfFortranObj,
     &                                  descript % xdmfMD % depthID)
      endif
      !
      ! set up metadata for each data type
      descript % xdmfMD % createdIDs = .false.
      descript % xdmfMD % positive = 'null'
      select case(descript % lun)
      case(63)
         descript % xdmfMD % variable_name = 'zeta'
         descript % xdmfMD % long_name = 'water surface elevation above geoid'
         descript % xdmfMD % standard_name = 'water_surface_elevation'
         descript % xdmfMD % coordinates = 'time y x'
         descript % xdmfMD % units = 'm'
      case(64)
         descript % xdmfMD % variable_name = 'vel'
         descript % xdmfMD % long_name = 'water column vertically averaged velocity'
         descript % xdmfMD % standard_name = 'water_velocity'
         descript % xdmfMD % coordinates = 'time y x'
         descript % xdmfMD % units = 'm s-1'
         descript % xdmfMD % positive = 'north/east'
      case(73)
         descript % xdmfMD % variable_name = 'pressure'
         descript % xdmfMD % long_name = 'air pressure at sea level'
         descript % xdmfMD % standard_name = 'air_pressure_at_sea_level'
         descript % xdmfMD % coordinates = 'time y x'
         descript % xdmfMD % units = 'meters of water'
      case(74)
         descript % xdmfMD % variable_name = 'wvel'
         descript % xdmfMD % long_name = 'wind_velocity'
         descript % xdmfMD % standard_name = 'wind'
         descript % xdmfMD % coordinates = 'time x y'
         descript % xdmfMD % units = 'm s-1'
         descript % xdmfMD % positive = 'north/east'
      case(93)
         descript % xdmfMD % variable_name = 'ice'
         descript % xdmfMD % long_name = 'ice coverage at at sea surface'
         descript % xdmfMD % standard_name = 'ice_pressure_at_sea_level'
         descript % xdmfMD % coordinates = 'time y x'
         descript % xdmfMD % units = 'percent'
      case default
         write(scratchMessage,'("Cannot write unit ",i0," files.")') descript % lun 
         call allMessage(ERROR, scratchMessage)
         call terminate()
      end select
      !
      ! set the time of this dataset in the XDMF file
      call XdmfSetTime(descript % xdmfMD % xdmfFortranObj, timeLoc)
      !
      ! set the metadata for this dataset
      !call writeMetaData(xdmfFortranObj, md)
      !
      select case(descript % num_items_per_record)
      case(1) ! scalar data
         attributeType = XDMF_ATTRIBUTE_TYPE_SCALAR
         if (mnproc.gt.1) then
            numValues = descript % num_fd_records
            attributeID = XdmfAddAttribute(
     &         descript % xdmfMD % xdmfFortranObj, 
     &         trim(descript % xdmfMD % variable_name)//CHAR(0), 
     &         XDMF_ATTRIBUTE_CENTER_NODE, attributeType, numValues, 
     &         XDMF_ARRAY_TYPE_FLOAT64, descript % array_g)
         else
            numValues = descript % num_records_this
            attributeID = XdmfAddAttribute(
     &         descript % xdmfMD % xdmfFortranObj, 
     &         trim(descript % xdmfMD % variable_name)//CHAR(0), 
     &         XDMF_ATTRIBUTE_CENTER_NODE, attributeType, numValues, 
     &         XDMF_ARRAY_TYPE_FLOAT64, descript % array)
         endif        
      case(2) ! vector data
         ! XDMF doesn't understand vectors with two components, so we
         ! have to allocate for three components, and set the z component
         ! of the vector to zero
         attributeType = XDMF_ATTRIBUTE_TYPE_VECTOR
         if (mnproc.gt.1) then
            numValues = 3*descript%num_fd_records
            if (numGrids.eq.0) then
               allocate(descript % xdmfMD 
     &             % data_array3(3,1:descript%num_fd_records))
            endif
            descript % xdmfMD % data_array3(1,:) = descript % array_g
            descript % xdmfMD % data_array3(2,:) = descript % array2_g
         else
            numValues = 3*descript%num_records_this
            if (numGrids.eq.0) then
               allocate(descript % xdmfMD 
     &             % data_array3(3,1:descript%num_records_this))
            endif
            descript % xdmfMD % data_array3(1,:) = descript % array
            descript % xdmfMD % data_array3(2,:) = descript % array2
         endif
         descript % xdmfMD % data_array3(3,:) = 0.d0
         attributeID = XdmfAddAttribute(descript % xdmfMD % xdmfFortranObj, 
     &      trim(descript % xdmfMD % variable_name)//CHAR(0), 
     &      XDMF_ATTRIBUTE_CENTER_NODE, attributeType, numValues, 
     &      XDMF_ARRAY_TYPE_FLOAT64, descript % xdmfMD % data_array3)
      case default
         write(scratchMessage,
     &     '("Cannot write data with ",i0," components.")') 
     &      descript % num_items_per_record
         call allMessage(ERROR,scratchMessage)
         call terminate()
      end select
      !
      ! set the time step of this dataset in the XDMF file
      write(timeStepString,'(i0)') it
      timeStepID = XdmfAddInformation(descript % xdmfMD % xdmfFortranObj,
     & 'IT'//CHAR(0), trim(timeStepString)//CHAR(0))
      ! 
      ! call the addGrid method; creates an unstructured mesh object with the
      ! specified name (2nd arg), then associates the geometry and topology 
      ! created above with this new unstructured mesh, also associates any
      ! informations or attributes with the new mesh, immediately writing
      ! it to the hdf5 file if the last argument is set to .true. 
      call XdmfAddGrid(descript % xdmfMD % xdmfFortranObj, 
     &   trim(agrid)//char(0), writeToHDF5)


      ! jgfdebug
      lastOutput = .false.
      timeStepsRemaining = nt - it
      if ( (nt - descript % endTimeStep) .gt.0 ) then 
         timeStepsRemaining = min(timeStepsRemaining, 
     &                            nt - descript % endTimeStep )
      endif
      write(6,'(a,i0)') 'timesteps remaining=',timeStepsRemaining
      if ( timeStepsRemaining .lt. descript % outputTimeStepIncrement) then
        lastOutput = .true.
      endif
      if ( lastOutput.eqv..true. ) then
         write(6,'(a)') 'finalizing output'
         call finalizeOutputFileXDMF(descript)
      endif
#endif
      !-----------------------------------------------------------------      
      end subroutine writeOutArrayXDMF
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !                   S U B R O U T I N E  
      !       F I N A L I Z E   O U T P U T   F I L E   X D M F
      !-----------------------------------------------------------------
      ! jgf51.21.24: Close up the grid collection and write xml.
      !-----------------------------------------------------------------
      subroutine finalizeOutputFileXDMF(descript)
      use global, only : OutputDataDescript_t
      implicit none
#ifdef ADCXDMF      
      include 'adcirc_Xdmf.f'
      type(OutputDataDescript_t) :: descript   
      !
      ! close this grid collection, writing it to the heavy data (HDF5) file
      ! if the value of writeToHDF5 is .true.; further grids cannot 
      ! be added to this collection
      call XdmfCloseGridCollection(descript % xdmfMD % xdmfFortranObj,
     &     writeToHDF5)
      !
      ! call the write method; arguments include the name of the xml file, 
      ! the max number of light data values, and whether to release memory
      ! after writing. This actually writes both the light and heavy data.
      write(6,'(a)') 'INFO: Writing data to file and releasing memory.'
      call XdmfWrite(descript % xdmfMD % xdmfFortranObj, 
     &   trim(descript % file_name)//'.xmf'//char(0), 
     &   lightDataLimit, release)
      !
      ! call the close method (deletes the XdmfFortran object)
      write(6,'(A)') 'INFO: Cleaning up and deleting XDMF object.'
      call xdmfClose(descript % xdmfMD % xdmfFortranObj)
#endif
      !----------------------------------------------------------------
      end subroutine finalizeOutputFileXDMF      
      !----------------------------------------------------------------
     
!----------------------------------------------------------------------
!       S U B R O U T I N E    W R I T E   M E T A   D A T A 
!----------------------------------------------------------------------
!     Write the metadata for the variable of interest.
!----------------------------------------------------------------------
      subroutine writeMetaData(descript)
      use global
      implicit none
#ifdef ADCXDMF
      include 'adcirc_Xdmf.f' 
      !
      type(OutputDataDescript_t) :: descript 
      !
      if (descript % xdmfMD % createdIDs.eqv..true.) then
         call XdmfAddPreviousInformation(
     &      descript % xdmfMD % xdmfFortranObj, 
     &      descript % xdmfMD % variable_name_id)
         call XdmfAddPreviousInformation(
     &      descript % xdmfMD % xdmfFortranObj, 
     &      descript % xdmfMD % long_name_id)
         call XdmfAddPreviousInformation(
     &      descript % xdmfMD % xdmfFortranObj, 
     &      descript % xdmfMD % standard_name_id)
         call XdmfAddPreviousInformation(
     &      descript % xdmfMD % xdmfFortranObj, 
     &      descript % xdmfMD % coordinates_id)
        call XdmfAddPreviousInformation(
     &   descript % xdmfMD % xdmfFortranObj, 
     &   descript % xdmfMD % units_id)
         if (trim(descript % xdmfMD % positive).ne.'null') then
            call XdmfAddPreviousInformation(
     &         descript % xdmfMD % xdmfFortranObj, 
     &         descript % xdmfMD % positive_id)
         endif
      else
         descript % xdmfMD % variable_name_id = 
     &      XdmfAddInformation(descript % xdmfMD % xdmfFortranObj, 
     &     'variable_name'//CHAR(0), 
     &     trim(descript % xdmfMD % variable_name)//CHAR(0))
         descript % xdmfMD % long_name_id = 
     &      XdmfAddInformation(descript % xdmfMD % xdmfFortranObj, 
     &     'long_name'//CHAR(0), 
     &     trim(descript % xdmfMD % long_name)//CHAR(0))
        descript % xdmfMD % standard_name_id = 
     &     XdmfAddInformation(descript % xdmfMD % xdmfFortranObj, 
     &    'standard_name'//CHAR(0), 
     &     trim(descript % xdmfMD % standard_name)//CHAR(0))
         descript % xdmfMD % coordinates_id = 
     &      XdmfAddInformation(descript % xdmfMD % xdmfFortranObj, 
     &     'coordinates'//CHAR(0), 
     &     trim(descript % xdmfMD % coordinates)//CHAR(0))
         descript % xdmfMD % units_id = 
     &      XdmfAddInformation(descript % xdmfMD % xdmfFortranObj, 
     &     'units'//CHAR(0), 
     &     trim(descript % xdmfMD % units)//CHAR(0))
!         if ( trim(descript % xdmfMD % md%positive).ne.'null') then
!    &      md%positive_id = XdmfAddInformation(descript % xdmfMD % xdmfFortranObj,'positive'//CHAR(0), 
!     &         trim(descript % xdmfMD % md%positive)//CHAR(0))
!         endif
         descript % xdmfMD % createdIDs = .true.
      endif
#endif
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
      subroutine addBoundaries(descript)
      use boundaries
      use global, only : OutputDataDescript_t
      implicit none
#ifdef ADCXDMF
      include 'adcirc_Xdmf.f' 
      !
      type(OutputDataDescript_t) :: descript
      character(1024) :: fluxBoundaryType ! character represenation of IBTYPE
      integer :: i, j  ! loop counter
      !
      write(6,'(A)') 'INFO: Adding boundary data.'
      do i=1, nope 
         ! jgftodo: this is hardcoded to 0 since this is the only type that ADCIRC
         ! supports, and most or all mesh files don't even have this value, although
         ! the documentation calls for it 
         elevationBoundaries(i)%informationID = 
     &     XdmfAddInformation(descript % xdmfMD % xdmfFortranObj, 
     &     'IBTYPEE'//CHAR(0), '0'//CHAR(0))
         elevationBoundaries(i)%setID = 
     &      XdmfAddSet(descript % xdmfMD % xdmfFortranObj, 
     &     'elevation_specified_boundary'//char(0), XDMF_SET_TYPE_NODE, 
     &     elevationBoundaries(i)%xdmf_nodes, nvdll(i), XDMF_ARRAY_TYPE_INT32)
      end do
      sfCount = 1
      efCount = 1
      ifCount = 1
      ifwpCount = 1 
      do i=1, nbou 
         write(fluxBoundaryType,'(i0)') ibtype_orig(i) 
         select case(ibtype_orig(i))
         case(0,1,2,10,11,12,20,21,22,30,32,52)
            simpleFluxBoundaries(sfCount)%informationID = 
     &         XdmfAddInformation(descript % xdmfMD % xdmfFortranObj, 
     &         'IBTYPE'//CHAR(0), 
     &         trim(fluxBoundaryType)//CHAR(0))
            simpleFluxBoundaries(sfCount)%setID = 
     &          XdmfAddSet(descript % xdmfMD % xdmfFortranObj, 
     &          'flux_specified_boundary'//char(0), XDMF_SET_TYPE_NODE, 
     &          simpleFluxBoundaries(sfCount)%xdmf_nodes(:), nvell(i), XDMF_ARRAY_TYPE_INT32)
            sfCount = sfCount + 1
         case(3,13,23)
            externalFluxBoundaries(efCount)%attributeIDs(1) = 
     &          XdmfAddAttribute(
     &          descript % xdmfMD % xdmfFortranObj, 'BARLANHT'//CHAR(0), 
     &          XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), 
     &          XDMF_ARRAY_TYPE_FLOAT64, externalFluxBoundaries(efCount)%barlanht)
            externalFluxBoundaries(efCount)%attributeIDs(2) = 
     &          XdmfAddAttribute(descript % xdmfMD % xdmfFortranObj, 
     &          'BARLANCFSP'//CHAR(0), 
     &          XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), 
     &          XDMF_ARRAY_TYPE_FLOAT64, externalFluxBoundaries(efCount)%barlancfsp)
            externalFluxBoundaries(efCount)%informationID = 
     &          XdmfAddInformation(descript % xdmfMD % xdmfFortranObj,
     &          'IBTYPE'//CHAR(0), 
     &          trim(fluxBoundaryType)//CHAR(0))
            externalFluxBoundaries(efCount)%setID = 
     &           XdmfAddSet(descript % xdmfMD % xdmfFortranObj, 
     &          'flux_specified_boundary'//char(0), XDMF_SET_TYPE_NODE, 
     &          externalFluxBoundaries(efCount)%xdmf_nodes(:), nvell(i), XDMF_ARRAY_TYPE_INT32)
            efCount = efCount + 1
         case(4,24)
            internalFluxBoundaries(ifCount)%attributeIDs(1) = 
     &          XdmfAddAttribute(descript % xdmfMD % xdmfFortranObj, 
     &          'IBCONN'//CHAR(0), 
     &          XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), 
     &          XDMF_ARRAY_TYPE_INT32, internalFluxBoundaries(ifCount)%xdmf_ibconn) ! zero offset
            internalFluxBoundaries(ifCount)%attributeIDs(2) = 
     &          XdmfAddAttribute(descript % xdmfMD % xdmfFortranObj, 
     &          'BARINHT'//CHAR(0), 
     &          XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), 
     &          XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundaries(ifCount)%barinht)
            internalFluxBoundaries(ifCount)%attributeIDs(3) = 
     &          XdmfAddAttribute(descript % xdmfMD % xdmfFortranObj, 
     &          'BARINCFSB'//CHAR(0), 
     &          XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), 
     &          XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundaries(ifCount)%barincfsb)
            internalFluxBoundaries(ifCount)%attributeIDs(4) = 
     &          XdmfAddAttribute(descript % xdmfMD % xdmfFortranObj, 
     &          'BARINCFSP'//CHAR(0), 
     &          XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i),
     &          XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundaries(ifCount)%barincfsp)
            internalFluxBoundaries(ifCount)%informationID = 
     &           XdmfAddInformation(descript % xdmfMD % xdmfFortranObj, 
     &          'IBTYPE'//CHAR(0), trim(fluxBoundaryType)//CHAR(0))
            internalFluxBoundaries(ifCount)%nodes(:) = 
     &         internalFluxBoundaries(ifCount)%nodes(:) - 1
            internalFluxBoundaries(ifCount)%setID = 
     &          XdmfAddSet(descript % xdmfMD % xdmfFortranObj, 
     &          'flux_specified_boundary'//char(0), XDMF_SET_TYPE_NODE, 
     &          internalFluxBoundaries(ifCount)%xdmf_nodes, nvell(i), XDMF_ARRAY_TYPE_INT32)
            ifCount = ifCount + 1
         case(5,25)
            internalFluxBoundariesWithPipes(ifwpCount)%attributeIDs(1) = 
     &          XdmfAddAttribute(descript % xdmfMD % xdmfFortranObj, 
     &          'IBCONN'//CHAR(0),
     &          XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), 
     &          XDMF_ARRAY_TYPE_INT32, 
     &          internalFluxBoundariesWithPipes(ifwpCount)%xdmf_ibconn(:)) ! zero offset
            internalFluxBoundariesWithPipes(ifwpCount)%attributeIDs(2) = 
     &          XdmfAddAttribute(descript % xdmfMD % xdmfFortranObj, 
     &          'BARINHT'//CHAR(0), 
     &          XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), 
     &          XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundariesWithPipes(ifwpCount)%barinht)
            internalFluxBoundariesWithPipes(ifwpCount)%attributeIDs(3) = 
     &          XdmfAddAttribute(descript % xdmfMD % xdmfFortranObj, 
     &          'BARINCFSB'//CHAR(0), 
     &          XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), 
     &          XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundariesWithPipes(ifwpCount)%barincfsb)
            internalFluxBoundariesWithPipes(ifwpCount)%attributeIDs(4) = 
     &          XdmfAddAttribute(descript % xdmfMD % xdmfFortranObj, 
     &          'BARINCFSP'//CHAR(0), 
     &          XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), 
     &          XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundariesWithPipes(ifwpCount)%barincfsp)   
            internalFluxBoundariesWithPipes(ifwpCount)%attributeIDs(5) = 
     &          XdmfAddAttribute(descript % xdmfMD % xdmfFortranObj, 
     &          'PIPEHT'//CHAR(0), 
     &          XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), 
     &          XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundariesWithPipes(ifwpCount)%pipeht)
            internalFluxBoundariesWithPipes(ifwpCount)%attributeIDs(6) = 
     &          XdmfAddAttribute(descript % xdmfMD % xdmfFortranObj, 
     &          'PIPECOEF'//CHAR(0), 
     &          XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), 
     &          XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundariesWithPipes(ifwpCount)%pipecoef)
            internalFluxBoundariesWithPipes(ifwpCount)%attributeIDs(7) = 
     &          XdmfAddAttribute(descript % xdmfMD % xdmfFortranObj, 
     &          'PIPEDIAM'//CHAR(0), 
     &          XDMF_ATTRIBUTE_CENTER_NODE, XDMF_ATTRIBUTE_TYPE_SCALAR, nvell(i), 
     &          XDMF_ARRAY_TYPE_FLOAT64, internalFluxBoundariesWithPipes(ifwpCount)%pipediam)
            internalFluxBoundariesWithPipes(ifwpCount)%informationID = 
     &          XdmfAddInformation(descript % xdmfMD % xdmfFortranObj, 
     &          'IBTYPE'//CHAR(0), 
     &          trim(fluxBoundaryType)//CHAR(0))
            internalFluxBoundariesWithPipes(ifwpCount)%setID = 
     &          XdmfAddSet(descript % xdmfMD % xdmfFortranObj,
     &          'flux_specified_boundary'//char(0), XDMF_SET_TYPE_NODE,
     &          internalFluxBoundariesWithPipes(ifwpCount)%xdmf_nodes, 
     &          nvell(i), XDMF_ARRAY_TYPE_INT32)
            ifwpCount = ifwpCount + 1
         case default
            write(6,
     &       '("ERROR: File contains IBTYPE=",i0," which is invalid.")'), ibtype_orig(i)
         end select
      end do
#endif
      !----------------------------------------------------------------------
      end subroutine addBoundaries
      !----------------------------------------------------------------------


      !----------------------------------------------------------------------
      ! S U B R O U T I N E   A D D   B O U N D A R Y  R E F E R E N C E S
      !----------------------------------------------------------------------
      ! Add boundary references to each grid in the collection.
      !----------------------------------------------------------------------
      subroutine addBoundaryReferences(descript)
      use boundaries
      use global, only : OutputDataDescript_t
      implicit none
#ifdef ADCXDMF      
      include 'adcirc_Xdmf.f' 
      !
      type(OutputDataDescript_t) :: descript
      character(len=256) :: fluxBoundaryType
      integer :: i, j  ! loop counter
      !
      do i=1, nope 
         elevationBoundaries(i)%informationID = 
     &        XdmfAddInformation(descript % xdmfMD % xdmfFortranObj, 
     &       'IBTYPEE'//CHAR(0), '0'//CHAR(0))
         elevationBoundaries(i)%setID = XdmfAddSet(
     &       descript % xdmfMD % xdmfFortranObj, 
     &       'elevation_specified_boundary'//char(0), XDMF_SET_TYPE_NODE, 
     &       elevationBoundaries(i)%xdmf_nodes, nvdll(i), XDMF_ARRAY_TYPE_INT32)
      end do
      sfCount = 1
      efCount = 1
      ifCount = 1
      ifwpCount = 1 
      do i=1, nbou 
         write(fluxBoundaryType,'(i0)') ibtype_orig(i)
         select case(ibtype_orig(i))
         case(0,1,2,10,11,12,20,21,22,30,52)  
            simpleFluxBoundaries(sfCount)%informationID = 
     &          XdmfAddInformation(descript % xdmfMD % xdmfFortranObj,
     &          'IBTYPE'//CHAR(0), 
     &          trim(fluxBoundaryType)//CHAR(0))
            simpleFluxBoundaries(sfCount)%setID = XdmfAddSet(
     &          descript % xdmfMD % xdmfFortranObj, 
     &          'flux_specified_boundary'//char(0), XDMF_SET_TYPE_NODE, 
     &          simpleFluxBoundaries(sfCount)%xdmf_nodes(:), nvell(i), XDMF_ARRAY_TYPE_INT32)
            sfCount = sfCount + 1
         case(3,13,23)
            do j=1, externalFluxBoundaries(efCount)%numAttributes
               call xdmfAddPreviousAttribute(descript % xdmfMD % xdmfFortranObj, 
     &            externalFluxBoundaries(efCount)%attributeIDs(j))
            end do
            externalFluxBoundaries(efCount)%informationID = 
     &          XdmfAddInformation(descript % xdmfMD % xdmfFortranObj, 
     &          'IBTYPE'//CHAR(0), 
     &          trim(fluxBoundaryType)//CHAR(0))
            externalFluxBoundaries(efCount)%setID = 
     &          XdmfAddSet(descript % xdmfMD % xdmfFortranObj, 
     &          'flux_specified_boundary'//char(0), XDMF_SET_TYPE_NODE, 
     &          externalFluxBoundaries(efCount)%xdmf_nodes(:), nvell(i), XDMF_ARRAY_TYPE_INT32)
            efCount = efCount + 1
         case(4,24)
            do j=1, internalFluxBoundaries(ifCount)%numAttributes
               call xdmfAddPreviousAttribute(descript % xdmfMD % xdmfFortranObj, 
     &                internalFluxBoundaries(ifCount)%attributeIDs(j))
            end do
            internalFluxBoundaries(ifCount)%informationID = 
     &         XdmfAddInformation(descript % xdmfMD % xdmfFortranObj, 
     &         'IBTYPE'//CHAR(0), trim(fluxBoundaryType)//CHAR(0))
            internalFluxBoundaries(ifCount)%setID = 
     &          XdmfAddSet(descript % xdmfMD % xdmfFortranObj, 
     &          'flux_specified_boundary'//char(0), XDMF_SET_TYPE_NODE, 
     &          internalFluxBoundaries(ifCount)%xdmf_nodes, nvell(i), XDMF_ARRAY_TYPE_INT32)
            ifCount = ifCount + 1
         case(5,25)
            do j=1, internalFluxBoundariesWithPipes(ifwpCount)%numAttributes
               call xdmfAddPreviousAttribute(descript % xdmfMD % xdmfFortranObj, 
     &                internalFluxBoundariesWithPipes(ifwpCount)%attributeIDs(j))
            end do
            internalFluxBoundariesWithPipes(ifwpCount)%informationID = 
     &          XdmfAddInformation(descript % xdmfMD % xdmfFortranObj, 
     &          'IBTYPE'//CHAR(0), 
     &          trim(fluxBoundaryType)//CHAR(0))
            internalFluxBoundariesWithPipes(ifwpCount)%setID = 
     &          XdmfAddSet(descript % xdmfMD % xdmfFortranObj,
     &          'flux_specified_boundary'//char(0), XDMF_SET_TYPE_NODE,
     &          internalFluxBoundariesWithPipes(ifwpCount)%xdmf_nodes, 
     &          nvell(i), XDMF_ARRAY_TYPE_INT32)
            ifwpCount = ifwpCount + 1
         case default
            write(6,'("ERROR: File contains IBTYPE=",i0," which is not a valid flux boundary type.")'), ibtype_orig(i)
         end select
      end do
#endif
      !----------------------------------------------------------------------
      end subroutine addBoundaryReferences
      !----------------------------------------------------------------------

      !---------------------------------------------------------------------
      ! S U B R O U T I N E    W R I T E    C O N T R O L   X D M F
      !---------------------------------------------------------------------
      subroutine writeControlXDMF(xdmfFortranObj)
      use control
      use boundaries, only : neta, nvel, nfluxf, lbcodei
      use mesh, only : slam0, sfea0
      implicit none
      include 'adcirc_Xdmf.f'
      integer*8, intent(in) :: xdmfFortranObj ! XDMF object to receive data
      character(len=80) :: wtimincLine
      character(len=80) :: rstimincLine
      character(len=80) :: cicetimincLine
      character(len=80) :: info
      character(len=80) :: varname
      character(len=80) :: stationType
      character(len=80) :: wtimincComment
      character(len=80) :: drampComment
      integer :: informationID
      integer :: controlnws ! nws as encoded with ncice and nrs parameters
      integer :: i 
      integer :: j
      !
      write(6,'(a)') 'INFO: Writing control parameters to XDMF.'
      informationID = XdmfAddInformation(xdmfFortranObj, 'rundes'//char(0), 
     &    trim(adjustl(rundes))//char(0))
      informationID = XdmfAddInformation(xdmfFortranObj, 'runid'//char(0), 
     &    trim(adjustl(runid))//char(0))
      write(info,fmt='(i0)') nfover
      informationID = XdmfAddInformation(xdmfFortranObj, 'nfover'//char(0), 
     &    trim(adjustl(info))//char(0))
      
      write(info,fmt='(i0)') nabout
      informationID = XdmfAddInformation(xdmfFortranObj, 'nabout'//char(0), 
     &    trim(adjustl(info))//char(0))
      write(info,fmt='(i0)') nscreen
      informationID = XdmfAddInformation(xdmfFortranObj, 'nscreen'//char(0), 
     &    trim(adjustl(info))//char(0))
      write(info,fmt='(i0)') ihot
      informationID = XdmfAddInformation(xdmfFortranObj, 'ihot'//char(0), 
     &    trim(adjustl(info))//char(0))
      write(info,fmt='(i0)') ics
      informationID = XdmfAddInformation(xdmfFortranObj, 'ics'//char(0), 
     &    trim(adjustl(info))//char(0))
      write(info,fmt='(i0)') im
      informationID = XdmfAddInformation(xdmfFortranObj, 'im'//char(0), 
     &    trim(adjustl(info))//char(0))
      if ((im.eq.21).or.(im.eq.31)) then
         write(info,fmt='(i0)') iden
         informationID = XdmfAddInformation(xdmfFortranObj, 'iden'//char(0), 
     &       trim(adjustl(info))//char(0))
      endif
      write(info,fmt='(i0)') nolibf
      informationID = XdmfAddInformation(xdmfFortranObj, 'nolibf'//char(0), 
     &    trim(adjustl(info))//char(0))
      write(info,fmt='(i0)') nolifa
      informationID = XdmfAddInformation(xdmfFortranObj, 'nolifa'//char(0), 
     &    trim(adjustl(info))//char(0))
      write(info,fmt='(i0)') nolica
      informationID = XdmfAddInformation(xdmfFortranObj, 'nolica'//char(0), 
     &    trim(adjustl(info))//char(0))
      write(info,fmt='(i0)') nolicat
      informationID = XdmfAddInformation(xdmfFortranObj, 'nolicat'//char(0), 
     &    trim(adjustl(info))//char(0))
      write(info,fmt='(i0)') nwp
      informationID = XdmfAddInformation(xdmfFortranObj, 'nwp'//char(0), 
     &    trim(adjustl(info))//char(0))
      
      if (nwp.ne.0) then
         do i=1, nwp
            write(info,'("attribute_name(",i0,")")') i
            informationID = XdmfAddInformation(xdmfFortranObj, 
     &          trim(adjustl(info))//char(0), trim(adjustl(attrName(i)))//char(0))
         end do
      endif
      write(info,fmt='(i0)') ncor
      informationID = XdmfAddInformation(xdmfFortranObj, 'ncor'//char(0), 
     &    trim(adjustl(info))//char(0))
      write(info,fmt='(i0)') ntip
      informationID = XdmfAddInformation(xdmfFortranObj, 'ntip'//char(0), 
     &    trim(adjustl(info))//char(0))
      controlnws = 1000 * ncice + 100 * nrs + abs(nws)
      if (nws.ne.0) then
         controlnws = controlnws * abs(nws)/nws ! to get the sign right
      endif
      write(info,fmt='(i0)') controlnws
      informationID = XdmfAddInformation(xdmfFortranObj, 'nws'//char(0), 
     &    trim(adjustl(info))//char(0))
      write(info,fmt='(i0)') nramp
      informationID = XdmfAddInformation(xdmfFortranObj, 'nramp'//char(0), 
     &    trim(adjustl(info))//char(0))
      write(info,fmt='(f10.5)') g
      informationID = XdmfAddInformation(xdmfFortranObj, 'g'//char(0), 
     &    trim(adjustl(info))//char(0))
      write(info,fmt='(f6.3)') tau0
      informationID = XdmfAddInformation(xdmfFortranObj, 'tau0'//char(0), 
     &    trim(adjustl(info))//char(0))
      if ((tau0.le.-5.0).and.(tau0.ge.-5.99)) then
         write(info,fmt='(f6.3,f6.3)') tau0fulldomainmin,tau0fulldomainmax
         informationID = XdmfAddInformation(xdmfFortranObj, 
     &       'tau0fulldomainmin,tau0fulldomainmax'//char(0),  
     &       trim(adjustl(info))//char(0))
      endif
      write(info,fmt='(f15.7)') dtdp
      informationID = XdmfAddInformation(xdmfFortranObj,'dtdp'//char(0), 
     &    trim(adjustl(info))//char(0))
      write(info,fmt='(f15.7)') statim
      informationID = XdmfAddInformation(xdmfFortranObj,'statim'//char(0), 
     &    trim(adjustl(info))//char(0))
      write(info,fmt='(f15.7)') reftim
      informationID = XdmfAddInformation(xdmfFortranObj,'reftim'//char(0), 
     &    trim(adjustl(info))//char(0))
      
      wtimincLine = ''
      wtimincComment = ''
      select case(abs(nws))
      case(0,1) 
         ! do nothing
      case(2,4,5,7,10,11,12,15,16)   
         wtimincComment = 'wtiminc'
      case(3)   
         wtimincComment = 'irefyr,irefmo,irefday,irefhr,irefmin,refsec'
         write(wtimincLine,*) irefyr,irefmo,irefday,irefhr,irefmin,refsec 
         informationID = XdmfAddInformation(xdmfFortranObj,trim(wtimincComment)//char(0),  
     &       trim(adjustl(wtimincLine))//char(0))
         wtimincComment = 'nwlat,nwlon,wlatmax,wlonmin,wlatinc,wloninc,wtiminc'
         write(wtimincLine,*) nwlat,nwlon,wlatmax,wlonmin,wlatinc,wloninc,wtiminc
      case(6)
         wtimincComment = 'nwlat,nwlon,wlatmax,wlonmin,wlatinc,wloninc,wtiminc'
         write(wtimincLine,*) nwlat,nwlon,wlatmax,wlonmin,wlatinc,wloninc,wtiminc
      case(8,19)
         wtimincComment = 'irefyr,irefmo,irefday,irefhr,stormnumber,bladj'
         write(wtimincLine,*) irefyr,irefmo,irefday,irefhr,stormnumber,bladj
      case(29)
         wtimincComment = 'irefyr,irefmo,irefday,irefhr,stormnumber,'//
     &    'bladj,wtiminc,purevortex,purebackground'
         write(wtimincLine,*) irefyr,irefmo,irefday,irefhr,stormnumber,
     &    bladj,wtiminc,purevortex, purebackground
      case default
         write(info,'(a,i0,a)') 'ERROR: The value of nws is ',nws,
     &     ' but this is not a valid value.'
         call exit(1)
      end select
      if (nrs.gt.0) then
          wtimincComment = trim(wtimincComment) // ',rstiminc'
          write(rstimincLine,*) rstiminc
          wtimincLine = trim(wtimincLine) // ' ' // trim(rstimincLine)
      endif
      if (ncice.gt.0) then
          wtimincComment = trim(wtimincComment) // ',cice_timinc'
          write(cicetimincLine,*) cice_timinc
          wtimincLine = trim(wtimincLine) // ' ' // trim(cicetimincLine)
      endif
      if ((nws.ne.0).or.(nrs.ne.0).or.(ncice.ne.0)) then
         informationID = XdmfAddInformation(xdmfFortranObj,trim(wtimincComment)//char(0), 
     &       trim(adjustl(wtimincLine))//char(0))
      endif
      
      
      write(info,fmt='(e20.10)') rnday
         informationID = XdmfAddInformation(xdmfFortranObj,'rnday'//char(0), 
     &       trim(adjustl(info))//char(0))
      !
      select case(nramp)
      case(0,1)
         write(info,fmt='(9f6.3)') dramp  
         drampComment = 'dramp'
      case(2) 
         write(info,fmt='(9f6.3)') 
     &       dramp,drampextflux,fluxsettlingtime
         drampComment = 'dramp,drampextflux,fluxsettlingtime'
      case(3) 
         write(info,fmt='(9f6.3)') 
     &       dramp,drampextflux,fluxsettlingtime,drampintflux
       drampComment = 'dramp,drampextflux,fluxsettlingtime,drampintflux'
      case(4) 
         write(info,fmt='(9f6.3)') 
     &       dramp,drampextflux,fluxsettlingtime,drampintflux, drampelev
         drampComment = 'dramp,drampextflux,fluxsettlingtime,'//
     &    'drampintflux, drampelev'
      case(5) 
         write(info,fmt='(9f6.3)') 
     &       dramp,drampextflux,fluxsettlingtime,drampintflux,drampelev,dramptip
         drampComment = 'dramp,drampextflux,fluxsettlingtime,drampintflux,drampelev,dramptip'
      case(6) 
         write(info,fmt='(9f6.3)') 
     &        dramp,drampextflux,fluxsettlingtime,drampintflux,drampelev,dramptip,drampmete
         drampComment = 'dramp,drampextflux,fluxsettlingtime,drampintflux,drampelev,dramptip,drampmete'
      case(7) 
         write(info,fmt='(9f6.3)') 
     &       dramp,drampextflux,fluxsettlingtime,drampintflux, 
     &              drampelev,dramptip,drampmete,drampwrad
       drampComment = 'dramp,drampextflux,fluxsettlingtime,drampintflux,' // 
     &              'drampelev,dramptip,drampmete,drampwrad'
      case(8) 
         write(info,fmt='(9f6.3)') 
     &          dramp,drampextflux,fluxsettlingtime,drampintflux, 
     &              drampelev,dramptip,drampmete,drampwrad,dunrampmete
      drampcomment = 'dramp,drampextflux,fluxsettlingtime,drampintflux,' // 
     &              'drampelev,dramptip,drampmete,drampwrad,dunrampmete'
      case default
         write(info,'("ERROR: The value of nramp is ",i2," but this is not a valid value.")') nramp
         call exit(1)
      end select
      informationID = XdmfAddInformation(xdmfFortranObj,trim(drampComment)//char(0),  
     &    trim(adjustl(info))//char(0))
      !
      write(info,fmt='(3f6.3)') a00,b00,c00
      informationID = XdmfAddInformation(xdmfFortranObj,'a00,b00,c00'//char(0), 
     &    trim(adjustl(info))//char(0))
      if (nolifa.eq.2) then
         write(info,fmt='(f15.7, i2, i2, f6.3)') 
     &       h0,nodedrymin,nodewetmin,velmin
         informationID = XdmfAddInformation(xdmfFortranObj, 
     &       'h0,nodedrymin,nodewetmin,velmin'//char(0), 
     &       trim(adjustl(info))//char(0))
      else
         write(info,fmt='(f15.7)') h0
         informationID = XdmfAddInformation(xdmfFortranObj,'h0'//char(0),  
     &       trim(adjustl(info))//char(0))
      endif
      write(info,fmt='(2f15.7)') slam0,sfea0
         informationID = XdmfAddInformation(xdmfFortranObj,'slam0,sfea0'//char(0), 
     &       trim(adjustl(info))//char(0))
      select case(nolibf)
      case(0)
         write(info,fmt='(f15.7)') tau
         informationID = XdmfAddInformation(xdmfFortranObj,'tau'//char(0), 
     &       trim(adjustl(info))//char(0))
      case(1)
         write(info,fmt='(f15.7)') cf
         informationID = XdmfAddInformation(xdmfFortranObj,'cf'//char(0),  
     &       trim(adjustl(info))//char(0))
      case(2)
         write(info,fmt='(4f15.7)') cf,hbreak,ftheta,fgamma
         informationID = XdmfAddInformation(xdmfFortranObj,
     &      'cf,hbreak,ftheta,fgamma'//char(0), 
     &       trim(adjustl(info))//char(0))
      case default
         write(info,'(a,i0,a)') 'ERROR: The value of nolibf is ',nolibf,
     &     ' but this is not a valid value.'
         call exit(1)
      end select
      if (im.eq.10) then
         write(info,fmt='(2f15.7)') eslm,eslc
         informationID = XdmfAddInformation(xdmfFortranObj,'eslm,eslc'//char(0), 
     &       trim(adjustl(info))//char(0))
      else
        write(info,fmt='(f15.7)') eslm
        informationID = XdmfAddInformation(xdmfFortranObj,'eslm'//char(0), 
     &       trim(adjustl(info))//char(0))
      endif
      write(info,fmt='(f15.7)') cori
      informationID = XdmfAddInformation(xdmfFortranObj,'cori'//char(0), 
     &    trim(adjustl(info))//char(0))
      write(info,fmt='(i0)') ntif
      informationID = XdmfAddInformation(xdmfFortranObj,'ntif'//char(0),  
     &    trim(adjustl(info))//char(0))
      do i=1,ntif
         write(info,'(a)')  trim(tipotag(i))
         write(varname,'("tipotag(",i0,")")') i
         informationID = XdmfAddInformation(xdmfFortranObj, 
     &       trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))     
         write(info,fmt='(5f15.7)') tpk(i),amigt(i),etrf(i),fft(i),facet(i)
         write(varname,'(5(a,i0),a)') "tpk(",i,"),amigt(",i,
     &      "),etrf(",i,"),fft(",i,"),facet(",i,")" 
         informationID = XdmfAddInformation(xdmfFortranObj,  
     &       trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))
      end do
      
      write(info,fmt='(i0)') nbfr
      informationID = XdmfAddInformation(xdmfFortranObj,'nbfr'//char(0),  
     &    trim(adjustl(info))//char(0))
      do i=1, nbfr
         write(info,'(a)') trim(bountag(i))
         write(varname,'("boundtag(",i0,")")') i
         informationID = XdmfAddInformation(xdmfFortranObj, 
     &       trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))   
         write(info,fmt='(3f15.7)') amig(i),ff(i),face(i)
         write(varname,'("amig(",i0,"),ff(",i0,"),face(",i0,")")' ) 
     &        i, i, i
         informationID = XdmfAddInformation(xdmfFortranObj,  
     &       trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))   
      end do
      do i=1, nbfr
         write(info,'(a)') trim(nbfr_alpha(i))
         write(varname,'("nbfr_alpha(",i0,")")') i
         informationID = XdmfAddInformation(xdmfFortranObj, 
     &       trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))   
         do j=1,neta
            write(varname,'("emo(",i0,",",i0,"),efa(",i0,",",i0,")")' ) 
     &        i, j, i, j
            write(info,fmt='(2f15.7)') emo(i,j),efa(i,j)
               informationID = XdmfAddInformation(xdmfFortranObj,  
     &             trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))
         end do
      end do
      write(info,fmt='(f15.7)') anginn
      informationID = XdmfAddInformation(xdmfFortranObj, 
     &    'anginn'//char(0), trim(adjustl(info))//char(0))
      if ((nfluxf.ne.0).and.(nffr.ne.0)) then
         write(info,fmt='(i0)') nffr
         informationID = XdmfAddInformation(xdmfFortranObj,  
     &       'nffr'//char(0), trim(adjustl(info))//char(0))
         if (nffr.gt.0) then
            do i=1,nffr
               write(info,'(a)') trim(fbountag(i))
               write(varname,'("fboundtag(",i0,")")') i
               informationID = XdmfAddInformation(xdmfFortranObj,  
     &             trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))
               write(info,fmt='(3f6.3)') 
     &             famig(i),fff(i),fface(i)
             write(varname,'("famig(",i0,"),fff(",i0,"),fface(",i0,")")' ) 
     &             i, i, i
               informationID = XdmfAddInformation(xdmfFortranObj, 
     &             trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))
            end do 
            do i=1,nffr
               write(info,'(a)') trim(nffr_alpha(i))
               write(varname,'("nffr_alpha(",i0,")")') i
               informationID = XdmfAddInformation(xdmfFortranObj,  
     &             trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))
               do j=1,nvel
                  select case(lbcodei(j))
                  case(2,12,22,52)
                     write(info,fmt='(2f6.3)') qnam(i,j), qnph(i,j)
           write(varname,'("qnam(",i0,",",i0,"),qnph(",i0,",",i0,")")' ) 
     &                   i, j, i, j
                     informationID = XdmfAddInformation(xdmfFortranObj,  
     &                   trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))
                  case default
                     cycle
                  end select
               end do
            end do
         endif
      endif
      
      write(info,fmt='(i0,1x,2f6.3,1x,i0)') noute,toutse,toutfe,nspoole
      informationID = XdmfAddInformation(xdmfFortranObj,
     &   'noute,toutse,toutfe,nspoole'//char(0), 
     &    trim(adjustl(info))//char(0))
      write(info,fmt='(i0)') nstae
      informationID = XdmfAddInformation(xdmfFortranObj,'nstae'//char(0), 
     &    trim(adjustl(info))//char(0))
      stationType = 'sea_surface_height_above_geoid'
      
      call writeStationsXDMF(nstae, xel, yel, stationType, xdmfFortranObj)
      
      write(info,fmt='(i0,1x,2f15.7,1x,i0)') noutv,toutsv,toutfv,nspoolv
      informationID = XdmfAddInformation(xdmfFortranObj,
     &   'noutv,toutsv,toutfv,nspoolv'//char(0), 
     &    trim(adjustl(info))//char(0))
      write(info,fmt='(i0)') nstav
      informationID = XdmfAddInformation(xdmfFortranObj,'nstav'//char(0), 
     &    trim(adjustl(info))//char(0))
      stationType = 'water_column_vertically_averaged_velocity'
      call writeStationsXDMF(nstav, xev, yev, stationType, xdmfFortranObj)
      if (im.eq.10) then
         write(info,fmt='(i0,1x,2f15.7,1x,i0)') noutc,toutsc,toutfc,nspoolc
         informationID = XdmfAddInformation(xdmfFortranObj,
     &    'noutc,toutsc,toutfc,nspoolc'//char(0), trim(adjustl(info))//char(0))
         write(info,fmt='(i0)') nstac
         informationID = XdmfAddInformation(xdmfFortranObj,'nstac'//char(0), 
     &       trim(adjustl(info))//char(0))
         stationType = 'water_column_vertically_averaged_concentration'   
         call writeStationsXDMF(nstav, xev, yev, stationType, xdmfFortranObj)
      endif
      if (nws.ne.0) then
         write(info,fmt='(i0,1x,2f15.7,1x,i0)') noutm,toutsm,toutfm,nspoolm
         informationID = XdmfAddInformation(xdmfFortranObj,
     &      'noutm,toutsm,toutfm,nspoolm'//char(0), trim(adjustl(info))//char(0))   
         write(info,fmt='(i0)') nstam
         informationID = XdmfAddInformation(xdmfFortranObj,'nstam'//char(0),  
     &       trim(adjustl(info))//char(0))   
         stationType = 
     &      'atmospheric_pressure_and_wind_velocity_at_sea_level'   
         call writeStationsXDMF(nstam, xem, yem, stationType, xdmfFortranObj)
      endif
      write(info,fmt='(i0,1x,2f15.7,1x,i0)') noutge,toutsge,toutfge,nspoolge
      informationID = XdmfAddInformation(xdmfFortranObj,
     &  'noutge,toutsge,toutfge,nspoolge'//char(0), trim(adjustl(info))//char(0))
      write(info,fmt='(i0,1x,2f15.7,1x,i0)') noutgv,toutsgv,toutfgv,nspoolgv
      informationID = XdmfAddInformation(xdmfFortranObj,
     &   'noutgv,toutsgv,toutfgv,nspoolgv'//char(0), trim(adjustl(info))//char(0))
      if (im.eq.10) then
         write(info,fmt='(i0,1x,2f15.7,1x,i0)') noutgc,toutsgc,toutfgc,nspoolgc
         informationID = XdmfAddInformation(xdmfFortranObj,
     &    'noutgc,toutsgc,toutfgc,nspoolgc'//char(0), trim(adjustl(info))//char(0))
      endif
      if (nws.ne.0) then
         write(info,fmt='(i0,1x,2f15.7,1x,i0)') noutgw,toutsgw,toutfgw,nspoolgw
         informationID = XdmfAddInformation(xdmfFortranObj,
     &    'noutgw,toutsgw,toutfgw,nspoolgw'//char(0), trim(adjustl(info))//char(0))   
      endif
      
      write(info,fmt='(i0)') nfreq
      informationID = XdmfAddInformation(xdmfFortranObj,'nfreq'//char(0), 
     &    trim(adjustl(info))//char(0))
      if (nfreq.ne.0) then
         do i=1,nfreq
            write(info,'(a10)') namefr(i)
            write(varname,'("namefr(",i0,")")') i
               informationID = XdmfAddInformation(xdmfFortranObj,  
     &             trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))      
            write(info,fmt='(3f15.7)') hafreq(i),haff(i),haface(i)
            write(varname,
     &         '("hafreq(",i0,"),haff(",i0,"),haface(",i0,")")' ) 
     &          i, i, i
            informationID = XdmfAddInformation(xdmfFortranObj, 
     &          trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))
         end do
      endif
      write(info,fmt='(2f15.7,1x,i0,1x,f15.7)') thas,thaf,nhainc,fmv
      informationID = XdmfAddInformation(xdmfFortranObj,
     &    'thas,thaf,nhainc,fmv' 
     &    //char(0), trim(adjustl(info))//char(0))
      write(info,fmt='(i0,1x,i0,1x,i0,1x,i0)') nhase,nhasv,nhage,nhagv
      informationID = XdmfAddInformation(xdmfFortranObj,
     &    'nhase,nhasv,nhage,nhagv' 
     &    //char(0), trim(adjustl(info))//char(0))
      write(info,fmt='(i0,1x,i0)') nhstar,nhsinc
      informationID = XdmfAddInformation(xdmfFortranObj,
     &    'nhstar,nhsinc' 
     &    //char(0), trim(adjustl(info))//char(0))
      write(info,fmt='(i0,1x,i0,f15.7,1x,i0)') ititer,isldia,convcr,itmax
      informationID = XdmfAddInformation(xdmfFortranObj,
     &    'ititer,isldia,convcr,itmax' 
     &    //char(0), trim(adjustl(info))//char(0))
      !
      ! read metadata if netcdf or xdmf format was specified in the output 
      ! identifiers
      if (readMetaData.eqv..true.) then 
         informationID = XdmfAddInformation(xdmfFortranObj,'title'//char(0), 
     &       trim(adjustl(title))//char(0))   
         informationID = XdmfAddInformation(xdmfFortranObj,'institution'//char(0), 
     &       trim(adjustl(institution))//char(0))
         informationID = XdmfAddInformation(xdmfFortranObj,'source'//char(0), 
     &       trim(adjustl(source))//char(0))   
         informationID = XdmfAddInformation(xdmfFortranObj,'history'//char(0), 
     &       trim(adjustl(history))//char(0))   
         informationID = XdmfAddInformation(xdmfFortranObj,'references'//char(0), 
     &       trim(adjustl(references))//char(0))   
         informationID = XdmfAddInformation(xdmfFortranObj,'comments'//char(0), 
     &       trim(adjustl(comments))//char(0))   
         informationID = XdmfAddInformation(xdmfFortranObj,'host'//char(0), 
     &       trim(adjustl(host))//char(0))   
         informationID = XdmfAddInformation(xdmfFortranObj,'convention'//char(0), 
     &       trim(adjustl(convention))//char(0))   
         informationID = XdmfAddInformation(xdmfFortranObj,'contact'//char(0), 
     &       trim(adjustl(contact))//char(0))   
         informationID = XdmfAddInformation(xdmfFortranObj,'base_date'//char(0), 
     &       trim(adjustl(base_date))//char(0))
      endif 
      !---------------------------------------------------------------------
      end subroutine writeControlXDMF
      !---------------------------------------------------------------------
      
      !---------------------------------------------------------------------
      ! S U B R O U T I N E    W R I T E   S T A T I O N S   X D M F
      !---------------------------------------------------------------------
      subroutine writeStationsXDMF(nsta, stax, stay, stationType, xdmfFortranObj)
      implicit none
      include 'adcirc_Xdmf.f'
      integer, intent(in) :: nsta
      real(sz), intent(in) :: stax(:)
      real(sz), intent(in) :: stay(:)
      character(len=80), intent(in) :: stationType
      integer*8, intent(in) :: xdmfFortranObj ! XDMF object to receive data
      integer :: informationID ! information index, used as ID to insert array
      !
      ! add the station locations as an array related to the station information object
      informationID = XdmfAddInformation(xdmfFortranObj,
     &    'station longitude'//char(0),  
     &    trim(adjustl(stationType))//char(0))
      ! arguments are: fortran obj, information index, values, numValues, arrayType
      call XdmfAddInformationArray(xdmfFortranObj, informationID, stax, nsta, 
     &    XDMF_ARRAY_TYPE_FLOAT64)
      informationID = XdmfAddInformation(xdmfFortranObj,
     &    'station latitude'//char(0),  
     &    trim(adjustl(stationType))//char(0))
      call XdmfAddInformationArray(xdmfFortranObj, informationID, stay, nsta, 
     &    XDMF_ARRAY_TYPE_FLOAT64)
      !---------------------------------------------------------------------
      end subroutine writeStationsXDMF
      !---------------------------------------------------------------------

      !-----------------------------------------------------------------
      end module xdmfio
      !-----------------------------------------------------------------
#endif
