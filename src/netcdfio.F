#ifdef ADCNETCDF
      !=================================================================
      !=================================================================
      !=================================================================
      !      =====                                           =====
      !      =====            MODULE NetCDFIO                =====
      !      =====                                           =====
      !=================================================================
      !=================================================================
      !=================================================================

      !=================================================================
      ! This module provides a NetCDF I/O capability for ADCIRC.
      !
      ! Revision history:
      !
      ! Date      Programmer                       Description of change
      ! ----      ----------                       ---------------------
      ! 03/30/07  Cristina Forbes, PSGS @ UNC-IMS  Wrote original code
      ! 03/30/08  Cristina Forbes, PSGS @ UNC-IMS  Modified code for
      !                                            globalio & hotstart
      !                                            from binary files
      ! 09/30/08 Cristina Forbes @ UNC-IMS         Modified metadata
      ! 10/15/08 Cristina Forbes @ UNC-IMS         Continued modifying metadata
      ! 5/21/08  Cristina Forbes @ UNC-IMS         Fixed hotstart write67
      !                                            seg-fault &  define portion
      !                                            for grids with no specified
      !                                            boundary forcing segments or nodes
      ! 10/20/09 Chris Massey @ USACE-ERDC-CHL     changed reserved word "count" to
      ! v49.01                                     kount to avoid conflicts.
      ! 07-08/10 Jason Fleming                     complete reorganization
      !                                            for greater modularity,
      !                                            flexibility, extensibility,
      !                                            and maintainability
      ! v51.20.06 Chris Massey @ USACE-ERDC-CHL    added time information
      !                                            for when max/min occurs
      ! 02.20.2018 William Pringle at UND          Changed NBVV and NBDV to
      !                                            vectors to avoid large netcdf
      !                                            files. Added in fort.51-54.nc
      !=================================================================
      MODULE NETCDFIO
      USE SIZES, ONLY : SZ, OFF, ASCII, BINARY, SPARSE_ASCII, NETCDF3,
     &   NETCDF4, XDMF
      USE GLOBAL, ONLY : DEBUG, ECHO, INFO, WARNING, ERROR,
     &   screenMessage, logMessage, allMessage, setMessageSource,
     &   unsetMessageSource, scratchMessage
     &
      USE NETCDF        ! enables the fortran 90 interface to netcdf
      IMPLICIT NONE
      include 'netcdf.inc' ! enables the fortran 77 interface to netcdf
C
      REAL(8), PARAMETER ::  doubleval=-99999.d0
      INTEGER, PARAMETER :: intFillValue=-99999
C
      TYPE, PRIVATE :: meshStructure
         LOGICAL :: initialized = .false.
         INTEGER :: X_id          ! x-coordinate or longitude
         INTEGER :: Y_id          ! y-coordinate or latitude
         INTEGER :: sigma_id      ! relative depth of 3D layer
         INTEGER :: DEPTH_id      ! distance from geoid
         INTEGER :: ELE_id        ! elements in grid
         INTEGER :: nbdvnc_id       ! nodes on elev spec boundary seg
         INTEGER :: nbvvnc_id       ! nodes on normal flow boundary seg
         INTEGER :: nvdllnc_id      ! num nodes on elev boundary seg
         INTEGER :: ibtypenc_id   ! discharge boundary type
         INTEGER :: ibtypeenc_id  ! elevation boundary type
         INTEGER :: nvellnc_id    ! nodes on norm flow spec boundary seg
         INTEGER :: max_nvdllnc_id ! max num nodes on any elev boundary seg
         INTEGER :: max_nvellnc_id ! max num nodes on any discharge boundary seg
         INTEGER :: slam0nc_id
         INTEGER :: sfea0nc_id
         INTEGER :: netanc_id
         INTEGER :: nvelnc_id
         INTEGER :: mesh_id        ! Mesh Definition (Corbitt)
!        Dimension ids
         INTEGER :: num_nodes_dim_id
         INTEGER :: num_elems_dim_id
         INTEGER :: nface_dim_id
         INTEGER :: nopenc_dim_id   ! num elev spec boundary forcing segs
         INTEGER :: nbounc_dim_id  ! number of normal flow specified boundary segment
         INTEGER :: ibtypenc_dim_id
         INTEGER :: ibtypeenc_dim_id  ! elevation boundary types
         INTEGER :: max_nvdllnc_dim_id
         INTEGER :: max_nvellnc_dim_id
         INTEGER :: num_v_nodes_dim_id
         INTEGER :: mesh_dim_id       ! Mesh Definition (Corbitt)
         !WJP 02.20.2018: for one dimensional boundary
         INTEGER :: nvelnc_dim_id
         INTEGER :: netanc_dim_id
!        Dimension lengths
         INTEGER :: num_nodes
         INTEGER :: num_elems
         INTEGER :: nface_len
         INTEGER :: nopenc
         INTEGER :: nbounc
         INTEGER :: max_nvdllnc
         INTEGER :: max_nvellnc
         INTEGER :: num_v_nodes
!        Rank (number of dimensions) for each variable
!        Variable shapes
         INTEGER :: x_dims(1)
         INTEGER :: y_dims(1)
         INTEGER :: sigma_dims(1)
         INTEGER :: depth_dims(1)
         INTEGER :: ele_dims(2)
         INTEGER :: nvdll_dims(1)
         INTEGER :: nbounc_dims(1)
         INTEGER :: ibtypenc_dims(1)
         INTEGER :: ibtypeenc_dims(1)
         INTEGER :: nvellnc_dims(1)
         INTEGER :: nvdllnc_dims(1)
         INTEGER :: nopenc_dims(1)
         INTEGER :: nbdvnc_dims(2)
         INTEGER :: nbvvnc_dims(2)
         INTEGER :: mesh_dims(1)       ! Mesh Definition (Corbitt)
         REAL(SZ), ALLOCATABLE :: xnc(:) ! x coordinate or longitude
         REAL(SZ), ALLOCATABLE :: ync(:) ! y coordinate or latitude
         !WJP 02.20.2018: changed to one dimensional boundary
         INTEGER, ALLOCATABLE ::  nbvvnc(:) ! boundary array
         INTEGER, ALLOCATABLE ::  nbdvnc(:) ! boundary array
         INTEGER, ALLOCATABLE ::  nvellnc(:)  ! boundary array
         INTEGER, ALLOCATABLE ::  nvdllnc(:)  ! boundary array
         INTEGER, ALLOCATABLE ::  ibtypenc(:)     ! boundary array
         INTEGER, ALLOCATABLE ::  ibtypeenc(:)
         INTEGER, ALLOCATABLE ::  nmnc(:,:)
         INTEGER, ALLOCATABLE ::  element(:,:)
         INTEGER :: netanc
         INTEGER :: nvelnc
         !WJP 02.20.2018: added fort.53-54 harmonic outputs
         INTEGER :: const_id
         INTEGER :: constl_dim_id
      END TYPE meshStructure

      type(meshStructure), private, save, target :: adcircMesh

      TYPE, PRIVATE :: fileData
         INTEGER :: record_counter
         INTEGER :: ncformat   ! netcdf file format to create
         LOGICAL :: createFile ! .true. if a new netCDF file must be created p
         CHARACTER(len=1024) :: filename
         LOGICAL :: fileFound  ! .true. if the netCDF file is present
      END Type fileData

      TYPE, PRIVATE :: timeData
         LOGICAL :: initialized = .false.
         INTEGER :: timenc_len = 1  ! number of time slices to write
         INTEGER :: timenc_dim_id
         INTEGER :: timenc_id
         INTEGER :: timenc_dims(1)
         REAL(8), ALLOCATABLE :: timenc(:)
      END TYPE timeData

      TYPE, PRIVATE :: stationData
         INTEGER :: ncid              ! the id of its netcdf file
         INTEGER :: num_stations      ! total number of stations
         INTEGER :: num_v_nodes
         INTEGER :: num_sta_dim_id
         INTEGER :: num_v_nodes_dim_id
         CHARACTER(50) :: varname(3)  ! name(s) of this data in netcdf file
         INTEGER, ALLOCATABLE :: name_lengths(:) ! lengths of station names
         CHARACTER(50), POINTER :: statnames(:)
         INTEGER :: slen_dim_id
         INTEGER :: station_id
         INTEGER :: station_dims(2)
         INTEGER :: x_id          ! station x-coordinate or longitude
         INTEGER :: y_id          ! station y-coordinate or latitude
         INTEGER :: x_dims(1)
         INTEGER :: y_dims(1)
         REAL(SZ), ALLOCATABLE :: x(:) ! x coordinate or longitude
         REAL(SZ), ALLOCATABLE :: y(:) ! y coordinate or latitude
         INTEGER :: station_data_id
         INTEGER :: u_station_data_id
         INTEGER :: v_station_data_id
         INTEGER :: w_station_data_id
         INTEGER :: station_data_dims(2)
         INTEGER :: station_data_dims_3D(3)
         !WJP 02.20.2018: added fort.51-52 harmonic outputs
         INTEGER :: const_id
         INTEGER :: constl_dim_id
         INTEGER :: station_ha_data_id
         INTEGER :: ha_u_station_data_id
         INTEGER :: ha_v_station_data_id
         INTEGER :: station_hg_data_id
         INTEGER :: hg_u_station_data_id
         INTEGER :: hg_v_station_data_id
         type(timeData) :: myTime
         type(fileData) :: myFile
      END TYPE stationData
C
      type(stationData), private, save :: elevSta ! elev stations (fort.61)
      type(stationData), private, save :: prSta   ! pressure stations (fort.71)
C---GML 20200712
      type(stationData), private, save :: IceSta   ! ice stations (fort.91)
C
      type(stationData), private, save :: velSta  ! velocity stations (fort.62)
      type(stationData), private, save :: wVelSta ! wind vel stations (fort.72)
      type(stationData), private, save :: densityStations3D    ! for fort.41
      type(stationData), private, save :: velocityStations3D   ! for fort.42
      type(stationData), private, save :: turbulenceStations3D ! for fort.43
      !WJP 02.20.2018: added fort.51-52 harmonic outputs
      type(stationData), private, save :: HAelevSta ! har elev stations (fort.51)
      type(stationData), private, save :: HAvelSta ! har vel stations (fort.52)

      !jgf: Water level correction data
      type(stationData), private, save :: dynamiccorrectionSta ! for dynamicwaterlevelcorrection.61
C
      TYPE, PRIVATE :: nodalData
         INTEGER ncid              ! the id of its netcdf file
         REAL(SZ) :: initial_value ! array will be initialized to this
         INTEGER :: nodal_data_id
         !WJP 02.20.2018: added fort.53-54 harmonic outputs
         INTEGER :: ha_nodal_data_id
         INTEGER :: hg_nodal_data_id
         INTEGER :: u_ha_nodal_data_id
         INTEGER :: u_hg_nodal_data_id
         INTEGER :: v_ha_nodal_data_id
         INTEGER :: v_hg_nodal_data_id
         INTEGER :: u_nodal_data_id
         INTEGER :: v_nodal_data_id
         INTEGER :: w_nodal_data_id
         INTEGER :: max_nodal_data_id
         INTEGER :: time_max_nodal_data_id        !tcm v51.20.06 added
         INTEGER :: nodal_data_dims(2)
         INTEGER :: nodal_data_dims_max(1)
         INTEGER :: nodal_data_dims_3D(3)
         CHARACTER(50) :: varnames(6) ! name(s) of this data in netcdf file
         type(timeData) :: myTime
         type(fileData) :: myFile
         type(meshStructure), pointer :: myMesh
      END TYPE nodalData
C
      type(nodalData), private, save :: elev       ! for fort.63
      type(nodalData), private, save :: pr         ! for fort.73
      type(nodalData), private, save :: wDrag      ! for windDragOut
      type(nodalData), private, save :: currentVel ! for fort.64
      type(nodalData), private, save :: windVel    ! for fort.74
      type(nodalData), private, save :: weirElev   ! for fort.77
      type(nodalData), private, save :: tau0nc     ! for fort.90
      type(nodalData), private, save, target :: EtaMax     ! for maxele.63
      type(nodalData), private, save, target :: PrMin      ! for minpr.63
      type(nodalData), private, save, target :: UMax       ! for maxvel.63
      type(nodalData), private, save, target :: WVMax      ! for maxwvel.63
      type(nodalData), private, save, target :: RSMax      ! for rads_max.63
      type(nodalData), private, save :: rads       ! for rads.64
C++ GML 20200528
      type(nodalData), private, save :: CICEAF     ! for fort.93
C++
#if defined CSWAN || defined ADCSWAN
      type(nodalData), private, save :: sw_hs      ! for swan_HS.63
      type(nodalData), private, save :: sw_dir     ! for swan_DIR.63
      type(nodalData), private, save :: sw_tm01    ! for swan_TM01.63
      type(nodalData), private, save :: sw_tps     ! for swan_TPS.63
      type(nodalData), private, save :: sw_wind    ! for swan_WIND.64
      type(nodalData), private, save :: sw_tm02    ! for swan_TM02.63
      type(nodalData), private, save :: sw_tmm10   ! for swan_tmm10.63
      type(nodalData), private, save :: sw_hs_max     ! for swan_HS_max.63
      type(nodalData), private, save :: sw_dir_max    ! for swan_DIR_max.63
      type(nodalData), private, save :: sw_tm01_max   ! for swan_TM01_max.63
      type(nodalData), private, save :: sw_tps_max    ! for swan_TPS_max.63
      type(nodalData), private, save :: sw_wind_max   ! for swan_WIND_max.63
      type(nodalData), private, save :: sw_tm02_max   ! for swan_TM02_max.63
      type(nodalData), private, save :: sw_tmm10_max  ! for swan_TMM10_max.63
#endif
      type(nodalData), private, save :: density3D    ! for fort.44
      type(nodalData), private, save :: velocity3D   ! for fort.45
      type(nodalData), private, save :: turbulence3D ! for fort.46
      type(nodalData), private, save :: futureSurfaceTemperature ! for fort.47
      !WJP 02.20.2018: added fort.53-54 harmonic outputs
      type(nodalData), private, save :: HAelev     ! for fort.53
      type(nodalData), private, save :: HAvel      ! for fort.54

      !jgf52.08.03: inundation data
      type(nodalData), private, save :: coldDry  ! for initiallydry.63
      type(nodalData), private, save, target :: maxInDep ! for maxinundepth.63
      type(nodalData), private, save, target :: inTime   ! for inundationtime.63
      type(nodalData), private, save :: eRisInun ! for endrisinginun.63
      type(nodalData), private, save, target :: evrDry ! for everdried.63

      !jgf: Water level correction data
      type(nodalData), private, save :: dynamiccorrection ! for dynamicwaterlevelcorrection.63
C
      TYPE, PRIVATE :: hotstartData
         INTEGER lun               ! ADCIRC's logical unit number
         INTEGER ncid              ! the id of its netcdf file
         INTEGER :: FileFmtMajorFile !
         INTEGER :: FileFmtMinorFile
         INTEGER :: FileFmtRevFile
         type(fileData) :: myFile
         type(timeData) :: myTime
         type(meshStructure), pointer :: myMesh
C        2D simulation state
         type(nodalData) :: zeta1
         type(nodalData) :: zeta2
         type(nodalData) :: zetad
         type(nodalData) :: vel
         type(nodalData) :: ch1
         type(nodalData) :: nodecodenc
         type(nodalData) :: noffnc
C        3D simulation state
         type(nodalData) :: duu
         type(nodalData) :: duv
         type(nodalData) :: dvv
         type(nodalData) :: uu
         type(nodalData) :: vv
         type(nodalData) :: bsx
         type(nodalData) :: bsy
         type(nodalData) :: density3D
         type(nodalData) :: velocity3D
         type(nodalData) :: turbulence3D
C        harmonic analysis components
         INTEGER :: namefr_len = 10
         INTEGER :: namefr_len_dim_id
         INTEGER :: namefr_dims(2)
         INTEGER :: mnharf_dim_id
         INTEGER :: load_vector_dim_id  ! 2x the number of frequencies
         INTEGER :: component_dims(1)
         INTEGER :: ha_dims(2)
         INTEGER :: hafreq_id
         INTEGER :: haff_id
         INTEGER :: haface_id
         INTEGER :: ha_id
         INTEGER :: namefr_id
C        harmonic analysis load vectors
         type(nodalData) :: gloelv
         type(nodalData) :: glovellv
         type(stationData) :: staelv
         type(stationData) :: stavellv
C        harmonic analysis means and variance calculations
         type(nodalData) :: xvelav
         type(nodalData) :: yvelav
         type(nodalData) :: xvelva
         type(nodalData) :: yvelva
         type(nodalData) :: elav
         type(nodalData) :: elva
      END TYPE hotstartData
C
      type(hotstartData), private, save, target :: hs67 ! for fort.67
      type(hotstartData), private, save, target :: hs68 ! for fort.68
      type(hotstartData), private, pointer :: hs   ! current hs file
      !
      ! lists of fulldomain node and element numbers that map to the
      ! nodes and elements of this subdomain in parallel
      INTEGER, ALLOCATABLE, TARGET :: fullDomainNodeList(:)
      INTEGER, ALLOCATABLE, TARGET :: fullDomainElementList(:)
      ! lists of fulldomain elevation and velocity station numbers that
      ! map to this subdomain's stations in parallel
      INTEGER, ALLOCATABLE, TARGET :: fullDomainElevationStationList(:)
      INTEGER, ALLOCATABLE, TARGET :: fullDomainVelocityStationList(:)
      ! used to specify the array of indexes to be used in a mapping
      INTEGER, POINTER :: fullDomainIndexList(:)
      ! usually initialized by readNetCDFHotstart(); otherwise in
      ! readAndMapToSubdomainMinMaxNetCDF()
      LOGICAL :: fullDomainIndexListsInitialized = .false.
C
      CHARACTER(1024) :: att_text ! scratch variable for attribute text
C
C     ZCobell:
      logical :: writerReadMetaData=.FALSE. !...variable to allow writers to read
                                            !   NetCDF metadata from ADCPREP processed
                                            !   files on first pass

      private initStationFile, initNodalDataFile,createNetCDFOutputFile,
     &   updateMetaData, initNetCDFCoord, defineMeshVariables,
     &   putMeshVariables, writeStationData,
     &   writeNodalData, setRecordCounterAndStoreTime,
     &   defineCoordinateAttributes, defineMetaData,
     &   defineTimeAttributes, getDimensions,
     &   mapFullDomainToSubdomain, mapFullDomainToSubdomainMByNP,
     &   mapFullDomainToSubdomainNPByM,
     &   putUnitsAttribute, metalength
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      CONTAINS
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   I N I T  N E T C D F  O U T P U T  F I L E
C-----------------------------------------------------------------------
C     jgf49.17.02 Allocates memory for NetCDF operations.
C-----------------------------------------------------------------------
      SUBROUTINE initNetCDFOutputFile(descript1, reterr, descript2,
     &   descript3)
      USE GLOBAL, ONLY : OutputDataDescript_t
#ifdef CMPI
      USE MESSENGER, ONLY : msg_fini
#endif
      IMPLICIT NONE
      TYPE(OutputDataDescript_t) :: descript1
      TYPE(OutputDataDescript_t), OPTIONAL :: descript2
      TYPE(OutputDataDescript_t), OPTIONAL :: descript3
      LOGICAL, intent(out) :: reterr
C
      call setMessageSource("initNetCDFOutputFile")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      reterr = .false.
      ! Don't allocate or initialize anything if this is a pure output file,
      ! (not an input/output file like max/min files) and the format is
      ! specified as something other than netcdf.
      IF ( ( (ABS(descript1 % specifier).ne.NETCDF3).and.
     &       (ABS(descript1 % specifier).ne.NETCDF4).and.
     &       (descript1 % readMaxMin .eqv..false.) ) .or.
           ! jgf52.08.14: If this is an input/output file (like a max/min
           ! file), we may still need to initialize the file so we can read
           ! the data from a previous run. But not if the file format specifier
           ! is OFF (i.e., 0) because won't be writing the file, so we don't
           ! need to initialize or read it.
     &       (descript1 % specifier.eq.0) ) then
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         RETURN
      ENDIF
C
      SELECT CASE(descript1 % lun)
         CASE(41)
            CALL initStationFile(densityStations3D,
     &           descript1, reterr)
         CASE(42)
            CALL initStationFile(velocityStations3D,
     &           descript1, reterr)
         CASE(43)
            CALL initStationFile(turbulenceStations3D,
     &           descript1, reterr)
         CASE(44)
            CALL initNodalDataFile(density3D,
     &           descript1, reterr, descript2, descript3)
         CASE(45)
            CALL initNodalDataFile(velocity3D,
     &           descript1, reterr, descript2, descript3)
         CASE(46)
            CALL initNodalDataFile(turbulence3D,
     &           descript1, reterr, descript2, descript3)
         CASE(47)
            CALL initNodalDataFile(futureSurfaceTemperature,
     &           descript1, reterr)
         !WJP 02.20.2018: added fort.51-54 harmonic outputs
         CASE(51)
            CALL initStationFile(HAelevSta, descript1, reterr)
         CASE(52)
            CALL initStationFile(HAvelSta, descript1, reterr)
         CASE(53)
            CALL initNodalDataFile(HAelev, descript1, reterr)
         CASE(54)
            CALL initNodalDataFile(HAvel, descript1, reterr)
         CASE(61)
            CALL initStationFile(elevSta, descript1, reterr)
         CASE(62)
            CALL initStationFile(velSta, descript1, reterr)
         CASE(63)
            CALL initNodalDataFile(elev, descript1, reterr)
         CASE(64)
            CALL initNodalDataFile(currentVel, descript1, reterr)
         CASE(71)
            CALL initStationFile(prSta, descript1, reterr)
         CASE(72)
            CALL initStationFile(wVelSta, descript1, reterr)
         CASE(73)
            CALL initNodalDataFile(pr, descript1, reterr)
         CASE(74)
            CALL initNodalDataFile(windVel, descript1, reterr)
C++ GML 20200528
         CASE (91)
            CALL initStationFile(IceSta, descript1, reterr)
         CASE (93)
            CALL initNodalDataFile(CICEAF, descript1, reterr)
C++
         CASE(173)
            CALL initNodalDataFile(wDrag, descript1, reterr)
         CASE(90)
            CALL initNodalDataFile(tau0nc, descript1, reterr)
         CASE(77)
            CALL initNodalDataFile(weirElev, descript1, reterr)
         CASE(164)
            CALL initNodalDataFile(rads, descript1, reterr)
#if defined CSWAN || defined ADCSWAN
         CASE(301)
            CALL initNodalDataFile(sw_hs, descript1, reterr)
         CASE(302)
            CALL initNodalDataFile(sw_dir, descript1, reterr)
         CASE(303)
            CALL initNodalDataFile(sw_tm01, descript1, reterr)
         CASE(304)
            CALL initNodalDataFile(sw_tps, descript1, reterr)
         CASE(305)
            CALL initNodalDataFile(sw_wind, descript1, reterr)
         CASE(306)
            CALL initNodalDataFile(sw_tm02, descript1, reterr)
         CASE(307)
            CALL initNodalDataFile(sw_tmm10, descript1, reterr)
#endif
         CASE(311)
            CALL initNodalDataFile(EtaMax, descript1, reterr)
         CASE(312)
            CALL initNodalDataFile(UMax, descript1, reterr)
         CASE(313)
            CALL initNodalDataFile(PrMin, descript1, reterr)
         CASE(314)
            CALL initNodalDataFile(WVMax, descript1, reterr)
         CASE(315)
            CALL initNodalDataFile(RSMax, descript1, reterr)
#if defined CSWAN || defined ADCSWAN
         CASE(316)
            CALL initNodalDataFile(sw_hs_max, descript1, reterr)
         CASE(317)
            CALL initNodalDataFile(sw_dir_max, descript1, reterr)
         CASE(318)
            CALL initNodalDataFile(sw_tm01_max, descript1, reterr)
         CASE(319)
            CALL initNodalDataFile(sw_tps_max, descript1, reterr)
         CASE(320)
            CALL initNodalDataFile(sw_wind_max, descript1, reterr)
         CASE(321)
            CALL initNodalDataFile(sw_tm02_max, descript1, reterr)
         CASE(322)
            CALL initNodalDataFile(sw_tmm10_max, descript1, reterr)
#endif
         CASE(400)
            CALL initNodalDataFile(inTime, descript1, reterr)
         CASE(401)
            CALL initNodalDataFile(maxInDep, descript1, reterr)
         CASE(402)
            CALL initNodalDataFile(coldDry, descript1, reterr)
         CASE(403)
            CALL initNodalDataFile(eRisInun, descript1, reterr)
         CASE(404)
            CALL initNodalDataFile(evrDry, descript1, reterr)
         CASE(108)
            CALL initNodalDataFile(dynamiccorrection, descript1, reterr)
         CASE(109)
            CALL initStationFile(dynamiccorrectionSta, descript1, reterr)

         CASE DEFAULT
            write(scratchMessage,
     &       '("No netCDF for files with unit number ",i0,".")') descript1 % lun
            call allMessage(ERROR,scratchMessage)
      END SELECT
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE initNetCDFOutputFile
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E    I N I T  S T A T I O N   F I L E
C-----------------------------------------------------------------------
C     jgf49.17.02 Sets up netCDF variables and allocates memory for a
C     station.
C-----------------------------------------------------------------------
      SUBROUTINE initStationFile(sta, descript1, reterr)
      USE SIZES, ONLY : MNPROC
      USE GLOBAL, ONLY : SNAMLEN, OutputDataDescript_t, RAD2DEG,
     &                   IHOT, NWS, STATNAME, STATNAMEV,
     &                   STATNAMEM, IDEN
      USE GLOBAL_3DVS, ONLY : STATNAMED, STATNAMEV3D, STATNAMET
      USE MESH, ONLY : ICS
#ifdef CMPI
      USE MESSENGER, ONLY : msg_fini
#endif
      USE HARM, ONLY: NAMEFR
      IMPLICIT NONE
C
      type(stationData), intent(inout) :: sta
      type(OutputDataDescript_t), intent(inout) :: descript1
      LOGICAL, intent(out) :: reterr
C
      INTEGER i,j
      INTEGER iret  ! success or failure of the netcdf call

C
!     date_string variables for time attribute
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer values(8)
      integer :: stationLuns3D(3) = (/ 41, 42, 43 /)
      integer :: stationLunsHA(2) = (/ 51, 52 /)
C
      call setMessageSource("initStationFile")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      reterr = .false.
C
C     jgf50.13: if netcdf output was requested, but there are no stations,
C     don't initialize the file; just return.
      IF ( descript1 % num_fd_records.eq.0 ) THEN
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         RETURN
      ENDIF
C
      sta % num_stations = descript1 % num_fd_records
      ALLOCATE(sta%myTime%timenc(sta%myTime%timenc_len))
C
      if (any(stationLuns3D.eq.descript1 % lun)
     &    .or.any(stationLunsHA.eq.descript1 % lun)) then
         sta%num_v_nodes = descript1 % num_items_per_record
      ENDIF
C
C     Memory allocation for the station
      ALLOCATE(sta%x(sta%num_stations))
      ALLOCATE(sta%y(sta%num_stations))
C
C     Initialize netCDF output file, creating a new one if necessary.
      CALL createNetCDFOutputFile(sta%ncid, sta%myFile, sta%myTime,
     &                            descript1, reterr)
C
C     if we didn't need to create a file, update metadata and return
      !jgf52.08.22: Only update metadata if running in serial; in
      !parallel the metadata was updated by adcprep
      IF ( sta%myFile%createFile.eqv..false.) THEN
         if ( (reterr.eqv..false.).and.(mnproc.eq.1) ) then
            CALL updateMetaData(sta%ncid,sta%myFile)
         endif
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         RETURN
      ENDIF
C
C     Set coordinates of each station, converting to degrees if we
C     are in spherical coordinates.
      IF(ICS.NE.1) THEN
         sta%x = descript1%x_coord * RAD2DEG
         sta%y = descript1%y_coord * RAD2DEG
      ELSE
         sta%x = descript1%x_coord
         sta%y = descript1%y_coord
      ENDIF
C
C     Create station dimension and station name dimension
      iret = nf90_def_dim(sta%ncid, 'station',
     &       sta%num_stations, sta%num_sta_dim_id)
      CALL check_err(iret)
      iret = nf90_def_dim(sta%ncid, 'namelen', SNAMLEN, sta%slen_dim_id)
      CALL check_err(iret)
      IF (any(stationLuns3D.eq.descript1%lun)) THEN
         iret = nf90_def_dim(sta%ncid, 'num_v_nodes', sta%num_v_nodes,
     &       sta%num_v_nodes_dim_id)
         sta%station_data_dims_3D(1) = sta%num_sta_dim_id
         sta%station_data_dims_3D(2) = sta%num_v_nodes_dim_id
         sta%station_data_dims_3D(3) = sta%myTime%timenc_dim_id
      !WJP 02.20.2018 Writing in constituent information for harmonic file
      ELSEIF (any(stationLunsHA.eq.descript1 % lun)) THEN
          call defineHarmonicAnalysisParametersInNetcdfFile(descript1,stationdat=sta)
      ENDIF
      CALL check_err(iret)
!
!     Define stations name
      sta%station_dims(1) = sta%slen_dim_id
      sta%station_dims(2) = sta%num_sta_dim_id

      iret = nf90_def_var(sta%ncid, 'station_name', NF90_CHAR,
     &    sta%station_dims, sta%station_id)
      CALL check_err(iret)
C
C     Define station locations
      sta%x_dims(1) = sta%num_sta_dim_id
      iret = nf90_def_var(sta%ncid, 'x', NF90_DOUBLE,
     &       sta%x_dims, sta%x_id)
      CALL check_err(iret)
      sta%y_dims(1) = sta%num_sta_dim_id
      iret = nf90_def_var(sta%ncid, 'y', NF90_DOUBLE,
     &       sta%y_dims, sta%y_id)
      CALL check_err(iret)
C
C     Set coordinates as representing either latitude or longitude,
C     or Cartesian x and y, depending on the value of ICS.
      CALL defineCoordinateAttributes(sta%ncid, sta%x_id, sta%y_id)
C
C     Fill in labels and populate variables as appropriate for the
C     different types of data in the station files. The labels and
C     units will also vary according to the coordinate system ADCIRC
C     is using (spherical or cartesian, according to the value of ICS)
C     as well as the units system (english or si according to the value of g).
      SELECT CASE(descript1 % lun)
         CASE(41) !       F O R T . 4 1
            iret = nf90_def_var(sta%ncid, 'sigmat', NF90_DOUBLE,
     &             sta%station_data_dims_3D, sta%u_station_data_id)
            CALL check_err(iret)
            IF ((IDEN.eq.2).or.(IDEN.eq.4)) THEN
               iret = nf90_def_var(sta%ncid, 'salinity', NF90_DOUBLE,
     &             sta%station_data_dims_3D, sta%v_station_data_id)
               CALL check_err(iret)
            ENDIF
            IF ((IDEN.eq.3).or.(IDEN.eq.4)) THEN
               iret = nf90_def_var(sta%ncid, 'temperature', NF90_DOUBLE,
     &             sta%station_data_dims_3D, sta%w_station_data_id)
               CALL check_err(iret)
            ENDIF
            ! sigma t
            att_text =
     &         'station water column vertically varying density'
            iret = nf90_put_att(sta%ncid,
     &         sta%u_station_data_id, 'long_name', trim(att_text))
            CALL check_err(iret)
            att_text = 'station_density_vertically_varying'
            iret = nf90_put_att(sta%ncid,
     &             sta%u_station_data_id, 'standard_name',
     &             trim(att_text))
            CALL check_err(iret)
            CALL putUnitsAttribute(sta%ncid,sta%u_station_data_id,
     &      "kg/m^3", "n/a")
            ! salinity
            IF ((IDEN.eq.2).or.(IDEN.eq.4)) THEN
               att_text =
     &            'station water column vertically varying salinity'
               iret = nf90_put_att(sta%ncid,
     &            sta%v_station_data_id, 'long_name',trim(att_text))
               CALL check_err(iret)
               att_text =
     &            'station_water_salinity_vertically_varying'
               iret = nf90_put_att(sta%ncid,
     &            sta%v_station_data_id, 'standard_name',
     &            trim(att_text))
               CALL check_err(iret)
               CALL putUnitsAttribute(sta%ncid,sta%v_station_data_id,
     &            "PSU","n/a")
            ENDIF
            ! temperature
            IF ((IDEN.eq.3).or.(IDEN.eq.4)) THEN
               att_text =
     &            'station water column vertically varying temperature'
               iret = nf90_put_att(sta%ncid,
     &            sta%w_station_data_id, 'long_name',
     &            trim(att_text))
               CALL check_err(iret)
               att_text =
     &            'station_water_temperature_vertically_varying'
               iret = nf90_put_att(sta%ncid,
     &            sta%w_station_data_id, 'standard_name',
     &            trim(att_text))
               CALL check_err(iret)
               CALL putUnitsAttribute(sta%ncid,sta%w_station_data_id,
     &            "Celsius","Fahrenheit")
            ENDIF
            sta%statnames => STATNAMED
C
         CASE(42) !       F O R T . 4 2
            iret = nf90_def_var(sta%ncid, 'u-vel3D', NF90_DOUBLE,
     &             sta%station_data_dims_3D, sta%u_station_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(sta%ncid, 'v-vel3D', NF90_DOUBLE,
     &             sta%station_data_dims_3D, sta%v_station_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(sta%ncid, 'w-vel3D', NF90_DOUBLE,
     &             sta%station_data_dims_3D,  sta%w_station_data_id)
            CALL check_err(iret)
            if (ics.ne.1) then
               ! u
               att_text =
     &    'station water column vertically varying east/west velocity'
               iret = nf90_put_att(sta%ncid,
     &         sta%u_station_data_id, 'long_name', trim(att_text))
               CALL check_err(iret)
               att_text =
     &     'station_eastward_water_velocity_vertically_varying'
               iret = nf90_put_att(sta%ncid,
     &             sta%u_station_data_id,
     &             'standard_name', trim(att_text))
               CALL check_err(iret)
               iret = nf90_put_att(sta%ncid,
     &             sta%u_station_data_id, 'positive', 'east')
               CALL check_err(iret)
               ! v
               att_text =
     &  'station water column vertically varying north/south velocity'
               iret = nf90_put_att(sta%ncid,
     &            sta%v_station_data_id, 'long_name',
     &            trim(att_text))
               CALL check_err(iret)
               att_text =
     &            'station_northward_water_velocity_vertically_varying'
               iret = nf90_put_att(sta%ncid,
     &            sta%v_station_data_id, 'standard_name',
     &            trim(att_text))
               CALL check_err(iret)
               iret = nf90_put_att(sta%ncid,
     &         sta%v_station_data_id, 'positive', 'north')
               CALL check_err(iret)
               ! w
               att_text =
     &  'station water column vertically varying up/down velocity'
               iret = nf90_put_att(sta%ncid,
     &            sta%w_station_data_id, 'long_name',
     &            trim(att_text))
               CALL check_err(iret)
               att_text =
     &            'station_upward_water_velocity_vertically_varying'
               iret = nf90_put_att(sta%ncid,
     &            sta%w_station_data_id, 'standard_name',
     &            trim(att_text))
               CALL check_err(iret)
               iret = nf90_put_att(sta%ncid,
     &         sta%w_station_data_id, 'positive', 'up')
               CALL check_err(iret)
            else
               ! u
               att_text =
     & 'station water column vertically varying velocity in x-direction'
               iret = nf90_put_att(sta%ncid,
     &            sta%u_station_data_id,'long_name', trim(att_text))
               CALL check_err(iret)
               att_text = 'station_x_water_velocity_vertically_varying'
               iret = nf90_put_att(sta%ncid,
     &            sta%u_station_data_id, 'standard_name',
     &            trim(att_text))
               CALL check_err(iret)
               iret = nf90_put_att(sta%ncid,
     &            sta%u_station_data_id, 'positive', 'right')
               CALL check_err(iret)
               ! v
               att_text =
     & 'station water column vertically varying velocity in y-direction'
               iret = nf90_put_att(sta%ncid,
     &            sta%v_station_data_id, 'long_name', trim(att_text))
               CALL check_err(iret)
               att_text = 'station_y_water_velocity_vertically_varying'
               iret = nf90_put_att(sta%ncid,
     &            sta%v_station_data_id, 'standard_name',
     &            trim(att_text))
               CALL check_err(iret)
               iret = nf90_put_att(sta%ncid,
     &            sta%v_station_data_id, 'positive',
     &            '90 degrees counterclockwise from x water velocity')
               CALL check_err(iret)
               ! w
               att_text = 'station water column '
     &             //'vertically varying velocity in z-direction'
               iret = nf90_put_att(sta%ncid,
     &            sta%w_station_data_id, 'long_name',trim(att_text))
               CALL check_err(iret)
               att_text = 'station_z_water_velocity_vertically_varying'
               iret = nf90_put_att(sta%ncid,
     &            sta%w_station_data_id, 'standard_name',
     &            trim(att_text))
               CALL check_err(iret)
               iret = nf90_put_att(sta%ncid,
     &            sta%w_station_data_id, 'positive',
     &            '90 degrees counterclockwise from x water velocity')
               CALL check_err(iret)
            endif
            CALL putUnitsAttribute(sta%ncid,sta%u_station_data_id,
     &      'm s-1', 'ft s-1')
            CALL putUnitsAttribute(sta%ncid,sta%v_station_data_id,
     &      'm s-1','ft s-1')
            CALL putUnitsAttribute(sta%ncid,sta%w_station_data_id,
     &      'm s-1','ft s-1')
            sta%statnames => STATNAMEV3D

         CASE(43) !       F O R T . 4 3
            iret = nf90_def_var(sta%ncid, 'q20', NF90_DOUBLE,
     &          sta%station_data_dims_3D, sta%u_station_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(sta%ncid, 'l', NF90_DOUBLE,
     &          sta%station_data_dims_3D, sta%v_station_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(sta%ncid, 'ev', NF90_DOUBLE,
     &          sta%station_data_dims_3D, sta%w_station_data_id)
            CALL check_err(iret)
            ! q20
            att_text =
     &      'station water column vertically varying turbulent'//
     &      'kinetic energy'
            iret = nf90_put_att(sta%ncid,
     &         sta%u_station_data_id, 'long_name', trim(att_text))
            CALL check_err(iret)
            att_text =
     &         'station_turbulent_kinetic_energy_vertically_varying'
            iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &             'standard_name', trim(att_text))
            CALL check_err(iret)
            ! l
            att_text =
     &         'station water column vertically varying mixing length'
            iret = nf90_put_att(sta%ncid,
     &            sta%v_station_data_id, 'long_name',
     &            trim(att_text))
            CALL check_err(iret)
            att_text =
     &            'station_water_mixing_length_vertically_varying'
               iret = nf90_put_att(sta%ncid,
     &            sta%v_station_data_id, 'standard_name',
     &            trim(att_text))
               CALL check_err(iret)
            !  ev
            att_text =
     &        'station water column vertically varying eddy viscosity'
               iret = nf90_put_att(sta%ncid,
     &            sta%w_station_data_id, 'long_name',
     &            trim(att_text))
               CALL check_err(iret)
               att_text =
     &            'station_water_eddy_viscosity_vertically_varying'
               iret = nf90_put_att(sta%ncid,
     &            sta%w_station_data_id, 'standard_name',
     &            trim(att_text))
               CALL check_err(iret)
C            CALL putUnitsAttribute(sta%ncid,sta%u_station_data_id,
C     &      'm s-1', 'ft s-1') !TODO: jgf49.48.01: units for q20??
            CALL putUnitsAttribute(sta%ncid,sta%v_station_data_id,
     &      "meters","n/a")
            CALL putUnitsAttribute(sta%ncid,sta%w_station_data_id,
     &      "m^2/s","n/a")
            sta%statnames => STATNAMET

         CASE(51) !       F O R T . 5 1        WJP 02.20.2018

            sta%station_data_dims(2) = sta%num_sta_dim_id
            sta%station_data_dims(1) = sta%num_v_nodes_dim_id
!           Put the amplitude stuff in
            iret = nf90_def_var(sta%ncid, 'amp', NF90_DOUBLE,
     &             sta%station_data_dims, sta%station_ha_data_id)
            CALL check_err(iret)
!           Define water surface elevation attributes
            iret = nf90_put_att(sta%ncid, sta%station_ha_data_id,
     &             'long_name',
     &             'station amplitude of tidal harmonic constituents')
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%station_ha_data_id,
     &             'standard_name', 'sta_hc_amp')
            CALL check_err(iret)
            CALL putUnitsAttribute(sta%ncid, sta%station_ha_data_id,
     &           'm', 'ft')
            iret = nf90_put_att(sta%ncid, sta%station_ha_data_id,
     &         '_FillValue',doubleval)
            CALL check_err(iret)
!           Put the phase stuff in
            iret = nf90_def_var(sta%ncid, 'phs', NF90_DOUBLE,
     &             sta%station_data_dims, sta%station_hg_data_id)
            CALL check_err(iret)
!           Define water surface elevation attributes
            iret = nf90_put_att(sta%ncid, sta%station_hg_data_id,
     &             'long_name',
     &             'station phase of tidal harmonic constituents')
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%station_hg_data_id,
     &             'standard_name', 'sta_hc_phs')
            CALL check_err(iret)
            CALL putUnitsAttribute(sta%ncid, sta%station_hg_data_id,
     &           'deg','deg')
            iret = nf90_put_att(sta%ncid, sta%station_hg_data_id,
     &         '_FillValue',doubleval)
            CALL check_err(iret)

            sta%statnames => STATNAME

         CASE(52) !       F O R T . 5 2        WJP 02.20.2018

            sta%station_data_dims(2) = sta%num_sta_dim_id
            sta%station_data_dims(1) = sta%num_v_nodes_dim_id
            iret = nf90_def_var(sta%ncid, 'u-vel-amp', NF90_DOUBLE,
     &             sta%station_data_dims, sta%ha_u_station_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(sta%ncid, 'u-vel-phs', NF90_DOUBLE,
     &             sta%station_data_dims, sta%hg_u_station_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(sta%ncid, 'v-vel-amp', NF90_DOUBLE,
     &             sta%station_data_dims, sta%ha_v_station_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(sta%ncid, 'v-vel-phs', NF90_DOUBLE,
     &             sta%station_data_dims, sta%hg_v_station_data_id)
            CALL check_err(iret)
               ! u amp
            iret = nf90_put_att(sta%ncid, sta%ha_u_station_data_id,
     &      'long_name',
     &  'station tidal harmonic amplitude of east/west velocity')
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%ha_u_station_data_id,
     &                   'standard_name',
     &       'sta_hc_amp_east_vel')
            CALL check_err(iret)
            CALL putUnitsAttribute(sta%ncid,sta%ha_u_station_data_id,
     &           'm/s', 'ft/s')
            iret = nf90_put_att(sta%ncid, sta%ha_u_station_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%ha_u_station_data_id,
     &             'dry_Value', doubleval)
            CALL check_err(iret)
               ! u phs
            iret = nf90_put_att(sta%ncid, sta%hg_u_station_data_id,
     &         'long_name',
     &  'station tidal harmonic phase of east/west velocity')
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%hg_u_station_data_id,
     &                   'standard_name',
     &       'sta_hc_phs_east_vel')
            CALL check_err(iret)
            CALL putUnitsAttribute(sta%ncid,sta%hg_u_station_data_id,
     &           'deg','deg')
            iret = nf90_put_att(sta%ncid, sta%hg_u_station_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%hg_u_station_data_id,
     &             'dry_Value', doubleval)
            CALL check_err(iret)
               ! v amp
            iret = nf90_put_att(sta%ncid,sta%ha_v_station_data_id,
     &         'long_name',
     &  'station tidal harmonic amplitude of north/south velocity')
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%ha_v_station_data_id,
     &         'standard_name',
     &         'sta_hc_amp_north_vel')
            CALL check_err(iret)
            CALL putUnitsAttribute(sta%ncid,sta%ha_v_station_data_id,
     &           'm/s','ft/s')
            iret = nf90_put_att(sta%ncid, sta%ha_v_station_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%ha_v_station_data_id,
     &             'dry_Value', doubleval)
            CALL check_err(iret)
               ! v phs
            iret = nf90_put_att(sta%ncid,sta%hg_v_station_data_id,
     &         'long_name',
     &  'station tidal harmonic phase of north/south velocity')
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%hg_v_station_data_id,
     &         'standard_name',
     &         'sta_hc_phs_north_vel')
            CALL check_err(iret)
            CALL putUnitsAttribute(sta%ncid,sta%hg_v_station_data_id,
     &           'deg','deg')
            iret = nf90_put_att(sta%ncid, sta%hg_v_station_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%hg_v_station_data_id,
     &             'dry_Value', doubleval)
            CALL check_err(iret)

            sta%statnames => STATNAMEV

         CASE(61) !       F O R T . 6 1

            sta%station_data_dims(1) = sta%num_sta_dim_id
            sta%station_data_dims(2) = sta%myTime%timenc_dim_id
            iret = nf90_def_var(sta%ncid, 'zeta', NF90_DOUBLE,
     &             sta%station_data_dims, sta%station_data_id)
            CALL check_err(iret)
!           Define water surface elevation attributes
            iret = nf90_put_att(sta%ncid, sta%station_data_id,
     &             'long_name', 'water surface elevation above geoid')
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%station_data_id,
     &             'standard_name', 'sea_surface_height_above_geoid')
            CALL check_err(iret)
            CALL putUnitsAttribute(sta%ncid, sta%station_data_id,
     &           'm', 'ft')
            iret = nf90_put_att(sta%ncid, sta%station_data_id,
     &         '_FillValue',doubleval)
            CALL check_err(iret)
!            iret = nf90_put_att(sta%ncid, sta%station_data_id,
!     &             'positive', 'up')
!            CALL check_err(iret)
            sta%statnames => STATNAME

         CASE(62) !       F O R T . 6 2

            sta%station_data_dims(1) = sta%num_sta_dim_id
            sta%station_data_dims(2) = sta%myTime%timenc_dim_id
            iret = nf90_def_var(sta%ncid, 'u-vel', NF90_DOUBLE,
     &             sta%station_data_dims, sta%u_station_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(sta%ncid, 'v-vel', NF90_DOUBLE,
     &             sta%station_data_dims, sta%v_station_data_id)
            CALL check_err(iret)
            if (ics.ne.1) then
               ! u
               iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &         'long_name',
     &    'station water column vertically averaged east/west velocity')
               CALL check_err(iret)

               iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &                   'standard_name',
     &       'station_eastward_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &                'positive', 'east')
               CALL check_err(iret)
               ! v
               iret = nf90_put_att(sta%ncid,sta%v_station_data_id,
     &         'long_name',
     & 'station water column vertically averaged north/south velocity')
               CALL check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &         'standard_name',
     &         'station_northward_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &         'positive', 'north')
               CALL check_err(iret)
            else
               ! u
               iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &               'long_name', 'station water column vertically '
     &               //'averaged velocity in x-direction')
               CALL check_err(iret)

               iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &             'standard_name',
     &             'station_x_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &                'positive', 'right')
               CALL check_err(iret)
               ! v
               iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &                'long_name',
     &            'station water column vertically averaged velocity '
     &            //'in y-direction')
               CALL check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &         'standard_name',
     &         'station_y_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &         'positive',
     &         '90 degrees counterclockwise from x water velocity')
               CALL check_err(iret)
            endif

            CALL putUnitsAttribute(sta%ncid, sta%u_station_data_id,
     &      'm s-1', 'ft s-1')
            iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &             'dry_Value', doubleval)
            CALL check_err(iret)
            CALL putUnitsAttribute(sta%ncid, sta%v_station_data_id,
     &      'm s-1','ft s-1')
            iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &      '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &      'dry_Value', doubleval)
            CALL check_err(iret)
            sta%statnames => STATNAMEV

         CASE(71) !         F  O R T . 7 1

            sta%station_data_dims(1) = sta%num_sta_dim_id
            sta%station_data_dims(2) = sta%myTime%timenc_dim_id

            iret = nf90_def_var(sta%ncid, 'pressure', NF90_DOUBLE,
     &         sta%station_data_dims, sta%station_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%station_data_id,
     &         'long_name', 'station air pressure at sea level')
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%station_data_id,
     &         'standard_name', 'station_air_pressure_at_sea_level')
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%station_data_id,
     &         'units', 'meters of water')
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%station_data_id,
     &         '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%station_data_id,
     &         'positive', 'up')
            CALL check_err(iret)
            sta%statnames => STATNAMEM

         CASE(72)  !      F O R T . 7 2

            sta%station_data_dims(1) = sta%num_sta_dim_id
            sta%station_data_dims(2) = sta%myTime%timenc_dim_id
            iret = nf90_def_var(sta%ncid, 'windx', NF90_DOUBLE,
     &         sta%station_data_dims, sta%u_station_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(sta%ncid, 'windy', NF90_DOUBLE,
     &         sta%station_data_dims, sta%v_station_data_id)
            CALL check_err(iret)
C
            if (ics.ne.1) then
               select case(abs(nws))
                  case(1,2)
                    iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &                      'long_name', 'station e/w wind stress')
                    CALL check_err(iret)
                    iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &                     'standard_name',
     &                     'station_eastward_surface_wind_stress')
                    CALL check_err(iret)
                    iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &                     'positive', 'east')
                    iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &                     'long_name', 'station n/s wind stress')
                    CALL check_err(iret)
                    iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &                    'standard_name',
     &                    'station_northward_surface_wind_stress')
                    CALL check_err(iret)
                    iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &                     'positive', 'north')
                    CALL check_err(iret)
                  case default
                    iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &                     'long_name', 'station e/w wind velocity')
                    CALL check_err(iret)
                    iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &                     'standard_name',
     &                     'station_eastward_wind')
                    CALL check_err(iret)
                    iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &                     'positive', 'east')
                    iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &                     'long_name',
     &                     'station n/s wind velocity')
                    CALL check_err(iret)
                    iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &                     'standard_name', 'station_northward_wind')
                    CALL check_err(iret)
                    iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &                     'positive', 'north')
               end select
            else
               select case(abs(nws))
                  case(1,2)
                    iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &                     'long_name',
     &                     'station wind stress in x-direction')
                    CALL check_err(iret)
                    iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &                     'standard_name',
     &                     'station_x_surface_wind_stress')
                    CALL check_err(iret)
                    iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &                     'positive', 'right')
                    iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &                     'long_name',
     &                     'station wind stress in y-direction')
                    CALL check_err(iret)
                    iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &                     'standard_name',
     &                     'station_y_surface_wind_stress')
                   CALL check_err(iret)
                    iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &                    'positive',
     & '90 degrees counterclockwise from wind velocity in x-direction')
                  case default
                    iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &                      'long_name',
     &                      'station wind velocity in x-direction')
                    CALL check_err(iret)
                    iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &                     'standard_name', 'station_x_wind')
                    CALL check_err(iret)
                    iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &                     'positive', 'right')
                    iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &                     'long_name',
     &                     'station wind velocity in y-direction')
                    CALL check_err(iret)
                    iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &                     'standard_name',
     &                     'station_y_wind')
                    CALL check_err(iret)
                    iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &                     'positive',
     &  '90 degrees counterclockwise from wind velocity in x-direction')
                     CALL check_err(iret)
               end select
            endif

            if (abs(nws).gt.2.and.abs(nws).lt.100) then
               call putUnitsAttribute(sta%ncid, sta%u_station_data_id,
     &              'm s-1', 'ft s-1')
               call putUnitsAttribute(sta%ncid, sta%v_station_data_id,
     &              'm s-1', 'ft s-1')
            else
               call putUnitsAttribute(sta%ncid, sta%u_station_data_id,
     &              'm2 s-2', 'ft s-2')
               call putUnitsAttribute(sta%ncid, sta%u_station_data_id,
     &              'm2 s-2', 'ft s-2')
            endif

            iret = nf90_put_att(sta%ncid, sta%u_station_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%v_station_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            sta%statnames => STATNAMEM

C++ GML 20200605 
         CASE(91) !         F  O R T . 9 1

            sta%station_data_dims(1) = sta%num_sta_dim_id
            sta%station_data_dims(2) = sta%myTime%timenc_dim_id

            iret = nf90_def_var(sta%ncid, 'iceaf', NF90_DOUBLE,
     &         sta%station_data_dims, sta%station_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%station_data_id,
     &         'long_name', 'station sea-ice af')
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%station_data_id,
     &         'standard_name', 'station_sea-ice area fraction')
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%station_data_id,
     &         'units', 'unitless')
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%station_data_id,
     &         '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%station_data_id,
     &         'positive', 'up')
            CALL check_err(iret)
            sta%statnames => STATNAMEM
c ++ GML

         CASE(109) !       O F F S E T . 6 1

            sta%station_data_dims(1) = sta%num_sta_dim_id
            sta%station_data_dims(2) = sta%myTime%timenc_dim_id
            iret = nf90_def_var(sta%ncid, 'dynamicWaterlevelCorrection', NF90_DOUBLE,
     &             sta%station_data_dims, sta%station_data_id)
            CALL check_err(iret)
!           Define water surface elevation attributes
            iret = nf90_put_att(sta%ncid, sta%station_data_id,
     &         'long_name',
     &         'dynamic water surface correction above water level')
            CALL check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%station_data_id,
     &         'standard_name',
     &         'dynamic_sea_surface_correction_above_water_level')
            CALL check_err(iret)
            CALL putUnitsAttribute(sta%ncid, sta%station_data_id,
     &           'm', 'ft')
            iret = nf90_put_att(sta%ncid, sta%station_data_id,
     &         '_FillValue',doubleval)
            CALL check_err(iret)
            sta%statnames => STATNAME


         CASE DEFAULT
            write(scratchMessage,
     &       '("No netCDF for station files with unit number ",i0,".")') descript1 % lun
            call allMessage(ERROR,scratchMessage)
      END SELECT
      iret = nf90_put_att(sta%ncid,sta%u_station_data_id,
     &      '_FillValue', doubleval)
      CALL check_err(iret)
      iret = nf90_put_att(sta%ncid,sta%u_station_data_id,
     &      'dry_Value', doubleval)
      CALL check_err(iret)
      iret = nf90_put_att(sta%ncid,sta%v_station_data_id,
     &      '_FillValue',  doubleval)
      CALL check_err(iret)
      iret = nf90_put_att(sta%ncid,sta%v_station_data_id,
     &      'dry_Value', doubleval)
      CALL check_err(iret)
      iret = nf90_put_att(sta%ncid,sta%w_station_data_id,
     &      '_FillValue',  doubleval)
      CALL check_err(iret)
      iret = nf90_put_att(sta%ncid,sta%w_station_data_id,
     &      'dry_Value',  doubleval)
      CALL check_err(iret)
C
C     jgf50.44: Automatically turn on compression if we are using the
C     netcdf4 file format.
#ifdef NETCDF_CAN_DEFLATE
      IF (abs(descript1%specifier).eq.5) THEN
         SELECT CASE(descript1 % lun)
         CASE(41)
            iret = nf_def_var_deflate(sta%ncid, sta%u_station_data_id,
     &          1, 1, 2)
            CALL check_err(iret)
            IF ((IDEN.eq.2).or.(IDEN.eq.4)) THEN
               iret = nf_def_var_deflate(sta%ncid,
     &                sta%v_station_data_id, 1, 1, 2)
               CALL check_err(iret)
            ENDIF
            IF ((IDEN.eq.3).or.(IDEN.eq.4)) THEN
               iret = nf_def_var_deflate(sta%ncid,
     &                sta%w_station_data_id, 1, 1, 2)
               CALL check_err(iret)
            ENDIF
         CASE(42,43)
            iret = nf_def_var_deflate(sta%ncid, sta%u_station_data_id,
     &          1, 1, 2)
            CALL check_err(iret)
            iret = nf_def_var_deflate(sta%ncid,
     &             sta%v_station_data_id, 1, 1, 2)
            CALL check_err(iret)
            iret = nf_def_var_deflate(sta%ncid,
     &             sta%w_station_data_id, 1, 1, 2)
            CALL check_err(iret)
         CASE(51)
            CALL check_err(nf90_def_var_deflate(sta%ncid, sta%station_ha_data_id,
     &             1, 1, 2))
            CALL check_err(nf90_def_var_deflate(sta%ncid, sta%station_hg_data_id,
     &             1, 1, 2))
         CASE(52)
            CALL check_err(nf90_def_var_deflate(sta%ncid, sta%ha_u_station_data_id,
     &             1, 1, 2))
            CALL check_err(nf90_def_var_deflate(sta%ncid, sta%hg_u_station_data_id,
     &             1, 1, 2))
            CALL check_err(nf90_def_var_deflate(sta%ncid, sta%ha_v_station_data_id,
     &             1, 1, 2))
            CALL check_err(nf90_def_var_deflate(sta%ncid, sta%hg_v_station_data_id,
     &             1, 1, 2))
         CASE(61,71,91,109) ! GML added 91 20210727
            iret = nf90_def_var_deflate(sta%ncid, sta%station_data_id,
     &             1, 1, 2)
            CALL check_err(iret)
         CASE(62,72)
            iret = nf90_def_var_deflate(sta%ncid, sta%u_station_data_id,
     &             1, 1, 2)
            CALL check_err(iret)
            iret = nf90_def_var_deflate(sta%ncid, sta%v_station_data_id,
     &             1, 1, 2)
            CALL check_err(iret)
         CASE DEFAULT ! should be unreachable
            write(scratchMessage,
     &       '("No netCDF for station files with unit number ",i0,".")') descript1 % lun
            call allMessage(ERROR,scratchMessage)
         END SELECT
      ENDIF
#endif
!
!     Define station names and code attributes
      iret = nf90_put_att(sta%ncid,sta%station_id,'long_name',
     &       'station name')
      CALL check_err(iret)
!
!     Define time attributes
      CALL defineTimeAttributes(sta%ncid, sta%myTime)
C
C     define metadata and selected fort.15 parameters in netcdf file
      CALL defineMetaData(sta%ncid)
!
!     Leave define mode
      iret = nf90_enddef(sta%ncid)
      CALL check_err(iret)
!
!     Store station name(s)
      iret = nf90_put_var(sta%ncid, sta%station_id, sta%statnames(:) )
      CALL check_err(iret)
C
C     Store station locations
      iret = nf90_put_var(sta%ncid, sta%x_id, sta%x)
      CALL check_err(iret)
      iret = nf90_put_var(sta%ncid, sta%y_id, sta%y)
      CALL check_err(iret)

      !WJP 02.20.2018 Store const name(s)
      IF (any(stationLunsHA.eq.descript1 % lun)) THEN
        CALL addHarmonicAnalysisParametersToNetcdfFile(sta%ncid)
      ENDIF

C
C     now close the initialized netcdf file
      iret = nf90_close(sta%ncid)
      CALL check_err(iret)
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE initStationFile
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E    I N I T  N O D A L  D A T A   F I L E
C-----------------------------------------------------------------------
C     jgf49.17.02 Sets up netCDF variables and allocates memory for
C     full domain data.
C-----------------------------------------------------------------------
      SUBROUTINE initNodalDataFile(dat, descript1, reterr,
     &   descript2, descript3)
      USE SIZES, ONLY : MNWPROC, MNPROC
      USE GLOBAL, ONLY : OutputDataDescript_t, NWS, C3D, IDEN, NE_G, NP_G
      USE GLOBAL_3DVS, ONLY : NFEN
      USE MESH, ONLY : ICS
      USE HARM, ONLY: NAMEFR
      IMPLICIT NONE
C
      type(nodalData), intent(inout) :: dat
      type(OutputDataDescript_t), intent(inout) :: descript1
      type(OutputDataDescript_t), intent(inout), optional :: descript2
      type(OutputDataDescript_t), intent(inout), optional :: descript3
      LOGICAL, intent(out) :: reterr
C
      INTEGER n    ! node counter
      INTEGER iret ! success or failure of netcdf call
      integer :: nodalLunsHA(2) = (/ 53, 54 /)
      REAL(SZ),ALLOCATABLE :: defaultValue(:)
      CHARACTER(len=1024) :: att_text  ! metadata
C
      call setMessageSource("initNodalDataFile")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      reterr = .false.
      dat%myMesh => adcircMesh
      IF (dat%myMesh%initialized .eqv..false.) THEN
         dat%myMesh%num_nodes = NP_G
         dat%myMesh%num_elems = NE_G
         IF (C3D.eqv..true.) THEN
            dat%myMesh%num_v_nodes = NFEN
         ENDIF
         dat%myMesh%nface_len = 3
      ENDIF
      ALLOCATE(dat%myTime%timenc(dat%myTime%timenc_len))
      ALLOCATE(defaultValue(dat%myMesh%num_nodes))
      defaultValue(:) = descript1%initial_value
      dat%myTime%initialized = .true.
C
C     Initialize netCDF output file, creating a new one if necessary.
      CALL createNetCDFOutputFile(dat%ncid, dat%myFile, dat%myTime,
     &                            descript1, reterr)
      !
      ! if we didn't need to create a file, update metadata and return
      IF (dat%myFile%createFile.eqv..false.) THEN
         IF ( (reterr.eqv..false.) .and.
     &        ( (abs(descript1%specifier).eq.NETCDF3) .or.
     &          (abs(descript1%specifier).eq.NETCDF4) ) ) then
            !
            ! jgf52.08.25: This subroutine is called in the case where
            ! a non-NetCDF min/max file will be written, but a NetCDF min/max
            ! file will be read. If this is the case, then don't update metadata
            ! in the netcdf file.
            !
Cobell...   The writers don't read the fort.15, so they need to
C           assume ADCPREP has placed the correct info on the first
C           pass, then call updateMetaData.
            IF ((MNWPROC.GT.0).AND.
     &           (writerReadMetaData.eqv..false.)) THEN
                CALL ReadMetaData(dat%ncid,dat%myFile)
                writerReadMetaData=.TRUE.
                dat%myMesh%num_nodes = NP_G
                dat%myMesh%num_elems = NE_G
            ENDIF
            ! jgf52.08.22: Only update meta data in serial, in parallel
            ! this was updated by adcprep
            if (mnproc.eq.1) then
               CALL updateMetaData(dat%ncid,dat%myFile)
            endif
         ENDIF
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         RETURN
      ENDIF
C
      IF (dat%myMesh%initialized.eqv..false.) THEN
         CALL initNetCDFCoord(dat%myMesh)
      ENDIF
      CALL defineMeshVariables(dat%ncid, dat%myMesh, dat%myFile)

      !WJP 02.20.2018 Writing in constituent information for harmonic file
      IF (any(nodalLunsHA.eq.descript1 % lun)) THEN
         CALL defineHarmonicAnalysisParametersinNetcdfFile(descript1,nodaldat=dat)
      ENDIF
c
      SELECT CASE(descript1 % lun)
C
         CASE(44)
            dat%nodal_data_dims_3D(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims_3D(2) = dat%myMesh%num_v_nodes_dim_id
            dat%nodal_data_dims_3D(3) = dat%myTime%timenc_dim_id
            iret = nf90_def_var(dat%ncid, 'sigmat', NF90_DOUBLE,
     &             dat%nodal_data_dims_3D, dat%u_nodal_data_id)
            CALL check_err(iret)
            IF ((IDEN.eq.2).or.(IDEN.eq.4)) THEN
               iret = nf90_def_var(dat%ncid, 'salinity', NF90_DOUBLE,
     &             dat%nodal_data_dims_3D, dat%v_nodal_data_id)
               CALL check_err(iret)
            ENDIF
            IF ((IDEN.eq.3).or.(IDEN.eq.4)) THEN
               iret = nf90_def_var(dat%ncid, 'temperature', NF90_DOUBLE,
     &             dat%nodal_data_dims_3D, dat%w_nodal_data_id)
               CALL check_err(iret)
            ENDIF
            ! sigma t
            att_text = "water column vertically varying density"
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &            'long_name', trim(att_text))
            CALL check_err(iret)
            att_text = "water_density_vertically_varying"
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &             'standard_name', trim(att_text))
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%u_nodal_data_id,
     &           "kg/m^3", "n/a")
            ! salinity
            IF ((IDEN.eq.2).or.(IDEN.eq.4)) THEN
               att_text = "water column vertically varying salinity"
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &            'long_name',trim(att_text))
               CALL check_err(iret)
               att_text = "water_salinity_vertically_varying"
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &             'standard_name', trim(att_text))
               CALL check_err(iret)
               CALL putUnitsAttribute(dat%ncid, dat%v_nodal_data_id,
     &           "PSU", "PSU")
            ENDIF
            ! temperature
            IF ((IDEN.eq.3).or.(IDEN.eq.4)) THEN
               att_text = "water column vertically varying temperature"
               iret = nf90_put_att(dat%ncid, dat%w_nodal_data_id,
     &            'long_name',trim(att_text))
               CALL check_err(iret)
               att_text = "water_temperature_vertically_varying"
               iret = nf90_put_att(dat%ncid, dat%w_nodal_data_id,
     &             'standard_name',trim(att_text))
               CALL check_err(iret)
               CALL putUnitsAttribute(dat%ncid, dat%w_nodal_data_id,
     &           "Celsius", "Fahrenheit")
            ENDIF
         CASE(45)
            dat%nodal_data_dims_3D(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims_3D(2) = dat%myMesh%num_v_nodes_dim_id
            dat%nodal_data_dims_3D(3) = dat%myTime%timenc_dim_id
            iret = nf90_def_var(dat%ncid, 'u-vel3D', NF90_DOUBLE,
     &             dat%nodal_data_dims_3D, dat%u_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(dat%ncid, 'v-vel3D', NF90_DOUBLE,
     &             dat%nodal_data_dims_3D, dat%v_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(dat%ncid, 'w-vel3D', NF90_DOUBLE,
     &             dat%nodal_data_dims_3D, dat%w_nodal_data_id)
            CALL check_err(iret)
            if (ics.ne.1) then
               ! u
               att_text =
     &         'water column vertically varying east/west velocity'
               iret = nf90_put_att(dat%ncid,
     &         dat%u_nodal_data_id, 'long_name', trim(att_text))
               CALL check_err(iret)
               att_text = 'eastward_water_velocity_vertically_varying'
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &             'standard_name', trim(att_text))
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid,
     &             dat%u_nodal_data_id, 'positive','east')
               CALL check_err(iret)
               ! v
               att_text =
     &         'water column vertically varying north/south velocity'
               iret = nf90_put_att(dat%ncid,
     &            dat%v_nodal_data_id, 'long_name',trim(att_text))
               CALL check_err(iret)
               att_text = 'northward_water_velocity_vertically_varying'
               iret = nf90_put_att(dat%ncid,
     &            dat%v_nodal_data_id, 'standard_name',trim(att_text))
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid,
     &            dat%v_nodal_data_id, 'positive', 'north')
               CALL check_err(iret)
               ! w
               att_text =
     &         'water column vertically varying up/down velocity'
               iret = nf90_put_att(dat%ncid, dat%w_nodal_data_id,
     &            'long_name',trim(att_text))
               CALL check_err(iret)
               att_text =
     &            'upward_water_velocity_vertically_varying'
               iret = nf90_put_att(dat%ncid, dat%w_nodal_data_id,
     &           'standard_name', trim(att_text))
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid,
     &           dat%w_nodal_data_id, 'positive', 'up')
               CALL check_err(iret)
            else
               ! u
               att_text
     & = 'water column vertically varying velocity in x-direction'
               iret = nf90_put_att(dat%ncid,
     &            dat%u_nodal_data_id,'long_name', trim(att_text))
               CALL check_err(iret)
               att_text = 'x_water_velocity_vertically_varying'
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &             'standard_name', trim(att_text))
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid,
     &            dat%u_nodal_data_id, 'positive', 'right')
               CALL check_err(iret)
               ! v
               att_text =
     & 'water column vertically varying velocity in y-direction'
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &           'long_name', trim(att_text))
               CALL check_err(iret)
               att_text = 'station_y_water_velocity_vertically_varying'
               iret = nf90_put_att(dat%ncid,
     &            dat%v_nodal_data_id, 'standard_name', trim(att_text))
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &           'positive', '90 degrees counterclockwise from x'//
     &            'water velocity')
               CALL check_err(iret)
               ! w
               att_text = 'water column vertically '
     &          //'varying velocity in z-direction'
               iret = nf90_put_att(dat%ncid,
     &            dat%w_nodal_data_id, 'long_name', trim(att_text))
               CALL check_err(iret)
               att_text = 'station_z_water_velocity_vertically_varying'
               iret = nf90_put_att(dat%ncid,
     &            dat%w_nodal_data_id, 'standard_name', trim(att_text))
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid,
     &            dat%w_nodal_data_id, 'positive',
     &            '90 degrees counterclockwise from x water velocity')
               CALL check_err(iret)
            endif
            CALL putUnitsAttribute(dat%ncid,dat%u_nodal_data_id,
     &         'm s-1', 'ft s-1')
            CALL putUnitsAttribute(dat%ncid,dat%v_nodal_data_id,
     &         'm s-1', 'ft s-1')
            CALL putUnitsAttribute(dat%ncid,dat%w_nodal_data_id,
     &         'm s-1', 'ft s-1')

         CASE(46) !       F O R T . 4 6
            dat%nodal_data_dims_3D(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims_3D(2) = dat%myMesh%num_v_nodes_dim_id
            dat%nodal_data_dims_3D(3) = dat%myTime%timenc_dim_id
            iret = nf90_def_var(dat%ncid, 'q20', NF90_DOUBLE,
     &          dat%nodal_data_dims_3D, dat%u_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(dat%ncid, 'l', NF90_DOUBLE,
     &          dat%nodal_data_dims_3D, dat%v_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(dat%ncid, 'ev', NF90_DOUBLE,
     &          dat%nodal_data_dims_3D, dat%w_nodal_data_id)
            CALL check_err(iret)
            ! q20
            att_text =
     &  'water column vertically varying turbulent kinetic energy'
            iret = nf90_put_att(dat%ncid,
     &         dat%u_nodal_data_id, 'long_name', trim(att_text))
            CALL check_err(iret)
            att_text =
     &         'turbulent_kinetic_energy_vertically_varying'
            iret = nf90_put_att(dat%ncid,
     &             dat%u_nodal_data_id,'standard_name',trim(att_text))
            CALL check_err(iret)
            ! l
            att_text =
     &         'water column vertically varying mixing length'
            iret = nf90_put_att(dat%ncid,
     &            dat%v_nodal_data_id, 'long_name',trim(att_text))
            CALL check_err(iret)
            att_text =
     &            'water_mixing_length_vertically_varying'
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &           'standard_name', trim(att_text))
               CALL check_err(iret)
            !  ev
            att_text =
     &        'water column vertically varying eddy viscosity'
               iret = nf90_put_att(dat%ncid,
     &            dat%w_nodal_data_id, 'long_name', trim(att_text))
               CALL check_err(iret)
               att_text = 'water_eddy_viscosity_vertically_varying'
               iret = nf90_put_att(dat%ncid,
     &            dat%w_nodal_data_id, 'standard_name', trim(att_text))
               CALL check_err(iret)
C            CALL putUnitsAttribute(dat%ncid,dat%u_nodal_data_id,
C     &      'm s-2', 'ft s-1') !TODO: jgf49.48.01: units for q20??
            CALL putUnitsAttribute(dat%ncid,dat%v_nodal_data_id,
     &      "meters","n/a")
            CALL putUnitsAttribute(dat%ncid,dat%w_nodal_data_id,
     &      "m^2/s","n/a")
C
         CASE(47)
C
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            dat%varnames(1) = 'qsurfkp1'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &         NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
            CALL check_err(iret)
            att_text = "sea surface temperature at the k+1 time level"
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'long_name', trim(att_text))
            CALL check_err(iret)
            att_text = "future sea surface temperature"
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'standard_name', trim(att_text))
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%nodal_data_id,
     &           "Celsius", "Fahrenheit")
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
C        WJP 02.20.2018: added fort.53-54 harmonic outputs
         CASE(53)
C
            dat%nodal_data_dims(2) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(1) = dat%myMesh%num_v_nodes_dim_id
!           Amplitude
            dat%varnames(1) = 'amp'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &         NF90_DOUBLE, dat%nodal_data_dims, dat%ha_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%ha_nodal_data_id,
     &            'long_name',
     &      'amplitude of tidal harmonic constituent elevation')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%ha_nodal_data_id,
     &             'standard_name', 'hc_elev_amp')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%ha_nodal_data_id,
     &             'coordinates', 'const y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%ha_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%ha_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%ha_nodal_data_id,
     &           'm', 'ft')
            iret = nf90_put_att(dat%ncid, dat%ha_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
!           Phase
            dat%varnames(1) = 'phs'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &         NF90_DOUBLE, dat%nodal_data_dims, dat%hg_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%hg_nodal_data_id,
     &            'long_name',
     &      'phase of tidal harmonic constituent elevation')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%hg_nodal_data_id,
     &             'standard_name', 'hc_elev_phs')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%hg_nodal_data_id,
     &             'coordinates', 'const y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%hg_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%hg_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%hg_nodal_data_id,
     &           'deg', 'deg')
            iret = nf90_put_att(dat%ncid, dat%hg_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
C
         CASE(54)
C
            dat%nodal_data_dims(2) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(1) = dat%myMesh%num_v_nodes_dim_id
!           u Amplitude
            dat%varnames(1) = 'u_amp'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &         NF90_DOUBLE, dat%nodal_data_dims, dat%u_ha_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_ha_nodal_data_id,
     &            'long_name',
     &      'amplitude of tidal harmonic constituent E-W velocity')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_ha_nodal_data_id,
     &             'standard_name', 'hc_u_amp')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_ha_nodal_data_id,
     &             'coordinates', 'const y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_ha_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_ha_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%u_ha_nodal_data_id,
     &           'm/s', 'ft/s')
            iret = nf90_put_att(dat%ncid, dat%u_ha_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
!           u Phase
            dat%varnames(1) = 'u_phs'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &         NF90_DOUBLE, dat%nodal_data_dims, dat%u_hg_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_hg_nodal_data_id,
     &            'long_name',
     &      'phase of tidal harmonic constituent E-W velocity')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_hg_nodal_data_id,
     &             'standard_name', 'hc_u_phs')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_hg_nodal_data_id,
     &             'coordinates', 'const y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_hg_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_hg_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%u_hg_nodal_data_id,
     &           'deg', 'deg')
            iret = nf90_put_att(dat%ncid, dat%u_hg_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
!           v amplitude
            dat%varnames(1) = 'v_amp'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &         NF90_DOUBLE, dat%nodal_data_dims, dat%v_ha_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_ha_nodal_data_id,
     &            'long_name',
     &      'amplitude of tidal harmonic constituent N-S velocity')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_ha_nodal_data_id,
     &             'standard_name', 'hc_v_amp')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_ha_nodal_data_id,
     &             'coordinates', 'const y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_ha_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_ha_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%v_ha_nodal_data_id,
     &           'm/s', 'ft/s')
            iret = nf90_put_att(dat%ncid, dat%v_ha_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
!           v Phase
            dat%varnames(1) = 'v_phs'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &         NF90_DOUBLE, dat%nodal_data_dims, dat%v_hg_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_hg_nodal_data_id,
     &            'long_name',
     &      'phase of tidal harmonic constituent N-S velocity')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_hg_nodal_data_id,
     &             'standard_name', 'hc_v_phs')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_hg_nodal_data_id,
     &             'coordinates', 'const y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_hg_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_hg_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%v_hg_nodal_data_id,
     &           'deg', 'deg')
            iret = nf90_put_att(dat%ncid, dat%v_hg_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
C
         CASE(63)
C
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            dat%varnames(1) = 'zeta'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &         NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'long_name', 'water surface elevation above geoid')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'standard_name', 'sea_surface_height_above_geoid')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%nodal_data_id,
     &           'm', 'ft')
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
C
         CASE(90)
C
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            dat%varnames(1) = 'tau0'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &         NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'long_name', 'primitive weighting in continuity equation')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'standard_name', 'primitive_weighting_in_continuity_equation')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%nodal_data_id,
     &           '1', '1')
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
C
         CASE(311)
C
            dat%varnames(1) = 'zeta_max'
            dat%varnames(2) = 'time_of_zeta_max'
            dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &      NF90_DOUBLE, dat%nodal_data_dims_max, dat%max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(2)),
     &      NF90_DOUBLE, dat%nodal_data_dims_max, dat%time_max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'long_name', 'maximum water surface elevation'//
     &            'above geoid')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &            'long_name', 'time of maximum water surface elevation'//
     &            'above geoid')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'standard_name', 'maximum_sea_surface_height_'//
     &             'above_geoid')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'standard_name', 'time_of_maximum_sea_surface_height_'//
     &             'above_geoid')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%max_nodal_data_id,
     &           'm', 'ft')
            CALL putUnitsAttribute(dat%ncid, dat%time_max_nodal_data_id,
     &           'sec', 'sec')
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)

C
         CASE(64)
C
            dat%varnames(1) = 'u-vel'
            dat%varnames(2) = 'v-vel'
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &          NF90_DOUBLE, dat%nodal_data_dims, dat%u_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(2)),
     &         NF90_DOUBLE, dat%nodal_data_dims, dat%v_nodal_data_id)
            CALL check_err(iret)
            if (ics.ne.1) then
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                'long_name',
     &   'water column vertically averaged east/west velocity')
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                'standard_name', 'eastward_water_velocity')
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                'positive', 'east')
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                'long_name',
     &    'water column vertically averaged north/south velocity')
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &               'standard_name','northward_water_velocity')
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                'positive', 'north')
               CALL check_err(iret)
            else
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                'long_name',
     &  'water column vertically averaged velocity in x-direction')
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &               'standard_name',
     &               'x_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                'positive', 'right')
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                'long_name',
     &       'water column vertically averaged velocity in y-direction')
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &               'standard_name',
     &               'y_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &             'positive',
     &             '90 degrees counterclockwise from x water velocity')
               CALL check_err(iret)
            endif
            CALL putUnitsAttribute(dat%ncid, dat%u_nodal_data_id,
     &           'm s-1', 'ft s-1')
            CALL putUnitsAttribute(dat%ncid, dat%v_nodal_data_id,
     &           'm s-1', 'ft s-1')
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &           '_FillValue',  doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &           '_FillValue',  doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &           'dry_Value',  doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &            'dry_Value',  doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &           'coordinates',  'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &           'coordinates',  'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &           'location',  'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &           'location',  'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &           'mesh',  'adcirc_mesh')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &           'mesh',  'adcirc_mesh')
            CALL check_err(iret)

         CASE(312)
            dat%varnames(1) = 'vel_max'
            dat%varnames(2) = 'time_of_vel_max'
            dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &             NF90_DOUBLE, dat%nodal_data_dims_max,
     &             dat%max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(2)),
     &             NF90_DOUBLE, dat%nodal_data_dims_max,
     &             dat%time_max_nodal_data_id)
            CALL check_err(iret)

            if (ics.ne.1) then
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &                'long_name',
     &         'maximum water column vertically averaged velocity')
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &                'long_name',
     &       'time of maximum water column vertically averaged velocity')
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &         'standard_name', 'maximum_water_velocity')
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &         'standard_name', 'time_of_maximum_water_velocity')
               CALL check_err(iret)
            else
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &                'long_name',
     &                'maximum water column vertically averaged'//
     &                'velocity')
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &                'long_name',
     &                'time of maximum water column vertically averaged'//
     &                'velocity')
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &               'standard_name',
     &               'maximum_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &               'standard_name',
     &               'time_of_maximum_water_velocity_depth_averaged')
               CALL check_err(iret)
            endif
            CALL putUnitsAttribute(dat%ncid, dat%max_nodal_data_id,
     &           'm s-1', 'ft s-1')
            CALL putUnitsAttribute(dat%ncid, dat%time_max_nodal_data_id,
     &           'sec', 'sec')
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &           '_FillValue',  doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &           '_FillValue',  doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &           'dry_Value',  doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &           'dry_Value',  doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &           'coordinates',  'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &           'coordinates',  'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &           'location',  'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &           'location',  'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &           'mesh',  'adcirc_mesh')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &           'mesh',  'adcirc_mesh')
            CALL check_err(iret)
C
         CASE(173)
C
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            dat%varnames(1) = 'winddrag'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &         NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'long_name', 'wind drag coefficient at sea level')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'standard_name', 'wind drag coefficient')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%nodal_data_id,
     &           'unitless', 'unitless')
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
C
         CASE(73)
C
            dat%varnames(1) = 'pressure'
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            iret = nf90_def_var(dat%ncid,
     &              trim(dat%varnames(1)),NF90_DOUBLE,
     &              dat%nodal_data_dims,
     &              dat%nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'long_name','air pressure at sea level')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'standard_name', 'air_pressure_at_sea_level')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id, 'units',
     &             'meters of water')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'mesh', 'adcirc_mesh')
            CALL check_err(iret)
C
         CASE(313)
C
            dat%varnames(1) = 'pressure_min'
            dat%varnames(2) = 'time_of_pressure_min'
            dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
            iret = nf90_def_var(dat%ncid,
     &              trim(dat%varnames(1)),NF90_DOUBLE,
     &              dat%nodal_data_dims_max,
     &              dat%max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(dat%ncid,
     &              trim(dat%varnames(2)),NF90_DOUBLE,
     &              dat%nodal_data_dims_max,
     &              dat%time_max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'long_name','minimum air pressure at sea level')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'long_name','time of minimum air pressure at sea level')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'standard_name', 'minimum_air_pressure_at_sea_level')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'standard_name',
     &         'time_of_minimum_air_pressure_at_sea_level')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,'units',
     &             'meters of water')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,'units',
     &             'seconds')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &            '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &            'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &            'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &            'mesh', 'adcirc_mesh')
            CALL check_err(iret)
C
         CASE(74)
C
            dat%varnames(1) = 'windx'
            dat%varnames(2) = 'windy'
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &      NF90_DOUBLE, dat%nodal_data_dims,dat%u_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(2)),
     &             NF90_DOUBLE, dat%nodal_data_dims,dat%v_nodal_data_id)
            CALL check_err(iret)
            if (ics.ne.1) then
               select case(abs(nws))
                  case(1,2)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                      'long_name',
     &                       'e/w wind stress')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                      'standard_name',
     &                      'eastward_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                      'positive', 'east')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                       'long_name',
     &                       'n/s wind stress')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                      'standard_name',
     &                      'northward_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                      'positive', 'north')
                     CALL check_err(iret)
                  case default
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &               'long_name', 'e/w wind velocity')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                      'standard_name', 'eastward_wind')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                     'positive', 'east')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                     'long_name', 'n/s wind velocity')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                       'standard_name', 'northward_wind')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                      'positive', 'north')
               end select
            else
               select case(abs(nws))
                  case(1,2)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                       'long_name','wind stress in x-direction')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                      'standard_name','x_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                      'positive', 'right')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                       'long_name','wind stress in y-direction')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                       'standard_name','y_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                     'positive',
     & '90 degrees counterclockwise from wind velocity in x-direction')
                     CALL check_err(iret)
                  case default
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                      'long_name','wind velocity in x-direction')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                      'standard_name','x_wind')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                     'positive', 'right')

                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                       'long_name','wind velocity in y-direction')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                      'standard_name', 'y_wind')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                      'positive',
     & '90 degrees counterclockwise from wind velocity in x-direction')
                     CALL check_err(iret)
               end select
            endif
            select case(abs(nws))
               case(1,2)
               call putUnitsAttribute(dat%ncid, dat%u_nodal_data_id,
     &              'm2 s-2', 'lbf s-2')
               call putUnitsAttribute(dat%ncid, dat%v_nodal_data_id,
     &              'm2 s-2', 'lbf s-2')
               case default
               call putUnitsAttribute(dat%ncid, dat%u_nodal_data_id,
     &              'm s-1', 'ft s-1')
               call putUnitsAttribute(dat%ncid, dat%v_nodal_data_id,
     &              'm s-1', 'ft s-1')
            end select
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &             'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &             'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
C
         CASE(314)
C
            dat%varnames(1) = 'wind_max'
            dat%varnames(2) = 'time_of_wind_max'
            dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &             NF90_DOUBLE, dat%nodal_data_dims_max,
     &             dat%max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(2)),
     &             NF90_DOUBLE, dat%nodal_data_dims_max,
     &             dat%time_max_nodal_data_id)
            CALL check_err(iret)
            if (ics.ne.1) then
               select case(abs(nws))
                  case(1,2)
                     iret=nf90_put_att(dat%ncid,dat%max_nodal_data_id,
     &                      'long_name',
     &                       'maximum wind stress')
                     CALL check_err(iret)
                     iret=nf90_put_att(dat%ncid,dat%time_max_nodal_data_id,
     &                      'long_name',
     &                       'time of maximum wind stress')
                     CALL check_err(iret)
                     iret=nf90_put_att(dat%ncid,dat%max_nodal_data_id,
     &                      'standard_name',
     &                      'maximum_surface_wind_stress')
                     CALL check_err(iret)
                     iret=nf90_put_att(dat%ncid,dat%time_max_nodal_data_id,
     &                      'standard_name',
     &                      'time_of_maximum_surface_wind_stress')
                     CALL check_err(iret)
                  case default
                     iret=nf90_put_att(dat%ncid,dat%max_nodal_data_id,
     &               'long_name', 'maximum wind velocity')
                     CALL check_err(iret)
                     iret=nf90_put_att(dat%ncid,dat%time_max_nodal_data_id,
     &               'long_name', 'time of maximum wind velocity')
                     CALL check_err(iret)
                     iret=nf90_put_att(dat%ncid,dat%max_nodal_data_id,
     &                      'standard_name', 'maximum_wind')
                     CALL check_err(iret)
                     iret=nf90_put_att(dat%ncid,dat%time_max_nodal_data_id,
     &                      'standard_name', 'time_of_maximum_wind')
                     CALL check_err(iret)
               end select
            else
               select case(abs(nws))
                  case(1,2)
                     iret=nf90_put_att(dat%ncid,dat%max_nodal_data_id,
     &                    'long_name','maximum wind stress')
                     CALL check_err(iret)
                     iret=nf90_put_att(dat%ncid,dat%time_max_nodal_data_id,
     &                    'long_name','time of maximum wind stress')
                     CALL check_err(iret)
                     iret=nf90_put_att(dat%ncid,dat%max_nodal_data_id,
     &                    'standard_name','maximum_surface_wind_stress')
                     CALL check_err(iret)
                     iret=nf90_put_att(dat%ncid,dat%time_max_nodal_data_id,
     &                    'standard_name',
     &                    'time_of_maximum_surface_wind_stress')
                     CALL check_err(iret)
                  case default
                  iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &                   'long_name','maximum wind velocity')
                  CALL check_err(iret)
                  iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &                   'long_name','time of maximum wind velocity')
                  CALL check_err(iret)
                  iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &                   'standard_name','maximum_wind')
                  CALL check_err(iret)
                  iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &                   'standard_name','time_of_maximum_wind')
                  CALL check_err(iret)
               end select
            endif
            select case(abs(nws))
               case(1,2)
               call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id,
     &              'm2 s-2', 'lbf ft-2')
               call putUnitsAttribute(dat%ncid, dat%time_max_nodal_data_id,
     &              'sec', 'sec')
               case default
               call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id,
     &              'm s-1', 'ft s-1')
               call putUnitsAttribute(dat%ncid, dat%time_max_nodal_data_id,
     &              'sec', 'sec')
            end select
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
C
         CASE(315)
C
            dat%varnames(1) = 'radstress_max'
            dat%varnames(2) = 'time_of_radstress_max'
            dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &            NF90_DOUBLE,dat%nodal_data_dims_max,
     &            dat%max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(2)),
     &            NF90_DOUBLE,dat%nodal_data_dims_max,
     &            dat%time_max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'long_name', 'maximum radiation stress gradient')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'long_name',
     &             'time of maximum radiation stress gradient')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'standard_name', 'maximum_radiation_stress')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'standard_name', 'time_of_maximum_radiation_stress')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%max_nodal_data_id,
     &           'm-2 s-2', 'lbf ft-2')
            CALL putUnitsAttribute(dat%ncid, dat%time_max_nodal_data_id,
     &           'sec', 'sec')
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)

        CASE(77)
            dat%varnames(1) = 'weir_dz'
            dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &        NF90_DOUBLE,dat%nodal_data_dims_max,dat%nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &           'long_name',
     &           'elevation change in weir boundary condition')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &           'standard_name','weir_elevation_change')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%nodal_data_id,
     &           'm', 'ft')
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &           'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             '_FillValue', 0.0d0)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)

C++ GML 20200528
         CASE(93)
            dat%varnames(1) = 'iceaf' ;
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id

            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &             NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
            CALL check_err(iret)

            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &           'long_name', 'sea-ice af')
            CALL check_err(iret)

            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &           'standard_name', 'sea-ice area fraction')
            CALL check_err(iret)

            !
            call putUnitsAttribute(dat%ncid, dat%nodal_data_id,
     &                 'unitless', 'unitless')
      
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)

            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'coordinates', 'time y x')
            CALL check_err(iret)

            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)


            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
C++

         CASE(164)
C
            dat%varnames(1) = 'radstress_x'
            dat%varnames(2) = 'radstress_y'
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &      NF90_DOUBLE, dat%nodal_data_dims, dat%u_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(2)),
     &             NF90_DOUBLE, dat%nodal_data_dims,dat%v_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &            'long_name', 'radiation stress gradient x component')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &            'long_name', 'radiation stress gradient y component')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &             'standard_name', 'radiation_stress_gradient_x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &             'standard_name', 'radiation_stress_gradient_y')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%u_nodal_data_id,
     &           'm-2 s-2', 'ft-2 s-2')
            CALL putUnitsAttribute(dat%ncid, dat%v_nodal_data_id,
     &           'm-2 s-2', 'ft-2 s-2')
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &             'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &             'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
C

#if defined CSWAN || defined ADCSWAN
C
C
         CASE(301)
C
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            dat%varnames(1) = 'swan_HS'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &         NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'long_name', 'significant wave height')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'standard_name',
     &             'sea_surface_wave_significant_height')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%nodal_data_id,
     &           'm', 'ft')
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'mesh', 'adcirc_mesh')
            CALL check_err(iret)
C
         CASE(316)
C
            dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
            dat%varnames(1) = 'swan_HS_max'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &      NF90_DOUBLE, dat%nodal_data_dims_max, dat%max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'long_name', 'maximum significant wave height')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'standard_name',
     &             'maximum_sea_surface_wave_significant_height')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%max_nodal_data_id,
     &           'm', 'ft')
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'mesh', 'adcirc_mesh')
            CALL check_err(iret)
C
         CASE(302)
C
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            dat%varnames(1) = 'swan_DIR'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &         NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'long_name', 'mean wave direction')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'standard_name', 'sea_surface_wave_to_direction')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%nodal_data_id,
     &           'degrees', 'degrees_CW_from_East')
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'mesh', 'adcirc_mesh')
            CALL check_err(iret)
C
         CASE(317)
C
            dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
            dat%varnames(1) = 'swan_DIR_max'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &             NF90_DOUBLE, dat%nodal_data_dims_max,
     &             dat%max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'long_name', 'maximum mean wave direction')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &      'standard_name', 'maximum_sea_surface_wave_to_direction')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%max_nodal_data_id,
     &           'degrees', 'degrees_CW_from_East')
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'mesh', 'adcirc_mesh')
            CALL check_err(iret)
C
         CASE(303)
C
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            dat%varnames(1) = 'swan_TM01'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &         NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'long_name', 'mean absolute wave period')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'standard_name',
     &      'sea_surface_wave_mean_period_from_variance_spectral'//
     &      '_density_first_frequency_moment')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%nodal_data_id,
     &           's', 's')
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'mesh', 'adcirc_mesh')
            CALL check_err(iret)
C
         CASE(318)
C
            dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
            dat%varnames(1) = 'swan_TM01_max'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &             NF90_DOUBLE, dat%nodal_data_dims_max,
     &             dat%max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'long_name', 'maximum TM01 mean wave period')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'standard_name',
     &      'maximum_sea_surface_wave_mean_period_from_variance_'//
     &      'spectral_density_first_frequency_moment')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%max_nodal_data_id,
     &           's', 's')
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'mesh', 'adcirc_mesh')
            CALL check_err(iret)
C
         CASE(304)
C
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            dat%varnames(1) = 'swan_TPS'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &         NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'long_name', 'smoothed peak period')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'standard_name','sea_surface_wave_period_'//
     &             'at_variance_spectral_density_maximum')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%nodal_data_id,
     &           's', 's')
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'mesh', 'adcirc_mesh')
            CALL check_err(iret)
C
         CASE(319)
C
            dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
            dat%varnames(1) = 'swan_TPS_max'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &      NF90_DOUBLE, dat%nodal_data_dims_max, dat%max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'long_name', 'maximum smoothed peak period')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'standard_name','maximum_sea_surface_wave_period'//
     &             '_at_variance_spectral_density_maximum')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%max_nodal_data_id,
     &           's', 's')
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'mesh', 'adcirc_mesh')
            CALL check_err(iret)
C
         CASE(305)
C
            dat%varnames(1) = 'swan_windx'
            dat%varnames(2) = 'swan_windy'
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &      NF90_DOUBLE,dat%nodal_data_dims,dat%u_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(2)),
     &             NF90_DOUBLE, dat%nodal_data_dims,dat%v_nodal_data_id)
            CALL check_err(iret)
            if (ics.ne.1) then
               select case(abs(nws))
                  case(1,2)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                      'long_name',
     &                       'e/w wind stress')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                      'standard_name',
     &                      'eastward_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                      'positive', 'east')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                       'long_name',
     &                       'n/s wind stress')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                      'standard_name',
     &                      'northward_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                      'positive', 'north')
                     CALL check_err(iret)
                  case default
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &               'long_name', 'e/w wind velocity')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                      'standard_name', 'eastward_wind')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                     'positive', 'east')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                     'long_name', 'n/s wind velocity')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                       'standard_name', 'northward_wind')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                      'positive', 'north')
               end select
            else
               select case(abs(nws))
                  case(1,2)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                       'long_name','wind stress in x-direction')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                      'standard_name','x_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                      'positive', 'right')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                       'long_name','wind stress in y-direction')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                       'standard_name','y_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                     'positive',
     & '90 degrees counterclockwise from wind velocity in x-direction')
                  case default
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                      'long_name','wind velocity in x-direction')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                      'standard_name','x_wind')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &                     'positive', 'right')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                       'long_name','wind velocity in y-direction')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                      'standard_name', 'y_wind')
                     CALL check_err(iret)
                     iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &                      'positive',
     & '90 degrees counterclockwise from wind velocity in x-direction')
                     CALL check_err(iret)
               end select
            endif

            select case(abs(nws))
               case(1,2)
               call putUnitsAttribute(dat%ncid, dat%u_nodal_data_id,
     &                 'm2 s-2', 'ft s-2')
               call putUnitsAttribute(dat%ncid, dat%v_nodal_data_id,
     &                 'm2 s-2', 'ft s-2')
               case default
               call putUnitsAttribute(dat%ncid, dat%u_nodal_data_id,
     &                 'm s-1', 'ft s-1')
               call putUnitsAttribute(dat%ncid, dat%v_nodal_data_id,
     &                 'm s-1', 'ft s-1')
            end select
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &             'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &             'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
C
         CASE(320)
C
            dat%varnames(1) = 'swan_wind_max'
            dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &      NF90_DOUBLE,dat%nodal_data_dims_max,dat%max_nodal_data_id)
            CALL check_err(iret)
            if (ics.ne.1) then
               select case(abs(nws))
                  case(1,2)
                  iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &                  'long_name',
     &                  'maximum wind stress')
                  CALL check_err(iret)
                  iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &                   'standard_name',
     &                   'maximum_surface_wind_stress')
                  CALL check_err(iret)
                  case default
                  iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &           'long_name', 'maximum wind velocity')
                  CALL check_err(iret)
                  iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &                   'standard_name', 'maximum_wind')
                  CALL check_err(iret)
               end select
            else
               select case(abs(nws))
                  case(1,2)
                  iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'long_name','maximum wind stress')
                  CALL check_err(iret)
                  iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'standard_name','maximum_surface_wind_stress')
                  CALL check_err(iret)
                  case default
                  iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'long_name','maximum wind velocity')
                  CALL check_err(iret)
                  iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'standard_name','maximum_wind')
                  CALL check_err(iret)
               end select
            endif
            select case(abs(nws))
               case(1,2)
               call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id,
     &                 'm2 s-2', 'ft s-2')
               case default
               call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id,
     &                 'm s-1', 'ft s-1')
            end select
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
C
         CASE(306)
C
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            dat%varnames(1) = 'swan_TM02'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &         NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'long_name', 'mean absoloute zero crossing period')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'standard_name',
     &      'sea_surface_wave_mean_period_from_variance_spectral'//
     &      '_density_second_frequency_moment')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%nodal_data_id,
     &           's', 's')
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'mesh', 'adcirc_mesh')
            CALL check_err(iret)
C
         CASE(321)
C
            dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
            dat%varnames(1) = 'swan_TM02_max'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &      NF90_DOUBLE, dat%nodal_data_dims_max, dat%max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'long_name', 'maximum TM02 mean wave period')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'standard_name',
     &      'maximum_sea_surface_wave_mean_period_from_variance_'//
     &      'spectral_density_second_frequency_moment')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%max_nodal_data_id,
     &           's', 's')
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'mesh', 'adcirc_mesh')
            CALL check_err(iret)
C
         CASE(307)
C
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            dat%varnames(1) = 'swan_TMM10'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &         NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'long_name', 'mean absolute wave period')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'standard_name',
     &      'sea_surface_wave_mean_period_from_variance_spectral'//
     &      '_density_inverse_frequency_moment')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%nodal_data_id,
     &           's', 's')
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'mesh', 'adcirc_mesh')
            CALL check_err(iret)
C
         CASE(322)
C
            dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
            dat%varnames(1) = 'swan_TMM10_max'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &      NF90_DOUBLE, dat%nodal_data_dims_max, dat%max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'long_name', 'maximum TMM10 mean wave period')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'standard_name',
     &      'maximum_sea_surface_wave_mean_period_from_variance_'//
     &      'spectral_density_inverse_frequency_moment')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%max_nodal_data_id,
     &           's', 's')
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'mesh', 'adcirc_mesh')
            CALL check_err(iret)
#endif
         ! inundationtime.63
         CASE(400)
C
            dat%varnames(1) = 'inun_time'
            dat%varnames(2) = 'onset_inun_time'
            dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &      NF90_DOUBLE, dat%nodal_data_dims_max, dat%max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(2)),
     &      NF90_DOUBLE, dat%nodal_data_dims_max, dat%time_max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &      'long_name', 'total time of inundation beyond the threshold')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &            'long_name',
     &            'time of onset of inundation beyond the threshold')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'standard_name',
     &             'total_time_of_inundation_beyond_the_threshold')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'standard_name',
     &             'time_of_onset_of_inundation_beyond_the_threshold')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%max_nodal_data_id,
     &           'm', 'ft')
            CALL putUnitsAttribute(dat%ncid, dat%time_max_nodal_data_id,
     &           'sec', 'sec')
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)

         ! maxinundepth.63
         CASE(401)
C
            dat%varnames(1) = 'inun_max'
            dat%varnames(2) = 'time_of_inun_max'
            dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &      NF90_DOUBLE, dat%nodal_data_dims_max, dat%max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(2)),
     &      NF90_DOUBLE, dat%nodal_data_dims_max, dat%time_max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'long_name', 'maximum water inundation depth'//
     &            'above geoid')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &            'long_name', 'time of maximum water inundation depth'//
     &            'above geoid')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'standard_name', 'maximum_water_inundation_depth')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'standard_name',
     &             'time_of_maximum_water_inundation_depth')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%max_nodal_data_id,
     &           'm', 'ft')
            CALL putUnitsAttribute(dat%ncid, dat%time_max_nodal_data_id,
     &           'sec', 'sec')
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)

        ! initiallydry.63
        CASE(402)
C
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            dat%varnames(1) = 'initiallydry'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &         NF90_INT, dat%nodal_data_dims, dat%nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'long_name',
     &            'dry nodes at cold start')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'standard_name', 'initially_dry')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)

            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             '_FillValue', intFillValue)
            CALL check_err(iret)

        ! endrisinginun.63
        CASE(403)
C
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            dat%varnames(1) = 'endrisinginun'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &         NF90_INT, dat%nodal_data_dims, dat%nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &            'long_name',
     &            'rising inundation at the end of the simulation')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'standard_name', 'end_rising_inundation')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             '_FillValue', intFillValue)
            CALL check_err(iret)
C
         ! everdried.63
         CASE(404)
C
            dat%varnames(1) = 'everdried'
            dat%varnames(2) = 'time_of_everdried'
            dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &      NF90_DOUBLE, dat%nodal_data_dims_max, dat%max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(2)),
     &      NF90_DOUBLE, dat%nodal_data_dims_max, dat%time_max_nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &            'long_name', 'ever dried')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &            'long_name',
     &            'time of most recent drying occurrence')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'standard_name', 'ever_dried')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'standard_name',
     &             'time_of_most_recent_drying_occurrence')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'coordinates', 'y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%max_nodal_data_id,
     &           'unitless', 'unitless')
            CALL putUnitsAttribute(dat%ncid, dat%time_max_nodal_data_id,
     &           'sec', 'sec')
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)

         !jgf: dynamicWaterlevelCorrection.63
         CASE(108)
C
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            dat%varnames(1) = 'dynamicWaterlevelCorrection'
            iret = nf90_def_var(dat%ncid,trim(dat%varnames(1)),
     &         NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &         'long_name',
     &         'dynamic water surface correction above water level')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &         'standard_name',
     &         'dynamic_sea_surface_correction_above_water_level')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'coordinates', 'time y x')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'location', 'node')
            CALL check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             'mesh', 'adcirc_mesh')
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%nodal_data_id,
     &           'm', 'ft')
            iret = nf90_put_att(dat%ncid, dat%nodal_data_id,
     &             '_FillValue', doubleval)
            CALL check_err(iret)

         CASE DEFAULT
            write(scratchMessage,
     &       '("No netCDF for files with unit number ",i0,".")') descript1 % lun
            call allMessage(ERROR,scratchMessage)
      END SELECT

      IF (IDEN.ne.44) THEN
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id,
     &      '_FillValue',  doubleval)
         CALL check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id,
     &      '_FillValue',  doubleval)
         CALL check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%w_nodal_data_id,
     &      '_FillValue',  doubleval)
         CALL check_err(iret)
      ENDIF


C
C     jgf50.44: Automatically turn on compression if we are using the
C     netcdf4 file format.
#ifdef NETCDF_CAN_DEFLATE
      IF (abs(descript1%specifier).eq.NETCDF4) THEN
         SELECT CASE(descript1 % lun)
         CASE(44)
            iret = nf_def_var_deflate(dat%ncid, dat%u_nodal_data_id,
     &          1, 1, 2)
            CALL check_err(iret)
            IF ((IDEN.eq.2).or.(IDEN.eq.4)) THEN
               iret = nf_def_var_deflate(dat%ncid, dat%v_nodal_data_id,
     &                1, 1, 2)
               CALL check_err(iret)
            ENDIF
            IF ((IDEN.eq.3).or.(IDEN.eq.4)) THEN
               iret = nf_def_var_deflate(dat%ncid, dat%w_nodal_data_id,
     &                1, 1, 2)
               CALL check_err(iret)
            ENDIF
         CASE(45,46)
            iret = nf_def_var_deflate(dat%ncid, dat%u_nodal_data_id,
     &          1, 1, 2)
            CALL check_err(iret)
            iret = nf_def_var_deflate(dat%ncid, dat%v_nodal_data_id,
     &             1, 1, 2)
            CALL check_err(iret)
            iret = nf_def_var_deflate(dat%ncid, dat%w_nodal_data_id,
     &             1, 1, 2)
            CALL check_err(iret)
         CASE(47,63,73,77,90,93,173,301:304,306,307,402,403,108) ! GML added 93 20210727
            iret = nf_def_var_deflate(dat%ncid, dat%nodal_data_id,
     &             1, 1, 2)
            CALL check_err(iret)
         !WJP 02.20.2018: added fort.53-54 harmonic outputs
         CASE(53)
            iret = nf_def_var_deflate(dat%ncid, dat%ha_nodal_data_id,
     &             1, 1, 2)
            CALL check_err(iret)
            iret = nf_def_var_deflate(dat%ncid, dat%hg_nodal_data_id,
     &             1, 1, 2)
            CALL check_err(iret)
         CASE(54)
            iret = nf_def_var_deflate(dat%ncid, dat%u_ha_nodal_data_id,
     &             1, 1, 2)
            CALL check_err(iret)
            iret = nf_def_var_deflate(dat%ncid, dat%u_hg_nodal_data_id,
     &             1, 1, 2)
            CALL check_err(iret)
            iret = nf_def_var_deflate(dat%ncid, dat%v_ha_nodal_data_id,
     &             1, 1, 2)
            CALL check_err(iret)
            iret = nf_def_var_deflate(dat%ncid, dat%v_hg_nodal_data_id,
     &             1, 1, 2)
            CALL check_err(iret)
         CASE(64,74,164,305)
            iret = nf_def_var_deflate(dat%ncid, dat%u_nodal_data_id,
     &             1, 1, 2)
            CALL check_err(iret)
            iret = nf_def_var_deflate(dat%ncid, dat%v_nodal_data_id,
     &             1, 1, 2)
            CALL check_err(iret)
         CASE(311:315,400,401,404)  !adcirc max/min values
            iret = nf_def_var_deflate(dat%ncid, dat%max_nodal_data_id,
     &             1, 1, 2)
            CALL check_err(iret)
            iret = nf_def_var_deflate(dat%ncid, dat%time_max_nodal_data_id,
     &             1, 1, 2)
            CALL check_err(iret)
         CASE(316:322)  !swan max/min values
            iret = nf_def_var_deflate(dat%ncid, dat%max_nodal_data_id,
     &             1, 1, 2)
            CALL check_err(iret)
         CASE DEFAULT ! should be unreachable
            write(scratchMessage,
     &       '("No netCDF for files with unit number ",i0,".")') descript1 % lun
            call allMessage(ERROR,scratchMessage)
         END SELECT
      ENDIF
#endif
C
C     RJW added 9/13/2010
C     to include time atributes in global data files
C     Define time attributes
      CALL defineTimeAttributes(dat%ncid, dat%myTime)
C
C     define metadata and selected fort.15 parameters in netcdf file
      CALL defineMetaData(dat%ncid)
C
C     Leave define mode
      iret = nf90_enddef(dat%ncid)
      CALL check_err(iret)
C
C     write mesh to netcdf file
      CALL putMeshVariables(dat%ncid,dat%myMesh)
      !
      ! jgf52.08.11: Initialize max/min values b/c they will be read in
      ! at the start of a hotstart run. The square of the initialized
      ! values will be compared with the sum of the squares of
      ! vector values in the live running code. Uninitialized values
      ! of -99999 become large and positive when squared, meaning that
      ! we will never be able set a new peak in the running code if
      ! we leave the max/min values uninitialized.
      IF ( (dat%myFile%createFile.eqv..true.).and.
     &     (descript1%readMaxMin.eqv..true.)) THEN
          iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id,defaultValue)
          CALL check_err(iret)
      ENDIF
C
      IF (any(nodalLunsHA.eq.descript1 % lun)) THEN
        CALL addHarmonicAnalysisParametersToNetcdfFile(dat%ncid)
      ENDIF
C     now close the initialized netcdf file
      iret = nf90_close(dat%ncid)
      CALL check_err(iret)

      call logMessage(INFO,'Initialized "'//trim(dat%myFile%filename)//
     &   '" file.')

#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE initNodalDataFile
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     SUBROUTINE defineHarmonicAnalysisParametersInNetcdfFile
C-----------------------------------------------------------------------
C     Defines the variables that carry frequency, nodal factor, and
C     equilibrium arg in the netcdf output file
C-----------------------------------------------------------------------
      SUBROUTINE defineHarmonicAnalysisParametersInNetcdfFile(descript1,
     &                                                        nodaldat,
     &                                                        stationdat)
        USE GLOBAL, ONLY : OutputDataDescript_t
        IMPLICIT NONE
        type(nodalData), intent(inout),optional :: nodaldat
        type(stationData), intent(inout),optional :: stationdat
        type(OutputDataDescript_t), intent(inout) :: descript1
        INTEGER :: iret,dmy

        IF(PRESENT(nodalDat))THEN
         nodalDat%myMesh%num_v_nodes = descript1 % num_items_per_record
         CALL check_err(nf90_def_dim(nodalDat%ncid, 'constlen', 10,
     &        nodalDat%myMesh%constl_dim_id))
         CALL check_err(nf90_def_dim(nodalDat%ncid, 'num_const',
     &        nodalDat%myMesh%num_v_nodes, nodalDat%myMesh%num_v_nodes_dim_id))
!        Put the constituent names in
         nodalDat%nodal_data_dims(1) = nodalDat%myMesh%constl_dim_id
         nodalDat%nodal_data_dims(2) = nodalDat%myMesh%num_v_nodes_dim_id
         iret = nf90_def_var(nodalDat%ncid, 'const', NF90_CHAR,
     &          nodalDat%nodal_data_dims, nodalDat%myMesh%const_id)
         CALL check_err(iret)
         CALL check_err(nf90_put_att(nodalDat%ncid, nodalDat%myMesh%const_id,
     &        'long_name',
     &        'names of the tidal harmonic constituents'))
         CALL check_err(nf90_def_var(nodalDat%ncid,'frequency',
     &          NF90_DOUBLE,nodalDat%myMesh%num_v_nodes_dim_id,dmy))
         CALL check_err(nf90_put_att(nodalDat%ncid,dmy,
     &          'long_name','frequency of harmonic constituents'))
         CALL check_err(nf90_put_att(nodalDat%ncid,dmy,'units','rad/s'))
         CALL check_err(nf90_def_var(nodalDat%ncid,'nodal_factor',
     &          NF90_DOUBLE,nodalDat%myMesh%num_v_nodes_dim_id,dmy))
         CALL check_err(nf90_put_att(nodalDat%ncid,dmy,'long_name',
     &          'nodal factor of harmonic constituents'))
         CALL check_err(nf90_def_var(nodalDat%ncid,'equilibrium_argument',
     &          NF90_DOUBLE,nodalDat%myMesh%num_v_nodes_dim_id,dmy))
         CALL check_err(nf90_put_att(nodalDat%ncid,dmy,'long_name',
     &          'equilibrium argument of harmonic constituents'))
         CALL check_err(nf90_put_att(nodalDat%ncid,dmy,'units','deg'))
        ELSEIF(PRESENT(stationDat))THEN
         CALL check_err(nf90_def_dim(stationDat%ncid, 'constlen', 10,
     &        stationDat%constl_dim_id))
         CALL check_err(nf90_def_dim(stationDat%ncid, 'num_const',
     &        stationDat%num_v_nodes, stationDat%num_v_nodes_dim_id))
!        Put the constituent names in
         stationDat%station_dims(1) = stationDat%constl_dim_id
         stationDat%station_dims(2) = stationDat%num_v_nodes_dim_id
         iret = nf90_def_var(stationDat%ncid, 'const', NF90_CHAR,
     &          stationDat%station_dims, stationDat%const_id)
         CALL check_err(iret)
         CALL check_err(nf90_put_att(stationDat%ncid, stationDat%const_id,
     &        'long_name',
     &        'names of the tidal harmonic constituents'))
         CALL check_err(nf90_def_var(stationDat%ncid,'frequency',
     &          NF90_DOUBLE,stationDat%num_v_nodes_dim_id,dmy))
         CALL check_err(nf90_put_att(stationDat%ncid,dmy,
     &          'long_name','frequency of harmonic constituents'))
         CALL check_err(nf90_put_att(stationDat%ncid,dmy,'units','rad/s'))
         CALL check_err(nf90_def_var(stationDat%ncid,'nodal_factor',
     &          NF90_DOUBLE,stationDat%num_v_nodes_dim_id,dmy))
         CALL check_err(nf90_put_att(stationDat%ncid,dmy,'long_name',
     &          'nodal factor of harmonic constituents'))
         CALL check_err(nf90_def_var(stationDat%ncid,'equilibrium_argument',
     &          NF90_DOUBLE,stationDat%num_v_nodes_dim_id,dmy))
         CALL check_err(nf90_put_att(stationDat%ncid,dmy,'long_name',
     &          'equilibrium argument of harmonic constituents'))
         CALL check_err(nf90_put_att(stationDat%ncid,dmy,'units','deg'))
        ENDIF

        RETURN

C-----------------------------------------------------------------------
      END SUBROUTINE defineHarmonicAnalysisParametersInNetcdfFile
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     SUBROUTINE addHarmonicAnalysisParametersToNetcdfFile
C-----------------------------------------------------------------------
C     Adds the variable data for frequency, nodal factor, and
C     equilibrium arg in the netcdf output file
C-----------------------------------------------------------------------
      SUBROUTINE addHarmonicAnalysisParametersToNetcdfFile(ncid)
         USE GLOBAL, ONLY : OutputDataDescript_t
         USE HARM, ONLY: NAMEFR,HAFREQ,HAFF,HAFACE
         IMPLICIT NONE
         INTEGER,INTENT(IN) :: ncid
         INTEGER :: iret,dmy

         CALL check_err(nf90_inq_varid(ncid,"const",dmy))
         CALL check_err(nf90_put_var(ncid, dmy, NAMEFR(:) ))
         CALL check_err(nf90_inq_varid(ncid,"frequency",dmy))
         CALL check_err(nf90_put_var(ncid,dmy,HAFREQ))
         CALL check_err(nf90_inq_varid(ncid,"nodal_factor",dmy))
         CALL check_err(nf90_put_var(ncid,dmy,HAFF))
         CALL check_err(nf90_inq_varid(ncid,"equilibrium_argument",dmy))
         CALL check_err(nf90_put_var(ncid,dmy,HAFACE))

         RETURN

C-----------------------------------------------------------------------
      END SUBROUTINE addHarmonicAnalysisParametersToNetcdfFile
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E  C R E A T E  N E T C D F  O U T P U T  F I L E
C-----------------------------------------------------------------------
C     jgf49.22 Creates a new netcdf output file if needed.
C-----------------------------------------------------------------------
      SUBROUTINE createNetCDFOutputFile(ncid, myFile, myTime,
     &                                  descript, ret_error)
      USE SIZES, ONLY : globaldir
      USE GLOBAL, ONLY : OutputDataDescript_t, IHOT, DEBUG, ERROR,
     &                   screenUnit, allMessage, scratchMessage
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      IMPLICIT NONE
      INTEGER, intent(out) :: ncid
      TYPE(fileData), intent(inout) :: myFile
      TYPE(timeData), intent(inout) :: myTime
      TYPE(OutputDataDescript_t), intent(inout) :: descript
      LOGICAL, intent(out) :: ret_error
C
      INTEGER iret
C
      call setMessageSource("createNetCDFOutputFile")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ret_error = .false.
      myFile%createFile = .false.
      myFile%fileFound = .false.

      ! create file name
      select case(descript%lun)
      case(67)
         myFile%filename = trim(globaldir)//'/fort.67.nc'
      case(68)
         myFile%filename = trim(globaldir)//'/fort.68.nc'
      case default
         ! trim(globaldir) has already been prepended to
         ! descript%file_name
         myFile%filename = trim(descript%file_name) // '.nc'
      end select
      !jgfdebug
      !write(6,*) 'myFile%filename=',trim(myFile%filename),
      !&  ' descript%file_name=',trim(descript%file_name)
C
C     jgf49.17.02: Simplified the criteria for creating a new netCDF
C     output file: coldstart, overwrite upon hotstart, or output
C     that does not already exist. These criteria do not apply to
C     netcdf hotstart files; we will always write those when called
C     upon to create them.
      INQUIRE(FILE=myFile%FILENAME,EXIST=myFile%fileFound)
      IF ((descript%lun.ne.67).and.(descript%lun.ne.68)) THEN
         IF ((IHOT.EQ.0).OR.(descript%specifier.lt.0).OR.
     &      (myFile%fileFound.eqv..false.)) THEN
#ifdef CMPI
            ! jgf49.31 when this subroutine is called by ADCIRC running in
            ! parallel, it should never create a new file, since that
            !  is the job of adcprep ... the file cannot be created
            ! here as a last resort since none of the processors have
            ! access to the full domain mesh and control files, whose
            ! data must also be stored in the netcdf output file
            IF (myFile%fileFound.eqv..false.) THEN
               write(scratchMessage,'(a,a,a)')
     &            "The NetCDF output file '",TRIM(myFile%FILENAME),
     &            "' was not found. "//
     &            "It should have been created by adcprep."
               call allMessage(ERROR,scratchMessage)
               ret_error = .true.
            ENDIF
#else
            ! these lines are executed by serial adcirc and adcprep
            !
            ! jgf52.08.25: It may be a non-NetCDF min/max will be written,
            ! but a netcdf min/max file of the same type may be read. In this
            ! case, don't create a  netcdf min/max file of this type.
            if ( (abs(descript%specifier).eq.NETCDF3) .or.
     &           (abs(descript%specifier).eq.NETCDF4) ) then
               myFile%createFile = .true.
               myFile%record_counter = 0
            endif
#endif
         ENDIF
      ELSE
         ! these lines are executed to create netcdf hotstart files
#ifdef CMPI
         IF (myFile%fileFound.eqv..true.) THEN
            call screenMessage(DEBUG,
     &      "Hotstart file was created by adcprep.")
            myFile%createFile = .false.
         ELSE
             call screenMessage(ERROR,"Hotstart file is missing.")
             call screenMessage(ERROR,
     &          "It should have been created by adcprep.")
             ret_error = .true.
         ENDIF
#else
         myFile%createFile = .true.
         myFile%record_counter = 0
#endif
      ENDIF
C
C     RETURN if we don't need to create a file.
      IF (myFile%createFile.eqv..false.) THEN
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         RETURN
      ENDIF
C
C     set the file format (netcdf3, or netcdf4/hdf5) ... the data model will
C     be classic (netcdf3) in any case
      select case(abs(descript%specifier))
      case(3,367,368)
         myFile%ncformat = NF90_CLOBBER
#ifdef HAVE_NETCDF4
      case(5,567,568)
         myFile%ncformat = ior(NF90_CLASSIC_MODEL,NF90_HDF5)
#else
      case(5,567,568)
         write(scratchMessage,'(A)')
     &    "This ADCIRC executable was compiled with the "
     &    //trim(nf90_inq_libvers())//" netcdf library."
         call allMessage(INFO,scratchMessage)
         write(scratchMessage,'(A,I3,A)') "File format specifier '",
     &      descript%specifier,"' requires NetCDF version 4."
         call allMessage(ERROR,scratchMessage)
         write(scratchMessage,'(A)')
     &    "It also requires the setting of NETCDF=enable and "
     &    //"NETCDF4=enable on the make command line."
         call allMessage(ERROR,scratchMessage)
         write(scratchMessage,'(A)')
     &    "You must recompile ADCIRC to use NetCDF4 formatting."
         call allMessage(ERROR,scratchMessage)
         ret_error = .true.
#endif
      case default
         write(scratchMessage,'(A,I3,A)') "File format specifier '",
     &      descript%specifier,"' is not valid."
         call allMessage(ERROR,scratchMessage)
         ret_error = .true.
      end select

      if ( ret_error.eqv..false. ) then
         iret = nf90_create(myFile%FILENAME, myFile%ncformat, ncid)
         CALL check_err(iret)
         ! Define time
         iret = nf90_def_dim(ncid, 'time', nf_unlimited,
     &       myTime%timenc_dim_id)
         CALL check_err(iret)
         myTime%timenc_dims(1) = myTime%timenc_dim_id
         iret = nf90_def_var(ncid, 'time', NF90_DOUBLE,
     &       myTime%timenc_dims, myTime%timenc_id)
         CALL check_err(iret)
      endif
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE createNetCDFOutputFile
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   I N I T   N E T C D F   C O O R D
C-----------------------------------------------------------------------
C     jgf49.17.02 Initializes NetCDF coordinates.
C-----------------------------------------------------------------------
      SUBROUTINE initNetCDFCoord(myMesh)
      USE GLOBAL, ONLY : NP_G, NE_G, RAD2DEG
      USE MESH, ONLY : X, Y, SLAM, SFEA, ICS, NM,
     &  SLAM0, SFEA0
      USE BOUNDARIES, ONLY : NBOU, NVEL, NOPE, NBVV, NVDLL, NBDV, NVELL,
     &  NETA, IBTYPEE, IBTYPE
      IMPLICIT NONE
      type(meshStructure), intent(inout) :: myMesh
      INTEGER :: i, j, k, kj  ! array indices
      INTEGER :: max_seg_nodes  ! size of longest list of nodes on a bdry seg
C
      call setMessageSource("initNetCDFCoord")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      myMesh%nopenc = nope
      myMesh%nbounc = nbou
      myMesh%netanc = neta
      myMesh%nvelnc = nvel
      myMesh%max_nvdllnc = maxval(nvdll) ! dimension of seg with most nodes
      myMesh%max_nvellnc = maxval(nvell) !WJP*2 dimension of seg with most nodes
C
      ALLOCATE(myMesh%xnc(NP_G))
      ALLOCATE(myMesh%ync(NP_G))
      ALLOCATE(myMesh%nvdllnc(myMesh%nopenc))
      ALLOCATE(myMesh%ibtypeenc(myMesh%nopenc))
      ALLOCATE(myMesh%ibtypenc(myMesh%nbounc))
      ALLOCATE(myMesh%nvellnc(myMesh%nbounc))
!     WJP to make nbdvnc arrays smaller:
      ALLOCATE(myMesh%nbdvnc(myMesh%netanc)) !(myMesh%nopenc,myMesh%max_nvdllnc))
C      ALLOCATE(myMesh%nbdvnc(myMesh%nopenc,myMesh%netanc))
!     WJP to make nbvvnc arrays smaller
      ALLOCATE(myMesh%nbvvnc(myMesh%nvelnc)) !(myMesh%nbounc,myMesh%max_nvellnc))
C      ALLOCATE(myMesh%nbvvnc(myMesh%nbounc,myMesh%nvelnc))
      ALLOCATE(myMesh%element(myMesh%nface_len,NE_G))
      ALLOCATE(myMesh%nmnc(NE_G,myMesh%nface_len))
C
C     Store nodal coordinates
      if (ics.eq.1) then
         myMesh%xnc = X
         myMesh%ync = Y
      else
         myMesh%xnc = SLAM * RAD2DEG ! convert back to degrees
         myMesh%ync = SFEA * RAD2DEG
      endif
!
!     elevation specified boundary forcing segments
      myMesh%nvdllnc = 0
      myMesh%nbdvnc = 0
      kj = 0 ! WJP initializing linear index
      DO k=1,myMesh%nopenc
         myMesh%nvdllnc(k) = nvdll(k)
         myMesh%ibtypeenc(k) = ibtypee(k)
         DO j=1,myMesh%nvdllnc(k)
            ! WJP counting the linear index
            kj = kj + 1
            myMesh%nbdvnc(kj) = nbdv(k,j)
            !myMesh%nbdvnc(k,j) = nbdv(k,j)
         END DO
      END DO
!
!     normal flow (discharge) specified boundary segments
C
      myMesh%nvellnc = 0
      myMesh%nbvvnc = 0
      kj = 0 !WJP initializing linear index
      DO k=1,myMesh%nbounc
         myMesh%nvellnc(k) = nvell(k)
         myMesh%ibtypenc(k) = ibtype(k)
         DO J=1,myMesh%nvellnc(k)
            ! WJP counting the linear index
            kj = kj + 1
            myMesh%nbvvnc(kj) = nbvv(k,j)
            !myMesh%nbvvnc(k,j)=nbvv(k,j)
         END DO
      END DO
!
      myMesh%nmnc=NM
!     Switch order in array for NETCDF
C       write(*,*)  myMesh%num_elems
C       write(*,*)  myMesh%nface_len
      do i=1, NE_G
         do j=1, myMesh%nface_len
            myMesh%element(j,i) = myMesh%nmnc(i,j)
         end do
      end do
!       write(99,*) myMesh%nmnc
!       write(98,*) NM
!       write(97,*) myMesh%element
!      stop
C
      myMesh%initialized = .true.
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE initNetCDFCoord
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E   W R I T E  O U T  A R R A Y  N E T C D F
C-----------------------------------------------------------------------
C     jgf48.03 This subroutine was created from cf's code in timestep.F to
C     write output files in NetCDF format.
C-----------------------------------------------------------------------
      SUBROUTINE writeOutArrayNetCDF(lun, timesec, it, descript1,
     &   descript2, descript3, descript4 )
      USE GLOBAL, ONLY : OutputDataDescript_t, ScreenUnit
C
      IMPLICIT NONE
C     args
      INTEGER, intent(in) :: lun ! logical unit number of file to write to
      REAL(8), intent(in) :: timesec ! seconds since cold start
      INTEGER, intent(in) :: it   ! number of time steps since cold start
      type(OutputDataDescript_t), intent(in) :: descript1 !describes output data
      type(OutputDataDescript_t), intent(in), optional :: descript2 !describes output data
      type(OutputDataDescript_t), intent(in), optional :: descript3 !describes output data
      type(OutputDataDescript_t), intent(in), optional :: descript4 !describes output data
C
      call setMessageSource("writeOutArrayNetCDF")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      SELECT CASE(lun)
         CASE(41)
            CALL writeStationData(densityStations3D, lun, descript1,
     &           timesec, descript2, descript3)
         CASE(42)
            CALL writeStationData(velocityStations3D, lun, descript1,
     &           timesec, descript2, descript3)
         CASE(43)
            CALL writeStationData(turbulenceStations3D, lun, descript1,
     &           timesec, descript2, descript3)
         CASE(44)
            CALL writeNodalData(density3D, lun, descript1,
     &           timesec, descript2, descript3)
         CASE(45)
            CALL writeNodalData(velocity3D, lun, descript1,
     &           timesec, descript2, descript3)
         CASE(46)
            CALL writeNodalData(turbulence3D, lun, descript1,
     &           timesec, descript2, descript3)
         CASE(47)
            CALL writeNodalData(futureSurfaceTemperature, lun,
     &           descript1, timesec)
         !WJP 02.20.2018: added fort.51-54 harmonic outputs
         CASE(51)
            CALL writeStationData(HAelevSta, lun, descript1, timesec,
     &           descript2)
         CASE(52)
            CALL writeStationData(HAvelSta, lun, descript1, timesec,
     &           descript2,descript3,descript4)
         CASE(53)
            CALL writeNodalData(HAelev, lun, descript1, timesec,
     &           descript2)
         CASE(54)
            CALL writeNodalData(HAVel, lun, descript1, timesec,
     &           descript2,descript3,descript4)
         CASE(61)
            CALL writeStationData(elevSta, lun, descript1, timesec)
         CASE(62)
            CALL writeStationData(velSta, lun, descript1, timesec)
         CASE(63)
            CALL writeNodalData(elev, lun, descript1, timesec)
         CASE(64)
            CALL writeNodalData(currentVel, lun, descript1, timesec)
         CASE(71)
            CALL writeStationData(prSta, lun, descript1, timesec)
         CASE(72)
            CALL writeStationData(wVelSta, lun, descript1, timesec)
         CASE(73)
            CALL writeNodalData(pr, lun, descript1, timesec)
         CASE(74)
            CALL writeNodalData(windVel, lun, descript1, timesec)
C++ GML 20200528 
         CASE(91)
            CALL writeStationData(IceSta, lun, descript1, timesec)
         CASE(93)
            CALL writeNodalData(CICEAF, lun, descript1, timesec)
C++
         CASE(173)
            CALL writeNodalData(wDrag, lun, descript1, timesec)
         CASE(77)
            CALL writeNodalData(weirElev, lun, descript1, timesec)
         CASE(90)
            CALL writeNodalData(tau0nc, lun, descript1, timesec)
         CASE(164)
            CALL writeNodalData(rads, lun, descript1, timesec)
#if defined CSWAN || defined ADCSWAN
         CASE(301)
            CALL writeNodalData(sw_hs, lun, descript1, timesec)
         CASE(302)
            CALL writeNodalData(sw_dir, lun, descript1, timesec)
         CASE(303)
            CALL writeNodalData(sw_tm01, lun, descript1, timesec)
         CASE(304)
            CALL writeNodalData(sw_tps, lun, descript1, timesec)
         CASE(305)
            CALL writeNodalData(sw_wind, lun, descript1, timesec)
         CASE(306)
            CALL writeNodalData(sw_tm02, lun, descript1, timesec)
         CASE(307)
            CALL writeNodalData(sw_tmm10, lun, descript1, timesec)
#endif
         CASE(311)
            CALL writeNodalData(EtaMax, lun, descript1, timesec)
         CASE(312)
            CALL writeNodalData(UMax, lun, descript1, timesec)
         CASE(313)
            CALL writeNodalData(PrMin, lun, descript1, timesec)
         CASE(314)
            CALL writeNodalData(WVMax, lun, descript1, timesec)
         CASE(315)
            CALL writeNodalData(RSMax, lun, descript1, timesec)
#if defined CSWAN || defined ADCSWAN
         CASE(316)
            CALL writeNodalData(sw_hs_max, lun, descript1, timesec)
         CASE(317)
            CALL writeNodalData(sw_dir_max, lun, descript1, timesec)
         CASE(318)
            CALL writeNodalData(sw_tm01_max, lun, descript1, timesec)
         CASE(319)
            CALL writeNodalData(sw_tps_max, lun, descript1, timesec)
         CASE(320)
            CALL writeNodalData(sw_wind_max, lun, descript1, timesec)
         CASE(321)
            CALL writeNodalData(sw_tm02_max, lun, descript1, timesec)
         CASE(322)
            CALL writeNodalData(sw_tmm10_max, lun, descript1, timesec)
#endif
         CASE(400)
            CALL writeNodalData(inTime, lun, descript1, timesec)
         CASE(401)
            CALL writeNodalData(maxInDep, lun, descript1, timesec)
         CASE(402)
            CALL writeNodalData(coldDry, lun, descript1, timesec)
         CASE(403)
            CALL writeNodalData(eRisInun, lun, descript1, timesec)
         CASE(404)
            CALL writeNodalData(evrDry, lun, descript1, timesec)

         CASE(108)
            CALL writeNodalData(dynamiccorrection, lun, descript1, timesec)
         CASE(109)
            CALL writeStationData(dynamiccorrectionSta, lun, descript1, timesec)

         CASE DEFAULT
            write(scratchMessage,
     &       '("No netCDF for files with unit number ",i0,".")') lun
            call allMessage(ERROR,scratchMessage)
      END SELECT
C

#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE writeOutArrayNetCDF
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E    W R I T E  S T A T I O N  D A T A
C-----------------------------------------------------------------------
C     jgf49.17.02 Writes data to station file.
C-----------------------------------------------------------------------
      SUBROUTINE writeStationData(sta, lun, descript1, timesec,
     &   descript2, descript3, descript4)
      USE SIZES, ONLY : MNPROC
      USE GLOBAL, ONLY : SNAMLEN, OutputDataDescript_t, scratchMessage,
     &  IDEN
      IMPLICIT NONE
C
      type(stationData), intent(inout) :: sta
      INTEGER, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: descript1
      type(OutputDataDescript_t), intent(in), optional :: descript2
      type(OutputDataDescript_t), intent(in), optional :: descript3
      type(OutputDataDescript_t), intent(in), optional :: descript4
      REAL(8), intent(in) :: timesec
C
      REAL(8), allocatable :: storedTimes(:) ! array of time values in file
      LOGICAL :: timeFound ! true if current time is in array of stored times
      INTEGER i,j
      INTEGER counti(1), starti(1), n
      INTEGER kount(2), start(2)
      INTEGER kount3D(4), start3D(3)
      INTEGER iret  ! success or failure of netcdf call
C
      call setMessageSource("writeStationData")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
C     jgf50.13: if netcdf output was requested, but there are no stations,
C     don't write to the file (it doesn't exist); just return.
      IF ( descript1 % num_fd_records.eq.0 ) THEN
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         RETURN
      ENDIF

      iret = nf90_open(sta%myFile%FILENAME, NF_WRITE, sta%ncid)
      CALL check_err(iret)
      sta%myTime%timenc(sta%myTime%timenc_len)=timesec

      CALL setRecordCounterAndStoreTime(sta%ncid, sta%myFile,
     &   sta%myTime)
      !
      ! get number of vertical nodes for 3D stations
      IF ((lun.eq.41).or.(lun.eq.42).or.(lun.eq.43)) THEN
         iret=nf90_inq_dimid(sta%ncid,"num_v_nodes",sta%num_v_nodes_dim_id)
         CALL check_err(iret)
         iret=nf90_inquire_dimension(sta%ncid, sta%num_v_nodes_dim_id,
     &        len=sta%num_v_nodes)
         CALL check_err(iret)
         !
         ! Set up the 3D netcdf data extents
         kount3D(1)=sta%num_stations
         kount3D(2)=sta%num_v_nodes
         kount3D(3)=sta%myTime%timenc_len
         start3D(1)=1
         start3D(2)=1
         start3D(3)=sta%myFile%record_counter
      ELSE
         ! Set up the 2D netcdf data extents
         kount(1)=sta%num_stations
         kount(2)=sta%myTime%timenc_len
         start(1)=1
         start(2)=sta%myFile%record_counter
      ENDIF
C
      SELECT CASE(lun)

      CASE(41)
         iret=nf90_inq_varid(sta%ncid, "sigmat", sta%u_station_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id,
     &          descript1%array2D, start3D, kount3D)
         ELSE
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id,
     &          descript1%array2D_g, start3D, kount3D)
         ENDIF
         CALL check_err(iret)

         IF ((iden.eq.2).or.(iden.eq.4)) THEN
            iret=nf90_inq_varid(sta%ncid, "salinity",
     &         sta%v_station_data_id)
            CALL check_err(iret)
            IF (MNPROC.eq.1) THEN
               iret = nf90_put_var(sta%ncid, sta%v_station_data_id,
     &             descript2%array2D, start3D, kount3D)
            ELSE
               iret = nf90_put_var(sta%ncid, sta%v_station_data_id,
     &             descript2%array2D_g, start3D, kount3D)
            ENDIF

            CALL check_err(iret)
         ENDIF
         IF ((iden.eq.3).or.(iden.eq.4)) THEN
            iret=nf90_inq_varid(sta%ncid, "temperature",
     &          sta%w_station_data_id)
            CALL check_err(iret)
            IF (MNPROC.eq.1) THEN
               iret = nf90_put_var(sta%ncid, sta%w_station_data_id,
     &            descript3%array2D, start3D, kount3D)
            ELSE
               iret = nf90_put_var(sta%ncid, sta%w_station_data_id,
     &            descript3%array2D_g, start3D, kount3D)
            ENDIF
            CALL check_err(iret)
         ENDIF
      CASE(42)
         iret=nf90_inq_varid(sta%ncid, "u-vel3D", sta%u_station_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(sta%ncid, "v-vel3D", sta%v_station_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(sta%ncid, "w-vel3D", sta%w_station_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id,
     &         descript1%array2D, start3D, kount3D)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%v_station_data_id,
     &         descript2%array2D, start3D, kount3D)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%w_station_data_id,
     &         descript3%array2D, start3D, kount3D)
            CALL check_err(iret)
         ELSE
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id,
     &         descript1%array2D_g, start3D, kount3D)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%v_station_data_id,
     &         descript2%array2D_g, start3D, kount3D)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%w_station_data_id,
     &         descript3%array2D_g, start3D, kount3D)
            CALL check_err(iret)
         ENDIF

      CASE(43)
         iret=nf90_inq_varid(sta%ncid, "q20", sta%u_station_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(sta%ncid, "l", sta%v_station_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(sta%ncid, "ev", sta%w_station_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id,
     &         descript1%array2D, start3D, kount3D)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%v_station_data_id,
     &         descript2%array2D, start3D, kount3D)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%w_station_data_id,
     &         descript3%array2D, start3D, kount3D)
            CALL check_err(iret)
         ELSE
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id,
     &         descript1%array2D_g, start3D, kount3D)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%v_station_data_id,
     &         descript2%array2D_g, start3D, kount3D)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%w_station_data_id,
     &         descript3%array2D_g, start3D, kount3D)
            CALL check_err(iret)
         ENDIF
      ! WJP 02.20.2018 Adding fort.51-52.nc capabilities
       CASE(51)
         iret=nf90_inq_varid(sta%ncid, "amp", sta%station_ha_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(sta%ncid, "phs", sta%station_hg_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(sta%ncid, sta%station_ha_data_id,
     &           descript1%array2D)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%station_hg_data_id,
     &           descript2%array2D)
            CALL check_err(iret)
         ELSE
            iret = nf90_put_var(sta%ncid, sta%station_ha_data_id,
     &         descript1%array2D_g)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%station_hg_data_id,
     &         descript2%array2D_g)
            CALL check_err(iret)
         ENDIF

      CASE(52)
         iret=nf90_inq_varid(sta%ncid, "u-vel-amp",
     &        sta%ha_u_station_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(sta%ncid, "u-vel-phs",
     &        sta%hg_u_station_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(sta%ncid, "v-vel-amp",
     &        sta%ha_v_station_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(sta%ncid, "v-vel-phs",
     &        sta%hg_v_station_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(sta%ncid, sta%ha_u_station_data_id,
     &         descript1%array2D)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%hg_u_station_data_id,
     &         descript2%array2D)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%ha_v_station_data_id,
     &         descript3%array2D)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%hg_v_station_data_id,
     &         descript4%array2D)
            CALL check_err(iret)
         ELSE
            iret = nf90_put_var(sta%ncid, sta%ha_u_station_data_id,
     &         descript1%array2D_g)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%hg_u_station_data_id,
     &         descript2%array2D_g)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%ha_v_station_data_id,
     &         descript3%array2D_g)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%hg_v_station_data_id,
     &         descript4%array2D_g)
            CALL check_err(iret)
         ENDIF

      CASE(61)
         iret=nf90_inq_varid(sta%ncid, "zeta", sta%station_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(sta%ncid, sta%station_data_id,
     &          descript1%array, start, kount)
         ELSE
            iret = nf90_put_var(sta%ncid, sta%station_data_id,
     &          descript1%array_g, start, kount)
         ENDIF
         CALL check_err(iret)
      CASE(62)
         iret=nf90_inq_varid(sta%ncid, "u-vel", sta%u_station_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(sta%ncid, "v-vel", sta%v_station_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id,
     &          descript1%array, start, kount)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%v_station_data_id,
     &          descript1%array2, start, kount)
         ELSE
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id,
     &          descript1%array_g, start, kount)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%v_station_data_id,
     &          descript1%array2_g, start, kount)
         ENDIF
      CASE(71)
         iret=nf90_inq_varid(sta%ncid, "pressure", sta%station_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(sta%ncid, sta%station_data_id,
     &          descript1%array, start, kount)
         ELSE
            iret = nf90_put_var(sta%ncid, sta%station_data_id,
     &          descript1%array_g, start, kount)
         ENDIF
         CALL check_err(iret)
      CASE(72)
           iret=nf90_inq_varid(sta%ncid,"windx",sta%u_station_data_id)
           CALL check_err(iret)
           iret=nf90_inq_varid(sta%ncid,"windy",sta%v_station_data_id)
           CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id,
     &          descript1%array, start, kount)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%v_station_data_id,
     &          descript1%array2, start, kount)
         ELSE
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id,
     &          descript1%array_g, start, kount)
            CALL check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%v_station_data_id,
     &          descript1%array2_g, start, kount)
         ENDIF

C++ GML 20200605
      CASE(91)
         iret=nf90_inq_varid(sta%ncid, "iceaf", sta%station_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(sta%ncid, sta%station_data_id,
     &          descript1%array, start, kount)
         ELSE
            iret = nf90_put_var(sta%ncid, sta%station_data_id,
     &          descript1%array_g, start, kount)
         ENDIF
         CALL check_err(iret)
C++

      CASE(109)
         iret=nf90_inq_varid(sta%ncid, "dynamicWaterlevelCorrection", sta%station_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(sta%ncid, sta%station_data_id,
     &          descript1%array, start, kount)
         ELSE
            iret = nf90_put_var(sta%ncid, sta%station_data_id,
     &          descript1%array_g, start, kount)
         ENDIF
         CALL check_err(iret)


      CASE DEFAULT
            write(scratchMessage,
     &       '("No netCDF for station files with unit number ",i0,".")') lun
            call allMessage(ERROR,scratchMessage)
      END SELECT
!
!     Close netCDF file
9999  iret = nf90_close(sta%ncid)
      CALL check_err(iret)

#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      end subroutine writeStationData
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E    W R I T E  N O D A L   D A T A
C-----------------------------------------------------------------------
C     jgf49.17.02 Writes data to a full domain file.
C-----------------------------------------------------------------------
      SUBROUTINE writeNodalData(dat, lun, descript1, timesec,
     &   descript2, descript3, descript4 )
      USE SIZES, ONLY : MNPROC, MYPROC
      USE GLOBAL, ONLY : SNAMLEN, OutputDataDescript_t, NODECODE,
     &                   scratchMessage, IDEN

      IMPLICIT NONE
C
      type(nodalData), intent(inout) :: dat
      INTEGER, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: descript1
      type(OutputDataDescript_t), intent(in), optional :: descript2
      type(OutputDataDescript_t), intent(in), optional :: descript3
      type(OutputDataDescript_t), intent(in), optional :: descript4
      REAL(8), intent(in) :: timesec
C
      INTEGER i,j, ii
      INTEGER counti(1), starti(1), n
      INTEGER kount(2), start(2)
      INTEGER kount3D(3), start3D(3)
      INTEGER iret ! success or failure of the netcdf call
      INTEGER num_values ! number of values in the array
      INTEGER,ALLOCATABLE  :: tempIntArray(:)
      REAL(SZ),ALLOCATABLE :: tempArray(:)
C
      call setMessageSource("writeNodalData")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      write(scratchMessage,'("Rank ",i0," writing to ",a,".")') myproc,
     &   trim(dat%myFile%filename)
      call allMessage(DEBUG,scratchMessage)

      iret = nf90_open(dat%myFile%FILENAME, NF_WRITE, dat%ncid)
      CALL check_err(iret)
      dat%myTime%timenc(dat%myTime%timenc_len)=timesec
C
      CALL setRecordCounterAndStoreTime(dat%ncid, dat%myFile,
     &   dat%myTime)
C
C     Set up the 2D netcdf data extents
      kount(1)=dat%myMesh%num_nodes
      kount(2)=dat%myTime%timenc_len
      start(1)=1
      start(2)=dat%myFile%record_counter

C
C     Set up the 3D netcdf data extents
      kount3D(1)=dat%myMesh%num_nodes
      kount3D(2)=dat%myMesh%num_v_nodes
      kount3D(3)=dat%myTime%timenc_len
      start3D(1)=1
      start3D(2)=1
      start3D(3)=dat%myFile%record_counter
C
C     Grab the data ids
      SELECT CASE(lun)
C
      CASE(44)
         iret=nf90_inq_varid(dat%ncid, "sigmat", dat%u_nodal_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(dat%ncid, dat%u_nodal_data_id,
     &          descript1%array2D, start3D, kount3D)
         ELSE
            iret = nf90_put_var(dat%ncid, dat%u_nodal_data_id,
     &          descript1%array2D_g, start3D, kount3D)
         ENDIF

         CALL check_err(iret)
         IF ((iden.eq.2).or.(iden.eq.4)) THEN
            iret=nf90_inq_varid(dat%ncid, "salinity",
     &          dat%v_nodal_data_id)
            CALL check_err(iret)
            IF (MNPROC.eq.1) THEN
               iret = nf90_put_var(dat%ncid, dat%v_nodal_data_id,
     &             descript2%array2D, start3D, kount3D)
            ELSE
               iret = nf90_put_var(dat%ncid, dat%v_nodal_data_id,
     &           descript2%array2D_g, start3D, kount3D)
            ENDIF
            CALL check_err(iret)
         ENDIF
         IF ((iden.eq.3).or.(iden.eq.4)) THEN
            iret=nf90_inq_varid(dat%ncid, "temperature",
     &         dat%w_nodal_data_id)
            CALL check_err(iret)
            IF (MNPROC.eq.1) THEN
               iret = nf90_put_var(dat%ncid, dat%w_nodal_data_id,
     &            descript3%array2D, start3D, kount3D)
            ELSE
               iret = nf90_put_var(dat%ncid, dat%w_nodal_data_id,
     &            descript3%array2D_g, start3D, kount3D)
            ENDIF
            CALL check_err(iret)
         ENDIF
      CASE(45)
         iret=nf90_inq_varid(dat%ncid, "u-vel3D", dat%u_nodal_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(dat%ncid, "v-vel3D", dat%v_nodal_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(dat%ncid, "w-vel3D", dat%w_nodal_data_id)
         CALL check_err(iret)
      CASE(46)
         iret=nf90_inq_varid(dat%ncid, "q20", dat%u_nodal_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(dat%ncid, "l", dat%v_nodal_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(dat%ncid, "ev", dat%w_nodal_data_id)
         CALL check_err(iret)
      CASE(47)
         iret=nf90_inq_varid(dat%ncid,"qsurfkp1",dat%nodal_data_id)
         CALL check_err(iret)
      ! WJP 02.20.2018 Adding capabilities for fort.53-54.nc
      CASE(53)
         iret=nf90_inq_varid(dat%ncid,"amp",dat%ha_nodal_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(dat%ncid,"phs",dat%hg_nodal_data_id)
         CALL check_err(iret)
      CASE(54)
         iret=nf90_inq_varid(dat%ncid, "u_amp", dat%ha_nodal_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(dat%ncid, "u_phs", dat%hg_nodal_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(dat%ncid, "v_amp", dat%v_ha_nodal_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(dat%ncid, "v_phs", dat%v_hg_nodal_data_id)
         CALL check_err(iret)
      CASE(63)
         iret=nf90_inq_varid(dat%ncid,"zeta",dat%nodal_data_id)
         CALL check_err(iret)
      CASE(64)
         iret=nf90_inq_varid(dat%ncid, "u-vel", dat%u_nodal_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(dat%ncid, "v-vel", dat%v_nodal_data_id)
         CALL check_err(iret)
      CASE(73)
         iret=nf90_inq_varid(dat%ncid, "pressure", dat%nodal_data_id)
         CALL check_err(iret)
      CASE(74)
         iret=nf90_inq_varid(dat%ncid,"windx",dat%u_nodal_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(dat%ncid,"windy",dat%v_nodal_data_id)
         CALL check_err(iret)
C++ GML 20200528
      CASE(93)
         iret=nf90_inq_varid(dat%ncid,"iceaf",dat%nodal_data_id)
         CALL check_err(iret)
C++
      CASE(173)
         iret=nf90_inq_varid(dat%ncid, "winddrag", dat%nodal_data_id)
         CALL check_err(iret)
      CASE(77)
         iret=nf90_inq_varid(dat%ncid,"weir_dz",dat%nodal_data_id)
         CALL check_err(iret)
      CASE(90)
         iret=nf90_inq_varid(dat%ncid, "tau0", dat%nodal_data_id)
         CALL check_err(iret)
      CASE(164)
        iret=nf90_inq_varid(dat%ncid, "radstress_x",dat%u_nodal_data_id)
        CALL check_err(iret)
        iret=nf90_inq_varid(dat%ncid, "radstress_y",dat%v_nodal_data_id)
        CALL check_err(iret)
#ifdef CSWAN
      CASE(301)
        iret=nf90_inq_varid(dat%ncid, "swan_HS", dat%nodal_data_id)
        CALL check_err(iret)
      CASE(302)
        iret=nf90_inq_varid(dat%ncid, "swan_DIR", dat%nodal_data_id)
        CALL check_err(iret)
      CASE(303)
        iret=nf90_inq_varid(dat%ncid, "swan_TM01", dat%nodal_data_id)
        CALL check_err(iret)
      CASE(304)
        iret=nf90_inq_varid(dat%ncid, "swan_TPS", dat%nodal_data_id)
        CALL check_err(iret)
      CASE(305)
        iret=nf90_inq_varid(dat%ncid, "swan_windx", dat%u_nodal_data_id)
        CALL check_err(iret)
        iret=nf90_inq_varid(dat%ncid, "swan_windy", dat%v_nodal_data_id)
        CALL check_err(iret)
      CASE(306)
        iret=nf90_inq_varid(dat%ncid, "swan_TM02", dat%nodal_data_id)
        CALL check_err(iret)
      CASE(307)
        iret=nf90_inq_varid(dat%ncid, "swan_TMM10", dat%nodal_data_id)
        CALL check_err(iret)
#endif

      CASE(311)
        iret=nf90_inq_varid(dat%ncid, "zeta_max", dat%max_nodal_data_id)
        CALL check_err(iret)
        iret=nf90_inq_varid(dat%ncid, "time_of_zeta_max",
     &        dat%time_max_nodal_data_id)
        CALL check_err(iret)
      CASE(312)
        iret=nf90_inq_varid(dat%ncid, "vel_max",
     &           dat%max_nodal_data_id)
        CALL check_err(iret)
        iret=nf90_inq_varid(dat%ncid, "time_of_vel_max",
     &           dat%time_max_nodal_data_id)
        CALL check_err(iret)
      CASE(313)
        iret=nf90_inq_varid(dat%ncid, "pressure_min",
     &          dat%max_nodal_data_id)
        CALL check_err(iret)
        iret=nf90_inq_varid(dat%ncid, "time_of_pressure_min",
     &          dat%time_max_nodal_data_id)
        CALL check_err(iret)
      CASE(314)
        iret=nf90_inq_varid(dat%ncid, "wind_max",
     &          dat%max_nodal_data_id)
        CALL check_err(iret)
        iret=nf90_inq_varid(dat%ncid, "time_of_wind_max",
     &          dat%time_max_nodal_data_id)
        CALL check_err(iret)
      CASE(315)
        iret=nf90_inq_varid(dat%ncid,"radstress_max",
     &           dat%max_nodal_data_id)
        CALL check_err(iret)
        iret=nf90_inq_varid(dat%ncid,"time_of_radstress_max",
     &           dat%time_max_nodal_data_id)
        CALL check_err(iret)
#ifdef CSWAN
      CASE(316)
        iret=nf90_inq_varid(dat%ncid, "swan_HS_max",
     &  dat%max_nodal_data_id)
        CALL check_err(iret)
      CASE(317)
        iret=nf90_inq_varid(dat%ncid, "swan_DIR_max",
     &  dat%max_nodal_data_id)
        CALL check_err(iret)
      CASE(318)
        iret=nf90_inq_varid(dat%ncid, "swan_TM01_max",
     &  dat%max_nodal_data_id)
        CALL check_err(iret)
      CASE(319)
        iret=nf90_inq_varid(dat%ncid, "swan_TPS_max",
     &  dat%max_nodal_data_id)
        CALL check_err(iret)
      CASE(320)
        iret=nf90_inq_varid(dat%ncid, "swan_wind_max",
     &  dat%max_nodal_data_id)
        CALL check_err(iret)
      CASE(321)
        iret=nf90_inq_varid(dat%ncid, "swan_TM02_max",
     &  dat%max_nodal_data_id)
        CALL check_err(iret)
      CASE(322)
        iret=nf90_inq_varid(dat%ncid, "swan_TMM10_max",
     &  dat%max_nodal_data_id)
        CALL check_err(iret)
#endif

      ! inundationtime.63
      CASE(400)
        iret=nf90_inq_varid(dat%ncid, "inun_time",
     &          dat%max_nodal_data_id)
        CALL check_err(iret)
        iret=nf90_inq_varid(dat%ncid, "onset_inun_time",
     &          dat%time_max_nodal_data_id)
        CALL check_err(iret)
      ! maxinundepth.63
      CASE(401)
        iret=nf90_inq_varid(dat%ncid, "inun_max",
     &          dat%max_nodal_data_id)
        CALL check_err(iret)
        iret=nf90_inq_varid(dat%ncid, "time_of_inun_max",
     &          dat%time_max_nodal_data_id)
        CALL check_err(iret)
      ! initiallydry.63
      CASE(402)
        iret=nf90_inq_varid(dat%ncid, "initiallydry",dat%nodal_data_id)
        CALL check_err(iret)
      ! endrisinginun.63
      CASE(403)
        iret=nf90_inq_varid(dat%ncid, "endrisinginun",dat%nodal_data_id)
        CALL check_err(iret)
      CASE(404)
        iret=nf90_inq_varid(dat%ncid, "everdried",
     &          dat%max_nodal_data_id)
        CALL check_err(iret)
        iret=nf90_inq_varid(dat%ncid, "time_of_everdried",
     &          dat%time_max_nodal_data_id)
        CALL check_err(iret)

      CASE(108)
         iret=nf90_inq_varid(dat%ncid,"dynamicWaterlevelCorrection",dat%nodal_data_id)
         CALL check_err(iret)


      CASE DEFAULT
         write(scratchMessage,
     &       '("No netCDF for files with unit number ",i0,".")') lun
         call allMessage(ERROR,scratchMessage)
      END SELECT
C
C     Write the array values
      SELECT CASE(lun)
      CASE(45,46)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(dat%ncid, dat%u_nodal_data_id,
     &          descript1%array2D, start3D, kount3D)
            CALL check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%v_nodal_data_id,
     &          descript2%array2D, start3D, kount3D)
            CALL check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%w_nodal_data_id,
     &          descript3%array2D, start3D, kount3D)
            CALL check_err(iret)
         ELSE
            iret = nf90_put_var(dat%ncid, dat%u_nodal_data_id,
     &          descript1%array2D_g, start3D, kount3D)
            CALL check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%v_nodal_data_id,
     &          descript2%array2D_g, start3D, kount3D)
            CALL check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%w_nodal_data_id,
     &          descript3%array2D_g, start3D, kount3D)
            CALL check_err(iret)
         ENDIF
      !WJP 02.20.2018 Adding capabilities for fort.53-54.nc
      CASE(53,54)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(dat%ncid, dat%ha_nodal_data_id,
     &          descript1%array2D)
            CALL check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%hg_nodal_data_id,
     &          descript2%array2D)
            CALL check_err(iret)
            IF (lun.eq.54) THEN
               iret = nf90_put_var(dat%ncid, dat%v_ha_nodal_data_id,
     &         descript3%array2D)
               CALL check_err(iret)
               iret = nf90_put_var(dat%ncid, dat%v_hg_nodal_data_id,
     &         descript4%array2D)
               CALL check_err(iret)
            ENDIF
         ELSE
            iret = nf90_put_var(dat%ncid, dat%ha_nodal_data_id,
     &          descript1%array2D_g)
            CALL check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%hg_nodal_data_id,
     &          descript2%array2D_g)
            CALL check_err(iret)
            IF (lun.eq.54) THEN
               iret = nf90_put_var(dat%ncid, dat%v_ha_nodal_data_id,
     &         descript3%array2D_g)
               CALL check_err(iret)
               iret = nf90_put_var(dat%ncid, dat%v_hg_nodal_data_id,
     &         descript4%array2D_g)
               CALL check_err(iret)
            ENDIF
         ENDIF
      CASE(316,317,318,319,321,322)  !swan max/min
         IF (MNPROC.eq.1) THEN ! SERIAL
            IF ( descript1%ConsiderWetDry.eqv..TRUE.) THEN
               allocate(tempArray(1:dat%myMesh%num_nodes))
               tempArray = merge(descript1%array,
     &              spread(descript1%alternate_value,1,kount(1)),nodecode.gt.0)
               iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id,tempArray)
            ELSE
               iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id,
     &            descript1%array)
            ENDIF
            CALL check_err(iret)
         ELSE
            iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id,
     &            descript1%array_g)
            CALL check_err(iret)
         ENDIF
      CASE(311,313)  !adcirc max/min only
         IF (MNPROC.eq.1) THEN ! SERIAL
            IF ( descript1%ConsiderWetDry.eqv..TRUE.) THEN
               allocate(tempArray(1:dat%myMesh%num_nodes))
               tempArray = merge(descript1%array,
     &              spread(descript1%alternate_value,1,kount(1)),nodecode.gt.0)
               iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id,tempArray)
               CALL check_err(iret)
               tempArray = merge(descript1%array2,
     &              spread(descript1%alternate_value,1,kount(1)),nodecode.gt.0)
               iret = nf90_put_var(dat%ncid, dat%time_max_nodal_data_id,tempArray)
               CALL check_err(iret)
            ELSE
               iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id,
     &            descript1%array)
               CALL check_err(iret)
               iret = nf90_put_var(dat%ncid, dat%time_max_nodal_data_id,
     &            descript1%array2)
               CALL check_err(iret)
            ENDIF
         ELSE
            iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id,
     &            descript1%array_g)
            CALL check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%time_max_nodal_data_id,
     &            descript1%array2_g)
            CALL check_err(iret)
         ENDIF
      CASE(312,314,315,400,401,404)  !adcirc max/min
         IF (MNPROC.eq.1) THEN ! SERIAL
            IF ( descript1%ConsiderWetDry.eqv..TRUE.) THEN
               allocate(tempArray(1:dat%myMesh%num_nodes))
               tempArray = merge(descript1%array,
     &              spread(descript1%alternate_value,1,kount(1)),nodecode.gt.0)
               iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id,tempArray)
               CALL check_err(iret)
               tempArray = merge(descript1%array2,
     &              spread(descript1%alternate_value,1,kount(1)),nodecode.gt.0)
               iret = nf90_put_var(dat%ncid, dat%time_max_nodal_data_id,tempArray)
               CALL check_err(iret)
            ELSE
               iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id,
     &             descript1%array)
               CALL check_err(iret)
               iret = nf90_put_var(dat%ncid, dat%time_max_nodal_data_id,
     &             descript1%array2)
               CALL check_err(iret)
            ENDIF
         ELSE                  ! PARALLEL
            iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id,
     &            descript1%array_g)
            CALL check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%time_max_nodal_data_id,
     &            descript1%array2_g)
            CALL check_err(iret)
         ENDIF
      CASE(320)  !swan max/min
         IF (MNPROC.eq.1) THEN ! SERIAL
            IF ( descript1%ConsiderWetDry.eqv..TRUE.) THEN
               allocate(tempArray(1:dat%myMesh%num_nodes))
               tempArray = merge(descript1%array,
     &              spread(descript1%alternate_value,1,kount(1)),nodecode.gt.0)
               iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id,tempArray)
               CALL check_err(iret)
            ELSE
               iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id,
     &             descript1%array)
               CALL check_err(iret)
            ENDIF
         ELSE                  ! PARALLEL
            iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id,
     &            descript1%array_g)
            CALL check_err(iret)
         ENDIF
      CASE(47,63,73,77,90,93,173,301,302,303,304,306,307,108) ! GML added 93 20210727
         IF (MNPROC.eq.1) THEN ! SERIAL
            IF ( descript1%ConsiderWetDry.eqv..TRUE.) THEN
               allocate(tempArray(1:dat%myMesh%num_nodes))
               tempArray = merge(descript1%array,
     &              spread(descript1%alternate_value,1,kount(1)),nodecode.gt.0) 

               iret = nf90_put_var(dat%ncid, dat%nodal_data_id,tempArray,start,kount)
            ELSE
               iret = nf90_put_var(dat%ncid, dat%nodal_data_id,
     &            descript1%array, start, kount)
            ENDIF
            CALL check_err(iret)
         ELSE

!c------------------------------
!            print*, "check ......", descript1%alternate_value, kount(1),
!     &               trim(descript1%file_name) ; 
!            DO ii = 1, kount(1)
!                  ! if ( isnan(descript1%array_g(ii))  ) then
!                  if ( ii == 2557 .OR. ii == 2573 ) then 
!                     print*, "check ..........", ii, descript1%array_g(ii), nodecode(ii) ; 
!                  endif
!            ENDDO
!            descript1%array_g(100) = (0.0/0.0);  ! Should get NaN
!            ! descript1%array_g(100) = descript1%array_g(100)*0 ; 
!            print*, "check ..........", 100, descript1%array_g(100), nodecode(100) ; 
!c------------------------------            

            iret = nf90_put_var(dat%ncid, dat%nodal_data_id,
     &            descript1%array_g, start, kount)


            CALL check_err(iret)
         ENDIF
      CASE(64,74,164,305)
         IF (MNPROC.eq.1) THEN ! SERIAL
            IF ( descript1%ConsiderWetDry.eqv..TRUE.) THEN
               allocate(tempArray(1:dat%myMesh%num_nodes))
               tempArray = merge(descript1%array,
     &              spread(descript1%alternate_value,1,kount(1)),nodecode.gt.0)
               iret = nf90_put_var(dat%ncid, dat%u_nodal_data_id,tempArray,start,kount)
               CALL check_err(iret)
               tempArray = merge(descript1%array2,
     &              spread(descript1%alternate_value,1,kount(1)),nodecode.gt.0)
               iret = nf90_put_var(dat%ncid, dat%v_nodal_data_id,tempArray,start,kount)
               CALL check_err(iret)
            ELSE
               iret = nf90_put_var(dat%ncid, dat%u_nodal_data_id,
     &             descript1%array, start, kount)
               CALL check_err(iret)
               iret = nf90_put_var(dat%ncid, dat%v_nodal_data_id,
     &             descript1%array2, start, kount)
               CALL check_err(iret)
            ENDIF
         ELSE                  ! PARALLEL
            iret = nf90_put_var(dat%ncid, dat%u_nodal_data_id,
     &            descript1%array_g, start, kount)
            CALL check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%v_nodal_data_id,
     &            descript1%array2_g, start, kount)
            CALL check_err(iret)
         ENDIF
      CASE(402,403)  ! integers at nodes (e.g., initiallydry.63, endrisinginun.63)
         IF (MNPROC.eq.1) THEN ! SERIAL
            IF ( descript1%ConsiderWetDry.eqv..TRUE.) THEN
               allocate(tempIntArray(1:dat%myMesh%num_nodes))
               tempIntArray = merge(descript1%iarray,
     &              spread(int(descript1%alternate_value),1,kount(1)),nodecode.gt.0)
               iret = nf90_put_var(dat%ncid, dat%nodal_data_id,tempIntArray,start,kount)
            ELSE
               iret = nf90_put_var(dat%ncid, dat%nodal_data_id,
     &            descript1%iarray, start, kount)
            ENDIF
            CALL check_err(iret)
         ELSE
            iret = nf90_put_var(dat%ncid, dat%nodal_data_id,
     &            descript1%iarray_g, start, kount)
            CALL check_err(iret)
         ENDIF
      CASE DEFAULT
            write(scratchMessage,
     &       '("No netCDF for files with unit number ",i0,".")') lun
            call allMessage(ERROR,scratchMessage)
      END SELECT
!
!     Close netCDF file
9999  iret = nf90_close(dat%ncid)
      CALL check_err(iret)

#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      end subroutine writeNodalData
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E
C      S E T   R E C O R D   C O U N T E R   A N D   S T O R E   T I M E
C-----------------------------------------------------------------------
C     jgf49.48.08 Compares the current ADCIRC simulation time with the
C     array of output times in the file, and if the simulation time
C     is before the end of the file, it sets the record counter to the
C     right place within the existing data. Data that occur after the
C     inserted data will remain, due to the inability of netcdf to
C     delete data from files.
C-----------------------------------------------------------------------
      SUBROUTINE setRecordCounterAndStoreTime(ncid, f, t)
      USE GLOBAL, ONLY : scratchMessage, scratchFormat
      IMPLICIT NONE
C
      INTEGER, intent(in) :: ncid
      type(fileData), intent(inout) :: f
      type(timeData), intent(inout) :: t
C
      REAL(8), allocatable :: storedTimes(:) ! array of time values in file
      LOGICAL :: timeFound ! true if current time is in array of stored times
      INTEGER :: ndim      ! number of dimensions in the netcdf file
      INTEGER :: nvar      ! number of variables in the netcdf file
      INTEGER :: natt      ! number of attributes in the netcdf file
      INTEGER :: counti(1), starti(1)
      INTEGER :: iret  ! success or failure of netcdf call
      INTEGER :: i     ! loop counter
C
      call setMessageSource("setRecordCounterAndStoreTime")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
!
!     Inquire time variable
      iret=nf90_inquire(ncid, ndim, nvar, natt, t%timenc_dim_id)
      CALL check_err(iret)
      iret=nf90_inquire_dimension(ncid, t%timenc_dim_id,
     &    len=f%record_counter)
      CALL check_err(iret)
      iret=nf90_inq_varid(ncid, "time", t%timenc_id)
      CALL check_err(iret)
C
C     Determine the relationship between the current simulation time
C     and the time array stored in the netcdf file. Set the record
C     counter based on this relationship.
      IF (f%record_counter.ne.0) THEN
         allocate(storedTimes(f%record_counter))
         iret=nf90_get_var(ncid, t%timenc_id, storedTimes)
         CALL check_err(iret)
         timeFound = .false.
         DO i=1,f%record_counter
            IF ((t%timenc(1).lt.storedTimes(i)).or.
     &         (abs(t%timenc(1)-storedTimes(i)).lt.1.0d-10)) THEN
               timeFound = .true.
               EXIT
            ENDIF
         ENDDO
         IF (timeFound.eqv..false.) THEN
            ! Increment the record counter so that we can store data at the
            ! next location in the netcdf file (i.e., all of the times
            ! in the netcdf file were found to be earlier than the current
            ! adcirc simulation time).
            f%record_counter = f%record_counter + 1
         ELSE
            ! jgf49.48.08: set the counter at the index that reflects the
            ! current time within the netcdf file (or is between two times
            ! found in the netcdf file).
            ! WARNING: all subsequent data will remain in the file, we
            ! are just overwriting it ... if we don't overwrite all of it,
            ! the pre-existing data will still be there, which is probably
            ! not what the user intended ... but apparently there is no
            ! way to delete data from netcdf files:
            ! http://www.unidata.ucar.edu/support/help/MailArchives/netcdf/msg02367.html
            scratchFormat =
     &      '("Overwriting pre-existing data in netcdf file ",a,'//
     &      '" for time=",f17.8,". ' //
     &      'Subsequent data in netcdf file remain unchanged.")'
            write(scratchMessage,scratchFormat)
     &          trim(f%FILENAME), t%timenc(1)
            CALL allMessage(INFO,scratchMessage)
            f%record_counter = i
         ENDIF
         deallocate(storedTimes)
      ELSE
         ! set the counter at 1 so we can record our first time value
         f%record_counter = 1
      ENDIF
!
!     Store simulation time in netcdf file
      starti(1)=f%record_counter
      counti(1)=t%timenc_len
      iret = nf90_put_var(ncid, t%timenc_id, t%timenc, starti, counti)
      CALL check_err(iret)

#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      end subroutine setRecordCounterAndStoreTime
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E    I N I T   N E T C D F   H O T S T A R T
C-----------------------------------------------------------------------
C     jgf49.35 Sets up netCDF variables and writes mesh data into netcdf
C     hotstart file.
C-----------------------------------------------------------------------
      SUBROUTINE initNetCDFHotstart(lun, Elev1Descript,
     &   Elev2Descript, VelDescript, CH1Descript, EtaDiscDescript,
     &   NodeCodeDescript, NOFFDescript, ncerror)
      USE GLOBAL, ONLY : OutputDataDescript_t, G, IM, IMHS, C3D, NE_G, NP_G
      USE MESH, ONLY : ICS
      IMPLICIT NONE
C
      INTEGER, intent(in) :: lun
      type(OutputDataDescript_t), intent(inout) :: Elev1Descript
      type(OutputDataDescript_t), intent(in) :: Elev2Descript
      type(OutputDataDescript_t), intent(in) :: VelDescript
      type(OutputDataDescript_t), intent(in) :: CH1Descript
      type(OutputDataDescript_t), intent(in) :: EtaDiscDescript
      type(OutputDataDescript_t), intent(in) :: NodeCodeDescript
      type(OutputDataDescript_t), intent(in) :: NOFFDescript
C
      LOGICAL ncerror
C
      INTEGER iret ! success or failure of the netcdf call

      INTEGER :: tempid
C
      call setMessageSource("initNetCDFHotstart")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ncerror = .false.
C
C     Point to the hotstart file we want to work on.
      IF (lun.eq.67) THEN
         hs => hs67
      ELSE
         hs => hs68
      ENDIF
      hs%myMesh => adcircMesh
      IF (hs%myMesh%initialized.eqv..false.) THEN
         hs%myMesh%num_nodes = NP_G
         hs%myMesh%num_elems = NE_G
         hs%myMesh%nface_len = 3
      ENDIF
      IF (hs%myTime%initialized.eqv..false.) THEN
         ALLOCATE(hs%myTime%timenc(hs%myTime%timenc_len))
         hs%myTime%initialized = .true.
      ENDIF
C
C     Initialize netCDF hotstart file, creating a new one
      Elev1Descript % lun = lun
      CALL createNetCDFOutputFile(hs%ncid, hs%myFile, hs%myTime,
     &                            Elev1Descript, ncerror)
      ! return an error flag to the calling routine if something went
      ! wrong when we tried to create the netcdf file
      if ( ncerror.eqv..true. ) then
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         return
      endif
      IF (hs%myMesh%initialized.eqv..false.) THEN
         hs%myMesh%num_nodes = NP_G
         hs%myMesh%num_elems = NE_G
         hs%myMesh%nface_len = 3
         CALL initNetCDFCoord(hs%myMesh)
      ENDIF
      CALL defineMeshVariables(hs%ncid, hs%myMesh, hs%myFile)
C
C     Z E T A 1
      hs%zeta1%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%zeta1%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'zeta1', NF90_DOUBLE,
     &   hs%zeta1%nodal_data_dims, hs%zeta1%nodal_data_id)
      CALL check_err(iret)

      iret = nf90_put_att(hs%ncid, hs%zeta1%nodal_data_id,
     &   'long_name', 'water surface elevation at previous time step')
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%zeta1%nodal_data_id,
     &   'standard_name',
     &      'water_surface_elevation_at_previous_time step')
      CALL check_err(iret)
      CALL putUnitsAttribute(hs%ncid, hs%zeta1%nodal_data_id,
     &        'm', 'ft')
      iret = nf90_put_att(hs%ncid, hs%zeta1%nodal_data_id,
     &         '_FillValue', doubleval)
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%zeta1%nodal_data_id,
     &      'positive', 'up')
      CALL check_err(iret)
C
C     Z E T A 2
      hs%zeta2%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%zeta2%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'zeta2', NF90_DOUBLE,
     &       hs%zeta2%nodal_data_dims, hs%zeta2%nodal_data_id)
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%zeta2%nodal_data_id,
     &      'long_name', 'water surface elevation at current time step')
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%zeta2%nodal_data_id,
     &      'standard_name',
     &      'water_surface_elevation_at_current_time_step')
      CALL check_err(iret)
      CALL putUnitsAttribute(hs%ncid, hs%zeta2%nodal_data_id,
     &        'm', 'ft')
      iret = nf90_put_att(hs%ncid, hs%zeta2%nodal_data_id,
     &        '_FillValue', doubleval)
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%zeta2%nodal_data_id,
     &      'positive', 'up')
      CALL check_err(iret)
C
C     Z E T A D
      hs%zetad%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%zetad%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'zetad', NF90_DOUBLE,
     &       hs%zetad%nodal_data_dims, hs%zetad%nodal_data_id)
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%zetad%nodal_data_id,
     &      'long_name',
     &      'water elevation at flux specified boundary')
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%zetad%nodal_data_id,
     &      'standard_name',
     &      'water_elevation_at_flux_specified_boundary')
      CALL check_err(iret)
      CALL putUnitsAttribute(hs%ncid, hs%zetad%nodal_data_id,
     &        'm', 'ft')
      iret = nf90_put_att(hs%ncid, hs%zetad%nodal_data_id,
     &        '_FillValue', doubleval)
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%zetad%nodal_data_id,'positive',
     &        'up')
      CALL check_err(iret)
C
C     U V E L
      hs%vel%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%vel%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'u-vel', NF90_DOUBLE,
     &       hs%vel%nodal_data_dims,  hs%vel%u_nodal_data_id)
      CALL check_err(iret)
      if (ics.ne.1) then
         iret = nf90_put_att(hs%ncid, hs%vel%u_nodal_data_id,
     &          'long_name', 'vertically averaged e/w velocity')
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%vel%u_nodal_data_id,
     &          'positive', 'east')
      else
         iret = nf90_put_att(hs%ncid, hs%vel%u_nodal_data_id,
     &       'long_name','vertically averaged velocity in x-direction')
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%vel%u_nodal_data_id,
     &          'positive', 'right')
      endif
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%vel%u_nodal_data_id,
     &       'standard_name', 'u_velocity')
      CALL check_err(iret)
      CALL putUnitsAttribute(hs%ncid, hs%vel%u_nodal_data_id,
     &        'm s-1', 'ft s-1')
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%vel%u_nodal_data_id,
     &        '_FillValue', doubleval)
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%vel%u_nodal_data_id,
     &      'dry_Value', doubleval)
      CALL check_err(iret)
C
C     V V E L
      iret = nf90_def_var(hs%ncid, 'v-vel', NF90_DOUBLE,
     &       hs%vel%nodal_data_dims, hs%vel%v_nodal_data_id)
      CALL check_err(iret)
      if (ics.ne.1) then
         iret = nf90_put_att(hs%ncid, hs%vel%v_nodal_data_id,
     &           'long_name', 'vertically averaged n/s velocity')
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%vel%v_nodal_data_id,
     &           'positive', 'north')
      else
         iret = nf90_put_att(hs%ncid, hs%vel%v_nodal_data_id,
     &       'long_name','vertically averaged velocity in y-direction')
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%vel%v_nodal_data_id,
     &  'positive', '90 degrees counterclockwise from x water velocity')
      endif
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%vel%v_nodal_data_id,
     &            'standard_name','v_velocity')
      CALL check_err(iret)
      CALL putUnitsAttribute(hs%ncid, hs%vel%v_nodal_data_id,
     &        'm s-1', 'ft s-1')
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%vel%v_nodal_data_id,
     &            '_FillValue', doubleval)
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%vel%v_nodal_data_id,
     &            'dry_Value', doubleval)
      CALL check_err(iret)
C
C     C H 1
      IF ((IM.EQ.10).OR.(IMHS.EQ.10)) THEN
         hs%ch1%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
         hs%ch1%nodal_data_dims(2) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'ch1', NF90_DOUBLE,
     &       hs%ch1%nodal_data_dims, hs%ch1%nodal_data_id)
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%ch1%nodal_data_id,
     &      'long_name', 'concentration')
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%ch1%nodal_data_id,
     &      'standard_name', 'concentration')
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%ch1%nodal_data_id,
     &        '_FillValue', doubleval)
         CALL check_err(iret)
      ENDIF
C
C     N O D E C O D E
      hs%nodecodenc%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%nodecodenc%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'nodecode', NF90_INT,
     &    hs%nodecodenc%nodal_data_dims, hs%nodecodenc%nodal_data_id)
      CALL check_err(iret)
      att_text = "wet or dry state of node where 1 indicates that the"
     &  //" node is wet and 0 indicates that the node is dry"
      iret = nf90_put_att(hs%ncid, hs%nodecodenc%nodal_data_id,
     &      'long_name', trim(att_text))
      CALL check_err(iret)
      att_text = "node_wet_or_dry"
      iret = nf90_put_att(hs%ncid, hs%nodecodenc%nodal_data_id,
     &      'standard_name', trim(att_text))
      CALL check_err(iret)
C
C     N O F F
      hs%noffnc%nodal_data_dims(1) = hs%myMesh%num_elems_dim_id
      hs%noffnc%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'noff', NF90_INT,
     &    hs%noffnc%nodal_data_dims, hs%noffnc%nodal_data_id)
      CALL check_err(iret)
      att_text = "wet or dry state of element where 1 indicates that"
     & //" the element is wet and 0 indicates that it is dry"
      iret = nf90_put_att(hs%ncid, hs%noffnc%nodal_data_id,
     &      'long_name', trim(att_text))
      CALL check_err(iret)
      att_text = "element_wet_or_dry"
      iret = nf90_put_att(hs%ncid, hs%noffnc%nodal_data_id,
     &      'standard_name', trim(att_text))
      CALL check_err(iret)
C
C     jgf50.44: Automatically turn on compression if we are using the
C     netcdf4 file format.
#ifdef NETCDF_CAN_DEFLATE
      IF ( (Elev1Descript%specifier.eq.5).or.
     &     (Elev1Descript%specifier.eq.567).or.
     &     (Elev1Descript%specifier.eq.568) ) THEN
         iret = nf_def_var_deflate(hs%ncid, hs%zeta1%nodal_data_id,
     &             1, 1, 2)
         CALL check_err(iret)
         iret = nf_def_var_deflate(hs%ncid, hs%zeta2%nodal_data_id,
     &             1, 1, 2)
         CALL check_err(iret)
         iret = nf_def_var_deflate(hs%ncid, hs%zetad%nodal_data_id,
     &             1, 1, 2)
         CALL check_err(iret)
         iret = nf_def_var_deflate(hs%ncid, hs%vel%u_nodal_data_id,
     &             1, 1, 2)
         CALL check_err(iret)
         iret = nf_def_var_deflate(hs%ncid, hs%vel%v_nodal_data_id,
     &             1, 1, 2)
         CALL check_err(iret)
         IF ((IM.EQ.10).OR.(IMHS.EQ.10)) THEN
            iret = nf_def_var_deflate(hs%ncid, hs%ch1%nodal_data_id,
     &             1, 1, 2)
            CALL check_err(iret)
         ENDIF
         iret = nf_def_var_deflate(hs%ncid, hs%nodecodenc%nodal_data_id,
     &             1, 1, 2)
         CALL check_err(iret)
         iret = nf_def_var_deflate(hs%ncid, hs%noffnc%nodal_data_id,
     &             1, 1, 2)
         CALL check_err(iret)
      ENDIF
#endif
!
!     Define hotstart parameters
      CALL check_err(iret)
      iret = nf90_def_var(hs%ncid, 'imhs', NF90_INT, varid=tempid)
      CALL check_err(iret)
      att_text = 'model_type'
      iret = nf90_put_att(hs%ncid, tempid, 'long_name', 'model_type')
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &       'model_type')
       CALL check_err(iret)
C
      iret = nf90_def_var(hs%ncid, 'iths', NF90_INT, varid=tempid)
      CALL check_err(iret)
      att_text =
     &     'model time step number since the beginning of the model run'
      iret = nf90_put_att(hs%ncid, tempid, 'long_name',trim(att_text))
      CALL check_err(iret)
      att_text = 'model_time_step'
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &       trim(att_text))
      CALL check_err(iret)
C
      iret = nf90_def_var(hs%ncid, 'iestp', NF90_INT, varid=tempid)
      CALL check_err(iret)
      att_text = 'line number (for ASCII output) or record number'
     &  //' (for binary output) of the most recent entry in the'
     &  //' elevation time series at specified elevation recording'
     &  //' stations output file'
      iret = nf90_put_att(hs%ncid, tempid, 'long_name',trim(att_text))
      CALL check_err(iret)
      att_text = 'line/record_number_of_last_entry_in_elev_rec_stations'
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &       trim(att_text))
      CALL check_err(iret)
c
      iret = nf90_def_var(hs%ncid, 'nscoue', NF90_INT, varid=tempid)
      CALL check_err(iret)
      att_text = 'time step counter to determine when the'
     & //' next entry will be written to the elevation time series at'
     & //' specified elevation recording Stations output file'
      iret = nf90_put_att(hs%ncid, tempid, 'long_name', trim(att_text))

      CALL check_err(iret)
      att_text = 'time_step_counter_for_next_entry_elev_rec_stations'
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &       trim(att_text))
      CALL check_err(iret)
C
      iret = nf90_def_var(hs%ncid, 'ivstp', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'long_name',
     &                   ' line number (for ASCII output) or record'
     & //' number (for binary output) of the most recent entry in the'
     & //' depth-averaged velocity time series at specified velocity'
     & //' recording stations output file')
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &               'record_number_of_last_entry_in_vel_rec_stations')
      CALL check_err(iret)
C
      iret = nf90_def_var(hs%ncid, 'nscouv', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'long_name',
     &                   'time step counter to determine when the'
     & //' next entry will be written to the depth-averaged velocity'
     & //' time series at specified velocity recording stations output'
     & //' file.')
      CALL check_err(iret)

      iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &           'time_step_counter_for_next_entry_vel_rec_stations')
      CALL check_err(iret)
C
      IF ((IM.EQ.10).OR.(IMHS.EQ.10)) THEN
         iret = nf90_def_var(hs%ncid, 'icstp', NF90_INT, varid=tempid)
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, tempid, 'long_name',
     &                       'line number (for ASCII output) or record'
     & //' number (for binary output) of the most recent entry in the'
     & //' scalar concentration time series at specified concentration'
     & //' recording stations output file')
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &               'record_number_of_last_entry_in_conc_rec_stations')
         CALL check_err(iret)
C
         iret = nf90_def_var(hs%ncid, 'nscouc', NF90_INT, varid=tempid)
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, tempid, 'long_name',
     &                'time step counter to determine when the next '
     & //'entry will be written to the scalar concentration time series'
     & //' at specified concentration recording stations output file')
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &           'time_step_counter_for_next_entry_conc_rec_stations')
         CALL check_err(iret)
      ENDIF
C
      iret = nf90_def_var(hs%ncid, 'ipstp', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'long_name',
     &                       'line number (for ASCII output) or record'
     & //' number (for binary output) of the most recent entry in the'
     & //' atmospheric pressure time series at specified meteorological'
     & //' recording stations')
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &       'record_number_of_last_entry_of_atm_press_at_rec_stations')
      CALL check_err(iret)
C
      iret = nf90_def_var(hs%ncid, 'iwstp', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'long_name',
     &                      'line number (for ASCII output) or record'
     & //' number (for binary output) of the most recent entry in the'
     & //' wind velocity time series at specified meteorological'
     & //' recording stations')
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &       'record_number_of_last_entry_of_wind_vel_at_rec_stations')
      CALL check_err(iret)
C
      iret = nf90_def_var(hs%ncid, 'nscoum', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'long_name',
     &                       'time step counter to determine when the'
     & //' next entry will be written to the atmospheric pressure time'
     & //' series at specified meteorological recording stations and'
     & //' wind velocity time series at specified meteorological'
     & //' recording stations output files')
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &   'time_step_counter_of_atm_press_and_wind_vel_at_rec_stations')
      CALL check_err(iret)
C
      iret = nf90_def_var(hs%ncid, 'igep', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'long_name',
     &                       'line number (for ASCII output) or record'
     & //' number (for binary output) of the most recent entry in the'
     & //' elevation time series at all nodes in the model grid output'
     & //' file')
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &       'record_number_of_last_entry_of_elev_at_model_nodes')
      CALL check_err(iret)
C
      iret = nf90_def_var(hs%ncid, 'nscouge', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'long_name',
     &                      'time step counter to determine when the'
     & //' next entry will be written to the elevation time series at'
     & //' all nodes in the model grid output file')
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &   'time_step_counter_of_elev_at_model_nodes')
      CALL check_err(iret)
C
      iret = nf90_def_var(hs%ncid, 'igvp', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'long_name',
     &                       'line number (for ASCII output) or record'
     & //' number (for binary output) of the most recent entry in the'
     & //' depth-averaged velocity time series at all nodes in the'
     & //' model grid output file')
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &       'record_number_of_last_entry_of_vel_at_model_nodes')
      CALL check_err(iret)
C
      iret = nf90_def_var(hs%ncid, 'nscougv', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'long_name',
     &                      'time step counter to determine when the'
     & //' next entry will be written to the depth-averaged velocity'
     & //' time series at all nodes in the model grid output file')
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &   'time_step_counter_of_vel_at_model_nodes')
      CALL check_err(iret)
C
      IF ((IM.EQ.10).OR.(IMHS.EQ.10)) THEN
         iret = nf90_def_var(hs%ncid, 'igcp', NF90_INT, varid=tempid)
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, tempid, 'long_name',
     &                     'line number (for ASCII output) or record'
     & //' number (for binary output) of the most recent entry in the'
     & //' scalar concentration time series at all nodes in the model'
     & //' grid output file')
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &       'record_number_of_last_entry_of_conc_at_model_nodes')
         CALL check_err(iret)
C
         iret = nf90_def_var(hs%ncid, 'nscougc', NF90_INT, varid=tempid)
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, tempid, 'long_name',
     &                     'time step counter to determine when the'
     & //' next entry will be written to the scalar concentration time'
     & //' series at All Nodes in the model grid output file')
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &   'time_step_counter_of_conc_at_model_nodes')
         CALL check_err(iret)
      ENDIF
C
      iret = nf90_def_var(hs%ncid, 'igpp', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'long_name',
     &                       'line number (for ASCII output) or record'
     & //' number (for binary output) of the most recent entry in the'
     & //' atmospheric pressure time series at all nodes in the model'
     & //' grid output file')
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &       'record_number_of_last_entry_of_atm_press_at_model_nodes')
      CALL check_err(iret)
C
      iret = nf90_def_var(hs%ncid, 'igwp', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'long_name',
     &                      'line number (for ASCII output) or record'
     & //' number (for binary output) of the most recent entry in the'
     & //' wind Stress or velocity time series at all nodes in the'
     & //' model grid output file')
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &       'record_number_of_last_entry_of_wind_vel_at_model_nodes')
      CALL check_err(iret)
C
      iret = nf90_def_var(hs%ncid, 'nscougw', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'long_name',
     &                       'time step counter to determine when the'
     & //' next entry will be written to the atmospheric pressure time'
     & //' series at all nodes in the model grid and wind stress or'
     & //' velocity time series at all nodes in the model grid output'
     & //' files')
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name',
     &   'time_step_counter_of_atm_press_and_wind_vel_at_model_nodes')
      CALL check_err(iret)
C
C     define time attributes
      CALL defineTimeAttributes(hs%ncid, hs%myTime)
C
C     define metadata and selected fort.15 parameters in netcdf file
      CALL defineMetaData(hs%ncid)
C
C     Leave define mode
      iret = nf90_enddef(hs%ncid)
      CALL check_err(iret)
C
C     write mesh to netcdf file
      CALL putMeshVariables(hs%ncid,hs%myMesh)
C
C     now close the initialized netcdf file
      iret = nf90_close(hs%ncid)
      CALL check_err(iret)
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE initNetCDFHotstart
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E
C        I N I T   N E T C D F   H O T S T A R T   H A R M O N I C
C-----------------------------------------------------------------------
C     jgf49.35 Sets up netCDF variables for hotstarting harmonic analysis.
C-----------------------------------------------------------------------
      SUBROUTINE initNetCDFHotstartHarmonic(lun, GLOELVDescript,
     &   STAELVDescript, GLOULVDescript, GLOVLVDescript, STAULVDescript,
     &   STAVLVDescript, err)
      USE SIZES, ONLY : MNHARF
      USE HARM, ONLY : NHASE, NHASV, NHAGE, NHAGV
      USE GLOBAL, ONLY : OutputDataDescript_t
      IMPLICIT NONE
C
      INTEGER, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: GLOELVDescript
      type(OutputDataDescript_t), intent(in) :: STAELVDescript
      type(OutputDataDescript_t), intent(in) :: GLOULVDescript
      type(OutputDataDescript_t), intent(in) :: GLOVLVDescript
      type(OutputDataDescript_t), intent(in) :: STAULVDescript
      type(OutputDataDescript_t), intent(in) :: STAVLVDescript
C
      LOGICAL err
C
      INTEGER iret ! success or failure of the netcdf call
      CHARACTER(1024) :: att_text
      INTEGER :: varid
C
      call setMessageSource("initNetCDFHotstartHarmonic")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      err = .false.
C
C     Point to the hotstart file we want to work on.
      IF (lun.eq.67) THEN
         hs => hs67
      ELSE
         hs => hs68
      ENDIF
C
C     Open existing NetCDF file
      iret = nf90_open(hs%myFile%FILENAME, NF_WRITE, hs%ncid)
      CALL check_err(iret)
C
C     Enter "redefine" mode
      iret = NF90_REDEF(hs%ncid)
      CALL check_err(iret)
!
!     Inquire variables
      CALL getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)
C
C     Define dimensions
      iret = nf90_def_dim(hs%ncid,'mnharf',MNHARF,
     &       hs%mnharf_dim_id)
      CALL check_err(iret)
      iret = nf90_def_dim(hs%ncid,'mnharfx2',(MNHARF*2),
     &       hs%load_vector_dim_id)
      CALL check_err(iret)
C
C     Create station dimension and station name dimension
      iret = nf90_def_dim(hs%ncid, 'namefrlen', hs%namefr_len,
     &   hs%namefr_len_dim_id)
      CALL check_err(iret)
!
!     Define harmonic analysis frequency names array
      hs%namefr_dims(1) = hs%namefr_len_dim_id
      hs%namefr_dims(2) = hs%mnharf_dim_id
      iret = nf90_def_var(hs%ncid, 'namefr', NF_CHAR,
     &   hs%namefr_dims, hs%namefr_id)
      CALL check_err(iret)
C
C     harmonic analysis components
      hs%component_dims(1) = hs%mnharf_dim_id
      iret = nf90_def_var(hs%ncid, 'hafreq', NF90_DOUBLE,
     &       hs%component_dims, hs%hafreq_id)
      CALL check_err(iret)
      att_text = "frequencies (rad/s) of harmonic analysis constituents"
      iret = nf90_put_att(hs%ncid, hs%hafreq_id,
     &   'long_name', trim(att_text))
      CALL check_err(iret)
      att_text = "harmonic analysis frequencies (rad/s)"
      iret = nf90_put_att(hs%ncid, hs%hafreq_id,
     &   'standard_name', trim(att_text))
C
      iret = nf90_def_var(hs%ncid, 'haff', NF90_DOUBLE,
     &       hs%component_dims, hs%haff_id)
      CALL check_err(iret)
      att_text = "nodal factors of harmonic analysis constituents"
      iret = nf90_put_att(hs%ncid, hs%haff_id,
     &   'long_name', trim(att_text))
      CALL check_err(iret)
      att_text = "harmonic analysis nodal factors"
      iret = nf90_put_att(hs%ncid, hs%haff_id,
     &   'standard_name', trim(att_text))
C
      iret = nf90_def_var(hs%ncid, 'haface', NF90_DOUBLE,
     &        hs%component_dims, hs%haface_id)
      CALL check_err(iret)
      att_text =
     &  "equilibrium arguments (degrees) of harmonic analysis "
     & //"constituents"
      iret = nf90_put_att(hs%ncid, hs%haface_id,
     &   'long_name',trim(att_text))
      CALL check_err(iret)
      att_text = "equilibrium arguments (degrees)"
      iret = nf90_put_att(hs%ncid, hs%haface_id,
     &   'standard_name', trim(att_text))
C
      hs%ha_dims(1) = hs%load_vector_dim_id
      hs%ha_dims(2) = hs%load_vector_dim_id
      iret = nf90_def_var(hs%ncid, 'ha',NF90_DOUBLE,hs%ha_dims,hs%ha_id)
      CALL check_err(iret)
      att_text = "left hand side matrix for harmonic analysis"
      iret = nf90_put_att(hs%ncid, hs%ha_id,'long_name', trim(att_text))
      CALL check_err(iret)
      att_text = "LHS for harmonic analysis"
      iret = nf90_put_att(hs%ncid, hs%ha_id,
     &   'standard_name', trim(att_text))
C
C     global elevation load vector
      IF (NHAGE.ne.0) THEN
         hs%gloelv%nodal_data_dims_3D(1) = hs%load_vector_dim_id
         hs%gloelv%nodal_data_dims_3D(2) = hs%myMesh%num_nodes_dim_id
         hs%gloelv%nodal_data_dims_3D(3) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'gloelv', NF90_DOUBLE,
     &       hs%gloelv%nodal_data_dims_3D, hs%gloelv%nodal_data_id)
         CALL check_err(iret)
         att_text = "full domain elevation load vector at each node"
         iret = nf90_put_att(hs%ncid, hs%gloelv%nodal_data_id,
     &      'long_name', trim(att_text))
         CALL check_err(iret)
         att_text = "full domain elevation load vector"
         iret = nf90_put_att(hs%ncid, hs%gloelv%nodal_data_id,
     &      'standard_name', trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%gloelv%nodal_data_id,
     &         '_FillValue', doubleval)
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%gloelv%nodal_data_id,
     &        'positive', 'up')
         CALL check_err(iret)
      ENDIF
C
C     global velocity load vectors
      IF (NHAGV.ne.0) THEN
         hs%glovellv%nodal_data_dims_3D(1) = hs%load_vector_dim_id
         hs%glovellv%nodal_data_dims_3D(2) = hs%myMesh%num_nodes_dim_id
         hs%glovellv%nodal_data_dims_3D(3) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'gloulv', NF90_DOUBLE,
     &       hs%glovellv%nodal_data_dims_3D, hs%glovellv%u_nodal_data_id)
         CALL check_err(iret)
         att_text = "full domain u velocity load vector at each node"
         iret = nf90_put_att(hs%ncid, hs%glovellv%u_nodal_data_id,
     &      'long_name', trim(att_text))
         CALL check_err(iret)
         att_text = "full domain u velocity load vector"
         iret = nf90_put_att(hs%ncid, hs%glovellv%u_nodal_data_id,
     &      'standard_name', trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%glovellv%u_nodal_data_id,
     &         '_FillValue', doubleval)
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%glovellv%u_nodal_data_id,
     &          'positive', 'east')
         CALL check_err(iret)
C
         iret = nf90_def_var(hs%ncid, 'glovlv', NF90_DOUBLE,
     &       hs%glovellv%nodal_data_dims_3D,hs%glovellv%v_nodal_data_id)
         CALL check_err(iret)
         att_text = "full domain v velocity load vector at each node"
         iret = nf90_put_att(hs%ncid, hs%glovellv%v_nodal_data_id,
     &      'long_name', trim(att_text))
         CALL check_err(iret)
         att_text = "full domain v velocity load vector"
         iret = nf90_put_att(hs%ncid, hs%glovellv%v_nodal_data_id,
     &      'standard_name', trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%glovellv%v_nodal_data_id,
     &         '_FillValue', doubleval)
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%glovellv%v_nodal_data_id,
     &           'positive', 'north')
         CALL check_err(iret)
      ENDIF
C
C     station elevation load vector
      IF (NHASE.ne.0) THEN
         hs%staelv%num_stations = STAELVDescript%num_fd_records
         iret = nf90_def_dim(hs%ncid, 'elevstation',
     &       hs%staelv%num_stations, hs%staelv%num_sta_dim_id)
         CALL check_err(iret)
         hs%staelv%station_data_dims_3D(1) = hs%load_vector_dim_id
         hs%staelv%station_data_dims_3D(2) = hs%staelv%num_sta_dim_id
         hs%staelv%station_data_dims_3D(3) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'staelv', NF90_DOUBLE,
     &       hs%staelv%station_data_dims_3D, hs%staelv%station_data_id)
         CALL check_err(iret)
         att_text = "elevation load vector at each elevation station"
         iret = nf90_put_att(hs%ncid, hs%staelv%station_data_id,
     &      'long_name',trim(att_text))
         CALL check_err(iret)
         att_text = "station elevation load vector"
         iret = nf90_put_att(hs%ncid, hs%staelv%station_data_id,
     &         'standard_name',trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%staelv%station_data_id,
     &         '_FillValue', doubleval)
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%staelv%station_data_id,
     &      'positive', 'up')
         CALL check_err(iret)
      ENDIF
C
C     station velocity load vectors
      IF (NHASV.ne.0) THEN
         hs%stavellv%num_stations = STAULVDescript%num_fd_records
         iret = nf90_def_dim(hs%ncid, 'velstation',
     &       hs%stavellv%num_stations, hs%stavellv%num_sta_dim_id)
         CALL check_err(iret)
C        define dimension
         hs%stavellv%station_data_dims_3D(1) = hs%load_vector_dim_id
         hs%stavellv%station_data_dims_3D(2) =hs%stavellv%num_sta_dim_id
         hs%stavellv%station_data_dims_3D(3) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'staulv', NF90_DOUBLE,
     &       hs%stavellv%station_data_dims_3D,
     &       hs%stavellv%u_station_data_id)
         CALL check_err(iret)
         att_text = "u velocity load vector at each velocity station"
         iret = nf90_put_att(hs%ncid, hs%stavellv%u_station_data_id,
     &      'long_name',trim(att_text))
         CALL check_err(iret)
         att_text = "station u velocity load vector"
         iret = nf90_put_att(hs%ncid, hs%stavellv%u_station_data_id,
     &      'standard_name', trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid,
     &         hs%stavellv%u_station_data_id,
     &         '_FillValue', doubleval)
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%stavellv%u_station_data_id,
     &          'positive', 'east')
         CALL check_err(iret)
C
         iret = nf90_def_var(hs%ncid, 'stavlv', NF90_DOUBLE,
     &       hs%stavellv%station_data_dims_3D,
     &       hs%stavellv%v_station_data_id)
         CALL check_err(iret)
         att_text = "v velocity load vector at each velocity station"
         iret = nf90_put_att(hs%ncid, hs%stavellv%v_station_data_id,
     &      'long_name', trim(att_text))
         CALL check_err(iret)
         att_text = "station v velocity load vector"
         iret = nf90_put_att(hs%ncid, hs%stavellv%v_station_data_id,
     &      'standard_name', trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%stavellv%v_station_data_id,
     &         '_FillValue', doubleval)
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%stavellv%v_station_data_id,
     &           'positive','north')
         CALL check_err(iret)
      ENDIF
C
C     jgf50.44: Automatically turn on compression if we are using the
C     netcdf4 file format.
#ifdef NETCDF_CAN_DEFLATE
      IF ( (GLOELVDescript%specifier.eq.5).or.
     &     (GLOELVDescript%specifier.eq.567).or.
     &     (GLOELVDescript%specifier.eq.568) ) THEN
         iret = nf_def_var_deflate(hs%ncid, hs%ha_id,
     &             1, 1, 2)
         CALL check_err(iret)
         IF (NHAGE.ne.0) THEN
            iret = nf90_def_var_deflate(hs%ncid,
     &            hs%gloelv%nodal_data_id, 1, 1, 2)
            CALL check_err(iret)
         ENDIF
         IF (NHAGV.ne.0) THEN
            iret = nf_def_var_deflate(hs%ncid,
     &            hs%glovellv%u_nodal_data_id, 1, 1, 2)
            CALL check_err(iret)
            iret = nf_def_var_deflate(hs%ncid,
     &            hs%glovellv%v_nodal_data_id, 1, 1, 2)
            CALL check_err(iret)
         ENDIF
         IF (NHASE.ne.0) THEN
            iret = nf_def_var_deflate(hs%ncid,
     &            hs%staelv%station_data_id, 1, 1, 2)
            CALL check_err(iret)
         ENDIF
         IF (NHASV.ne.0) THEN
            iret = nf_def_var_deflate(hs%ncid,
     &            hs%stavellv%u_station_data_id, 1, 1, 2)
            CALL check_err(iret)
            iret = nf_def_var_deflate(hs%ncid,
     &            hs%stavellv%v_station_data_id, 1, 1, 2)
            CALL check_err(iret)
         ENDIF
      ENDIF
#endif
!
!     Define harmonic analysis parameters
      call defineParameterWithText(hs%ncid, 'icha', NF90_INT,
     &     "harmonic analysis spool counter",
     &     "harmonic analysis spool counter")
C
      call defineParameterWithText(hs%ncid, 'nz', NF90_INT,
     &     "set to 0 if a steady harmonic component is included",
     &     "indicator of steady harmonic component")
C
      call defineParameterWithText(hs%ncid, 'nf', NF90_INT,
     &  "set to 1 if a steady harmonic component is included",
     &  "steady harmonic component number")
C
      call defineParameterWithText(hs%ncid, 'mm', NF90_INT,
     & "2x the number of harmonic frequencies plus any steady component",
     & "2x harmonic frequencies [plus 1]")
C
      call defineParameterWithText(hs%ncid, 'nstae', NF90_INT,
     &  "number of elevation recording stations for harmonic analysis",
     &  "number of elevation recording stations")
C
      call defineParameterWithText(hs%ncid, 'nstav', NF90_INT,
     &  "number of velocity recording stations for harmonic analysis",
     &  "number of velocity recording stations")
C
       call defineParameterWithText(hs%ncid, 'nhase', NF90_INT,
     &  "indicator for perfomance and formatting of harmonic analysis "
     & //"of elevation station data",
     &  "elevation station harmonic analysis indicator")
C
       call defineParameterWithText(hs%ncid, 'nhasv', NF90_INT,
     &  "indicator for perfomance and formatting of harmonic analysis "
     & //"of velocity station data",
     &  "velocity station harmonic analysis indicator")
C
       call defineParameterWithText(hs%ncid, 'nhage', NF90_INT,
     &  "indicator for perfomance and formatting of harmonic analysis "
     &  //"of full domain elevation data (at every node)",
     &  "full domain elevation harmonic analysis indicator")
C
       call defineParameterWithText(hs%ncid, 'nhagv', NF90_INT,
     &  "indicator for perfomance and formatting of harmonic analysis "
     &  //"of full domain velocity data (at every node)",
     &  "full domain velocity harmonic analysis indicator")
C
       call defineParameterWithText(hs%ncid, 'icall', NF90_INT,
     &  "number of subroutine calls to update load vectors and left "
     & //"matrix for harmonic analysis",
     &  "number of calls to update harmonic analysis")
C
       call defineParameterWithText(hs%ncid, 'nfreq', NF90_INT,
     &  "number of frequencies under consideration in harmonic analysis"
     & //" not including a steady component, if any",
     &  "number of frequencies for harmonic analysis")
C
       call defineParameterWithText(hs%ncid, 'timeud', NF90_DOUBLE,
     &  "ADCIRC time at the most recent update of the load vectors for "
     &  //"harmonic analysis",
     &  "update time for load vectors")
C
       call defineParameterWithText(hs%ncid, 'itud', NF90_INT,
     &  "ADCIRC time step at the most recent update of the load vectors"
     &  //" for harmonic analysis",
     &  "update time step for load vectors")
C
C
C     Leave define mode
      iret = nf90_enddef(hs%ncid)
      CALL check_err(iret)
C
C     now close the initialized netcdf file
      iret = nf90_close(hs%ncid)
      CALL check_err(iret)
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE initNetCDFHotstartHarmonic
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E
C        I N I T   N E T C D F   H O T S T A R T   H A R M O N I C
C                                  M E A N S   V A R I A N C E S
C-----------------------------------------------------------------------
C     jgf49.43.14 Sets up netCDF variables for hotstarting harmonic analysis
C     means and variances calculations.
C-----------------------------------------------------------------------
      SUBROUTINE initNetCDFHotstartHarmonicMeansVariances(lun,
     &   ELAVDescript, ELVADescript, XVELAVDescript, YVELAVDescript,
     &   XVELVADescript, YVELVADescript, reterror)
      USE SIZES, ONLY : MNHARF
      USE HARM, ONLY : NHASE, NHASV, NHAGE, NHAGV
      USE GLOBAL, ONLY : OutputDataDescript_t
      IMPLICIT NONE
C
      INTEGER, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: ELAVDescript
      type(OutputDataDescript_t), intent(in) :: ELVADescript
      type(OutputDataDescript_t), intent(in) :: XVELAVDescript
      type(OutputDataDescript_t), intent(in) :: YVELAVDescript
      type(OutputDataDescript_t), intent(in) :: XVELVADescript
      type(OutputDataDescript_t), intent(in) :: YVELVADescript
C
      INTEGER iret ! success or failure of the netcdf call
      CHARACTER(1024) :: att_text
      INTEGER :: varid
      LOGICAL :: reterror
C
      call setMessageSource("initNetCDFHotstartHarmonicMeansVariances")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      reterror = .false.
C
C     Point to the hotstart file we want to work on. Memory allocation
C     was already done for means and variances by initNetCDFHotstartHarmonic.
      IF (lun.eq.67) THEN
         hs => hs67
      ELSE
         hs => hs68
      ENDIF
C
C     Open existing NetCDF file
      iret = nf90_open(hs%myFile%FILENAME, NF_WRITE, hs%ncid)
      CALL check_err(iret)
C
C     Enter "redefine" mode
      iret = nf90_redef(hs%ncid)
      CALL check_err(iret)
!
!     Inquire variables
      CALL getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)
C
C     elevation
      IF (NHAGE.ne.0) THEN
         ! ELAV
         hs%elav%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
         hs%elav%nodal_data_dims(2) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'elav', NF90_DOUBLE,
     &       hs%elav%nodal_data_dims, hs%elav%nodal_data_id)
         CALL check_err(iret)
         att_text = "sum of elevations computed by ADCIRC, at every "
     &   //"node in the model grid, over all time steps since harmonic "
     &   //"analysis means and variance checking has begun"
         iret = nf90_put_att(hs%ncid, hs%elav%nodal_data_id,
     &      'long_name',trim(att_text))
         CALL check_err(iret)
         att_text = "sum of elevations"
         iret = nf90_put_att(hs%ncid, hs%elav%nodal_data_id,
     &      'standard_name', trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%elav%nodal_data_id,
     &         '_FillValue', doubleval)
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%elav%nodal_data_id,
     &        'positive', 'up')
         CALL check_err(iret)
         ! ELVA
         hs%elva%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
         hs%elva%nodal_data_dims(2) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'elva', NF90_DOUBLE,
     &       hs%elva%nodal_data_dims, hs%elva%nodal_data_id)
         CALL check_err(iret)
         att_text = "sum of squares of elevations computed by ADCIRC, "
     & //"at every node in the model grid, over all time steps since "
     & //"harmonic analysis means and variance checking has begun"
         iret = nf90_put_att(hs%ncid, hs%elva%nodal_data_id,
     &      'long_name', trim(att_text))
         CALL check_err(iret)
         att_text = "sum of squares of elevations"
         iret = nf90_put_att(hs%ncid, hs%elva%nodal_data_id,
     &      'standard_name', trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%elva%nodal_data_id,
     &         '_FillValue', doubleval)
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%elva%nodal_data_id,
     &        'positive', 'up')
         CALL check_err(iret)
      ENDIF
C
C     global velocity load vectors
      IF (NHAGV.ne.0) THEN
         ! XVELAV
         hs%xvelav%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
         hs%xvelav%nodal_data_dims(2) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'xvelav', NF90_DOUBLE,
     &       hs%xvelav%nodal_data_dims,
     &       hs%xvelav%nodal_data_id)
         CALL check_err(iret)
         att_text = "sum of depth-averaged u velocities computed by "
     & //"ADCIRC, at every node in the model grid, over all time steps "
     & //"since harmonic analysis means and variance checking has begun"
         iret = nf90_put_att(hs%ncid, hs%xvelav%nodal_data_id,
     &      'long_name', trim(att_text))
         CALL check_err(iret)
         att_text = "sum of depth averaged u velocities"
         iret = nf90_put_att(hs%ncid, hs%xvelav%nodal_data_id,
     &      'standard_name', trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%xvelav%nodal_data_id,
     &         '_FillValue', doubleval)
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%xvelav%nodal_data_id,
     &        'positive', 'east')
         CALL check_err(iret)
         ! YVELAV
         hs%yvelav%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
         hs%yvelav%nodal_data_dims(2) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'yvelav', NF90_DOUBLE,
     &       hs%yvelav%nodal_data_dims, hs%yvelav%nodal_data_id)
         CALL check_err(iret)
         att_text = "sum of depth-averaged v velocities computed by "
     & //"ADCIRC, at every node in the model grid, over all time steps "
     & //"since harmonic analysis means and variance checking has begun"
         iret = nf90_put_att(hs%ncid, hs%yvelav%nodal_data_id,
     &      'long_name', trim(att_text))
         CALL check_err(iret)
         att_text = "sum of depth averaged v velocities"
         iret = nf90_put_att(hs%ncid, hs%yvelav%nodal_data_id,
     &      'standard_name', trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%yvelav%nodal_data_id,
     &         '_FillValue', doubleval)
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%yvelav%nodal_data_id,
     &        'positive', 'north')
         CALL check_err(iret)

         ! XVELVA
         hs%xvelva%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
         hs%xvelva%nodal_data_dims(2) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'xvelva', NF90_DOUBLE,
     &       hs%xvelva%nodal_data_dims,
     &       hs%xvelva%nodal_data_id)
         CALL check_err(iret)
         att_text = "sum of squares of depth averaged u velocities "
     & //"computed by ADCIRC, at every node in the model grid, over all"
     & //" time steps since harmonic analysis means and variance "
     & //"checking has begun"
         iret = nf90_put_att(hs%ncid, hs%xvelva%nodal_data_id,
     &      'long_name',trim(att_text))
         CALL check_err(iret)
         att_text = "sum of squares of depth averaged u velocities"
         iret = nf90_put_att(hs%ncid, hs%xvelva%nodal_data_id,
     &      'standard_name', trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%xvelva%nodal_data_id,
     &         '_FillValue',  doubleval)
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%xvelva%nodal_data_id,
     &        'positive',  'east')
         CALL check_err(iret)
         ! YVELVA
         hs%yvelva%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
         hs%yvelva%nodal_data_dims(2) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'yvelva', NF90_DOUBLE,
     &       hs%yvelva%nodal_data_dims, hs%yvelva%nodal_data_id)
         CALL check_err(iret)
         att_text = "sum of squares of depth averaged v velocities "
     &  //"computed by ADCIRC, at every node in the model grid, over "
     &  //"all time steps since harmonic analysis means and variance "
     &  //"checking has begun"
         iret = nf90_put_att(hs%ncid, hs%yvelva%nodal_data_id,
     &      'long_name', trim(att_text))
         CALL check_err(iret)
         att_text = "sum of squares of depth averaged v velocities"
         iret = nf90_put_att(hs%ncid, hs%yvelva%nodal_data_id,
     &      'standard_name', trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%yvelva%nodal_data_id,
     &         '_FillValue', doubleval)
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%yvelva%nodal_data_id,
     &        'positive', 'up')
         CALL check_err(iret)
      ENDIF

C
C     jgf50.44: Automatically turn on compression if we are using the
C     netcdf4 file format.
#ifdef NETCDF_CAN_DEFLATE
      IF ( (ELAVDescript%specifier.eq.5).or.
     &     (ELAVDescript%specifier.eq.567).or.
     &     (ELAVDescript%specifier.eq.568) ) THEN
         IF (NHAGE.ne.0) THEN
            iret = nf_def_var_deflate(hs%ncid,
     &            hs%elav%nodal_data_id, 1, 1, 2)
            CALL check_err(iret)
            iret = nf_def_var_deflate(hs%ncid,
     &            hs%elva%nodal_data_id, 1, 1, 2)
            CALL check_err(iret)
         ENDIF
         IF (NHAGV.ne.0) THEN
            iret = nf_def_var_deflate(hs%ncid,
     &            hs%xvelav%nodal_data_id, 1, 1, 2)
            CALL check_err(iret)
            iret = nf_def_var_deflate(hs%ncid,
     &            hs%xvelva%nodal_data_id, 1, 1, 2)
            CALL check_err(iret)
            iret = nf_def_var_deflate(hs%ncid,
     &            hs%yvelav%nodal_data_id, 1, 1, 2)
            CALL check_err(iret)
            iret = nf_def_var_deflate(hs%ncid,
     &            hs%yvelva%nodal_data_id, 1, 1, 2)
            CALL check_err(iret)
         ENDIF
      ENDIF
#endif
!
!     Define harmonic analysis parameters
      call defineParameterWithText(hs%ncid, 'ntsteps', NF90_INT,
     &     "number of time steps since start of means and variance",
     &     "number of time steps since start of means and variance")
C
C
C     Leave define mode
      iret = nf90_enddef(hs%ncid)
      CALL check_err(iret)
C
C     now close the initialized netcdf file
      iret = nf90_close(hs%ncid)
      CALL check_err(iret)
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE initNetCDFHotstartHarmonicMeansVariances
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E
C        I N I T   N E T C D F   H O T S T A R T   3 D
C-----------------------------------------------------------------------
C     jgf49.49.02 Sets up netCDF variables for hotstarting 3D.
C-----------------------------------------------------------------------
      SUBROUTINE initNetCDFHotstart3D(lun, netcdf_format)
      USE SIZES, ONLY : MNPROC
      USE GLOBAL, ONLY : OutputDataDescript_t, IDEN, scratchMessage
      USE GLOBAL_3DVS, ONLY : NFEN
      IMPLICIT NONE
C
      INTEGER, intent(in) :: lun
      INTEGER, intent(in) :: netcdf_format ! whether netcdf3 or netcdf4
                                           ! format (netcdf4=hdf5)
                                           ! classic data model in any case
C
      LOGICAL err
C
      INTEGER iret ! success or failure of the netcdf call
      CHARACTER(1024) :: att_text
      INTEGER :: tempid
C
      call setMessageSource("initNetCDFHotstart3D")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      err = .false.
C
C     Point to the hotstart file we want to work on.
      IF (lun.eq.67) THEN
         hs => hs67
      ELSE
         hs => hs68
      ENDIF
C
C     nodes and elements
      hs%myMesh => adcircMesh
C     jgf50.06: If this is called during a parallel run, and all other
C     fulldomain 3D output is turned off, the number of vertical nodes will not
C     have been set.
      IF (MNPROC.gt.1) THEN
         hs%myMesh%num_v_nodes = NFEN
      ENDIF
C
C     Open existing NetCDF file
      iret = nf90_open(hs%myFile%FILENAME, NF_WRITE, hs%ncid)
      CALL check_err(iret)
C
C     Enter "redefine" mode
      iret = nf90_redef(hs%ncid)
      CALL check_err(iret)
!
!     Inquire variables
      CALL getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)
C
C     DUU
      hs%duu%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%duu%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'duu', NF90_DOUBLE,
     &    hs%duu%nodal_data_dims, hs%duu%nodal_data_id)
      CALL check_err(iret)
      att_text = "velocity dispersion term"
      iret = nf90_put_att(hs%ncid, hs%duu%nodal_data_id,
     &   'long_name', trim(att_text))
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%duu%nodal_data_id,
     &      '_FillValue', doubleval)
      CALL check_err(iret)
C
C     DUV
      hs%duv%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%duv%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'duv', NF90_DOUBLE,
     &    hs%duv%nodal_data_dims, hs%duv%nodal_data_id)
      CALL check_err(iret)
      att_text = "velocity dispersion term"
      iret = nf90_put_att(hs%ncid, hs%duv%nodal_data_id,
     &   'long_name', trim(att_text))
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%duv%nodal_data_id,
     &      '_FillValue', doubleval)
      CALL check_err(iret)
C
C     DVV
      hs%dvv%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%dvv%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'dvv', NF90_DOUBLE,
     &    hs%dvv%nodal_data_dims, hs%dvv%nodal_data_id)
      CALL check_err(iret)
      att_text = "velocity dispersion term"
      iret = nf90_put_att(hs%ncid, hs%dvv%nodal_data_id,
     &   'long_name', trim(att_text))
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%dvv%nodal_data_id,
     &      '_FillValue',  doubleval)
      CALL check_err(iret)
C
C     UU
      hs%uu%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%uu%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'uu', NF90_DOUBLE,
     &    hs%uu%nodal_data_dims, hs%uu%nodal_data_id)
      CALL check_err(iret)
      att_text = "vertically averaged velocity in east direction"
      iret = nf90_put_att(hs%ncid, hs%uu%nodal_data_id,
     &   'long_name',trim(att_text))
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%uu%nodal_data_id,
     &      '_FillValue', doubleval)
      CALL check_err(iret)
C
C     VV
      hs%vv%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%vv%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'vv', NF90_DOUBLE,
     &    hs%vv%nodal_data_dims, hs%vv%nodal_data_id)
      CALL check_err(iret)
      att_text = "vertically averaged velocity in north direction"
      iret = nf90_put_att(hs%ncid, hs%vv%nodal_data_id,
     &   'long_name', trim(att_text))
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%vv%nodal_data_id,
     &      '_FillValue', doubleval)
      CALL check_err(iret)
C
C     BSX
      hs%bsx%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%bsx%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'bsx', NF90_DOUBLE,
     &    hs%bsx%nodal_data_dims, hs%bsx%nodal_data_id)
      CALL check_err(iret)
      att_text = "bottom stress in east direction"
      iret = nf90_put_att(hs%ncid, hs%bsx%nodal_data_id,
     &   'long_name', trim(att_text))
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%bsx%nodal_data_id,
     &      '_FillValue',doubleval)
      CALL check_err(iret)
C
C     BSY
      hs%bsy%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%bsy%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'bsy', NF90_DOUBLE,
     &    hs%bsy%nodal_data_dims, hs%bsy%nodal_data_id)
      CALL check_err(iret)
      att_text = "bottom stress in north direction"
      iret = nf90_put_att(hs%ncid, hs%bsy%nodal_data_id,
     &   'long_name',trim(att_text))
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%bsy%nodal_data_id,
     &      '_FillValue', doubleval)
      CALL check_err(iret)
C
C     3D DENSITY
      IF (IDEN.ne.0) THEN
         hs%density3D%nodal_data_dims_3D(1) = hs%myMesh%num_nodes_dim_id
         hs%density3D%nodal_data_dims_3D(2)=hs%myMesh%num_v_nodes_dim_id
         hs%density3D%nodal_data_dims_3D(3) = hs%myTime%timenc_dim_id
      ENDIF
      IF (IDEN.eq.1) THEN
         iret = nf90_def_var(hs%ncid, 'sigt', NF90_DOUBLE,
     &     hs%density3D%nodal_data_dims_3D,hs%density3D%u_nodal_data_id)
         CALL check_err(iret)
         att_text = "sigma t density"
         iret = nf90_put_att(hs%ncid, hs%density3D%u_nodal_data_id,
     &      'long_name',trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%density3D%u_nodal_data_id,
     &         '_FillValue', doubleval)
         CALL check_err(iret)
      ENDIF
      IF ((IDEN.eq.2).or.(IDEN.eq.4)) THEN
         iret = nf90_def_var(hs%ncid, 'salinity', NF90_DOUBLE,
     &       hs%density3D%nodal_data_dims_3D,
     &       hs%density3D%v_nodal_data_id)
         CALL check_err(iret)
         att_text = "salinity"
         iret = nf90_put_att(hs%ncid, hs%density3D%v_nodal_data_id,
     &      'long_name',trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%density3D%v_nodal_data_id,
     &         '_FillValue', doubleval)
         CALL check_err(iret)
      ENDIF
      IF ((IDEN.eq.3).or.(IDEN.eq.4)) THEN
         iret = nf90_def_var(hs%ncid, 'temperature', NF90_DOUBLE,
     &       hs%density3D%nodal_data_dims_3D,
     &       hs%density3D%w_nodal_data_id)
         CALL check_err(iret)
         att_text = "salinity"
         iret = nf90_put_att(hs%ncid, hs%density3D%w_nodal_data_id,
     &      'long_name', trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%density3D%w_nodal_data_id,
     &         '_FillValue', doubleval)
         CALL check_err(iret)
      ENDIF
C
C     3D VELOCITY
      hs%velocity3D%nodal_data_dims_3D(1) = hs%myMesh%num_nodes_dim_id
      hs%velocity3D%nodal_data_dims_3D(2) = hs%myMesh%num_v_nodes_dim_id
      hs%velocity3D%nodal_data_dims_3D(3) = hs%myTime%timenc_dim_id
C     u-vel3D
      iret = nf90_def_var(hs%ncid, 'u-vel3D', NF90_DOUBLE,
     &    hs%velocity3D%nodal_data_dims_3D,
     &    hs%velocity3D%u_nodal_data_id)
      CALL check_err(iret)
      att_text = "3D fulldomain velocity in east direction"
      iret = nf90_put_att(hs%ncid, hs%velocity3D%u_nodal_data_id,
     &   'long_name', trim(att_text))
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%velocity3D%u_nodal_data_id,
     &      '_FillValue', doubleval)
      CALL check_err(iret)
C     v-vel3D
      iret = nf90_def_var(hs%ncid, 'v-vel3D', NF90_DOUBLE,
     &    hs%velocity3D%nodal_data_dims_3D,
     &    hs%velocity3D%v_nodal_data_id)
      CALL check_err(iret)
      att_text = "3D fulldomain velocity in north direction"
      iret = nf90_put_att(hs%ncid, hs%velocity3D%v_nodal_data_id,
     &   'long_name',trim(att_text))
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%velocity3D%v_nodal_data_id,
     &      '_FillValue', doubleval)
      CALL check_err(iret)
C     w-vel3D
      iret = nf90_def_var(hs%ncid, 'w-vel3D', NF90_DOUBLE,
     &    hs%velocity3D%nodal_data_dims_3D,
     &    hs%velocity3D%w_nodal_data_id)
      CALL check_err(iret)
      att_text = "3D full domain velocity in the vertical direction"
      iret = nf90_put_att(hs%ncid, hs%velocity3D%w_nodal_data_id,
     &   'long_name', trim(att_text))
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%velocity3D%w_nodal_data_id,
     &      '_FillValue', doubleval)
      CALL check_err(iret)
C
C     3D TURBULENCE
      hs%turbulence3D%nodal_data_dims_3D(1) = hs%myMesh%num_nodes_dim_id
      hs%turbulence3D%nodal_data_dims_3D(2)=hs%myMesh%num_v_nodes_dim_id
      hs%turbulence3D%nodal_data_dims_3D(3) = hs%myTime%timenc_dim_id
C     Q20
      iret = nf90_def_var(hs%ncid, 'q20', NF90_DOUBLE,
     &    hs%turbulence3D%nodal_data_dims_3D,
     &    hs%turbulence3D%u_nodal_data_id)
      CALL check_err(iret)
      att_text = "3D fulldomain turbulence kinetic energy"
      iret = nf90_put_att(hs%ncid, hs%turbulence3D%u_nodal_data_id,
     &   'long_name',trim(att_text))
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%turbulence3D%u_nodal_data_id,
     &      '_FillValue', doubleval)
      CALL check_err(iret)
C     L
      iret = nf90_def_var(hs%ncid, 'l', NF90_DOUBLE,
     &    hs%turbulence3D%nodal_data_dims_3D,
     &    hs%turbulence3D%v_nodal_data_id)
      CALL check_err(iret)
      att_text = "3D fulldomain turbulence length scale"
      iret = nf90_put_att(hs%ncid, hs%turbulence3D%v_nodal_data_id,
     &   'long_name', trim(att_text))
      CALL check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%turbulence3D%v_nodal_data_id,
     &      '_FillValue', doubleval)
      CALL check_err(iret)
C
C     jgf50.44: Automatically turn on compression if we are using the
C     netcdf4 file format.
#ifdef NETCDF_CAN_DEFLATE
      IF ( (netcdf_format.eq.5).or.
     &     (netcdf_format.eq.567).or.
     &     (netcdf_format.eq.568) ) THEN
         iret = nf_def_var_deflate(hs%ncid, hs%duu%nodal_data_id,
     &             1, 1, 2)
         CALL check_err(iret)
         iret = nf_def_var_deflate(hs%ncid, hs%dvv%nodal_data_id,
     &             1, 1, 2)
         CALL check_err(iret)
         iret = nf_def_var_deflate(hs%ncid, hs%uu%nodal_data_id,
     &             1, 1, 2)
         CALL check_err(iret)
         iret = nf_def_var_deflate(hs%ncid, hs%vv%nodal_data_id,
     &             1, 1, 2)
         CALL check_err(iret)
         iret = nf_def_var_deflate(hs%ncid, hs%bsx%nodal_data_id,
     &             1, 1, 2)
         CALL check_err(iret)
         iret = nf_def_var_deflate(hs%ncid, hs%bsy%nodal_data_id,
     &             1, 1, 2)
         CALL check_err(iret)
         IF (IDEN.eq.1) THEN
            iret = nf_def_var_deflate(hs%ncid,
     &               hs%density3D%u_nodal_data_id, 1, 1, 2)
            CALL check_err(iret)
         ENDIF
         IF ((IDEN.eq.2).or.(IDEN.eq.4)) THEN
            iret = nf_def_var_deflate(hs%ncid,
     &           hs%density3D%v_nodal_data_id, 1, 1, 2)
            CALL check_err(iret)
         ENDIF
         IF ((IDEN.eq.3).or.(IDEN.eq.4)) THEN
            iret = nf_def_var_deflate(hs%ncid,
     &           hs%density3D%w_nodal_data_id, 1, 1, 2)
            CALL check_err(iret)
         ENDIF
         iret = nf_def_var_deflate(hs%ncid,
     &          hs%velocity3D%u_nodal_data_id, 1, 1, 2)
         CALL check_err(iret)
         iret = nf_def_var_deflate(hs%ncid,
     &          hs%velocity3D%v_nodal_data_id, 1, 1, 2)
         CALL check_err(iret)
         iret = nf_def_var_deflate(hs%ncid,
     &          hs%velocity3D%w_nodal_data_id, 1, 1, 2)
         CALL check_err(iret)
         iret = nf_def_var_deflate(hs%ncid,
     &          hs%turbulence3D%u_nodal_data_id, 1, 1, 2)
         CALL check_err(iret)
         iret = nf_def_var_deflate(hs%ncid,
     &          hs%turbulence3D%v_nodal_data_id, 1, 1, 2)
         CALL check_err(iret)

      ENDIF
#endif
C
      iret = nf90_def_var(hs%ncid, 'n3dsd', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_def_var(hs%ncid, 'i3dsdrec', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_def_var(hs%ncid, 'n3dsv', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_def_var(hs%ncid, 'i3dsvrec', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_def_var(hs%ncid, 'n3dst', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_def_var(hs%ncid, 'i3dstrec', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_def_var(hs%ncid, 'n3dgd', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_def_var(hs%ncid, 'i3dgdrec', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_def_var(hs%ncid, 'n3dgv', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_def_var(hs%ncid, 'i3dgvrec', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_def_var(hs%ncid, 'n3dgt', NF90_INT, varid=tempid)
      CALL check_err(iret)
      iret = nf90_def_var(hs%ncid, 'i3dgtrec', NF90_INT, varid=tempid)
      CALL check_err(iret)
C
C
C     Leave define mode
      iret = nf90_enddef(hs%ncid)
      CALL check_err(iret)
C
C     now close the initialized netcdf file
      iret = nf90_close(hs%ncid)
      CALL check_err(iret)
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE initNetCDFHotstart3D
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E
C        D E F I N E   P A R A M E T E R   W I T H   T E X T
C-----------------------------------------------------------------------
C     jgf49.44 Defines a variable in the netcdf file and associates
C     attribute text with it.
C-----------------------------------------------------------------------
      SUBROUTINE defineParameterWithText(ncid, param, varType,
     &                                  longName, standardName)
      IMPLICIT NONE
      INTEGER, intent(in) :: ncid
      CHARACTER(len=*), intent(in) :: param
      INTEGER, intent(in) :: varType      ! netcdf-defined data type
      CHARACTER(len=*), intent(in) :: longName
      CHARACTER(len=*), intent(in) :: standardName
C
      INTEGER tempid  ! variable id for attaching to text
      INTEGER iret   ! netcdf err indicator
C
      call setMessageSource("defineParameterWithText")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      iret = nf90_def_var(ncid, param, varType, varid=tempid)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, tempid, 'long_name', trim(longName))
      CALL check_err(iret)
      iret = nf90_put_att(ncid, tempid, 'standard_name',
     &       trim(standardName))
      CALL check_err(iret)
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE defineParameterWithText
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E    W R I T E   N E T C D F   H O T S T A R T
C-----------------------------------------------------------------------
C     jgf49.17.02 Writes data to hotstart file.
C-----------------------------------------------------------------------
      SUBROUTINE writeNetCDFHotstart(lun, Elev1Descript,
     &        Elev2Descript, VelDescript, CH1Descript, EtaDiscDescript,
     &        NodeCodeDescript, NOFFDescript, timesec, it)
      USE SIZES, ONLY : globaldir, mnproc
      USE GLOBAL, ONLY : OutputDataDescript_t,
     &                   im, iestp, nscoue, ivstp, imhs,
     &                   nscouv, icstp, nscouc, ipstp, iwstp, nscoum,
     &                   igep, nscouge, igvp, nscougv, igcp, nscougc,
     &                   igpp, igwp, nscougw
      IMPLICIT NONE
C
      INTEGER, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: Elev1Descript
      type(OutputDataDescript_t), intent(in) :: Elev2Descript
      type(OutputDataDescript_t), intent(in) :: VelDescript
      type(OutputDataDescript_t), intent(in) :: CH1Descript
      type(OutputDataDescript_t), intent(in) :: EtaDiscDescript
      type(OutputDataDescript_t), intent(in) :: NodeCodeDescript
      type(OutputDataDescript_t), intent(in) :: NOFFDescript
      REAL(8), intent(in) :: timesec
      INTEGER, intent(in) :: it   ! current ADCIRC time step
C
      INTEGER i,j
      INTEGER counti(1), starti(1), n
      INTEGER kount(2), start(2) ! for nodally based data
      INTEGER elekount(2) ! for elementally based data
      INTEGER iret ! success or failure of the netcdf call
      INTEGER tempid
      CHARACTER(len=10) :: fext
C
      call setMessageSource("writeNetCDFHotstart")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
C     Point to the hotstart file we want to work on
      IF (lun.eq.67) THEN
         hs => hs67
      ELSE
         hs => hs68
      ENDIF
      hs%myMesh => adcircMesh
C
      ! create file name
      write(fext,'(i0)') lun
      hs%myFile%filename = trim(globaldir)//'/fort.'//trim(fext)//'.nc'
      call logMessage(INFO,'Opening hotstart file "'
     &    //trim(hs%myFile%filename)//'" for writing.')
C
      iret = nf90_open(hs%myFile%FILENAME, NF_WRITE, hs%ncid)
      CALL check_err(iret)
!
!     Inquire variables
      CALL getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)
      IF (hs%myTime%initialized.eqv..false.) THEN
         ALLOCATE(hs%myTime%timenc(hs%myTime%timenc_len))
         hs%myTime%initialized = .true.
      ENDIF
C
C
C     Don't increment the record counter so that we can store data at the
C     next location in the netcdf file -- the hotstart file is only intended
C     to have a single snapshot of data in it.
C      hs%myFile%record_counter = hs%myFile%record_counter + 1
      hs%myFile%record_counter = 1
!
!     Store time
      iret=nf90_inq_varid(hs%ncid,"time",hs%myTime%timenc_id)
      starti(1)=hs%myFile%record_counter
      counti(1)=hs%myTime%timenc_len
      hs%myTime%timenc(hs%myTime%timenc_len)=timesec
      iret = nf90_put_var(hs%ncid, hs%myTime%timenc_id,hs%myTime%timenc,
     &       starti, counti)
      CALL check_err(iret)

C
      kount(1)=hs%myMesh%num_nodes
      kount(2)=hs%myTime%timenc_len
      elekount(1)=hs%myMesh%num_elems
      elekount(2)=kount(2)
      start(1)=1
      start(2)=hs%myFile%record_counter
C
C     Get the NetCDF IDs of the relevant variables from the file
      iret=nf90_inq_varid(hs%ncid,"zeta1",hs%zeta1%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"zeta2",hs%zeta2%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"zetad",hs%zetad%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"u-vel",hs%vel%u_nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"v-vel",hs%vel%v_nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"nodecode",hs%nodecodenc%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"noff",hs%noffnc%nodal_data_id)
      CALL check_err(iret)
C
C     Write the nodal data to the netcdf file
      IF (MNPROC.eq.1) THEN
         iret = nf90_put_var(hs%ncid, hs%zeta1%nodal_data_id,
     &       Elev1Descript%array, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%zeta2%nodal_data_id,
     &      Elev2Descript%array, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%zetad%nodal_data_id,
     &       EtaDiscDescript%array, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%vel%u_nodal_data_id,
     &       VelDescript%array, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%vel%v_nodal_data_id,
     &      VelDescript%array2, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%nodecodenc%nodal_data_id,
     &      NodeCodeDescript%iarray, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%noffnc%nodal_data_id,
     &       NOFFDescript%iarray, start, elekount)
         CALL check_err(iret)
      ELSE
         iret = nf90_put_var(hs%ncid, hs%zeta1%nodal_data_id,
     &       Elev1Descript%array_g, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%zeta2%nodal_data_id,
     &      Elev2Descript%array_g, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%zetad%nodal_data_id,
     &       EtaDiscDescript%array_g, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%vel%u_nodal_data_id,
     &       VelDescript%array_g, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%vel%v_nodal_data_id,
     &      VelDescript%array2_g, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%nodecodenc%nodal_data_id,
     &      NodeCodeDescript%iarray_g, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%noffnc%nodal_data_id,
     &       NOFFDescript%iarray_g, start, elekount)
         CALL check_err(iret)
      ENDIF

C
C     Get each variable ID for the model parameters in the netcdf file
C     and immediately write the parameter value to that variable ID before
C     going on to the next one.
      iret = nf90_inq_varid(hs%ncid,"imhs",tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, im)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"iths",tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, it)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"iestp",tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, iestp)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nscoue",tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nscoue)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"ivstp",tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, ivstp)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nscouv",tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nscouv)
      CALL check_err(iret)
      IF ((IM.EQ.10).OR.(IMHS.EQ.10)) THEN
         iret = nf90_inq_varid(hs%ncid,"icstp",tempid)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, tempid, icstp)
         CALL check_err(iret)
         iret = nf90_inq_varid(hs%ncid,"nscouc",tempid)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, tempid, nscouc)
         CALL check_err(iret)
      ENDIF
      iret = nf90_inq_varid(hs%ncid,"ipstp",tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, ipstp)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"iwstp",tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, iwstp)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nscoum",tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nscoum)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"igep",tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, igep)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nscouge",tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nscouge)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"igvp",tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, igvp)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nscougv",tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nscougv)
      CALL check_err(iret)
      IF ((IM.EQ.10).OR.(IMHS.EQ.10)) THEN
         iret = nf90_inq_varid(hs%ncid,"igcp",tempid)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, tempid, igcp)
         CALL check_err(iret)
         iret = nf90_inq_varid(hs%ncid,"nscougc",tempid)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, tempid, nscougc)
         CALL check_err(iret)
      ENDIF
      iret = nf90_inq_varid(hs%ncid,"igpp",tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, igpp)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"igwp",tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, igwp)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nscougw",tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nscougw)
      CALL check_err(iret)
C
C     now close the netcdf file
      iret = nf90_close(hs%ncid)
      CALL check_err(iret)
      !
      ! jgf52.17: If we started with this hotstart file, the
      ! metadata still reflects the fort.15 from the previous run,
      ! instead of this one, so we need to update it.
      CALL updateMetaData(hs%ncid,hs%myFile)

C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE writeNetCDFHotstart
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E
C        W R I T E   N E T C D F   H O T S T A R T   H A R M O N I C
C-----------------------------------------------------------------------
C     jgf49.44.03 Writes harmonic analysis data to hotstart file.
C-----------------------------------------------------------------------
      SUBROUTINE writeNetCDFHotstartHarmonic(lun,
     &            GLOELVDescript, STAELVDescript,
     &            GLOULVDescript, GLOVLVDescript,
     &            STAULVDescript, STAVLVDescript)
      USE SIZES, ONLY : MNHARF, MNPROC, globaldir
      USE GLOBAL, ONLY : OutputDataDescript_t, NSTAE_G, NSTAV_G
      USE HARM, ONLY : nz, nf, mm, nhase, nhasv, nhage, nhagv, icall,
     &                 nfreq, timeud, itud, namefr, hafreq, haff,
     &                 haface, ha, icha
      IMPLICIT NONE
C
      INTEGER, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: GLOELVDescript
      type(OutputDataDescript_t), intent(in) :: STAELVDescript
      type(OutputDataDescript_t), intent(in) :: GLOULVDescript
      type(OutputDataDescript_t), intent(in) :: GLOVLVDescript
      type(OutputDataDescript_t), intent(in) :: STAULVDescript
      type(OutputDataDescript_t), intent(in) :: STAVLVDescript
C
      INTEGER i,j
      INTEGER kount(3), start(3) ! for nodally based data
      INTEGER hakount(2), hastart(2) ! for lhs
      INTEGER iret ! success or failure of the netcdf call
      INTEGER tempid
C
      call setMessageSource("writeNetCDFHotstartHarmonic")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
C     Point to the hotstart file we want to work on
      IF (lun.eq.67) THEN
         hs => hs67
      ELSE
         hs => hs68
      ENDIF
C
      iret = nf90_open(hs%myFile%FILENAME, NF_WRITE, hs%ncid)
      CALL check_err(iret)
!
!     Inquire variables (time)
      CALL getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)

      IF (NHASE.ne.0) THEN
C        elevation station dimension
         iret=nf90_inq_dimid(hs%ncid,"elevstation",
     &        hs%staelv%num_sta_dim_id)
         CALL check_err(iret)
         iret=nf90_inquire_dimension(hs%ncid,hs%staelv%num_sta_dim_id,
     &        len=hs%staelv%num_stations)
         CALL check_err(iret)
      ENDIF
      IF (NHASV.ne.0) THEN
C        velocity station dimension
         iret=nf90_inq_dimid(hs%ncid,"velstation",
     &        hs%stavellv%num_sta_dim_id)
         CALL check_err(iret)
         iret=nf90_inquire_dimension(hs%ncid,hs%stavellv%num_sta_dim_id,
     &        len=hs%stavellv%num_stations)
         CALL check_err(iret)
      ENDIF
C
C
C     Don't increment the record counter so that we can store data at the
C     next location in the netcdf file -- the hotstart file is only intended
C     to have a single snapshot of data in it.
      hs%myFile%record_counter = 1
C
      kount(1)=MNHARF*2             ! for load vector data
      kount(2)=hs%myMesh%num_nodes  ! for nodal data
      kount(3)=hs%myTime%timenc_len
      start(1)=1
      start(2)=1
      start(3)=hs%myFile%record_counter
C
C     Get the NetCDF IDs of the relevant variables from the file
      IF (NHAGE.ne.0) THEN
         iret=nf90_inq_varid(hs%ncid,"gloelv",hs%gloelv%nodal_data_id)
         CALL check_err(iret)
      ENDIF
      IF (NHASE.ne.0) THEN
         iret=nf90_inq_varid(hs%ncid,"staelv",hs%staelv%station_data_id)
         CALL check_err(iret)
      ENDIF
      IF (NHAGV.ne.0) THEN
       iret=nf90_inq_varid(hs%ncid,"gloulv",hs%glovellv%u_nodal_data_id)
       CALL check_err(iret)
       iret=nf90_inq_varid(hs%ncid,"glovlv",hs%glovellv%v_nodal_data_id)
       CALL check_err(iret)
      ENDIF
      IF (NHASV.ne.0) THEN
         iret=nf90_inq_varid(hs%ncid,"staulv",
     &                     hs%stavellv%u_station_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(hs%ncid,"stavlv",
     &                     hs%stavellv%v_station_data_id)
         CALL check_err(iret)
      ENDIF
C
C     Write the nodal data to the netcdf file
      IF (MNPROC.eq.1) THEN
         IF (NHAGE.ne.0) THEN
            kount(2)=hs%myMesh%num_nodes  ! for nodal data
            iret = nf90_put_var(hs%ncid, hs%gloelv%nodal_data_id,
     &          GLOELVDescript%array2D, start, kount)
            CALL check_err(iret)
         ENDIF
         IF (NHASE.ne.0) THEN
            kount(2) = hs%staelv%num_stations ! for elevation stations
            iret = nf90_put_var(hs%ncid, hs%staelv%station_data_id,
     &         STAELVDescript%array2D, start, kount)
            CALL check_err(iret)
         ENDIF
         IF (NHAGV.ne.0) THEN
            kount(2)=hs%myMesh%num_nodes  ! for nodal data
            iret = nf90_put_var(hs%ncid, hs%glovellv%u_nodal_data_id,
     &          GLOULVDescript%array2D, start, kount)
            CALL check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%glovellv%v_nodal_data_id,
     &          GLOVLVDescript%array2D, start, kount)
            CALL check_err(iret)
         ENDIF
         IF (NHASV.ne.0) THEN
            kount(2) = hs%stavellv%num_stations ! for velocity stations
            iret = nf90_put_var(hs%ncid, hs%stavellv%u_station_data_id,
     &          STAULVDescript%array2D, start, kount)
            CALL check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%stavellv%v_station_data_id,
     &          STAVLVDescript%array2D, start, kount)
            CALL check_err(iret)
         ENDIF
      ELSE
         IF (NHAGE.ne.0) THEN
            kount(2)=hs%myMesh%num_nodes  ! for nodal data
            iret = nf90_put_var(hs%ncid, hs%gloelv%nodal_data_id,
     &          GLOELVDescript%array2D_g, start, kount)
            CALL check_err(iret)
         ENDIF
         IF (NHASE.ne.0) THEN
            kount(2) = hs%staelv%num_stations ! for elevation stations
            iret = nf90_put_var(hs%ncid, hs%staelv%station_data_id,
     &         STAELVDescript%array2D_g, start, kount)
            CALL check_err(iret)
         ENDIF
         IF (NHAGV.ne.0) THEN
            kount(2)=hs%myMesh%num_nodes  ! for nodal data
            iret = nf90_put_var(hs%ncid, hs%glovellv%u_nodal_data_id,
     &          GLOULVDescript%array2D_g, start, kount)
            CALL check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%glovellv%v_nodal_data_id,
     &          GLOVLVDescript%array2D_g, start, kount)
            CALL check_err(iret)
         ENDIF
         IF (NHASV.ne.0) THEN
            kount(2) = hs%stavellv%num_stations ! for velocity stations
            iret = nf90_put_var(hs%ncid, hs%stavellv%u_station_data_id,
     &          STAULVDescript%array2D_g, start, kount)
            CALL check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%stavellv%v_station_data_id,
     &          STAVLVDescript%array2D_g, start, kount)
            CALL check_err(iret)
         ENDIF
      ENDIF
C
C     Get each variable ID for the model parameters in the netcdf file
C     and immediately write the parameter value to that variable ID before
C     going on to the next one.
      iret = nf90_inq_varid(hs%ncid,"icha", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, icha)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nz", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nz)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nf", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nf)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"mm", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, mm)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nstae", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nstae_g)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nstav", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nstav_g)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nhase", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nhase)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nhasv", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nhasv)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nhage", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nhage)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nhagv", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nhagv)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"icall", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, icall)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nfreq", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nfreq)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"timeud", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, timeud)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"itud", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, itud)
      CALL check_err(iret)
C
C     left hand side
      hakount(1) = 2*MNHARF
      hakount(2) = 2*MNHARF
      hastart(1) = 1
      hastart(2) = 1
      iret = nf90_inq_varid(hs%ncid,"ha",hs%ha_id)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, hs%ha_id, ha, hastart, hakount)
      CALL check_err(iret)
C     frequency names
      iret = nf90_inq_varid(hs%ncid,"namefr",hs%namefr_id)
      CALL check_err(iret)
      do i=1,mnharf
         start(1)=1
         start(2)=i
         kount(1)= len(namefr(i))
         kount(2)=1
         iret = nf90_put_var(hs%ncid, hs%namefr_id, namefr(i),
     &       start, kount)
         CALL check_err(iret)
      end do
C     harmonic constituents
      start(1) = 1
      start(2) = 1
      kount(1) = MNHARF ! for constituents
      kount(2) = 1
      iret = nf90_inq_varid(hs%ncid,"hafreq",hs%hafreq_id)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, hs%hafreq_id, hafreq, start, kount)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"haff",hs%haff_id)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, hs%haff_id, haff, start, kount)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"haface",hs%haface_id)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, hs%haface_id, haface, start, kount)
      CALL check_err(iret)
C
C     now close the netcdf file
      iret = nf90_close(hs%ncid)
      CALL check_err(iret)
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE writeNetCDFHotstartHarmonic
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E
C        W R I T E   N E T C D F   H O T S T A R T   H A R M O N I C
C                                     M E A N S   V A R I A N C E S
C-----------------------------------------------------------------------
C     jgf49.43.14 Writes harmonic analysis means and variance data to
C     hotstart file.
C-----------------------------------------------------------------------
      SUBROUTINE writeNetCDFHotstartHarmonicMeansVariances(lun,
     &            ELAVDescript, ELVADescript,
     &            XVELAVDescript, YVELAVDescript,
     &            XVELVADescript, YVELVADescript)
      USE SIZES, ONLY : MNHARF, MNPROC, globaldir
      USE GLOBAL, ONLY : OutputDataDescript_t, NSTAE_G, NSTAV_G
      USE HARM, ONLY : nhase, nhasv, nhage, nhagv, ntsteps
      IMPLICIT NONE
C
      INTEGER, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: ELAVDescript
      type(OutputDataDescript_t), intent(in) :: ELVADescript
      type(OutputDataDescript_t), intent(in) :: XVELAVDescript
      type(OutputDataDescript_t), intent(in) :: YVELAVDescript
      type(OutputDataDescript_t), intent(in) :: XVELVADescript
      type(OutputDataDescript_t), intent(in) :: YVELVADescript
C
      INTEGER i,j
      INTEGER kount(2), start(2) ! for nodally based data
      INTEGER iret ! success or failure of the netcdf call
      INTEGER tempid
C
      call setMessageSource("writeNetCDFHotstartHarmonicMeansVariances")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
C     Point to the hotstart file we want to work on
      IF (lun.eq.67) THEN
         hs => hs67
      ELSE
         hs => hs68
      ENDIF
C
      iret = nf90_open(hs%myFile%FILENAME, NF_WRITE, hs%ncid)
      CALL check_err(iret)
!
!     Inquire variables
!     time dimension
      CALL getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)
C
C     Don't increment the record counter so that we can store data at the
C     next location in the netcdf file -- the hotstart file is only intended
C     to have a single snapshot of data in it.
      hs%myFile%record_counter = 1
C
      kount(1)=hs%myMesh%num_nodes  ! for nodal data
      kount(2)=hs%myTime%timenc_len
      start(1)=1
      start(2)=hs%myFile%record_counter
C
C     Get the NetCDF IDs of the relevant variables from the file
      IF (NHAGE.ne.0) THEN
         iret=nf90_inq_varid(hs%ncid,"elav",hs%elav%nodal_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(hs%ncid,"elva",hs%elva%nodal_data_id)
         CALL check_err(iret)
      ENDIF
      IF (NHAGV.ne.0) THEN
         iret=nf90_inq_varid(hs%ncid,"xvelav",hs%xvelav%nodal_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(hs%ncid,"yvelav",hs%yvelav%nodal_data_id)
         CALL check_err(iret)
C
         iret=nf90_inq_varid(hs%ncid,"xvelva",hs%xvelva%nodal_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(hs%ncid,"yvelva",hs%yvelva%nodal_data_id)
         CALL check_err(iret)
      ENDIF
C
C     Write the nodal data to the netcdf file
      IF (MNPROC.eq.1) THEN
         IF (NHAGE.ne.0) THEN
            iret = nf90_put_var(hs%ncid, hs%elav%nodal_data_id,
     &          ELAVDescript%array, start, kount)
            CALL check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%elva%nodal_data_id,
     &          ELVADescript%array, start, kount)
            CALL check_err(iret)
         ENDIF
         IF (NHAGV.ne.0) THEN
            iret = nf90_put_var(hs%ncid, hs%xvelav%nodal_data_id,
     &       XVELAVDescript%array, start, kount)
            CALL check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%yvelav%nodal_data_id,
     &       YVELAVDescript%array, start, kount)
            CALL check_err(iret)
C
            iret = nf90_put_var(hs%ncid, hs%xvelva%nodal_data_id,
     &        XVELVADescript%array, start, kount)
            CALL check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%yvelva%nodal_data_id,
     &        YVELVADescript%array, start, kount)
            CALL check_err(iret)
         ENDIF
      ELSE
         IF (NHAGE.ne.0) THEN
            iret = nf90_put_var(hs%ncid, hs%elav%nodal_data_id,
     &          ELAVDescript%array_g, start, kount)
            CALL check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%elva%nodal_data_id,
     &          ELVADescript%array_g, start, kount)
            CALL check_err(iret)
         ENDIF
         IF (NHAGV.ne.0) THEN
            iret = nf90_put_var(hs%ncid, hs%xvelav%nodal_data_id,
     &       XVELAVDescript%array_g, start, kount)
            CALL check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%yvelav%nodal_data_id,
     &       YVELAVDescript%array_g, start, kount)
            CALL check_err(iret)
C
            iret = nf90_put_var(hs%ncid, hs%xvelva%nodal_data_id,
     &        XVELVADescript%array_g, start, kount)
            CALL check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%yvelva%nodal_data_id,
     &        YVELVADescript%array_g, start, kount)
            CALL check_err(iret)
         ENDIF
      ENDIF

C
C     Get each variable ID for the model parameters in the netcdf file
C     and immediately write the parameter value to that variable ID before
C     going on to the next one.
      iret = nf90_inq_varid(hs%ncid,"ntsteps", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, ntsteps)
      CALL check_err(iret)
C
C     now close the netcdf file
      iret = nf90_close(hs%ncid)
      CALL check_err(iret)
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE writeNetCDFHotstartHarmonicMeansVariances
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E
C        W R I T E   N E T C D F   H O T S T A R T   3 D
C-----------------------------------------------------------------------
C     jgf49.49.02 Writes 3D data to hotstart file. Only includes arrays
C     of rank 1 (dimensioned by number of nodes, NP). The 3D variables
C     (i.e., of rank 2) will be written in a subsequent subroutine call.
C-----------------------------------------------------------------------
      SUBROUTINE writeNetCDFHotstart3D(lun,DUUDescript,
     &      DUVDescript, DVVDescript, UUDescript, VVDescript,
     &      BSXDescript, BSYDescript)
      USE SIZES, ONLY : MNPROC, globaldir
      USE GLOBAL, ONLY : OutputDataDescript_t
      USE GLOBAL_3DVS, ONLY : n3dsd, i3dsdrec, n3dsv, i3dsvrec,
     &   n3dst, i3dstrec, n3dgd, i3dgdrec, n3dgv, i3dgvrec,
     &   n3dgt, i3dgtrec
      IMPLICIT NONE
C
      INTEGER, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: DUUDescript
      type(OutputDataDescript_t), intent(in) :: DUVDescript
      type(OutputDataDescript_t), intent(in) :: DVVDescript
      type(OutputDataDescript_t), intent(in) :: UUDescript
      type(OutputDataDescript_t), intent(in) :: VVDescript
      type(OutputDataDescript_t), intent(in) :: BSXDescript
      type(OutputDataDescript_t), intent(in) :: BSYDescript
C
      INTEGER i,j
      INTEGER kount(2), start(2) ! for nodally based data
      INTEGER iret ! success or failure of the netcdf call
      INTEGER tempid
C
      call setMessageSource("writeNetCDFHotstart3D")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
C     Point to the hotstart file we want to work on
      IF (lun.eq.67) THEN
         hs => hs67
      ELSE
         hs => hs68
      ENDIF
C
      iret = nf90_open(hs%myFile%FILENAME, NF_WRITE, hs%ncid)
      CALL check_err(iret)
!
!     Inquire variables
!     time dimension
      CALL getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)
C
C
C     Don't increment the record counter so that we can store data at the
C     next location in the netcdf file -- the hotstart file is only intended
C     to have a single snapshot of data in it.
      hs%myFile%record_counter = 1
C
      kount(1)=hs%myMesh%num_nodes  ! for nodal data
      kount(2)=hs%myTime%timenc_len
      start(1)=1
      start(2)=hs%myFile%record_counter
C
C     Get the NetCDF IDs of the relevant variables from the file
      iret=nf90_inq_varid(hs%ncid,"duu",hs%duu%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"duv",hs%duv%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"dvv",hs%dvv%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"uu",hs%uu%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"vv",hs%vv%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"bsx",hs%bsx%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"bsy",hs%bsy%nodal_data_id)
      CALL check_err(iret)
C
C     Write the nodal data to the netcdf file
      IF (MNPROC.eq.1) THEN
         iret = nf90_put_var(hs%ncid, hs%duu%nodal_data_id,
     &       DUUDescript%array, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%duv%nodal_data_id,
     &       DUVDescript%array, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%dvv%nodal_data_id,
     &       DVVDescript%array, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%uu%nodal_data_id,
     &       UUDescript%array, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%vv%nodal_data_id,
     &       VVDescript%array, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%bsx%nodal_data_id,
     &       BSXDescript%array, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%bsy%nodal_data_id,
     &       BSYDescript%array, start, kount)
         CALL check_err(iret)
      ELSE
         iret = nf90_put_var(hs%ncid, hs%duu%nodal_data_id,
     &       DUUDescript%array_g, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%duv%nodal_data_id,
     &       DUVDescript%array_g, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%dvv%nodal_data_id,
     &       DVVDescript%array_g, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%uu%nodal_data_id,
     &       UUDescript%array_g, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%vv%nodal_data_id,
     &       VVDescript%array_g, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%bsx%nodal_data_id,
     &       BSXDescript%array_g, start, kount)
         CALL check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%bsy%nodal_data_id,
     &       BSYDescript%array_g, start, kount)
         CALL check_err(iret)
      ENDIF
C
C     Get each variable ID for the model parameters in the netcdf file
C     and immediately write the parameter value to that variable ID before
C     going on to the next one.
      iret = nf90_inq_varid(hs%ncid,"n3dsd", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, n3dsd)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"i3dsdrec", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, i3dsdrec)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"n3dsv", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, n3dsv)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"i3dsvrec", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, i3dsvrec)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"n3dst", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, n3dst)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"i3dstrec", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, i3dstrec)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"n3dgd", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, n3dgd)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"i3dgdrec", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, i3dgdrec)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"n3dgv", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, n3dgv)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"i3dgvrec", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, i3dgvrec)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"n3dgt", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, n3dgt)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"i3dgtrec", tempid)
      CALL check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, i3dgtrec)
      CALL check_err(iret)
C
C     now close the netcdf file
      iret = nf90_close(hs%ncid)
      CALL check_err(iret)
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE writeNetCDFHotstart3D
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E
C        W R I T E   N E T C D F   H O T S T A R T   3 D   V A R
C-----------------------------------------------------------------------
C     jgf49.49.02 Writes 3D data to hotstart file for a single variable
C     of rank 2.
C-----------------------------------------------------------------------
      SUBROUTINE writeNetCDFHotstart3DVar(lun, descript)
      USE SIZES, ONLY : MNPROC, globaldir
      USE GLOBAL, ONLY : OutputDataDescript_t
      IMPLICIT NONE
C
      INTEGER, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: descript
C
      INTEGER i,j
      INTEGER kount(3), start(3) ! for nodally based data
      INTEGER iret ! success or failure of the netcdf call
      INTEGER tempid
C
      call setMessageSource("writeNetCDFHotstart3DVar")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
C     Point to the hotstart file we want to work on
      IF (lun.eq.67) THEN
         hs => hs67
      ELSE
         hs => hs68
      ENDIF
C
      iret = nf90_open(hs%myFile%FILENAME, NF_WRITE, hs%ncid)
      CALL check_err(iret)
!
!     Inquire variables
!     time dimension
      CALL getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)
!     vertical dimension
      iret=nf90_inq_dimid(hs%ncid,"num_v_nodes",
     &        hs%myMesh%num_v_nodes_dim_id)
      CALL check_err(iret)
      iret=nf90_inquire_dimension(hs%ncid,hs%myMesh%num_v_nodes_dim_id,
     &     len=hs%myMesh%num_v_nodes)
      CALL check_err(iret)
C
C
C     Don't increment the record counter so that we can store data at the
C     next location in the netcdf file -- the hotstart file is only intended
C     to have a single snapshot of data in it.
      hs%myFile%record_counter = 1
C
      kount(1)=hs%myMesh%num_nodes
      kount(2)=hs%myMesh%num_v_nodes
      kount(3)=hs%myTime%timenc_len
      start(1)=1
      start(2)=1
      start(3)=hs%myFile%record_counter
C
C     Get the NetCDF ID of the relevant variable from the file
      SELECT CASE(trim(descript%field_name))
      CASE("SigmaT")
         iret=nf90_inq_varid(hs%ncid,"sigt",hs%density3D%u_nodal_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(hs%ncid, hs%density3D%u_nodal_data_id,
     &          descript%array2D, start, kount)
         ELSE
            iret = nf90_put_var(hs%ncid, hs%density3D%u_nodal_data_id,
     &          descript%array2D_g, start, kount)
         ENDIF
         CALL check_err(iret)
      CASE("Salinity")
         iret=nf90_inq_varid(hs%ncid,"salinity",
     &      hs%density3D%v_nodal_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(hs%ncid, hs%density3D%v_nodal_data_id,
     &         descript%array2D, start, kount)
          ELSE
            iret = nf90_put_var(hs%ncid, hs%density3D%v_nodal_data_id,
     &         descript%array2D_g, start, kount)
         ENDIF
         CALL check_err(iret)
      CASE("Temperature")
         iret = nf90_inq_varid(hs%ncid,"temperature",
     &      hs%density3D%nodal_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(hs%ncid,
     &          hs%density3D%w_nodal_data_id,
     &          descript%array2D, start, kount)
         ELSE
            iret = nf90_put_var(hs%ncid,
     &          hs%density3D%w_nodal_data_id,
     &          descript%array2D_g, start, kount)
         ENDIF
         CALL check_err(iret)
      CASE("u-vel3D")
         iret=nf90_inq_varid(hs%ncid,"u-vel3D",
     &      hs%velocity3D%u_nodal_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(hs%ncid, hs%velocity3D%u_nodal_data_id,
     &          descript%array2D, start, kount)
         ELSE
            iret = nf90_put_var(hs%ncid, hs%velocity3D%u_nodal_data_id,
     &          descript%array2D_g, start, kount)
         ENDIF
         CALL check_err(iret)
      CASE("v-vel3D")
         iret=nf90_inq_varid(hs%ncid,"v-vel3D",
     &       hs%velocity3D%v_nodal_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(hs%ncid,
     &          hs%velocity3D%v_nodal_data_id, descript%array2D,
     &          start, kount)
         ELSE
            iret = nf90_put_var(hs%ncid,
     &          hs%velocity3D%v_nodal_data_id, descript%array2D_g,
     &          start, kount)
         ENDIF
         CALL check_err(iret)
      CASE("w-vel3D")
         iret = nf90_inq_varid(hs%ncid, "w-vel3D",
     &       hs%velocity3D%w_nodal_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(hs%ncid, hs%velocity3D%w_nodal_data_id,
     &         descript%array2D, start, kount)
          ELSE
            iret = nf90_put_var(hs%ncid, hs%velocity3D%w_nodal_data_id,
     &         descript%array2D_g, start, kount)
         ENDIF
         CALL check_err(iret)
      CASE("q20")
         iret = nf90_inq_varid(hs%ncid, "q20",
     &       hs%turbulence3D%u_nodal_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(hs%ncid,
     &         hs%turbulence3D%u_nodal_data_id,
     &         descript%array2D, start, kount)
         ELSE
            iret = nf90_put_var(hs%ncid,
     &         hs%turbulence3D%u_nodal_data_id,
     &         descript%array2D_g, start, kount)
         ENDIF
         CALL check_err(iret)
      CASE("l")
         iret=nf90_inq_varid(hs%ncid,"l",hs%turbulence3D%v_nodal_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_put_var(hs%ncid,
     &          hs%turbulence3D%v_nodal_data_id,
     &          descript%array2D, start, kount)
         ELSE
            iret = nf90_put_var(hs%ncid,
     &          hs%turbulence3D%v_nodal_data_id,
     &          descript%array2D_g, start, kount)
         ENDIF
         CALL check_err(iret)
      CASE DEFAULT

      END SELECT
C
C     now close the netcdf file
      iret = nf90_close(hs%ncid)
      CALL check_err(iret)
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE writeNetCDFHotstart3DVar
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E    R E A D   N E T C D F   H O T S T A R T
C-----------------------------------------------------------------------
C     jgf49.17.02 Reads data from the hotstart file.
C-----------------------------------------------------------------------
      SUBROUTINE readNetCDFHotstart(lun, timeLoc)
      USE VERSION, ONLY : FileFmtRev, FileFmtMinor, FileFmtMajor
      USE SIZES, ONLY : globaldir, mnproc, myproc
      USE GLOBAL, ONLY : OutputDataDescript_t, imhs, iths, iestp, nscoue, ivstp,
     &                   nscouv, icstp, nscouc, ipstp, iwstp, nscoum,
     &                   igep, nscouge, igvp, nscougv, igcp, nscougc,
     &                   igpp, igwp, nscougw, ETA1, ETA2, EtaDisc,
     &                   UU2, VV2, NNODECODE, NOFF,
     &                   IM, NP_G, NE_G, NODES_LG, IMAP_EL_LG
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      IMPLICIT NONE
C
      INTEGER, intent(in) :: lun
      REAL(8), intent(out) :: timeLoc
C

      INTEGER i,j
      INTEGER counti(1), starti(1), n
      INTEGER kount(2), start(2)
      INTEGER iret ! success or failure of the netcdf call
      INTEGER :: tempid

      INTEGER :: sd_node_number    ! subdomain index to map
      INTEGER :: sd_element_number ! subdomain index to map
      CHARACTER(len=10) :: fext
C
      call setMessageSource("readNetCDFHotstart")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
C     Point to the hotstart file we want to work on
      IF (lun.eq.67) THEN
         hs => hs67
      ELSE
         hs => hs68
      ENDIF
C
      hs%myMesh => adcircMesh
      hs%myFile%fileFound = .false.
      ! create file name
      write(fext,'(i0)') lun
      hs%myFile%filename = trim(globaldir)//'/fort.'//trim(fext)//'.nc'
      call logMessage(INFO,'Opening hot start file "'
     &               //trim(hs%myFile%filename)//'" for reading.')
C
C     Open fulldomain file
      INQUIRE(FILE=hs%myFile%FILENAME,EXIST=hs%myFile%fileFound)
      IF (hs%myFile%fileFound.eqv..false.) THEN
         write(scratchMessage,'(a,a,a)') 'The file ',
     &       trim(adjustl(hs%myFile%FILENAME)),
     &       ' was not found; ADCIRC terminating.'
         call allMessage(ERROR,scratchMessage)
#ifdef CMPI
            CALL MSG_FINI()
#endif
         CALL EXIT(1)
      ELSE
         iret = nf90_open(hs%myFile%FILENAME, NF90_NOWRITE, hs%ncid)
         CALL check_err(iret)
      ENDIF
!
!     Inquire variables
      call getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)
C
C     nodes and elements
      hs%myMesh => adcircMesh
      IF (hs%myMesh%initialized.eqv..false.) THEN
         hs%myMesh%num_nodes = NP_G
         hs%myMesh%num_elems = NE_G
         hs%myMesh%nface_len = 3
      ENDIF
      IF (hs%myTime%initialized.eqv..false.) THEN
         ALLOCATE(hs%myTime%timenc(hs%myTime%timenc_len))
         hs%myTime%initialized = .true.
      ENDIF
C
      hs%myFile%record_counter = 1
      kount(1)=hs%myMesh%num_nodes
      kount(2)=hs%myTime%timenc_len
      start(1)=1
      start(2)=hs%myFile%record_counter
!
!     Get time
      iret=nf90_inq_varid(hs%ncid, "time", hs%myTime%timenc_id)
      CALL check_err(iret)
      starti(1)=hs%myFile%record_counter
      counti(1)=hs%myTime%timenc_len
      iret = nf90_get_var(hs%ncid, hs%myTime%timenc_id,
     &    hs%myTime%timenc, starti, counti)
      CALL check_err(iret)
      ! set timeLoc in hstart.F to the current time in the hotstart file
      timeLoc = hs%myTime%timenc(hs%myFile%record_counter)
      !
      ! get array variable ids
      iret=nf90_inq_varid(hs%ncid, "zeta1", hs%zeta1%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid, "zeta2", hs%zeta2%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid, "zetad", hs%zetad%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid, "u-vel", hs%vel%u_nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid, "v-vel", hs%vel%v_nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"nodecode",hs%nodecodenc%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid, "noff", hs%noffnc%nodal_data_id)
      CALL check_err(iret)
      ! serial
      IF (MNPROC.eq.1) THEN
C        Elev1
         iret = nf90_get_var(hs%ncid, hs%zeta1%nodal_data_id,
     &     eta1, start, kount)
         CALL check_err(iret)
C        Elev2
         iret = nf90_get_var(hs%ncid, hs%zeta2%nodal_data_id,
     &       eta2, start, kount)
         CALL check_err(iret)
C        EtaDisc
         iret = nf90_get_var(hs%ncid, hs%zetad%nodal_data_id,
     &      EtaDisc, start, kount)
         CALL check_err(iret)
C        Vel
         iret = nf90_get_var(hs%ncid, hs%vel%u_nodal_data_id,
     &      uu2, start, kount)
         CALL check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%vel%v_nodal_data_id,
     &      vv2, start, kount)
         CALL check_err(iret)
C        NodeCode
         iret = nf90_get_var(hs%ncid, hs%nodecodenc%nodal_data_id,
     &       nnodecode, start, kount)
         CALL check_err(iret)
C        NOFF
         start(1)=1
         kount(1)=hs%myMesh%num_elems
         iret = nf90_get_var(hs%ncid, hs%noffnc%nodal_data_id,
     &       noff, start, kount)
         CALL check_err(iret)

      ELSE       ! parallel
C
         ! Create the fulldomain node and element index lists for this subdomain
         call createFullDomainIndexLists()
         !
         ! get fulldomain data and map the data to this subdomain
         fullDomainIndexList => fullDomainNodeList
         call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes,
     &       hs%zeta1%nodal_data_id, subdomain_reals=eta1)
         call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes,
     &       hs%zeta2%nodal_data_id, subdomain_reals=eta2)
         call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes,
     &       hs%zetad%nodal_data_id, subdomain_reals=EtaDisc)
         call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes,
     &       hs%vel%u_nodal_data_id, subdomain_reals=uu2)
         call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes,
     &       hs%vel%v_nodal_data_id, subdomain_reals=vv2)
         call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes,
     &       hs%nodecodenc%nodal_data_id, subdomain_ints=nnodecode)
         fullDomainIndexList => fullDomainElementList
         call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_elems,
     &       hs%noffnc%nodal_data_id, subdomain_ints=noff)
      ENDIF
C
C     Read in model parameters to ADCIRC variables
      iret = nf90_inq_varid(hs%ncid,"imhs", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, imhs)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"iths", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, iths)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"iestp", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, iestp)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nscoue", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, nscoue)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"ivstp", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, ivstp)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nscouv", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, nscouv)
      CALL check_err(iret)
      IF ((IM.EQ.10).OR.(IMHS.EQ.10)) THEN
         iret = nf90_inq_varid(hs%ncid,"icstp", tempid)
         CALL check_err(iret)
         iret = nf90_get_var(hs%ncid, tempid, icstp)
         CALL check_err(iret)
         iret = nf90_inq_varid(hs%ncid,"nscouc", tempid)
         CALL check_err(iret)
         iret = nf90_get_var(hs%ncid, tempid, nscouc)
         CALL check_err(iret)
      ENDIF
      iret = nf90_inq_varid(hs%ncid,"ipstp", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, ipstp)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"iwstp", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, iwstp)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nscoum", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, nscoum)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"igep", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, igep)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nscouge", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, nscouge)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"igvp", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, igvp)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nscougv", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, nscougv)
      CALL check_err(iret)
      IF ((IM.EQ.10).OR.(IMHS.EQ.10)) THEN
         iret = nf90_inq_varid(hs%ncid,"igcp", tempid)
         CALL check_err(iret)
         iret = nf90_get_var(hs%ncid, tempid, igcp)
         CALL check_err(iret)
         iret = nf90_inq_varid(hs%ncid,"nscougc", tempid)
         CALL check_err(iret)
         iret = nf90_get_var(hs%ncid, tempid, nscougc)
         CALL check_err(iret)
      ENDIF
      iret = nf90_inq_varid(hs%ncid,"igpp", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, igpp)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"igwp", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, igwp)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nscougw", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, nscougw)
      CALL check_err(iret)
C
      DEALLOCATE(hs%myTime%timenc)
      hs%myTime%initialized = .false. !..zc50.93 - reinitialize if we write a hot start later
C
C     now close the netcdf file
      iret = nf90_close(hs%ncid)
      CALL check_err(iret)
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE readNetCDFHotstart
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C                    S U B R O U T I N E
C                 R E A D   A N D   M A P   T O
C         S U B D O M A I N   M I N   M A X   N E T C D F
C-----------------------------------------------------------------------
C     jgf52.08.03: Read in the min/max values from a previous run
C     at hotstart. In order to read min/max data from netCDF, we also
C     have to be writing that min/max data in netCDF. TODO: relax this
C     limitation.
C-----------------------------------------------------------------------
      subroutine readAndMapToSubdomainMaxMinNetCDF(descript, timeloc)
      use sizes, only : mnproc
      use global, only : OutputDataDescript_t, allMessage, logMessage, nfover
      implicit none
      type(OutputDataDescript_t), intent(inout) :: descript
      real(8), intent(in) :: timeloc ! adcirc time in seconds since cold start
C
      type(nodalData), pointer :: dat
      integer :: i,j
      integer :: counti(1), starti(1)
      integer :: kount(2), start(2)
      integer :: iret ! success or failure of the netcdf call
      logical :: timestampsFound   !.true. if file contains time of occurrence dataset
      character(len=1024) :: vn    ! min/max variable name in netcdf file
      character(len=1025) :: tvn   ! time of min/max occurrence variable name in netcdf

      call setMessageSource("readAndMapToSubdomainMaxMinNetCDF")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      select case(descript % lun)
      case(311) ! maxele
         dat => EtaMax
         vn = 'zeta_max'
         tvn = 'time_of_zeta_max'
      case(312) ! maxvel
         dat => UMax
         vn = 'vel_max'
         tvn = 'time_of_vel_max'
      case(313) ! minpr
         dat => PrMin
         vn = 'pressure_min'
         tvn = 'time_of_pressure_min'
      case(314) ! maxwvel
         dat => WVMax
         vn = 'wind_max'
         tvn = 'time_of_wind_max'
      case(315) ! maxrs
         dat => RSMax
         vn = 'radstress_max'
         tvn = 'time_of_radstress_max'
      case(400) ! inundationtime
         dat => inTime
         vn = 'inun_time'
         tvn = 'onset_inun_time'
      case(401) ! maxinundepth
         dat => maxInDep
         vn = 'inun_max'
         tvn = 'time_of_inun_max'
      case(404) ! everdried
         dat => evrDry
         vn = 'everdried'
         tvn = 'time_of_everdried'
      case default
         call allMessage(ERROR,'Files of the type "'
     &  //trim(descript % file_name)//'.nc" cannot be read by ADCIRC.')
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         RETURN ! RETURN RETURN
      end select

      dat%myFile%filename = trim(descript % file_name)//'.nc'
      call logMessage(INFO,'Opening file "'
     &               //trim(dat%myFile%filename)//'" for reading.')
C
C     Open fulldomain file
      inquire(file=dat%myfile%filename,exist=dat%myfile%filefound)
      if (dat%myFile%fileFound.eqv..false.) then
         call allMessage(INFO,'The file "'//trim(dat%myFile%FILENAME)//
     &       '" was not found.')
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         RETURN ! RETURN RETURN
      else
         ! the file was found, let's open it
         iret = nf90_open(dat%myFile%FILENAME, NF90_NOWRITE, dat%ncid)
         call check_err(iret)
      endif
C
C     nodes and elements
      dat%myMesh => adcircMesh
!
!     Inquire variables
      call getDimensions(dat%ncid, dat%myTime, dat%myMesh, dat%myFile)
      !
      ! jgf52.08.21: If the netcdf min/max file contains no time records,
      ! it was just created by adcprep, and there are no data to load.
      if ( dat%myFile%record_counter.eq.0 ) then
         call allMessage(INFO,'The file '//trim(dat%myFile%filename)//
     &      'contains no data, '//
     &      'so the min/max record will be started anew.')
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         !
         ! now close the netcdf file
         iret = nf90_close(dat%ncid)
         CALL check_err(iret)
         return ! EARLY RETURN
      else
         write(scratchMessage,3333) trim(descript % file_name)
3333    format('Values from ',(a),
     &       ' will be reflected from the solution prior to this hotstart.')
         call allMessage(INFO,scratchMessage)
      endif
      IF (dat%myTime%initialized.eqv..false.) THEN
         ALLOCATE(dat%myTime%timenc(dat%myTime%timenc_len))
         dat%myTime%initialized = .true.
      ENDIF
      !
      ! jgf52.08.20: Read time and reject min/max data that correspond
      ! to a future time. This can happen when an analyst hotstarts a run,
      ! allows it to complete, then makes some adjustments to the input and
      ! then tries to re-run it while forgetting to remove or replace the
      ! min/max file from the previous attempt.
      !
      ! It is hard to imagine a scenario where someone would want to keep
      ! the min/max record from the previous attempt under these circumstances.
      ! As a result, I am logging it as an error but allowing the run to
      ! continue if non-fatal override is enabled, since it would be
      ! annoying to have ADCIRC bomb out every time the analyst mistakenly
      ! leaves a min/max file in place.
      iret=nf90_inq_varid(dat%ncid, "time", dat%myTime%timenc_id)
      CALL check_err(iret)
      starti(1)=1  ! min/max files have one element in the time array
      counti(1)=1
      iret = nf90_get_var(dat%ncid, dat%myTime%timenc_id,
     &    dat%myTime%timenc, starti, counti)
      CALL check_err(iret)
      if (  dat%myTime%timenc(1) .gt. timeloc ) then
         call allMessage(ERROR,'Max/min file '//trim(dat%myFile%filename)//
     &     ' contains data written after the hotstart time. '//
     &   ' It may have been produced using different input parameters, '//
     &     'perhaps during a previously attempted hot start run. '//
     &     ' Its values will not be read.')
         if (nfover.eq.0) then
            call allMessage(ERROR,
     &       'Execution terminated due to invalid '//
     &        trim(dat%myFile%filename)//' file.')
             call netcdfTerminate()
         else
            call allMessage(INFO,'The record for '//
     &        trim(dat%myFile%filename)//
     &        ' will be started anew and exection will continue.')
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
             call allMessage(DEBUG,"Return.")
#endif
             call unsetMessageSource()
             !
             ! now close the netcdf file
             iret = nf90_close(dat%ncid)
             CALL check_err(iret)
             return ! EARLY RETURN
         endif
      endif

      dat%myFile%record_counter = 1
      kount(1)=dat%myMesh%num_nodes
      kount(2)=dat%myTime%timenc_len
      start(1)=1
      start(2)=dat%myFile%record_counter
      !
      ! get array variable ids
      !call logMessage(DEBUG,'netcdf var is '//trim(vn))
      iret=nf90_inq_varid(dat%ncid,trim(vn),dat%max_nodal_data_id)
      CALL check_err(iret)
      timestampsFound = .true.
      iret=nf90_inq_varid(dat%ncid,trim(tvn),dat%time_max_nodal_data_id)
      if ( iret.ne.NF90_NOERR ) then
         call logMessage(INFO,'The file "'//trim(dat%myFile%filename)//
     &  '" does not contain time of occurrence for the min/max data. '//
     &  'As a result, time of occurrence will be initialized to -99999.')
         timestampsFound = .false.
      endif
      !
      ! serial
      if (mnproc.eq.1) then

         if (descript % isInteger.eqv..true.) then
            ! integer
            iret = nf90_get_var(dat%ncid, dat%max_nodal_data_id,
     &             descript % iarray, start, kount)
            call check_err(iret)
         else
            ! float
            iret = nf90_get_var(dat%ncid, dat%max_nodal_data_id,
     &             descript % array, start, kount)
            call check_err(iret)
         endif
         ! time of occurrence
         if (timestampsFound.eqv..true.) then
            iret = nf90_get_var(dat%ncid, dat%time_max_nodal_data_id,
     &             descript % array2, start, kount)
            call check_err(iret)
         endif
C
      else       ! parallel
C
         ! Create the fulldomain node and element index lists for this subdomain
         ! if it has not already happened for the reading of the netcdf
         ! hotstart file.
         if (fullDomainIndexListsInitialized.eqv..false.) then
            call createFullDomainIndexLists()
         endif
         !
         ! get fulldomain data and map the data to this subdomain
         fullDomainIndexList => fullDomainNodeList
         if ( descript % isInteger.eqv..true.) then
            ! integer
            call mapFullDomainToSubdomain(dat%ncid, dat%myMesh%num_nodes,
     &          dat%max_nodal_data_id, subdomain_ints=descript%iarray)
         else
            ! float
            call mapFullDomainToSubdomain(dat%ncid, dat%myMesh%num_nodes,
     &          dat%max_nodal_data_id, subdomain_reals=descript%array)
         endif
         ! time of occurrence
         if (timestampsFound.eqv..true.) then
            call mapFullDomainToSubdomain(dat%ncid, dat%myMesh%num_nodes,
     &       dat%time_max_nodal_data_id, subdomain_reals=descript%array2)
         endif
      endif
C
C     now close the netcdf file
      iret = nf90_close(dat%ncid)
      CALL check_err(iret)

#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      end subroutine readAndMapToSubdomainMaxMinNetCDF
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E
C        C R E A T E   F U L L D O M A I N   I N D E X   L I S T S
C-----------------------------------------------------------------------
C     jgf52.08.03 Create fulldomain node and element index lists for
C     this subdomain. This subroutine uses code originally found in
C     readNetCDFHotstart; these lists need to be created for reading
C     min/max files as well as hotstart files. If a netCDF hotstart
C     file was already read, these index lists will have already been
C     initialized, so we wouldn't have to do it again.
C-----------------------------------------------------------------------
      subroutine createFullDomainIndexLists()
      use mesh, only : ne, np
      use global, only : nodes_lg, imap_el_lg
      implicit none
      integer :: sd_element_number
      integer :: sd_node_number
      call setMessageSource("createFullDomainIndexLists")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ! make a list of full domain nodes that correspond to this
      ! subdomain's nodes
      allocate(fullDomainNodeList(np))
      ! loop over subdomain indexes to form a list of corresponding
      ! fulldomain indexes
      forall (sd_node_number=1:np)
         ! get the corresponding fulldomain indexes
         fullDomainNodeList(sd_node_number)
     &         = ABS(nodes_lg(sd_node_number))
      end forall
      ! make a list of full domain elements that correspond to this
      ! subdomain's elements
      allocate(fullDomainElementList(ne))
      ! loop over subdomain indexes to form a list of corresponding
      ! fulldomain indexes
      forall (sd_element_number=1:ne)
         ! get the corresponding fulldomain indexes
         fullDomainElementList(sd_element_number)
     &         = abs(imap_el_lg(sd_element_number))
      end forall
      fullDomainIndexListsInitialized = .true.

#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      end subroutine createFullDomainIndexLists
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E
C        M A P   F U L L D O M A I N   T O   S U B D O M A I N
C-----------------------------------------------------------------------
C     jgf50.60.03 Maps full domain data to subdomain.
C-----------------------------------------------------------------------
      SUBROUTINE  mapFullDomainToSubdomain(ncid, fd_array_size,
     &       data_id, subdomain_reals, subdomain_ints)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: ncid ! file id to pull data from
      INTEGER, INTENT(in) :: fd_array_size ! highest index in fulldomain array
      INTEGER, INTENT(in) :: data_id ! netcdf variable id in file
      REAL(sz), OPTIONAL, INTENT(out) :: subdomain_reals(:) ! we need
      INTEGER, OPTIONAL, INTENT(out) :: subdomain_ints(:)  ! we need
C
      REAL(sz), ALLOCATABLE :: work_reals(:) ! holds fulldomain data
      INTEGER, ALLOCATABLE :: work_ints(:)   ! holds fulldomain data
      INTEGER iret ! success or failure of the netcdf call
C
      call setMessageSource("mapFullDomainToSubDomain")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      ! grab array of full domain data from netcdf file and pull out the
      ! data that is needed for this subdomain
      if (present(subdomain_reals)) then
         allocate(work_reals(fd_array_size))
         iret = nf90_get_var(ncid, data_id, work_reals)
         call check_err(iret)
         subdomain_reals(:) = work_reals(fullDomainIndexList)
         deallocate(work_reals)
      else
         allocate(work_ints(fd_array_size))
         iret = nf90_get_var(ncid, data_id, work_ints)
         call check_err(iret)
         subdomain_ints(:) = work_ints(fullDomainIndexList)
         deallocate(work_ints)
      endif
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE mapFulldomainToSubdomain
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E
C     M A P   F U L L D O M A I N   T O   S U B D O M A I N  M   B Y   NP
C
C-----------------------------------------------------------------------
C     jgf50.60.04 Maps mutlidimensional data from full domain data to
C     subdomain; specifically, this subroutine is useful for harmonic
C     analysis data.
C-----------------------------------------------------------------------
      SUBROUTINE  mapFullDomainToSubdomainMByNP(ncid, m, n,
     &       data_id, subdomain_reals)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: ncid ! file id to pull data from
      INTEGER, INTENT(in) :: m    ! non-nodal dimension
      INTEGER, INTENT(in) :: n    ! number of horizontal nodes
      INTEGER, INTENT(in) :: data_id ! netcdf variable id in file
      REAL(sz), OPTIONAL, INTENT(out) :: subdomain_reals(:,:) ! we need
C
      REAL(sz), ALLOCATABLE :: work_reals(:,:) ! holds fulldomain data
      INTEGER iret ! success or failure of the netcdf call
      INTEGER i    ! loop counter
C
      call setMessageSource("mapFullDomainToSubDomainMbyNP")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ! grab array of full domain data from netcdf file and pull out the
      ! data that is needed for this subdomain
      allocate(work_reals(m,n))
      iret = nf90_get_var(ncid, data_id, work_reals)
      call check_err(iret)
      do i=1,size(subdomain_reals(1,:))
         subdomain_reals(:,i) = work_reals(:,fullDomainIndexList(i))
      end do
      deallocate(work_reals)
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE mapFulldomainToSubdomainMByNP
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E
C     M A P   F U L L D O M A I N   T O   S U B D O M A I N  NP   B Y   M
C
C-----------------------------------------------------------------------
C     jgf50.60.04 Maps mutlidimensional data from full domain data to
C     subdomain; specifically, this subroutine is useful for 3D data.
C-----------------------------------------------------------------------
      SUBROUTINE  mapFullDomainToSubdomainNPByM(ncid, m, n,
     &       data_id, subdomain_reals)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: ncid ! file id to pull data from
      INTEGER, INTENT(in) :: m    ! non-nodal dimension
      INTEGER, INTENT(in) :: n    ! number of horizontal nodes
      INTEGER, INTENT(in) :: data_id ! netcdf variable id in file
      REAL(sz), OPTIONAL, INTENT(out) :: subdomain_reals(:,:) ! we need
C
      REAL(sz), ALLOCATABLE :: work_reals(:,:) ! holds fulldomain data
      INTEGER iret ! success or failure of the netcdf call
C
      call setMessageSource("mapFullDomainToSubDomainNPbyM")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ! grab array of full domain data from netcdf file and pull out the
      ! data that is needed for this subdomain
      allocate(work_reals(n,m))
      iret = nf90_get_var(ncid, data_id, work_reals)
      call check_err(iret)
      subdomain_reals(:,:) = work_reals(fullDomainIndexList,:)
      deallocate(work_reals)
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE mapFulldomainToSubdomainNPByM
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   G E T   D I M E N S I O N S
C-----------------------------------------------------------------------
C     jgf50.60.04 Pulls the number of nodes, number of elements, and
C     number of time records from the netcdf file, as well as the
C     dimension ids of those dimensions.
C-----------------------------------------------------------------------
      SUBROUTINE getDimensions(ncid, time_struct, mesh_struct,
     &    file_struct)
      IMPLICIT NONE
      integer :: ncid  ! file id
      type(timeData) :: time_struct
      type(fileData) :: file_struct
      type(meshStructure) :: mesh_struct
C
      integer :: iret ! success or failure of netcdf call
      integer :: ndim      ! number of dimensions in the netcdf file
      integer :: nvar      ! number of variables in the netcdf file
      integer :: natt      ! number of attributes in the netcdf file
C
      call setMessageSource("getDimensions")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
!
!     Inquire variables
      time_struct%timenc_dim_id = 0
      iret=nf90_inquire(ncid,ndim,nvar,natt,time_struct%timenc_dim_id)
      CALL check_err(iret)
       iret=nf90_inquire_dimension(ncid, time_struct%timenc_dim_id,
     &        len=file_struct%record_counter)
       CALL check_err(iret)
C
      iret=nf90_inq_dimid(ncid, "node", mesh_struct%num_nodes_dim_id)
      CALL check_err(iret)
      iret=nf90_inquire_dimension(ncid, mesh_struct%num_nodes_dim_id,
     &     len=mesh_struct%num_nodes)
      CALL check_err(iret)
C
      iret=nf90_inq_dimid(ncid, "nele", mesh_struct%num_elems_dim_id)
      CALL check_err(iret)
      iret=nf90_inquire_dimension(ncid, mesh_struct%num_elems_dim_id,
     &     len=mesh_struct%num_elems)
      CALL check_err(iret)
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE getDimensions
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E
C        R E A D   N E T C D F   H O T S T A R T   H A R M O N I C
C-----------------------------------------------------------------------
C     jgf49.44.11 Reads harmonic analysis data from the hotstart file.
C-----------------------------------------------------------------------
      SUBROUTINE readNetCDFHotstartHarmonic(lun)
      USE SIZES, ONLY : MNPROC, globaldir, myproc
      USE GLOBAL, ONLY : OutputDataDescript_t, NODES_LG,
     &                   IMAP_STAE_LG, NSTAE, IMAP_STAV_LG, NSTAV
      USE HARM, ONLY : GLOELV, STAELV, GLOULV, GLOVLV, STAULV, STAVLV,
     &                 NHASE, NHASV, NHAGE, NHAGV, MNHARF, ICHA, INZ,
     &                 INZ, INF, IMM, INSTAE, INSTAV, INHASE, INHASV,
     &                 INHASE, INHAGE, INHAGV, ICALL, INFREQ, TIMEUD,
     &                 ITUD, HA, INAMEFR, INP, IFF, IFACE, IFREQ
      USE MESH, ONLY : NP
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      IMPLICIT NONE
C
      INTEGER, intent(in) :: lun
C
      INTEGER i,j
      INTEGER counti(1), starti(1), n
      INTEGER kount(3), start(3)
      INTEGER hakount(2), hastart(2)
      INTEGER sd_station_number
      INTEGER iret ! success or failure of the netcdf call
      INTEGER tempid
      type(OutputDataDescript_t) :: descript
C
      call setMessageSource("readNetCDFHotstartHarmonic")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
C     Point to the hotstart file we want to work on
      IF (lun.eq.67) THEN
         hs => hs67
      ELSE
         hs => hs68
      ENDIF
C
      hs%myMesh => adcircMesh
      hs%myFile%fileFound = .false.
C
C     Open fulldomain file
      INQUIRE(FILE=hs%myFile%FILENAME,EXIST=hs%myFile%fileFound)
      IF (hs%myFile%fileFound.eqv..false.) THEN
         WRITE(*,*) "ERROR: The file ",hs%myFile%FILENAME,
     &      " was not found; ADCIRC terminating."
#ifdef CMPI
            CALL MSG_FINI()
#endif
         CALL EXIT(1)
      ELSE
         iret = nf90_open(hs%myFile%FILENAME, NF_NOWRITE, hs%ncid)
         CALL check_err(iret)
      ENDIF
!
!     Inquire variables
      call getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)
C
C     elevation station dimension
      IF (NHASE.ne.0) THEN
         iret=nf90_inq_dimid(hs%ncid,"elevstation",
     &        hs%staelv%num_sta_dim_id)
         CALL check_err(iret)
         iret=nf90_inquire_dimension(hs%ncid, hs%staelv%num_sta_dim_id,
     &        len=hs%staelv%num_stations)
         CALL check_err(iret)
         IF (MNPROC.gt.1) THEN
            ! make a list of full domain elevation stations that
            ! correspond to this subdomain's stations
        allocate(fullDomainElevationStationList(hs%staelv%num_stations))
            ! loop over subdomain indexes to form a list of corresponding
            ! fulldomain indexes
            forall (sd_station_number=1:nstae)
               ! get the corresponding fulldomain indexes
               fullDomainElevationStationList(sd_station_number)
     &          = ABS(imap_stae_lg(sd_station_number))
            end forall
         ENDIF
      ENDIF
C     velocity station dimension
      IF (NHASV.ne.0) THEN
         iret=nf90_inq_dimid(hs%ncid,"velstation",
     &        hs%stavellv%num_sta_dim_id)
         CALL check_err(iret)
         iret=nf90_inquire_dimension(hs%ncid,hs%stavellv%num_sta_dim_id,
     &        len=hs%stavellv%num_stations)
         CALL check_err(iret)
         IF (MNPROC.gt.1) THEN
            ! make a list of full domain velocity stations that
            ! correspond to this subdomain's velocity stations
            allocate(fullDomainVelocityStationList
     &          (hs%stavellv%num_stations))
            ! loop over subdomain indexes to form a list of corresponding
            ! fulldomain indexes
            forall (sd_station_number=1:nstav)
               ! get the corresponding fulldomain indexes
               fullDomainVelocityStationList(sd_station_number)
     &            = ABS(imap_stav_lg(sd_station_number))
            end forall
         ENDIF
      ENDIF
C
C     Point to the hotstart file we want to work on.
      IF (lun.eq.67) THEN
         hs => hs67
      ELSE
         hs => hs68
      ENDIF
C
      hs%myFile%record_counter = 1
C
      kount(1)=MNHARF*2             ! for load vector data
      kount(2)=hs%myMesh%num_nodes  ! for nodal data
      kount(3)=hs%myTime%timenc_len
      start(1)=1
      start(2)=1
      start(3)=hs%myFile%record_counter
C
C     Read in fulldomain load vector data
C     GLOELV - full domain elevation
      fullDomainIndexList => fullDomainNodeList
      IF (NHAGE.ne.0) THEN
         iret=nf90_inq_varid(hs%ncid,"gloelv",hs%gloelv%nodal_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_get_var(hs%ncid, hs%gloelv%nodal_data_id,
     &           GLOELV, start, kount)
            CALL check_err(iret)
         ELSE
            call mapFullDomainToSubdomainMByNP(hs%ncid,
     &         2*MNHARF, hs%myMesh%num_nodes,
     &         hs%gloelv%nodal_data_id, subdomain_reals=GLOELV)
         ENDIF
      ENDIF
C     GLOULV - fulldomain u velocity
C     GLOVLV - fulldomain v velocity
      IF (NHAGV.ne.0) THEN
       iret=nf90_inq_varid(hs%ncid,"gloulv",hs%glovellv%u_nodal_data_id)
         CALL check_err(iret)
       iret=nf90_inq_varid(hs%ncid,"glovlv",hs%glovellv%v_nodal_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_get_var(hs%ncid, hs%glovellv%u_nodal_data_id,
     &          GLOULV, start, kount)
            CALL check_err(iret)
            iret = nf90_get_var(hs%ncid, hs%glovellv%v_nodal_data_id,
     &          GLOVLV, start, kount)
            CALL check_err(iret)
         ELSE
            call mapFullDomainToSubdomainMByNP(hs%ncid,
     &         2*MNHARF, hs%myMesh%num_nodes,
     &         hs%glovellv%u_nodal_data_id, subdomain_reals=GLOULV)
            call mapFullDomainToSubdomainMByNP(hs%ncid,
     &         2*MNHARF, hs%myMesh%num_nodes,
     &         hs%glovellv%v_nodal_data_id, subdomain_reals=GLOVLV)
         ENDIF
      ENDIF
C     STAELV - station elevation
      IF (NHASE.ne.0) THEN
         fullDomainIndexList => fullDomainElevationStationList
         iret=nf90_inq_varid(hs%ncid,"staelv",hs%staelv%station_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            kount(2) = hs%staelv%num_stations ! for elevation stations
            iret = nf90_get_var(hs%ncid, hs%staelv%station_data_id,
     &          STAELV, start, kount)
            CALL check_err(iret)
         ELSE
         ! jgf52.30.12: Fix segfault by checking number of stations
         ! (provided by Damrongsak Wirasaet)
            IF ( NSTAE > 0 ) THEN
               call mapFullDomainToSubdomainMByNP(hs%ncid,
     &            2*MNHARF, hs%staelv%num_stations,
     &            hs%staelv%station_data_id, subdomain_reals=STAELV)
            ENDIF
         ENDIF
      ENDIF
C     STAULV/STAVLV - station u and v velocity

      IF (NHASV.ne.0) THEN

         fullDomainIndexList => fullDomainVelocityStationList
         kount(2) = hs%stavellv%num_stations ! for velocity stations
         iret=nf90_inq_varid(hs%ncid,"staulv",
     &        hs%stavellv%u_station_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(hs%ncid,"stavlv",
     &        hs%stavellv%v_station_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_get_var(hs%ncid, hs%stavellv%u_station_data_id,
     &          STAULV, start, kount)
            CALL check_err(iret)
            iret = nf90_get_var(hs%ncid, hs%stavellv%v_station_data_id,
     &          STAVLV, start, kount)
            CALL check_err(iret)
         ELSE
            ! jgf52.30.12: Fix segfault by checking number of stations
            ! (provided by Damrongsak Wirasaet)
            IF ( NSTAV > 0 ) THEN
               call mapFullDomainToSubdomainMByNP(hs%ncid,
     &            2*MNHARF, hs%stavellv%num_stations,
     &            hs%stavellv%u_station_data_id, subdomain_reals=STAULV)
               call mapFullDomainToSubdomainMByNP(hs%ncid,
     &            2*MNHARF, hs%stavellv%num_stations,
     &            hs%stavellv%v_station_data_id, subdomain_reals=STAVLV)
            ENDIF
         ENDIF
      ENDIF
C
C     Read in model parameters to ADCIRC variables
      iret = nf90_inq_varid(hs%ncid,"icha", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, icha)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nz", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, inz)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nf", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, inf)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"mm", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, imm)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nstae", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, instae)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nstav", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, instav)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nhase", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, inhase)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nhasv", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, inhasv)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nhage", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, inhage)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nhagv", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, inhagv)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"icall", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, icall)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"nfreq", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, infreq)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"timeud", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, timeud)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"itud", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, itud)
      CALL check_err(iret)
C
C     Variables that are used to check that harmonic analysis data match
C     the simulation that is reading this hotstart file.
      inp = hs%myMesh%num_nodes
C
C     left hand side
      hakount(1) = 2*MNHARF
      hakount(2) = 2*MNHARF
      hastart(1) = 1
      hastart(2) = 1
      iret = nf90_inq_varid(hs%ncid,"ha",hs%ha_id)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, hs%ha_id, ha, hastart, hakount)
      CALL check_err(iret)
C     frequency names
      iret = nf90_inq_varid(hs%ncid,"namefr",hs%namefr_id)
      CALL check_err(iret)
      do i=1,mnharf
         start(1)=1
         start(2)=i
         kount(1)= len(inamefr(i))
         kount(2)=1
        iret = nf90_get_var(hs%ncid,hs%namefr_id,inamefr(i),start,kount)
         CALL check_err(iret)
      end do
C     harmonic constituents
      start(1) = 1
      start(2) = 1
      kount(1) = MNHARF ! for constituents
      kount(2) = 1
      iret = nf90_inq_varid(hs%ncid,"hafreq",hs%hafreq_id)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, hs%hafreq_id, ifreq, start, kount)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"haff",hs%haff_id)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, hs%haff_id, iff, start, kount)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"haface",hs%haface_id)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, hs%haface_id, iface, start, kount)
      CALL check_err(iret)
C
C     now close the netcdf file
      iret = nf90_close(hs%ncid)
      CALL check_err(iret)
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE readNetCDFHotstartHarmonic
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E
C        R E A D   N E T C D F   H O T S T A R T   H A R M O N I C
C                                   M E A N S   V A R I A N C E S
C-----------------------------------------------------------------------
C     jgf49.44.11 Reads harmonic analysis data from the hotstart file.
C-----------------------------------------------------------------------
      SUBROUTINE readNetCDFHotstartHarmonicMeansVariances(lun)
      USE SIZES, ONLY : MNPROC, globaldir
      USE GLOBAL, ONLY : OutputDataDescript_t, NODES_LG
      USE HARM, ONLY : ELAV, ELVA, XVELAV, YVELAV, XVELVA, YVELVA,
     &                 NTSTEPS, NHAGE, NHAGV
      USE MESH, ONLY : NP
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      IMPLICIT NONE
C
      INTEGER, intent(in) :: lun
C
      INTEGER i,j
      INTEGER kount(2), start(2)
      INTEGER iret ! success or failure of the netcdf call
      INTEGER tempid
      type(OutputDataDescript_t) :: descript
C
      call setMessageSource("readNetCDFHotstartHarmonicMeansVariances")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
C     Point to the hotstart file we want to work on
      IF (lun.eq.67) THEN
         hs => hs67
      ELSE
         hs => hs68
      ENDIF
C
      hs%myMesh => adcircMesh
      hs%myFile%fileFound = .false.
C
C     Open fulldomain file
      INQUIRE(FILE=hs%myFile%FILENAME,EXIST=hs%myFile%fileFound)
      IF (hs%myFile%fileFound.eqv..false.) THEN
         WRITE(*,*) "ERROR: The file ",hs%myFile%FILENAME,
     &      " was not found; ADCIRC terminating."
#ifdef CMPI
            CALL MSG_FINI()
#endif
         CALL EXIT(1)
      ELSE
         iret = nf90_open(hs%myFile%FILENAME, NF_NOWRITE, hs%ncid)
         CALL check_err(iret)
      ENDIF
!
!     Inquire variables
      call getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)
C
      hs%myFile%record_counter = 1
C
      kount(1)=hs%myMesh%num_nodes  ! for nodal data
      kount(2)=hs%myTime%timenc_len
      start(1)=1
      start(2)=hs%myFile%record_counter
C
      fullDomainIndexList => fullDomainNodeList
      IF (NHAGE.ne.0) THEN
         iret=nf90_inq_varid(hs%ncid,"elav",hs%elav%nodal_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(hs%ncid,"elva",hs%elva%nodal_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_get_var(hs%ncid, hs%elav%nodal_data_id,
     &           ELAV, start, kount)
            CALL check_err(iret)
            iret = nf90_get_var(hs%ncid, hs%elva%nodal_data_id,
     &           ELVA, start, kount)
         CALL check_err(iret)
         ELSE
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes,
     &         hs%elav%nodal_data_id, subdomain_reals=ELAV)
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes,
     &         hs%elva%nodal_data_id, subdomain_reals=ELVA)
         ENDIF
      ENDIF
      IF (NHAGV.ne.0) THEN
         iret=nf90_inq_varid(hs%ncid,"xvelav",hs%xvelav%nodal_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(hs%ncid,"yvelav",hs%yvelav%nodal_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(hs%ncid,"xvelva",hs%xvelva%nodal_data_id)
         CALL check_err(iret)
         iret=nf90_inq_varid(hs%ncid,"yvelva",hs%yvelva%nodal_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_get_var(hs%ncid, hs%xvelav%nodal_data_id,
     &          XVELAV, start, kount)
            CALL check_err(iret)
            iret = nf90_get_var(hs%ncid, hs%yvelav%nodal_data_id,
     &         YVELAV, start, kount)
            CALL check_err(iret)
            iret = nf90_get_var(hs%ncid, hs%xvelva%nodal_data_id,
     &          XVELVA, start, kount)
            CALL check_err(iret)
            iret = nf90_get_var(hs%ncid, hs%yvelva%nodal_data_id,
     &          YVELVA, start, kount)
            CALL check_err(iret)
         ELSE
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes,
     &         hs%xvelav%nodal_data_id, subdomain_reals=XVELAV)
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes,
     &         hs%xvelva%nodal_data_id, subdomain_reals=XVELVA)
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes,
     &         hs%yvelav%nodal_data_id, subdomain_reals=YVELAV)
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes,
     &         hs%yvelva%nodal_data_id, subdomain_reals=YVELVA)
         ENDIF
      ENDIF
C
C     Read in model parameters to ADCIRC variables
      iret = nf90_inq_varid(hs%ncid,"ntsteps",tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, ntsteps)
      CALL check_err(iret)
C
C     now close the netcdf file
      iret = nf90_close(hs%ncid)
      CALL check_err(iret)
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE readNetCDFHotstartHarmonicMeansVariances
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E
C        R E A D   N E T C D F   H O T S T A R T   3 D
C-----------------------------------------------------------------------
C     jgf49.48.03 Reads 3D data from the hotstart file.
C-----------------------------------------------------------------------
      SUBROUTINE readNetCDFHotstart3D(lun)
      USE SIZES, ONLY : MNPROC, globaldir
      USE GLOBAL, ONLY : OutputDataDescript_t, IDEN,
     &   duu1, duv1, dvv1, uu2, vv2, bsx1, bsy1, NP_G
      USE GLOBAL_3DVS, ONLY : q, iy, nfen, sigt, sal, temp, wz, q20,
     &   l, n3dsd, i3dsdrec, n3dsv, i3dsvrec, n3dst, i3dstrec, n3dgd,
     &   i3dgdrec, n3dgv, i3dgvrec, n3dgt, i3dgtrec
      USE MESH, ONLY : NP
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      IMPLICIT NONE
C
      INTEGER, intent(in) :: lun
C
      INTEGER i,j
      INTEGER kount(2), start(2)
      REAL(SZ), ALLOCATABLE :: data2D(:,:)
      INTEGER kount3D(3), start3D(3)
      REAL(SZ), ALLOCATABLE :: data3D(:,:,:)
      REAL(SZ), ALLOCATABLE :: rp(:,:) ! real part of Q (subdomain), i.e. u-vel
      REAL(SZ), ALLOCATABLE :: ip(:,:) ! imag part of Q (subdomain), i.e. v-vel
      INTEGER iret ! success or failure of the netcdf call
      INTEGER tempid
C
      call setMessageSource("readNetCDFHotstart3D")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
C     Point to the hotstart file we want to work on
      IF (lun.eq.67) THEN
         hs => hs67
      ELSE
         hs => hs68
      ENDIF
C
      hs%myMesh => adcircMesh
      hs%myFile%fileFound = .false.
C
C     Open fulldomain file
      INQUIRE(FILE=hs%myFile%FILENAME,EXIST=hs%myFile%fileFound)
      IF (hs%myFile%fileFound.eqv..false.) THEN
         WRITE(*,*) "ERROR: The file ",hs%myFile%FILENAME,
     &      " was not found; ADCIRC terminating."
#ifdef CMPI
            CALL MSG_FINI()
#endif
         CALL EXIT(1)
      ELSE
         iret = nf90_open(hs%myFile%FILENAME, NF_NOWRITE, hs%ncid)
         CALL check_err(iret)
      ENDIF
!
!     Inquire variables
      CALL getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)
C     vertical node dimension
      iret=nf90_inq_dimid(hs%ncid, "num_v_nodes",
     &   hs%myMesh%num_v_nodes_dim_id)
      CALL check_err(iret)
      iret=nf90_inquire_dimension(hs%ncid, hs%myMesh%num_v_nodes_dim_id,
     &     len=hs%myMesh%num_v_nodes)
      CALL check_err(iret)
C
      hs%myFile%record_counter = 1
C
      kount(1)=hs%myMesh%num_nodes  ! for nodal data
      kount(2)=hs%myTime%timenc_len
      start(1)=1
      start(2)=hs%myFile%record_counter
C
      kount3D(1)=hs%myMesh%num_nodes
      kount3D(2)=hs%myMesh%num_v_nodes  ! for 3D data
      kount3D(3)=hs%myTime%timenc_len
      start3D(1)=1
      start3D(2)=1
      start3D(3)=hs%myFile%record_counter

      iret=nf90_inq_varid(hs%ncid,"duu",hs%duu%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"duv",hs%duv%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"dvv",hs%dvv%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"uu",hs%uu%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"vv",hs%vv%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"bsx",hs%bsx%nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"bsy",hs%bsy%nodal_data_id)
      CALL check_err(iret)
      IF (MNPROC.eq.1) THEN
         iret = nf90_get_var(hs%ncid, hs%duu%nodal_data_id,
     &       DUU1, start, kount)
         CALL check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%duv%nodal_data_id,
     &       DUV1, start, kount)
         CALL check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%dvv%nodal_data_id,
     &       DVV1, start, kount)
         CALL check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%uu%nodal_data_id,
     &       UU2, start, kount)
         CALL check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%vv%nodal_data_id,
     &       VV2, start, kount)
         CALL check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%bsx%nodal_data_id,
     &       BSX1, start, kount)
         CALL check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%bsy%nodal_data_id,
     &       BSY1, start, kount)
         CALL check_err(iret)
      ELSE
         fullDomainIndexList => fullDomainNodeList
         call mapFullDomainToSubdomain(hs%ncid,
     &      hs%myMesh%num_nodes, hs%duu%nodal_data_id,
     &      subdomain_reals=DUU1)
         call mapFullDomainToSubdomain(hs%ncid,
     &      hs%myMesh%num_nodes, hs%duv%nodal_data_id,
     &      subdomain_reals=DUV1)
         call mapFullDomainToSubdomain(hs%ncid,
     &      hs%myMesh%num_nodes, hs%dvv%nodal_data_id,
     &      subdomain_reals=DVV1)
         call mapFullDomainToSubdomain(hs%ncid,
     &      hs%myMesh%num_nodes, hs%uu%nodal_data_id,
     &      subdomain_reals=UU2)
         call mapFullDomainToSubdomain(hs%ncid,
     &      hs%myMesh%num_nodes, hs%vv%nodal_data_id,
     &      subdomain_reals=VV2)
         call mapFullDomainToSubdomain(hs%ncid,
     &      hs%myMesh%num_nodes, hs%bsx%nodal_data_id,
     &      subdomain_reals=BSX1)
         call mapFullDomainToSubdomain(hs%ncid,
     &      hs%myMesh%num_nodes, hs%bsy%nodal_data_id,
     &      subdomain_reals=BSY1)
      ENDIF
      !
      ! 3D Density
      IF (ABS(IDEN).eq.1) THEN
        iret=nf90_inq_varid(hs%ncid,"sigt",hs%density3D%u_nodal_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_get_var(hs%ncid,
     &         hs%density3D%u_nodal_data_id, SIGT, start3D, kount3D)
            CALL check_err(iret)
         ELSE
            call mapFullDomainToSubdomainNPByM(hs%ncid,
     &         nfen, hs%myMesh%num_nodes, hs%density3D%u_nodal_data_id,
     &         subdomain_reals=SIGT)
         ENDIF
      ENDIF
      IF ((ABS(IDEN).eq.2).or.(ABS(IDEN).eq.4)) THEN
         iret=nf90_inq_varid(hs%ncid,"salinity",
     &      hs%density3D%v_nodal_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_get_var(hs%ncid, hs%density3D%v_nodal_data_id,
     &          SAL, start3D, kount3D)
            CALL check_err(iret)
         ELSE
            call mapFullDomainToSubdomainNPByM(hs%ncid,
     &         nfen, hs%myMesh%num_nodes, hs%density3D%v_nodal_data_id,
     &         subdomain_reals=SAL)
         ENDIF
      ENDIF
      IF ((ABS(IDEN).eq.3).or.(ABS(IDEN).eq.4)) THEN
         iret=nf90_inq_varid(hs%ncid,"temperature",
     &      hs%density3D%w_nodal_data_id)
         CALL check_err(iret)
         IF (MNPROC.eq.1) THEN
            iret = nf90_get_var(hs%ncid, hs%density3D%w_nodal_data_id,
     &         TEMP, start3D, kount3D)
            CALL check_err(iret)
         ELSE
            call mapFullDomainToSubdomainNPByM(hs%ncid,
     &         nfen, hs%myMesh%num_nodes, hs%density3D%w_nodal_data_id,
     &         subdomain_reals=TEMP)
         ENDIF
      ENDIF
      !
      ! 3D velocity
      iret=nf90_inq_varid(hs%ncid,"u-vel3D",hs%velocity3D%u_nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"v-vel3D",hs%velocity3D%v_nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"w-vel3D",hs%velocity3D%w_nodal_data_id)
      CALL check_err(iret)
      allocate(rp(np,nfen),ip(np,nfen))
      IF (MNPROC.eq.1) THEN
         iret = nf90_get_var(hs%ncid, hs%velocity3D%u_nodal_data_id,
     &      rp, start3D, kount3D)
         CALL check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%velocity3D%v_nodal_data_id,
     &      ip, start3D, kount3D)
         CALL check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%velocity3D%w_nodal_data_id,
     &      wz, start3D, kount3D)
         CALL check_err(iret)
      ELSE
         call mapFullDomainToSubdomainNPByM(hs%ncid,
     &      nfen, hs%myMesh%num_nodes, hs%velocity3D%u_nodal_data_id,
     &      subdomain_reals=rp)
         call mapFullDomainToSubdomainNPByM(hs%ncid,
     &      nfen, hs%myMesh%num_nodes, hs%velocity3D%v_nodal_data_id,
     &      subdomain_reals=ip)
         call mapFullDomainToSubdomainNPByM(hs%ncid,
     &      nfen, hs%myMesh%num_nodes, hs%velocity3D%w_nodal_data_id,
     &      subdomain_reals=wz)
      ENDIF
      q(:,:) = rp(:,:) + iy*ip(:,:) ! construct q from real and imaginary
      deallocate(rp,ip)
      !
      ! 3D turbulence
      iret=nf90_inq_varid(hs%ncid,"q20",hs%turbulence3D%u_nodal_data_id)
      CALL check_err(iret)
      iret=nf90_inq_varid(hs%ncid,"l",hs%turbulence3D%v_nodal_data_id)
      CALL check_err(iret)
      IF (MNPROC.eq.1) THEN
         iret = nf90_get_var(hs%ncid, hs%turbulence3D%u_nodal_data_id,
     &      q20, start3D, kount3D)
         CALL check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%turbulence3D%v_nodal_data_id,
     &      l, start3D, kount3D)
         CALL check_err(iret)
      ELSE
         call mapFullDomainToSubdomainNPByM(hs%ncid,
     &      nfen, hs%myMesh%num_nodes, hs%turbulence3D%u_nodal_data_id,
     &      subdomain_reals=q20)
         call mapFullDomainToSubdomainNPByM(hs%ncid,
     &      nfen, hs%myMesh%num_nodes, hs%velocity3D%v_nodal_data_id,
     &      subdomain_reals=l)
      ENDIF
C
C     Read in model parameters to ADCIRC variables
      iret = nf90_inq_varid(hs%ncid,"n3dsd", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, n3dsd)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"i3dsdrec", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, i3dsdrec)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"n3dsv", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, n3dsv)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"i3dsvrec", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, i3dsvrec)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"n3dst", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, n3dst)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"i3dstrec", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, i3dstrec)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"n3dgd", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, n3dgd)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"i3dgdrec", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, i3dgdrec)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"n3dgv", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, n3dgv)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"i3dgvrec", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, i3dgvrec)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"n3dgt", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, n3dgt)
      CALL check_err(iret)
      iret = nf90_inq_varid(hs%ncid,"i3dgtrec", tempid)
      CALL check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, i3dgtrec)
      CALL check_err(iret)
C
C     now close the netcdf file
      iret = nf90_close(hs%ncid)
      CALL check_err(iret)
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE readNetCDFHotstart3D
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   D E F I N E   M E S H   V A R I A B L E S
C-----------------------------------------------------------------------
C     jgf49.17.02 Defines data that are common to all netCDF files.
C-----------------------------------------------------------------------
      subroutine defineMeshVariables(ncid, myMesh, myFile)
      USE GLOBAL, ONLY : RAD2DEG, C3D, scratchMessage
      USE MESH, ONLY : SLAM0, SFEA0, AGRID
      USE GLOBAL_3DVS, ONLY : NFEN
      USE SIZES, ONLY : SZ
      USE NodalAttributes, ONLY : nwp, nolibf, tau0, cf, eslm
      USE VERSION
      IMPLICIT NONE
      INTEGER iret   ! Error status return
      type(meshStructure) :: myMesh
      INTEGER, intent(in) :: ncid   ! netCDF id
      type(fileData), intent(in) :: myFile ! file format info
      INTEGER i,j
      INTEGER counti(1), starti(1), n
      INTEGER kount(2), start(2)
      REAL    realval(1)
      CHARACTER(1024) :: att_text ! reused to hold attribute text
      INTEGER :: line ! line counter
!     -----------------
!      date_string variables for time attribute
!     -----------------
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer values(8)
C
      call setMessageSource("defineMeshVariables")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
!
!     Define dimensions
      iret = nf90_def_dim(ncid,'node',myMesh%num_nodes,
     &       myMesh%num_nodes_dim_id)
      CALL check_err(iret)
      iret = nf90_def_dim(ncid,'nele',myMesh%num_elems,
     &       myMesh%num_elems_dim_id)
      call check_err(iret)
      IF (C3D.eqv..true.) THEN
         myMesh%num_v_nodes = NFEN
         iret = nf90_def_dim(ncid,'num_v_nodes',myMesh%num_v_nodes,
     &       myMesh%num_v_nodes_dim_id)
         call check_err(iret)
         myMesh%sigma_dims(1) = myMesh%num_v_nodes_dim_id
         iret = nf90_def_var(ncid, 'sigma', NF90_DOUBLE,
     &       myMesh%sigma_dims, myMesh%sigma_id)
         CALL check_err(iret)
      ENDIF
      iret = nf90_def_dim(ncid, 'nvertex', 3, myMesh%nface_dim_id)
      call check_err(iret)
      if(myMesh%nopenc.ne.0) then
         iret = nf90_def_dim(ncid,'nope', myMesh%nopenc,
     &         myMesh%nopenc_dim_id)
         call check_err(iret)
!        WJP to make arrays smaller only have dimension of neta
         iret = nf90_def_dim(ncid,'neta', myMesh%netanc,
     &         myMesh%netanc_dim_id)
         call check_err(iret)
         iret = nf90_def_dim(ncid,'max_nvdll', myMesh%max_nvdllnc,
     &         myMesh%max_nvdllnc_dim_id)
         call check_err(iret)
      endif
      if(myMesh%nbounc.ne.0) then
         iret = nf90_def_dim(ncid, 'nbou',
     &          myMesh%nbounc, myMesh%nbounc_dim_id)
         call check_err(iret)
!        WJP to make arrays smaller only have dimension of nvel
         iret = nf90_def_dim(ncid,'nvel', myMesh%nvelnc,
     &         myMesh%nvelnc_dim_id)
         call check_err(iret)
         iret = nf90_def_dim(ncid, 'max_nvell',
     &          myMesh%max_nvellnc, myMesh%max_nvellnc_dim_id)
      endif
!Corbitt: Define ADCIRC-Mesh Dimension
      iret = nf90_def_dim(ncid,'mesh',1,myMesh%mesh_dim_id)
      call check_err(iret)
!
!     Define variables
!     Define X
      myMesh%X_dims(1) = myMesh%num_nodes_dim_id
      iret = nf90_def_var(ncid, 'x', NF90_DOUBLE,
     &       myMesh%X_dims, myMesh%X_id)
      CALL check_err(iret)
!     Define Y coordinate
      myMesh%Y_dims(1) = myMesh%num_nodes_dim_id
      iret = nf90_def_var(ncid, 'y', NF90_DOUBLE,
     &       myMesh%Y_dims, myMesh%Y_id)
      CALL check_err(iret)
!     Define elements
      myMesh%ELE_dims(1) = myMesh%nface_dim_id
      myMesh%ELE_dims(2) = myMesh%num_elems_dim_id
      iret = nf90_def_var(ncid, 'element',NF90_INT,
     &       myMesh%ELE_dims, myMesh%ELE_id)
      CALL check_err(iret)
!Corbitt: Define ADCIRC-Mesh Variable
      myMesh%MESH_dims(1) = myMesh%mesh_dim_id
      iret = nf90_def_var(ncid, 'adcirc_mesh',NF90_INT,
     &       myMesh%MESH_dims, myMesh%MESH_id)
      CALL check_err(iret)
      !
      ! jgf50.44: Turn on compression if this is a netcdf4-formatted file.
#ifdef NETCDF_CAN_DEFLATE
      IF (myFile%ncformat.eq.ior(NF_CLASSIC_MODEL,NF_NETCDF4)) THEN
         iret = nf90_def_var_deflate(ncid, myMesh%X_id, 1, 1, 2)
         CALL check_err(iret)
         iret = nf90_def_var_deflate(ncid, myMesh%Y_id, 1, 1, 2)
         CALL check_err(iret)
         iret = nf90_def_var_deflate(ncid, myMesh%ELE_id, 1, 1, 2)
         CALL check_err(iret)
      ENDIF
#endif
!
!     Define elevation specified boundary forcing segments information
      if (myMesh%nopenc.ne.0) then
         iret = nf90_def_var(ncid,'neta', NF90_INT,
     &       varid=myMesh%netanc_id)
         CALL check_err(iret)
         myMesh%nvdllnc_dims(1) = myMesh%nopenc_dim_id
         iret = nf90_def_var(ncid, 'nvdll',NF90_INT,
     &                  myMesh%nvdllnc_dims, myMesh%nvdllnc_id)
         CALL check_err(iret)
         iret = nf90_def_var(ncid, 'max_nvdll',NF90_INT,
     &                  varid=myMesh%max_nvdllnc_id)
         CALL check_err(iret)
         myMesh%ibtypeenc_dims(1) = myMesh%nopenc_dim_id
         iret = nf90_def_var(ncid, 'ibtypee',NF90_INT,
     &          myMesh%ibtypeenc_dims,
     &          myMesh%ibtypeenc_id)
         CALL check_err(iret)
!        WJP to make arrays smaller only have dimension of neta
         myMesh%nbdvnc_dims(1) = myMesh%netanc_dim_id
C         myMesh%nbdvnc_dims(1) = myMesh%nopenc_dim_id
C         myMesh%nbdvnc_dims(2) = myMesh%netanc_dim_id
C         myMesh%nbdvnc_dims(2) = myMesh%max_nvdllnc_dim_id
         iret = nf90_def_var(ncid, 'nbdv',NF90_INT,
     &          myMesh%nbdvnc_dims(1), myMesh%nbdvnc_id)
         CALL check_err(iret)
         !jgf50.44: Turn on compression if this is a netcdf4 formatted file.
#ifdef NETCDF_CAN_DEFLATE
         IF (myFile%ncformat.eq.ior(NF_CLASSIC_MODEL,NF_NETCDF4)) THEN
            iret = nf90_def_var_deflate(ncid, myMesh%nvdllnc_id,
     &             1, 1, 2)
            CALL check_err(iret)
            iret = nf90_def_var_deflate(ncid, myMesh%ibtypeenc_id,
     &             1, 1, 2)
            CALL check_err(iret)
            iret = nf90_def_var_deflate(ncid, myMesh%nbdvnc_id,
     &             1, 1, 2)
            CALL check_err(iret)
         ENDIF
#endif
      endif
!     Define normal flow boundary information
      if (myMesh%nbounc.ne.0) then
         iret = nf90_def_var(ncid,'nvel', NF90_INT,
     &      varid=myMesh%nvelnc_id)
         CALL check_err(iret)
         myMesh%nvellnc_dims(1) = myMesh%nbounc_dim_id
         iret = nf90_def_var(ncid, 'nvell',NF90_INT,
     &          myMesh%nvellnc_dims, myMesh%nvellnc_id)
         CALL check_err(iret)
         iret = nf90_def_var(ncid, 'max_nvell',NF90_INT,
     &                  varid=myMesh%max_nvellnc_id)
         CALL check_err(iret)
         myMesh%ibtypenc_dims(1) = myMesh%nbounc_dim_id
         iret = nf90_def_var(ncid, 'ibtype',NF90_INT,
     &          myMesh%ibtypenc_dims, myMesh%ibtypenc_id)
         CALL check_err(iret)
!        WJP to make arrays smaller only have dimension of nvel
         myMesh%nbvvnc_dims(1) = myMesh%nvelnc_dim_id
C         myMesh%nbvvnc_dims(1) = myMesh%nbounc_dim_id
C         myMesh%nbvvnc_dims(2) = myMesh%nvelnc_dim_id
C         myMesh%nbvvnc_dims(2) = myMesh%max_nvellnc_dim_id
         iret = nf90_def_var(ncid, 'nbvv',NF90_INT,
     &          myMesh%nbvvnc_dims(1), myMesh%nbvvnc_id)
         CALL check_err(iret)
         !jgf50.44: Turn on compression if this is a netcdf4 formatted file.
#ifdef NETCDF_CAN_DEFLATE
         IF (myFile%ncformat.eq.ior(NF_CLASSIC_MODEL,NF_NETCDF4)) THEN
            iret = nf90_def_var_deflate(ncid, myMesh%nvellnc_id,
     &             1, 1, 2)
            CALL check_err(iret)
            iret = nf90_def_var_deflate(ncid, myMesh%ibtypenc_id,
     &             1, 1, 2)
            CALL check_err(iret)
            iret = nf90_def_var_deflate(ncid, myMesh%nbvvnc_id,
     &             1, 1, 2)
            CALL check_err(iret)
         ENDIF
#endif
      endif
!     -------------------
!     Define Z coordinate
!     --------------------
      myMesh%DEPTH_dims(1) = myMesh%num_nodes_dim_id
      iret = nf90_def_var(ncid, 'depth',NF90_DOUBLE,
     &       myMesh%DEPTH_dims, myMesh%DEPTH_id)
      CALL check_err(iret)
#ifdef NETCDF_CAN_DEFLATE
      IF (myFile%ncformat.eq.ior(NF_CLASSIC_MODEL,NF_NETCDF4)) THEN
         iret = nf90_def_var_deflate(ncid, myMesh%DEPTH_id,
     &          1, 1, 2)
         CALL check_err(iret)
      ENDIF
#endif
C
C     Set coordinates as representing latitude or longitude, depending on
C     the value of ICS
      CALL defineCoordinateAttributes(ncid, myMesh%X_id, myMesh%Y_id)
!
!     Define depth attributes (Corbitt)
      iret = nf90_put_att(ncid, myMesh%DEPTH_id, 'long_name', 'distance
     &  below geoid')
      CALL check_err(iret)
      iret = nf90_put_att(ncid, myMesh%DEPTH_id, 'standard_name',
     & 'depth below geoid')
      CALL check_err(iret)
      iret = nf90_put_att(ncid, myMesh%DEPTH_id, 'coordinates',
     & 'time y x')
      CALL check_err(iret)
      iret = nf90_put_att(ncid, myMesh%DEPTH_id, 'location', 'node')
      CALL check_err(iret)
      iret = nf90_put_att(ncid, myMesh%DEPTH_id, 'mesh', 'adcirc_mesh')
      CALL check_err(iret)
      iret = nf90_put_att(ncid, myMesh%DEPTH_id, 'units', 'm')
      CALL check_err(iret)
!      iret = nf90_put_att(ncid, myMesh%DEPTH_id, 'positive', 'down')
!      CALL check_err(iret) Corbitt: Not CF Compliant

!     Define Element Attributes (Corbitt)
      iret = nf90_put_att(ncid, myMesh%ELE_id,'long_name', 'element')
      CALL check_err(iret)
!Corbitt 11/12/13 - Obsolete CF-UGRID Convention
!      iret = nf90_put_att(ncid, myMesh%ELE_id,'standard_name',
!     & 'face_node_connectivity')
!      CALL check_err(iret)
      iret = nf90_put_att(ncid, myMesh%ELE_id,'cf_role',
     & 'face_node_connectivity')
      CALL check_err(iret)
      iret = nf90_put_att(ncid, myMesh%ELE_id,'start_index',1)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, myMesh%ELE_id, 'units', 'nondimensional')
      CALL check_err(iret)
!
!     Define text attributes for boundary segments
      if (myMesh%nopenc.ne.0) then
C        NOPE
C         att_text = "number of elevation specified boundary &
C     &forcing segments"
C         iret = nf90_put_att(ncid, myMesh%nopenc_id, 'long_name',
C     &          len_trim(att_text), trim(att_text))
C         CALL check_err(iret)
C         iret = nf90_put_att(ncid, myMesh%nopenc_id, 'units', 14,
C     &                     'nondimensional')
C         CALL check_err(iret)
C        NETA
         att_text = "total number of elevation specified boundary nodes"
         iret = nf90_put_att(ncid, myMesh%netanc_id, 'long_name',
     &          trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(ncid, myMesh%netanc_id, 'units',
     &   'nondimensional')
         CALL check_err(iret)
C        NVDLL
         att_text = "number of nodes in each elevation specified "
     &   //"boundary segment"
         iret = nf90_put_att(ncid, myMesh%nvdllnc_id, 'long_name',
     &          trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(ncid, myMesh%nvdllnc_id, 'units',
     &   'nondimensional')
         CALL check_err(iret)
C        IBTYPEE
         att_text = "elevation boundary type"
         iret = nf90_put_att(ncid, myMesh%ibtypeenc_id, 'long_name',
     &       trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(ncid, myMesh%ibtypeenc_id, 'units',
     &       'nondimensional')
         CALL check_err(iret)
C        NBDV
         att_text = "node numbers on each elevation specified boundary "
     &   //"segment"
         iret = nf90_put_att(ncid, myMesh%nbdvnc_id, 'long_name',
     &          trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(ncid, myMesh%nbdvnc_id, 'units',
     &   'nondimensional')
         CALL check_err(iret)
      endif
      if (myMesh%nbounc.ne.0) then
C        NBOU
C         att_text = "number of normal flow (discharge) specified &
C     &boundary segments"
C         iret = nf90_put_att(ncid, myMesh%nbounc_id, 'long_name',
C     &          len_trim(att_text), trim(att_text))
C         CALL check_err(iret)
C         iret = nf90_put_att(ncid, myMesh%nbounc_id, 'units', 14,
C     &                      'nondimensional')
C         CALL check_err(iret)
C        NVEL
         att_text = "total number of normal flow specified boundary "
     &      //"nodes including both the front and back nodes on "
     &      //"internal barrier boundaries"
         iret = nf90_put_att(ncid, myMesh%nvelnc_id, 'long_name',
     &          trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(ncid, myMesh%nvelnc_id, 'units',
     &   'nondimensional')
         CALL check_err(iret)
C        IBTYPE
         att_text = "type of normal flow (discharge) boundary"
         iret = nf90_put_att(ncid, myMesh%ibtypenc_id, 'long_name',
     &       trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(ncid, myMesh%ibtypenc_id, 'units',
     &   'nondimensional')
         CALL check_err(iret)
C        NVELL
         att_text = 'number of nodes in each normal flow '
     &      //'specified boundary segment'
         iret = nf90_put_att(ncid, myMesh%nvellnc_id, 'long_name',
     &      trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(ncid, myMesh%nvellnc_id, 'units',
     &                      'nondimensional')
         CALL check_err(iret)
C        NBVV
         att_text = "node numbers on normal flow boundary segment"
         iret = nf90_put_att(ncid, myMesh%nbvvnc_id, 'long_name',
     &          trim(att_text))
         CALL check_err(iret)
         iret = nf90_put_att(ncid, myMesh%nbvvnc_id, 'units',
     &                      'nondimensional')
         CALL check_err(iret)
      endif
!Corbitt: Define ADCIRC-Mesh Attributes
         iret = nf90_put_att(ncid, myMesh%mesh_id, 'long_name',
     &                      'mesh_topology')
         CALL check_err(iret)
         iret = nf90_put_att(ncid, myMesh%mesh_id, 'cf_role',
     &                      'mesh_topology')
         CALL check_err(iret)
!Corbitt: 11/12/13 - Obsolete CF-UGRID Convention
!         iret = nf90_put_att(ncid, myMesh%mesh_id, 'standard_name',
!     &                      'mesh_topology')
!         CALL check_err(iret)
!         iret = nf90_put_att(ncid, myMesh%mesh_id, 'dimension',
!     &                      2)
!         CALL check_err(iret)
         iret = nf90_put_att(ncid, myMesh%mesh_id, 'topology_dimension',
     &                      2)
         CALL check_err(iret)
         iret = nf90_put_att(ncid, myMesh%mesh_id, 'node_coordinates',
     &                      'x y')
         CALL check_err(iret)
         iret = nf90_put_att(ncid, myMesh%mesh_id,
     &   'face_node_connectivity','element')
         CALL check_err(iret)

#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE defineMeshVariables
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E
C        D E F I N E  C O O R D I N A T E  A T T R I B U T E S
C-----------------------------------------------------------------------
C     jgf49.21 Defines coordinate attributes to identify coordinates as
C     either x and y in Cartesian length (feet or meters) or latitude
C     and longitude depending on the value of ICS.
C-----------------------------------------------------------------------
      subroutine defineCoordinateAttributes(ncid, xid, yid)
      USE MESH, ONLY : ICS
      IMPLICIT NONE
      INTEGER ncid
      INTEGER xid
      INTEGER yid
      INTEGER iret
C
      call setMessageSource("defineCoordinateAttributes")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      if(ics.ne.1) then  ! this indicates spherical coordinates

C        Define longitude attributes
         iret = nf90_put_att(ncid, xid, 'long_name','longitude')
         CALL check_err(iret)
         iret = nf90_put_att(ncid, xid,'standard_name','longitude')
         CALL check_err(iret)
         iret = nf90_put_att(ncid, xid,'units', 'degrees_east')
         CALL check_err(iret)
         iret = nf90_put_att(ncid, xid, 'positive', 'east')
         CALL check_err(iret)
!        Define latitude attributes
         iret = nf90_put_att(ncid, yid, 'long_name','latitude')
         CALL check_err(iret)
         iret = nf90_put_att(ncid, yid,'standard_name','latitude')
         CALL check_err(iret)
         iret = nf90_put_att(ncid, yid, 'units', 'degrees_north')
         CALL check_err(iret)
         iret = nf90_put_att(ncid, yid, 'positive','north')
         CALL check_err(iret)
      else   ! must be using Cartesian (x,y) coordinates
!        Define x-coordinate attributes
         iret = nf90_put_att(ncid, xid,'long_name',
     &   'Cartesian coordinate x')
         CALL check_err(iret)
         iret = nf90_put_att(ncid, xid, 'standard_name','x_coordinate')
         CALL check_err(iret)
!        determine variable units
         CALL putUnitsAttribute(ncid, xid, 'm', 'ft')
         iret = nf90_put_att(ncid, xid, 'positive', 'right')
         CALL check_err(iret)
!        Define y-coordinate attributes
         iret = nf90_put_att(ncid, yid,'long_name',
     &   'Cartesian coordinate y')
         CALL check_err(iret)
         iret = nf90_put_att(ncid, yid, 'standard_name', 'y_coordinate')
         CALL check_err(iret)
         CALL putUnitsAttribute(ncid, yid, 'm', 'ft')
         iret = nf90_put_att(ncid, yid, 'positive',
     &          '90 degrees counterclockwise from x')
         CALL check_err(iret)
      endif

#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE defineCoordinateAttributes
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   P U T   M E S H   V A R I A B L E S
C-----------------------------------------------------------------------
C     jgf49.17.02 Writes data that are common to all netCDF files
C     to the specified file.
C-----------------------------------------------------------------------
      subroutine putMeshVariables(ncid, myMesh)
      USE SIZES, ONLY : SZ
      USE MESH, ONLY : DP
      USE GLOBAL, ONLY : C3D
      USE GLOBAL_3DVS, ONLY : SIGMA
      USE NodalAttributes, ONLY : nwp, nolibf, tau0, cf, eslm
      USE VERSION
      IMPLICIT NONE
      INTEGER, intent(in) :: ncid
      type(meshStructure), intent(inout) :: myMesh
      INTEGER iret   ! Error status return
C
      INTEGER counti(1), starti(1), n
      INTEGER kount(2), start(2)
      REAL    realval(1)
C
      call setMessageSource("putMeshVariables")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
C     Store nodal coordinates
      iret = nf90_put_var(ncid, myMesh%X_id, myMesh%xnc)
      CALL check_err(iret)
      iret = nf90_put_var(ncid, myMesh%Y_id, myMesh%ync)
      CALL check_err(iret)
      IF (C3D.eqv..true.) THEN
         iret = nf90_put_var(ncid, myMesh%sigma_id, sigma)
         CALL check_err(iret)
      ENDIF
C
!     Store depth
      iret = nf90_put_var(ncid, myMesh%DEPTH_id, DP)
      CALL check_err(iret)
!     Store elements
      kount(1)=myMesh%nface_len
      kount(2)=myMesh%num_elems
      start(1)=1
      start(2)=1
      iret=nf90_put_var(ncid, myMesh%ele_id,
     & myMesh%element,start,kount)
      call check_err(iret)

!     Store elevation boundary information
      if(myMesh%nopenc.ne.0) then
         iret = nf90_put_var(ncid, myMesh%netanc_id, myMesh%netanc)
         CALL check_err(iret)
         iret = nf90_put_var(ncid, myMesh%nvdllnc_id, myMesh%nvdllnc)
         CALL check_err(iret)
         iret = nf90_put_var(ncid, myMesh%max_nvdllnc_id,
     &          myMesh%max_nvdllnc)
         CALL check_err(iret)
         iret = nf90_put_var(ncid,myMesh%ibtypeenc_id,
     &          myMesh%ibtypeenc)
         CALL check_err(iret)
         iret = nf90_put_var(ncid, myMesh%nbdvnc_id ,myMesh%nbdvnc)
         CALL check_err(iret)
      endif
!     Store normal flow boundary information
      if(myMesh%nbounc.ne.0) then
         iret = nf90_put_var(ncid, myMesh%nvelnc_id, myMesh%nvelnc)
         CALL check_err(iret)
         iret = nf90_put_var(ncid,myMesh%ibtypenc_id, myMesh%ibtypenc)
         CALL check_err(iret)
         iret = nf90_put_var(ncid, myMesh%nvellnc_id, myMesh%nvellnc)
         CALL check_err(iret)
         iret = nf90_put_var(ncid, myMesh%max_nvellnc_id,
     &          myMesh%max_nvellnc)
         CALL check_err(iret)
         iret = nf90_put_var(ncid, myMesh%nbvvnc_id, myMesh%nbvvnc)
         CALL check_err(iret)
      endif

#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE putMeshVariables
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   D E F I N E   T I M E   A T T R I B U T E S
C-----------------------------------------------------------------------
C     jgf50.38 Writes time data that are common to all netCDF files
C     to the specified file.
C-----------------------------------------------------------------------
      SUBROUTINE defineTimeAttributes(ncid, myTime)
      USE GLOBAL, ONLY : base_date
      IMPLICIT NONE
      INTEGER, intent(in) :: ncid  ! file id of netcdf file
      type(TimeData), intent(in) :: myTime ! time data for this netcdf file
C
      character date_string*60
      integer :: iret    ! error status of netcdf call
C
      iret = nf90_put_att(ncid, myTime%timenc_id,
     &       'long_name','model time')
      CALL check_err(iret)
      iret = nf90_put_att(ncid, myTime%timenc_id,
     &       'standard_name', 'time')
      CALL check_err(iret)
      date_string = 'seconds since '//adjustl(trim(base_date))
      iret = nf90_put_att(ncid, myTime%timenc_id, 'units',
     &        trim(date_string))
      CALL check_err(iret)
      iret = nf90_put_att(ncid, myTime%timenc_id,
     &       'base_date', adjustl(trim(base_date)))
      CALL check_err(iret)
C-----------------------------------------------------------------------
      END SUBROUTINE defineTimeAttributes
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   D E F I N E   M E T A D A T A
C-----------------------------------------------------------------------
C     jgf49.17.02 Writes data that are common to all netCDF files
C     to the specified file.
C-----------------------------------------------------------------------
      subroutine defineMetaData(ncid)
      USE VERSION
      USE MESH, ONLY : ICS, AGRID, SLAM0, SFEA0
      USE GLOBAL, ONLY : RUNDES, title, institution, source,
     &                   history, references, comments, host,
     &                   convention, contact, dtdp, ihot,
     &                   nolifa, nolica, nolicat,
     &                   ncor, ntip, nws, nramp, statim,
     &                   reftim, rnday, dramp, h0,
     &                   cori, ntif, nbfr, RAD2DEG, C3D, inunThresh,
     &                   inundationOutput, runid
      USE GLOBAL_3DVS, ONLY : iden, islip, kp, z0s, z0b, theta1, theta2,
     &   ievc, evmin, evcon, alp1, alp2, alp3, igc, nlsd, nvsd, nltd,
     &   nvtd, alp4
      USE GWCE, ONLY : a00, b00, c00
      USE NodalAttributes, ONLY: nolibf, nwp, tau0, cf, eslm
      IMPLICIT NONE
      INTEGER :: ncid
      INTEGER :: iret ! success or failure of the netcdf call
C
      REAL(8) SLAM0DEG
      REAL(8) SFEA0DEG
!     date_string variables for time attribute
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer values(8)
C
      call setMessageSource("defineMetaData")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
!
!     Convert back to degrees ... the original input is in degrees,
C     but this gets converted to radians immediately and unfortunately
C     the values that were read in get overwritten ... need to go
C     back to degrees to write them back out
      SLAM0DEG=SLAM0*RAD2DEG
      SFEA0DEG=SFEA0*RAD2DEG
!     -----------------
!     Global attributes
!     -----------------
      iret = nf90_put_att(ncid, NF_GLOBAL, 'model', 'ADCIRC')
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NF_GLOBAL, 'version', ADC_VERSION)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NF_GLOBAL, 'grid_type', 'Triangular')
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NF_GLOBAL, 'description',
     &       trim(adjustl(rundes)))
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NF_GLOBAL, 'agrid',
     &        trim(adjustl(agrid)))
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NF_GLOBAL, 'rundes',
     &       trim(adjustl(rundes)))
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NF_GLOBAL, 'runid',
     &       trim(adjustl(runid)))
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NF_GLOBAL, 'title',
     &        trim(adjustl(title)))
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NF_GLOBAL, 'institution',
     &        trim(adjustl(institution)))
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NF_GLOBAL, 'source',
     &       trim(adjustl(source)))
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NF_GLOBAL, 'history',
     &       trim(adjustl(history)))
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NF_GLOBAL, 'references',
     &       trim(adjustl(references)))
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NF_GLOBAL, 'comments',
     &       trim(adjustl(comments)))
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NF_GLOBAL, 'host',
     &       trim(adjustl(host)))
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NF_GLOBAL, 'convention',
     &       trim(adjustl(convention)))
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NF_GLOBAL, 'Conventions',
     &       'UGRID-0.9.0')
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NF_GLOBAL, 'contact',
     &       trim(adjustl(contact)))
      CALL check_err(iret)
      CALL date_and_time(now_date,big_ben,zone,values)
      WRITE(date_string,71) values(1),values(2),values(3),
     &      values(5),values(6),values(7),   (values(4))/60
 71   FORMAT(I4,'-',I2.2,'-',i2.2,' ',i2,':',i2.2,':',i2.2,' '
     &  ,i3.2,':00')
      iret = nf90_put_att(ncid, NF_GLOBAL,'creation_date',
     &       trim(date_string))
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NF_GLOBAL,'modification_date',
     &       trim(date_string))
      CALL check_err(iret)

!     -------------------------------------------
!     writing global attributes from fort.15 file
!     -------------------------------------------
      iret = nf90_put_att(ncid, NF_GLOBAL, 'fort.15',
     &                      '==== Input File Parameters (below) ====')
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'dt', dtdp)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'ihot', ihot)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'ics', ics)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'nolibf', nolibf)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'nolifa', nolifa)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'nolica', nolica)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'nolicat', nolicat)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'nwp', nwp)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'ncor', ncor)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'ntip', ntip)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'nws', nws)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'nramp', nramp)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'tau0', tau0)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'statim', statim)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'reftim', reftim)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'rnday', rnday)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'dramp', dramp)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'a00', a00)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'b00', b00)

      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'c00', c00)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'h0', h0)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'slam0', slam0deg)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'sfea0', sfea0deg)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'cf', cf)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'eslm', eslm)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'cori', cori)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'ntif', ntif)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'nbfr', nbfr)
      CALL check_err(iret)
      IF (C3D.eqv..true.) THEN
         iret = nf90_put_att(ncid, NCGLOBAL,'iden', iden)
         CALL check_err(iret)
         iret = nf90_put_att(ncid, NCGLOBAL,'islip', islip)
         CALL check_err(iret)
         iret = nf90_put_att(ncid, NCGLOBAL,'kp', kp)
         CALL check_err(iret)
         iret = nf90_put_att(ncid, NCGLOBAL,'z0s', z0s)
         CALL check_err(iret)
         iret = nf90_put_att(ncid, NCGLOBAL,'z0b', z0b)
         CALL check_err(iret)
         iret = nf90_put_att(ncid, NCGLOBAL,'alp1', alp1)
         CALL check_err(iret)
         iret = nf90_put_att(ncid, NCGLOBAL,'alp2', alp2)
         CALL check_err(iret)
         iret = nf90_put_att(ncid, NCGLOBAL,'alp3', alp3)
         CALL check_err(iret)
         iret = nf90_put_att(ncid, NCGLOBAL,'igc', igc)
         CALL check_err(iret)
         iret = nf90_put_att(ncid, NCGLOBAL,'ievc', ievc)
         CALL check_err(iret)
         iret = nf90_put_att(ncid, NCGLOBAL,'evmin', evmin)
         CALL check_err(iret)
         iret = nf90_put_att(ncid, NCGLOBAL,'evcon', evcon)
         CALL check_err(iret)
         IF ((ievc.eq.50).or.(ievc.eq.51)) THEN
            iret = nf90_put_att(ncid, NCGLOBAL,'theta1', theta1)
            CALL check_err(iret)
            iret = nf90_put_att(ncid, NCGLOBAL,'theta2', theta2)
            CALL check_err(iret)
         ENDIF
         IF (iden.gt.0) THEN
            iret = nf90_put_att(ncid, NCGLOBAL,'nlsd', nlsd)
            CALL check_err(iret)
            iret = nf90_put_att(ncid, NCGLOBAL,'nvsd', nvsd)
            CALL check_err(iret)
            iret = nf90_put_att(ncid, NCGLOBAL,'nltd', nltd)
            CALL check_err(iret)
            iret = nf90_put_att(ncid, NCGLOBAL,'nvtd', nvtd)
            CALL check_err(iret)
            iret = nf90_put_att(ncid, NCGLOBAL,'alp4', alp4)
            CALL check_err(iret)
         ENDIF
      ENDIF
      if (inundationOutput.eqv..true.) then
         iret = nf90_put_att(ncid, NCGLOBAL,'inunThresh', inunThresh)
         CALL check_err(iret)
      endif

#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE defineMetaData
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     F U N C T I O N     M E T A   L E N G T H
C-----------------------------------------------------------------------
C     jgf49.29 Gets the length of the metadata line, depending on the
C     presence of a "!" in the line; if a "!" is present, it is used to
C     terminate the metadata ... if not, then the whole line is used (up
C     to 80 characters, or as declared in global.F).
C-----------------------------------------------------------------------
      function metalength(string)
      integer :: metalength
      character(*), intent(in) :: string
      metalength = index(string,"!") ! use the "!" as terminator if present
      if (metalength.eq.0) then
         ! there is no embedded "!" in the metadata line -- use the full line
         metalength = len_trim(string)
      else
         ! trim space between end of metadata and embedded "!"
         metalength = len_trim(string(1:metalength-1))
      endif
C-----------------------------------------------------------------------
      end function metalength
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   U P D A T E   M E T A   D A T A
C-----------------------------------------------------------------------
C     jgf49.17.02 Updates data that are common to all netCDF files
C     in the specified file.
C-----------------------------------------------------------------------
      subroutine updateMetaData(ncid,myFile)
      USE SIZES, ONLY : SZ
      USE GLOBAL, ONLY : ihot, nramp, rnday, dramp, h0, rundes, runid
      USE NodalAttributes, ONLY : nwp, nolibf, tau0, cf, eslm
      USE GWCE, ONLY : a00, b00, c00
      USE VERSION
      IMPLICIT NONE
      INTEGER :: ncid
      type(fileData), intent(inout) :: myFile
C
      INTEGER iret   ! Error status return
!     date_string variables for time attribute
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer values(8)
C
      call setMessageSource("updateMetaData")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C     Open existing NetCDF file
      iret = nf90_open(myFile%FILENAME, NF_WRITE, ncid)
      CALL check_err(iret)
C
      iret = nf90_redef (ncid)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'rundes', trim(adjustl(rundes)))
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'runid', trim(adjustl(runid)))
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'ihot', ihot)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'nramp', nramp)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'tau0', tau0)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'rnday', rnday)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'dramp', dramp)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'a00', a00)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'b00', b00)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'c00', c00)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'h0', h0)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'cf', cf)
      CALL check_err(iret)
      iret = nf90_put_att(ncid, NCGLOBAL, 'eslm', eslm)
      CALL check_err(iret)
      CALL date_and_time(now_date,big_ben,zone,values)
      WRITE(date_string,71) values(1),values(2),values(3),
     &      values(5),values(6),values(7),   (values(4))/60
 71   FORMAT(I4,'-',I2.2,'-',i2.2,' ',i2,':',i2.2,':',i2.2,' '
     &  ,i3.2,':00')
      iret = nf90_put_att(ncid, NF_GLOBAL,'modification_date',
     &       date_string)
      CALL check_err(iret)
      iret = nf90_enddef (ncid)
      CALL check_err(iret)
C     now close the updated netcdf file
      iret = nf90_close(ncid)
      CALL check_err(iret)
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE updateMetaData
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E  R E A D  M E T A  D A T A
C-----------------------------------------------------------------------
C     zc50.74.05 - Give the writer processors a way to read the metadata
C     they need to update in the NetCDF file.
C-----------------------------------------------------------------------
      SUBROUTINE ReadMetaData(ncid,myFile)
C-----------------------------------------------------------------------
        USE SIZES,ONLY:SZ
        USE GLOBAL,ONLY:ihot,nramp,rnday,dramp,h0,np_g,ne_g,rundes,runid
        USE NodalAttributes,ONLY:nwp,nolibf,tau0,cf,eslm
        USE VERSION
        USE GWCE, ONLY : a00,b00,c00
        IMPLICIT NONE
        INTEGER :: NCID,dimid_node,dimid_nele
        TYPE(fileData),INTENT(INOUT) :: myFile
        INTEGER :: iret

        call setMessageSource("readMetaData")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif
        CALL CHECK_ERR(nf90_open(myFile%FILENAME,NF_NOWRITE,NCID))
        CALL CHECK_ERR(nf90_get_att(NCID,NCGLOBAL,'rundes',rundes))
        CALL CHECK_ERR(nf90_get_att(NCID,NCGLOBAL,'runid',runid))
        CALL CHECK_ERR(nf90_get_att(NCID,NCGLOBAL,'ihot',ihot))
        CALL CHECK_ERR(nf90_get_att(NCID,NCGLOBAL,'nramp',nramp))
        CALL CHECK_ERR(nf90_get_att(NCID,NCGLOBAL,'tau0',tau0))
        CALL CHECK_ERR(nf90_get_att(NCID,NCGLOBAL,'rnday',rnday))
        CALL CHECK_ERR(nf90_get_att(NCID,NCGLOBAL,'dramp',dramp))
        CALL CHECK_ERR(nf90_get_att(NCID,NCGLOBAL,'a00',a00))
        CALL CHECK_ERR(nf90_get_att(NCID,NCGLOBAL,'b00',b00))
        CALL CHECK_ERR(nf90_get_att(NCID,NCGLOBAL,'c00',c00))
        CALL CHECK_ERR(nf90_get_att(NCID,NCGLOBAL,'h0',h0))
        CALL CHECK_ERR(nf90_get_att(NCID,NCGLOBAL,'cf',cf))
        CALL CHECK_ERR(nf90_get_att(NCID,NCGLOBAL,'eslm',eslm))
        CALL CHECK_ERR(nf90_inq_dimid(NCID,"node",dimid_node))
        CALL CHECK_ERR(nf90_inq_dimid(NCID,"nele",dimid_nele))
        CALL CHECK_ERR(nf90_inquire_dimension(ncid,dimid_node,len=NP_G))
        CALL CHECK_ERR(nf90_inquire_dimension(ncid,dimid_nele,len=NE_G))
        CALL CHECK_ERR(nf90_close(NCID))
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
        CALL allMessage(DEBUG,"Enter.")
#endif
        CALL unsetMessageSource()
        RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   P U T   U N I T S   A T T R I B U T E
C-----------------------------------------------------------------------
C     jgf49.17.02 Puts the right units label based on whether ADCIRC was
C     run with English units or SI units.
C-----------------------------------------------------------------------
      subroutine putUnitsAttribute(ncid, var_id, metric, english)
      USE GLOBAL, ONLY : G
      IMPLICIT NONE
      INTEGER ncid
      INTEGER var_id
      INTEGER iret  ! success or failure of netcdf call
      CHARACTER(*) metric
      CHARACTER(*) english
C
      call setMessageSource("putUnitsAttribute")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      if (G.LT.11.D0) then
          iret = nf90_put_att(ncid, var_id, 'units', metric)
      else
          iret = nf90_put_att(ncid, var_id, 'units', english)
      endif
      CALL check_err(iret)

#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C
C-----------------------------------------------------------------------
      end subroutine putUnitsAttribute
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   C H E C K  _  E R R
C-----------------------------------------------------------------------
C     jgf49.17.02 Checks the return value from netCDF calls; if there
C     was an error, it writes the error message to the screen and to the
C     fort.16 file.
C-----------------------------------------------------------------------
      subroutine check_err(iret)
      USE SIZES, ONLY : SZ, myproc
      USE GLOBAL, ONLY : screenUnit
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      IMPLICIT NONE
      INTEGER, intent(in) :: iret

      call setMessageSource("check_err")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      if (iret .ne. NF_NOERR) then
         call allMessage(ERROR,nf90_strerror(iret))
         call netcdfTerminate()
      endif
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      end subroutine check_err
C-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!     S U B R O U T I N E   N E T C D F   T E R M I N A T E
!-----------------------------------------------------------------------
      SUBROUTINE netcdfTerminate(NO_MPI_FINALIZE)
#ifdef CMPI
      USE MESSENGER
#endif
      USE GLOBAL, ONLY : setMessageSource, unsetMessageSource,
     &   allMessage, DEBUG, ECHO, INFO, WARNING, ERROR
      IMPLICIT NONE
      LOGICAL, OPTIONAL :: NO_MPI_FINALIZE
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      REAL, ALLOCATABLE :: dummy(:)
#endif
C
      call setMessageSource("netcdfTerminate")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      call allMessage(INFO,"ADCIRC Terminating.")

#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
         ! intentionally create a segmentation fault so that we can get
         ! a stack trace to determine the line number of the netcdf call
         ! that went bad ... this assumes that the code was compiled with
         ! debugging symbols, bounds checking, and stack trace turned on.
         dummy(1) = 99.9d0
#endif

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
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      call unsetMessageSource()
!-----------------------------------------------------------------------
      END SUBROUTINE netcdfTerminate
!-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E   S E T   A D C I R C   P A R A M E T E R S
C-----------------------------------------------------------------------
C     jgf49.31 Called by adcprep to populate the global.F data normally
C     collected by read_input.F in adcirc. This is somewhat
C     convoluted, since the data are collected by adcprep in read14()
C     and read15() and used to populate pre_global.F, and then we
C     call this subroutine in the netcdf module to populate the global
C     and nodalattributes modules with the same data.
C
C     This twisted logic is the result of the unfortunate overlap between the
C     pre_global and global modules, among other things. Once adcprep
C     becomes integrated into adcirc, this silly subroutine will no longer
C     be needed. At the very least, adcprep should be made to populate
C     global, not pre_global with these data, in other words, both
C     adcirc and adcprep should use the global module.
C-----------------------------------------------------------------------
      SUBROUTINE setADCIRCParameters(
     &   base_date_p, NE_G_p,
     &   NBOU_p, NVEL_p, NOPE_p, NP_G_p, SLAM0_p, SFEA0_p, NBVV_p,
     &   NVDLL_p, NBDV_p, NVELL_p, X_p, Y_p,
     &   IBTYPE_p, IBTYPEE_p, SLAM_p, SFEA_p,
     &   NODECODE_p, G_p, FileFmtRev_p, FileFmtMinor_p,
     &   FileFmtMajor_p, im_p, iestp_p, nscoue_p, ivstp_p, nscouv_p,
     &   icstp_p, nscouc_p, ipstp_p, iwstp_p, nscoum_p, igep_p,
     &   nscouge_p, igvp_p, nscougv_p, igcp_p, nscougc_p, igpp_p,
     &   igwp_p, nscougw_p, NM_p, DP_p, RUNDES_p, AGRID_p, title_p,
     &   institution_p, source_p, history_p, references_p, comments_p,
     &   host_p, convention_p, contact_p, dtdp_p, ihot_p, ics_p,
     &   nolifa_p, nolica_p, nolicat_p, ncor_p, ntip_p, nws_p, nramp_p,
     &   statim_p, reftim_p, rnday_p, dramp_p, a00_p, b00_p, c00_p,
     &   h0_p, cori_p, ntif_p, nbfr_p, myProc_p, screenUnit_p, nolibf_p,
     &   nwp_p, tau0_p, cf_p, eslm_p, nstae_p, nstav_p, nstam_p, neta_p,
     &   nabout_p, nscreen_p,
     &   nfen_p,iden_p, islip_p, kp_p, z0s_p, z0b_p, theta1_p, theta2_p,
     &   ievc_p, evmin_p, evcon_p, alp1_p, alp2_p, alp3_p, igc_p,nlsd_p,
     &   nvsd_p, nltd_p, nvtd_p, alp4_p, C3D_p, RUNID_p)
      USE SIZES, ONLY : myproc
      USE MESH, ONLY : X, Y, SLAM, SFEA, NM, DP, ICS,
     &   SLAM0, SFEA0, AGRID
      USE BOUNDARIES, ONLY : NBOU, NVEL, NOPE, NBVV, NVDLL, NBDV, NVELL,
     &   IBTYPE, IBTYPEE, NETA
      USE VERSION, ONLY : FileFmtRev, FileFmtMinor, FileFmtMajor
      USE GLOBAL, ONLY : base_date, SNAMLEN, RAD2DEG, NODECODE, G,
     &   NP_G, NE_G, im, iestp, nscoue, ivstp, nscouv,
     &   icstp, nscouc, ipstp, iwstp, nscoum, igep, nscouge, igvp,
     &   nscougv, igcp, nscougc, igpp, igwp, nscougw,
     &   RUNDES, title, institution, source, history,
     &   references, comments, host, convention, contact, dtdp, ihot,
     &   nolifa, nolica, nolicat, ncor, ntip, nws, nramp, statim,
     &   reftim, rnday, dramp, h0, cori, ntif, nbfr,
     &   screenUnit, nabout, nscreen, C3D, runid
      USE GLOBAL_3DVS, ONLY :
     &   nfen, n3dsd, i3dsdrec, n3dsv, i3dsvrec,
     &   n3dst, i3dstrec, n3dgd, i3dgdrec, n3dgv, i3dgvrec,
     &   n3dgt, i3dgtrec,
     &   iden, islip, kp, z0s, z0b, theta1, theta2,
     &   ievc, evmin, evcon, alp1, alp2, alp3, igc, nlsd, nvsd, nltd,
     &   nvtd, alp4
      USE GWCE, ONLY : a00, b00, c00
      USE NodalAttributes, ONLY : nolibf, nwp, tau0, cf, eslm
      IMPLICIT NONE
C
C     Declare the argument variables coming in from adcprep.
      CHARACTER(80), intent(in) :: base_date_p
      INTEGER, intent(in) :: NE_G_p
      INTEGER, intent(in) :: NBOU_p
      INTEGER, intent(in) :: NVEL_p
      INTEGER, intent(in) :: NOPE_p
      INTEGER, intent(in) :: NP_G_p
      REAL(8), intent(in) :: SLAM0_p
      REAL(8), intent(in) :: SFEA0_p
      INTEGER, intent(in) :: NBVV_p(:,:)
      INTEGER, intent(in) :: NVDLL_p(:)
      INTEGER, intent(in) :: NBDV_p(:,:)
      INTEGER, intent(in) :: NVELL_p(:)
      REAL(8), intent(in) :: X_p(:)
      REAL(8), intent(in) :: Y_p(:)
      INTEGER, intent(in) :: IBTYPE_p(:)
      INTEGER, intent(in) :: IBTYPEE_p(:)
      REAL(8), intent(in) :: SLAM_p(:)
      REAL(8), intent(in) :: SFEA_p(:)
      INTEGER, intent(in) :: NODECODE_p(:)
      REAL(SZ), intent(in) :: G_p
      INTEGER, intent(in) :: FileFmtRev_p
      INTEGER, intent(in) :: FileFmtMinor_p
      INTEGER, intent(in) :: FileFmtMajor_p
      INTEGER, intent(in) :: im_p
      INTEGER, intent(in) :: iestp_p
      INTEGER, intent(in) :: nscoue_p
      INTEGER, intent(in) :: ivstp_p
      INTEGER, intent(in) :: nscouv_p
      INTEGER, intent(in) :: icstp_p
      INTEGER, intent(in) :: nscouc_p
      INTEGER, intent(in) :: ipstp_p
      INTEGER, intent(in) :: iwstp_p
      INTEGER, intent(in) :: nscoum_p
      INTEGER, intent(in) :: igep_p
      INTEGER, intent(in) :: nscouge_p
      INTEGER, intent(in) :: igvp_p
      INTEGER, intent(in) :: nscougv_p
      INTEGER, intent(in) :: igcp_p
      INTEGER, intent(in) :: nscougc_p
      INTEGER, intent(in) :: igpp_p
      INTEGER, intent(in) :: igwp_p
      INTEGER, intent(in) :: nscougw_p
      INTEGER, intent(in) :: NM_p(:,:)
      REAL(SZ), intent(in) :: DP_p(:)
      CHARACTER(80), intent(in) :: RUNDES_p
      CHARACTER(80), intent(in) :: RUNID_p
      CHARACTER(80), intent(in) :: AGRID_p
      CHARACTER(80), intent(in) :: title_p
      CHARACTER(80), intent(in) :: institution_p
      CHARACTER(80), intent(in) :: source_p
      CHARACTER(80), intent(in) :: history_p
      CHARACTER(80), intent(in) :: references_p
      CHARACTER(80), intent(in) :: comments_p
      CHARACTER(80), intent(in) :: host_p
      CHARACTER(80), intent(in) :: convention_p
      CHARACTER(80), intent(in) :: contact_p
      REAL(8), intent(in) :: dtdp_p
      INTEGER, intent(in) :: ihot_p
      INTEGER, intent(in) :: ics_p
      INTEGER, intent(in) :: nolifa_p
      INTEGER, intent(in) :: nolica_p
      INTEGER, intent(in) :: nolicat_p
      INTEGER, intent(in) :: ncor_p
      INTEGER, intent(in) :: ntip_p
      INTEGER, intent(in) :: nws_p
      INTEGER, intent(in) :: nramp_p
      REAL(8), intent(in) :: statim_p
      REAL(8), intent(in) :: reftim_p
      REAL(SZ), intent(in) :: rnday_p
      REAL(SZ), intent(in) :: dramp_p
      REAL(SZ), intent(in) :: a00_p
      REAL(SZ), intent(in) :: b00_p
      REAL(SZ), intent(in) :: c00_p
      REAL(SZ), intent(in) :: h0_p
      REAL(SZ), intent(in) :: cori_p
      INTEGER, intent(in) :: ntif_p
      INTEGER, intent(in) :: nbfr_p
      INTEGER, intent(in) :: myProc_p
      INTEGER, intent(in) :: screenUnit_p
      INTEGER, intent(in) :: nolibf_p
      INTEGER, intent(in) :: nwp_p
      REAL(SZ), intent(in) :: tau0_p
      REAL(SZ), intent(in) :: cf_p
      REAL(SZ), intent(in) :: eslm_p
      INTEGER, intent(in) :: nstae_p
      INTEGER, intent(in) :: nstav_p
      INTEGER, intent(in) :: nstam_p
      INTEGER, intent(in) :: neta_p
      INTEGER, intent(in) :: nabout_p
      INTEGER, intent(in) :: nscreen_p
C
      INTEGER, intent(in) :: nfen_p
      INTEGER, intent(in) :: iden_p
      INTEGER, intent(in) :: islip_p
      REAL(SZ), intent(in) :: kp_p
      REAL(SZ), intent(in) :: z0s_p
      REAL(SZ), intent(in) :: z0b_p
      REAL(SZ), intent(in) :: theta1_p
      REAL(SZ), intent(in) :: theta2_p
      INTEGER, intent(in) :: ievc_p
      REAL(SZ), intent(in) :: evmin_p
      REAL(SZ), intent(in) :: evcon_p
      REAL(SZ), intent(in) :: alp1_p
      REAL(SZ), intent(in) :: alp2_p
      REAL(SZ), intent(in) :: alp3_p
      INTEGER, intent(in) :: igc_p
      REAL(SZ), intent(in) :: nlsd_p
      REAL(SZ), intent(in) :: nvsd_p
      REAL(SZ), intent(in) :: nltd_p
      REAL(SZ), intent(in) :: nvtd_p
      REAL(SZ), intent(in) :: alp4_p
      LOGICAL, intent(in) :: C3D_p
C
      call setMessageSource("setADCIRCParameters")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      base_date = base_date_p
      NE_G = NE_G_p
      NBOU = NBOU_p
      NVEL = NVEL_p
      NOPE = NOPE_p
      NP_G = NP_G_p
      SLAM0 = SLAM0_p
      SFEA0 = SFEA0_p
      ALLOCATE(NBVV(NBOU_p,0:NVEL_p))
      NBVV = NBVV_p
      ALLOCATE(NVDLL(NOPE_p))
      NVDLL = NVDLL_p
      ALLOCATE(NBDV(NOPE_p,NETA_p))
      NBDV = NBDV_p
      ALLOCATE(NVELL(NBOU_p))
      NVELL = NVELL_p
      ALLOCATE(X(NP_G_p))
      X = X_p
      ALLOCATE(Y(NP_G_p))
      Y = Y_p
      ALLOCATE(IBTYPEE(NOPE_p))
      IBTYPEE = IBTYPEE_p
      ALLOCATE(IBTYPE(NBOU_p))
      IBTYPE = IBTYPE_p
      ALLOCATE(SLAM(NP_G_p))
      SLAM = SLAM_p
      ALLOCATE(SFEA(NP_G_p))
      SFEA = SFEA_p
      ALLOCATE(NODECODE(NP_G_p))
      NODECODE = NODECODE_p
      G = G_p
      FileFmtRev = FileFmtRev_p
      FileFmtMinor = FileFmtMinor_p
      FileFmtMajor = FileFmtMajor_p
      im = im_p
      iestp = iestp_p
      nscoue = nscoue_p
      ivstp = ivstp_p
      nscouv = nscouv_p
      icstp = icstp_p
      nscouc = nscouc_p
      ipstp = ipstp_p
      iwstp = iwstp_p
      nscoum = nscoum_p
      igep = igep_p
      nscouge = nscouge_p
      igvp = igvp_p
      nscougv = nscougv_p
      igcp = igcp_p
      nscougc = nscougc_p
      igpp = igpp_p
      igwp = igwp_p
      nscougw = nscougw_p
      ALLOCATE(NM(NE_G_p,3))
      NM = NM_p
      ALLOCATE(DP(NP_G_p))
      DP = DP_p
      RUNDES = RUNDES_p
      RUNID = RUNID_p
      AGRID = AGRID_p
      title = title_p
      institution = institution_p
      source = source_p
      history = history_p
      references = references_p
      comments = comments_p
      host = host_p
      convention = convention_p
      contact = contact_p
      dtdp = dtdp_p
      ihot = ihot_p
      ics = ics_p
      nolifa = nolifa_p
      nolica = nolica_p
      nolicat = nolicat_p
      ncor = ncor_p
      ntip = ntip_p
      nws = nws_p
      nramp = nramp_p
      statim = statim_p
      reftim = reftim_p
      rnday = rnday_p
      dramp = dramp_p
      a00 = a00_p
      b00 = b00_p
      c00 =  c00_p
      h0 = h0_p
      cori = cori_p
      ntif = ntif_p
      nbfr = nbfr_p
      myProc = myProc_p
      screenUnit = screenUnit_p
      nolibf = nolibf_p
      nwp = nwp_p
      tau0 = tau0_p
      cf = cf_p
      eslm = eslm_p
      neta = neta_p
      nabout = nabout_p
      nscreen = nscreen_p
C
      nfen = nfen_p
      iden = iden_p
      islip = islip_p
      kp = kp_p
      z0s = z0s_p
      z0b = z0b_p
      theta1 = theta1_p
      theta2 = theta2_p
      ievc = ievc_p
      evmin = evmin_p
      evcon = evcon_p
      alp1 = alp1_p
      alp2 = alp2_p
      alp3 = alp3_p
      igc = igc_p
      nlsd = nlsd_p
      nvsd = nvsd_p
      nltd = nltd_p
      nvtd = nvtd_p
      alp4 = alp4_p
      C3D = C3D_p
C
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE setADCIRCParameters
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   F R E E    N E T C D F   C O O R D
C-----------------------------------------------------------------------
C     jgf49.43.05 Frees memory allocated for mesh and boundary data
C     in NetCDF files.
C-----------------------------------------------------------------------
      SUBROUTINE freeNetCDFCoord()
      IMPLICIT NONE
C
      call setMessageSource("freeNetCDFCoord.")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      IF (adcircMesh%initialized.eqv..true.) THEN
         DEALLOCATE(adcircMesh%xnc)
         DEALLOCATE(adcircMesh%ync)
         DEALLOCATE(adcircMesh%nvdllnc)
         DEALLOCATE(adcircMesh%ibtypeenc)
         DEALLOCATE(adcircMesh%ibtypenc)
         DEALLOCATE(adcircMesh%nvellnc)
         DEALLOCATE(adcircMesh%nbdvnc)
         DEALLOCATE(adcircMesh%nbvvnc)
         DEALLOCATE(adcircMesh%element)
         DEALLOCATE(adcircMesh%nmnc)
         adcircMesh%initialized = .false.
      ENDIF
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE freeNetCDFCoord
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      END MODULE NETCDFIO
C-----------------------------------------------------------------------
#endif
