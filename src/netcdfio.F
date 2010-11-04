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
      ! 
      !=================================================================
      MODULE NETCDFIO
      USE SIZES, ONLY : SZ
      IMPLICIT NONE
      include 'netcdf.inc'
C
      DOUBLE PRECISION, PARAMETER ::  doubleval(1)=-99999.d0
C
      CHARACTER*35, PARAMETER ::
     &    VAR61_LN='water surface elevation above geoid'
      CHARACTER*23, PARAMETER :: VAR61_ST='water_surface_elevation'

      CHARACTER*29, PARAMETER :: VAR62_LN='station velocity'
      CHARACTER*25, PARAMETER :: VAR62_ST='station_velocity'

      CHARACTER*35, PARAMETER ::
     &   VAR63_LN='water surface elevation above geoid'
      CHARACTER*23, PARAMETER :: VAR63_ST='water_surface_elevation'

      CHARACTER*33, PARAMETER ::
     &   VAR71_LN='station air pressure at sea level'
      CHARACTER*33, PARAMETER ::
     &   VAR71_ST='station_air_pressure_at_sea_level'
      CHARACTER*15, PARAMETER :: VAR71_UN='meters of water'

      CHARACTER*25, PARAMETER :: VAR73_LN='air pressure at sea level'
      CHARACTER*25, PARAMETER :: VAR73_ST='air_pressure_at_sea_level'
      CHARACTER*15, PARAMETER :: VAR73_UN='meters of water'
      
      CHARACTER*14 , PARAMETER :: VAR_NOND='nondimensional'
      CHARACTER*1  , PARAMETER :: VAR_SECS='s'
      CHARACTER*1  , PARAMETER :: VAR_METS='m'
      CHARACTER*2  , PARAMETER :: VAR_FEET='ft'
      CHARACTER*2  , PARAMETER :: VAR_UP='up'
      CHARACTER*34 , PARAMETER ::
     &   VAR_UPY='90 degrees counterclockwise from x'
      CHARACTER*49 , PARAMETER ::
     &   VAR_UPV='90 degrees counterclockwise from x water velocity'
      CHARACTER*61 , PARAMETER :: VAR_UPW=
     &   '90 degrees counterclockwise from wind velocity in x-direction'
      CHARACTER*59 , PARAMETER :: VAR_UPWS=
     &   '90 degrees counterclockwise from wind stress in x-direction'
      CHARACTER*4  , PARAMETER :: VAR_DW='down'
      CHARACTER*5  , PARAMETER :: VAR_RG='right'
      CHARACTER*4  , PARAMETER :: VAR_EAST='east'
      CHARACTER*5  , PARAMETER :: VAR_NORTH='north'
      CHARACTER*5  , PARAMETER :: VAR_MPS='m s-1'    !'meters/second' (13)
      CHARACTER*6  , PARAMETER :: VAR_FPS='ft s-1'   !'feet/second' (11)
      CHARACTER*6  , PARAMETER :: VAR_MPS_2='m2 s-2' !'meters**2/second**2' (19)
      CHARACTER*7  , PARAMETER :: VAR_FPS_2='ft2 s-2'!'feet**2/second**2'(17)
      CHARACTER*23 , PARAMETER :: VAR_TUN='seconds since base_date'
      CHARACTER*12 , PARAMETER :: VAR_DEGE='degrees_east'
      CHARACTER*13 , PARAMETER :: VAR_DEGN='degrees_north'

      CHARACTER*19 , PARAMETER :: VAR_L_D='distance from geoid'
      CHARACTER*19 , PARAMETER :: VAR_S_D='depth'
      CHARACTER*52 , PARAMETER :: VAR_L_MIN_D=
     &   'minimum distance from geoid to water surface in grid'
      CHARACTER*52 , PARAMETER :: VAR_L_MAX_D=
     &   'maximum distance from geoid to water surface in grid'
      CHARACTER*13 , PARAMETER :: VAR_S_MIN_D='minimum_depth'
      CHARACTER*13 , PARAMETER :: VAR_S_MAX_D='maximum_depth'

      CHARACTER*22 , PARAMETER :: VAR_L_X='Cartesian coordinate x'
      CHARACTER*12 , PARAMETER :: VAR_S_X='x_coordinate'
      CHARACTER*22 , PARAMETER :: VAR_L_Y='Cartesian coordinate y'
      CHARACTER*12 , PARAMETER :: VAR_S_Y='y_coordinate'
      
      TYPE, PRIVATE :: meshStructure
         INTEGER :: X_id          ! x-coordinate or longitude
         INTEGER :: Y_id          ! y-coordinate or latitude
         INTEGER :: DEPTH_id      ! distance from geoid
         INTEGER :: ELE_id        ! elements in grid
         INTEGER :: nbdvnc_id       ! nodes on elev spec boundary seg
         INTEGER :: nbvvnc_id       ! nodes on normal flow boundary seg
         INTEGER :: nvdllnc_id      ! num nodes on elev boundary seg
         INTEGER :: nbounc_id
         INTEGER :: ibtypenc_id   ! discharge boundary type
         INTEGER :: ibtypeenc_id  ! elevation boundary type
         INTEGER :: nvellnc_id    ! nodes on norm flow spec boundary seg
         INTEGER :: nopenc_id  
         INTEGER :: slam0nc_id
         INTEGER :: sfea0nc_id
         INTEGER :: netanc_id
         INTEGER :: nvelnc_id
!        Dimension ids
         INTEGER :: num_nodes_dim_id
         INTEGER :: num_elems_dim_id
         INTEGER :: nface_dim_id
         INTEGER :: nopenc_dim_id   ! num elev spec boundary forcing segs
         INTEGER :: netanc_dim_id   ! total number of elevation specified boundary nodes
         INTEGER :: nbounc_dim_id  ! number of normal flow specified boundary segment
         INTEGER :: nvelnc_dim_id  ! total number of normal flow specified boundary nodes
         INTEGER :: ibtypenc_dim_id
         INTEGER :: ibtypeenc_dim_id  ! elevation boundary types
!        Dimension lengths
         INTEGER  num_nodes
         INTEGER  num_elems
         INTEGER  nface_len
         INTEGER  nopenc
         INTEGER  netanc
         INTEGER  nbounc
         INTEGER  nvelnc
!        Rank (number of dimensions) for each variable
         INTEGER :: depth_rank = 1
         INTEGER :: ele_rank = 2
         INTEGER :: x_rank = 1
         INTEGER :: y_rank = 1 
         INTEGER :: nbounc_rank = 1
         INTEGER :: nvdll_rank = 1
         INTEGER :: ibtypenc_rank = 1
         INTEGER :: ibtypeenc_rank = 1
         INTEGER :: nvellnc_rank = 1 
         INTEGER :: nvdllnc_rank = 1
         INTEGER :: nopenc_rank = 1
         INTEGER :: nbdvnc_rank = 2
         INTEGER :: nbvvnc_rank = 2
!        Variable shapes
         INTEGER  x_dims(1)
         INTEGER  y_dims(1)
         INTEGER  depth_dims(1)
         INTEGER  ele_dims(2)
         INTEGER  nvdll_dims(1)
         INTEGER  nbounc_dims(1)
         INTEGER  nvelnc_dims(1)
         INTEGER  ibtypenc_dims(1)
         INTEGER  ibtypeenc_dims(1)
         INTEGER  nvellnc_dims(1)
         INTEGER  nvdllnc_dims(1)
         INTEGER  nopenc_dims(1)
         INTEGER  netanc_dims(1)
         INTEGER  nbdvnc_dims(2)
         INTEGER  nbvvnc_dims(2)
         REAL(SZ), ALLOCATABLE :: xnc(:) ! x coordinate or longitude
         REAL(SZ), ALLOCATABLE :: ync(:) ! y coordinate or latitude  
         INTEGER, ALLOCATABLE ::  nbvvnc(:,:) ! boundary array
         INTEGER, ALLOCATABLE ::  nbdvnc(:,:) ! boundary array
         INTEGER, ALLOCATABLE ::  nvellnc(:)  ! boundary array
         INTEGER, ALLOCATABLE ::  nvdllnc(:)  ! boundary array
         INTEGER, ALLOCATABLE ::  nodecodenc(:)
         INTEGER, ALLOCATABLE ::  noffnc(:)     ! boundary array
         INTEGER, ALLOCATABLE ::  ibtypenc(:)     ! boundary array
         INTEGER, ALLOCATABLE ::  ibtypeenc(:)
         INTEGER, ALLOCATABLE ::  nmnc(:,:) 
         INTEGER, ALLOCATABLE ::  element(:,:)
      END TYPE meshStructure
 
      TYPE, PRIVATE :: fileData
         INTEGER record_counter
         LOGICAL createFile ! .true. if a new netCDF file must be created p
         CHARACTER*5 FILEFORT
         CHARACTER*10 FILENAME
         CHARACTER*3 FILEEXT
         LOGICAL fileFound  ! .true. if the netCDF file is present
      END Type fileData
 
      TYPE, PRIVATE :: timeData
         INTEGER :: timenc_len = 1  ! number of time slices to write
         INTEGER :: timenc_dim_id
         INTEGER :: timenc_rank = 1
         INTEGER :: timenc_id
         INTEGER :: timenc_dims(1)
         REAL(8), ALLOCATABLE :: timenc(:)
      END TYPE timeData 

      TYPE, PRIVATE :: stationData
         INTEGER ncid              ! the id of its netcdf file
         INTEGER num_stations      ! total number of stations 
         INTEGER :: num_sta_dim_id
         INTEGER, ALLOCATABLE :: name_lengths(:) ! lengths of station names 
         CHARACTER(50), POINTER :: statnames(:)
         INTEGER :: slen_dim_id
         INTEGER :: scode_dim_id
         INTEGER :: station_rank = 2
         INTEGER :: station_id
         INTEGER :: scode_rank = 2
         INTEGER :: station_dims(2)
         INTEGER :: scode_dims(2)
         INTEGER :: x_id          ! station x-coordinate or longitude
         INTEGER :: y_id          ! station y-coordinate or latitude
         INTEGER :: x_dims(1)
         INTEGER :: y_dims(1)
         INTEGER :: x_rank = 1
         INTEGER :: y_rank = 1 
         REAL(SZ), ALLOCATABLE :: station_data(:,:) ! for scalar data 
         REAL(SZ), ALLOCATABLE :: u_station_data(:,:) ! for vector component 1
         REAL(SZ), ALLOCATABLE :: v_station_data(:,:) ! for vector component 2 
         REAL(SZ), ALLOCATABLE :: x(:) ! x coordinate or longitude
         REAL(SZ), ALLOCATABLE :: y(:) ! y coordinate or latitude  
         INTEGER :: station_data_id
         INTEGER :: station_data_rank = 2
         INTEGER :: u_station_data_id
         INTEGER :: v_station_data_id
         INTEGER :: station_data_dims(2)
         type(timeData) :: myTime
         type(fileData) :: myFile
      END TYPE stationData

      type(stationData), private, save :: elevSta ! elev stations (fort.61)
      type(stationData), private, save :: prSta   ! pressure stations (fort.71)
      type(stationData), private, save :: velSta  ! velocity stations (fort.62)
      type(stationData), private, save :: wVelSta ! wind vel stations (fort.72)

      TYPE, PRIVATE :: nodalData
         INTEGER ncid              ! the id of its netcdf file
         REAL(SZ) :: initial_value ! array will be initialized to this
         INTEGER :: int_initial_value ! int array will be initialized to this
         REAL(SZ), ALLOCATABLE :: nodal_data(:,:)   ! for scalar data 
         INTEGER, ALLOCATABLE :: nodal_int(:,:)    ! for scalar integer data 
         REAL(SZ), ALLOCATABLE :: u_nodal_data(:,:) ! for 1st vector component
         REAL(SZ), ALLOCATABLE :: v_nodal_data(:,:) ! for 2nd vector component
         INTEGER :: nodal_data_id
         INTEGER :: nodal_data_rank = 2
         INTEGER :: u_nodal_data_id
         INTEGER :: v_nodal_data_id
         INTEGER :: nodal_data_dims(2)
         type(timeData) :: myTime
         type(fileData) :: myFile
         type(meshStructure) :: myMesh
      END TYPE nodalData

      type(nodalData), private, save :: elev       ! for fort.63
      type(nodalData), private, save :: pr         ! for fort.73
      type(nodalData), private, save :: currentVel ! for fort.64
      type(nodalData), private, save :: windVel    ! for fort.74

      TYPE, PRIVATE :: hotstartData
         INTEGER lun               ! ADCIRC's logical unit number
         INTEGER ncid              ! the id of its netcdf file
         INTEGER varid
         INTEGER :: FileFmtMajorFile
         INTEGER :: FileFmtMinorFile
         INTEGER :: FileFmtRevFile
         type(nodalData) :: zeta1
         type(nodalData) :: zeta2
         type(nodalData) :: zetad
         type(nodalData) :: vel
         type(nodalData) :: ch1
         type(nodalData) :: nodecodenc
         type(nodalData) :: noffnc      
         type(fileData) :: myFile
         type(timeData) :: myTime
         type(meshStructure) :: myMesh
      END TYPE hotstartData   
C
      type(hotstartData), private, save, target :: hs67 ! for fort.67
      type(hotstartData), private, save, target :: hs68 ! for fort.68
      type(hotstartData), pointer :: hs   ! indicates the current hs file    
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      CONTAINS
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   A L L O C   N E T C D F   A R R A Y S  
C-----------------------------------------------------------------------
C     jgf49.17.02 Allocates memory for NetCDF operations.
C-----------------------------------------------------------------------
      SUBROUTINE allocNetCDFArrays(lun, descript, error)
      USE GLOBAL, ONLY : OutputDataDescript_t
#ifdef CMPI
      USE MESSENGER, ONLY : msg_fini
#endif 
      IMPLICIT NONE
      INTEGER, intent(in) :: lun
      TYPE(OutputDataDescript_t) :: descript
      LOGICAL error
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering allocNetCDFArrays."
#endif
C 
      ! don't allocate or initialize anything if this output type is not
      ! using netcdf
      IF (ABS(descript % specifier).ne.3) THEN
#ifdef NETCDF_DEBUG
         write(*,*) "DEBUG: Returning from allocNetCDFArrays."
#endif
         RETURN  
      ENDIF
C
      error = .false.
      SELECT CASE(lun)
    
         CASE(61)
            CALL initStation(elevSta, lun, descript, error)
         CASE(62) 
            CALL initStation(velSta, lun, descript, error)
         CASE(63)
            CALL initNodalData(elev, lun, descript, error)
         CASE(64) 
            CALL initNodalData(currentVel, lun, descript, error)
         CASE(71) 
            CALL initStation(prSta, lun, descript, error)
         CASE(72) 
            CALL initStation(wVelSta, lun, descript, error)
         CASE(73)
            CALL initNodalData(pr, lun, descript, error)
         CASE(74) 
            CALL initNodalData(windVel, lun, descript, error)
         CASE DEFAULT
             WRITE(*,*) 'ERROR: No netCDF for this filetype.' 

      END SELECT    
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from allocNetCDFArrays."
#endif
C-----------------------------------------------------------------------
      END SUBROUTINE allocNetCDFArrays
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E    I N I T  S T A T I O N   
C-----------------------------------------------------------------------
C     jgf49.17.02 Sets up netCDF variables and allocates memory for a 
C     station.
C-----------------------------------------------------------------------
      SUBROUTINE initStation(sta, lun, descript, error)
      USE GLOBAL, ONLY : SNAMLEN, OutputDataDescript_t, ICS, RAD2DEG,
     &                   IHOT, NWS, base_date, STATNAME, STATNAMEV,
     &                   STATNAMEM
#ifdef CMPI
      USE MESSENGER, ONLY : msg_fini
#endif
      IMPLICIT NONE
C
      type(stationData), intent(inout) :: sta
      INTEGER, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: descript
C
      INTEGER i,j
      INTEGER counti(1), starti(1), n
      INTEGER kount(2), start(2)
      INTEGER iret  ! success or failure of the netcdf call
      LOGICAL error
C
!     date_string variables for time attribute
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer values(8)
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering initStation."
#endif
      error = .false.
      sta % num_stations = descript % num_fd_records
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: There are ",sta%num_stations," stations."
#endif
      ALLOCATE(sta%myTime%timenc(sta%myTime%timenc_len))
C      
C     Memory allocation for the station type, according to the number
C     of stations and time
      select case(descript%num_items_per_record)

         case(1)
            ALLOCATE(sta%station_data(sta%num_stations,
     &              sta%myTime%timenc_len))
            DO n=1, sta%num_stations
               sta%station_data(n,sta%myTime%timenc_len) 
     &             = descript%initial_value
            END DO
         case(2)
            ALLOCATE(sta%u_station_data(sta%num_stations,
     &                                  sta%myTime%timenc_len))
            ALLOCATE(sta%v_station_data(sta%num_stations,
     &                                  sta%myTime%timenc_len))
            DO n=1, sta % num_stations
               sta%u_station_data(n,sta%myTime%timenc_len) 
     &              = descript%initial_value
               sta%v_station_data(n,sta%myTime%timenc_len) 
     &              = descript%initial_value
            END DO
         case default
            WRITE(6,*) "ERROR: The station has ",
     &        descript%num_items_per_record," components; this is not ",
     &        "yet supported in ADCIRC's NetCDF."       

      end select
C
C     Initialize netCDF output file, creating a new one if necessary.
      CALL createNetCDFOutputFile(sta%ncid, sta%myFile, lun, sta%myTime,
     &                            descript, error)
C     
C     if we didn't need to create a file, update metadata and return
      IF (sta%myFile%createFile.eqv..false.) THEN
         IF (error.eqv..false.) THEN
            CALL updateMetaData(sta%ncid,sta%myFile)
         ENDIF
#ifdef NETCDF_DEBUG
         write(*,*) "DEBUG: Returning from initStation."
#endif
         RETURN
      ENDIF     
C
C     Set coordinates of each station, converting to degrees if we 
C     are in spherical coordinates. 
      ALLOCATE(sta%x(sta%num_stations))
      ALLOCATE(sta%y(sta%num_stations))
      IF(ICS.EQ.2) THEN
         sta%x = descript%x_coord * RAD2DEG
         sta%y = descript%y_coord * RAD2DEG
      ELSE
         sta%x = descript%x_coord
         sta%y = descript%y_coord
      ENDIF
C 
C     Create station dimension and station name dimension
      iret = nf_def_dim(sta%ncid, 'station',
     &       sta%num_stations, sta%num_sta_dim_id)
      CALL check_err(iret)
      iret = nf_def_dim(sta%ncid, 'namelen', SNAMLEN, sta%slen_dim_id)
      CALL check_err(iret)
!     
!     Define stations name
      sta%station_dims(1) = sta%slen_dim_id
      sta%station_dims(2) = sta%num_sta_dim_id
      iret = nf_def_var(sta%ncid, 'station_name', NF_CHAR, sta%station_rank, 
     &   sta%station_dims, sta%station_id)
      CALL check_err(iret)
C
C     Define station locations
      sta%x_dims(1) = sta%num_sta_dim_id
      iret = nf_def_var(sta%ncid, 'x', NF_DOUBLE, sta%x_rank,
     &       sta%x_dims, sta%x_id)
      sta%y_dims(1) = sta%num_sta_dim_id
      iret = nf_def_var(sta%ncid, 'y', NF_DOUBLE, sta%y_rank,
     &       sta%y_dims, sta%y_id)
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
      SELECT CASE(lun)

         CASE(61) !       F O R T . 6 1
   
            sta%station_data_dims(1) = sta%num_sta_dim_id
            sta%station_data_dims(2) = sta%myTime%timenc_dim_id
            iret = nf_def_var(sta%ncid, 'zeta', NF_DOUBLE,
     &             sta%station_data_rank, sta%station_data_dims,
     &             sta%station_data_id)
            CALL check_err(iret)
!           Define water surface elevation attributes
            iret = nf_put_att_text(sta%ncid, sta%station_data_id,
     &             'long_name',35, VAR61_LN)
            CALL check_err(iret)
            iret = nf_put_att_text(sta%ncid, sta%station_data_id,
     &             'standard_name', 23, VAR61_ST)
            CALL check_err(iret)
            CALL putUnitsAttribute(sta%ncid, sta%station_data_id,
     &           VAR_METS, VAR_FEET)
            iret = nf_put_att_DOUBLE(sta%ncid, sta%station_data_id,
     &              '_FillValue', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            iret = nf_put_att_text(sta%ncid, sta%station_data_id,
     &             'positive', 2, VAR_UP)
            CALL check_err(iret)
            sta%statnames => STATNAME 
 
         CASE(62) !       F O R T . 6 2 
 
            sta%station_data_dims(1) = sta%num_sta_dim_id
            sta%station_data_dims(2) = sta%myTime%timenc_dim_id
            iret = nf_def_var(sta%ncid, 'u-vel', NF_DOUBLE, 
     &             sta%station_data_rank, sta%station_data_dims,
     &             sta%u_station_data_id)
            CALL check_err(iret)
            iret = nf_def_var(sta%ncid, 'v-vel', NF_DOUBLE,
     &             sta%station_data_rank, sta%station_data_dims,
     &             sta%v_station_data_id)
            CALL check_err(iret)
            if (ics.eq.2) then
               ! u
               iret = nf_put_att_text(sta%ncid, sta%u_station_data_id,
     &         'long_name', 59, 
     &    'station water column vertically averaged east/west velocity')
               CALL check_err(iret)

               iret = nf_put_att_text(sta%ncid, sta%u_station_data_id,
     &                   'standard_name', 46, 
     &       'station_eastward_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf_put_att_text(sta%ncid, sta%u_station_data_id,
     &                'positive', 4, VAR_EAST)
               CALL check_err(iret)
               ! v
               iret = nf_put_att_text(sta%ncid,sta%v_station_data_id,
     &         'long_name', 61, 
     &'station water column vertically averaged north/south velocity')
               CALL check_err(iret)
               iret = nf_put_att_text(sta%ncid, sta%v_station_data_id,
     &         'standard_name', 47, 
     &         'station_northward_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf_put_att_text(sta%ncid, sta%v_station_data_id,
     &         'positive', 5, VAR_NORTH)
               CALL check_err(iret)
            else
               ! u
               iret = nf_put_att_text(sta%ncid, sta%u_station_data_id,
     &               'long_name', 64, 
     &            'station water column vertically averaged velocity in 
     &x-direction')
               CALL check_err(iret)

               iret = nf_put_att_text(sta%ncid, sta%u_station_data_id,
     &                   'standard_name', 39, 
     &       'station_x_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf_put_att_text(sta%ncid, sta%u_station_data_id,
     &                'positive', 5, VAR_RG)
               CALL check_err(iret)
               ! v
               iret = nf_put_att_text(sta%ncid, sta%v_station_data_id,
     &                'long_name', 64, 
     &            'station water column vertically averaged velocity in 
     &y-direction')
               CALL check_err(iret)
               iret = nf_put_att_text(sta%ncid, sta%v_station_data_id,
     &         'standard_name', 39, 
     &         'station_y_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf_put_att_text(sta%ncid, sta%v_station_data_id,
     &         'positive', 49, VAR_UPV)
               CALL check_err(iret)
            endif

            CALL putUnitsAttribute(sta%ncid, sta%u_station_data_id, 
     &      VAR_MPS, VAR_FPS)
            iret = nf_put_att_double(sta%ncid, sta%u_station_data_id,
     &             '_FillValue', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            iret = nf_put_att_double(sta%ncid, sta%u_station_data_id,
     &             'dry_Value', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            CALL putUnitsAttribute(sta%ncid, sta%v_station_data_id,
     &      VAR_MPS,VAR_FPS)
            iret = nf_put_att_DOUBLE(sta%ncid, sta%v_station_data_id,
     &      '_FillValue', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            iret = nf_put_att_DOUBLE(sta%ncid, sta%v_station_data_id,
     &      'dry_Value', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            sta%statnames => STATNAMEV

         CASE(71) !         F  O R T . 7 1 

            sta%station_data_dims(1) = sta%num_sta_dim_id
            sta%station_data_dims(2) = sta%myTime%timenc_dim_id
      
            iret = nf_def_var(sta%ncid, 'pressure', NF_DOUBLE,
     &      sta%station_data_rank, sta%station_data_dims, sta%station_data_id)
            CALL check_err(iret)
            iret = nf_put_att_text(sta%ncid, sta%station_data_id,
     &      'long_name',33, VAR71_LN)
            CALL check_err(iret)
            iret = nf_put_att_text(sta%ncid, sta%station_data_id,
     &      'standard_name', 25, VAR71_ST)
            CALL check_err(iret)
            iret = nf_put_att_text(sta%ncid, sta%station_data_id,
     &      'units', 15, VAR71_UN)
            CALL check_err(iret)
            iret = nf_put_att_DOUBLE(sta%ncid, sta%station_data_id,
     &      '_FillValue', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            iret = nf_put_att_text(sta%ncid, sta%station_data_id,
     &      'positive', 2, VAR_UP)
            CALL check_err(iret)
            sta%statnames => STATNAMEM 
          
         CASE(72)  !      F O R T . 7 2 
      
            sta%station_data_dims(1) = sta%num_sta_dim_id
            sta%station_data_dims(2) = sta%myTime%timenc_dim_id
            iret = nf_def_var(sta%ncid, 'windx', NF_DOUBLE,
     &      sta%station_data_rank, sta%station_data_dims,
     &      sta%u_station_data_id)
            CALL check_err(iret)
            iret = nf_def_var(sta%ncid, 'windy', NF_DOUBLE,
     &      sta%station_data_rank, sta%station_data_dims, 
     &      sta%v_station_data_id)
            CALL check_err(iret)
C
            if (ics.eq.2) then
               select case(abs(nws)) 
                  case(1,2)
                     iret = nf_put_att_text(sta%ncid, sta%u_station_data_id,
     &                      'long_name', 23, 'station e/w wind stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%u_station_data_id,
     &                      'standard_name', 36,
     &                      'station_eastward_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%u_station_data_id,
     &                      'positive', 4, VAR_EAST)
                     iret = nf_put_att_text(sta%ncid, sta%v_station_data_id,
     &                      'long_name', 23, 'station n/s wind stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%v_station_data_id,
     &                     'standard_name', 37, 
     &                     'station_northward_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%v_station_data_id,
     &                      'positive', 5, VAR_NORTH)
                     CALL check_err(iret)
                  case default       
                     iret = nf_put_att_text(sta%ncid, sta%u_station_data_id,
     &                      'long_name', 25,'station e/w wind velocity')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%u_station_data_id,
     &                      'standard_name', 21, 
     &                      'station_eastward_wind')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%u_station_data_id,
     &                      'positive', 4, VAR_EAST)
                     iret = nf_put_att_text(sta%ncid, sta%v_station_data_id,
     &                      'long_name', 25, 
     &                      'station n/s wind velocity')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%v_station_data_id,
     &                      'standard_name',22,'station_northward_wind')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%v_station_data_id,
     &                      'positive', 5,VAR_NORTH)
               end select
            else
               select case(abs(nws))
                  case(1,2)
                     iret = nf_put_att_text(sta%ncid, sta%u_station_data_id,
     &                      'long_name', 33, 
     &                      'station wind stress in x-direction')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%u_station_data_id,
     &                      'standard_name', 29, 
     &                      'station_x_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%u_station_data_id,
     &                      'positive', 5, VAR_RG)
                     iret = nf_put_att_text(sta%ncid, sta%v_station_data_id,
     &                      'long_name', 33, 
     &                      'station wind stress in y-direction')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%v_station_data_id,
     &                      'standard_name', 29, 
     &                      'station_y_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%v_station_data_id,
     &                     'positive', 59, VAR_UPWS)
                  case default 
                     iret = nf_put_att_text(sta%ncid, sta%u_station_data_id,
     &                       'long_name', 36, 
     &                       'station wind velocity in x-direction')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%u_station_data_id,
     &                      'standard_name', 36, 'station_x_wind')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%u_station_data_id,
     &                      'positive', 5, VAR_RG)
                     iret = nf_put_att_text(sta%ncid, sta%v_station_data_id,
     &                      'long_name', 36, 
     &                      'station wind velocity in y-direction')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%v_station_data_id,
     &                      'standard_name', 36, 
     &                      'station_y_wind')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%v_station_data_id,
     &                      'positive', 61, VAR_UPW)
                     CALL check_err(iret)
               end select
            endif

            if (abs(nws).gt.2.and.abs(nws).lt.100) then
               call putUnitsAttribute(sta%ncid, sta%u_station_data_id,
     &              VAR_MPS, VAR_FPS)
               call putUnitsAttribute(sta%ncid, sta%v_station_data_id,
     &              VAR_MPS, VAR_FPS)
            else
               call putUnitsAttribute(sta%ncid, sta%u_station_data_id,
     &              VAR_MPS_2, VAR_FPS_2)
               call putUnitsAttribute(sta%ncid, sta%u_station_data_id,
     &              VAR_MPS_2, VAR_FPS_2)
            endif

            iret = nf_put_att_double(sta%ncid, sta%u_station_data_id,
     &             '_FillValue', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            iret = nf_put_att_DOUBLE(sta%ncid, sta%v_station_data_id,
     &             '_FillValue', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            sta%statnames => STATNAMEM !jgfdebug
    
         CASE DEFAULT
            WRITE(*,*) "ERROR: The unit number ",lun," cannot be ",
     &                "initialized as a NetCDF file."
                 
          
      END SELECT
!  
!     Define station names and code attributes
      iret = nf_put_att_text(sta%ncid,sta%station_id,'long_name',12,
     &       'station name')
      CALL check_err(iret)
!
!     Define time attributes
      iret = nf_put_att_text(sta%ncid, sta%myTime%timenc_id,
     &       'long_name',10,'model time')
      CALL check_err(iret)
      iret = nf_put_att_text(sta%ncid,sta%myTime%timenc_id,
     &       'standard_name', 4, 'time')
      CALL check_err(iret)
      iret = nf_put_att_text(sta%ncid, sta%myTime%timenc_id, 'units', 23, 
     .       'seconds since base_date')
      CALL check_err(iret)
      iret = nf_put_att_text(sta%ncid, sta%myTime%timenc_id,
     &       'base_date', metalength(adjustl(base_date)),
     &        adjustl(base_date))
      CALL check_err(iret)

      CALL date_and_time(now_date,big_ben,zone,values)
      WRITE(date_string,71) values(1),values(2),values(3),
     &      values(5),values(6),values(7),   (values(4))/60
 71   FORMAT(I4,'-',I2.2,'-',i2.2,' ',i2,':',i2.2,':',i2.2,' '
     &  ,i3.2,':00')
      iret = nf_put_att_text(sta%ncid, NF_GLOBAL,'creation_date',
     &       len_trim(date_string),date_string)
      CALL check_err(iret)
      iret = nf_put_att_text(sta%ncid, NF_GLOBAL,'modification_date',
     &       len_trim(date_string),date_string)
      CALL check_err(iret)
C
C     define metadata and selected fort.15 parameters in netcdf file
      CALL defineMetaData(sta%ncid)
!
!     Leave define mode
      iret = nf_enddef(sta%ncid)
      CALL check_err(iret)
!   
!     Store station name
      do i=1,sta%num_stations
         start(1)=1
         start(2)=i
         kount(1)= LEN(TRIM(sta%statnames(i)))
         kount(2)=1
         iret = nf_put_vara_text(sta%ncid,sta%station_id,start,kount,
     &   sta%statnames(i) )
         CALL check_err(iret)
      end do
C
C     Store station locations
      iret = nf_put_var_double(sta%ncid, sta%x_id, sta%x)
      CALL check_err(iret)
      iret = nf_put_var_double(sta%ncid, sta%y_id, sta%y)
      CALL check_err(iret)
C
C     now close the initialized netcdf file
      iret = nf_close(sta%ncid)
      CALL check_err(iret)
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from initStation."
#endif
C-----------------------------------------------------------------------
      END SUBROUTINE initStation
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E    I N I T  N O D A L  D A T A  
C-----------------------------------------------------------------------
C     jgf49.17.02 Sets up netCDF variables and allocates memory for  
C     full domain data.
C-----------------------------------------------------------------------
      SUBROUTINE initNodalData(dat, lun, descript, error)
      USE GLOBAL, ONLY : OutputDataDescript_t, ICS, NWS, NE_G,base_date
      IMPLICIT NONE
C
      type(nodalData), intent(inout) :: dat
      INTEGER, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: descript
      LOGICAL :: error
C
      INTEGER n    ! node counter
      INTEGER iret ! success or failure of netcdf call
C 
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering initNodalData."
#endif
      error = .false.
      dat%myMesh%num_nodes = descript % num_fd_records
      dat%myMesh%num_elems = NE_G
      dat%myMesh%nface_len = 3
      ALLOCATE(dat%myTime%timenc(dat%myTime%timenc_len))
     
      SELECT CASE(descript%num_items_per_record)
      CASE(1)
         ALLOCATE(dat%nodal_data(dat%myMesh%num_nodes,
     &       dat%myTime%timenc_len))
C        initialization
         DO n=1, dat%myMesh%num_nodes
            dat%nodal_data(n,dat%myTime%timenc_len)
     &         = descript%initial_value
         END DO
      CASE(2)
         ALLOCATE(dat%u_nodal_data(dat%myMesh%num_nodes,
     &       dat%myTime%timenc_len))
         ALLOCATE(dat%v_nodal_data(dat%myMesh%num_nodes,
     &       dat%myTime%timenc_len))
C        initialization
         DO n=1, dat%myMesh%num_nodes
            dat%u_nodal_data(n,dat%myTime%timenc_len) 
     &         = descript%initial_value
            dat%v_nodal_data(n,dat%myTime%timenc_len) 
     &         = descript%initial_value
         END DO
      CASE DEFAULT
         WRITE(*,*) "ADCIRC's NetCDF cannot write ",
     &     descript%num_items_per_record, " dimensional arrays."
      END SELECT
C
C     Initialize netCDF output file, creating a new one if necessary.
      CALL createNetCDFOutputFile(dat%ncid, dat%myFile, lun, dat%myTime,
     &                            descript, error)
C     
C     if we didn't need to create a file, update metadata and return
      IF (dat%myFile%createFile.eqv..false.) THEN
         IF (error.eqv..false.) THEN
            CALL updateMetaData(dat%ncid,dat%myFile)
         ENDIF
#ifdef NETCDF_DEBUG
         write(*,*) "DEBUG: Returning from initNodalData."
#endif
         RETURN
      ENDIF     
C
      CALL initNetCDFCoord(dat%myMesh)
      CALL defineMeshVariables(dat%ncid, dat%myMesh)
c
      SELECT CASE(lun)
C
         CASE(63)
C
            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            iret = nf_def_var(dat%ncid, 'zeta', NF_DOUBLE,
     &             dat%nodal_data_rank, dat%nodal_data_dims,
     &             dat%nodal_data_id)
            CALL check_err(iret)     
            iret = nf_put_att_text(dat%ncid, dat%nodal_data_id,
     &            'long_name',35, VAR63_LN)
            CALL check_err(iret)
            iret = nf_put_att_text(dat%ncid, dat%nodal_data_id,
     &             'standard_name', 23, VAR63_ST)
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%nodal_data_id,
     &           VAR_METS, VAR_FEET)  
            iret = nf_put_att_DOUBLE(dat%ncid, dat%nodal_data_id,
     &             '_FillValue', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            iret = nf_put_att_text(dat%ncid, dat%nodal_data_id,
     &             'positive', 2, VAR_UP)
            CALL check_err(iret)

         CASE(64)

            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            iret = nf_def_var(dat%ncid, 'u-vel', NF_DOUBLE,
     &             dat%nodal_data_rank, dat%nodal_data_dims,
     &             dat%u_nodal_data_id)
            CALL check_err(iret)
            iret = nf_def_var(dat%ncid, 'v-vel', NF_DOUBLE,
     &             dat%nodal_data_rank, dat%nodal_data_dims,
     &             dat%v_nodal_data_id)
            CALL check_err(iret)
            if (ics.eq.2) then
               iret = nf_put_att_text(dat%ncid, dat%u_nodal_data_id,
     &                'long_name', 51, 
     &   'water column vertically averaged east/west velocity')
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%u_nodal_data_id,
     &                'standard_name', 38, 
     &                'eastward_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%u_nodal_data_id,
     &                'positive', 4, VAR_EAST)
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%v_nodal_data_id,
     &                'long_name', 53, 
     &    'water column vertically averaged north/south velocity')
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%v_nodal_data_id,
     &               'standard_name', 39, 
     &                      'northward_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%v_nodal_data_id,
     &                'positive', 5, VAR_NORTH)
               CALL check_err(iret)
            else
               iret = nf_put_att_text(dat%ncid, dat%u_nodal_data_id,
     &                'long_name', 56, 
     &  'water column vertically averaged velocity in x-direction')
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%u_nodal_data_id, 
     &               'standard_name', 31, 
     &               'x_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%u_nodal_data_id,
     &                'positive', 5, VAR_RG)
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%v_nodal_data_id,
     &                'long_name', 56, 
     &       'water column vertically averaged velocity in y-direction')
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%v_nodal_data_id, 
     &               'standard_name', 31, 
     &               'y_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%v_nodal_data_id,
     &                'positive', 49, VAR_UPV)
               CALL check_err(iret)
            endif
            CALL putUnitsAttribute(dat%ncid, dat%u_nodal_data_id,
     &           VAR_MPS, VAR_FPS)
            CALL putUnitsAttribute(dat%ncid, dat%v_nodal_data_id,
     &           VAR_MPS, VAR_FPS)
            iret = nf_put_att_double(dat%ncid, dat%u_nodal_data_id,
     &           '_FillValue', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
      
            iret = nf_put_att_double(dat%ncid, dat%u_nodal_data_id,
     &           'dry_Value', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            iret = nf_put_att_double(dat%ncid, dat%v_nodal_data_id,
     &           '_FillValue', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
      
            iret = nf_put_att_double(dat%ncid, dat%v_nodal_data_id,
     &            'dry_Value', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)

         CASE(73)

            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id     
            iret = nf_def_var(dat%ncid,'pressure',NF_DOUBLE,
     &              dat%nodal_data_rank, dat%nodal_data_dims,
     &              dat%nodal_data_id) 
            CALL check_err(iret)
!           Define sea surface pressure attributes
            iret = nf_put_att_text(dat%ncid, dat%nodal_data_id,
     &             'long_name',25, VAR73_LN)
            CALL check_err(iret)
            iret = nf_put_att_text(dat%ncid, dat%nodal_data_id,
     &             'standard_name', 25, VAR73_ST)
            CALL check_err(iret)
            iret = nf_put_att_text(dat%ncid, dat%nodal_data_id, 'units',
     &             15, VAR73_UN)
            CALL check_err(iret)
            iret = nf_put_att_DOUBLE(dat%ncid, dat%nodal_data_id,
     &            '_FillValue', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            iret = nf_put_att_text(dat%ncid, dat%nodal_data_id, 
     &             'positive', 2, VAR_UP)
            CALL check_err(iret)

         CASE(74)

            dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
            dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
            iret = nf_def_var(dat%ncid, 'windx', NF_DOUBLE,
     &             dat%nodal_data_rank, dat%nodal_data_dims,
     &             dat%u_nodal_data_id)
            CALL check_err(iret)
            iret = nf_def_var(dat%ncid, 'windy', NF_DOUBLE,
     &             dat%nodal_data_rank, dat%nodal_data_dims,
     &             dat%v_nodal_data_id)
            CALL check_err(iret)
            if (ics.eq.2) then
               select case(abs(nws))
                  case(1,2)
                     iret = nf_put_att_text(dat%ncid, dat%u_nodal_data_id,
     &                      'long_name', 15, 
     &                       'e/w wind stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%u_nodal_data_id,
     &                      'standard_name', 28, 
     &                      'eastward_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%u_nodal_data_id,
     &                      'positive', 4, VAR_EAST)
                     iret = nf_put_att_text(dat%ncid, dat%v_nodal_data_id,
     &                       'long_name', 15, 
     &                       'n/s wind stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%v_nodal_data_id,
     &                      'standard_name', 29, 
     &                      'northward_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%v_nodal_data_id,
     &                      'positive', 5, VAR_NORTH)
                     CALL check_err(iret)

                  case default
                     iret = nf_put_att_text(dat%ncid, dat%u_nodal_data_id,
     &               'long_name', 17, 
     &              'e/w wind velocity')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%u_nodal_data_id,
     &                      'standard_name', 13, 
     &                      'eastward_wind')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%u_nodal_data_id,
     &                     'positive', 4, VAR_EAST)
                     iret = nf_put_att_text(dat%ncid, dat%v_nodal_data_id,
     &                     'long_name', 17, 
     &       'n/s wind velocity')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%v_nodal_data_id,
     &                       'standard_name', 14, 
     &       'northward_wind')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%v_nodal_data_id,
     &                      'positive', 5, VAR_NORTH)
               end select
            else
               select case(abs(nws))
                  case(1,2)
                     iret = nf_put_att_text(dat%ncid, dat%u_nodal_data_id,
     &                       'long_name', 26, 
     &                       'wind stress in x-direction')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%u_nodal_data_id,
     &                      'standard_name', 21, 
     &                      'x_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%u_nodal_data_id,
     &                      'positive', 5, VAR_RG)
                     CALL check_err(iret)                   
                     iret = nf_put_att_text(dat%ncid, dat%v_nodal_data_id,
     &                       'long_name', 26, 
     &                       'wind stress in y-direction')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%v_nodal_data_id,
     &                       'standard_name', 21, 
     &                       'y_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%v_nodal_data_id,
     &                     'positive', 59, VAR_UPWS)
                  case default
                     iret = nf_put_att_text(dat%ncid, dat%u_nodal_data_id,
     &                      'long_name', 28,
     &                      'wind velocity in x-direction')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%u_nodal_data_id,
     &                      'standard_name', 28, 
     &                      'x_wind')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%u_nodal_data_id,
     &                     'positive', 5, VAR_RG)

                     iret = nf_put_att_text(dat%ncid, dat%v_nodal_data_id,
     &                       'long_name', 28, 
     &                       'wind velocity in y-direction')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%v_nodal_data_id,
     &                      'standard_name', 28, 
     &                      'y_wind')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%v_nodal_data_id,
     &                      'positive', 61, VAR_UPW)
                     CALL check_err(iret)
               end select
            endif

            select case(abs(nws))
               case(1,2)
                  call putUnitsAttribute(dat%ncid, dat%u_nodal_data_id,
     &                 VAR_MPS_2, VAR_FPS_2)
                  call putUnitsAttribute(dat%ncid, dat%v_nodal_data_id,
     &                 VAR_MPS_2, VAR_FPS_2)
               case default
                  call putUnitsAttribute(dat%ncid, dat%u_nodal_data_id,
     &                 VAR_MPS, VAR_FPS)
                  call putUnitsAttribute(dat%ncid, dat%v_nodal_data_id,
     &                 VAR_MPS, VAR_FPS)
            end select
            iret = nf_put_att_double(dat%ncid, dat%u_nodal_data_id,
     &             '_FillValue', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            iret = nf_put_att_DOUBLE(dat%ncid, dat%v_nodal_data_id,
     &             '_FillValue', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
         CASE DEFAULT

      END SELECT
C
C     RJW added 9/13/2010
C     to include time atributes in global data files
C     Define time attributes
      iret = nf_put_att_text(dat%ncid, dat%myTime%timenc_id,
     &       'long_name',10,'model time')
      CALL check_err(iret)
      iret = nf_put_att_text(dat%ncid,dat%myTime%timenc_id,
     &       'standard_name', 4, 'time')
      CALL check_err(iret)
      iret = nf_put_att_text(dat%ncid, dat%myTime%timenc_id, 'units', 
     &       23, 'seconds since base_date')
      CALL check_err(iret)
      iret = nf_put_att_text(dat%ncid, dat%myTime%timenc_id,
     &       'base_date', metalength(adjustl(base_date)),
     &        adjustl(base_date))
      CALL check_err(iret)
C
C     define metadata and selected fort.15 parameters in netcdf file
      CALL defineMetaData(dat%ncid)
C
C     Leave define mode
      iret = nf_enddef(dat%ncid)
      CALL check_err(iret)
C
C     write mesh to netcdf file
      CALL putMeshVariables(dat%ncid,dat%myMesh)
C
C     now close the initialized netcdf file
      iret = nf_close(dat%ncid)
      CALL check_err(iret)

#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from initNodalData."
#endif
C-----------------------------------------------------------------------
      END SUBROUTINE initNodalData
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E  S E T  N E T C D F  F I L E  N A M E  
C-----------------------------------------------------------------------
C     jgf49.41 Sets the name of the netcdf file. 
C-----------------------------------------------------------------------
      SUBROUTINE setNetCDFFileName(lun, myFile)
      IMPLICIT NONE
      INTEGER, intent(in) :: lun ! logical unit number
      TYPE(fileData) :: myFile
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering setNetCDFFileName."
#endif
      myFile%FILEFORT = 'fort.'
      myFile%FILEEXT = '.nc'
      WRITE(myFile%FILENAME,'(A5,I2,A3)') 
     &      myFile%FILEFORT, lun, myFile%FILEEXT
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from setNetCDFFileName."
#endif
C-----------------------------------------------------------------------
      END SUBROUTINE setNetCDFFileName
C-----------------------------------------------------------------------



C-----------------------------------------------------------------------
C     S U B R O U T I N E  C R E A T E  N E T C D F  O U T P U T  F I L E 
C-----------------------------------------------------------------------
C     jgf49.22 Creates a new netcdf output file if needed. 
C-----------------------------------------------------------------------
      SUBROUTINE createNetCDFOutputFile(ncid, myFile, lun, myTime,
     &                                  descript, error)
      USE GLOBAL, ONLY : OutputDataDescript_t, IHOT, screenUnit
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      IMPLICIT NONE
      INTEGER ncid
      TYPE(fileData) :: myFile
      INTEGER :: lun
      TYPE(timeData) :: myTime
      TYPE(OutputDataDescript_t) :: descript
      LOGICAL :: error
C
      INTEGER iret
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering createNetCDFOutputFile."
#endif
      error = .false.
      myFile%createFile = .false.
      myFile%fileFound = .false.
      CALL setNetCDFFileName(lun,myFile)
C
C     jgf49.17.02: Simplified the criteria for creating a new netCDF
C     output file: coldstart, overwrite upon hotstart, or output
C     that does not already exist. These criteria do not apply to 
C     netcdf hotstart files; we will always write those when called
C     upon to create them.
      INQUIRE(FILE=myFile%FILENAME,EXIST=myFile%fileFound)
      IF ((lun.ne.67).and.(lun.ne.68)) THEN
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
               write(screenUnit,*) 
     &           "ERROR: The NetCDF output file '",myFile%FILENAME,
     &           "' was not found. It should have been created by 
     &adcprep."
               write(screenUnit,*) "ERROR: ADCIRC execution terminated."
               error = .true.
            ENDIF
#else
            ! these lines are executed by serial adcirc and adcprep
            myFile%createFile = .true.
            myFile%record_counter = 1
#endif
         ENDIF
      ELSE
         ! these lines are executed to create netcdf hotstart files
         ! don't overwrite existing hotstart files
         IF (myFile%fileFound.eqv..true.) THEN
            write(*,*) "INFO: Hotstart file ",lun,
     &        " already exists."
           myFile%createFile = .false.
         ELSE
            myFile%createFile = .true.
            myFile%record_counter = 1 
         ENDIF        
      ENDIF
C     
C     RETURN if we don't need to create a file.
      IF (myFile%createFile.eqv..false.) THEN
#ifdef NETCDF_DEBUG
         write(*,*) "DEBUG: Returning from createNetCDFOutputFile."
#endif
         RETURN
      ENDIF     
C
      iret = nf_create(myFile%FILENAME, NF_CLOBBER, ncid)
      CALL check_err(iret)
C
C       
!     Define time
      iret = nf_def_dim(ncid, 'time', nf_unlimited,
     &       myTime%timenc_dim_id)
      CALL check_err(iret)
      myTime%timenc_dims(1) = myTime%timenc_dim_id
      iret = nf_def_var(ncid, 'time', NF_DOUBLE, myTime%timenc_rank,
     &       myTime%timenc_dims, myTime%timenc_id)
      CALL check_err(iret)

#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from createNetCDFOutputFile."
#endif
C-----------------------------------------------------------------------
      END SUBROUTINE createNetCDFOutputFile
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   I N I T   N E T C D F   C O O R D  
C-----------------------------------------------------------------------
C     jgf49.17.02 Initializes NetCDF coordinates. 
C-----------------------------------------------------------------------
      SUBROUTINE initNetCDFCoord(myMesh)
      USE GLOBAL, ONLY : RAD2DEG, NBOU, NVEL, NOPE, NP_G, SLAM0, SFEA0,
     &                   NBVV, NVDLL, NBDV, NVELL, NETA, X, Y, SLAM,
     &                   SFEA, ICS, IBTYPEE, IBTYPE
      IMPLICIT NONE
      type(meshStructure), intent(inout) :: myMesh
      INTEGER :: i, j, k  ! array indices
      INTEGER :: max_seg_nodes  ! size of longest list of nodes on a bdry seg
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering initNetCDFCoord."
#endif
C
      ALLOCATE(myMesh%xnc(NP_G))
      ALLOCATE(myMesh%ync(NP_G))
      ALLOCATE(myMesh%nvdllnc(nope))
      ALLOCATE(myMesh%ibtypeenc(nope))
      ALLOCATE(myMesh%ibtypenc(nbou))
      ALLOCATE(myMesh%nvellnc(nbou))
      max_seg_nodes = maxval(nvdll) ! dimension of seg with most nodes
      ALLOCATE(myMesh%nbdvnc(nope,max_seg_nodes))
      max_seg_nodes = maxval(nvell) ! dimension of seg with most nodes
      ALLOCATE(myMesh%nbvvnc(nbou,max_seg_nodes))
      ALLOCATE(myMesh%element(myMesh%nface_len,myMesh%num_elems))
      ALLOCATE(myMesh%nmnc(myMesh%num_elems,myMesh%nface_len))
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
      myMesh%nopenc = nope
      myMesh%netanc = neta
      DO k=1,nope
         myMesh%nvdllnc(k) = nvdll(k)
         myMesh%ibtypeenc(k) = ibtypee(k)
         DO j=1,nvdll(k)
            myMesh%nbdvnc(k,j) = nbdv(k,j)
         END DO  
      END DO
!
!     normal flow (discharge) specified boundary segments
      myMesh%nbounc = nbou
      myMesh%nvelnc = nvel
C
      DO k=1,nbou
         myMesh%nvellnc(k) = nvell(k)
         myMesh%ibtypenc(k) = ibtype(k)
         DO J=1,nvell(k)
            myMesh%nbvvnc(k,j)=nbvv(k,j)
         END DO
      END DO
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from initNetCDFCoord."
#endif
C-----------------------------------------------------------------------
      END SUBROUTINE initNetCDFCoord
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E   W R I T E  O U T  A R R A Y  N E T C D F 
C-----------------------------------------------------------------------
C     jgf48.03 This subroutine was created from cf's code in timestep.F to
C     write output files in NetCDF format.
C-----------------------------------------------------------------------
      SUBROUTINE writeOutArrayNetCDF(lun, time, it, descript)
      USE GLOBAL, ONLY : OutputDataDescript_t, ScreenUnit
C
      IMPLICIT NONE
C     args
      INTEGER, intent(in) :: lun ! logical unit number of file to write to
      REAL(8), intent(in) :: time ! seconds since cold start
      INTEGER, intent(in) :: it   ! number of time steps since cold start
      type(OutputDataDescript_t), intent(in) :: descript !describes output data
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering writeOutArrayNetCDF."
#endif
C
      SELECT CASE(lun)
C    
         CASE(61)
            CALL writeStation(elevSta, lun, descript, time)
         CASE(62) 
            CALL writeStation(velSta, lun, descript, time)
         CASE(63)
            CALL writeNodalData(elev, lun, descript, time)
         CASE(64) 
            CALL writeNodalData(currentVel, lun, descript, time)
         CASE(71) 
            CALL writeStation(prSta, lun, descript, time)
         CASE(72) 
            CALL writeStation(wVelSta, lun, descript, time)
         CASE(73)
            CALL writeNodalData(pr, lun, descript, time)
         CASE(74) 
            CALL writeNodalData(windVel, lun, descript, time)
         CASE DEFAULT
             WRITE(ScreenUnit,*) 'ERROR: No netCDF for this filetype.' 
      END SELECT
C

#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from writeOutArrayNetCDF."
#endif
C-----------------------------------------------------------------------
      END SUBROUTINE writeOutArrayNetCDF
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E    W R I T E  S T A T I O N   
C-----------------------------------------------------------------------
C     jgf49.17.02 Writes data to station file.
C-----------------------------------------------------------------------
      SUBROUTINE writeStation(sta, lun, descript, time)
      USE SIZES, ONLY : MNPROC
      USE GLOBAL, ONLY : SNAMLEN, OutputDataDescript_t
      IMPLICIT NONE
C
      type(stationData), intent(inout) :: sta
      INTEGER, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: descript
      REAL(8), intent(in) :: time
C
      INTEGER i,j
      INTEGER counti(1), starti(1), n
      INTEGER kount(2), start(2)
      INTEGER iret  ! success or failure of netcdf call
      INTEGER recid ! netcdf id of the record counter information
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering writeStation."
#endif
      iret = nf_open(sta%myFile%FILENAME, NF_WRITE, sta%ncid)
      CALL check_err(iret)   
!   
!     Inquire variables
      iret=nf_inq_unlimdim(sta%ncid,recid)
      CALL check_err(iret)   
      iret=nf_inq_dimlen(sta%ncid,recid,sta%myFile%record_counter)
      CALL check_err(iret)   
C
C     Increment the record counter so that we can store data at the
C     next location in the netcdf file
      sta%myFile%record_counter = sta%myFile%record_counter + 1       
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: sta%myFile%record_counter = ",
     &           sta%myFile%record_counter
#endif
!    
!     Store time
      iret=nf_inq_varid(sta%ncid, "time", sta%myTime%timenc_id)
      starti(1)=sta%myFile%record_counter
      counti(1)=sta%myTime%timenc_len
      sta%myTime%timenc(sta%myTime%timenc_len)=time
      iret = nf_put_vara_DOUBLE(sta%ncid, sta%myTime%timenc_id,
     &       starti, counti, sta%myTime%timenc)
      CALL check_err(iret)
C
      SELECT CASE(lun)

      CASE(61)
         iret=nf_inq_varid(sta%ncid, "zeta", sta%station_data_id)
         CALL check_err(iret)
      CASE(62)
         iret=nf_inq_varid(sta%ncid, "u-vel", sta%u_station_data_id)
         CALL check_err(iret)
         iret=nf_inq_varid(sta%ncid, "v-vel", sta%v_station_data_id)
         CALL check_err(iret)
      CASE(71)
         iret=nf_inq_varid(sta%ncid, "pressure", sta%station_data_id)
         CALL check_err(iret)
      CASE(72)  
           iret=nf_inq_varid(sta%ncid,"windx",sta%u_station_data_id)
           CALL check_err(iret)
           iret=nf_inq_varid(sta%ncid,"windy",sta%v_station_data_id)
           CALL check_err(iret)

      CASE DEFAULT
           WRITE(*,*) "ERROR: The unit number ",lun," is not ",
     &                "available in NetCDF."
          
      END SELECT
C
C     Set up the netcdf data extents
      kount(1)=sta%num_stations
      kount(2)=sta%myTime%timenc_len
      start(1)=1
      start(2)=sta%myFile%record_counter
C
      SELECT CASE(descript%num_items_per_record) 
      CASE(1)
         IF (MNPROC.eq.1) THEN
            sta%station_data(:,kount(2)) = descript%array 
         ELSE
            sta%station_data(:,kount(2)) = descript%array_g
         ENDIF
#ifdef NETCDF_DEBUG
         DO i=1,sta%num_stations
            write(*,*) "DEBUG: station_data(",i,")=",
     &                  sta%station_data(i,kount(2))
         ENDDO
#endif
         iret = nf_put_vara_DOUBLE(sta%ncid, sta%station_data_id,
     &          start, kount, sta%station_data)
         CALL check_err(iret)
      CASE(2)
         IF (MNPROC.eq.1) THEN
            sta%u_station_data(:,kount(2)) = descript%array
            sta%v_station_data(:,kount(2)) = descript%array2 
         ELSE
            sta%u_station_data(:,kount(2)) = descript%array_g
            sta%v_station_data(:,kount(2)) = descript%array2_g
         ENDIF
         iret = nf_put_vara_DOUBLE(sta%ncid, sta%u_station_data_id,
     &          start, kount, sta%u_station_data)
         CALL check_err(iret)
         iret = nf_put_vara_DOUBLE(sta%ncid, sta%v_station_data_id,
     &          start, kount, sta%v_station_data)
         CALL check_err(iret)
      CASE DEFAULT
         WRITE(*,*) "ADCIRC's NetCDF cannot write ",
     &         descript%num_items_per_record, " dimensional arrays."
      END SELECT
!
!     Close netCDF file
9999  iret = nf_close(sta%ncid)
      CALL check_err(iret)

#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from writeStation."
#endif
C-----------------------------------------------------------------------
      end subroutine writeStation
C-----------------------------------------------------------------------     


C-----------------------------------------------------------------------
C     S U B R O U T I N E    W R I T E  N O D A L   D A T A   
C-----------------------------------------------------------------------
C     jgf49.17.02 Writes data to a full domain file.
C-----------------------------------------------------------------------
      SUBROUTINE writeNodalData(dat, lun, descript, time)
      USE SIZES, ONLY : MNPROC
      USE GLOBAL, ONLY : SNAMLEN, OutputDataDescript_t, NODECODE
      IMPLICIT NONE
C
      type(nodalData), intent(inout) :: dat
      INTEGER, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: descript
      REAL(8), intent(in) :: time
C
      INTEGER i,j
      INTEGER counti(1), starti(1), n
      INTEGER kount(2), start(2)
      INTEGER iret ! success or failure of the netcdf call
      INTEGER recid ! netcdf id of the record counter info
      INTEGER num_values ! number of values in the array
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering writeNodalData."
#endif
C
      iret = nf_open(dat%myFile%FILENAME, NF_WRITE, dat%ncid)
      CALL check_err(iret)   
!   
!     Inquire variables
      iret=nf_inq_unlimdim(dat%ncid,recid)
      CALL check_err(iret)   
C
C
      iret=nf_inq_dimlen(dat%ncid,recid,dat%myFile%record_counter)
      CALL check_err(iret)   
C
C     Increment the record counter so that we can store data at the
C     next location in the netcdf file
      dat%myFile%record_counter = dat%myFile%record_counter + 1       
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: dat%myFile%record_counter = ",
     &            dat%myFile%record_counter
#endif
!     
!     Store time
      iret=nf_inq_varid(dat%ncid,"time",dat%myTime%timenc_id)
      starti(1)=dat%myFile%record_counter
      counti(1)=dat%myTime%timenc_len
      dat%myTime%timenc(dat%myTime%timenc_len)=time
      iret = nf_put_vara_DOUBLE(dat%ncid, dat%myTime%timenc_id,
     &       starti,counti, dat%myTime%timenc)
      CALL check_err(iret)
C
      SELECT CASE(lun)

      CASE(63)
         iret=nf_inq_varid(dat%ncid,"zeta",dat%nodal_data_id)
         CALL check_err(iret)
      CASE(64)
         iret=nf_inq_varid(dat%ncid, "u-vel", dat%u_nodal_data_id)
         CALL check_err(iret)
         iret=nf_inq_varid(dat%ncid, "v-vel", dat%v_nodal_data_id)
         CALL check_err(iret)
      CASE(73)
         iret=nf_inq_varid(dat%ncid, "pressure", dat%nodal_data_id)
         CALL check_err(iret)
      CASE(74)
         iret=nf_inq_varid(dat%ncid,"windx",dat%u_nodal_data_id)
         CALL check_err(iret)
         iret=nf_inq_varid(dat%ncid,"windy",dat%v_nodal_data_id)
         CALL check_err(iret)
      CASE DEFAULT
         WRITE(*,*) "ERROR: NetCDF output is not available for the ",
     &         "ADCIRC unit number ",lun,"."
      END SELECT
C
C     Set up the netcdf data extents
      kount(2)=dat%myTime%timenc_len
      start(1)=1
      start(2)=dat%myFile%record_counter
C     
C     Set up the array values that need to be written
      IF (MNPROC.eq.1) THEN ! SERIAL
         kount(1)=descript%num_records_this
         IF ( descript%ConsiderWetDry.eqv..TRUE.) THEN
            SELECT CASE(descript%num_items_per_record) 
            CASE(1)
               WHERE (NODECODE.EQ.0) 
                  dat%nodal_data(:,kount(2)) = descript%alternate_value
               ELSEWHERE 
                  dat%nodal_data(:,kount(2)) = descript%array
               END WHERE
            CASE(2)
               WHERE (NODECODE.EQ.0) 
                  dat%u_nodal_data(:,kount(2))=descript%alternate_value
                  dat%v_nodal_data(:,kount(2))=descript%alternate_value
               ELSEWHERE 
                  dat%u_nodal_data(:,kount(2))=descript%array
                  dat%v_nodal_data(:,kount(2))=descript%array
               END WHERE
            CASE DEFAULT
               WRITE(*,*) "ADCIRC's NetCDF cannot write ",
     &         descript%num_items_per_record, " dimensional arrays."
            END SELECT
         ELSE 
            SELECT CASE(descript%num_items_per_record) 
            CASE(1)
               dat%nodal_data(:,kount(2)) = descript%array 
            CASE(2)
               dat%u_nodal_data(:,kount(2)) = descript%array
               dat%v_nodal_data(:,kount(2)) = descript%array2 
            CASE DEFAULT
               WRITE(*,*) "ERROR: ADCIRC's NetCDF cannot write ",
     &         descript%num_items_per_record, " dimensional arrays."
            END SELECT
         ENDIF
      ELSE                  ! PARALLEL
         ! the globalio routine deals with inserting alternate values for 
         ! dry nodes, so we don't need to worry about it here (ie, in parallel)
         kount(1)=descript%num_fd_records
         SELECT CASE(descript%num_items_per_record) 
         CASE(1)
            dat%nodal_data(:,kount(2)) = descript%array_g
         CASE(2)
            dat%u_nodal_data(:,kount(2)) = descript%array_g
            dat%v_nodal_data(:,kount(2)) = descript%array2_g
         CASE DEFAULT
            WRITE(*,*) "ERROR: ADCIRC's NetCDF cannot write ",
     &      descript%num_items_per_record, " dimensional arrays."
         END SELECT
      ENDIF
C            
C     Write the values
      SELECT CASE(descript%num_items_per_record)
      CASE(1)
         iret = nf_put_vara_DOUBLE(dat%ncid, dat%nodal_data_id, start,
     &       kount, dat%nodal_data)
         CALL check_err(iret)
      CASE(2)
         iret=nf_put_vara_DOUBLE(dat%ncid, dat%u_nodal_data_id, start,
     &        kount, dat%u_nodal_data)
         CALL check_err(iret)
         iret = nf_put_vara_DOUBLE(dat%ncid, dat%v_nodal_data_id, start,
     &          kount, dat%v_nodal_data)
         CALL check_err(iret)
      CASE DEFAULT
         WRITE(*,*) "ADCIRC's NetCDF cannot write ",
     &   descript%num_items_per_record, " dimensional arrays."
      END SELECT
!
!     Close netCDF file
9999  iret = nf_close(dat%ncid)
      CALL check_err(iret)

#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from writeNodalData."
#endif
C-----------------------------------------------------------------------
      end subroutine writeNodalData
C-----------------------------------------------------------------------     


C-----------------------------------------------------------------------
C     S U B R O U T I N E    I N I T   N E T C D F   H O T S T A R T
C-----------------------------------------------------------------------
C     jgf49.35 Sets up netCDF variables and writes mesh data into netcdf  
C     hotstart file.
C-----------------------------------------------------------------------
      SUBROUTINE initNetCDFHotstart(lun, Elev1Descript,
     &   Elev2Descript, VelDescript, CH1Descript, EtaDiscDescript,
     &   NodeCodeDescript, NOFFDescript, error) 
      USE GLOBAL, ONLY : OutputDataDescript_t, ICS, G, IM, NE_G, NP_G
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
C
      LOGICAL error
C
      INTEGER iret ! success or failure of the netcdf call     
      CHARACTER(1024) :: att_text
      INTEGER :: varid
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering initNetCDFHotstart."
#endif
      error = .false.
C
C     Point to the correct hotstart block and then allocate 
C     the netcdf data arrays.
      CALL allocateNetCDFHotstart(lun)
C
C     initialize array values with implicit loop
      hs%zeta1%nodal_data = Elev1Descript%initial_value
      hs%zeta2%nodal_data = Elev2Descript%initial_value
      hs%zetad%nodal_data = EtaDiscDescript%initial_value
      hs%vel%u_nodal_data = VelDescript%initial_value
      hs%vel%v_nodal_data = VelDescript%initial_value
      hs%ch1%nodal_data = CH1Descript%initial_value
      hs%nodecodenc%nodal_int = NodeCodeDescript%int_initial_value
      hs%noffnc%nodal_int = NOFFDescript%int_initial_value
C
C     Initialize netCDF hotstart file, creating a new one 
      CALL createNetCDFOutputFile(hs%ncid, hs%myFile, lun, hs%myTime,
     &                            Elev1Descript, error)
      CALL initNetCDFCoord(hs%myMesh)
      CALL defineMeshVariables(hs%ncid, hs%myMesh)
C
C     Z E T A 1
      hs%zeta1%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%zeta1%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf_def_var(hs%ncid, 'zeta1', NF_DOUBLE,
     &       hs%zeta1%nodal_data_rank, hs%zeta1%nodal_data_dims,
     &       hs%zeta1%nodal_data_id)
      CALL check_err(iret)     
      iret = nf_put_att_text(hs%ncid, hs%zeta1%nodal_data_id,
     &      'long_name',45, 
     &      'water surface elevation at previous time step')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%zeta1%nodal_data_id,
     &      'standard_name', 45, 
     &      'water_surface_elevation_at_previous_time step')
      CALL check_err(iret)
      CALL putUnitsAttribute(hs%ncid, hs%zeta1%nodal_data_id,
     &        VAR_METS, VAR_FEET)  
      iret = nf_put_att_DOUBLE(hs%ncid, hs%zeta1%nodal_data_id,
     &         '_FillValue', NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%zeta1%nodal_data_id, 
     &      'positive', 2, VAR_UP)
      CALL check_err(iret)
C
C     Z E T A 2
      hs%zeta2%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%zeta2%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf_def_var(hs%ncid, 'zeta2', NF_DOUBLE,
     &       hs%zeta2%nodal_data_rank, hs%zeta2%nodal_data_dims,
     &       hs%zeta2%nodal_data_id)
      CALL check_err(iret)     
      iret = nf_put_att_text(hs%ncid, hs%zeta2%nodal_data_id,
     &      'long_name',44, 
     &      'water surface elevation at current time step')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%zeta2%nodal_data_id,
     &      'standard_name', 44, 
     &      'water_surface_elevation_at_current_time_step')
      CALL check_err(iret)
      CALL putUnitsAttribute(hs%ncid, hs%zeta2%nodal_data_id,
     &        VAR_METS, VAR_FEET)  
      iret = nf_put_att_DOUBLE(hs%ncid, hs%zeta2%nodal_data_id,
     &        '_FillValue', NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%zeta2%nodal_data_id,
     &      'positive', 2, VAR_UP)
      CALL check_err(iret)
C
C     Z E T A D
      hs%zetad%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%zetad%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf_def_var(hs%ncid, 'zetad', NF_DOUBLE,
     &       hs%zetad%nodal_data_rank, hs%zetad%nodal_data_dims,
     &       hs%zetad%nodal_data_id)
      CALL check_err(iret)     
      iret = nf_put_att_text(hs%ncid, hs%zetad%nodal_data_id,
     &      'long_name',42,       
     &      'water elevation at flux specified boundary')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%zetad%nodal_data_id,
     &      'standard_name', 42, 
     &      'water_elevation_at_flux_specified_boundary')
      CALL check_err(iret)
      CALL putUnitsAttribute(hs%ncid, hs%zetad%nodal_data_id,
     &        VAR_METS, VAR_FEET)  
      iret = nf_put_att_DOUBLE(hs%ncid, hs%zetad%nodal_data_id,
     &        '_FillValue', NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%zetad%nodal_data_id,'positive',
     &       2, VAR_UP)
      CALL check_err(iret)
C
C     U V E L
      hs%vel%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%vel%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf_def_var(hs%ncid, 'u-vel', NF_DOUBLE,
     &       hs%vel%nodal_data_rank, hs%vel%nodal_data_dims,
     &       hs%vel%nodal_data_id)
      CALL check_err(iret)
      if (ics.eq.2) then
         iret = nf_put_att_text(hs%ncid, hs%vel%u_nodal_data_id, 
     &          'long_name', 32, 
     &          'vertically averaged e/w velocity')
         CALL check_err(iret)
         iret = nf_put_att_text(hs%ncid, hs%vel%u_nodal_data_id,
     &          'positive', 4, VAR_EAST)
      else
         iret = nf_put_att_text(hs%ncid, hs%vel%u_nodal_data_id,
     &          'long_name', 43,
     &          'vertically averaged velocity in x-direction')
         CALL check_err(iret)
         iret = nf_put_att_text(hs%ncid, hs%vel%u_nodal_data_id,
     &          'positive', 5, VAR_RG)
      endif
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%vel%u_nodal_data_id, 
     &       'standard_name', 10, 'u_velocity')
      CALL check_err(iret)
      CALL putUnitsAttribute(hs%ncid, hs%vel%u_nodal_data_id,
     &        VAR_MPS, VAR_FPS)  
      CALL check_err(iret)
      iret = nf_put_att_double(hs%ncid, hs%vel%u_nodal_data_id,
     &        '_FillValue', NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(hs%ncid, hs%vel%u_nodal_data_id,
     &      'dry_Value', NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
C
C     V V E L   
      iret = nf_def_var(hs%ncid, 'v-vel', NF_DOUBLE,
     &       hs%vel%nodal_data_rank, hs%vel%nodal_data_dims,
     &       hs%vel%nodal_data_id)
      CALL check_err(iret)     
      if (ics.eq.2) then
         iret = nf_put_att_text(hs%ncid, hs%vel%v_nodal_data_id,
     &           'long_name', 32, 
     &           'vertically averaged n/s velocity')
         CALL check_err(iret)
         iret = nf_put_att_text(hs%ncid, hs%vel%v_nodal_data_id,
     &           'positive', 5, VAR_NORTH)
      else
         iret = nf_put_att_text(hs%ncid, hs%vel%v_nodal_data_id,
     &           'long_name', 43, 
     &           'vertically averaged velocity in y-direction')
         CALL check_err(iret)
         iret = nf_put_att_text(hs%ncid, hs%vel%v_nodal_data_id,
     &           'positive', 49, VAR_UPV)
      endif
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%vel%v_nodal_data_id,
     &            'standard_name', 10, 'v_velocity')
      CALL check_err(iret)
      CALL putUnitsAttribute(hs%ncid, hs%vel%v_nodal_data_id,
     &        VAR_MPS, VAR_FPS)  
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(hs%ncid, hs%vel%v_nodal_data_id,
     &            '_FillValue', NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(hs%ncid, hs%vel%v_nodal_data_id,
     &            'dry_Value', NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
C
C     C H 1
      IF (IM.EQ.10) THEN
         hs%ch1%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
         hs%ch1%nodal_data_dims(2) = hs%myTime%timenc_dim_id
         iret = nf_def_var(hs%ncid, 'ch1', NF_DOUBLE,
     &       hs%ch1%nodal_data_rank, hs%ch1%nodal_data_dims,
     &       hs%ch1%nodal_data_id)
         CALL check_err(iret)     
         att_text = "concentration"
         iret = nf_put_att_text(hs%ncid, hs%ch1%nodal_data_id,
     &      'long_name',len_trim(att_text), trim(att_text))
         CALL check_err(iret)
         att_text = "concentration" 
         iret = nf_put_att_text(hs%ncid, hs%ch1%nodal_data_id,
     &      'standard_name', len_trim(att_text), trim(att_text))
         CALL check_err(iret)
         iret = nf_put_att_DOUBLE(hs%ncid, hs%ch1%nodal_data_id,
     &        '_FillValue', NF_DOUBLE, 1, doubleval)
         CALL check_err(iret)
      ENDIF
C
C     N O D E C O D E
      hs%nodecodenc%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%nodecodenc%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf_def_var(hs%ncid, 'nodecode', NF_INT,
     &       hs%nodecodenc%nodal_data_rank, 
     &       hs%nodecodenc%nodal_data_dims,
     &       hs%nodecodenc%nodal_data_id)
      CALL check_err(iret)     
      att_text = "wet or dry state of node where 1 indicates that the 
     &node is wet and 0 indicates that the node is dry"
      iret = nf_put_att_text(hs%ncid, hs%nodecodenc%nodal_data_id,
     &      'long_name',len_trim(att_text), trim(att_text))
      CALL check_err(iret)
      att_text = "node_wet_or_dry"
      iret = nf_put_att_text(hs%ncid, hs%nodecodenc%nodal_data_id,
     &      'standard_name', len_trim(att_text), trim(att_text))
      CALL check_err(iret)
C
C     N O F F
      hs%noffnc%nodal_data_dims(1) = hs%myMesh%num_elems_dim_id
      hs%noffnc%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf_def_var(hs%ncid, 'noff', NF_INT,
     &       hs%noffnc%nodal_data_rank, 
     &       hs%noffnc%nodal_data_dims,
     &       hs%noffnc%nodal_data_id)
      CALL check_err(iret)  
      att_text = "wet or dry state of element where 1 indicates that  
     &the element is wet and 0 indicates that it is dry"
      iret = nf_put_att_text(hs%ncid, hs%nodecodenc%nodal_data_id,
     &      'long_name',len_trim(att_text), trim(att_text))
      CALL check_err(iret)
      att_text = "element_wet_or_dry"
      iret = nf_put_att_text(hs%ncid, hs%nodecodenc%nodal_data_id,
     &      'standard_name', len_trim(att_text), trim(att_text))
      CALL check_err(iret)
!
!     Define hotstart parameters
      CALL check_err(iret)
      iret = nf_def_var(hs%ncid, 'imhs', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 18, 
     &                       'model type')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 10,
     &                'model_type')
      CALL check_err(iret)
C
      iret = nf_def_var(hs%ncid, 'iths', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 59, 
     &       'model time step number since the beginning
     & of the model run')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 15,
     &       'model_time_step')
      CALL check_err(iret)
C
      iret = nf_def_var(hs%ncid, 'iestp', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 176, 
     &                       'line number (for ASCII output) or record
     & number (for binary output) of the most recent entry in the
     & elevation time series at specified elevation recording stations
     & output file')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 54,
     &      'line/record_number_of_last_entry_in_elev_rec_stations')
      CALL check_err(iret)
c
      iret = nf_def_var(hs%ncid, 'nscoue', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 150, 
     &                       'time step counter to determine when the
     & next entry will be written to the elevation time series at
     & specified elevation recording Stations output file')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 15,
     &           'time_step_counter_for_next_entry_elev_rec_stations)')
      CALL check_err(iret)
C
      iret = nf_def_var(hs%ncid, 'ivstp', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 190, 
     &                       ' line number (for ASCII output) or record 
     & number (for binary output) of the most recent entry in the
     & depth-averaged velocity time series at specified velocity 
     & recording stations output file')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 48,
     &               'record_number_of_last_entry_in_vel_rec_stations')
      CALL check_err(iret)
C
      iret = nf_def_var(hs%ncid, 'nscouv', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 165, 
     &                       'time step counter to determine when the
     & next entry will be written to the Depth-averaged Velocity Time
     & series at specified velocity recording stations output file.')
      CALL check_err(iret)

      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 50,
     &           'time_step_counter_for_next_entry_vel_rec_stations')
      CALL check_err(iret)
C
      iret = nf_def_var(hs%ncid, 'icstp', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 190, 
     &                       'line number (for ASCII output) or record
     & number (for binary output) of the most recent entry in the scalar
     & concentration time series at specified concentration recording
     & stations output file')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 49,
     &               'record_number_of_last_entry_in_conc_rec_stations')
      CALL check_err(iret)
C
      iret = nf_def_var(hs%ncid, 'nscouc', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 165, 
     &                       'time step counter to determine when the
     & next entry will be written to the scalar concentration Time 
     &series at specified concentration Recording Stations output file')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 51,
     &           'time_step_counter_for_next_entry_conc_rec_stations')
      CALL check_err(iret)
C
      iret = nf_def_var(hs%ncid, 'ipstp', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 192, 
     &                       'line number (for ASCII output) or record
     & number (for binary output) of the most recent entry in the 
     &atmospheric pressure time series at specified meteorological 
     &recording stations')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 57,
     &       'record_number_of_last_entry_of_atm_press_at_rec_stations')
      CALL check_err(iret)
C
      iret = nf_def_var(hs%ncid, 'iwstp', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 190, 
     &                       'line number (for ASCII output) or record
     & number (for binary output) of the most recent entry in the Wind 
     & Velocity time series at specified meteorological recording 
     &stations')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 56,
     &       'record_number_of_last_entry_of_wind_vel_at_rec_stations')
      CALL check_err(iret)
C
      iret = nf_def_var(hs%ncid, 'nscoum', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 244, 
     &                       'time step counter to determine when the
     & next entry will be written to the atmospheric pressure time 
     &series at specified meteorological recording stations and wind 
     &velocity time series at specified meteorological recording 
     &stations output files')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 60,
     &   'time_step_counter_of_atm_press_and_wind_vel_at_rec_stations')
      CALL check_err(iret)
C
      iret = nf_def_var(hs%ncid, 'igep', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 165, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &Elevation time series at All Nodes in the Model Grid output file')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 51,
     &       'record_number_of_last_entry_of_elev_at_model_nodes')
      CALL check_err(iret)
C
      iret = nf_def_var(hs%ncid, 'nscouge', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 139, 
     &                       'time step counter to determine when the
     & next entry will be written to the  Elevation time series at All 
     &Nodes in the Model Grid output file')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 41,
     &   'time_step_counter_of_elev_at_model_nodes')
      CALL check_err(iret)
C
      iret = nf_def_var(hs%ncid, 'igvp', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 179, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &Depth-averaged velocity time series at All Nodes in the Model 
     &Grid output file')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 50,
     &       'record_number_of_last_entry_of_vel_at_model_nodes')
      CALL check_err(iret)
C
      iret = nf_def_var(hs%ncid, 'nscougv', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 153, 
     &                       'time step counter to determine when the 
     &next entry will be written to the Depth-averaged velocity time 
     &series at All Nodes in the Model Grid output file')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 40,
     &   'time_step_counter_of_vel_at_model_nodes')
      CALL check_err(iret)
C
      iret = nf_def_var(hs%ncid, 'igcp', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 176, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &scalar Concentration time series at All Nodes in the Model 
     &Grid output file')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 51,
     &       'record_number_of_last_entry_of_conc_at_model_nodes')
      CALL check_err(iret)
C
      iret = nf_def_var(hs%ncid, 'nscougc', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 150, 
     &                       'time step counter to determine when the 
     &next entry will be written to the scalar Concentration time 
     &series at All Nodes in the Model Grid output file')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 41,
     &   'time_step_counter_of_conc_at_model_nodes')
      CALL check_err(iret)
C
      iret = nf_def_var(hs%ncid, 'igpp', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 176, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &atmospheric pressure time series at All Nodes in the Model 
     &Grid output file')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 61,
     &       'record_number_of_last_entry_of_atm_press_at_model_nodes')
      CALL check_err(iret)
C
      iret = nf_def_var(hs%ncid, 'igwp', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 179, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &Wind Stress or velocity time series at All Nodes in the Model 
     &Grid output file')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 55,
     &       'record_number_of_last_entry_of_wind_vel_at_model_nodes')
      CALL check_err(iret)
C
      iret = nf_def_var(hs%ncid, 'nscougw', NF_INT, 0, 1, varid)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'long_name', 222, 
     &                       'time step counter to determine when the 
     &next entry will be written to the atmospheric pressure time  
     &series at All Nodes in the Model Grid and Wind Stress or velocity
     &time series at All Nodes in the Model Grid output files')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, varid, 'standard_name', 59,
     &   'time_step_counter_of_atm_press_and_wind_vel_at_model_nodes')
      CALL check_err(iret)
C
C     define metadata and selected fort.15 parameters in netcdf file
      CALL defineMetaData(hs%ncid)
C
C     Leave define mode
      iret = nf_enddef(hs%ncid)
      CALL check_err(iret)
C
C     write mesh to netcdf file
      CALL putMeshVariables(hs%ncid,hs%myMesh)
C
C     now close the initialized netcdf file
      iret = nf_close(hs%ncid)
      CALL check_err(iret)
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from initNetCDFHotstart."
#endif
C-----------------------------------------------------------------------
      END SUBROUTINE initNetCDFHotstart
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   A L L O C A T E   N E T C D F   H O T S T A R T
C-----------------------------------------------------------------------
C     jgf49.41 Allocates memory for netcdf hotstart data.
C-----------------------------------------------------------------------
      SUBROUTINE allocateNetCDFHotstart(lun)
      USE GLOBAL, ONLY : NE_G, NP_G, IM, IMHS
      IMPLICIT NONE
      INTEGER, intent(in) :: lun ! unit number of hotstart file (67 or 68)
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering allocateNetCDFHotstart."
#endif
C
C     Point to the hotstart file we want to work on.
      IF (lun.eq.67) THEN
         hs => hs67
      ELSE 
         hs => hs68
      ENDIF
C
C     nodes and elements
      hs%myMesh%num_nodes = NP_G
      hs%myMesh%num_elems = NE_G
      hs%myMesh%nface_len = 3
      ALLOCATE(hs%myTime%timenc(hs%myTime%timenc_len))
C
C     Elev1
      ALLOCATE(hs%zeta1%nodal_data(hs%myMesh%num_nodes,
     &       hs%myTime%timenc_len))
C     Elev2
      ALLOCATE(hs%zeta2%nodal_data(hs%myMesh%num_nodes,
     &       hs%myTime%timenc_len))
C     EtaDisc
      ALLOCATE(hs%zetad%nodal_data(hs%myMesh%num_nodes,
     &       hs%myTime%timenc_len))
C     Vel
      ALLOCATE(hs%vel%u_nodal_data(hs%myMesh%num_nodes,
     &       hs%myTime%timenc_len))
      ALLOCATE(hs%vel%v_nodal_data(hs%myMesh%num_nodes,
     &       hs%myTime%timenc_len))
C     CH1
      IF ((IM.EQ.10).OR.(IMHS.EQ.10)) THEN
         ALLOCATE(hs%ch1%nodal_data(hs%myMesh%num_nodes,
     &          hs%myTime%timenc_len))
      ENDIF
C     NodeCode 
      ALLOCATE(hs%nodecodenc%nodal_int(hs%myMesh%num_nodes,
     &       hs%myTime%timenc_len))
C     NOFF 
      ALLOCATE(hs%noffnc%nodal_int(hs%myMesh%num_elems,
     &       hs%myTime%timenc_len))
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from allocateNetCDFHotstart."
#endif
C-----------------------------------------------------------------------
      END SUBROUTINE allocateNetCDFHotstart
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E    W R I T E   N E T C D F   H O T S T A R T   
C-----------------------------------------------------------------------
C     jgf49.17.02 Writes data to hotstart file.
C-----------------------------------------------------------------------
      SUBROUTINE writeNetCDFHotstart(lun, Elev1Descript,
     &        Elev2Descript, VelDescript, CH1Descript, EtaDiscDescript,
     &        NodeCodeDescript, NOFFDescript, time, it)
      USE GLOBAL, ONLY : OutputDataDescript_t, MNPROC, 
     &                   im, iestp, nscoue, ivstp,
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
      REAL(8), intent(in) :: time
      INTEGER, intent(in) :: it   ! current ADCIRC time step
C
      INTEGER i,j
      INTEGER counti(1), starti(1), n
      INTEGER kount(2), start(2) ! for nodally based data
      INTEGER elekount(2) ! for elementally based data
      INTEGER iret ! success or failure of the netcdf call
      INTEGER recid    
      INTEGER varid
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering writeNetCDFHotstart."
#endif
C
C     Point to the hotstart file we want to work on
      IF (lun.eq.67) THEN
         hs => hs67
      ELSE 
         hs => hs68
      ENDIF
C 
      CALL setNetCDFFileName(lun,hs%myFile)
C
      iret = nf_open(hs%myFile%FILENAME, NF_WRITE, hs%ncid)
      CALL check_err(iret)   
!   
!     Inquire variables
      iret=nf_inq_unlimdim(hs%ncid,recid)
      CALL check_err(iret)   
      iret=nf_inq_dimlen(hs%ncid,recid,hs%myFile%record_counter)
      CALL check_err(iret)   
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: hs%myFile%record_counter = ",
     &            hs%myFile%record_counter
#endif
C
C     Don't increment the record counter so that we can store data at the
C     next location in the netcdf file -- the hotstart file is only intended
C     to have a single snapshot of data in it.
C      hs%myFile%record_counter = hs%myFile%record_counter + 1       
      hs%myFile%record_counter = 1       
!     
!     Store time
      iret=nf_inq_varid(hs%ncid,"time",hs%myTime%timenc_id)
      starti(1)=hs%myFile%record_counter
      counti(1)=hs%myTime%timenc_len
      hs%myTime%timenc(hs%myTime%timenc_len)=time
      iret = nf_put_vara_DOUBLE(hs%ncid, hs%myTime%timenc_id,
     &       starti,counti, hs%myTime%timenc)
      CALL check_err(iret)
C
      kount(1)=hs%myMesh%num_nodes
      kount(2)=hs%myTime%timenc_len
      elekount(1)=hs%myMesh%num_elems
      elekount(2)=kount(2)
      start(1)=1
      start(2)=hs%myFile%record_counter
      hs%myTime%timenc(hs%myTime%timenc_len)=time
C
C     Get the NetCDF IDs of the relevant variables from the file
      iret=nf_inq_varid(hs%ncid,"zeta1",hs%zeta1%nodal_data_id)
      CALL check_err(iret)
      iret=nf_inq_varid(hs%ncid,"zeta2",hs%zeta2%nodal_data_id)
      CALL check_err(iret)
      iret=nf_inq_varid(hs%ncid,"zetad",hs%zetad%nodal_data_id)
      CALL check_err(iret)
      iret=nf_inq_varid(hs%ncid,"u-vel",hs%vel%u_nodal_data_id)
      CALL check_err(iret)
      iret=nf_inq_varid(hs%ncid,"v-vel",hs%vel%v_nodal_data_id)
      CALL check_err(iret)
      iret=nf_inq_varid(hs%ncid,"nodecode",hs%nodecodenc%nodal_data_id)
      CALL check_err(iret)
      iret=nf_inq_varid(hs%ncid,"noff",hs%noffnc%nodal_data_id)
      CALL check_err(iret)
C
C     Place the nodal data into the netcdf arrays
      IF (MNPROC.eq.1) THEN
         hs%zeta1%nodal_data(:,kount(2)) = Elev1Descript%array 
         hs%zeta2%nodal_data(:,kount(2)) = Elev2Descript%array 
         hs%zetad%nodal_data(:,kount(2)) = EtaDiscDescript%array 
         hs%vel%u_nodal_data(:,kount(2)) = VelDescript%array 
         hs%vel%v_nodal_data(:,kount(2)) = VelDescript%array2 
         hs%nodecodenc%nodal_int(:,kount(2)) = NodeCodeDescript%iarray 
         hs%noffnc%nodal_int(:,elekount(2)) = NOFFDescript%iarray 
      ELSE
         hs%zeta1%nodal_data(:,kount(2)) = Elev1Descript%array_g 
         hs%zeta2%nodal_data(:,kount(2)) = Elev2Descript%array_g 
         hs%zetad%nodal_data(:,kount(2)) = EtaDiscDescript%array_g 
         hs%vel%u_nodal_data(:,kount(2)) = VelDescript%array_g
         hs%vel%v_nodal_data(:,kount(2)) = VelDescript%array2_g
         hs%nodecodenc%nodal_int(:,kount(2)) =
     &                                  NodeCodeDescript%iarray_g 
         hs%noffnc%nodal_int(:,elekount(2)) = NOFFDescript%iarray_g 
      ENDIF
C
C     Write the nodal data to the netcdf file 
      iret = nf_put_vara_DOUBLE(hs%ncid, hs%zeta1%nodal_data_id, start,
     &       kount, hs%zeta1%nodal_data)
      CALL check_err(iret)
      iret = nf_put_vara_DOUBLE(hs%ncid, hs%zeta2%nodal_data_id, start,
     &       kount, hs%zeta2%nodal_data)
      CALL check_err(iret)
      iret = nf_put_vara_DOUBLE(hs%ncid, hs%zetad%nodal_data_id, start,
     &       kount, hs%zetad%nodal_data)
      CALL check_err(iret)
      iret = nf_put_vara_DOUBLE(hs%ncid, hs%vel%u_nodal_data_id, start,
     &       kount, hs%vel%u_nodal_data)
      CALL check_err(iret)
      iret = nf_put_vara_DOUBLE(hs%ncid, hs%vel%v_nodal_data_id, start,
     &       kount, hs%vel%v_nodal_data)
      CALL check_err(iret)
      iret = nf_put_vara_INT(hs%ncid, hs%nodecodenc%nodal_data_id,
     &      start, kount, hs%nodecodenc%nodal_int)
      CALL check_err(iret)
      iret = nf_put_vara_INT(hs%ncid, hs%noffnc%nodal_data_id,
     &       start, elekount, hs%noffnc%nodal_int)
      CALL check_err(iret)
C
C     Get each variable ID for the model parameters in the netcdf file
C     and immediately write the parameter value to that variable ID before
C     going on to the next one. 
      iret = nf_inq_varid(hs%ncid,"imhs",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, im)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"iths",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, it)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"iestp",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, iestp)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"nscoue",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, nscoue)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"ivstp",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, ivstp)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"nscouv",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, nscouv)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"icstp",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, icstp)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"nscouc",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, nscouc)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"ipstp",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, ipstp)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"iwstp",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, iwstp)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"nscoum",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, nscoum)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"igep",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, igep)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"nscouge",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, nscouge)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"igvp",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, igvp)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"nscougv",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, nscougv)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"igcp",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, igcp)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"nscougc",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, nscougc)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"igpp",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, igpp)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"igwp",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, igwp)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"nscougw",varid)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, varid, nscougw)
      CALL check_err(iret)
C
C     now close the netcdf file
      iret = nf_close(hs%ncid)
      CALL check_err(iret)
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from writeNetCDFHotstart."
#endif
C-----------------------------------------------------------------------
      END SUBROUTINE writeNetCDFHotstart
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E    R E A D   N E T C D F   H O T S T A R T   
C-----------------------------------------------------------------------
C     jgf49.17.02 Reads data from the hotstart file.
C-----------------------------------------------------------------------
      SUBROUTINE readNetCDFHotstart(lun)
      USE GLOBAL, ONLY : OutputDataDescript_t, FileFmtRev, FileFmtMinor,
     &                   FileFmtMajor, imhs, iths, iestp, nscoue, ivstp,
     &                   nscouv, icstp, nscouc, ipstp, iwstp, nscoum,
     &                   igep, nscouge, igvp, nscougv, igcp, nscougc,
     &                   igpp, igwp, nscougw, ETA1, ETA2, EtaDisc, 
     &                   UU2, VV2, NNODECODE, NOFF, NODES_LG,
     &                   IMAP_EL_LG, MNPROC, MYPROC, NE, NP
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      IMPLICIT NONE
C
      INTEGER, intent(in) :: lun
C      
      INTEGER i,j
      INTEGER counti(1), starti(1), n
      INTEGER kount(2), start(2)
      INTEGER iret ! success or failure of the netcdf call
      INTEGER recid
      INTEGER fd_node_number
      INTEGER sd_node_number
      INTEGER fd_ele_number
      INTEGER sd_ele_number
      INTEGER varid
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering readNetCDFHotstart."
#endif
C
C     Point to the hotstart file we want to work on
      IF (lun.eq.67) THEN
         hs => hs67
      ELSE 
         hs => hs68
      ENDIF
C
      hs%myFile%fileFound = .false.
      CALL setNetCDFFileName(lun,hs%myFile)
C
C     Open fulldomain file
      INQUIRE(FILE=hs%myFile%FILENAME,EXIST=hs%myFile%fileFound)
      IF (hs%myFile%fileFound.eqv..false.) THEN
         WRITE(*,*) "ERROR: The file ",hs%myFile%FILENAME,
     &      " was not found; ADCIRC terminating."
#ifdef CMPI
            CALL MSG_FINI()
#endif
         STOP
      ELSE
         iret = nf_open(hs%myFile%FILENAME, NF_NOWRITE, hs%ncid)
         CALL check_err(iret) 
      ENDIF
! 
!     Inquire variables
      iret=nf_inq_unlimdim(hs%ncid,recid)
      CALL check_err(iret)   
      iret=nf_inq_dimlen(hs%ncid,recid,hs%myFile%record_counter)
      CALL check_err(iret)   
C
      iret=nf_inq_dimid(hs%ncid, "node", hs%myMesh%num_nodes_dim_id)
      CALL check_err(iret)   
      iret=nf_inq_dimlen(hs%ncid, hs%myMesh%num_nodes_dim_id,
     &     hs%myMesh%num_nodes)
      CALL check_err(iret)   
      iret=nf_inq_dimid(hs%ncid, "nele", hs%myMesh%num_elems_dim_id)
      CALL check_err(iret)   
      iret=nf_inq_dimlen(hs%ncid, hs%myMesh%num_elems_dim_id,
     &     hs%myMesh%num_elems)
      CALL check_err(iret)   
C
      CALL allocateNetCDFHotstart(lun)
C
      hs%myFile%record_counter = 1
      kount(1)=hs%myMesh%num_nodes
      kount(2)=hs%myTime%timenc_len
      start(1)=1
      start(2)=hs%myFile%record_counter
!   
!     Get time
      iret=nf_inq_varid(hs%ncid, "time", hs%myTime%timenc_id)
      CALL check_err(iret)
      starti(1)=hs%myFile%record_counter
      counti(1)=hs%myTime%timenc_len
      iret = nf_get_vara_DOUBLE(hs%ncid, hs%myTime%timenc_id, starti, 
     &       counti, hs%myTime%timenc)
      CALL check_err(iret)
C
C     Read in fulldomain data
C     Elev1
      iret=nf_inq_varid(hs%ncid, "zeta1", hs%zeta1%nodal_data_id)
      CALL check_err(iret)
      iret = nf_get_vara_DOUBLE(hs%ncid, hs%zeta1%nodal_data_id, start,
     &       kount, hs%zeta1%nodal_data)
      CALL check_err(iret)
C     Elev2
      iret=nf_inq_varid(hs%ncid, "zeta2", hs%zeta2%nodal_data_id)
      CALL check_err(iret)
      iret = nf_get_vara_DOUBLE(hs%ncid, hs%zeta2%nodal_data_id, start,
     &       kount, hs%zeta2%nodal_data)
      CALL check_err(iret)
C     EtaDisc
      iret=nf_inq_varid(hs%ncid, "zetad", hs%zetad%nodal_data_id)
      CALL check_err(iret)
      iret = nf_get_vara_DOUBLE(hs%ncid, hs%zetad%nodal_data_id, start,
     &       kount, hs%zetad%nodal_data)
      CALL check_err(iret)
C     Vel
      iret=nf_inq_varid(hs%ncid, "u-vel", hs%vel%u_nodal_data_id)
      CALL check_err(iret)
      iret = nf_get_vara_DOUBLE(hs%ncid, hs%vel%u_nodal_data_id, start,
     &       kount, hs%vel%u_nodal_data)
      CALL check_err(iret)
C
      iret=nf_inq_varid(hs%ncid, "v-vel", hs%vel%v_nodal_data_id)
      CALL check_err(iret)
      iret = nf_get_vara_DOUBLE(hs%ncid, hs%vel%v_nodal_data_id, start,
     &       kount, hs%vel%v_nodal_data)
      CALL check_err(iret)
C     NodeCode
      iret=nf_inq_varid(hs%ncid, "nodecode",hs%nodecodenc%nodal_data_id)
      CALL check_err(iret)
      iret = nf_get_vara_INT(hs%ncid, hs%nodecodenc%nodal_data_id, 
     &       start, kount, hs%nodecodenc%nodal_int)
      CALL check_err(iret)
C     NOFF
      start(1)=1
      kount(1)=hs%myMesh%num_elems
      iret=nf_inq_varid(hs%ncid, "noff", hs%noffnc%nodal_data_id)
      CALL check_err(iret)
      iret = nf_get_vara_INT(hs%ncid, hs%noffnc%nodal_data_id, start,
     &       kount, hs%noffnc%nodal_int)
      CALL check_err(iret)
C
C     Place nodal data in ADCIRC arrays 

      ! serial
      IF (MNPROC.eq.1) THEN 
         ETA1(:) = hs%zeta1%nodal_data(:,kount(2))
         ETA2(:) = hs%zeta2%nodal_data(:,kount(2))
         UU2(:) = hs%vel%u_nodal_data(:,kount(2))
         VV2(:) = hs%vel%v_nodal_data(:,kount(2))
         EtaDisc(:) = hs%zetad%nodal_data(:,kount(2))
         NNODECODE(:) = hs%nodecodenc%nodal_int(:,kount(2))
         NOFF(:) = hs%noffnc%nodal_int(:,kount(2))
      ! parallel
      ELSE 
         ! loop over subdomain nodes
         DO sd_node_number=1,np
            ! get the corresponding fulldomain node number 
            fd_node_number = ABS(NODES_LG(sd_node_number)) 
            ! fill in the subdomain arrays with the corresponding
            ! fulldomain values
            ETA1(sd_node_number) = 
     &          hs%zeta1%nodal_data(fd_node_number,kount(2))
            ETA2(sd_node_number) = 
     &          hs%zeta2%nodal_data(fd_node_number,kount(2))
            UU2(sd_node_number) = 
     &          hs%vel%u_nodal_data(fd_node_number,kount(2))
            VV2(sd_node_number) = 
     &          hs%vel%v_nodal_data(fd_node_number,kount(2))
            EtaDisc(sd_node_number) = 
     &          hs%zetad%nodal_data(fd_node_number,kount(2))
            NNODECODE(sd_node_number) = 
     &          hs%nodecodenc%nodal_int(fd_node_number,kount(2))
         ENDDO   
         ! loop over subdomain elements
         DO sd_ele_number=1,ne
            fd_ele_number = ABS(IMAP_EL_LG(sd_ele_number))
            NOFF(sd_ele_number) = 
     &          hs%noffnc%nodal_int(fd_ele_number,kount(2))
         ENDDO
      ENDIF
C
C     Read in model parameters to ADCIRC variables
      iret = nf_inq_varid(hs%ncid,"imhs",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, imhs)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"iths",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, iths)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"iestp",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, iestp)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"nscoue",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, nscoue)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"ivstp",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, ivstp)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"nscouv",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, nscouv)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"icstp",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, icstp)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"nscouc",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, nscouc)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"ipstp",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, ipstp)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"iwstp",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, iwstp)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"nscoum",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, nscoum)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"igep",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, igep)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"nscouge",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, nscouge)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"igvp",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, igvp)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"nscougv",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, nscougv)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"igcp",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, igcp)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"nscougc",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, nscougc)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"igpp",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, igpp)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"igwp",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, igwp)
      CALL check_err(iret)
      iret = nf_inq_varid(hs%ncid,"nscougw",varid)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, varid, nscougw)
      CALL check_err(iret)
C  
      DEALLOCATE(hs%zeta1%nodal_data)
      DEALLOCATE(hs%zeta2%nodal_data)
      DEALLOCATE(hs%zetad%nodal_data)
      DEALLOCATE(hs%vel%u_nodal_data)
      DEALLOCATE(hs%vel%v_nodal_data)
      IF (IMHS.EQ.10) THEN
         DEALLOCATE(hs%CH1%nodal_data)
      ENDIF
      DEALLOCATE(hs%nodecodenc%nodal_int)
      DEALLOCATE(hs%noffnc%nodal_int)
      DEALLOCATE(hs%myTime%timenc)
C
C     now close the netcdf file
      iret = nf_close(hs%ncid)
      CALL check_err(iret)
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from readNetCDFHotstart."
#endif
C-----------------------------------------------------------------------
      END SUBROUTINE readNetCDFHotstart
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   D E F I N E   M E S H   V A R I A B L E S 
C-----------------------------------------------------------------------
C     jgf49.17.02 Defines data that are common to all netCDF files.
C-----------------------------------------------------------------------
      subroutine defineMeshVariables(ncid, myMesh)
      USE GLOBAL, ONLY : RAD2DEG, SLAM0, SFEA0, NM
      USE SIZES, ONLY : SZ
      USE NodalAttributes, ONLY : nwp, nolibf, tau0, cf, eslm
      USE VERSION
      IMPLICIT NONE
      INTEGER iret   ! Error status return
      type(meshStructure) :: myMesh
      INTEGER, intent(in) :: ncid   ! netCDF id
      INTEGER i,j
      INTEGER counti(1), starti(1), n
      INTEGER kount(2), start(2)
      REAL    realval(1)
      CHARACTER(1024) :: att_text ! reused to hold attribute text
!     -----------------
!      date_string variables for time attribute
!     -----------------
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer values(8)
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering defineMeshVariables."
#endif
      myMesh%nmnc=NM
!     Switch order in array for NETCDF
       write(*,*)  myMesh%num_elems
       write(*,*)  myMesh%nface_len
      do i=1, myMesh%num_elems
         do j=1, myMesh%nface_len
            myMesh%element(j,i) = myMesh%nmnc(i,j)
         end do
      end do
!       write(99,*) myMesh%nmnc
!       write(98,*) NM
!       write(97,*) myMesh%element
!      stop
! 
!     Define dimensions
      iret = nf_def_dim(ncid,'node',myMesh%num_nodes,
     &       myMesh%num_nodes_dim_id)
      CALL check_err(iret)
      iret = nf_def_dim(ncid,'nele',myMesh%num_elems,
     &       myMesh%num_elems_dim_id)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'nvertex', 3, myMesh%nface_dim_id)
      call check_err(iret)
      if(myMesh%nopenc.ne.0) then
         iret = nf_def_dim(ncid,'nope', myMesh%nopenc,
     &         myMesh%nopenc_dim_id)
         call check_err(iret)
         iret = nf_def_dim(ncid,'neta', myMesh%netanc,
     &         myMesh%netanc_dim_id)
         call check_err(iret)
      endif
      if(myMesh%nbounc.ne.0) then
         iret = nf_def_dim(ncid, 'nbou',
     &          myMesh%nbounc, myMesh%nbounc_dim_id)
         call check_err(iret)
         iret = nf_def_dim(ncid, 'nvel',
     &          myMesh%nvelnc, myMesh%nvelnc_dim_id)
         call check_err(iret)
      endif
! 
!     Define variables
!     Define X
      myMesh%X_dims(1) = myMesh%num_nodes_dim_id
      iret = nf_def_var(ncid, 'x', NF_DOUBLE,
     &       myMesh%X_rank, myMesh%X_dims, myMesh%X_id)
      CALL check_err(iret)
!     Define Y coordinate
      myMesh%Y_dims(1) = myMesh%num_nodes_dim_id
      iret = nf_def_var(ncid, 'y', NF_DOUBLE,
     &       myMesh%Y_rank, myMesh%Y_dims, myMesh%Y_id)
      CALL check_err(iret)
!     Define elements
      myMesh%ELE_dims(1) = myMesh%nface_dim_id
      myMesh%ELE_dims(2) = myMesh%num_elems_dim_id
      iret = nf_def_var(ncid, 'element',NF_INT,
     &       myMesh%ELE_rank, myMesh%ELE_dims, myMesh%ELE_id)
      CALL check_err(iret)
!     Define elevation specified boundary forcing segments information
      if (myMesh%nopenc.ne.0) then
C         iret = nf_def_var(ncid,'nope', NF_INT, 1, 1, myMesh%nopenc_id) 
C         CALL check_err(iret)
C         iret = nf_def_var(ncid,'neta', NF_INT, 1, 1, myMesh%netanc_id) 
C         CALL check_err(iret)
         myMesh%nvdllnc_dims(1) = myMesh%nopenc_dim_id      
         iret = nf_def_var(ncid, 'nvdll',NF_INT, myMesh%nvdllnc_rank, 
     &                  myMesh%nvdllnc_dims, myMesh%nvdllnc_id)
         CALL check_err(iret)
         myMesh%ibtypeenc_dims(1) = myMesh%nopenc_dim_id      
         iret = nf_def_var(ncid, 'ibtypee',NF_INT, 
     &          myMesh%ibtypeenc_rank, myMesh%ibtypeenc_dims,
     &          myMesh%ibtypeenc_id)
         CALL check_err(iret)
         myMesh%nbdvnc_dims(1) = myMesh%nopenc_dim_id
         myMesh%nbdvnc_dims(2) = myMesh%netanc_dim_id      
         iret = nf_def_var(ncid, 'nbdv',NF_INT, myMesh%nbdvnc_rank,
     &          myMesh%nbdvnc_dims, myMesh%nbdvnc_id)
         CALL check_err(iret)
      endif
!     Define normal flow boundary information
      if (myMesh%nbounc.ne.0) then
C         iret = nf_def_var(ncid,'nbou', NF_INT, 1, 1, myMesh%nbounc_id) 
C         CALL check_err(iret)
C         iret = nf_def_var(ncid,'nvel', NF_INT, 1, 1, myMesh%nvelnc_id) 
C         CALL check_err(iret)
         myMesh%nvellnc_dims(1) = myMesh%nbounc_dim_id      
         iret = nf_def_var(ncid, 'nvell',NF_INT, myMesh%nvellnc_rank, 
     &          myMesh%nvellnc_dims, myMesh%nvellnc_id)
         CALL check_err(iret)
         myMesh%ibtypenc_dims(1) = myMesh%nbounc_dim_id      
         iret = nf_def_var(ncid, 'ibtype',NF_INT, myMesh%ibtypenc_rank, 
     &          myMesh%ibtypenc_dims, myMesh%ibtypenc_id)
         CALL check_err(iret)
         myMesh%nbvvnc_dims(1) = myMesh%nbounc_dim_id      
         myMesh%nbvvnc_dims(2) = myMesh%nvelnc_dim_id      
         iret = nf_def_var(ncid, 'nbvv',NF_INT, myMesh%nbvvnc_rank,
     &          myMesh%nbvvnc_dims, myMesh%nbvvnc_id)
         CALL check_err(iret)
      endif
!     -------------------
!     Define Z coordinate
!     --------------------
      myMesh%DEPTH_dims(1) = myMesh%num_nodes_dim_id
      iret = nf_def_var(ncid, 'depth',NF_DOUBLE, myMesh%DEPTH_rank, 
     &       myMesh%DEPTH_dims, myMesh%DEPTH_id)
      CALL check_err(iret)
C
C     Set coordinates as representing latitude or longitude, depending on
C     the value of ICS
      CALL defineCoordinateAttributes(ncid, myMesh%X_id, myMesh%Y_id)
!   
!     Define depth attributes
      iret = nf_put_att_text(ncid, myMesh%DEPTH_id, 'long_name', 19, VAR_L_D)
      iret = nf_put_att_text(ncid, myMesh%DEPTH_id, 'standard_name', 5,VAR_S_D)
      iret = nf_put_att_text(ncid, myMesh%DEPTH_id, 'units', 1, VAR_METS)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, myMesh%DEPTH_id, 'positive', 4, VAR_DW)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, myMesh%ELE_id,
     &   'long_name', 7, 'element')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, myMesh%ELE_id, 'units', 14, VAR_NOND)
      CALL check_err(iret)
!     
!     Define text attributes for boundary segments
      if (myMesh%nopenc.ne.0) then
C        NOPE
         att_text = "number of elevation specified boundary &
     &              forcing segments"
         iret = nf_put_att_text(ncid, myMesh%nopenc_id, 'long_name',
     &          len_trim(att_text), trim(att_text)) 
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, myMesh%nopenc_id, 'units', 14, 
     &                     VAR_NOND)
         CALL check_err(iret)
C        NETA
         att_text = "total number of elevation specified  
     &boundary nodes"
         iret = nf_put_att_text(ncid, myMesh%netanc_id, 'long_name',
     &          len_trim(att_text), trim(att_text)) 
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, myMesh%netanc_id, 'units', 14, 
     &                     VAR_NOND)
         CALL check_err(iret)
C        NVDLL
         att_text = "number of nodes in each elevation specified 
     &boundary segment"
         iret = nf_put_att_text(ncid, myMesh%nvdllnc_id, 'long_name',
     &          len_trim(att_text), trim(att_text))
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, myMesh%nvdllnc_id, 'units', 14, 
     &                      VAR_NOND)
         CALL check_err(iret)
C        IBTYPEE
         att_text = "elevation boundary type"
         iret = nf_put_att_text(ncid, myMesh%ibtypeenc_id, 'long_name',
     &          len_trim(att_text), trim(att_text))
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, myMesh%ibtypeenc_id, 'units', 14, 
     &                      VAR_NOND)
         CALL check_err(iret)
C        NBDV
         att_text = "node numbers on each elevation specified boundary 
     &segment"
         iret = nf_put_att_text(ncid, myMesh%nbdvnc_id, 'long_name',
     &          len_trim(att_text), trim(att_text))
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, myMesh%nbdvnc_id, 'units', 14, 
     &                      VAR_NOND)
         CALL check_err(iret)
      endif
      if (myMesh%nbounc.ne.0) then
C        NBOU
         att_text = "number of normal flow (discharge) specified 
     &boundary segments"
         iret = nf_put_att_text(ncid, myMesh%nbounc_id, 'long_name',
     &          len_trim(att_text), trim(att_text))
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, myMesh%nbounc_id, 'units', 14, 
     &                      VAR_NOND)
         CALL check_err(iret)
C        NVEL
         att_text = "total number of normal flow specified boundary 
     &nodes including both the front and back nodes on internal barrier
     & boundaries"
         iret = nf_put_att_text(ncid, myMesh%nvelnc_id, 'long_name',
     &          len_trim(att_text), trim(att_text))
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, myMesh%nvelnc_id, 'units', 14, 
     &                      VAR_NOND)
         CALL check_err(iret)
C        IBTYPE
         att_text = "type of normal flow (discharge) boundary"
         iret = nf_put_att_text(ncid, myMesh%ibtypenc_id, 'long_name', 
     &       len_trim(att_text), trim(att_text))
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, myMesh%ibtypenc_id, 'units', 14, 
     &                      VAR_NOND)
         CALL check_err(iret)
C        NVELL
         iret = nf_put_att_text(ncid, myMesh%nvellnc_id, 'long_name',62,
     & 'number of nodes in each normal flow specified boundary segment')
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, myMesh%nvellnc_id, 'units', 14, 
     &                      VAR_NOND)
         CALL check_err(iret)
C        NBVV
         att_text = "node numbers on normal flow boundary segment"
         iret = nf_put_att_text(ncid, myMesh%nbvvnc_id, 'long_name',
     &          len_trim(att_text), trim(att_text))
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, myMesh%nbvvnc_id, 'units', 14, 
     &                      VAR_NOND)
         CALL check_err(iret)
      endif
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from defineMeshVariables."
#endif
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
      USE GLOBAL, ONLY : ICS
      IMPLICIT NONE
      INTEGER ncid
      INTEGER xid
      INTEGER yid
      INTEGER iret
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering defineCoordinateAttributes."
#endif
C
      if(ics.eq.2) then  ! this indicates spherical coordinates
      
C        Define longitude attributes
         iret = nf_put_att_text(ncid, xid, 'long_name',9,'longitude')
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, xid,'standard_name', 9,
     &       'longitude')
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, xid,'units', 12, VAR_DEGE)
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, xid, 'positive', 4, VAR_EAST)
         CALL check_err(iret)
!        Define latitude attributes
         iret = nf_put_att_text(ncid, yid, 'long_name', 8, 'latitude')
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, yid,
     &      'standard_name', 8, 'latitude')
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, yid, 'units', 13, VAR_DEGN)
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, yid, 'positive', 5, VAR_NORTH)
         CALL check_err(iret)
      else   ! must be using Cartesian (x,y) coordinates
!        Define x-coordinate attributes
         iret = nf_put_att_text(ncid, xid, 'long_name', 22, VAR_L_X)
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, xid, 'standard_name',
     &                          12, VAR_S_X)
         CALL check_err(iret)
!        determine variable units
         CALL putUnitsAttribute(ncid, xid, VAR_METS, VAR_FEET)
         iret = nf_put_att_text(ncid, xid, 'positive', 5,  VAR_RG)
         CALL check_err(iret)
!        Define y-coordinate attributes
         iret = nf_put_att_text(ncid, yid, 'long_name', 22, VAR_L_Y)
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, yid, 'standard_name', 
     &                          12, VAR_S_Y)
         CALL check_err(iret)
         CALL putUnitsAttribute(ncid, yid, VAR_METS, VAR_FEET)
         iret = nf_put_att_text(ncid, yid, 'positive', 34, VAR_UPY)
         CALL check_err(iret)
      endif

#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from defineCoordinateAttributes."
#endif
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
      USE GLOBAL, ONLY : DP
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
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering putMeshVariables."
#endif
C
C     Store nodal coordinates
      iret = nf_put_var_double(ncid, myMesh%X_id, myMesh%xnc)
      CALL check_err(iret)
      iret = nf_put_var_double(ncid, myMesh%Y_id, myMesh%ync)
      CALL check_err(iret)
C
!     Store depth
      iret = nf_put_var_double(ncid, myMesh%DEPTH_id, DP)
      CALL check_err(iret)
!     Store elements
      kount(1)=myMesh%nface_len
      kount(2)=myMesh%num_elems
      start(1)=1
      start(2)=1
      iret=nf_put_vara_int(ncid,myMesh%ele_id,start,kount,myMesh%element)
      call check_err(iret)
!     Store elevation boundary information
      if(myMesh%nopenc.ne.0) then
C         iret = nf_put_var_int(ncid, myMesh%nopenc_id, myMesh%nopenc)
C         CALL check_err(iret)
C         iret = nf_put_var_int(ncid, myMesh%netanc_id ,myMesh%netanc)
         CALL check_err(iret)
         iret = nf_put_var_int(ncid, myMesh%nvdllnc_id ,myMesh%nvdllnc)
         CALL check_err(iret)
         iret = nf_put_var_int(ncid,myMesh%ibtypeenc_id,
     &          myMesh%ibtypeenc)
         CALL check_err(iret)
         iret = nf_put_var_int(ncid, myMesh%nbdvnc_id ,myMesh%nbdvnc)
         CALL check_err(iret)
      endif
!     Store normal flow boundary information
      if(myMesh%nbounc.ne.0) then
C         iret = nf_put_var_int(ncid, myMesh%nbounc_id, myMesh%nbounc)
C         CALL check_err(iret)
C         iret = nf_put_var_int(ncid, myMesh%nvelnc_id, myMesh%nvelnc)
C         CALL check_err(iret)
         iret = nf_put_var_int(ncid,myMesh%ibtypenc_id, myMesh%ibtypenc)
         CALL check_err(iret)
         iret = nf_put_var_int(ncid, myMesh%nvellnc_id, myMesh%nvellnc)
         CALL check_err(iret)
         iret = nf_put_var_int(ncid, myMesh%nbvvnc_id, myMesh%nbvvnc)
         CALL check_err(iret)
      endif

#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from putMeshVariables."
#endif
C-----------------------------------------------------------------------
      END SUBROUTINE putMeshVariables 
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   D E F I N E   M E T A D A T A 
C-----------------------------------------------------------------------
C     jgf49.17.02 Writes data that are common to all netCDF files
C     to the specified file.
C-----------------------------------------------------------------------
      subroutine defineMetaData(ncid)
      USE VERSION
      USE GLOBAL, ONLY : RUNDES, AGRID, title, institution, source,
     &                   history, references, comments, host,
     &                   convention, contact, dtdp, ihot, 
     &                   ics, nolifa, nolica, nolicat, 
     &                   ncor, ntip, nws, nramp, statim,
     &                   reftim, rnday, dramp, a00, b00, c00, h0,
     &                   cori, ntif, nbfr, RAD2DEG, SLAM0, SFEA0
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
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering defineMetaData."
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
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'model', 6, 'ADCIRC')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'version', 
     &       LEN(TRIM(ADC_VERSION)), ADC_VERSION)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'grid_type', 10, 
     &                      'Triangular')
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'description',
     &       metalength(trim(adjustl(rundes))), trim(adjustl(rundes)))
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'grid',
     &       metalength(trim(adjustl(agrid))), trim(adjustl(agrid)))   
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'title',
     &       metalength(trim(adjustl(title))), trim(adjustl(title)))
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'institution', 
     &       metalength(trim(adjustl(institution))),
     &       trim(adjustl(institution)))
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'source', 
     &       metalength(trim(adjustl(source))), trim(adjustl(source))) 
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'history', 
     &       metalength(trim(adjustl(history))), trim(adjustl(history)))
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'references', 
     &       metalength(trim(adjustl(references))),
     &       trim(adjustl(references)))
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'comments', 
     &       metalength(trim(adjustl(comments))),
     &       trim(adjustl(comments)))
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'host', 
     &       metalength(trim(adjustl(host))), trim(adjustl(host)))    
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'convention', 
     &       metalength(trim(adjustl(convention))),
     &       trim(adjustl(convention)))
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'contact', 
     &       metalength(trim(adjustl(contact))), trim(adjustl(contact)))
      CALL check_err(iret)
      CALL date_and_time(now_date,big_ben,zone,values)
      WRITE(date_string,71) values(1),values(2),values(3),
     &      values(5),values(6),values(7),   (values(4))/60
 71   FORMAT(I4,'-',I2.2,'-',i2.2,' ',i2,':',i2.2,':',i2.2,' '
     &  ,i3.2,':00')
      iret = nf_put_att_text(ncid, NF_GLOBAL,'creation_date',
     &       len_trim(date_string),date_string)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL,'modification_date',
     &       len_trim(date_string),date_string)
      CALL check_err(iret)
!     -------------------------------------------
!     writing global attributes from fort.15 file 
!     -------------------------------------------
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'fort.15', 40, 
     &                      '==== Input File Parameters (below) ====')
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'dt', NF_DOUBLE, 1, dtdp)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ihot', NF_INT, 1, ihot)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ics', NF_INT, 1, ics)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolibf', NF_INT, 1, nolibf)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolifa', NF_INT, 1, nolifa)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolica', NF_INT, 1, nolica)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolicat', NF_INT, 1, 
     &                      nolicat)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nwp', NF_INT, 1, nwp)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ncor', NF_INT, 1, ncor)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ntip', NF_INT, 1, ntip)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nws', NF_INT, 1, nws)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nramp', NF_INT, 1, nramp)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'tau0', NF_DOUBLE, 1, 
     &                         tau0)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'statim', NF_DOUBLE, 1, 
     &                         statim)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'reftim', NF_DOUBLE, 1, 
     &                         reftim)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'rnday', NF_DOUBLE, 1,
     &                         rnday)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'dramp', NF_DOUBLE, 1, 
     &                         dramp)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'a00', NF_DOUBLE, 1, 
     &                         a00)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'b00', NF_DOUBLE, 1, 
     &                         b00)

      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'c00', NF_DOUBLE, 1, 
     &                         c00)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'h0', NF_DOUBLE, 1, 
     &                         h0)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'slam0', NF_DOUBLE, 1, 
     &                         slam0deg)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'sfea0', NF_DOUBLE, 1, 
     &                         sfea0deg)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'cf', NF_DOUBLE, 1, 
     &                         cf)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'eslm', NF_DOUBLE, 1, 
     &                         eslm)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'cori', NF_DOUBLE, 1, 
     &                         cori)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ntif', NF_INT, 1, ntif)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nbfr', NF_INT, 1, nbfr)
      CALL check_err(iret)
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from defineMetaData."
#endif
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
      USE GLOBAL, ONLY : ihot, nramp, rnday, dramp, a00, b00, c00, h0
      USE NodalAttributes, ONLY : nwp, nolibf, tau0, cf, eslm
      USE VERSION
      IMPLICIT NONE
      INTEGER, intent(in) :: ncid 
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
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering updateMetaData."
#endif
C     Open existing NetCDF file
      iret = nf_open(myFile%FILENAME, NF_WRITE, ncid)
      CALL check_err(iret)   
C
      iret = NF_REDEF (ncid)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ihot', NF_INT, 1, ihot)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nramp', NF_INT, 1,nramp)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'tau0', NF_DOUBLE, 1, 
     &                      tau0)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'rnday', NF_DOUBLE, 1,
     &                      rnday)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'dramp', NF_DOUBLE, 1, 
     &                      dramp)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'a00', NF_DOUBLE, 1, 
     &                      a00)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'b00', NF_DOUBLE, 1, 
     &                      b00)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'c00', NF_DOUBLE, 1, 
     &                      c00)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'h0', NF_DOUBLE, 1, 
     &                      h0)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'cf', NF_DOUBLE, 1, 
     &                      cf)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'eslm', NF_DOUBLE, 1, 
     &                      eslm)
      CALL check_err(iret)
      CALL date_and_time(now_date,big_ben,zone,values)
      WRITE(date_string,71) values(1),values(2),values(3),
     &      values(5),values(6),values(7),   (values(4))/60
 71   FORMAT(I4,'-',I2.2,'-',i2.2,' ',i2,':',i2.2,':',i2.2,' '
     &  ,i3.2,':00')
      iret = nf_put_att_text(ncid, NF_GLOBAL,'modification_date',
     &       len_trim(date_string),date_string)
      CALL check_err(iret)
      iret = NF_ENDDEF (ncid)
      CALL check_err(iret)
C     now close the updated netcdf file
      iret = nf_close(ncid)
      CALL check_err(iret)
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from updateMetaData."
#endif
C-----------------------------------------------------------------------
      END SUBROUTINE updateMetaData 
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
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering putUnitsAttribute."
#endif
C
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, var_id, 'units',
     &           len_trim(metric), trim(metric))
      else
          iret = nf_put_att_text(ncid, var_id, 'units',
     &           len_trim(english), trim(english)) 
      endif
      CALL check_err(iret)

#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from putUnitsAttribute."
#endif
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
      USE SIZES, ONLY : SZ
      USE GLOBAL, ONLY : myProc, screenUnit
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      IMPLICIT NONE 
      INTEGER, intent(in) :: iret
#ifdef NETCDF_DEBUG
      REAL, ALLOCATABLE :: dummy(:)
#endif
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering check_err."
#endif
      if (iret .ne. NF_NOERR) then
         write(ScreenUnit,*) "ERROR: ",nf_strerror(iret)
         write(ScreenUnit,*) "ERROR: ADCIRC execution terminated."
         write(16,*) "ERROR: ",nf_strerror(iret)
         write(16,*) "ERROR: ADCIRC execution terminated."
#ifdef CMPI
         call msg_fini()
#endif
#ifdef NETCDF_DEBUG
         ! intentionally create a segmentation fault so that we can get
         ! a stack trace to determine the line number of the netcdf call
         ! that went bad ... this assumes that the code was compiled with 
         ! debugging symbols, bounds checking, and stack trace turned on. 
         dummy(1) = 99.9d0
#endif
         stop
      endif
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from check_err."
#endif
C-----------------------------------------------------------------------
      end subroutine check_err
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   S E T   A D C I R C   P A R A M E T E R S      
C-----------------------------------------------------------------------
C     jgf49.31 Called by adcprep to populate the global.F data normally
C     collected by read_input.F in adcirc. This is somewhat 
C     convoluted, since the data are collected by adcprep in read14()
C     and read15() and used to populate in pre_global.F, and then we
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
     &   base_date_p, STATNAME_p, STATNAMEV_p, STATNAMEM_p, NE_G_p, 
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
     &   nwp_p, tau0_p, cf_p, eslm_p, nstae_p, nstav_p, nstam_p, neta_p)
      USE GLOBAL, ONLY : 
     &   base_date, STATNAME, STATNAMEV, STATNAMEM, NE_G, NBOU,
     &   NVEL, NOPE, NP_G, SLAM0, SFEA0, NBVV, NVDLL, NBDV, NVELL, X, Y,
     &   IBTYPE, IBTYPEE,
     &   SLAM, SFEA, SNAMLEN, RAD2DEG, NODECODE, G, FileFmtRev,
     &   FileFmtMinor, FileFmtMajor, im, iestp, nscoue, ivstp, nscouv,
     &   icstp, nscouc, ipstp, iwstp, nscoum, igep, nscouge, igvp,
     &   nscougv, igcp, nscougc, igpp, igwp, nscougw, NM,
     &   DP, RUNDES, AGRID, title, institution, source, history,
     &   references, comments, host, convention, contact, dtdp, ihot,
     &   ics, nolifa, nolica, nolicat, ncor, ntip, nws, nramp, statim,
     &   reftim, rnday, dramp, a00, b00, c00, h0, cori, ntif, nbfr,
     &   myProc, screenUnit, neta
      USE NodalAttributes, ONLY : nolibf, nwp, tau0, cf, eslm
      IMPLICIT NONE
C
C     Declare the argument variables coming in from adcprep.
      CHARACTER(80), intent(in) :: base_date_p
      CHARACTER(50), intent(in) :: STATNAME_p(:)
      CHARACTER(50), intent(in) :: STATNAMEV_p(:)
      CHARACTER(50), intent(in) :: STATNAMEM_p(:)
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
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Entering setADCIRCParameters."
#endif
C
      base_date = base_date_p
      ALLOCATE(STATNAME(NSTAE_p))
      ALLOCATE(STATNAMEV(NSTAV_p))
      ALLOCATE(STATNAMEM(NSTAM_p))
      STATNAME = STATNAME_p
      STATNAMEV = STATNAMEV_p
      STATNAMEM = STATNAMEM_p
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
C
#ifdef NETCDF_DEBUG
      write(*,*) "DEBUG: Returning from setADCIRCParameters."
#endif
C-----------------------------------------------------------------------
      END SUBROUTINE setADCIRCParameters
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      END MODULE NETCDFIO
C-----------------------------------------------------------------------
