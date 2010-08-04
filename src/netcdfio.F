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
      ! 
      !=================================================================
      MODULE NETCDFIO
      USE SIZES, ONLY : SZ
      IMPLICIT NONE
      include 'netcdf.inc'
C
      DOUBLE PRECISION, PARAMETER ::  doubleval(1)=-99999.d0

      INTEGER ncmode !jgf: 0 for adcprep, 1 for serial and 2 for parallel?
      
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
         INTEGER :: nbdv_id       ! nodes on elev spec boundary seg
         INTEGER :: nbvv_id       ! nodes on normal flow boundary seg
         INTEGER :: nvdll_id      ! num nodes on elev boundary seg
         INTEGER :: nbounc_id
         INTEGER :: ibtypenc_id   ! boundary type
         INTEGER :: nvellnc_id    ! nodes on norm flow spec boundary seg
         INTEGER :: nvdllnc_id
         INTEGER :: nopenc_id  
!        Dimension ids
         INTEGER :: num_nodes_dim_id
         INTEGER :: num_elems_dim_id
         INTEGER :: nface_dim_id
         INTEGER :: nope_dim_id   ! num elev spec boundary forcing segs
         INTEGER :: neta_dim_id   ! total number of elevation specified boundary nodes
         INTEGER :: nbou_dim_id  ! number of normal flow specified boundary segment
         INTEGER :: nvel_dim_id  ! total number of normal flow specified boundary nodes
         INTEGER :: nvelnc_dim_id
         INTEGER :: netanc_dim_id
!        Dimension lengths
         INTEGER  num_nodes
         INTEGER  num_elems
         INTEGER  nface_len
         INTEGER  nope_len
         INTEGER  neta_len
         INTEGER  nbou_len
         INTEGER  nvel_len
         INTEGER  nvelnc_len
         INTEGER  netanc_len
!        Rank (number of dimensions) for each variable
         INTEGER :: depth_rank = 1
         INTEGER :: ele_rank = 2
         INTEGER :: x_rank = 1
         INTEGER :: y_rank = 1 
         INTEGER :: nbounc_rank = 1
         INTEGER :: nvdll_rank = 1
         INTEGER :: ibtypenc_rank = 1
         INTEGER :: nvellnc_rank = 1 
         INTEGER :: nvdllnc_rank = 1
         INTEGER :: nopenc_rank = 1
         INTEGER :: nbdv_rank = 2
         INTEGER :: nbvv_rank = 2
!        Variable shapes
         INTEGER  x_dims(1)
         INTEGER  y_dims(1)
         INTEGER  depth_dims(1)
         INTEGER  ele_dims(2)
         INTEGER  nvdll_dims(1)
         INTEGER  nbounc_dims(1)
         INTEGER  ibtypenc_dims(1)
         INTEGER  nvellnc_dims(1)
         INTEGER  nvdllnc_dims(1)
         INTEGER  nopenc_dims(1)
         INTEGER  nbdv_dims(2)
         INTEGER  nbvv_dims(2)
         REAL(SZ), ALLOCATABLE :: xnc(:) ! x coordinate or longitude
         REAL(SZ), ALLOCATABLE :: ync(:) ! y coordinate or latitude  
         INTEGER, ALLOCATABLE ::  nbvvnc(:,:) ! boundary array
         INTEGER, ALLOCATABLE ::  nbounc(:)   ! boundary array
         INTEGER, ALLOCATABLE ::  nvellnc(:)  ! boundary array
         INTEGER, ALLOCATABLE ::  nopenc(:)   ! boundary array
         INTEGER, ALLOCATABLE ::  nvdllnc(:)  ! boundary array
         INTEGER, ALLOCATABLE ::  nodecodenc(:) ! boundary array
         INTEGER, ALLOCATABLE ::  noffnc(:)     ! boundary array
         INTEGER, ALLOCATABLE ::  ibtypenc(:)     ! boundary array
         INTEGER, ALLOCATABLE ::  nm(:,:) 
         INTEGER, ALLOCATABLE ::  element(:,:)
         INTEGER nvelnc
         INTEGER netanc
      END TYPE meshStructure
 
      TYPE, PRIVATE :: fileData
         INTEGER record_counter
         LOGICAL createFile ! .true. if a new netCDF file must be created p
         CHARACTER*5 FILEFORT
         CHARACTER*10 FILENAME
         CHARACTER*3 FILEEXT
         LOGICAL fileFound  ! .true. if the netCDF file is present
      END Type fileData
 
      TYPE, PRIVATE :: stationData
         PRIVATE
         INTEGER lun               ! ADCIRC's logical unit number
         INTEGER ncid              ! the id of its netcdf file
         INTEGER num_stations      ! total number of stations 
         INTEGER, ALLOCATABLE :: name_lengths(:) ! lengths of station names 
         CHARACTER(50), POINTER :: statnames(:,:)
         REAL(SZ), ALLOCATABLE :: station_data(:,:) ! for scalar data 
         REAL(SZ), ALLOCATABLE :: u_station_data(:,:) ! for vector component 1
         REAL(SZ), ALLOCATABLE :: v_station_data(:,:) ! for vector component 2 
         REAL(SZ), ALLOCATABLE :: sta_x(:) ! x coordinate or longitude
         REAL(SZ), ALLOCATABLE :: sta_y(:) ! y coordinate or latitude  
         INTEGER :: timenc_len       ! number of time records  
         INTEGER :: timenc_dim_id
         INTEGER :: time_dim_id
         INTEGER :: timenc_rank = 1
         INTEGER :: timenc_id
         REAL(8), ALLOCATABLE :: timenc(:)
         INTEGER :: slen_dim_id
         INTEGER :: scode_dim_id
         INTEGER :: num_sta_dim_id
         INTEGER :: station_rank = 2
         INTEGER :: station_id
         INTEGER :: scode_rank = 2
         INTEGER :: zeta_id
         INTEGER :: zeta_rank = 2
         INTEGER :: u_id
         INTEGER :: u_rank = 2
         INTEGER :: v_id
         INTEGER :: v_rank = 2
         INTEGER :: timenc_dims(1)
         INTEGER :: station_dims(2)
         INTEGER :: scode_dims(2)
         INTEGER :: zeta_dims(2)
         INTEGER :: u_dims(2)
         INTEGER :: v_dims(2)

         type(fileData) :: myFile
      END TYPE stationData

      type(stationData), save :: elevSta ! for elevation stations (fort.61)
      type(stationData), save :: prSta   ! for pressure stations (fort.71)      
      type(stationData), save :: velSta  ! for velocity stations (fort.62)
      type(stationData), save :: wVelSta ! for wind velocity stations (fort.72)

      TYPE, PRIVATE :: nodalData
         PRIVATE
         INTEGER lun               ! ADCIRC's logical unit number
         INTEGER ncid              ! the id of its netcdf file
         INTEGER num_nodes         ! total number of stations 
         INTEGER num_nodes_dim_id
         REAL(SZ) :: initial_value ! array will be initialized to this
         REAL(SZ), ALLOCATABLE :: nodal_data(:,:)   ! for scalar data 
         REAL(SZ), ALLOCATABLE :: u_nodal_data(:,:) ! for 1st vector component
         REAL(SZ), ALLOCATABLE :: v_nodal_data(:,:) ! for 2nd vector component
         INTEGER :: timenc_len       ! number of time records  
         INTEGER :: timenc_dim_id
         INTEGER :: time_dim_id
         INTEGER :: timenc_rank = 1
         INTEGER :: timenc_id
         REAL(8), ALLOCATABLE :: timenc(:)
         INTEGER :: zeta_id
         INTEGER :: zeta_rank = 2
         INTEGER :: u_id
         INTEGER :: u_rank = 2
         INTEGER :: v_id
         INTEGER :: v_rank = 2
         INTEGER :: timenc_dims(1)
         INTEGER :: zeta_dims(2)
         INTEGER :: u_dims(2)
         INTEGER :: v_dims(2)
         type(fileData) :: myFile
         type(meshStructure) :: myMesh
      END TYPE nodalData

      type(nodalData), save :: elev       ! for fort.63
      type(nodalData), save :: pr         ! for fort.73
      type(nodalData), save :: currentVel ! for fort.64
      type(nodalData), save :: windVel    ! for fort.74

      TYPE, PRIVATE :: hotstartData
         INTEGER lun               ! ADCIRC's logical unit number
         INTEGER ncid              ! the id of its netcdf file
         INTEGER num_nodes         ! total number of stations 
         INTEGER num_nodes_dim_id
         INTEGER num_elems         ! total number of stations 
         INTEGER num_elems_dim_id
         INTEGER u_id
         INTEGER v_id
         INTEGER zeta1_id
         INTEGER zeta2_id
         INTEGER zetad_id
         INTEGER nodecodenc_id
         INTEGER noffnc_id
         INTEGER imhs_id
         INTEGER iths_id
         INTEGER iestp_id
         INTEGER nscoue_id
         INTEGER ivstp_id
         INTEGER nscouv_id
         INTEGER icstp_id
         INTEGER nscouc_id
         INTEGER ipstp_id
         INTEGER iwstp_id
         INTEGER nscoum_id
         INTEGER igep_id
         INTEGER nscouge_id
         INTEGER igvp_id
         INTEGER nscougv_id
         INTEGER igcp_id
         INTEGER nscougc_id
         INTEGER igpp_id
         INTEGER igwp_id
         INTEGER nscougw_id
         INTEGER :: u_rank = 2
         INTEGER :: v_rank = 2
         INTEGER :: zeta1_rank = 2
         INTEGER :: zeta2_rank = 2
         INTEGER :: zetad_rank = 2
         INTEGER :: nodecodenc_rank = 1
         INTEGER :: noffnc_rank = 1
         INTEGER :: imhs_rank = 0
         INTEGER :: iths_rank = 0 
         INTEGER :: IESTP_rank = 0
         INTEGER :: NSCOUE_rank = 0
         INTEGER :: IVSTP_rank = 0
         INTEGER :: NSCOUV_rank = 0
         INTEGER :: ICSTP_rank = 0
         INTEGER :: NSCOUC_rank = 0
         INTEGER :: IPSTP_rank = 0
         INTEGER :: IWSTP_rank = 0
         INTEGER :: NSCOUM_rank = 0
         INTEGER :: IGEP_rank = 0
         INTEGER :: NSCOUGE_rank = 0
         INTEGER :: IGVP_rank = 0
         INTEGER :: NSCOUGV_rank = 0
         INTEGER :: IGCP_rank = 0
         INTEGER :: NSCOUGC_rank = 0
         INTEGER :: IGPP_rank = 0
         INTEGER :: IGWP_rank = 0
         INTEGER :: NSCOUGW_rank = 0
         INTEGER  u_dims(2)
         INTEGER  v_dims(2)
         INTEGER  zeta1_dims(2)
         INTEGER  zeta2_dims(2)
         INTEGER  zetad_dims(2)
         INTEGER  nodecodenc_dims(1)
         INTEGER  noffnc_dims(1)
         REAL(SZ), ALLOCATABLE :: u(:,:)  
         REAL(SZ), ALLOCATABLE :: v(:,:)   
         REAL(SZ), ALLOCATABLE :: zeta1(:,:)  
         REAL(SZ), ALLOCATABLE :: zeta2(:,:)  
         REAL(SZ), ALLOCATABLE :: zetad(:,:)  
         INTEGER, ALLOCATABLE :: nodecodenc(:,:)  
         INTEGER, ALLOCATABLE :: noffnc(:,:)  
         INTEGER :: timenc_len       ! number of time records  
         INTEGER :: timenc_dim_id
         INTEGER :: timenc_rank = 1
         INTEGER :: timenc_id
         REAL(8), ALLOCATABLE :: timenc(:)
         INTEGER :: timenc_dims(1)
         type(fileData) :: myFile
         INTEGER :: FileFmtMajorFile
         INTEGER :: FileFmtMinorFile
         INTEGER :: FileFmtRevFile
         type(meshStructure) :: myMesh
      END TYPE hotstartData   
C
      type(hotstartData), save :: hs  ! for fort.67 or fort.68
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      CONTAINS
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E   S E T   N E T C D F   M O D E  
C-----------------------------------------------------------------------
C     jgf49.17.02 Allocates memory for NetCDF operations.
C-----------------------------------------------------------------------
      SUBROUTINE setNetCDFMode(mnproc)
      IMPLICIT NONE
      INTEGER, intent(in) :: mnproc
C
      ncmode = 1     ! reasonable default
      IF (mnproc.eq.0) THEN
         ncmode = 0  ! for adcprep?
      ENDIF
      IF (mnproc.eq.1) THEN
         ncmode = 1  ! for serial
      ENDIF
      IF (mnproc.gt.1) THEN
         ncmode = 2  ! for parallel
      ENDIF          
C-----------------------------------------------------------------------
      END SUBROUTINE setNetCDFMode
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   A L L O C   N E T C D F   A R R A Y S  
C-----------------------------------------------------------------------
C     jgf49.17.02 Allocates memory for NetCDF operations.
C-----------------------------------------------------------------------
      SUBROUTINE allocNetCDFArrays(lun, descript)
      USE GLOBAL, ONLY : OutputDataDescript_t
      IMPLICIT NONE
      INTEGER, intent(in) :: lun
      TYPE(OutputDataDescript_t) :: descript
C 
      ! don't allocate or initialize anything if this output type is not
      ! using netcdf
      IF (ABS(descript % specifier).ne.3) THEN
         RETURN  
      ENDIF
C
      SELECT CASE(lun)
    
         CASE(61)
            CALL initStation(elevSta, lun, descript)
         CASE(62) 
            CALL initStation(velSta, lun, descript)
         CASE(63)
            CALL initNodalData(elev, lun, descript)
         CASE(64) 
            CALL initNodalData(currentVel, lun, descript)
         CASE(67,68)
            ! do nothing, this is a binary hotstart file
         CASE(-67,-68)
            CALL initHotstart(lun, descript)
         CASE(71) 
            CALL initStation(prSta, lun, descript)
         CASE(72) 
            CALL initStation(wVelSta, lun, descript)
         CASE(73)
            CALL initNodalData(pr, lun, descript)
         CASE(74) 
            CALL initNodalData(windVel, lun, descript)
         CASE DEFAULT
             WRITE(*,*) 'ERROR: No netCDF for this filetype.' 

      END SELECT    
C
C-----------------------------------------------------------------------
      END SUBROUTINE allocNetCDFArrays
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E    I N I T  S T A T I O N   
C-----------------------------------------------------------------------
C     jgf49.17.02 Sets up netCDF variables and allocates memory for a 
C     station.
C-----------------------------------------------------------------------
      SUBROUTINE initStation(sta, lun, descript)
      USE GLOBAL, ONLY : SNAMLEN, OutputDataDescript_t, ICS, RAD2DEG,
     &                   IHOT, NWS, base_date, STATNAME, STATNAMEV,
     &                   STATNAMEM
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
C
!     date_string variables for time attribute
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer values(8)
C
      sta % lun = lun
      sta % num_stations = descript % num_fd_records
      
      select case(descript%num_items_per_record)

         case(1)
            ALLOCATE(sta%station_data(sta%num_stations, sta%timenc_len))
            DO n=1, sta%num_stations
               sta%station_data(n,sta%timenc_len) = descript%initial_value
            END DO
         case(2)
            ALLOCATE(sta%u_station_data(sta%num_stations, sta%timenc_len))
            ALLOCATE(sta%v_station_data(sta%num_stations, sta%timenc_len))
            DO n=1, sta % num_stations
               sta%u_station_data(n,sta%timenc_len) = descript%initial_value
               sta%v_station_data(n,sta%timenc_len) = descript%initial_value
            END DO
         case default
            WRITE(6,*) "ERROR: The station has ",
     &        descript%num_items_per_record," components; this is not ",
     &        "yet supported in ADCIRC's NetCDF."       

      end select
C
C     Set coordinates of each station, converting to degrees if we 
C     are in spherical coordinates. 
      ALLOCATE(sta%sta_x(sta % num_stations))
      ALLOCATE(sta%sta_y(sta % num_stations))
      IF(ICS.EQ.2) THEN
         DO I=1, sta % num_stations
            sta%sta_x(i)=descript%x_coord(i) * RAD2DEG
            sta%sta_y(i)=descript%y_coord(i) * RAD2DEG
         END DO
      ELSE
         DO I=1,sta%num_stations
            sta%sta_x(i)=descript%x_coord(i)
            sta%sta_y(i)=descript%y_coord(i)
         END DO
      ENDIF
C
C     Initialize netCDF output file, creating a new one if necessary.
      sta%myFile%record_counter = 0
      sta%myFile%createFile = .false.
      sta%myFile%fileFound = .false.
      sta%myFile%FILEFORT = 'fort.'
      sta%myFile%FILEEXT = '.nc'
      WRITE(sta%myFile%FILENAME,'(A5,I2,A3)') 
     &      sta%myFile%FILEFORT, lun, sta%myFile%FILEEXT
C
C     jgf49.17.02: Simplified the criteria for creating a new netCDF
C     file: coldstart, overwrite upon hotstart, or output
C     that does not already exist.
      INQUIRE(FILE=sta%myFile%FILENAME,EXIST=sta%myFile%fileFound)
      IF ((IHOT.EQ.0).OR.(descript%specifier.lt.0).OR.
     &   (sta%myFile%fileFound.eqv..false.)) THEN
         sta%myFile%createFile = .true.
      ENDIF
C     
C     RETURN if we don't need to create a file.
      IF (sta%myFile%createFile.eqv..false.) THEN
         RETURN
      ENDIF     

      iret = nf_create(sta%myFile%FILENAME, NF_CLOBBER, sta%ncid)
      CALL check_err(iret)
C
      CALL putMetaData(sta%ncid)
C       
!     Define time
      iret = nf_def_dim(sta%ncid, 'time',nf_unlimited,sta%timenc_dim_id)
      CALL check_err(iret)
      sta%timenc_dims(1) = sta%timenc_dim_id
      iret = nf_def_var(sta%ncid, 'time', NF_DOUBLE, sta%timenc_rank,
     &       sta%timenc_dims, sta%timenc_id)
      CALL check_err(iret)
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
C     Fill in labels and populate variables as appropriate for the
C     different types of data in the station files. The labels and 
C     units will also vary according to the coordinate system ADCIRC
C     is using (spherical or cartesian, according to the value of ICS)
C     as well as the units system (english or si according to the value of g).
      SELECT CASE(lun)

         CASE(61) !       F O R T . 6 1
   
            sta%zeta_dims(1) = sta%num_sta_dim_id
            sta%zeta_dims(2) = sta%timenc_dim_id
            iret = nf_def_var(sta%ncid, 'zeta', NF_DOUBLE,
     &             sta%zeta_rank, sta%zeta_dims, sta%zeta_id)
            CALL check_err(iret)
!           Define water surface elevation attributes
            iret = nf_put_att_text(sta%ncid, sta%zeta_id,
     &             'long_name',35, VAR61_LN)
            CALL check_err(iret)
            iret = nf_put_att_text(sta%ncid, sta%zeta_id,
     &             'standard_name', 23, VAR61_ST)
            CALL check_err(iret)
            CALL putUnitsAttribute(sta%ncid, sta%zeta_id, VAR_METS, VAR_FEET)
            iret = nf_put_att_DOUBLE(sta%ncid, sta%zeta_id, '_FillValue',
     &             NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            iret = nf_put_att_text(sta%ncid, sta%zeta_id,'positive', 2, VAR_UP)
            CALL check_err(iret)
            sta%statnames => STATNAME
 
         CASE(62) !       F O R T . 6 2 
 
            sta%u_dims(1) = sta%num_sta_dim_id
            sta%u_dims(2) = sta%timenc_dim_id
            iret = nf_def_var(sta%ncid, 'u-vel', NF_DOUBLE, 
     &             sta%u_rank, sta%u_dims, sta%u_id)
            CALL check_err(iret)
            sta%v_dims(1) = sta%num_sta_dim_id
            sta%v_dims(2) = sta%timenc_dim_id
            iret = nf_def_var(sta%ncid, 'v-vel', NF_DOUBLE,
     &             sta%v_rank, sta%v_dims, sta%v_id)
            CALL check_err(iret)
            if (ics.eq.2) then
               ! u
               iret = nf_put_att_text(sta%ncid, sta%u_id, 'long_name', 59, 
     &    'station water column vertically averaged east/west velocity')
               CALL check_err(iret)

               iret = nf_put_att_text(sta%ncid, sta%u_id,
     &                   'standard_name', 46, 
     &       'station_eastward_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf_put_att_text(sta%ncid, sta%u_id, 'positive',
     &                4, VAR_EAST)
               CALL check_err(iret)
               ! v
               iret = nf_put_att_text(sta%ncid,sta%v_id, 'long_name', 61, 
     &'station water column vertically averaged north/south velocity')
               CALL check_err(iret)
               iret = nf_put_att_text(sta%ncid,sta%v_id, 'standard_name', 47, 
     &         'station_northward_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf_put_att_text(sta%ncid, sta%v_id,
     &         'positive', 5, VAR_NORTH)
               CALL check_err(iret)
            else
               ! u
               iret = nf_put_att_text(sta%ncid, sta%u_id,
     &               'long_name', 64, 
     &            'station water column vertically averaged velocity in 
     &x-direction')
               CALL check_err(iret)

               iret = nf_put_att_text(sta%ncid, sta%u_id,
     &                   'standard_name', 39, 
     &       'station_x_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf_put_att_text(sta%ncid,sta%u_id,
     &                'positive', 5, VAR_RG)
               CALL check_err(iret)
               ! v
               iret = nf_put_att_text(sta%ncid, sta%v_id,
     &                'long_name', 64, 
     &            'station water column vertically averaged velocity in 
     &y-direction')
               CALL check_err(iret)
               iret = nf_put_att_text(sta%ncid, sta%v_id,
     &         'standard_name', 39, 
     &         'station_y_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf_put_att_text(sta%ncid, sta%v_id, 'positive',
     &         49, VAR_UPV)
               CALL check_err(iret)
            endif

            CALL putUnitsAttribute(sta%ncid, sta%u_id,VAR_MPS, VAR_FPS)
            iret = nf_put_att_double(sta%ncid, sta%u_id,
     &             '_FillValue', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            iret = nf_put_att_double(sta%ncid, sta%u_id,
     &             'dry_Value', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            CALL putUnitsAttribute(sta%ncid, sta%v_id,VAR_MPS,VAR_FPS)
            iret = nf_put_att_DOUBLE(sta%ncid, sta%v_id,
     &      '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
            CALL check_err(iret)
            iret = nf_put_att_DOUBLE(sta%ncid, sta%v_id,
     &      'dry_Value', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            sta%statnames => STATNAMEV

         CASE(71) !         F  O R T . 7 1 

            sta%zeta_dims(1) = sta%num_sta_dim_id
            sta%zeta_dims(2) = sta%timenc_dim_id
      
            iret = nf_def_var(sta%ncid, 'pressure', NF_DOUBLE,
     &      sta%zeta_rank, sta%zeta_dims, sta%zeta_id)
            CALL check_err(iret)
            iret = nf_put_att_text(sta%ncid, sta%zeta_id,
     &      'long_name',33, VAR71_LN)
            CALL check_err(iret)
            iret = nf_put_att_text(sta%ncid, sta%zeta_id,
     &      'standard_name', 25, VAR71_ST)
            CALL check_err(iret)
            iret = nf_put_att_text(sta%ncid, sta%zeta_id,
     &      'units', 15, VAR71_UN)
            CALL check_err(iret)
            iret = nf_put_att_DOUBLE(sta%ncid, sta%zeta_id,
     &      '_FillValue', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            iret = nf_put_att_text(sta%ncid, sta%zeta_id,
     &      'positive', 2, VAR_UP)
            CALL check_err(iret)
            sta%statnames => STATNAMEV
          
         CASE(72)  !      F O R T . 7 2 
      
            sta%u_dims(1) = sta%num_sta_dim_id
            sta%u_dims(2) = sta%timenc_dim_id
            iret = nf_def_var(sta%ncid, 'windx', NF_DOUBLE,
     &      sta%u_rank, sta%u_dims, sta%u_id)
            CALL check_err(iret)
            sta%v_dims(1) = sta%num_sta_dim_id
            sta%v_dims(2) = sta%timenc_dim_id
            iret = nf_def_var(sta%ncid, 'windy', NF_DOUBLE,
     &      sta%v_rank, sta%v_dims, sta%v_id)
            CALL check_err(iret)
C
            if (ics.eq.2) then
               select case(abs(nws)) 
                  case(1,2)
                     iret = nf_put_att_text(sta%ncid, sta%u_id,
     &                      'long_name', 23, 'station e/w wind stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid,sta%u_id,
     &                      'standard_name', 36,
     &                      'station_eastward_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%u_id,
     &                      'positive', 4, VAR_EAST)
                     iret = nf_put_att_text(sta%ncid, sta%v_id,
     &                      'long_name', 23, 'station n/s wind stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%v_id,
     &                     'standard_name', 37, 
     &                     'station_northward_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%v_id,
     &                      'positive', 5, VAR_NORTH)
                     CALL check_err(iret)
                  case default       
                     iret = nf_put_att_text(sta%ncid, sta%u_id,
     &                      'long_name', 25,'station e/w wind velocity')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%u_id,
     &                      'standard_name', 21, 
     &                      'station_eastward_wind')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%u_id,
     &                      'positive', 4, VAR_EAST)
                     iret = nf_put_att_text(sta%ncid, sta%v_id,
     &                      'long_name', 25, 
     &                      'station n/s wind velocity')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%v_id,
     &                      'standard_name',22,'station_northward_wind')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%v_id,
     &                      'positive', 5,VAR_NORTH)
               end select
            else
               select case(abs(nws))
                  case(1,2)
                     iret = nf_put_att_text(sta%ncid, sta%u_id,
     &                      'long_name', 33, 
     &                      'station wind stress in x-direction')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%u_id,
     &                      'standard_name', 29, 
     &                      'station_x_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%u_id,
     &                      'positive', 5, VAR_RG)
                     iret = nf_put_att_text(sta%ncid, sta%v_id,
     &                      'long_name', 33, 
     &                      'station wind stress in y-direction')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%v_id,
     &                      'standard_name', 29, 
     &                      'station_y_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%v_id,
     &                     'positive', 59, VAR_UPWS)
                  case default 
                     iret = nf_put_att_text(sta%ncid, sta%u_id,
     &                       'long_name', 36, 
     &                       'station wind velocity in x-direction')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%u_id,
     &                      'standard_name', 36, 'station_x_wind')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%u_id,
     &                      'positive', 5, VAR_RG)
                     iret = nf_put_att_text(sta%ncid, sta%v_id,
     &                      'long_name', 36, 
     &                      'station wind velocity in y-direction')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%v_id,
     &                      'standard_name', 36, 
     &                      'station_y_wind')
                     CALL check_err(iret)
                     iret = nf_put_att_text(sta%ncid, sta%v_id,
     &                      'positive', 61, VAR_UPW)
                     CALL check_err(iret)
               end select
            endif

            if (abs(nws).gt.2.and.abs(nws).lt.100) then
               call putUnitsAttribute(sta%ncid, sta%u_id, VAR_MPS, VAR_FPS)
               call putUnitsAttribute(sta%ncid, sta%v_id, VAR_MPS, VAR_FPS)
            else
               call putUnitsAttribute(sta%ncid, sta%u_id, VAR_MPS_2, VAR_FPS_2)
               call putUnitsAttribute(sta%ncid, sta%u_id, VAR_MPS_2, VAR_FPS_2)
            endif

            iret = nf_put_att_double(sta%ncid, sta%u_id, '_FillValue',
     &             NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            iret = nf_put_att_DOUBLE(sta%ncid, sta%v_id, '_FillValue',
     &             NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            sta%statnames => STATNAMEM
    
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
      iret = nf_put_att_text(sta%ncid,sta%timenc_id,
     &       'long_name',10,'model time')
      CALL check_err(iret)
      iret = nf_put_att_text(sta%ncid,sta%timenc_id,
     &       'standard_name', 4, 'time')
      CALL check_err(iret)
      iret = nf_put_att_text(sta%ncid, sta%timenc_id, 'units', 23, 
     .       'seconds since base_date')
      CALL check_err(iret)
      iret = nf_put_att_text(sta%ncid, sta%timenc_id, 'base_date',
     .       (LEN(TRIM(base_date(1:INDEX(base_date,"!")-1)))), base_date)
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
!
!     Leave define mode
      iret = nf_enddef(sta%ncid)
      CALL check_err(iret)
!   
!     Store station name
      do i=1,sta%num_stations
         start(1)=1
         start(2)=i
         kount(1)= LEN(TRIM(sta%statnames(i,SNAMLEN)))
         kount(2)=1
         iret = nf_put_vara_text(sta%ncid,sta%station_id,start,kount,
     &   sta%statnames(i,SNAMLEN) )
         CALL check_err(iret)
      end do
C-----------------------------------------------------------------------
      END SUBROUTINE initStation
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E    I N I T  N O D A L  D A T A  
C-----------------------------------------------------------------------
C     jgf49.17.02 Sets up netCDF variables and allocates memory for  
C     full domain data.
C-----------------------------------------------------------------------
      SUBROUTINE initNodalData(dat, lun, descript)
      USE GLOBAL, ONLY : OutputDataDescript_t, ICS, NWS
      IMPLICIT NONE
C
      type(nodalData), intent(inout) :: dat
      INTEGER, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: descript
C
      INTEGER n    ! node counter
      INTEGER iret ! success or failure of netcdf call
C 
      dat % lun = lun
      dat % num_nodes = descript % num_fd_records
      dat % initial_value = descript % initial_value
      ALLOCATE(dat % nodal_data(dat % num_nodes, dat % timenc_len))
C
C     initialization
      DO n=1, dat % num_nodes
         dat % nodal_data(n,dat%timenc_len) = descript % initial_value
      END DO
C
      CALL initNetCDFCoord(dat%myMesh)
      CALL defineMeshVariables(dat%ncid, dat%myMesh)

      SELECT CASE(lun)

         CASE(63)

            dat%zeta_dims(1) = dat%num_nodes_dim_id
            dat%zeta_dims(2) = dat%timenc_dim_id
            iret = nf_def_var(dat%ncid, 'zeta', NF_DOUBLE,
     &             dat%zeta_rank, dat%zeta_dims, dat%zeta_id)
            CALL check_err(iret)     
            iret = nf_put_att_text(dat%ncid, dat%zeta_id,
     &            'long_name',35, VAR63_LN)
            CALL check_err(iret)
            iret = nf_put_att_text(dat%ncid, dat%zeta_id,
     &             'standard_name', 23, VAR63_ST)
            CALL check_err(iret)
            CALL putUnitsAttribute(dat%ncid, dat%zeta_id,
     &           VAR_METS, VAR_FEET)  
            iret = nf_put_att_DOUBLE(dat%ncid, dat%zeta_id,
     &             '_FillValue', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            iret = nf_put_att_text(dat%ncid, dat%zeta_id,
     &             'positive', 2, VAR_UP)
            CALL check_err(iret)

         CASE(64)

            dat%u_dims(1) = dat%num_nodes_dim_id
            dat%u_dims(2) = dat%timenc_dim_id
            iret = nf_def_var(dat%ncid, 'u-vel', NF_DOUBLE,
     &             dat%u_rank, dat%u_dims, dat%u_id)
            CALL check_err(iret)

            dat%v_dims(1) = dat%num_nodes_dim_id
            dat%v_dims(2) = dat%timenc_dim_id
            iret = nf_def_var(dat%ncid, 'v-vel', NF_DOUBLE,
     &             dat%v_rank, dat%v_dims, dat%v_id)
            CALL check_err(iret)
            if (ics.eq.2) then
               iret = nf_put_att_text(dat%ncid, dat%u_id, 'long_name', 51, 
     &   'water column vertically averaged east/west velocity')
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%u_id,
     &                'standard_name', 38, 
     &                'eastward_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%u_id, 'positive', 4, VAR_EAST)
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%v_id, 'long_name', 53, 
     &    'water column vertically averaged north/south velocity')
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%v_id,'standard_name', 39, 
     &                      'northward_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%v_id, 'positive', 5, VAR_NORTH)
               CALL check_err(iret)
            else
               iret = nf_put_att_text(dat%ncid, dat%u_id, 'long_name', 56, 
     &  'water column vertically averaged velocity in x-direction')
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%u_id,'standard_name', 31, 
     &               'x_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%u_id, 'positive', 5, VAR_RG)
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%v_id, 'long_name', 56, 
     &       'water column vertically averaged velocity in y-direction')
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%v_id,'standard_name', 31, 
     &                      'y_water_velocity_depth_averaged')
               CALL check_err(iret)
               iret = nf_put_att_text(dat%ncid, dat%v_id, 'positive', 49, VAR_UPV)
               CALL check_err(iret)
            endif
            CALL putUnitsAttribute(dat%ncid, dat%u_id, VAR_MPS, VAR_FPS)
            CALL putUnitsAttribute(dat%ncid, dat%v_id, VAR_MPS, VAR_FPS)
            iret = nf_put_att_double(dat%ncid, dat%u_id, '_FillValue',
     &             NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
      
            iret = nf_put_att_double(dat%ncid, dat%u_id, 'dry_Value',
     &             NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            iret = nf_put_att_double(dat%ncid, dat%v_id, '_FillValue',
     &             NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
      
            iret = nf_put_att_double(dat%ncid, dat%v_id, 'dry_Value',
     &             NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)

         CASE(73)

            dat%zeta_dims(1) = dat%num_nodes_dim_id
            dat%zeta_dims(2) = dat%timenc_dim_id     
            iret = nf_def_var(dat%ncid,'pressure',NF_DOUBLE,
     &              dat%zeta_rank, dat%zeta_dims, dat%zeta_id) 
            CALL check_err(iret)
!           Define sea surface pressure attributes
            iret = nf_put_att_text(dat%ncid, dat%zeta_id,
     &             'long_name',25, VAR73_LN)
            CALL check_err(iret)
            iret = nf_put_att_text(dat%ncid, dat%zeta_id,
     &             'standard_name', 33, VAR73_ST)
            CALL check_err(iret)
            iret = nf_put_att_text(dat%ncid, dat%zeta_id, 'units', 15, VAR73_UN)
            CALL check_err(iret)
            iret = nf_put_att_DOUBLE(dat%ncid, dat%zeta_id,
     &            '_FillValue', NF_DOUBLE, 1, doubleval)
            CALL check_err(iret)
            iret = nf_put_att_text(dat%ncid, dat%zeta_id,'positive', 2, VAR_UP)
            CALL check_err(iret)

         CASE(74)

            dat%u_dims(1) = dat%num_nodes_dim_id
            dat%u_dims(2) = dat%timenc_dim_id
            iret = nf_def_var(dat%ncid, 'windx', NF_DOUBLE,
     &             dat%u_rank, dat%u_dims, dat%u_id)
            CALL check_err(iret)
            dat%v_dims(1) = dat%num_nodes_dim_id
            dat%v_dims(2) = dat%timenc_dim_id
            iret = nf_def_var(dat%ncid, 'windy', NF_DOUBLE,
     &             dat%v_rank, dat%v_dims, dat%v_id)
            CALL check_err(iret)
            if (ics.eq.2) then
               select case(abs(nws))
                  case(1,2)
                     iret = nf_put_att_text(dat%ncid, dat%u_id, 'long_name', 15, 
     &                       'e/w wind stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%u_id, 'standard_name', 28, 
     &                      'eastward_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%u_id, 'positive', 4, VAR_EAST)
                     iret = nf_put_att_text(dat%ncid, dat%v_id, 'long_name', 15, 
     &                       'n/s wind stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%v_id, 'standard_name', 29, 
     &                      'northward_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%v_id, 'positive', 5, VAR_NORTH)
                     CALL check_err(iret)

                  case default
                     iret = nf_put_att_text(dat%ncid, dat%u_id, 'long_name', 17, 
     &              'e/w wind velocity')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%u_id, 'standard_name', 13, 
     &                      'eastward_wind')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%u_id, 'positive', 4, VAR_EAST)
                     iret = nf_put_att_text(dat%ncid, dat%v_id, 'long_name', 17, 
     &       'n/s wind velocity')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%v_id, 'standard_name', 14, 
     &       'northward_wind')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%v_id, 'positive', 5, VAR_NORTH)
               end select
            else
               select case(abs(nws))
                  case(1,2)
                     iret = nf_put_att_text(dat%ncid, dat%u_id, 'long_name', 26, 
     &                       'wind stress in x-direction')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%u_id, 'standard_name', 21, 
     &                      'x_surface_wind_stress')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%u_id, 'positive', 5, VAR_RG)
                     CALL check_err(iret)                   
                     iret = nf_put_att_text(dat%ncid, dat%v_id, 'long_name', 26, 
     &                       'wind stress in y-direction')
              CALL check_err(iret)
              iret = nf_put_att_text(dat%ncid, dat%v_id, 'standard_name', 21, 
     &                      'y_surface_wind_stress')
             CALL check_err(iret)
             iret = nf_put_att_text(dat%ncid, dat%v_id, 'positive', 59, VAR_UPWS)
                  case default
                     iret = nf_put_att_text(dat%ncid, dat%u_id, 'long_name', 28, 
     &                       'wind velocity in x-direction')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%u_id, 'standard_name', 28, 
     &                      'x_wind')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%u_id, 'positive', 5, VAR_RG)

                     iret = nf_put_att_text(dat%ncid, dat%v_id, 'long_name', 28, 
     &                       'wind velocity in y-direction')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%v_id, 'standard_name', 28, 
     &                      'y_wind')
                     CALL check_err(iret)
                     iret = nf_put_att_text(dat%ncid, dat%v_id, 'positive', 61, VAR_UPW)
                     CALL check_err(iret)
               end select
            endif

            select case(abs(nws))
               case(1,2)
                  call putUnitsAttribute(dat%ncid, dat%u_id, VAR_MPS_2, VAR_FPS_2)
                  call putUnitsAttribute(dat%ncid, dat%v_id, VAR_MPS_2, VAR_FPS_2)
               case default
                  call putUnitsAttribute(dat%ncid, dat%u_id, VAR_MPS, VAR_FPS)
                  call putUnitsAttribute(dat%ncid, dat%v_id, VAR_MPS, VAR_FPS)
            end select
            iret = nf_put_att_double(dat%ncid, dat%u_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
            CALL check_err(iret)
            iret = nf_put_att_DOUBLE(dat%ncid, dat%v_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
         CASE DEFAULT

      END SELECT


C-----------------------------------------------------------------------
      END SUBROUTINE initNodalData
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     S U B R O U T I N E   I N I T   N E T C D F   C O O R D  
C-----------------------------------------------------------------------
C     jgf49.17.02 Initializes NetCDF coordinates. 
C-----------------------------------------------------------------------
      SUBROUTINE initNetCDFCoord(myMesh)
      USE GLOBAL, ONLY : RAD2DEG, NBOU, NVEL, NOPE, NP_G, SLAM0, SFEA0,
     &                   NBVV, NVDLL, NBDV, NVELL
      IMPLICIT NONE
      type(meshStructure), intent(inout) :: myMesh
      INTEGER :: i, j, k  ! array indices
      REAL(8) SLAM0DEG, SFEA0DEG
C
      ALLOCATE(myMesh%nbvvnc(nbou,nvel))
      ALLOCATE(myMesh%xnc(NP_G))
      ALLOCATE(myMesh%ync(NP_G))
      ALLOCATE(myMesh%nvellnc(nbou))
      ALLOCATE(myMesh%nvdllnc(nope))
      ALLOCATE(myMesh%ibtypenc(nbou))
      ALLOCATE(myMesh%nbounc(myMesh%nvelnc))
      ALLOCATE(myMesh%nopenc(myMesh%netanc))
      ALLOCATE(myMesh%element(myMesh%nface_len,myMesh%num_elems))
      ALLOCATE(myMesh%nm(myMesh%num_elems,myMesh%nface_len))

!     --------------------------------
!     Convert to degress for NETCDF I/O
!     --------------------------------
      SLAM0DEG=SLAM0
      SFEA0DEG=SFEA0
      IF (NCMODE.EQ.1) THEN
         SLAM0DEG=SLAM0*RAD2DEG
         SFEA0DEG=SFEA0*RAD2DEG
      ENDIF
C
      DO I=1,nbou
         DO J=1,nvel
            myMesh%nbvvnc(i,j)=nbvv(i,j)
         END DO
      END DO
!
!        elevation specified boundary forcing segments
      myMesh%netanc=0
      DO i=1,nope
         myMesh%nvdllnc(i)=nvdll(i)
         myMesh%netanc=myMesh%netanc+myMesh%nvdllnc(i)
      END DO
C
      k=0
      DO i=1,nope
         DO j=1,myMesh%nvdllnc(i)
            k=k+1
            myMesh%nopenc(k)=nbdv(i,j)
         END DO
      END DO
!
!        normal flow (discharge) specified boundary segments
      myMesh%nvelnc=0
      DO i=1,nbou
         myMesh%nvellnc(i)=nvell(i)
         myMesh%nvelnc=myMesh%nvelnc+myMesh%nvellnc(i)
      END DO
C
      k=0
      DO i=1,nbou
         DO j=1,myMesh%nvellnc(i)
            k=k+1
            myMesh%nbounc(k)=nbvv(i,j)
         END DO
      END DO
C
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
      USE SIZES
      USE GLOBAL
      USE NodalAttributes, ONLY :  nolibf, nwp, tau0, cf, eslm
C
      IMPLICIT NONE
C     args
      INTEGER, intent(in) :: lun ! logical unit number of file to write to
      REAL(8), intent(in) :: time ! seconds since cold start
      INTEGER, intent(in) :: it   ! number of time steps since cold start
      type(OutputDataDescript_t), intent(in) :: descript !describes output data
C     local vars
      INTEGER :: I,J,K           ! loop counter
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
         CASE(67,68) 
            CALL writeHotstart(lun, descript, time)
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
C-----------------------------------------------------------------------
      END SUBROUTINE writeOutArrayNetCDF
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E    W R I T E  S T A T I O N   
C-----------------------------------------------------------------------
C     jgf49.17.02 Writes data to station file.
C-----------------------------------------------------------------------
      SUBROUTINE writeStation(sta, lun, descript, time)
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
      sta%timenc(sta%timenc_len)=time
      iret = NF_OPEN(sta%myFile%FILENAME, NF_WRITE, sta%ncid)
      CALL check_err(iret)   
!   
!     Inquire variables
      iret=nf_inq_unlimdim(sta%ncid,recid)
      CALL check_err(iret)   
      iret=nf_inq_dimlen(sta%ncid,recid,sta%myFile%record_counter)
      CALL check_err(iret)   
!
!     Update Global Attributes
      CALL updateMetaData(sta%ncid)
      IF (ncmode.eq.0) go to 9999 ! i.e., if called by adcprep?

      kount(1)=sta%num_stations
      kount(2)=sta%timenc_len
      start(1)=1
      start(2)=sta%myFile%record_counter

      SELECT CASE(lun)

      CASE(61)
         iret=nf_inq_varid(sta%ncid, "zeta", sta%zeta_id)
         CALL check_err(iret)
         iret = nf_put_vara_DOUBLE(sta%ncid, sta%zeta_id, start, kount,
     &          sta%station_data)
         CALL check_err(iret)
      CASE(62)
         iret=nf_inq_varid(sta%ncid, "u-vel", sta%u_id)
         iret=nf_inq_varid(sta%ncid, "v-vel", sta%v_id)
         iret = nf_put_vara_DOUBLE(sta%ncid, sta%u_id, start, kount,
     &          sta%u_station_data)
         CALL check_err(iret)
         iret = nf_put_vara_DOUBLE(sta%ncid, sta%v_id, start, kount,
     &          sta%v_station_data)
         CALL check_err(iret)
      CASE(71)
         iret=nf_inq_varid(sta%ncid, "pressure", sta%zeta_id)
         CALL check_err(iret)
         iret = nf_put_vara_DOUBLE(sta%ncid, sta%zeta_id, start, kount,
     &          sta%station_data)
         CALL check_err(iret)
      CASE(72)  
           iret=nf_inq_varid(sta%ncid,"windx",sta%u_id)
           iret=nf_inq_varid(sta%ncid,"windy",sta%v_id)
           iret = nf_put_vara_DOUBLE(sta%ncid, sta%u_id, start, kount,
     &            sta%u_station_data)
           CALL check_err(iret)
           iret = nf_put_vara_DOUBLE(sta%ncid, sta%v_id, start, kount,
     &            sta%v_station_data)
           CALL check_err(iret)

      CASE DEFAULT
           WRITE(*,*) "ERROR: The unit number ",lun," is not ",
     &                "available in NetCDF."
          
      END SELECT
!    
!     Store time
      iret=nf_inq_varid(sta%ncid, "time", sta%timenc_id)
      starti(1)=sta%myFile%record_counter
      counti(1)=sta%timenc_len
      iret = nf_put_vara_DOUBLE(sta%ncid, sta%timenc_id, starti,counti,
     &       sta%timenc)
      CALL check_err(iret)
!
!     Close netCDF file
9999  iret = nf_close(sta%ncid)
      CALL check_err(iret)
C-----------------------------------------------------------------------
      end subroutine writeStation
C-----------------------------------------------------------------------     


C-----------------------------------------------------------------------
C     S U B R O U T I N E    W R I T E  N O D A L   D A T A   
C-----------------------------------------------------------------------
C     jgf49.17.02 Writes data to a full domain file.
C-----------------------------------------------------------------------
      SUBROUTINE writeNodalData(dat, lun, descript, time)
      USE GLOBAL, ONLY : SNAMLEN, OutputDataDescript_t
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
C
      dat%timenc(dat%timenc_len)=time
      iret = NF_OPEN(dat%myFile%FILENAME, NF_WRITE, dat%ncid)
      CALL check_err(iret)   
!   
!     Inquire variables
      iret=nf_inq_unlimdim(dat%ncid,recid)
      CALL check_err(iret)   
      iret=nf_inq_dimlen(dat%ncid,recid,dat%myFile%record_counter)
      CALL check_err(iret)   
!
!     Update Global Attributes
      CALL updateMetaData(dat%ncid)
      IF (ncmode.eq.0) go to 9999

      kount(1)=dat%num_nodes
      kount(2)=dat%timenc_len
      start(1)=1
      start(2)=dat%myFile%record_counter

      SELECT CASE(lun)

      CASE(63)
         iret=nf_inq_varid(dat%ncid,"zeta",dat%zeta_id)
         iret = nf_put_vara_DOUBLE(dat%ncid, dat%zeta_id, start, kount,
     &          dat%nodal_data)
         CALL check_err(iret)
      CASE(64)
         iret=nf_inq_varid(dat%ncid, "u-vel", dat%u_id)
         iret=nf_inq_varid(dat%ncid, "v-vel", dat%v_id)
         iret=nf_put_vara_DOUBLE(dat%ncid, dat%u_id, start, kount,
     &        dat%u_nodal_data)
         CALL check_err(iret)
         iret = nf_put_vara_DOUBLE(dat%ncid, dat%v_id, start, kount,
     &          dat%v_nodal_data)
         CALL check_err(iret)
      CASE(73)
         iret=nf_inq_varid(dat%ncid, "pressure", dat%zeta_id)
         CALL check_err(iret)
         iret=nf_put_vara_DOUBLE(dat%ncid, dat%zeta_id, start, kount,
     &        dat%nodal_data)
         CALL check_err(iret)
      CASE(74)
         iret=nf_inq_varid(dat%ncid,"windx",dat%u_id)
         iret=nf_inq_varid(dat%ncid,"windy",dat%v_id)
         iret=nf_put_vara_DOUBLE(dat%ncid, dat%u_id, start, kount,
     &        dat%u_nodal_data)
         CALL check_err(iret)
         iret = nf_put_vara_DOUBLE(dat%ncid, dat%v_id, start, kount,
     &          dat%v_nodal_data)
         CALL check_err(iret)
      CASE DEFAULT
 
      END SELECT
!    
!     Store time
      iret=nf_inq_varid(dat%ncid,"time",dat%timenc_id)
      starti(1)=dat%myFile%record_counter
      counti(1)=dat%timenc_len
      iret = nf_put_vara_DOUBLE(dat%ncid, dat%timenc_id, starti,counti,
     &       dat%timenc)
      CALL check_err(iret)
!
!     Close netCDF file
9999  iret = nf_close(dat%ncid)
      CALL check_err(iret)
C-----------------------------------------------------------------------
      end subroutine writeNodalData
C-----------------------------------------------------------------------     


C-----------------------------------------------------------------------
C     S U B R O U T I N E    I N I T  H O T S T A R T  
C-----------------------------------------------------------------------
C     jgf49.17.02 Sets up netCDF variables and allocates memory for  
C     hotstart file.
C-----------------------------------------------------------------------
      SUBROUTINE initHotstart(lun, descript)
      USE GLOBAL, ONLY : OutputDataDescript_t, ICS, G
      IMPLICIT NONE
C
      INTEGER, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: descript
C
      INTEGER iret ! success or failure of the netcdf call     
C
      hs%zeta1_dims(1) = hs%num_nodes_dim_id
      hs%zeta1_dims(2) = hs%timenc_dim_id
      iret = nf_def_var(hs%ncid, 'zeta1',NF_DOUBLE,
     &       hs%zeta1_rank, hs%zeta1_dims, hs%zeta1_id)
      CALL check_err(iret)

      hs%zeta2_dims(1) = hs%num_nodes_dim_id
      hs%zeta2_dims(2) = hs%timenc_dim_id
      iret = nf_def_var(hs%ncid, 'zeta2',NF_DOUBLE,
     &       hs%zeta2_rank, hs%zeta2_dims, hs%zeta2_id)
      CALL check_err(iret)

      hs%zetad_dims(1) = hs%num_nodes_dim_id
      hs%zetad_dims(2) = hs%timenc_dim_id
      iret = nf_def_var(hs%ncid, 'zeta1',NF_DOUBLE,
     &       hs%zetad_rank, hs%zetad_dims, hs%zetad_id)
      CALL check_err(iret)

      hs%u_dims(1) = hs%num_nodes_dim_id
      hs%u_dims(2) = hs%timenc_dim_id
      iret = nf_def_var(hs%ncid, 'u-vel', NF_DOUBLE, hs%u_rank, hs%u_dims, 
     &                  hs%u_id)
      CALL check_err(iret)

      hs%v_dims(1) = hs%num_nodes_dim_id
      hs%v_dims(2) = hs%timenc_dim_id
      iret = nf_def_var(hs%ncid, 'v-vel', NF_DOUBLE, hs%v_rank, hs%v_dims, 
     &                  hs%v_id)
      CALL check_err(iret)
!
!     Define hotstart parameters
      iret = nf_def_var(hs%ncid, 'nodecode', NF_INT, hs%nodecodenc_rank, 
     &                  hs%num_nodes_dim_id, hs%nodecodenc_id)
      CALL check_err(iret)
      iret = nf_def_var(hs%ncid, 'noff', NF_INT, hs%noffnc_rank, 
     &                 hs%num_elems_dim_id, hs%noffnc_id)
      CALL check_err(iret)
      iret = nf_def_var(hs%ncid, 'imhs', NF_INT, hs%imhs_rank, 0, 
     &                  hs%imhs_id)
      CALL check_err(iret)
      iret = nf_def_var(hs%ncid, 'iths', NF_INT, hs%iths_rank, 0, 
     &                  hs%iths_id)
      CALL check_err(iret)
      iret = nf_def_var(hs%ncid, 'iestp', NF_INT, hs%iestp_rank, 0, 
     &                  hs%iestp_id)
      CALL check_err(iret)
      iret = nf_def_var(hs%ncid, 'nscoue', NF_INT, hs%nscoue_rank, 0, 
     &                  hs%nscoue_id)
      CALL check_err(iret)

!         WRITE(6,*) "ivstp: "
      iret = nf_def_var(hs%ncid, 'ivstp', NF_INT, hs%ivstp_rank, 0, 
     &                  hs%ivstp_id)
      CALL check_err(iret)
      iret = nf_def_var(hs%ncid, 'nscouv', NF_INT, hs%nscouv_rank, 0, 
     &                  hs%nscouv_id)
      CALL check_err(iret)
      iret = nf_def_var(hs%ncid, 'icstp', NF_INT, hs%icstp_rank, 0, 
     &                  hs%icstp_id)
      CALL check_err(iret)
      iret = nf_def_var(hs%ncid, 'nscouc', NF_INT, hs%nscouc_rank, 0, 
     &                  hs%nscouc_id)
      CALL check_err(iret)

      iret = nf_def_var(hs%ncid, 'ipstp', NF_INT, hs%ipstp_rank, 0, 
     &                  hs%ipstp_id)
      CALL check_err(iret)
      iret = nf_def_var(hs%ncid, 'iwstp', NF_INT, hs%iwstp_rank, 0, 
     &                  hs%iwstp_id)
      CALL check_err(iret)
      iret = nf_def_var(hs%ncid, 'nscoum', NF_INT, hs%nscoum_rank, 0, 
     &                  hs%nscoum_id)
      CALL check_err(iret)
      iret = nf_def_var(hs%ncid, 'igep', NF_INT, hs%igep_rank, 0, 
     &                  hs%igep_id)
      CALL check_err(iret)
      
      iret = nf_def_var(hs%ncid, 'nscouge', NF_INT, hs%nscouge_rank, 0, 
     &                  hs%nscouge_id)
      CALL check_err(iret)
      iret = nf_def_var(hs%ncid, 'igvp', NF_INT, hs%igvp_rank, 0, 
     &                  hs%igvp_id)
      CALL check_err(iret)
      iret = nf_def_var(hs%ncid, 'nscougv', NF_INT, hs%nscougv_rank, 0, 
     &                  hs%nscougv_id)
      CALL check_err(iret)
      iret = nf_def_var(hs%ncid, 'igcp', NF_INT, hs%igcp_rank, 0, 
     &                  hs%igcp_id)
      CALL check_err(iret)
!         WRITE(6,*) "nscougc: "

      iret = nf_def_var(hs%ncid, 'nscougc', NF_INT, hs%nscougc_rank, 0, 
     &                  hs%nscougc_id)
      CALL check_err(iret)
      iret = nf_def_var(hs%ncid, 'igpp', NF_INT, hs%igpp_rank, 0, 
     &                  hs%igpp_id)
      CALL check_err(iret)
      iret = nf_def_var(hs%ncid, 'igwp', NF_INT, hs%igwp_rank, 0, 
     &                  hs%igwp_id)
      CALL check_err(iret)
      iret = nf_def_var(hs%ncid, 'nscougw', NF_INT, hs%nscougw_rank, 0, 
     &                  hs%nscougw_id)
      CALL check_err(iret)


       iret = nf_put_att_text(hs%ncid, hs%nodecodenc_id, 'long_name', 82, 
     &                       'node code indicating whether the node is
     & presently wet (active) or dry (inactive)')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%nodecodenc_id, 'standard_name', 18,
     &                'wet_dry_indicator')
      CALL check_err(iret)
!         WRITE(6,*) "Assign noffnc_id: "
      iret = nf_put_att_text(hs%ncid, hs%noffnc_id, 'long_name', 34, 
     &                       'Element based wetting/drying flag')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%noffnc_id, 'standard_name', 28,
     &                'element_wetting/drying_flag')
      CALL check_err(iret)

!         WRITE(6,*) "Assign imhs_id: "

      iret = nf_put_att_text(hs%ncid, hs%imhs_id, 'long_name', 18, 
     &                       'model type')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%imhs_id, 'standard_name', 10,
     &                'model_type')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%iths_id, 'long_name', 59, 
     &                       'model time step number since the beginning
     & of the model run')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%iths_id, 'standard_name', 15,
     &           'model_time_step')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%iestp_id, 'long_name', 176, 
     &                       'line number (for ASCII output) or record
     & number (for binary output) of the most recent entry in the
     & elevation time series at specified elevation recording stations
     & output file')
      CALL check_err(iret)
    
      iret = nf_put_att_text(hs%ncid, hs%iestp_id, 'standard_name', 54,
     &      'line/record_number_of_last_entry_in_elev_rec_stations')
      CALL check_err(iret)
      
      iret = nf_put_att_text(hs%ncid, hs%nscoue_id, 'long_name', 150, 
     &                       'time step counter to determine when the
     & next entry will be written to the elevation time series at
     & specified elevation recording Stations output file')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%nscoue_id, 'standard_name', 15,
     &           'time_step_counter_for_next_entry_elev_rec_stations)')
      CALL check_err(iret)
      
      iret = nf_put_att_text(hs%ncid, hs%ivstp_id, 'long_name', 190, 
     &                       ' line number (for ASCII output) or record 
     & number (for binary output) of the most recent entry in the
     & depth-averaged velocity time series at specified velocity 
     & recording stations output file')
      CALL check_err(iret)
      
      iret = nf_put_att_text(hs%ncid, hs%ivstp_id, 'standard_name', 48,
     &               'record_number_of_last_entry_in_vel_rec_stations')
      CALL check_err(iret)

       iret = nf_put_att_text(hs%ncid, hs%nscouv_id, 'long_name', 165, 
     &                       'time step counter to determine when the
     & next entry will be written to the Depth-averaged Velocity Time
     & series at specified velocity recording stations output file.')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%nscouv_id, 'standard_name', 50,
     &           'time_step_counter_for_next_entry_vel_rec_stations')
      CALL check_err(iret)

      iret = nf_put_att_text(hs%ncid, hs%icstp_id, 'long_name', 190, 
     &                       'line number (for ASCII output) or record
     & number (for binary output) of the most recent entry in the scalar
     & concentration time series at specified concentration recording
     & stations output file')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%icstp_id, 'standard_name', 49,
     &               'record_number_of_last_entry_in_conc_rec_stations')
      CALL check_err(iret)
     
      iret = nf_put_att_text(hs%ncid, hs%nscouc_id, 'long_name', 165, 
     &                       'time step counter to determine when the
     & next entry will be written to the scalar concentration Time 
     &series at specified concentration Recording Stations output file')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%nscouc_id, 'standard_name', 51,
     &           'time_step_counter_for_next_entry_conc_rec_stations')
      CALL check_err(iret)

!         WRITE(6,*) "Assign ipstp_id: "

      iret = nf_put_att_text(hs%ncid, hs%ipstp_id, 'long_name', 192, 
     &                       'line number (for ASCII output) or record
     & number (for binary output) of the most recent entry in the 
     &atmospheric pressure time series at specified meteorological 
     &recording stations')
      CALL check_err(iret)
     
      iret = nf_put_att_text(hs%ncid, hs%ipstp_id, 'standard_name', 57,
     &       'record_number_of_last_entry_of_atm_press_at_rec_stations')
      CALL check_err(iret)
     
!         WRITE(6,*) "Assign iwstp_id: "
      iret = nf_put_att_text(hs%ncid, hs%iwstp_id, 'long_name', 190, 
     &                       'line number (for ASCII output) or record
     & number (for binary output) of the most recent entry in the Wind 
     & Velocity time series at specified meteorological recording 
     &stations')
      CALL check_err(iret)

      iret = nf_put_att_text(hs%ncid, hs%iwstp_id, 'standard_name', 56,
     &       'record_number_of_last_entry_of_wind_vel_at_rec_stations')
      CALL check_err(iret)

      iret = nf_put_att_text(hs%ncid, hs%nscoum_id, 'long_name', 244, 
     &                       'time step counter to determine when the
     & next entry will be written to the atmospheric pressure time 
     &series at specified meteorological recording stations and wind 
     &velocity time series at specified meteorological recording 
     &stations output files')
      CALL check_err(iret)

!         WRITE(6,*) "Assign nscoum_id: "
      iret = nf_put_att_text(hs%ncid, hs%nscoum_id, 'standard_name', 60,
     &   'time_step_counter_of_atm_press_and_wind_vel_at_rec_stations')
      CALL check_err(iret)

      iret = nf_put_att_text(hs%ncid, hs%igep_id, 'long_name', 165, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &Elevation time series at All Nodes in the Model Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(hs%ncid, hs%igep_id, 'standard_name', 51,
     &       'record_number_of_last_entry_of_elev_at_model_nodes')
      CALL check_err(iret)

!         WRITE(6,*) "Assign nscouge_id: "
      iret = nf_put_att_text(hs%ncid, hs%nscouge_id, 'long_name', 139, 
     &                       'time step counter to determine when the
     & next entry will be written to the  Elevation time series at All 
     &Nodes in the Model Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(hs%ncid, hs%nscouge_id, 'standard_name', 41,
     &   'time_step_counter_of_elev_at_model_nodes')
      CALL check_err(iret)

      iret = nf_put_att_text(hs%ncid, hs%igvp_id, 'long_name', 179, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &Depth-averaged velocity time series at All Nodes in the Model 
     &Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(hs%ncid, hs%igvp_id, 'standard_name', 50,
     &       'record_number_of_last_entry_of_vel_at_model_nodes')
      CALL check_err(iret)

!         WRITE(6,*) "Assign nscougv_id: "
      iret = nf_put_att_text(hs%ncid, hs%nscougv_id, 'long_name', 153, 
     &                       'time step counter to determine when the 
     &next entry will be written to the Depth-averaged velocity time 
     &series at All Nodes in the Model Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(hs%ncid, hs%nscougv_id, 'standard_name', 40,
     &   'time_step_counter_of_vel_at_model_nodes')
      CALL check_err(iret)

      iret = nf_put_att_text(hs%ncid, hs%igcp_id, 'long_name', 176, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &scalar Concentration time series at All Nodes in the Model 
     &Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(hs%ncid, hs%igcp_id, 'standard_name', 51,
     &       'record_number_of_last_entry_of_conc_at_model_nodes')
      CALL check_err(iret)

!         WRITE(6,*) "Assign nscougc_id: "
      iret = nf_put_att_text(hs%ncid, hs%nscougc_id, 'long_name', 150, 
     &                       'time step counter to determine when the 
     &next entry will be written to the scalar Concentration time 
     &series at All Nodes in the Model Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(hs%ncid, hs%nscougc_id, 'standard_name', 41,
     &   'time_step_counter_of_conc_at_model_nodes')
      CALL check_err(iret)

      iret = nf_put_att_text(hs%ncid, hs%igpp_id, 'long_name', 176, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &atmospheric pressure time series at All Nodes in the Model 
     &Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(hs%ncid, hs%igpp_id, 'standard_name', 61,
     &       'record_number_of_last_entry_of_atm_press_at_model_nodes')
      CALL check_err(iret)

      iret = nf_put_att_text(hs%ncid, hs%igwp_id, 'long_name', 179, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &Wind Stress or velocity time series at All Nodes in the Model 
     &Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(hs%ncid, hs%igwp_id, 'standard_name', 55,
     &       'record_number_of_last_entry_of_wind_vel_at_model_nodes')
      CALL check_err(iret)

!         WRITE(6,*) "Assign nscougw_id: "
      iret = nf_put_att_text(hs%ncid, hs%nscougw_id, 'long_name', 222, 
     &                       'time step counter to determine when the 
     &next entry will be written to the atmospheric pressure time series 
     &at All Nodes in the Model Grid and Wind Stress or velocity time 
     &series at All Nodes in the Model Grid output files')
      CALL check_err(iret)

      iret = nf_put_att_text(hs%ncid, hs%nscougw_id, 'standard_name', 59,
     &   'time_step_counter_of_atm_press_and_wind_vel_at_model_nodes')
      CALL check_err(iret)


!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(hs%ncid, hs%u_id, 'long_name', 32, 
     &                       'vertically averaged e/w velocity')
       else
         iret = nf_put_att_text(hs%ncid, hs%u_id, 'long_name', 43, 
     &            'vertically averaged velocity in x-direction')
      endif
      CALL check_err(iret)

      iret = nf_put_att_text(hs%ncid, hs%u_id, 'standard_name', 10, 
     &                      'u_velocity')
      CALL check_err(iret)
!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
          iret = nf_put_att_text(hs%ncid, hs%u_id, 'units',5, VAR_MPS)
      else
          iret = nf_put_att_text(hs%ncid, hs%u_id, 'units',6, VAR_FPS)
      endif
      CALL check_err(iret)
!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(hs%ncid, hs%u_id, 'positive', 4, VAR_EAST)
      else
          iret = nf_put_att_text(hs%ncid, hs%u_id, 'positive', 5, VAR_RG)
      endif
      CALL check_err(iret)
      iret = nf_put_att_double(hs%ncid, hs%u_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
      
      iret = nf_put_att_double(hs%ncid, hs%u_id, 'dry_Value', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
!    -----------
!     Define v velocity
!    -----------
!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(hs%ncid, hs%v_id, 'long_name', 32, 
     &           'vertically averaged n/s velocity')
      else
          iret = nf_put_att_text(hs%ncid, hs%v_id, 'long_name', 43, 
     &           'vertically averaged velocity in y-direction')
      endif
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%v_id, 'standard_name', 10, 
     &                      'v_velocity')
      CALL check_err(iret)
!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
         iret = nf_put_att_text(hs%ncid, hs%v_id, 'units',5, VAR_MPS)
      else
         iret = nf_put_att_text(hs%ncid, hs%v_id, 'units',6, VAR_FPS)
      endif
      CALL check_err(iret)

!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(hs%ncid, hs%v_id, 'positive', 5, VAR_NORTH)
      else
          iret = nf_put_att_text(hs%ncid, hs%v_id, 'positive', 49, VAR_UPV)
      endif
      CALL check_err(iret)

      iret = nf_put_att_DOUBLE(hs%ncid, hs%v_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
!      doubleval(1) = -99999.
      iret = nf_put_att_DOUBLE(hs%ncid, hs%v_id, 'dry_Value', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)

!    --------------------------------------------
!     Define water surface elevation attributes
!    --------------------------------------------
      
      iret = nf_put_att_text(hs%ncid, hs%zeta1_id,'long_name',45, 
     &      'water surface elevation at previous time step')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%zeta1_id, 'standard_name', 45, 
     &      'water_surface_elevation_at_previous_time step')
      CALL check_err(iret)
      
      if (G.LT.11.D0) then
          iret = nf_put_att_text(hs%ncid, hs%zeta1_id, 'units', 1, VAR_METS)
          CALL check_err(iret)
      else
         iret = nf_put_att_text(hs%ncid, hs%zeta1_id, 'units', 2, VAR_FEET)
         CALL check_err(iret)
      endif      
      
      iret = nf_put_att_DOUBLE(hs%ncid, hs%zeta1_id, '_FillValue', NF_DOUBLE, 
     &                      1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%zeta1_id,'positive', 2, VAR_UP)
      CALL check_err(iret)
!    --------------------------------------------
!     Define water surface elevation attributes
!    --------------------------------------------
      
      iret = nf_put_att_text(hs%ncid, hs%zeta2_id,'long_name',44, 
     &      'water surface elevation at current time step')
                       
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%zeta2_id, 'standard_name', 44, 
     &      'water_surface_elevation_at_current_time_step')
      CALL check_err(iret)
      
      if (G.LT.11.D0) then
          iret = nf_put_att_text(hs%ncid, hs%zeta2_id, 'units', 1, VAR_METS)
          CALL check_err(iret)
      else
         iret = nf_put_att_text(hs%ncid, hs%zeta2_id, 'units', 2, VAR_FEET)
         CALL check_err(iret)
      endif      
      
      iret = nf_put_att_DOUBLE(hs%ncid, hs%zeta2_id, '_FillValue', NF_DOUBLE, 
     &                      1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%zeta2_id,'positive', 2, VAR_UP)
      CALL check_err(iret)
!    --------------------------------------------
!     Define water surface elevation attributes
!    --------------------------------------------

      iret = nf_put_att_text(hs%ncid, hs%zetad_id,'long_name',42,       
     &      'water elevation at flux specified boundary')
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%zetad_id, 'standard_name', 42, 
     &      'water_elevation_at_flux_specified_boundary')
      CALL check_err(iret)
      
      if (G.LT.11.D0) then
          iret = nf_put_att_text(hs%ncid, hs%zetad_id, 'units', 1, VAR_METS)
          CALL check_err(iret)
      else
         iret = nf_put_att_text(hs%ncid, hs%zetad_id, 'units', 2, VAR_FEET)
         CALL check_err(iret)
      endif      
      
      iret = nf_put_att_DOUBLE(hs%ncid, hs%zetad_id, '_FillValue', NF_DOUBLE, 
     &                      1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(hs%ncid, hs%zetad_id,'positive', 2, VAR_UP)
      CALL check_err(iret)

C-----------------------------------------------------------------------
      END SUBROUTINE initHotstart
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E    W R I T E  H O T S T A R T   
C-----------------------------------------------------------------------
C     jgf49.17.02 Writes data to hotstart file.
C-----------------------------------------------------------------------
      SUBROUTINE writeHotstart(lun, descript, time)
      USE GLOBAL
      IMPLICIT NONE
C
      INTEGER, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: descript
      REAL(8), intent(in) :: time
C
      INTEGER i,j
      INTEGER counti(1), starti(1), n
      INTEGER kount(2), start(2)
      INTEGER iret ! success or failure of the netcdf call
C
      hs%timenc(hs%timenc_len)=time
      iret=nf_inq_varid(hs%ncid,"zeta1",hs%zeta1_id)
      iret=nf_inq_varid(hs%ncid,"zeta2",hs%zeta2_id)
      iret=nf_inq_varid(hs%ncid,"zetad",hs%zetad_id)
      iret=nf_inq_varid(hs%ncid,"u-vel",hs%u_id)
      iret=nf_inq_varid(hs%ncid,"v-vel",hs%v_id)
      iret=nf_inq_varid(hs%ncid,"nodecode",hs%nodecodenc_id)
      iret=nf_inq_varid(hs%ncid,"noff",hs%noffnc_id)
      iret=nf_inq_varid(hs%ncid,"imhs",hs%imhs_id)
      iret=nf_inq_varid(hs%ncid,"iths",hs%iths_id)
      iret=nf_inq_varid(hs%ncid,"iestp",hs%iestp_id)
      iret=nf_inq_varid(hs%ncid,"nscoue",hs%nscoue_id)
      iret=nf_inq_varid(hs%ncid,"ivstp",hs%ivstp_id)
      iret=nf_inq_varid(hs%ncid,"nscouv",hs%nscouv_id)
      iret=nf_inq_varid(hs%ncid,"icstp",hs%icstp_id)
      iret=nf_inq_varid(hs%ncid,"nscouc",hs%nscouc_id)
      iret=nf_inq_varid(hs%ncid,"ipstp",hs%ipstp_id)
      iret=nf_inq_varid(hs%ncid,"iwstp",hs%iwstp_id)
      iret=nf_inq_varid(hs%ncid,"nscoum",hs%nscoum_id)
      iret=nf_inq_varid(hs%ncid,"igep",hs%igep_id)
      iret=nf_inq_varid(hs%ncid,"nscouge",hs%nscouge_id)
      iret=nf_inq_varid(hs%ncid,"igvp",hs%igvp_id)
      iret=nf_inq_varid(hs%ncid,"nscougv",hs%nscougv_id)
      iret=nf_inq_varid(hs%ncid,"igcp",hs%igcp_id)
      iret=nf_inq_varid(hs%ncid,"nscougc",hs%nscougc_id)
      iret=nf_inq_varid(hs%ncid,"igpp",hs%igpp_id)
      iret=nf_inq_varid(hs%ncid,"igwp",hs%igwp_id)
      iret=nf_inq_varid(hs%ncid,"nscougw",hs%nscougw_id)

      kount(1)=hs%num_nodes
      kount(2)=hs%timenc_len
      start(1)=1
      start(2)=hs%myFile%record_counter

      iret = nf_put_vara_DOUBLE(hs%ncid, hs%zeta1_id, start, kount, hs%zeta1)
      CALL check_err(iret)
      iret = nf_put_vara_DOUBLE(hs%ncid, hs%zeta2_id, start, kount, hs%zeta2)
      CALL check_err(iret)
      iret = nf_put_vara_DOUBLE(hs%ncid, hs%zetad_id, start, kount, hs%zetad)
      CALL check_err(iret)
      iret = nf_put_vara_DOUBLE(hs%ncid, hs%u_id, start, kount, hs%u)
      CALL check_err(iret)
      iret = nf_put_vara_DOUBLE(hs%ncid, hs%v_id, start, kount, hs%v)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%nodecodenc_id, hs%nodecodenc)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%noffnc_id, hs%noffnc)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%imhs_id, imhs)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%iths_id, iths)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%iestp_id, iestp)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%nscoue_id, nscoue)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%ivstp_id, ivstp)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%nscouv_id, nscouv)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%icstp_id, icstp)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%nscouc_id, nscouc)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%ipstp_id, ipstp)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%iwstp_id, iwstp)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%nscoum_id, nscoum)
      CALL check_err(iret)
!               WRITE(6,*) "=================Store igep", igep
      iret = nf_put_var_int(hs%ncid, hs%igep_id, igep)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%nscouge_id, nscouge)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%igvp_id, igvp)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%nscougv_id, nscougv)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%igcp_id, igcp)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%nscougc_id, nscougc)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%igpp_id, igpp)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%igwp_id, igwp)
      CALL check_err(iret)
      iret = nf_put_var_int(hs%ncid, hs%nscougw_id, nscougw)
      CALL check_err(iret)
C-----------------------------------------------------------------------
      END SUBROUTINE writeHotstart
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E    R E A D  H O T S T A R T   
C-----------------------------------------------------------------------
C     jgf49.17.02 Reads data from the hotstart file.
C-----------------------------------------------------------------------
      SUBROUTINE readHotstart(lun)
      USE GLOBAL, ONLY : OutputDataDescript_t, FileFmtRev, FileFmtMinor,
     &                   FileFmtMajor, imhs, iths, iestp, nscoue, ivstp,
     &                   nscouv, icstp, nscouc, ipstp, iwstp, nscoum,
     &                   igep, nscouge, igvp, nscougv, igcp, nscougc,
     &                   igpp, igwp, nscougw
      IMPLICIT NONE
C
      INTEGER, intent(in) :: lun
C      
      INTEGER i,j
      INTEGER counti(1), starti(1), n
      INTEGER kount(2), start(2)
      INTEGER iret ! success or failure of the netcdf call
C
      iret=nf_inq_dimid(hs%ncid, "node", hs%num_nodes_dim_id)
      CALL check_err(iret)   
      iret=nf_inq_dimlen(hs%ncid, hs%num_nodes_dim_id, hs%num_nodes)
      CALL check_err(iret)   
      iret=nf_inq_dimid(hs%ncid, "nele", hs%num_elems_dim_id)
      CALL check_err(iret)   
      iret=nf_inq_dimlen(hs%ncid, hs%num_elems_dim_id, hs%num_elems)
      CALL check_err(iret)   
C       IF(hs%nlen.NE.hs%node_len.OR.hs%nelen.NE.hs%elem_len) THEN
C         WRITE(6,887) hs%node_len,hs%elem_len,hs%nlen,hs%nelen
C887      FORMAT("===================================================",/,
C     &   "The grid input file you are reading is not",/, 
C     &   "consistent with the grid in the hotstart file.",/,
C     &   "Your grid file has the following parameters:",/,/,
C     &   "Number of nodes    = ",I10,/,
C     &   "Number of elements = ",I10,/,/,
C     &   "The hotstart file has the following parameters",/, 
C     &   "Number of nodes    = ",I10,/,
C     &   "Number of elements = ",I10,/,/,
C     &    "         ******* ADCIRC will stop! *******",/,
C     &    "====================================================")
C         STOP
C      ENDIF

      iret = nf_get_att_int(hs%ncid, NCGLOBAL, 'FileFmtMajor',  
     &                      hs%FileFmtMajorFile)
      CALL check_err(iret)
      iret = nf_get_att_int(hs%ncid, NCGLOBAL, 'FileFmtMinor',  
     &                      hs%FileFmtMinorFile)
      CALL check_err(iret)
      iret = nf_get_att_int(hs%ncid, NCGLOBAL, 'FileFmtRev', 
     &                      hs%FileFmtRevFile)
      CALL check_err(iret)
      
      IF(hs%FileFmtMajorFile.NE.FileFmtMajor.OR.
     &   hs%FileFmtMinorFile.NE.FileFmtMinor.OR.
     &   hs%FileFmtRevFile.NE.FileFmtRev)THEN
         WRITE(6,888) hs%FileFmtMajorFile,hs%FileFmtMinorFile,hs%FileFmtRevFile,
     & FileFmtMajor,FileFmtMinor,FileFmtRev
888      FORMAT("===================================================",/,
     &   "The input hotstart file you are reading is not",/, 
     &    "consistent with the version of the model you are running.",/,
     &    "Your file has the following file format parameters:",/,/,
     &    "FileFmtMajorFile = ",I3,/,
     &    "FileFmtMinorFile = ",I3,/,
     &    "FileFmtRevFile   = ",I3,/,/,
     &    "The current version of ADCIRC needs a hotstart file",/, 
     &    "with the following file format parameters:",/,/,
     &    "FileFmtMajor = ",I3,/,
     &    "FileFmtMinor = ",I3,/,
     &    "FileFmtRev   = ",I3,/,/,
     &    "         ******* ADCIRC will stop! *******",/,
     &    "====================================================")
         STOP
      ENDIF
        
      hs%myFile%record_counter = 1
      kount(1)=hs%num_nodes
      kount(2)=hs%timenc_len
      start(1)=1
      start(2)=hs%myFile%record_counter
!   
!     Store time
      starti(1)=hs%myFile%record_counter
      counti(1)=hs%timenc_len
      iret = nf_get_vara_DOUBLE(hs%ncid, hs%timenc_id, starti,counti,
     &       hs%timenc)
      CALL check_err(iret)
      iret = nf_get_vara_DOUBLE(hs%ncid, hs%zeta1_id, start, kount, hs%zeta1)
      CALL check_err(iret)
      iret = nf_get_vara_DOUBLE(hs%ncid, hs%zeta2_id, start, kount, hs%zeta2)
      CALL check_err(iret)
      iret = nf_get_vara_DOUBLE(hs%ncid, hs%zetad_id, start, kount, hs%zetad)
      CALL check_err(iret)
      iret = nf_get_vara_DOUBLE(hs%ncid, hs%u_id, start, kount, hs%u)
      CALL check_err(iret)
      iret = nf_get_vara_DOUBLE(hs%ncid, hs%v_id, start, kount, hs%v)
      CALL check_err(iret)
      iret = nf_get_vara_int(hs%ncid, hs%nodecodenc_id, start, kount, 
     &                       hs%nodecodenc)
      CALL check_err(iret)
      start(1)=1
      kount(1)=hs%num_elems
      iret = nf_get_vara_int(hs%ncid, hs%noffnc_id, start, kount, hs%noffnc)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, hs%imhs_id, imhs)
      CALL check_err(iret)
!      PRINT *,"NETCDFIO imhs = ", imhs
      iret = nf_get_var_int(hs%ncid, hs%iths_id, iths)
      CALL check_err(iret)
!      PRINT *,"NETCDFIO iths = ", iths
      iret = nf_get_var_int(hs%ncid, hs%iestp_id,iestp)
      CALL check_err(iret)
!      PRINT *,"NETCDFIO IESTP = ", iestp
      iret = nf_get_var_int(hs%ncid, hs%nscoue_id, nscoue)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, hs%ivstp_id, ivstp)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, hs%nscouv_id, nscouv)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, hs%icstp_id, icstp)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, hs%nscouc_id, nscouc)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, hs%ipstp_id, ipstp)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, hs%iwstp_id, iwstp)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, hs%nscoum_id,nscoum)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, hs%igep_id, igep)
      CALL check_err(iret)
!      PRINT *,"NETCDFIO igep = ", igep
      iret = nf_get_var_int(hs%ncid, hs%nscouge_id, nscouge)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, hs%igvp_id, igvp)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, hs%nscougv_id, nscougv)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, hs%igcp_id, igcp)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, hs%nscougc_id, nscougc)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, hs%igpp_id, igpp)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, hs%igwp_id, igwp)
      CALL check_err(iret)
      iret = nf_get_var_int(hs%ncid, hs%nscougw_id, nscougw)
      CALL check_err(iret)
C-----------------------------------------------------------------------
      END SUBROUTINE readHotstart
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   D E F I N E   M E S H   V A R I A B L E S 
C-----------------------------------------------------------------------
C     jgf49.17.02 Defines data that are common to all netCDF files.
C-----------------------------------------------------------------------
      subroutine defineMeshVariables(ncid, myMesh)
      USE GLOBAL, ONLY : RAD2DEG, SLAM0, SFEA0
      USE SIZES, ONLY : SZ
      USE NodalAttributes, ONLY : nwp, nolibf, tau0, cf, eslm
      USE VERSION
      IMPLICIT NONE
      INTEGER iret   ! Error status return
      type(meshStructure) :: myMesh
      INTEGER, intent(in) :: ncid   ! netCDF id
      INTEGER i,j
      REAL(8) SLAM0DEG, SFEA0DEG
      INTEGER counti(1), starti(1), n
      INTEGER kount(2), start(2)
      REAL    realval(1)
!     -----------------
!      date_string variables for time attribute
!     -----------------
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer values(8)
C
!     Switch order in array for NETCDF
      do i=1, myMesh%num_elems
         do j=1, myMesh%nface_len
            myMesh%element(j,i) = myMesh%nm(i,j)
         end do
      end do
!
!     Convert to degress for NETCDF I/O
      SLAM0DEG=SLAM0
      SFEA0DEG=SFEA0
      IF (NCMODE.EQ.1) THEN
         SLAM0DEG=SLAM0*RAD2DEG
         SFEA0DEG=SFEA0*RAD2DEG
      ENDIF
! 
!     Define dimensions
      iret = nf_def_dim(ncid,'node',myMesh%num_nodes,myMesh%num_nodes_dim_id)
      CALL check_err(iret)
      iret = nf_def_dim(ncid,'nele',myMesh%num_elems,myMesh%num_elems_dim_id)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'nvertex', 3, myMesh%nface_dim_id)
      call check_err(iret)
      if(myMesh%nope_len.ne.0) then
         iret = nf_def_dim(ncid,'nope',
     &             myMesh%nope_len,myMesh%nope_dim_id)
         call check_err(iret)
      endif
      if(myMesh%netanc_len.ne.0) then
         iret = nf_def_dim(ncid, 'neta', 
     &          myMesh%netanc_len, myMesh%neta_dim_id)
         call check_err(iret)
      endif
      if(myMesh%nbou_len.ne.0) then
         iret = nf_def_dim(ncid, 'nbou',
     &          myMesh%nbou_len, myMesh%nbou_dim_id)
         call check_err(iret)
      endif
      if(myMesh%nvel_len.ne.0) then
         iret = nf_def_dim(ncid, 'nvel',
     &          myMesh%nvel_len, myMesh%nvel_dim_id)
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
      if((myMesh%nope_len.ne.0).and.(myMesh%neta_len.ne.0)) then
         myMesh%nvdllnc_dims(1) = myMesh%nope_dim_id      
         iret = nf_def_var(ncid, 'nvdll',NF_INT, myMesh%nvdllnc_rank, 
     &                  myMesh%nvdllnc_dims, myMesh%nvdllnc_id)
         CALL check_err(iret)
         myMesh%nopenc_dims(1) = myMesh%neta_dim_id      
         iret = nf_def_var(ncid, 'nbdv',NF_INT, myMesh%nopenc_rank,
     &          myMesh%nopenc_dims, myMesh%nopenc_id)
         CALL check_err(iret)
      endif
!     Define normal flow boundary information
      if((myMesh%nbou_len.ne.0).and.(myMesh%nvel_len.ne.0)) then
         myMesh%nvellnc_dims(1) = myMesh%nbou_dim_id      
         iret = nf_def_var(ncid, 'nvell',NF_INT, myMesh%nvellnc_rank, 
     &                  myMesh%nvellnc_dims, myMesh%nvellnc_id)
         CALL check_err(iret)
         myMesh%ibtypenc_dims(1) = myMesh%nbou_dim_id      
         iret = nf_def_var(ncid, 'ibtype',NF_INT, myMesh%ibtypenc_rank, 
     &          myMesh%ibtypenc_dims, myMesh%ibtypenc_id)
         CALL check_err(iret)
         myMesh%nbounc_dims(1) = myMesh%nvel_dim_id      
         iret = nf_def_var(ncid, 'nbvv',NF_INT, myMesh%nbounc_rank,
     &          myMesh%nbounc_dims, myMesh%nbounc_id)
         CALL check_err(iret)
      endif
!     -------------------
!     Define Z coordinate
!     --------------------
      myMesh%DEPTH_dims(1) = myMesh%num_nodes_dim_id
      iret = nf_def_var(ncid, 'depth',NF_DOUBLE, myMesh%DEPTH_rank, 
     &       myMesh%DEPTH_dims, myMesh%DEPTH_id)
      CALL check_err(iret)
C-----------------------------------------------------------------------
      END SUBROUTINE defineMeshVariables
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E   P U T   M E S H   V A R I A B L E S 
C-----------------------------------------------------------------------
C     jgf49.17.02 Writes data that are common to all netCDF files
C     to the specified file.
C-----------------------------------------------------------------------
      subroutine putMeshVariables(ncid, myMesh)
      USE SIZES, ONLY : SZ
      USE GLOBAL, ONLY : ICS, X, Y, DP
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
!     Define longitude attributes
      if(ics.eq.2) then  ! this indicates spherical coordinates
      
         iret = nf_put_att_text(ncid, myMesh%X_id, 
     &      'long_name',9,'longitude')
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, myMesh%X_id,'standard_name', 9,
     &       'longitude')
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, myMesh%X_id,
     &       'units', 12, VAR_DEGE)
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, myMesh%X_id,
     &       'positive', 4, VAR_EAST)
         CALL check_err(iret)
!     Define latitude attributes
         iret = nf_put_att_text(ncid, myMesh%Y_id,
     &      'long_name', 8, 'latitude')
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, myMesh%Y_id,
     &      'standard_name', 8, 'latitude')
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, myMesh%Y_id,
     &      'units', 13, VAR_DEGN)
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, myMesh%Y_id,
     &      'positive', 5, VAR_NORTH)
         CALL check_err(iret)
      else   ! must be using Cartesian (x,y) coordinates
      
!     Define x-coordinate attributes
!     ------------------------------
         iret = nf_put_att_text(ncid, myMesh%X_id,
     &      'long_name', 22, VAR_L_X)
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, myMesh%X_id, 'standard_name',
     &                          12, VAR_S_X)
         CALL check_err(iret)
!    -------------------------
!     determine variable units
!    -------------------------
         CALL putUnitsAttribute(ncid, myMesh%Y_id, VAR_METS, VAR_FEET)

         iret = nf_put_att_text(ncid, myMesh%X_id, 'positive', 5,  VAR_RG)
         CALL check_err(iret)
!     ------------------------------
!     Define y-coordinate attributes
!     ------------------------------     
         iret = nf_put_att_text(ncid, myMesh%Y_id, 'long_name', 22, VAR_L_Y)
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, myMesh%Y_id, 'standard_name', 
     &                          12, VAR_S_Y)
         CALL check_err(iret)
          
         CALL putUnitsAttribute(ncid, myMesh%Y_id, VAR_METS, VAR_FEET)

         iret = nf_put_att_text(ncid, myMesh%Y_id, 'positive', 34, VAR_UPY)
         CALL check_err(iret)
      endif
!     ------------------------------
!     Define depth attributes
!     ------------------------------
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
!     ---------------------------------------------------------------
!     Define node numbers on normal flow boundary segments attributes
!     ----------------------------------------------------------------  
      iret = nf_put_att_text(ncid, myMesh%nopenc_id, 'long_name', 57, 
     &      'node numbers on all elevation specified boundary segments')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, myMesh%nopenc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, myMesh%nvdllnc_id, 'long_name', 60, 
     &'number of nodes in each elevation specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, myMesh%nvdllnc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, myMesh%nbounc_id, 'long_name', 58, 
     & 'node numbers in all normal flow specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, myMesh%nbounc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, myMesh%ibtypenc_id, 'long_name', 13, 
     &       'boundary type')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, myMesh%ibtypenc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, myMesh%nvellnc_id, 'long_name', 62, 
     &'number of nodes in each normal flow specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, myMesh%nvellnc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
!
!     Store X
      iret = nf_put_var_double(ncid, myMesh%X_id, X)
      CALL check_err(iret)
!     Store Y
      iret = nf_put_var_double(ncid, myMesh%Y_id, Y)
      CALL check_err(iret)
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
      if(myMesh%neta_len.ne.0) then
         iret = nf_put_var_int(ncid, myMesh%nopenc_id, myMesh%nopenc)
         CALL check_err(iret)
      endif
      if(myMesh%nope_len.ne.0) then
         iret = nf_put_var_int(ncid, myMesh%nvdllnc_id ,myMesh%nvdllnc)
         CALL check_err(iret)
      endif
!     Store normal flow boundary information
      if(myMesh%nvel_len.ne.0) then
         iret = nf_put_var_int(ncid, myMesh%nbounc_id, myMesh%nbounc)
         CALL check_err(iret)
      endif
      if(myMesh%nbou_len.ne.0) then
         iret = nf_put_var_int(ncid,myMesh%ibtypenc_id, myMesh%ibtypenc)
         CALL check_err(iret)
      endif
      if(myMesh%nbou_len.ne.0) then
         iret = nf_put_var_int(ncid, myMesh%nvellnc_id, myMesh%nvellnc)
         CALL check_err(iret)
      endif
C-----------------------------------------------------------------------
      END SUBROUTINE putMeshVariables 
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   P U T   M E T A D A T A 
C-----------------------------------------------------------------------
C     jgf49.17.02 Writes data that are common to all netCDF files
C     to the specified file.
C-----------------------------------------------------------------------
      subroutine putMetaData(ncid)
      USE VERSION
      USE GLOBAL, ONLY : RUNDES, AGRID, title, institution, source,
     &                   history, references, comments, host,
     &                   convention, contact, dtdp, ihot, 
     &                   ics, nolifa, nolica, nolicat, 
     &                   ncor, ntip, nws, nramp, statim,
     &                   reftim, rnday, dramp, a00, b00, c00, h0,
     &                   cori, ntif, nbfr
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
     &       LEN(TRIM(rundes(1:INDEX(rundes,"!")-1))), rundes)      
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'grid',
     &       LEN(TRIM(agrid)), agrid)       
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'title',
     &       LEN(TRIM(title(1:INDEX(title,"!")-1))), title)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'institution', 
     &       LEN(TRIM(institution(1:INDEX(institution,"!")-1))), 
     &       institution)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'source', 
     &       LEN(TRIM(source(1:INDEX(source,"!")-1))), source)     
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'history', 
     &       LEN(TRIM(history(1:INDEX(history,"!")-1))),  history)    
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'references', 
     &       LEN(TRIM(references(1:INDEX(references,"!")-1))), 
     &       references)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'comments', 
     &       LEN(TRIM(comments(1:INDEX(comments,"!")-1))), 
     &       comments)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'host', 
     &       LEN(TRIM(host(1:INDEX(host,"!")-1))),  host)    
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'convention', 
     &       LEN(TRIM(convention(1:INDEX(convention,"!")-1))), 
     &       convention)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'contact', 
     &       LEN(TRIM(contact(1:INDEX(contact,"!")-1))), contact)
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

C-----------------------------------------------------------------------
      END SUBROUTINE putMetaData 
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   U P D A T E   M E T A   D A T A 
C-----------------------------------------------------------------------
C     jgf49.17.02 Updates data that are common to all netCDF files
C     in the specified file.
C-----------------------------------------------------------------------
      subroutine updateMetaData(ncid)
      USE SIZES, ONLY : SZ
      USE GLOBAL, ONLY : ihot, nramp, rnday, dramp, a00, b00, c00, h0
      USE NodalAttributes, ONLY : nwp, nolibf, tau0, cf, eslm
      USE VERSION
      IMPLICIT NONE
      INTEGER, intent(in) :: ncid 
C
      INTEGER iret   ! Error status return
!     date_string variables for time attribute
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer values(8)

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
C-----------------------------------------------------------------------
      END SUBROUTINE updateMetaData 
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   P U T   U N I T   A T T R I B U T E S
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
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, var_id, 'units',
     &           len_trim(metric), trim(metric))
      else
          iret = nf_put_att_text(ncid, var_id, 'units',
     &           len_trim(english), trim(english)) 
      endif
      CALL check_err(iret)
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
      if (iret .ne. NF_NOERR) then
         write(ScreenUnit,*) "ERROR: NetCDF: ",nf_strerror(iret)
         write(ScreenUnit,*) "ERROR: ADCIRC execution terminated."
         write(16,*) "ERROR: NetCDF: ",nf_strerror(iret)
         write(16,*) "ERROR: ADCIRC execution terminated."
#ifdef CMPI
         call msg_fini()
#endif
         stop
      endif
C-----------------------------------------------------------------------
      end subroutine check_err
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      END MODULE NETCDFIO
C-----------------------------------------------------------------------
