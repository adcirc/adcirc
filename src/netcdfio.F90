!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2025 R.A. Luettich, Jr., J.J. Westerink
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------!
#include "logging_macros.h"
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
module NETCDFIO
   use SIZES, only: OFF, ASCII, SPARSE_ASCII, NETCDF3, NETCDF4, XDMF
   use NETCDF_ERROR, only: CHECK_ERR
   use mod_logging, only: DEBUG, ECHO, INFO, WARNING, ERROR, &
                          screenMessage, logMessage, allMessage, t_log_scope, init_log_scope

   use NETCDF, only: NF90_NOERR, NF90_NOWRITE, NF90_WRITE, NF90_CLOBBER, &
                     NF90_NOCLOBBER, NF90_CHAR, NF90_DOUBLE, NF90_INT, &
                     NF90_FLOAT, NF90_UBYTE, NF90_INT64, NF90_UNLIMITED, &
                     NF90_CLASSIC_MODEL, NF90_HDF5, NF90_NETCDF4, NF90_GLOBAL, &
                     nf90_create, nf90_open, nf90_close, nf90_inq_varid, &
                     nf90_def_var, nf90_put_var, nf90_get_var, nf90_def_dim, &
                     nf90_inq_dimid, nf90_put_att, nf90_get_att, nf90_redef, nf90_inquire, &
                     nf90_def_dim, nf90_strerror, nf90_inquire_dimension, nf90_enddef
#ifdef NETCDF_CAN_DEFLATE
   use NETCDF, only: nf90_def_var_deflate
#endif

   implicit none

   private

   real(8), parameter ::  doubleval = -99999.d0
   integer, parameter :: intFillValue = -99999

   type, private :: meshStructure
      logical :: initialized = .false.
      integer :: X_id ! x-coordinate or longitude
      integer :: Y_id ! y-coordinate or latitude
      integer :: sigma_id ! relative depth of 3D layer
      integer :: DEPTH_id ! distance from geoid
      integer :: ELE_id ! elements in grid
      integer :: nbdvnc_id ! nodes on elev spec boundary seg
      integer :: nbvvnc_id ! nodes on normal flow boundary seg
      integer :: ibconn_id ! connected internal boundary node
      integer :: barht_id ! barrier height
      integer :: barsp_id ! barrier supercritical flow coef
      integer :: barsb_id ! barrier subcritical flow coef
      integer :: pipeht_id ! barrier pipe height id
      integer :: pipecoef_id ! barrier pipe coefficient id
      integer :: pipediam_id ! barrier pipe diameter id
      integer :: nvdllnc_id ! num nodes on elev boundary seg
      integer :: ibtypenc_id ! discharge boundary type
      integer :: ibtypeenc_id ! elevation boundary type
      integer :: nvellnc_id ! nodes on norm flow spec boundary seg
      integer :: slam0nc_id
      integer :: sfea0nc_id
      integer :: netanc_id
      integer :: nvelnc_id
      integer :: mesh_id ! Mesh Definition (Corbitt)
!     Dimension ids
      integer :: num_nodes_dim_id
      integer :: num_elems_dim_id
      integer :: nface_dim_id
      integer :: nopenc_dim_id ! num elev spec boundary forcing segs
      integer :: nbounc_dim_id ! number of normal flow specified boundary segment
      integer :: ibtypenc_dim_id
      integer :: ibtypeenc_dim_id ! elevation boundary types
      integer :: max_nvdllnc_dim_id
      integer :: max_nvellnc_dim_id
      integer :: num_v_nodes_dim_id
      integer :: nvelnc_dim_id
      integer :: netanc_dim_id
!     Dimension lengths
      integer :: num_nodes
      integer :: num_elems
      integer :: nface_len
      integer :: nopenc
      integer :: nbounc
      integer :: max_nvdllnc
      integer :: max_nvellnc
      integer :: num_v_nodes
!     Rank (number of dimensions) for each variable
!     Variable shapes
      integer :: x_dims(1)
      integer :: y_dims(1)
      integer :: sigma_dims(1)
      integer :: depth_dims(1)
      integer :: ele_dims(2)
      integer :: nvdll_dims(1)
      integer :: nbounc_dims(1)
      integer :: ibtypenc_dims(1)
      integer :: ibtypeenc_dims(1)
      integer :: nvellnc_dims(1)
      integer :: nvdllnc_dims(1)
      integer :: nopenc_dims(1)
      integer :: nbdvnc_dims(2)
      integer :: nbvvnc_dims(2)
      real(8), allocatable :: xnc(:) ! x coordinate or longitude
      real(8), allocatable :: ync(:) ! y coordinate or latitude
      integer, allocatable ::  nbvvnc(:) ! boundary array
      integer, allocatable ::  ibconnnc(:) ! boundary array
      integer, allocatable ::  nbdvnc(:) ! boundary array
      integer, allocatable ::  nvellnc(:) ! boundary array
      integer, allocatable ::  nvdllnc(:) ! boundary array
      integer, allocatable ::  ibtypenc(:) ! boundary array
      integer, allocatable ::  ibtypeenc(:)
      integer, allocatable ::  nmnc(:, :)
      integer, allocatable ::  element(:, :)
      real(8), allocatable ::  barht(:)
      real(8), allocatable ::  barsp(:)
      real(8), allocatable ::  barsb(:)
      real(8), allocatable ::  pipeht(:)
      real(8), allocatable ::  pipediam(:)
      real(8), allocatable ::  pipecoef(:)
      integer :: netanc
      integer :: nvelnc
      integer :: const_id
      integer :: constl_dim_id
      logical :: hasWeirs
      logical :: hasInternalWeirs
      logical :: hasPipes
   end type meshStructure

   type(meshStructure), private, target :: adcircMesh

   type, private :: fileData
      integer :: record_counter
      integer :: ncformat ! netcdf file format to create
      logical :: createFile ! .true. if a new netCDF file must be created p
      character(len=1024) :: filename
      logical :: fileFound ! .true. if the netCDF file is present
   end type fileData

   type, private :: timeData
      logical :: initialized = .false.
      integer :: timenc_len = 1 ! number of time slices to write
      integer :: timenc_dim_id
      integer :: timenc_id
      integer :: timenc_dims(1)
      real(8), allocatable :: timenc(:)
   end type timeData

   type, private :: stationData
      integer :: ncid ! the id of its netcdf file
      integer :: num_stations ! total number of stations
      integer :: num_v_nodes
      integer :: num_sta_dim_id
      integer :: num_v_nodes_dim_id
      character(50) :: varname(3) ! name(s) of this data in netcdf file
      integer, allocatable :: name_lengths(:) ! lengths of station names
      character(50), pointer :: statnames(:)
      integer :: slen_dim_id
      integer :: station_id
      integer :: station_dims(2)
      integer :: x_id ! station x-coordinate or longitude
      integer :: y_id ! station y-coordinate or latitude
      integer :: x_dims(1)
      integer :: y_dims(1)
      real(8), allocatable :: x(:) ! x coordinate or longitude
      real(8), allocatable :: y(:) ! y coordinate or latitude
      integer :: station_data_id
      integer :: u_station_data_id
      integer :: v_station_data_id
      integer :: w_station_data_id
      integer :: station_data_dims(2)
      integer :: station_data_dims_3D(3)
      integer :: const_id
      integer :: constl_dim_id
      integer :: station_ha_data_id
      integer :: ha_u_station_data_id
      integer :: ha_v_station_data_id
      integer :: station_hg_data_id
      integer :: hg_u_station_data_id
      integer :: hg_v_station_data_id
      type(timeData) :: myTime
      type(fileData) :: myFile
   end type stationData

   type(stationData), private :: elevSta ! elev stations (fort.61)
   type(stationData), private :: prSta ! pressure stations (fort.71)
   type(stationData), private :: IceSta ! ice stations (fort.91)

   type(stationData), private :: velSta ! velocity stations (fort.62)
   type(stationData), private :: wVelSta ! wind vel stations (fort.72)
   type(stationData), private :: densityStations3D ! for fort.41
   type(stationData), private :: velocityStations3D ! for fort.42
   type(stationData), private :: turbulenceStations3D ! for fort.43
   type(stationData), private :: HAelevSta ! har elev stations (fort.51)
   type(stationData), private :: HAvelSta ! har vel stations (fort.52)
   type(stationData), private :: dynamiccorrectionSta ! for dynamicwaterlevelcorrection.61

   type, private :: nodalData
      integer :: ncid ! the id of its netcdf file
      real(8) :: initial_value ! array will be initialized to this
      integer :: nodal_data_id
      integer :: ha_nodal_data_id
      integer :: hg_nodal_data_id
      integer :: u_ha_nodal_data_id
      integer :: u_hg_nodal_data_id
      integer :: v_ha_nodal_data_id
      integer :: v_hg_nodal_data_id
      integer :: u_nodal_data_id
      integer :: v_nodal_data_id
      integer :: w_nodal_data_id
      integer :: max_nodal_data_id
      integer :: time_max_nodal_data_id !tcm v51.20.06 added
      integer :: nodal_data_dims(2)
      integer :: nodal_data_dims_max(1)
      integer :: nodal_data_dims_3D(3)
      character(50) :: varnames(6) ! name(s) of this data in netcdf file
      type(timeData) :: myTime
      type(fileData) :: myFile
      type(meshStructure), pointer :: myMesh
   end type nodalData

   type(nodalData), private :: elev ! for fort.63
   type(nodalData), private :: pr ! for fort.73
   type(nodalData), private :: wDrag ! for windDragOut
   type(nodalData), private :: currentVel ! for fort.64
   type(nodalData), private :: windVel ! for fort.74
   type(nodalData), private:: weirElev ! for fort.77
   type(nodalData), private:: tau0nc ! for fort.90
   type(nodalData), private, target :: EtaMax ! for maxele.63
   type(nodalData), private, target :: PrMin ! for minpr.63
   type(nodalData), private, target :: UMax ! for maxvel.63
   type(nodalData), private, target :: WVMax ! for maxwvel.63
   type(nodalData), private, target :: RSMax ! for rads_max.63
   type(nodalData), private :: rads ! for rads.64
   type(nodalData), private :: CICEAF ! for fort.93

   type(nodalData), private :: sw_hs ! for swan_HS.63
   type(nodalData), private :: sw_dir ! for swan_DIR.63
   type(nodalData), private :: sw_tm01 ! for swan_TM01.63
   type(nodalData), private :: sw_tps ! for swan_TPS.63
   type(nodalData), private :: sw_wind ! for swan_WIND.64
   type(nodalData), private :: sw_tm02 ! for swan_TM02.63
   type(nodalData), private :: sw_tmm10 ! for swan_tmm10.63
   type(nodalData), private :: sw_hs_max ! for swan_HS_max.63
   type(nodalData), private :: sw_dir_max ! for swan_DIR_max.63
   type(nodalData), private :: sw_tm01_max ! for swan_TM01_max.63
   type(nodalData), private :: sw_tps_max ! for swan_TPS_max.63
   type(nodalData), private :: sw_wind_max ! for swan_WIND_max.63
   type(nodalData), private :: sw_tm02_max ! for swan_TM02_max.63
   type(nodalData), private :: sw_tmm10_max ! for swan_TMM10_max.63

   type(nodalData), private :: density3D ! for fort.44
   type(nodalData), private :: velocity3D ! for fort.45
   type(nodalData), private :: turbulence3D ! for fort.46
   type(nodalData), private :: futureSurfaceTemperature ! for fort.47
   type(nodalData), private :: internal3D ! fort.48 variables
   type(nodalData), private :: HAelev ! for fort.53
   type(nodalData), private :: HAvel ! for fort.54

   type(nodalData), private:: coldDry ! for initiallydry.63
   type(nodalData), private, target :: maxInDep ! for maxinundepth.63
   type(nodalData), private, target :: inTime ! for inundationtime.63
   type(nodalData), private:: eRisInun ! for endrisinginun.63
   type(nodalData), private, target :: evrDry ! for everdried.63
   type(nodalData), private :: dynamiccorrection ! for dynamicwaterlevelcorrection.63

   type, private :: hotstartData
      integer :: lun ! ADCIRC's logical unit number
      integer :: ncid ! the id of its netcdf file
      integer :: FileFmtMajorFile !
      integer :: FileFmtMinorFile
      integer :: FileFmtRevFile
      type(fileData) :: myFile
      type(timeData) :: myTime
      type(meshStructure), pointer :: myMesh
!     2D simulation state
      type(nodalData) :: zeta1
      type(nodalData) :: zeta2
      type(nodalData) :: zetad
      type(nodalData) :: vel
      type(nodalData) :: ch1
      type(nodalData) :: nodecodenc
      type(nodalData) :: noffnc

      type(nodalData) :: htot1, htot2
      type(nodalData) :: rs1, rs2
      type(nodalData) :: swan_rs1, swan_rs2

!     3D simulation state
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
!     harmonic analysis components
      integer :: namefr_len = 10
      integer :: namefr_len_dim_id
      integer :: namefr_dims(2)
      integer :: mnharf_dim_id
      integer :: load_vector_dim_id ! 2x the number of frequencies
      integer :: component_dims(1)
      integer :: ha_dims(2)
      integer :: hafreq_id
      integer :: haff_id
      integer :: haface_id
      integer :: ha_id
      integer :: namefr_id
!     harmonic analysis load vectors
      type(nodalData) :: gloelv
      type(nodalData) :: glovellv
      type(stationData) :: staelv
      type(stationData) :: stavellv
!     harmonic analysis means and variance calculations
      type(nodalData) :: xvelav
      type(nodalData) :: yvelav
      type(nodalData) :: xvelva
      type(nodalData) :: yvelva
      type(nodalData) :: elav
      type(nodalData) :: elva
   end type hotstartData

   type(hotstartData), private, target :: hs67 ! for fort.67
   type(hotstartData), private, target :: hs68 ! for fort.68
   type(hotstartData), private, pointer :: hs ! current hs file

   ! lists of fulldomain node and element numbers that map to the
   ! nodes and elements of this subdomain in parallel
   integer, allocatable, target :: fullDomainNodeList(:)
   integer, allocatable, target :: fullDomainElementList(:)
   ! lists of fulldomain elevation and velocity station numbers that
   ! map to this subdomain's stations in parallel
   integer, allocatable, target :: fullDomainElevationStationList(:)
   integer, allocatable, target :: fullDomainVelocityStationList(:)
   ! used to specify the array of indexes to be used in a mapping
   integer, pointer :: fullDomainIndexList(:)
   ! usually initialized by readNetCDFHotstart(); otherwise in
   ! readAndMapToSubdomainMinMaxNetCDF()
   logical :: fullDomainIndexListsInitialized = .false.

   character(1024) :: att_text ! scratch variable for attribute text

   logical :: writerReadMetaData = .false. !...variable to allow writers to read
   !   NetCDF metadata from ADCPREP processed
   !   files on first pass

   public :: initNetCDFOutputFile, readNetCDFHotstart, &
             writeNetCDFHotstart, updateMetaData, &
             writeOutArrayNetCDF, initNetCDFHotstart, &
             initNetCDFHotstartHarmonic, initNetCDFHotstartHarmonicMeansVariances, &
             initNetCDFHotstart3D, writeNetCDFHotstartHarmonic, &
             writeNetCDFHotstartHarmonicMeansVariances, writeNetCDFHotstart3D, &
             writeNetCDFHotstart3DVar, readAndMapToSubdomainMaxMinNetCDF, &
             readNetCDFHotstartHarmonic, readNetCDFHotstartHarmonicMeansVariances, &
             readNetCDFHotstart3D, setADCIRCParameters, freeNetCDFCoord

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E   I N I T  N E T C D F  O U T P U T  F I L E
!-----------------------------------------------------------------------
!     jgf49.17.02 Allocates memory for NetCDF operations.
!-----------------------------------------------------------------------
   subroutine initNetCDFOutputFile(descript1, reterr)
      use GLOBAL, only: OutputDataDescript_t
      use GLOBAL_3DVS, only: outputFort48
      implicit none
      type(OutputDataDescript_t), intent(inout) :: descript1
      logical, intent(out) :: reterr
      character(1024) :: scratchMessage

      LOG_SCOPE_TRACED("initNetCDFOutputFile", NETCDFIO_TRACING)

      reterr = .false.
      ! Don't allocate or initialize anything if this is a pure output file,
      ! (not an input/output file like max/min files) and the format is
      ! specified as something other than netcdf.
      if (((abs(descript1%specifier) /= NETCDF3) .and. &
           (abs(descript1%specifier) /= NETCDF4) .and. &
           (descript1%readMaxMin .eqv. .false.)) .or. &
          (descript1%specifier == 0)) then
         return
      end if

      select case (descript1%lun)
      case (41)
         call initStationFile(densityStations3D, &
                              descript1, reterr)
      case (42)
         call initStationFile(velocityStations3D, &
                              descript1, reterr)
      case (43)
         call initStationFile(turbulenceStations3D, &
                              descript1, reterr)
      case (44)
         call initNodalDataFile(density3D, descript1, reterr)
      case (45)
         call initNodalDataFile(velocity3D, descript1, reterr)
      case (46)
         call initNodalDataFile(turbulence3D, descript1, reterr)
      case (47)
         call initNodalDataFile(futureSurfaceTemperature, &
                                descript1, reterr)
      case (48)
         if (outputFort48) then
            call initNodalDataFile(internal3D, descript1, reterr)
         end if
      case (51)
         call initStationFile(HAelevSta, descript1, reterr)
      case (52)
         call initStationFile(HAvelSta, descript1, reterr)
      case (53)
         call initNodalDataFile(HAelev, descript1, reterr)
      case (54)
         call initNodalDataFile(HAvel, descript1, reterr)
      case (61)
         call initStationFile(elevSta, descript1, reterr)
      case (62)
         call initStationFile(velSta, descript1, reterr)
      case (63)
         call initNodalDataFile(elev, descript1, reterr)
      case (64)
         call initNodalDataFile(currentVel, descript1, reterr)
      case (71)
         call initStationFile(prSta, descript1, reterr)
      case (72)
         call initStationFile(wVelSta, descript1, reterr)
      case (73)
         call initNodalDataFile(pr, descript1, reterr)
      case (74)
         call initNodalDataFile(windVel, descript1, reterr)
      case (91)
         call initStationFile(IceSta, descript1, reterr)
      case (93)
         call initNodalDataFile(CICEAF, descript1, reterr)
      case (173)
         call initNodalDataFile(wDrag, descript1, reterr)
      case (90)
         call initNodalDataFile(tau0nc, descript1, reterr)
      case (77)
         call initNodalDataFile(weirElev, descript1, reterr)
      case (164)
         call initNodalDataFile(rads, descript1, reterr)
      case (301)
         call initNodalDataFile(sw_hs, descript1, reterr)
      case (302)
         call initNodalDataFile(sw_dir, descript1, reterr)
      case (303)
         call initNodalDataFile(sw_tm01, descript1, reterr)
      case (304)
         call initNodalDataFile(sw_tps, descript1, reterr)
      case (305)
         call initNodalDataFile(sw_wind, descript1, reterr)
      case (306)
         call initNodalDataFile(sw_tm02, descript1, reterr)
      case (307)
         call initNodalDataFile(sw_tmm10, descript1, reterr)
      case (311)
         call initNodalDataFile(EtaMax, descript1, reterr)
      case (312)
         call initNodalDataFile(UMax, descript1, reterr)
      case (313)
         call initNodalDataFile(PrMin, descript1, reterr)
      case (314)
         call initNodalDataFile(WVMax, descript1, reterr)
      case (315)
         call initNodalDataFile(RSMax, descript1, reterr)
      case (316)
         call initNodalDataFile(sw_hs_max, descript1, reterr)
      case (317)
         call initNodalDataFile(sw_dir_max, descript1, reterr)
      case (318)
         call initNodalDataFile(sw_tm01_max, descript1, reterr)
      case (319)
         call initNodalDataFile(sw_tps_max, descript1, reterr)
      case (320)
         call initNodalDataFile(sw_wind_max, descript1, reterr)
      case (321)
         call initNodalDataFile(sw_tm02_max, descript1, reterr)
      case (322)
         call initNodalDataFile(sw_tmm10_max, descript1, reterr)
      case (400)
         call initNodalDataFile(inTime, descript1, reterr)
      case (401)
         call initNodalDataFile(maxInDep, descript1, reterr)
      case (402)
         call initNodalDataFile(coldDry, descript1, reterr)
      case (403)
         call initNodalDataFile(eRisInun, descript1, reterr)
      case (404)
         call initNodalDataFile(evrDry, descript1, reterr)
      case (108)
         call initNodalDataFile(dynamiccorrection, descript1, reterr)
      case (109)
         call initStationFile(dynamiccorrectionSta, descript1, reterr)
      case DEFAULT
         write (scratchMessage, &
                '("No netCDF for files with unit number ",i0,".")') descript1%lun
         call allMessage(ERROR, scratchMessage)
      end select

!-----------------------------------------------------------------------
   end subroutine initNetCDFOutputFile
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E    I N I T  S T A T I O N   F I L E
!-----------------------------------------------------------------------
!     jgf49.17.02 Sets up netCDF variables and allocates memory for a
!     station.
!-----------------------------------------------------------------------
   subroutine initStationFile(sta, descript1, reterr)
      use SIZES, only: MNPROC
      use ADC_CONSTANTS, only: RAD2DEG
      use GLOBAL, only: SNAMLEN, OutputDataDescript_t, &
                        NWS, STATNAME, STATNAMEV, &
                        STATNAMEM, IDEN
      use GLOBAL_3DVS, only: STATNAMED, STATNAMEV3D, STATNAMET
      use MESH, only: ICS
      implicit none

      type(stationData), intent(inout) :: sta
      type(OutputDataDescript_t), intent(inout) :: descript1
      logical, intent(out) :: reterr
      integer :: iret ! success or failure of the netcdf call
      character(1024) :: scratchMessage

!     date_string variables for time attribute
      integer, parameter :: stationLuns3D(3) = [41, 42, 43]
      integer, parameter :: stationLunsHA(2) = [51, 52]

      LOG_SCOPE_TRACED("initStationFile", NETCDFIO_TRACING)
      reterr = .false.

!     jgf50.13: if netcdf output was requested, but there are no stations,
!     don't initialize the file; just return.
      if (descript1%num_fd_records == 0) then
         return
      end if

      sta%num_stations = descript1%num_fd_records
      allocate (sta%myTime%timenc(sta%myTime%timenc_len))

      if (any(stationLuns3D == descript1%lun) &
          .or. any(stationLunsHA == descript1%lun)) then
         sta%num_v_nodes = descript1%num_items_per_record
      end if

!     Memory allocation for the station
      allocate (sta%x(sta%num_stations))
      allocate (sta%y(sta%num_stations))

!     Initialize netCDF output file, creating a new one if necessary.
      call createNetCDFOutputFile(sta%ncid, sta%myFile, sta%myTime, &
                                  descript1, reterr)

!     if we didn't need to create a file, update metadata and return
!     !jgf52.08.22: Only update metadata if running in serial; in
!     !parallel the metadata was updated by adcprep
      if (sta%myFile%createFile .eqv. .false.) then
         if ((reterr .eqv. .false.) .and. (mnproc == 1)) then
            call updateMetaData(sta%ncid, sta%myFile)
         end if
         return
      end if

!     Set coordinates of each station, converting to degrees if we
!     are in spherical coordinates.
      if (ICS /= 1) then
         sta%x = descript1%x_coord*RAD2DEG
         sta%y = descript1%y_coord*RAD2DEG
      else
         sta%x = descript1%x_coord
         sta%y = descript1%y_coord
      end if

!     Create station dimension and station name dimension
      iret = nf90_def_dim(sta%ncid, 'station', &
                          sta%num_stations, sta%num_sta_dim_id)
      call check_err(iret)
      iret = nf90_def_dim(sta%ncid, 'namelen', SNAMLEN, sta%slen_dim_id)
      call check_err(iret)
      if (any(stationLuns3D == descript1%lun)) then
         iret = nf90_def_dim(sta%ncid, 'num_v_nodes', sta%num_v_nodes, &
                             sta%num_v_nodes_dim_id)
         sta%station_data_dims_3D(1) = sta%num_sta_dim_id
         sta%station_data_dims_3D(2) = sta%num_v_nodes_dim_id
         sta%station_data_dims_3D(3) = sta%myTime%timenc_dim_id
!     !WJP 02.20.2018 Writing in constituent information for harmonic file
      elseif (any(stationLunsHA == descript1%lun)) then
         call defineHarmonicAnalysisParametersInNetcdfFile(descript1, stationdat=sta)
      end if
      call check_err(iret)

!     Define stations name
      sta%station_dims(1) = sta%slen_dim_id
      sta%station_dims(2) = sta%num_sta_dim_id

      iret = nf90_def_var(sta%ncid, 'station_name', NF90_CHAR, &
                          sta%station_dims, sta%station_id)
      call check_err(iret)

!     Define station locations
      sta%x_dims(1) = sta%num_sta_dim_id
      iret = nf90_def_var(sta%ncid, 'x', NF90_DOUBLE, &
                          sta%x_dims, sta%x_id)
      call check_err(iret)
      sta%y_dims(1) = sta%num_sta_dim_id
      iret = nf90_def_var(sta%ncid, 'y', NF90_DOUBLE, &
                          sta%y_dims, sta%y_id)
      call check_err(iret)

!     Set coordinates as representing either latitude or longitude,
!     or Cartesian x and y, depending on the value of ICS.
      call defineCoordinateAttributes(sta%ncid, sta%x_id, sta%y_id)

!     Fill in labels and populate variables as appropriate for the
!     different types of data in the station files. The labels and
!     units will also vary according to the coordinate system ADCIRC
!     is using (spherical or cartesian, according to the value of ICS)
!     as well as the units system (english or si according to the value of g).
      select case (descript1%lun)
      case (41) !       F O R T . 4 1
         iret = nf90_def_var(sta%ncid, 'sigmat', NF90_DOUBLE, &
                             sta%station_data_dims_3D, sta%u_station_data_id)
         call check_err(iret)
         if ((IDEN == 2) .or. (IDEN == 4)) then
            iret = nf90_def_var(sta%ncid, 'salinity', NF90_DOUBLE, &
                                sta%station_data_dims_3D, sta%v_station_data_id)
            call check_err(iret)
         end if
         if ((IDEN == 3) .or. (IDEN == 4)) then
            iret = nf90_def_var(sta%ncid, 'temperature', NF90_DOUBLE, &
                                sta%station_data_dims_3D, sta%w_station_data_id)
            call check_err(iret)
         end if
         ! sigma t
         att_text = &
            'station water column vertically varying density'
         iret = nf90_put_att(sta%ncid, &
                             sta%u_station_data_id, 'long_name', trim(att_text))
         call check_err(iret)
         att_text = 'station_density_vertically_varying'
         iret = nf90_put_att(sta%ncid, &
                             sta%u_station_data_id, 'standard_name', &
                             trim(att_text))
         call check_err(iret)
         call putUnitsAttribute(sta%ncid, sta%u_station_data_id, &
                                "kg/m^3", "n/a")
         ! salinity
         if ((IDEN == 2) .or. (IDEN == 4)) then
            att_text = &
               'station water column vertically varying salinity'
            iret = nf90_put_att(sta%ncid, &
                                sta%v_station_data_id, 'long_name', trim(att_text))
            call check_err(iret)
            att_text = &
               'station_water_salinity_vertically_varying'
            iret = nf90_put_att(sta%ncid, &
                                sta%v_station_data_id, 'standard_name', &
                                trim(att_text))
            call check_err(iret)
            call putUnitsAttribute(sta%ncid, sta%v_station_data_id, &
                                   "PSU", "n/a")
         end if
         ! temperature
         if ((IDEN == 3) .or. (IDEN == 4)) then
            att_text = &
               'station water column vertically varying temperature'
            iret = nf90_put_att(sta%ncid, &
                                sta%w_station_data_id, 'long_name', &
                                trim(att_text))
            call check_err(iret)
            att_text = &
               'station_water_temperature_vertically_varying'
            iret = nf90_put_att(sta%ncid, &
                                sta%w_station_data_id, 'standard_name', &
                                trim(att_text))
            call check_err(iret)
            call putUnitsAttribute(sta%ncid, sta%w_station_data_id, &
                                   "Celsius", "Fahrenheit")
         end if
         sta%statnames => STATNAMED

      case (42) !       F O R T . 4 2
         iret = nf90_def_var(sta%ncid, 'u-vel3D', NF90_DOUBLE, &
                             sta%station_data_dims_3D, sta%u_station_data_id)
         call check_err(iret)
         iret = nf90_def_var(sta%ncid, 'v-vel3D', NF90_DOUBLE, &
                             sta%station_data_dims_3D, sta%v_station_data_id)
         call check_err(iret)
         iret = nf90_def_var(sta%ncid, 'w-vel3D', NF90_DOUBLE, &
                             sta%station_data_dims_3D, sta%w_station_data_id)
         call check_err(iret)
         if (ics /= 1) then
            ! u
            att_text = &
               'station water column vertically varying east/west velocity'
            iret = nf90_put_att(sta%ncid, &
                                sta%u_station_data_id, 'long_name', trim(att_text))
            call check_err(iret)
            att_text = &
               'station_eastward_water_velocity_vertically_varying'
            iret = nf90_put_att(sta%ncid, &
                                sta%u_station_data_id, &
                                'standard_name', trim(att_text))
            call check_err(iret)
            iret = nf90_put_att(sta%ncid, &
                                sta%u_station_data_id, 'positive', 'east')
            call check_err(iret)
            ! v
            att_text = &
               'station water column vertically varying north/south velocity'
            iret = nf90_put_att(sta%ncid, &
                                sta%v_station_data_id, 'long_name', &
                                trim(att_text))
            call check_err(iret)
            att_text = &
               'station_northward_water_velocity_vertically_varying'
            iret = nf90_put_att(sta%ncid, &
                                sta%v_station_data_id, 'standard_name', &
                                trim(att_text))
            call check_err(iret)
            iret = nf90_put_att(sta%ncid, &
                                sta%v_station_data_id, 'positive', 'north')
            call check_err(iret)
            ! w
            att_text = &
               'station water column vertically varying up/down velocity'
            iret = nf90_put_att(sta%ncid, &
                                sta%w_station_data_id, 'long_name', &
                                trim(att_text))
            call check_err(iret)
            att_text = &
               'station_upward_water_velocity_vertically_varying'
            iret = nf90_put_att(sta%ncid, &
                                sta%w_station_data_id, 'standard_name', &
                                trim(att_text))
            call check_err(iret)
            iret = nf90_put_att(sta%ncid, &
                                sta%w_station_data_id, 'positive', 'up')
            call check_err(iret)
         else
            ! u
            att_text = &
               'station water column vertically varying velocity in x-direction'
            iret = nf90_put_att(sta%ncid, &
                                sta%u_station_data_id, 'long_name', trim(att_text))
            call check_err(iret)
            att_text = 'station_x_water_velocity_vertically_varying'
            iret = nf90_put_att(sta%ncid, &
                                sta%u_station_data_id, 'standard_name', &
                                trim(att_text))
            call check_err(iret)
            iret = nf90_put_att(sta%ncid, &
                                sta%u_station_data_id, 'positive', 'right')
            call check_err(iret)
            ! v
            att_text = &
               'station water column vertically varying velocity in y-direction'
            iret = nf90_put_att(sta%ncid, &
                                sta%v_station_data_id, 'long_name', trim(att_text))
            call check_err(iret)
            att_text = 'station_y_water_velocity_vertically_varying'
            iret = nf90_put_att(sta%ncid, &
                                sta%v_station_data_id, 'standard_name', &
                                trim(att_text))
            call check_err(iret)
            iret = nf90_put_att(sta%ncid, &
                                sta%v_station_data_id, 'positive', &
                                '90 degrees counterclockwise from x water velocity')
            call check_err(iret)
            ! w
            att_text = 'station water column ' &
                       //'vertically varying velocity in z-direction'
            iret = nf90_put_att(sta%ncid, &
                                sta%w_station_data_id, 'long_name', trim(att_text))
            call check_err(iret)
            att_text = 'station_z_water_velocity_vertically_varying'
            iret = nf90_put_att(sta%ncid, &
                                sta%w_station_data_id, 'standard_name', &
                                trim(att_text))
            call check_err(iret)
            iret = nf90_put_att(sta%ncid, &
                                sta%w_station_data_id, 'positive', &
                                '90 degrees counterclockwise from x water velocity')
            call check_err(iret)
         end if
         call putUnitsAttribute(sta%ncid, sta%u_station_data_id, &
                                'm s-1', 'ft s-1')
         call putUnitsAttribute(sta%ncid, sta%v_station_data_id, &
                                'm s-1', 'ft s-1')
         call putUnitsAttribute(sta%ncid, sta%w_station_data_id, &
                                'm s-1', 'ft s-1')
         sta%statnames => STATNAMEV3D

      case (43) !       F O R T . 4 3
         iret = nf90_def_var(sta%ncid, 'q20', NF90_DOUBLE, &
                             sta%station_data_dims_3D, sta%u_station_data_id)
         call check_err(iret)
         iret = nf90_def_var(sta%ncid, 'l', NF90_DOUBLE, &
                             sta%station_data_dims_3D, sta%v_station_data_id)
         call check_err(iret)
         iret = nf90_def_var(sta%ncid, 'ev', NF90_DOUBLE, &
                             sta%station_data_dims_3D, sta%w_station_data_id)
         call check_err(iret)
         ! q20
         att_text = &
            'station water column vertically varying turbulent'// &
            'kinetic energy'
         iret = nf90_put_att(sta%ncid, &
                             sta%u_station_data_id, 'long_name', trim(att_text))
         call check_err(iret)
         att_text = &
            'station_turbulent_kinetic_energy_vertically_varying'
         iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                             'standard_name', trim(att_text))
         call check_err(iret)
         ! l
         att_text = &
            'station water column vertically varying mixing length'
         iret = nf90_put_att(sta%ncid, &
                             sta%v_station_data_id, 'long_name', &
                             trim(att_text))
         call check_err(iret)
         att_text = &
            'station_water_mixing_length_vertically_varying'
         iret = nf90_put_att(sta%ncid, &
                             sta%v_station_data_id, 'standard_name', &
                             trim(att_text))
         call check_err(iret)
         !  ev
         att_text = &
            'station water column vertically varying eddy viscosity'
         iret = nf90_put_att(sta%ncid, &
                             sta%w_station_data_id, 'long_name', &
                             trim(att_text))
         call check_err(iret)
         att_text = &
            'station_water_eddy_viscosity_vertically_varying'
         iret = nf90_put_att(sta%ncid, &
                             sta%w_station_data_id, 'standard_name', &
                             trim(att_text))
         call check_err(iret)
!            CALL putUnitsAttribute(sta%ncid,sta%u_station_data_id,
!     &      'm s-1', 'ft s-1') !TODO: jgf49.48.01: units for q20??
         call putUnitsAttribute(sta%ncid, sta%v_station_data_id, &
                                "meters", "n/a")
         call putUnitsAttribute(sta%ncid, sta%w_station_data_id, &
                                "m^2/s", "n/a")
         sta%statnames => STATNAMET

      case (51) !       F O R T . 5 1        WJP 02.20.2018

         sta%station_data_dims(2) = sta%num_sta_dim_id
         sta%station_data_dims(1) = sta%num_v_nodes_dim_id
!        Put the amplitude stuff in
         iret = nf90_def_var(sta%ncid, 'amp', NF90_DOUBLE, &
                             sta%station_data_dims, sta%station_ha_data_id)
         call check_err(iret)
!        Define water surface elevation attributes
         iret = nf90_put_att(sta%ncid, sta%station_ha_data_id, &
                             'long_name', &
                             'station amplitude of tidal harmonic constituents')
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%station_ha_data_id, &
                             'standard_name', 'sta_hc_amp')
         call check_err(iret)
         call putUnitsAttribute(sta%ncid, sta%station_ha_data_id, &
                                'm', 'ft')
         iret = nf90_put_att(sta%ncid, sta%station_ha_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
!        Put the phase stuff in
         iret = nf90_def_var(sta%ncid, 'phs', NF90_DOUBLE, &
                             sta%station_data_dims, sta%station_hg_data_id)
         call check_err(iret)
!        Define water surface elevation attributes
         iret = nf90_put_att(sta%ncid, sta%station_hg_data_id, &
                             'long_name', &
                             'station phase of tidal harmonic constituents')
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%station_hg_data_id, &
                             'standard_name', 'sta_hc_phs')
         call check_err(iret)
         call putUnitsAttribute(sta%ncid, sta%station_hg_data_id, &
                                'deg', 'deg')
         iret = nf90_put_att(sta%ncid, sta%station_hg_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)

         sta%statnames => STATNAME

      case (52) !       F O R T . 5 2        WJP 02.20.2018

         sta%station_data_dims(2) = sta%num_sta_dim_id
         sta%station_data_dims(1) = sta%num_v_nodes_dim_id
         iret = nf90_def_var(sta%ncid, 'u-vel-amp', NF90_DOUBLE, &
                             sta%station_data_dims, sta%ha_u_station_data_id)
         call check_err(iret)
         iret = nf90_def_var(sta%ncid, 'u-vel-phs', NF90_DOUBLE, &
                             sta%station_data_dims, sta%hg_u_station_data_id)
         call check_err(iret)
         iret = nf90_def_var(sta%ncid, 'v-vel-amp', NF90_DOUBLE, &
                             sta%station_data_dims, sta%ha_v_station_data_id)
         call check_err(iret)
         iret = nf90_def_var(sta%ncid, 'v-vel-phs', NF90_DOUBLE, &
                             sta%station_data_dims, sta%hg_v_station_data_id)
         call check_err(iret)
!              ! u amp
         iret = nf90_put_att(sta%ncid, sta%ha_u_station_data_id, &
                             'long_name', &
                             'station tidal harmonic amplitude of east/west velocity')
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%ha_u_station_data_id, &
                             'standard_name', &
                             'sta_hc_amp_east_vel')
         call check_err(iret)
         call putUnitsAttribute(sta%ncid, sta%ha_u_station_data_id, &
                                'm/s', 'ft/s')
         iret = nf90_put_att(sta%ncid, sta%ha_u_station_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%ha_u_station_data_id, &
                             'dry_Value', doubleval)
         call check_err(iret)
!              ! u phs
         iret = nf90_put_att(sta%ncid, sta%hg_u_station_data_id, &
                             'long_name', &
                             'station tidal harmonic phase of east/west velocity')
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%hg_u_station_data_id, &
                             'standard_name', &
                             'sta_hc_phs_east_vel')
         call check_err(iret)
         call putUnitsAttribute(sta%ncid, sta%hg_u_station_data_id, &
                                'deg', 'deg')
         iret = nf90_put_att(sta%ncid, sta%hg_u_station_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%hg_u_station_data_id, &
                             'dry_Value', doubleval)
         call check_err(iret)
!              ! v amp
         iret = nf90_put_att(sta%ncid, sta%ha_v_station_data_id, &
                             'long_name', &
                             'station tidal harmonic amplitude of north/south velocity')
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%ha_v_station_data_id, &
                             'standard_name', &
                             'sta_hc_amp_north_vel')
         call check_err(iret)
         call putUnitsAttribute(sta%ncid, sta%ha_v_station_data_id, &
                                'm/s', 'ft/s')
         iret = nf90_put_att(sta%ncid, sta%ha_v_station_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%ha_v_station_data_id, &
                             'dry_Value', doubleval)
         call check_err(iret)
!              ! v phs
         iret = nf90_put_att(sta%ncid, sta%hg_v_station_data_id, &
                             'long_name', &
                             'station tidal harmonic phase of north/south velocity')
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%hg_v_station_data_id, &
                             'standard_name', &
                             'sta_hc_phs_north_vel')
         call check_err(iret)
         call putUnitsAttribute(sta%ncid, sta%hg_v_station_data_id, &
                                'deg', 'deg')
         iret = nf90_put_att(sta%ncid, sta%hg_v_station_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%hg_v_station_data_id, &
                             'dry_Value', doubleval)
         call check_err(iret)

         sta%statnames => STATNAMEV

      case (61) !       F O R T . 6 1

         sta%station_data_dims(1) = sta%num_sta_dim_id
         sta%station_data_dims(2) = sta%myTime%timenc_dim_id
         iret = nf90_def_var(sta%ncid, 'zeta', NF90_DOUBLE, &
                             sta%station_data_dims, sta%station_data_id)
         call check_err(iret)
!        Define water surface elevation attributes
         iret = nf90_put_att(sta%ncid, sta%station_data_id, &
                             'long_name', 'water surface elevation above geoid')
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%station_data_id, &
                             'standard_name', 'sea_surface_height_above_geoid')
         call check_err(iret)
         call putUnitsAttribute(sta%ncid, sta%station_data_id, &
                                'm', 'ft')
         iret = nf90_put_att(sta%ncid, sta%station_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
!            iret = nf90_put_att(sta%ncid, sta%station_data_id,
!     &             'positive', 'up')
!            CALL check_err(iret)
         sta%statnames => STATNAME

      case (62) !       F O R T . 6 2

         sta%station_data_dims(1) = sta%num_sta_dim_id
         sta%station_data_dims(2) = sta%myTime%timenc_dim_id
         iret = nf90_def_var(sta%ncid, 'u-vel', NF90_DOUBLE, &
                             sta%station_data_dims, sta%u_station_data_id)
         call check_err(iret)
         iret = nf90_def_var(sta%ncid, 'v-vel', NF90_DOUBLE, &
                             sta%station_data_dims, sta%v_station_data_id)
         call check_err(iret)
         if (ics /= 1) then
            ! u
            iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                                'long_name', &
                                'station water column vertically averaged east/west velocity')
            call check_err(iret)

            iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                                'standard_name', &
                                'station_eastward_water_velocity_depth_averaged')
            call check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                                'positive', 'east')
            call check_err(iret)
            ! v
            iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                                'long_name', &
                                'station water column vertically averaged north/south velocity')
            call check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                                'standard_name', &
                                'station_northward_water_velocity_depth_averaged')
            call check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                                'positive', 'north')
            call check_err(iret)
         else
            ! u
            iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                                'long_name', 'station water column vertically ' &
                                //'averaged velocity in x-direction')
            call check_err(iret)

            iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                                'standard_name', &
                                'station_x_water_velocity_depth_averaged')
            call check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                                'positive', 'right')
            call check_err(iret)
!              ! v
            iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                                'long_name', &
                                'station water column vertically averaged velocity ' &
                                //'in y-direction')
            call check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                                'standard_name', &
                                'station_y_water_velocity_depth_averaged')
            call check_err(iret)
            iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                                'positive', &
                                '90 degrees counterclockwise from x water velocity')
            call check_err(iret)
         end if

         call putUnitsAttribute(sta%ncid, sta%u_station_data_id, &
                                'm s-1', 'ft s-1')
         iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                             'dry_Value', doubleval)
         call check_err(iret)
         call putUnitsAttribute(sta%ncid, sta%v_station_data_id, &
                                'm s-1', 'ft s-1')
         iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                             'dry_Value', doubleval)
         call check_err(iret)
         sta%statnames => STATNAMEV

      case (71) !         F  O R T . 7 1

         sta%station_data_dims(1) = sta%num_sta_dim_id
         sta%station_data_dims(2) = sta%myTime%timenc_dim_id

         iret = nf90_def_var(sta%ncid, 'pressure', NF90_DOUBLE, &
                             sta%station_data_dims, sta%station_data_id)
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%station_data_id, &
                             'long_name', 'station air pressure at sea level')
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%station_data_id, &
                             'standard_name', 'station_air_pressure_at_sea_level')
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%station_data_id, &
                             'units', 'meters of water')
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%station_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%station_data_id, &
                             'positive', 'up')
         call check_err(iret)
         sta%statnames => STATNAMEM

      case (72) !      F O R T . 7 2

         sta%station_data_dims(1) = sta%num_sta_dim_id
         sta%station_data_dims(2) = sta%myTime%timenc_dim_id
         iret = nf90_def_var(sta%ncid, 'windx', NF90_DOUBLE, &
                             sta%station_data_dims, sta%u_station_data_id)
         call check_err(iret)
         iret = nf90_def_var(sta%ncid, 'windy', NF90_DOUBLE, &
                             sta%station_data_dims, sta%v_station_data_id)
         call check_err(iret)

         if (ics /= 1) then
            select case (abs(nws))
            case (1, 2)
               iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                                   'long_name', 'station wind stress u-component')
               call check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                                   'standard_name', &
                                   'station_eastward_surface_wind_stress')
               call check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                                   'positive', 'east')
               iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                                   'long_name', 'station wind stress v-component')
               call check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                                   'standard_name', &
                                   'station_northward_surface_wind_stress')
               call check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                                   'positive', 'north')
               call check_err(iret)
            case default
               iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                                   'long_name', 'station wind speed at 10m u-component')
               call check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                                   'standard_name', &
                                   'station_eastward_wind')
               call check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                                   'positive', 'east')
               iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                                   'long_name', &
                                   'station wind speed at 10m v-component')
               call check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                                   'standard_name', 'station_northward_wind')
               call check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                                   'positive', 'north')
            end select
         else
            select case (abs(nws))
            case (1, 2)
               iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                                   'long_name', &
                                   'station wind stress in x-direction')
               call check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                                   'standard_name', &
                                   'station_x_surface_wind_stress')
               call check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                                   'positive', 'right')
               iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                                   'long_name', &
                                   'station wind stress in y-direction')
               call check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                                   'standard_name', &
                                   'station_y_surface_wind_stress')
               call check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                                   'positive', &
                                   '90 degrees counterclockwise from wind velocity in x-direction')
            case default
               iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                                   'long_name', &
                                   'station wind velocity in x-direction')
               call check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                                   'standard_name', 'station_x_wind')
               call check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                                   'positive', 'right')
               iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                                   'long_name', &
                                   'station wind velocity in y-direction')
               call check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                                   'standard_name', &
                                   'station_y_wind')
               call check_err(iret)
               iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                                   'positive', &
                                   '90 degrees counterclockwise from wind velocity in x-direction')
               call check_err(iret)
            end select
         end if

         if (abs(nws) > 2 .and. abs(nws) < 100) then
            call putUnitsAttribute(sta%ncid, sta%u_station_data_id, &
                                   'm s-1', 'ft s-1')
            call putUnitsAttribute(sta%ncid, sta%v_station_data_id, &
                                   'm s-1', 'ft s-1')
         else
            call putUnitsAttribute(sta%ncid, sta%u_station_data_id, &
                                   'm2 s-2', 'ft s-2')
            call putUnitsAttribute(sta%ncid, sta%u_station_data_id, &
                                   'm2 s-2', 'ft s-2')
         end if

         iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         sta%statnames => STATNAMEM

      case (91) !         F  O R T . 9 1

         sta%station_data_dims(1) = sta%num_sta_dim_id
         sta%station_data_dims(2) = sta%myTime%timenc_dim_id

         iret = nf90_def_var(sta%ncid, 'iceaf', NF90_DOUBLE, &
                             sta%station_data_dims, sta%station_data_id)
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%station_data_id, &
                             'long_name', 'station sea-ice af')
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%station_data_id, &
                             'standard_name', 'station_sea-ice area fraction')
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%station_data_id, &
                             'units', 'unitless')
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%station_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%station_data_id, &
                             'positive', 'up')
         call check_err(iret)
         sta%statnames => STATNAMEM

      case (109) !       O F F S E T . 6 1

         sta%station_data_dims(1) = sta%num_sta_dim_id
         sta%station_data_dims(2) = sta%myTime%timenc_dim_id
         iret = nf90_def_var(sta%ncid, 'dynamicWaterlevelCorrection', NF90_DOUBLE, &
                             sta%station_data_dims, sta%station_data_id)
         call check_err(iret)
!        Define water surface elevation attributes
         iret = nf90_put_att(sta%ncid, sta%station_data_id, &
                             'long_name', &
                             'dynamic water surface correction above water level')
         call check_err(iret)
         iret = nf90_put_att(sta%ncid, sta%station_data_id, &
                             'standard_name', &
                             'dynamic_sea_surface_correction_above_water_level')
         call check_err(iret)
         call putUnitsAttribute(sta%ncid, sta%station_data_id, &
                                'm', 'ft')
         iret = nf90_put_att(sta%ncid, sta%station_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         sta%statnames => STATNAME

      case DEFAULT
         write (scratchMessage, &
                '("No netCDF for station files with unit number ",i0,".")') descript1%lun
         call allMessage(ERROR, scratchMessage)
      end select
      iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)
      iret = nf90_put_att(sta%ncid, sta%u_station_data_id, &
                          'dry_Value', doubleval)
      call check_err(iret)
      iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)
      iret = nf90_put_att(sta%ncid, sta%v_station_data_id, &
                          'dry_Value', doubleval)
      call check_err(iret)
      iret = nf90_put_att(sta%ncid, sta%w_station_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)
      iret = nf90_put_att(sta%ncid, sta%w_station_data_id, &
                          'dry_Value', doubleval)
      call check_err(iret)

!     jgf50.44: Automatically turn on compression if we are using the
!     netcdf4 file format.
#ifdef NETCDF_CAN_DEFLATE
      if (abs(descript1%specifier) == 5) then
         select case (descript1%lun)
         case (41)
            iret = nf90_def_var_deflate(sta%ncid, sta%u_station_data_id, &
                                        1, 1, 2)
            call check_err(iret)
            if ((IDEN == 2) .or. (IDEN == 4)) then
               iret = nf90_def_var_deflate(sta%ncid, &
                                           sta%v_station_data_id, 1, 1, 2)
               call check_err(iret)
            end if
            if ((IDEN == 3) .or. (IDEN == 4)) then
               iret = nf90_def_var_deflate(sta%ncid, &
                                           sta%w_station_data_id, 1, 1, 2)
               call check_err(iret)
            end if
         case (42, 43)
            iret = nf90_def_var_deflate(sta%ncid, sta%u_station_data_id, &
                                        1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(sta%ncid, &
                                        sta%v_station_data_id, 1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(sta%ncid, &
                                        sta%w_station_data_id, 1, 1, 2)
            call check_err(iret)
         case (51)
            call check_err(nf90_def_var_deflate(sta%ncid, sta%station_ha_data_id, &
                                                1, 1, 2))
            call check_err(nf90_def_var_deflate(sta%ncid, sta%station_hg_data_id, &
                                                1, 1, 2))
         case (52)
            call check_err(nf90_def_var_deflate(sta%ncid, sta%ha_u_station_data_id, &
                                                1, 1, 2))
            call check_err(nf90_def_var_deflate(sta%ncid, sta%hg_u_station_data_id, &
                                                1, 1, 2))
            call check_err(nf90_def_var_deflate(sta%ncid, sta%ha_v_station_data_id, &
                                                1, 1, 2))
            call check_err(nf90_def_var_deflate(sta%ncid, sta%hg_v_station_data_id, &
                                                1, 1, 2))
         case (61, 71, 91, 109) ! GML added 91 20210727
            iret = nf90_def_var_deflate(sta%ncid, sta%station_data_id, &
                                        1, 1, 2)
            call check_err(iret)
         case (62, 72)
            iret = nf90_def_var_deflate(sta%ncid, sta%u_station_data_id, &
                                        1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(sta%ncid, sta%v_station_data_id, &
                                        1, 1, 2)
            call check_err(iret)
         case DEFAULT ! should be unreachable
            write (scratchMessage, &
                   '("No netCDF for station files with unit number ",i0,".")') descript1%lun
            call allMessage(ERROR, scratchMessage)
         end select
      end if
#endif

!     Define station names and code attributes
      iret = nf90_put_att(sta%ncid, sta%station_id, 'long_name', &
                          'station name')
      call check_err(iret)

!     Define time attributes
      call defineTimeAttributes(sta%ncid, sta%myTime)

!     define metadata and selected fort.15 parameters in netcdf file
      call defineMetaData(sta%ncid)

!     Leave define mode
      iret = nf90_enddef(sta%ncid)
      call check_err(iret)

!     Store station name(s)
      iret = nf90_put_var(sta%ncid, sta%station_id, sta%statnames(:))
      call check_err(iret)

!     Store station locations
      iret = nf90_put_var(sta%ncid, sta%x_id, sta%x)
      call check_err(iret)
      iret = nf90_put_var(sta%ncid, sta%y_id, sta%y)
      call check_err(iret)

!     !WJP 02.20.2018 Store const name(s)
      if (any(stationLunsHA == descript1%lun)) then
         call addHarmonicAnalysisParametersToNetcdfFile(sta%ncid)
      end if

!     now close the initialized netcdf file
      iret = nf90_close(sta%ncid)
      call check_err(iret)
!-----------------------------------------------------------------------
   end subroutine initStationFile
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E    I N I T  N O D A L  D A T A   F I L E
!-----------------------------------------------------------------------
!     jgf49.17.02 Sets up netCDF variables and allocates memory for
!     full domain data.
!-----------------------------------------------------------------------
   subroutine initNodalDataFile(dat, descript1, reterr)
      use SIZES, only: MNWPROC, MNPROC, MYPROC
      use GLOBAL, only: OutputDataDescript_t, NWS, C3D, IDEN, NE_G, NP_G
      use GLOBAL_3DVS, only: NFEN
      use MESH, only: ICS
      implicit none

      type(nodalData), intent(inout) :: dat
      type(OutputDataDescript_t), intent(inout) :: descript1
      logical, intent(out) :: reterr

      integer :: iret ! success or failure of netcdf call
      integer, parameter :: nodalLunsHA(2) = [53, 54]
      real(8), allocatable :: defaultValue(:)
      character(len=1024) :: att_text ! metadata
      character(len=1024) :: scratchMessage

      LOG_SCOPE_TRACED("initNodalDataFile", NETCDFIO_TRACING)
      reterr = .false.
      dat%myMesh => adcircMesh
      if (dat%myMesh%initialized .eqv. .false.) then
         dat%myMesh%num_nodes = NP_G
         dat%myMesh%num_elems = NE_G
         if (C3D .eqv. .true.) then
            dat%myMesh%num_v_nodes = NFEN
         end if
         dat%myMesh%nface_len = 3
      end if
      allocate (dat%myTime%timenc(dat%myTime%timenc_len))
      allocate (defaultValue(dat%myMesh%num_nodes))
      defaultValue(:) = descript1%initial_value
      dat%myTime%initialized = .true.

!     Initialize netCDF output file, creating a new one if necessary.
      call createNetCDFOutputFile(dat%ncid, dat%myFile, dat%myTime, &
                                  descript1, reterr)

      ! if we didn't need to create a file, update metadata and return
      if (dat%myFile%createFile .eqv. .false.) then
         if ((reterr .eqv. .false.) .and. &
             ((abs(descript1%specifier) == NETCDF3) .or. &
              (abs(descript1%specifier) == NETCDF4))) then

            ! This subroutine is called in the case where
            ! a non-NetCDF min/max file will be written, but a NetCDF min/max
            ! file will be read. If this is the case, then don't update metadata
            ! in the netcdf file.
            ! The writers don't read the fort.15, so they need to
            ! assume ADCPREP has placed the correct info on the first
            ! pass, then call updateMetaData.
            if ((MNWPROC > 0) .and. &
                (writerReadMetaData .eqv. .false.) .and. &
                (MYPROC /= 0)) then
               call ReadMetaData(dat%ncid, dat%myFile)
               writerReadMetaData = .true.
               dat%myMesh%num_nodes = NP_G
               dat%myMesh%num_elems = NE_G
            end if
            ! jgf52.08.22: Only update meta data in serial, in parallel
            ! this was updated by adcprep
            if (mnproc == 1) then
               call updateMetaData(dat%ncid, dat%myFile)
            end if
         end if
         return
      end if

      if (dat%myMesh%initialized .eqv. .false.) then
         call initNetCDFCoord(dat%myMesh)
      end if
      call defineMeshVariables(dat%ncid, dat%myMesh, dat%myFile)

!     !WJP 02.20.2018 Writing in constituent information for harmonic file
      if (any(nodalLunsHA == descript1%lun)) then
         call defineHarmonicAnalysisParametersinNetcdfFile(descript1, nodaldat=dat)
      end if

      select case (descript1%lun)

      case (44)
         dat%nodal_data_dims_3D(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims_3D(2) = dat%myMesh%num_v_nodes_dim_id
         dat%nodal_data_dims_3D(3) = dat%myTime%timenc_dim_id
         iret = nf90_def_var(dat%ncid, 'sigmat', NF90_DOUBLE, &
                             dat%nodal_data_dims_3D, dat%u_nodal_data_id)
         call check_err(iret)
         if ((IDEN == 2) .or. (IDEN == 4)) then
            iret = nf90_def_var(dat%ncid, 'salinity', NF90_DOUBLE, &
                                dat%nodal_data_dims_3D, dat%v_nodal_data_id)
            call check_err(iret)
         end if
         if ((IDEN == 3) .or. (IDEN == 4)) then
            iret = nf90_def_var(dat%ncid, 'temperature', NF90_DOUBLE, &
                                dat%nodal_data_dims_3D, dat%w_nodal_data_id)
            call check_err(iret)
         end if
         ! sigma t
         att_text = "water column vertically varying density"
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             'long_name', trim(att_text))
         call check_err(iret)
         att_text = "water_density_vertically_varying"
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             'standard_name', trim(att_text))
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%u_nodal_data_id, &
                                "kg/m^3", "n/a")
         ! salinity
         if ((IDEN == 2) .or. (IDEN == 4)) then
            att_text = "water column vertically varying salinity"
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                'long_name', trim(att_text))
            call check_err(iret)
            att_text = "water_salinity_vertically_varying"
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                'standard_name', trim(att_text))
            call check_err(iret)
            call putUnitsAttribute(dat%ncid, dat%v_nodal_data_id, &
                                   "PSU", "PSU")
         end if
         ! temperature
         if ((IDEN == 3) .or. (IDEN == 4)) then
            att_text = "water column vertically varying temperature"
            iret = nf90_put_att(dat%ncid, dat%w_nodal_data_id, &
                                'long_name', trim(att_text))
            call check_err(iret)
            att_text = "water_temperature_vertically_varying"
            iret = nf90_put_att(dat%ncid, dat%w_nodal_data_id, &
                                'standard_name', trim(att_text))
            call check_err(iret)
            call putUnitsAttribute(dat%ncid, dat%w_nodal_data_id, &
                                   "Celsius", "Fahrenheit")
         end if
      case (45)
         dat%nodal_data_dims_3D(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims_3D(2) = dat%myMesh%num_v_nodes_dim_id
         dat%nodal_data_dims_3D(3) = dat%myTime%timenc_dim_id
         iret = nf90_def_var(dat%ncid, 'u-vel3D', NF90_DOUBLE, &
                             dat%nodal_data_dims_3D, dat%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_def_var(dat%ncid, 'v-vel3D', NF90_DOUBLE, &
                             dat%nodal_data_dims_3D, dat%v_nodal_data_id)
         call check_err(iret)
         iret = nf90_def_var(dat%ncid, 'w-vel3D', NF90_DOUBLE, &
                             dat%nodal_data_dims_3D, dat%w_nodal_data_id)
         call check_err(iret)
         if (ics /= 1) then
            ! u
            att_text = &
               'water column vertically varying east/west velocity'
            iret = nf90_put_att(dat%ncid, &
                                dat%u_nodal_data_id, 'long_name', trim(att_text))
            call check_err(iret)
            att_text = 'eastward_water_velocity_vertically_varying'
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                'standard_name', trim(att_text))
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, &
                                dat%u_nodal_data_id, 'positive', 'east')
            call check_err(iret)
            ! v
            att_text = &
               'water column vertically varying north/south velocity'
            iret = nf90_put_att(dat%ncid, &
                                dat%v_nodal_data_id, 'long_name', trim(att_text))
            call check_err(iret)
            att_text = 'northward_water_velocity_vertically_varying'
            iret = nf90_put_att(dat%ncid, &
                                dat%v_nodal_data_id, 'standard_name', trim(att_text))
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, &
                                dat%v_nodal_data_id, 'positive', 'north')
            call check_err(iret)
            ! w
            att_text = &
               'water column vertically varying up/down velocity'
            iret = nf90_put_att(dat%ncid, dat%w_nodal_data_id, &
                                'long_name', trim(att_text))
            call check_err(iret)
            att_text = &
               'upward_water_velocity_vertically_varying'
            iret = nf90_put_att(dat%ncid, dat%w_nodal_data_id, &
                                'standard_name', trim(att_text))
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, &
                                dat%w_nodal_data_id, 'positive', 'up')
            call check_err(iret)
         else
            ! u
            att_text &
               = 'water column vertically varying velocity in x-direction'
            iret = nf90_put_att(dat%ncid, &
                                dat%u_nodal_data_id, 'long_name', trim(att_text))
            call check_err(iret)
            att_text = 'x_water_velocity_vertically_varying'
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                'standard_name', trim(att_text))
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, &
                                dat%u_nodal_data_id, 'positive', 'right')
            call check_err(iret)
            ! v
            att_text = &
               'water column vertically varying velocity in y-direction'
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                'long_name', trim(att_text))
            call check_err(iret)
            att_text = 'station_y_water_velocity_vertically_varying'
            iret = nf90_put_att(dat%ncid, &
                                dat%v_nodal_data_id, 'standard_name', trim(att_text))
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                'positive', '90 degrees counterclockwise from x'// &
                                'water velocity')
            call check_err(iret)
            ! w
            att_text = 'water column vertically ' &
                       //'varying velocity in z-direction'
            iret = nf90_put_att(dat%ncid, &
                                dat%w_nodal_data_id, 'long_name', trim(att_text))
            call check_err(iret)
            att_text = 'station_z_water_velocity_vertically_varying'
            iret = nf90_put_att(dat%ncid, &
                                dat%w_nodal_data_id, 'standard_name', trim(att_text))
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, &
                                dat%w_nodal_data_id, 'positive', &
                                '90 degrees counterclockwise from x water velocity')
            call check_err(iret)
         end if
         call putUnitsAttribute(dat%ncid, dat%u_nodal_data_id, &
                                'm s-1', 'ft s-1')
         call putUnitsAttribute(dat%ncid, dat%v_nodal_data_id, &
                                'm s-1', 'ft s-1')
         call putUnitsAttribute(dat%ncid, dat%w_nodal_data_id, &
                                'm s-1', 'ft s-1')

      case (46) !       F O R T . 4 6
         dat%nodal_data_dims_3D(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims_3D(2) = dat%myMesh%num_v_nodes_dim_id
         dat%nodal_data_dims_3D(3) = dat%myTime%timenc_dim_id
         iret = nf90_def_var(dat%ncid, 'q20', NF90_DOUBLE, &
                             dat%nodal_data_dims_3D, dat%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_def_var(dat%ncid, 'l', NF90_DOUBLE, &
                             dat%nodal_data_dims_3D, dat%v_nodal_data_id)
         call check_err(iret)
         iret = nf90_def_var(dat%ncid, 'ev', NF90_DOUBLE, &
                             dat%nodal_data_dims_3D, dat%w_nodal_data_id)
         call check_err(iret)
         ! q20
         att_text = &
            'water column vertically varying turbulent kinetic energy'
         iret = nf90_put_att(dat%ncid, &
                             dat%u_nodal_data_id, 'long_name', trim(att_text))
         call check_err(iret)
         att_text = &
            'turbulent_kinetic_energy_vertically_varying'
         iret = nf90_put_att(dat%ncid, &
                             dat%u_nodal_data_id, 'standard_name', trim(att_text))
         call check_err(iret)
         ! l
         att_text = &
            'water column vertically varying mixing length'
         iret = nf90_put_att(dat%ncid, &
                             dat%v_nodal_data_id, 'long_name', trim(att_text))
         call check_err(iret)
         att_text = &
            'water_mixing_length_vertically_varying'
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             'standard_name', trim(att_text))
         call check_err(iret)
         !  ev
         att_text = &
            'water column vertically varying eddy viscosity'
         iret = nf90_put_att(dat%ncid, &
                             dat%w_nodal_data_id, 'long_name', trim(att_text))
         call check_err(iret)
         att_text = 'water_eddy_viscosity_vertically_varying'
         iret = nf90_put_att(dat%ncid, &
                             dat%w_nodal_data_id, 'standard_name', trim(att_text))
         call check_err(iret)
!            CALL putUnitsAttribute(dat%ncid,dat%u_nodal_data_id,
!     &      'm s-2', 'ft s-1') !TODO: jgf49.48.01: units for q20??
         call putUnitsAttribute(dat%ncid, dat%v_nodal_data_id, &
                                "meters", "n/a")
         call putUnitsAttribute(dat%ncid, dat%w_nodal_data_id, &
                                "m^2/s", "n/a")

      case (47)

         dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
         dat%varnames(1) = 'qsurfkp1'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
         call check_err(iret)
         att_text = "sea surface temperature at the k+1 time level"
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'long_name', trim(att_text))
         call check_err(iret)
         att_text = "future sea surface temperature"
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'standard_name', trim(att_text))
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%nodal_data_id, &
                                "Celsius", "Fahrenheit")
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)

      case (48)
         dat%nodal_data_dims_3D(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims_3D(2) = dat%myMesh%num_v_nodes_dim_id
         dat%nodal_data_dims_3D(3) = dat%myTime%timenc_dim_id
         iret = nf90_def_var(dat%ncid, 'bpgx', NF90_DOUBLE, &
                             dat%nodal_data_dims_3D, dat%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_def_var(dat%ncid, 'bpgy', NF90_DOUBLE, &
                             dat%nodal_data_dims_3D, dat%v_nodal_data_id)
         call check_err(iret)
         ! bpgx
         att_text = "u-component of baroclinic pressure gradient"
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             'long_name', trim(att_text))
         call check_err(iret)
         att_text = "u-component of baroclinic pressure gradient"
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             'standard_name', trim(att_text))
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%u_nodal_data_id, &
                                "n/a", "n/a")
         ! bpgy
         att_text = "v-component of baroclinic pressure gradient"
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             'long_name', trim(att_text))
         call check_err(iret)
         att_text = "v-component of baroclinic pressure gradient"
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             'standard_name', trim(att_text))
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%v_nodal_data_id, &
                                "n/a", "n/a")

      case (53)
         dat%nodal_data_dims(2) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(1) = dat%myMesh%num_v_nodes_dim_id
!        Amplitude
         dat%varnames(1) = 'amp'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%ha_nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%ha_nodal_data_id, &
                             'long_name', &
                             'amplitude of tidal harmonic constituent elevation')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%ha_nodal_data_id, &
                             'standard_name', 'hc_elev_amp')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%ha_nodal_data_id, &
                             'coordinates', 'const y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%ha_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%ha_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%ha_nodal_data_id, &
                                'm', 'ft')
         iret = nf90_put_att(dat%ncid, dat%ha_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
!        Phase
         dat%varnames(1) = 'phs'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%hg_nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%hg_nodal_data_id, &
                             'long_name', &
                             'phase of tidal harmonic constituent elevation')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%hg_nodal_data_id, &
                             'standard_name', 'hc_elev_phs')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%hg_nodal_data_id, &
                             'coordinates', 'const y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%hg_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%hg_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%hg_nodal_data_id, &
                                'deg', 'deg')
         iret = nf90_put_att(dat%ncid, dat%hg_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)

      case (54)
         dat%nodal_data_dims(2) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(1) = dat%myMesh%num_v_nodes_dim_id
!        u Amplitude
         dat%varnames(1) = 'u_amp'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%u_ha_nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_ha_nodal_data_id, &
                             'long_name', &
                             'amplitude of tidal harmonic constituent E-W velocity')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_ha_nodal_data_id, &
                             'standard_name', 'hc_u_amp')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_ha_nodal_data_id, &
                             'coordinates', 'const y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_ha_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_ha_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%u_ha_nodal_data_id, &
                                'm/s', 'ft/s')
         iret = nf90_put_att(dat%ncid, dat%u_ha_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
!        u Phase
         dat%varnames(1) = 'u_phs'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%u_hg_nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_hg_nodal_data_id, &
                             'long_name', &
                             'phase of tidal harmonic constituent E-W velocity')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_hg_nodal_data_id, &
                             'standard_name', 'hc_u_phs')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_hg_nodal_data_id, &
                             'coordinates', 'const y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_hg_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_hg_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%u_hg_nodal_data_id, &
                                'deg', 'deg')
         iret = nf90_put_att(dat%ncid, dat%u_hg_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
!        v amplitude
         dat%varnames(1) = 'v_amp'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%v_ha_nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_ha_nodal_data_id, &
                             'long_name', &
                             'amplitude of tidal harmonic constituent N-S velocity')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_ha_nodal_data_id, &
                             'standard_name', 'hc_v_amp')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_ha_nodal_data_id, &
                             'coordinates', 'const y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_ha_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_ha_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%v_ha_nodal_data_id, &
                                'm/s', 'ft/s')
         iret = nf90_put_att(dat%ncid, dat%v_ha_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
!        v Phase
         dat%varnames(1) = 'v_phs'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%v_hg_nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_hg_nodal_data_id, &
                             'long_name', &
                             'phase of tidal harmonic constituent N-S velocity')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_hg_nodal_data_id, &
                             'standard_name', 'hc_v_phs')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_hg_nodal_data_id, &
                             'coordinates', 'const y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_hg_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_hg_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%v_hg_nodal_data_id, &
                                'deg', 'deg')
         iret = nf90_put_att(dat%ncid, dat%v_hg_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)

      case (63)
         dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
         dat%varnames(1) = 'zeta'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'long_name', 'water surface elevation above geoid')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'standard_name', 'sea_surface_height_above_geoid')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%nodal_data_id, &
                                'm', 'ft')
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
      case (90)
         dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
         dat%varnames(1) = 'tau0'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'long_name', 'primitive weighting in continuity equation')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'standard_name', 'primitive_weighting_in_continuity_equation')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%nodal_data_id, &
                                '1', '1')
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
      case (311)
         dat%varnames(1) = 'zeta_max'
         dat%varnames(2) = 'time_of_zeta_max'
         dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(2)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, dat%time_max_nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'long_name', 'maximum water surface elevation'// &
                             'above geoid')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'long_name', 'time of maximum water surface elevation'// &
                             'above geoid')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'standard_name', 'maximum_sea_surface_height_'// &
                             'above_geoid')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'standard_name', 'time_of_maximum_sea_surface_height_'// &
                             'above_geoid')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id, &
                                'm', 'ft')
         call putUnitsAttribute(dat%ncid, dat%time_max_nodal_data_id, &
                                'sec', 'sec')
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)

      case (64)
         dat%varnames(1) = 'u-vel'
         dat%varnames(2) = 'v-vel'
         dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(2)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%v_nodal_data_id)
         call check_err(iret)
         if (ics /= 1) then
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                'long_name', &
                                'water column vertically averaged east/west velocity')
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                'standard_name', 'eastward_water_velocity')
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                'positive', 'east')
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                'long_name', &
                                'water column vertically averaged north/south velocity')
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                'standard_name', 'northward_water_velocity')
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                'positive', 'north')
            call check_err(iret)
         else
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                'long_name', &
                                'water column vertically averaged velocity in x-direction')
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                'standard_name', &
                                'x_water_velocity_depth_averaged')
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                'positive', 'right')
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                'long_name', &
                                'water column vertically averaged velocity in y-direction')
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                'standard_name', &
                                'y_water_velocity_depth_averaged')
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                'positive', &
                                '90 degrees counterclockwise from x water velocity')
            call check_err(iret)
         end if
         call putUnitsAttribute(dat%ncid, dat%u_nodal_data_id, &
                                'm s-1', 'ft s-1')
         call putUnitsAttribute(dat%ncid, dat%v_nodal_data_id, &
                                'm s-1', 'ft s-1')
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             'dry_Value', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             'dry_Value', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (312)
         dat%varnames(1) = 'vel_max'
         dat%varnames(2) = 'time_of_vel_max'
         dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, &
                             dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(2)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, &
                             dat%time_max_nodal_data_id)
         call check_err(iret)

         if (ics /= 1) then
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                'long_name', &
                                'maximum water column vertically averaged velocity')
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                                'long_name', &
                                'time of maximum water column vertically averaged velocity')
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                'standard_name', 'maximum_water_velocity')
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                                'standard_name', 'time_of_maximum_water_velocity')
            call check_err(iret)
         else
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                'long_name', &
                                'maximum water column vertically averaged'// &
                                'velocity')
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                                'long_name', &
                                'time of maximum water column vertically averaged'// &
                                'velocity')
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                'standard_name', &
                                'maximum_water_velocity_depth_averaged')
            call check_err(iret)
            iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                                'standard_name', &
                                'time_of_maximum_water_velocity_depth_averaged')
            call check_err(iret)
         end if
         call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id, &
                                'm s-1', 'ft s-1')
         call putUnitsAttribute(dat%ncid, dat%time_max_nodal_data_id, &
                                'sec', 'sec')
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'dry_Value', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'dry_Value', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (173)
         dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
         dat%varnames(1) = 'winddrag'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'long_name', 'wind drag coefficient at sea level')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'standard_name', 'wind drag coefficient')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%nodal_data_id, &
                                'unitless', 'unitless')
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)

      case (73)
         dat%varnames(1) = 'pressure'
         dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
         iret = nf90_def_var(dat%ncid, &
                             trim(dat%varnames(1)), NF90_DOUBLE, &
                             dat%nodal_data_dims, &
                             dat%nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'long_name', 'air pressure at sea level')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'standard_name', 'air_pressure_at_sea_level')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, 'units', &
                             'meters of water')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (313)
         dat%varnames(1) = 'pressure_min'
         dat%varnames(2) = 'time_of_pressure_min'
         dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
         iret = nf90_def_var(dat%ncid, &
                             trim(dat%varnames(1)), NF90_DOUBLE, &
                             dat%nodal_data_dims_max, &
                             dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_def_var(dat%ncid, &
                             trim(dat%varnames(2)), NF90_DOUBLE, &
                             dat%nodal_data_dims_max, &
                             dat%time_max_nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'long_name', 'minimum air pressure at sea level')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'long_name', 'time of minimum air pressure at sea level')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'standard_name', 'minimum_air_pressure_at_sea_level')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'standard_name', &
                             'time_of_minimum_air_pressure_at_sea_level')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, 'units', &
                             'meters of water')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, 'units', &
                             'seconds')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (74)
         dat%varnames(1) = 'windx'
         dat%varnames(2) = 'windy'
         dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(2)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%v_nodal_data_id)
         call check_err(iret)
         if (ics /= 1) then
            select case (abs(nws))
            case (1, 2)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'long_name', &
                                   'wind stress u-component')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'standard_name', &
                                   'eastward_surface_wind_stress')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'positive', 'east')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'long_name', &
                                   'wind stress v-component')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'standard_name', &
                                   'northward_surface_wind_stress')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'positive', 'north')
               call check_err(iret)
            case default
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'long_name', 'wind speed at 10m u-component')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'standard_name', 'eastward_wind')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'positive', 'east')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'long_name', 'wind speed at 10m v-component')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'standard_name', 'northward_wind')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'positive', 'north')
            end select
         else
            select case (abs(nws))
            case (1, 2)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'long_name', 'wind stress in x-direction')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'standard_name', 'x_surface_wind_stress')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'positive', 'right')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'long_name', 'wind stress in y-direction')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'standard_name', 'y_surface_wind_stress')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'positive', &
                                   '90 degrees counterclockwise from wind velocity in x-direction')
               call check_err(iret)
            case default
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'long_name', 'wind velocity in x-direction')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'standard_name', 'x_wind')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'positive', 'right')

               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'long_name', 'wind velocity in y-direction')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'standard_name', 'y_wind')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'positive', &
                                   '90 degrees counterclockwise from wind velocity in x-direction')
               call check_err(iret)
            end select
         end if
         select case (abs(nws))
         case (1, 2)
            call putUnitsAttribute(dat%ncid, dat%u_nodal_data_id, &
                                   'm2 s-2', 'lbf s-2')
            call putUnitsAttribute(dat%ncid, dat%v_nodal_data_id, &
                                   'm2 s-2', 'lbf s-2')
         case default
            call putUnitsAttribute(dat%ncid, dat%u_nodal_data_id, &
                                   'm s-1', 'ft s-1')
            call putUnitsAttribute(dat%ncid, dat%v_nodal_data_id, &
                                   'm s-1', 'ft s-1')
         end select
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (314)
         dat%varnames(1) = 'wind_max'
         dat%varnames(2) = 'time_of_wind_max'
         dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, &
                             dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(2)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, &
                             dat%time_max_nodal_data_id)
         call check_err(iret)
         if (ics /= 1) then
            select case (abs(nws))
            case (1, 2)
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                   'long_name', &
                                   'maximum wind stress')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                                   'long_name', &
                                   'time of maximum wind stress')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                   'standard_name', &
                                   'maximum_surface_wind_stress')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                                   'standard_name', &
                                   'time_of_maximum_surface_wind_stress')
               call check_err(iret)
            case default
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                   'long_name', 'maximum wind velocity')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                                   'long_name', 'time of maximum wind velocity')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                   'standard_name', 'maximum_wind')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                                   'standard_name', 'time_of_maximum_wind')
               call check_err(iret)
            end select
         else
            select case (abs(nws))
            case (1, 2)
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                   'long_name', 'maximum wind stress')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                                   'long_name', 'time of maximum wind stress')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                   'standard_name', 'maximum_surface_wind_stress')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                                   'standard_name', &
                                   'time_of_maximum_surface_wind_stress')
               call check_err(iret)
            case default
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                   'long_name', 'maximum wind velocity')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                                   'long_name', 'time of maximum wind velocity')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                   'standard_name', 'maximum_wind')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                                   'standard_name', 'time_of_maximum_wind')
               call check_err(iret)
            end select
         end if
         select case (abs(nws))
         case (1, 2)
            call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id, &
                                   'm2 s-2', 'lbf ft-2')
            call putUnitsAttribute(dat%ncid, dat%time_max_nodal_data_id, &
                                   'sec', 'sec')
         case default
            call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id, &
                                   'm s-1', 'ft s-1')
            call putUnitsAttribute(dat%ncid, dat%time_max_nodal_data_id, &
                                   'sec', 'sec')
         end select
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (315)
         dat%varnames(1) = 'radstress_max'
         dat%varnames(2) = 'time_of_radstress_max'
         dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, &
                             dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(2)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, &
                             dat%time_max_nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'long_name', 'maximum radiation stress gradient')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'long_name', &
                             'time of maximum radiation stress gradient')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'standard_name', 'maximum_radiation_stress')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'standard_name', 'time_of_maximum_radiation_stress')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id, &
                                'm-2 s-2', 'lbf ft-2')
         call putUnitsAttribute(dat%ncid, dat%time_max_nodal_data_id, &
                                'sec', 'sec')
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (77)
         dat%varnames(1) = 'weir_dz'
         dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, dat%nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'long_name', &
                             'elevation change in weir boundary condition')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'standard_name', 'weir_elevation_change')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%nodal_data_id, &
                                'm', 'ft')
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             '_FillValue', 0.0d0)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (93)
         dat%varnames(1) = 'iceaf'
         dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id

         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
         call check_err(iret)

         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'long_name', 'sea-ice af')
         call check_err(iret)

         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'standard_name', 'sea-ice area fraction')
         call check_err(iret)

         call putUnitsAttribute(dat%ncid, dat%nodal_data_id, &
                                'unitless', 'unitless')

         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)

         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)

         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'location', 'node')
         call check_err(iret)

         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (164)
         dat%varnames(1) = 'radstress_x'
         dat%varnames(2) = 'radstress_y'
         dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(2)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%v_nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             'long_name', 'radiation stress gradient x component')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             'long_name', 'radiation stress gradient y component')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             'standard_name', 'radiation_stress_gradient_x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             'standard_name', 'radiation_stress_gradient_y')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%u_nodal_data_id, &
                                'm-2 s-2', 'ft-2 s-2')
         call putUnitsAttribute(dat%ncid, dat%v_nodal_data_id, &
                                'm-2 s-2', 'ft-2 s-2')
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (301)
         dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
         dat%varnames(1) = 'swan_HS'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'long_name', 'significant wave height')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'standard_name', &
                             'sea_surface_wave_significant_height')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%nodal_data_id, &
                                'm', 'ft')
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (316)
         dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
         dat%varnames(1) = 'swan_HS_max'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'long_name', 'maximum significant wave height')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'standard_name', &
                             'maximum_sea_surface_wave_significant_height')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id, &
                                'm', 'ft')
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (302)
         dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
         dat%varnames(1) = 'swan_DIR'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'long_name', 'mean wave direction')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'standard_name', 'sea_surface_wave_to_direction')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%nodal_data_id, &
                                'degrees', 'degrees_CW_from_East')
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (317)
         dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
         dat%varnames(1) = 'swan_DIR_max'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, &
                             dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'long_name', 'maximum mean wave direction')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'standard_name', 'maximum_sea_surface_wave_to_direction')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id, &
                                'degrees', 'degrees_CW_from_East')
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (303)
         dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
         dat%varnames(1) = 'swan_TM01'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'long_name', 'mean absolute wave period')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'standard_name', &
                             'sea_surface_wave_mean_period_from_variance_spectral'// &
                             '_density_first_frequency_moment')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%nodal_data_id, &
                                's', 's')
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (318)
         dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
         dat%varnames(1) = 'swan_TM01_max'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, &
                             dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'long_name', 'maximum TM01 mean wave period')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'standard_name', &
                             'maximum_sea_surface_wave_mean_period_from_variance_'// &
                             'spectral_density_first_frequency_moment')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id, &
                                's', 's')
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (304)
         dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
         dat%varnames(1) = 'swan_TPS'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'long_name', 'smoothed peak period')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'standard_name', 'sea_surface_wave_period_'// &
                             'at_variance_spectral_density_maximum')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%nodal_data_id, &
                                's', 's')
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (319)
         dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
         dat%varnames(1) = 'swan_TPS_max'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'long_name', 'maximum smoothed peak period')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'standard_name', 'maximum_sea_surface_wave_period'// &
                             '_at_variance_spectral_density_maximum')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id, &
                                's', 's')
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (305)
         dat%varnames(1) = 'swan_windx'
         dat%varnames(2) = 'swan_windy'
         dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(2)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%v_nodal_data_id)
         call check_err(iret)
         if (ics /= 1) then
            select case (abs(nws))
            case (1, 2)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'long_name', &
                                   'wind stress u-component')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'standard_name', &
                                   'eastward_surface_wind_stress')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'positive', 'east')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'long_name', &
                                   'wind stress v-component')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'standard_name', &
                                   'northward_surface_wind_stress')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'positive', 'north')
               call check_err(iret)
            case default
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'long_name', 'wind speed at 10m u-component')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'standard_name', 'eastward_wind')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'positive', 'east')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'long_name', 'wind speed at 10m v-component')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'standard_name', 'northward_wind')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'positive', 'north')
            end select
         else
            select case (abs(nws))
            case (1, 2)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'long_name', 'wind stress in x-direction')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'standard_name', 'x_surface_wind_stress')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'positive', 'right')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'long_name', 'wind stress in y-direction')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'standard_name', 'y_surface_wind_stress')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'positive', &
                                   '90 degrees counterclockwise from wind velocity in x-direction')
            case default
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'long_name', 'wind velocity in x-direction')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'standard_name', 'x_wind')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                                   'positive', 'right')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'long_name', 'wind velocity in y-direction')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'standard_name', 'y_wind')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                                   'positive', &
                                   '90 degrees counterclockwise from wind velocity in x-direction')
               call check_err(iret)
            end select
         end if

         select case (abs(nws))
         case (1, 2)
            call putUnitsAttribute(dat%ncid, dat%u_nodal_data_id, &
                                   'm2 s-2', 'ft s-2')
            call putUnitsAttribute(dat%ncid, dat%v_nodal_data_id, &
                                   'm2 s-2', 'ft s-2')
         case default
            call putUnitsAttribute(dat%ncid, dat%u_nodal_data_id, &
                                   'm s-1', 'ft s-1')
            call putUnitsAttribute(dat%ncid, dat%v_nodal_data_id, &
                                   'm s-1', 'ft s-1')
         end select
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (320)
         dat%varnames(1) = 'swan_wind_max'
         dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, dat%max_nodal_data_id)
         call check_err(iret)
         if (ics /= 1) then
            select case (abs(nws))
            case (1, 2)
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                   'long_name', &
                                   'maximum wind stress')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                   'standard_name', &
                                   'maximum_surface_wind_stress')
               call check_err(iret)
            case default
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                   'long_name', 'maximum wind velocity')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                   'standard_name', 'maximum_wind')
               call check_err(iret)
            end select
         else
            select case (abs(nws))
            case (1, 2)
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                   'long_name', 'maximum wind stress')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                   'standard_name', 'maximum_surface_wind_stress')
               call check_err(iret)
            case default
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                   'long_name', 'maximum wind velocity')
               call check_err(iret)
               iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                                   'standard_name', 'maximum_wind')
               call check_err(iret)
            end select
         end if
         select case (abs(nws))
         case (1, 2)
            call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id, &
                                   'm2 s-2', 'ft s-2')
         case default
            call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id, &
                                   'm s-1', 'ft s-1')
         end select
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (306)
         dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
         dat%varnames(1) = 'swan_TM02'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'long_name', 'mean absoloute zero crossing period')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'standard_name', &
                             'sea_surface_wave_mean_period_from_variance_spectral'// &
                             '_density_second_frequency_moment')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%nodal_data_id, &
                                's', 's')
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (321)
         dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
         dat%varnames(1) = 'swan_TM02_max'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'long_name', 'maximum TM02 mean wave period')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'standard_name', &
                             'maximum_sea_surface_wave_mean_period_from_variance_'// &
                             'spectral_density_second_frequency_moment')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id, &
                                's', 's')
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (307)
         dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
         dat%varnames(1) = 'swan_TMM10'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'long_name', 'mean absolute wave period')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'standard_name', &
                             'sea_surface_wave_mean_period_from_variance_spectral'// &
                             '_density_inverse_frequency_moment')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%nodal_data_id, &
                                's', 's')
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (322)
         dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
         dat%varnames(1) = 'swan_TMM10_max'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'long_name', 'maximum TMM10 mean wave period')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'standard_name', &
                             'maximum_sea_surface_wave_mean_period_from_variance_'// &
                             'spectral_density_inverse_frequency_moment')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id, &
                                's', 's')
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

      case (400)
         dat%varnames(1) = 'inun_time'
         dat%varnames(2) = 'onset_inun_time'
         dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(2)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, dat%time_max_nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'long_name', 'total time of inundation beyond the threshold')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'long_name', &
                             'time of onset of inundation beyond the threshold')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'standard_name', &
                             'total_time_of_inundation_beyond_the_threshold')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'standard_name', &
                             'time_of_onset_of_inundation_beyond_the_threshold')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id, &
                                'm', 'ft')
         call putUnitsAttribute(dat%ncid, dat%time_max_nodal_data_id, &
                                'sec', 'sec')
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)

      case (401)
         dat%varnames(1) = 'inun_max'
         dat%varnames(2) = 'time_of_inun_max'
         dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(2)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, dat%time_max_nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'long_name', 'maximum water inundation depth'// &
                             'above geoid')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'long_name', 'time of maximum water inundation depth'// &
                             'above geoid')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'standard_name', 'maximum_water_inundation_depth')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'standard_name', &
                             'time_of_maximum_water_inundation_depth')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id, &
                                'm', 'ft')
         call putUnitsAttribute(dat%ncid, dat%time_max_nodal_data_id, &
                                'sec', 'sec')
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)

      case (402)
         dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
         dat%varnames(1) = 'initiallydry'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_INT, dat%nodal_data_dims, dat%nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'long_name', &
                             'dry nodes at cold start')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'standard_name', 'initially_dry')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)

         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             '_FillValue', intFillValue)
         call check_err(iret)

      case (403)
         dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
         dat%varnames(1) = 'endrisinginun'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_INT, dat%nodal_data_dims, dat%nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'long_name', &
                             'rising inundation at the end of the simulation')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'standard_name', 'end_rising_inundation')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             '_FillValue', intFillValue)
         call check_err(iret)

      case (404)
         dat%varnames(1) = 'everdried'
         dat%varnames(2) = 'time_of_everdried'
         dat%nodal_data_dims_max(1) = dat%myMesh%num_nodes_dim_id
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(2)), &
                             NF90_DOUBLE, dat%nodal_data_dims_max, dat%time_max_nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'long_name', 'ever dried')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'long_name', &
                             'time of most recent drying occurrence')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'standard_name', 'ever_dried')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'standard_name', &
                             'time_of_most_recent_drying_occurrence')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'coordinates', 'y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%max_nodal_data_id, &
                                'unitless', 'unitless')
         call putUnitsAttribute(dat%ncid, dat%time_max_nodal_data_id, &
                                'sec', 'sec')
         iret = nf90_put_att(dat%ncid, dat%max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%time_max_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)

      case (108)
         dat%nodal_data_dims(1) = dat%myMesh%num_nodes_dim_id
         dat%nodal_data_dims(2) = dat%myTime%timenc_dim_id
         dat%varnames(1) = 'dynamicWaterlevelCorrection'
         iret = nf90_def_var(dat%ncid, trim(dat%varnames(1)), &
                             NF90_DOUBLE, dat%nodal_data_dims, dat%nodal_data_id)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'long_name', &
                             'dynamic water surface correction above water level')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'standard_name', &
                             'dynamic_sea_surface_correction_above_water_level')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'coordinates', 'time y x')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'location', 'node')
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             'mesh', 'adcirc_mesh')
         call check_err(iret)
         call putUnitsAttribute(dat%ncid, dat%nodal_data_id, &
                                'm', 'ft')
         iret = nf90_put_att(dat%ncid, dat%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)

      case DEFAULT
         write (scratchMessage, &
                '("No netCDF for files with unit number ",i0,".")') descript1%lun
         call allMessage(ERROR, scratchMessage)
      end select

      if (IDEN /= 44) then
         iret = nf90_put_att(dat%ncid, dat%u_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%v_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(dat%ncid, dat%w_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
      end if

!     jgf50.44: Automatically turn on compression if we are using the
!     netcdf4 file format.
#ifdef NETCDF_CAN_DEFLATE
      if (abs(descript1%specifier) == NETCDF4) then
         select case (descript1%lun)
         case (44)
            iret = nf90_def_var_deflate(dat%ncid, dat%u_nodal_data_id, &
                                        1, 1, 2)
            call check_err(iret)
            if ((IDEN == 2) .or. (IDEN == 4)) then
               iret = nf90_def_var_deflate(dat%ncid, dat%v_nodal_data_id, &
                                           1, 1, 2)
               call check_err(iret)
            end if
            if ((IDEN == 3) .or. (IDEN == 4)) then
               iret = nf90_def_var_deflate(dat%ncid, dat%w_nodal_data_id, &
                                           1, 1, 2)
               call check_err(iret)
            end if
         case (45, 46)
            iret = nf90_def_var_deflate(dat%ncid, dat%u_nodal_data_id, &
                                        1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(dat%ncid, dat%v_nodal_data_id, &
                                        1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(dat%ncid, dat%w_nodal_data_id, &
                                        1, 1, 2)
            call check_err(iret)
         case (47, 63, 73, 77, 90, 93, 173, 301:304, 306, 307, 402, 403, 108) ! GML added 93 20210727
            iret = nf90_def_var_deflate(dat%ncid, dat%nodal_data_id, &
                                        1, 1, 2)
            call check_err(iret)
         case (48)
            iret = nf90_def_var_deflate(dat%ncid, dat%u_nodal_data_id, &
                                        1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(dat%ncid, dat%v_nodal_data_id, &
                                        1, 1, 2)
            call check_err(iret)
!        !WJP 02.20.2018: added fort.53-54 harmonic outputs
         case (53)
            iret = nf90_def_var_deflate(dat%ncid, dat%ha_nodal_data_id, &
                                        1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(dat%ncid, dat%hg_nodal_data_id, &
                                        1, 1, 2)
            call check_err(iret)
         case (54)
            iret = nf90_def_var_deflate(dat%ncid, dat%u_ha_nodal_data_id, &
                                        1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(dat%ncid, dat%u_hg_nodal_data_id, &
                                        1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(dat%ncid, dat%v_ha_nodal_data_id, &
                                        1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(dat%ncid, dat%v_hg_nodal_data_id, &
                                        1, 1, 2)
            call check_err(iret)
         case (64, 74, 164, 305)
            iret = nf90_def_var_deflate(dat%ncid, dat%u_nodal_data_id, &
                                        1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(dat%ncid, dat%v_nodal_data_id, &
                                        1, 1, 2)
            call check_err(iret)
         case (311:315, 400, 401, 404) !adcirc max/min values
            iret = nf90_def_var_deflate(dat%ncid, dat%max_nodal_data_id, &
                                        1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(dat%ncid, dat%time_max_nodal_data_id, &
                                        1, 1, 2)
            call check_err(iret)
         case (316:322) !swan max/min values
            iret = nf90_def_var_deflate(dat%ncid, dat%max_nodal_data_id, &
                                        1, 1, 2)
            call check_err(iret)
         case DEFAULT ! should be unreachable
            write (scratchMessage, &
                   '("No netCDF for files with unit number ",i0,".")') descript1%lun
            call allMessage(ERROR, scratchMessage)
         end select
      end if
#endif

!     RJW added 9/13/2010
!     to include time atributes in global data files
!     Define time attributes
      call defineTimeAttributes(dat%ncid, dat%myTime)

!     define metadata and selected fort.15 parameters in netcdf file
      call defineMetaData(dat%ncid)

!     Leave define mode
      iret = nf90_enddef(dat%ncid)
      call check_err(iret)

!     write mesh to netcdf file
      call putMeshVariables(dat%ncid, dat%myMesh)
      !
      ! jgf52.08.11: Initialize max/min values b/c they will be read in
      ! at the start of a hotstart run. The square of the initialized
      ! values will be compared with the sum of the squares of
      ! vector values in the live running code. Uninitialized values
      ! of -99999 become large and positive when squared, meaning that
      ! we will never be able set a new peak in the running code if
      ! we leave the max/min values uninitialized.
      if ((dat%myFile%createFile .eqv. .true.) .and. &
          (descript1%readMaxMin .eqv. .true.)) then
         iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id, defaultValue)
         call check_err(iret)
      end if

      if (any(nodalLunsHA == descript1%lun)) then
         call addHarmonicAnalysisParametersToNetcdfFile(dat%ncid)
      end if
!     now close the initialized netcdf file
      iret = nf90_close(dat%ncid)
      call check_err(iret)

      call logMessage(INFO, 'Initialized "'//trim(dat%myFile%filename)// &
                      '" file.')

!-----------------------------------------------------------------------
   end subroutine initNodalDataFile
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     SUBROUTINE defineHarmonicAnalysisParametersInNetcdfFile
!-----------------------------------------------------------------------
!     Defines the variables that carry frequency, nodal factor, and
!     equilibrium arg in the netcdf output file
!-----------------------------------------------------------------------
   subroutine defineHarmonicAnalysisParametersInNetcdfFile(descript1, &
                                                           nodaldat, &
                                                           stationdat)
      use GLOBAL, only: OutputDataDescript_t
      implicit none
      type(nodalData), intent(inout), optional :: nodaldat
      type(stationData), intent(inout), optional :: stationdat
      type(OutputDataDescript_t), intent(inout) :: descript1
      integer :: iret, dmy

      if (present(nodalDat)) then
         nodalDat%myMesh%num_v_nodes = descript1%num_items_per_record
         call check_err(nf90_def_dim(nodalDat%ncid, 'constlen', 10, &
                                     nodalDat%myMesh%constl_dim_id))
         call check_err(nf90_def_dim(nodalDat%ncid, 'num_const', &
                                     nodalDat%myMesh%num_v_nodes, nodalDat%myMesh%num_v_nodes_dim_id))
!        Put the constituent names in
         nodalDat%nodal_data_dims(1) = nodalDat%myMesh%constl_dim_id
         nodalDat%nodal_data_dims(2) = nodalDat%myMesh%num_v_nodes_dim_id
         iret = nf90_def_var(nodalDat%ncid, 'const', NF90_CHAR, &
                             nodalDat%nodal_data_dims, nodalDat%myMesh%const_id)
         call check_err(iret)
         call check_err(nf90_put_att(nodalDat%ncid, nodalDat%myMesh%const_id, &
                                     'long_name', &
                                     'names of the tidal harmonic constituents'))
         call check_err(nf90_def_var(nodalDat%ncid, 'frequency', &
                                     NF90_DOUBLE, nodalDat%myMesh%num_v_nodes_dim_id, dmy))
         call check_err(nf90_put_att(nodalDat%ncid, dmy, &
                                     'long_name', 'frequency of harmonic constituents'))
         call check_err(nf90_put_att(nodalDat%ncid, dmy, 'units', 'rad/s'))
         call check_err(nf90_def_var(nodalDat%ncid, 'nodal_factor', &
                                     NF90_DOUBLE, nodalDat%myMesh%num_v_nodes_dim_id, dmy))
         call check_err(nf90_put_att(nodalDat%ncid, dmy, 'long_name', &
                                     'nodal factor of harmonic constituents'))
         call check_err(nf90_def_var(nodalDat%ncid, 'equilibrium_argument', &
                                     NF90_DOUBLE, nodalDat%myMesh%num_v_nodes_dim_id, dmy))
         call check_err(nf90_put_att(nodalDat%ncid, dmy, 'long_name', &
                                     'equilibrium argument of harmonic constituents'))
         call check_err(nf90_put_att(nodalDat%ncid, dmy, 'units', 'deg'))
      elseif (present(stationDat)) then
         call check_err(nf90_def_dim(stationDat%ncid, 'constlen', 10, &
                                     stationDat%constl_dim_id))
         call check_err(nf90_def_dim(stationDat%ncid, 'num_const', &
                                     stationDat%num_v_nodes, stationDat%num_v_nodes_dim_id))
!        Put the constituent names in
         stationDat%station_dims(1) = stationDat%constl_dim_id
         stationDat%station_dims(2) = stationDat%num_v_nodes_dim_id
         iret = nf90_def_var(stationDat%ncid, 'const', NF90_CHAR, &
                             stationDat%station_dims, stationDat%const_id)
         call check_err(iret)
         call check_err(nf90_put_att(stationDat%ncid, stationDat%const_id, &
                                     'long_name', &
                                     'names of the tidal harmonic constituents'))
         call check_err(nf90_def_var(stationDat%ncid, 'frequency', &
                                     NF90_DOUBLE, stationDat%num_v_nodes_dim_id, dmy))
         call check_err(nf90_put_att(stationDat%ncid, dmy, &
                                     'long_name', 'frequency of harmonic constituents'))
         call check_err(nf90_put_att(stationDat%ncid, dmy, 'units', 'rad/s'))
         call check_err(nf90_def_var(stationDat%ncid, 'nodal_factor', &
                                     NF90_DOUBLE, stationDat%num_v_nodes_dim_id, dmy))
         call check_err(nf90_put_att(stationDat%ncid, dmy, 'long_name', &
                                     'nodal factor of harmonic constituents'))
         call check_err(nf90_def_var(stationDat%ncid, 'equilibrium_argument', &
                                     NF90_DOUBLE, stationDat%num_v_nodes_dim_id, dmy))
         call check_err(nf90_put_att(stationDat%ncid, dmy, 'long_name', &
                                     'equilibrium argument of harmonic constituents'))
         call check_err(nf90_put_att(stationDat%ncid, dmy, 'units', 'deg'))
      end if

!-----------------------------------------------------------------------
   end subroutine defineHarmonicAnalysisParametersInNetcdfFile
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     SUBROUTINE addHarmonicAnalysisParametersToNetcdfFile
!-----------------------------------------------------------------------
!     Adds the variable data for frequency, nodal factor, and
!     equilibrium arg in the netcdf output file
!-----------------------------------------------------------------------
   subroutine addHarmonicAnalysisParametersToNetcdfFile(ncid)
      use GLOBAL, only: OutputDataDescript_t
      use HARM, only: NAMEFR, HAFREQ, HAFF, HAFACE
      implicit none
      integer, intent(IN) :: ncid
      integer :: dmy

      call check_err(nf90_inq_varid(ncid, "const", dmy))
      call check_err(nf90_put_var(ncid, dmy, NAMEFR(:)))
      call check_err(nf90_inq_varid(ncid, "frequency", dmy))
      call check_err(nf90_put_var(ncid, dmy, HAFREQ))
      call check_err(nf90_inq_varid(ncid, "nodal_factor", dmy))
      call check_err(nf90_put_var(ncid, dmy, HAFF))
      call check_err(nf90_inq_varid(ncid, "equilibrium_argument", dmy))
      call check_err(nf90_put_var(ncid, dmy, HAFACE))

!-----------------------------------------------------------------------
   end subroutine addHarmonicAnalysisParametersToNetcdfFile
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E  C R E A T E  N E T C D F  O U T P U T  F I L E
!-----------------------------------------------------------------------
!     jgf49.22 Creates a new netcdf output file if needed.
!-----------------------------------------------------------------------
   subroutine createNetCDFOutputFile(ncid, myFile, myTime, &
                                     descript, ret_error, static_time_length)
      use SIZES, only: globaldir
      use GLOBAL, only: OutputDataDescript_t, IHOT
      use mod_logging, only: allMessage
#ifdef CMPI
      use mod_logging, only: debug
#endif
      implicit none
      integer, intent(out) :: ncid
      type(fileData), intent(inout) :: myFile
      type(timeData), intent(inout) :: myTime
      type(OutputDataDescript_t), intent(inout) :: descript
      logical, intent(out) :: ret_error
      integer, intent(in), optional :: static_time_length
      integer :: iret
      character(1024) :: scratchMessage

      LOG_SCOPE_TRACED("createNetCDFOutputFile", NETCDFIO_TRACING)
      ret_error = .false.
      myFile%createFile = .false.
      myFile%fileFound = .false.

      ! create file name
      select case (descript%lun)
      case (67)
         myFile%filename = trim(globaldir)//'/fort.67.nc'
      case (68)
         myFile%filename = trim(globaldir)//'/fort.68.nc'
      case default
!        ! trim(globaldir) has already been prepended to
!        ! descript%file_name
         myFile%filename = trim(descript%file_name)//'.nc'
      end select
!     !jgfdebug
!     !write(6,*) 'myFile%filename=',trim(myFile%filename),
!     !&  ' descript%file_name=',trim(descript%file_name)
!
!     jgf49.17.02: Simplified the criteria for creating a new netCDF
!     output file: coldstart, overwrite upon hotstart, or output
!     that does not already exist. These criteria do not apply to
!     netcdf hotstart files; we will always write those when called
!     upon to create them.
      inquire (FILE=myFile%FILENAME, EXIST=myFile%fileFound)
      if ((descript%lun /= 67) .and. (descript%lun /= 68)) then
         if ((IHOT == 0) .or. (descript%specifier < 0) .or. &
             (myFile%fileFound .eqv. .false.)) then
#ifdef CMPI
            ! jgf49.31 when this subroutine is called by ADCIRC running in
            ! parallel, it should never create a new file, since that
            !  is the job of adcprep ... the file cannot be created
            ! here as a last resort since none of the processors have
            ! access to the full domain mesh and control files, whose
            ! data must also be stored in the netcdf output file
            if (myFile%fileFound .eqv. .false.) then
               write (scratchMessage, '(a,a,a)') &
                  "The NetCDF output file '", trim(myFile%FILENAME), &
                  "' was not found. "// &
                  "It should have been created by adcprep."
               call allMessage(ERROR, scratchMessage)
               ret_error = .true.
            end if
#else
            ! these lines are executed by serial adcirc and adcprep
            !
            ! jgf52.08.25: It may be a non-NetCDF min/max will be written,
            ! but a netcdf min/max file of the same type may be read. In this
            ! case, don't create a  netcdf min/max file of this type.
            if ((abs(descript%specifier) == NETCDF3) .or. &
                (abs(descript%specifier) == NETCDF4)) then
               myFile%createFile = .true.
               myFile%record_counter = 0
            end if
#endif
         end if
      else
         ! these lines are executed to create netcdf hotstart files
#ifdef CMPI
         if (myFile%fileFound .eqv. .true.) then
            call screenMessage(DEBUG, &
                               "Hotstart file was created by adcprep.")
            myFile%createFile = .false.
         else
            call screenMessage(ERROR, "Hotstart file is missing.")
            call screenMessage(ERROR, &
                               "It should have been created by adcprep.")
            ret_error = .true.
         end if
#else
         myFile%createFile = .true.
         myFile%record_counter = 0
#endif
      end if

!     RETURN if we don't need to create a file.
      if (myFile%createFile .eqv. .false.) then
         return
      end if

!     set the file format (netcdf3, or netcdf4/hdf5) ... the data model will
!     be classic (netcdf3) in any case
      select case (abs(descript%specifier))
      case (3, 367, 368)
         myFile%ncformat = NF90_CLOBBER
#ifdef HAVE_NETCDF4
      case (5, 567, 568)
         myFile%ncformat = ior(NF90_CLASSIC_MODEL, NF90_HDF5)
#else
      case (5, 567, 568)
         write (scratchMessage, '(A)') &
            "This ADCIRC executable was compiled with the " &
            //trim(nf90_inq_libvers())//" netcdf library."
         call allMessage(INFO, scratchMessage)
         write (scratchMessage, '(A,I3,A)') "File format specifier '", &
            descript%specifier, "' requires NetCDF version 4."
         call allMessage(ERROR, scratchMessage)
         write (scratchMessage, '(A)') &
            "It also requires the setting of NETCDF=enable and " &
            //"NETCDF4=enable on the make command line."
         call allMessage(ERROR, scratchMessage)
         write (scratchMessage, '(A)') &
            "You must recompile ADCIRC to use NetCDF4 formatting."
         call allMessage(ERROR, scratchMessage)
         ret_error = .true.
#endif
      case default
         write (scratchMessage, '(A,I3,A)') "File format specifier '", &
            descript%specifier, "' is not valid."
         call allMessage(ERROR, scratchMessage)
         ret_error = .true.
      end select

      if (ret_error .eqv. .false.) then
         iret = nf90_create(myFile%FILENAME, myFile%ncformat, ncid)
         call check_err(iret)
         ! Define time
         if (present(static_time_length)) then
            iret = nf90_def_dim(ncid, 'time', static_time_length, &
                                myTime%timenc_dim_id)
         else
            iret = nf90_def_dim(ncid, 'time', nf90_unlimited, &
                                myTime%timenc_dim_id)
         end if
         call check_err(iret)
         myTime%timenc_dims(1) = myTime%timenc_dim_id
         iret = nf90_def_var(ncid, 'time', NF90_DOUBLE, &
                             myTime%timenc_dims, myTime%timenc_id)
         call check_err(iret)
      end if

!-----------------------------------------------------------------------
   end subroutine createNetCDFOutputFile
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E   I N I T   N E T C D F   C O O R D
!-----------------------------------------------------------------------
!     jgf49.17.02 Initializes NetCDF coordinates.
!-----------------------------------------------------------------------
   subroutine initNetCDFCoord(myMesh)
      use GLOBAL, only: NP_G, NE_G
      use ADC_CONSTANTS, only: RAD2DEG
      use MESH, only: X, Y, SLAM, SFEA, ICS, NM
      use BOUNDARIES, only: NBOU, NVEL, NOPE, NBVV, NVDLL, NBDV, NVELL, &
                            NETA, IBTYPEE, IBTYPE, IBCONNR, barlanhtr, barinhtr, barincfspr, &
                            barlancfspr, barincfsbr, pipehtr, pipecoefr, pipediamr
      implicit none
      type(meshStructure), intent(inout) :: myMesh
      integer :: i, j, k, kj ! array indices

      LOG_SCOPE_TRACED("initNetCDFCoord", NETCDFIO_TRACING)

      myMesh%nopenc = nope
      myMesh%nbounc = nbou
      myMesh%netanc = neta
      myMesh%nvelnc = nvel
      myMesh%max_nvdllnc = maxval(nvdll) ! dimension of seg with most nodes
      myMesh%max_nvellnc = maxval(nvell) !WJP*2 dimension of seg with most nodes

      allocate (myMesh%xnc(NP_G))
      allocate (myMesh%ync(NP_G))
      allocate (myMesh%nvdllnc(myMesh%nopenc))
      allocate (myMesh%ibtypeenc(myMesh%nopenc))
      allocate (myMesh%ibtypenc(myMesh%nbounc))
      allocate (myMesh%nvellnc(myMesh%nbounc))
!     WJP to make nbdvnc arrays smaller:
      allocate (myMesh%nbdvnc(myMesh%netanc)) !(myMesh%nopenc,myMesh%max_nvdllnc))
      allocate (myMesh%nbvvnc(myMesh%nvelnc)) !(myMesh%nbounc,myMesh%max_nvellnc))
      allocate (myMesh%element(myMesh%nface_len, NE_G))
      allocate (myMesh%ibconnnc(myMesh%nvelnc)) !(myMesh%nbounc,myMesh%max_nvellnc))
      allocate (myMesh%barht(myMesh%nvelnc))
      allocate (myMesh%barsp(myMesh%nvelnc))
      allocate (myMesh%barsb(myMesh%nvelnc))
      allocate (myMesh%pipecoef(myMesh%nvelnc))
      allocate (myMesh%pipediam(myMesh%nvelnc))
      allocate (myMesh%pipeht(myMesh%nvelnc))
      allocate (myMesh%nmnc(NE_G, myMesh%nface_len))

!     Store nodal coordinates
      if (ics == 1) then
         myMesh%xnc = X
         myMesh%ync = Y
      else
         myMesh%xnc = SLAM*RAD2DEG ! convert back to degrees
         myMesh%ync = SFEA*RAD2DEG
      end if

!     elevation specified boundary forcing segments
      myMesh%nvdllnc = 0
      myMesh%nbdvnc = 0
      kj = 0 ! WJP initializing linear index
      do k = 1, myMesh%nopenc
         myMesh%nvdllnc(k) = nvdll(k)
         myMesh%ibtypeenc(k) = ibtypee(k)
         do j = 1, myMesh%nvdllnc(k)
            ! WJP counting the linear index
            kj = kj + 1
            myMesh%nbdvnc(kj) = nbdv(k, j)
         end do
      end do

!     normal flow (discharge) specified boundary segments
      myMesh%nvellnc = 0
      myMesh%nbvvnc = 0
      myMesh%ibconnnc = 0
      myMesh%barht = -99999d0
      myMesh%barsp = -99999d0
      myMesh%barsb = -99999d0
      myMesh%pipeht = -99999d0
      myMesh%pipecoef = -99999d0
      myMesh%pipediam = -99999d0
      kj = 0 !WJP initializing linear index
      do k = 1, myMesh%nbounc
         myMesh%nvellnc(k) = nvell(k)
         myMesh%ibtypenc(k) = ibtype(k)

         ! Set the boundary type specific variables
         select case (ibtype(k))
         case (3, 13, 23)
            myMesh%hasWeirs = .true.
         case (4, 14, 24)
            myMesh%hasWeirs = .true.
            myMesh%hasInternalWeirs = .true.
         case (5, 25)
            myMesh%hasWeirs = .true.
            myMesh%hasInternalWeirs = .true.
            myMesh%hasPipes = .true.
         end select

         do J = 1, myMesh%nvellnc(k)
            kj = kj + 1
            myMesh%nbvvnc(kj) = nbvv(k, j)
            if (ibtype(k) == 3 .or. ibtype(k) == 13 .or. ibtype(k) == 23) then
               myMesh%barht(kj) = barlanhtr(k, j)
               myMesh%barsp(kj) = barlancfspr(k, j)
            elseif (ibtype(k) == 4 .or. ibtype(k) == 24) then
               myMesh%ibconnnc(kj) = ibconnr(k, j)
               myMesh%barht(kj) = barinhtr(k, j)
               myMesh%barsp(kj) = barincfspr(k, j)
               myMesh%barsb(kj) = barincfsbr(k, j)
            elseif (ibtype(k) == 5 .or. ibtype(k) == 25) then
               myMesh%ibconnnc(kj) = ibconnr(k, j)
               myMesh%barht(kj) = barinhtr(k, j)
               myMesh%barsp(kj) = barincfspr(k, j)
               myMesh%barsb(kj) = barincfsbr(k, j)
               myMesh%pipeht(kj) = pipehtr(k, j)
               myMesh%pipecoef(kj) = pipecoefr(k, j)
               myMesh%pipediam(kj) = pipediamr(k, j)
            end if
         end do
      end do

      myMesh%nmnc = NM
!     Switch order in array for NETCDF
      do i = 1, NE_G
         do j = 1, myMesh%nface_len
            myMesh%element(j, i) = myMesh%nmnc(i, j)
         end do
      end do

      myMesh%initialized = .true.
!-----------------------------------------------------------------------
   end subroutine initNetCDFCoord
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E   W R I T E  O U T  A R R A Y  N E T C D F
!-----------------------------------------------------------------------
!     jgf48.03 This subroutine was created from cf's code in timestep.F to
!     write output files in NetCDF format.
!-----------------------------------------------------------------------
   subroutine writeOutArrayNetCDF(lun, timesec, descript1, &
                                  descript2, descript3, descript4)
      use GLOBAL, only: OutputDataDescript_t
      use GLOBAL_3DVS, only: outputFort48

      implicit none

      integer, intent(in) :: lun ! logical unit number of file to write to
      real(8), intent(in) :: timesec ! seconds since cold start
      type(OutputDataDescript_t), intent(in) :: descript1 !describes output data
      type(OutputDataDescript_t), intent(in), optional :: descript2 !describes output data
      type(OutputDataDescript_t), intent(in), optional :: descript3 !describes output data
      type(OutputDataDescript_t), intent(in), optional :: descript4 !describes output data
      character(1024) :: scratchMessage

      LOG_SCOPE_TRACED("writeOutArrayNetCDF", NETCDFIO_TRACING)

      select case (lun)
      case (41)
         call writeStationData(densityStations3D, lun, descript1, &
                               timesec, descript2, descript3)
      case (42)
         call writeStationData(velocityStations3D, lun, descript1, &
                               timesec, descript2, descript3)
      case (43)
         call writeStationData(turbulenceStations3D, lun, descript1, &
                               timesec, descript2, descript3)
      case (44)
         call writeNodalData(density3D, lun, descript1, &
                             timesec, descript2, descript3)
      case (45)
         call writeNodalData(velocity3D, lun, descript1, &
                             timesec, descript2, descript3)
      case (46)
         call writeNodalData(turbulence3D, lun, descript1, &
                             timesec, descript2, descript3)
      case (47)
         call writeNodalData(futureSurfaceTemperature, lun, &
                             descript1, timesec)
      case (48)
         if (outputFort48) then
            call writeNodalData(internal3D, lun, descript1, &
                                timesec, descript2, descript3)
         end if
      case (51)
         call writeStationData(HAelevSta, lun, descript1, timesec, &
                               descript2)
      case (52)
         call writeStationData(HAvelSta, lun, descript1, timesec, &
                               descript2, descript3, descript4)
      case (53)
         call writeNodalData(HAelev, lun, descript1, timesec, &
                             descript2)
      case (54)
         call writeNodalData(HAVel, lun, descript1, timesec, &
                             descript2, descript3, descript4)
      case (61)
         call writeStationData(elevSta, lun, descript1, timesec)
      case (62)
         call writeStationData(velSta, lun, descript1, timesec)
      case (63)
         call writeNodalData(elev, lun, descript1, timesec)
      case (64)
         call writeNodalData(currentVel, lun, descript1, timesec)
      case (71)
         call writeStationData(prSta, lun, descript1, timesec)
      case (72)
         call writeStationData(wVelSta, lun, descript1, timesec)
      case (73)
         call writeNodalData(pr, lun, descript1, timesec)
      case (74)
         call writeNodalData(windVel, lun, descript1, timesec)
      case (91)
         call writeStationData(IceSta, lun, descript1, timesec)
      case (93)
         call writeNodalData(CICEAF, lun, descript1, timesec)
      case (173)
         call writeNodalData(wDrag, lun, descript1, timesec)
      case (77)
         call writeNodalData(weirElev, lun, descript1, timesec)
      case (90)
         call writeNodalData(tau0nc, lun, descript1, timesec)
      case (164)
         call writeNodalData(rads, lun, descript1, timesec)
      case (301)
         call writeNodalData(sw_hs, lun, descript1, timesec)
      case (302)
         call writeNodalData(sw_dir, lun, descript1, timesec)
      case (303)
         call writeNodalData(sw_tm01, lun, descript1, timesec)
      case (304)
         call writeNodalData(sw_tps, lun, descript1, timesec)
      case (305)
         call writeNodalData(sw_wind, lun, descript1, timesec)
      case (306)
         call writeNodalData(sw_tm02, lun, descript1, timesec)
      case (307)
         call writeNodalData(sw_tmm10, lun, descript1, timesec)
      case (311)
         call writeNodalData(EtaMax, lun, descript1, timesec)
      case (312)
         call writeNodalData(UMax, lun, descript1, timesec)
      case (313)
         call writeNodalData(PrMin, lun, descript1, timesec)
      case (314)
         call writeNodalData(WVMax, lun, descript1, timesec)
      case (315)
         call writeNodalData(RSMax, lun, descript1, timesec)
      case (316)
         call writeNodalData(sw_hs_max, lun, descript1, timesec)
      case (317)
         call writeNodalData(sw_dir_max, lun, descript1, timesec)
      case (318)
         call writeNodalData(sw_tm01_max, lun, descript1, timesec)
      case (319)
         call writeNodalData(sw_tps_max, lun, descript1, timesec)
      case (320)
         call writeNodalData(sw_wind_max, lun, descript1, timesec)
      case (321)
         call writeNodalData(sw_tm02_max, lun, descript1, timesec)
      case (322)
         call writeNodalData(sw_tmm10_max, lun, descript1, timesec)
      case (400)
         call writeNodalData(inTime, lun, descript1, timesec)
      case (401)
         call writeNodalData(maxInDep, lun, descript1, timesec)
      case (402)
         call writeNodalData(coldDry, lun, descript1, timesec)
      case (403)
         call writeNodalData(eRisInun, lun, descript1, timesec)
      case (404)
         call writeNodalData(evrDry, lun, descript1, timesec)
      case (108)
         call writeNodalData(dynamiccorrection, lun, descript1, timesec)
      case (109)
         call writeStationData(dynamiccorrectionSta, lun, descript1, timesec)
      case DEFAULT
         write (scratchMessage, &
                '("No netCDF for files with unit number ",i0,".")') lun
         call allMessage(ERROR, scratchMessage)
      end select

!-----------------------------------------------------------------------
   end subroutine writeOutArrayNetCDF
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E    W R I T E  S T A T I O N  D A T A
!-----------------------------------------------------------------------
!     jgf49.17.02 Writes data to station file.
!-----------------------------------------------------------------------
   subroutine writeStationData(sta, lun, descript1, timesec, &
                               descript2, descript3, descript4)
      use SIZES, only: MNPROC
      use GLOBAL, only: OutputDataDescript_t, IDEN
      implicit none

      type(stationData), intent(inout) :: sta
      integer, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: descript1
      type(OutputDataDescript_t), intent(in), optional :: descript2
      type(OutputDataDescript_t), intent(in), optional :: descript3
      type(OutputDataDescript_t), intent(in), optional :: descript4
      real(8), intent(in) :: timesec

      integer :: kount(2), start(2)
      integer :: kount3D(4), start3D(3)
      integer :: iret ! success or failure of netcdf call
      character(1024) :: scratchMessage

      LOG_SCOPE_TRACED("writeStationData", NETCDFIO_TRACING)

!     jgf50.13: if netcdf output was requested, but there are no stations,
!     don't write to the file (it doesn't exist); just return.
      if (descript1%num_fd_records == 0) then
         return
      end if

      iret = nf90_open(sta%myFile%FILENAME, NF90_WRITE, sta%ncid)
      call check_err(iret)
      sta%myTime%timenc(sta%myTime%timenc_len) = timesec

      call setRecordCounterAndStoreTime(sta%ncid, sta%myFile, &
                                        sta%myTime)

      ! get number of vertical nodes for 3D stations
      if ((lun == 41) .or. (lun == 42) .or. (lun == 43)) then
         iret = nf90_inq_dimid(sta%ncid, "num_v_nodes", sta%num_v_nodes_dim_id)
         call check_err(iret)
         iret = nf90_inquire_dimension(sta%ncid, sta%num_v_nodes_dim_id, &
                                       len=sta%num_v_nodes)
         call check_err(iret)

         ! Set up the 3D netcdf data extents
         kount3D(1) = sta%num_stations
         kount3D(2) = sta%num_v_nodes
         kount3D(3) = sta%myTime%timenc_len
         start3D(1) = 1
         start3D(2) = 1
         start3D(3) = sta%myFile%record_counter
      else
         ! Set up the 2D netcdf data extents
         kount(1) = sta%num_stations
         kount(2) = sta%myTime%timenc_len
         start(1) = 1
         start(2) = sta%myFile%record_counter
      end if

      select case (lun)

      case (41)
         iret = nf90_inq_varid(sta%ncid, "sigmat", sta%u_station_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id, &
                                descript1%array2D, start3D, kount3D)
         else
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id, &
                                descript1%array2D_g, start3D, kount3D)
         end if
         call check_err(iret)

         if ((iden == 2) .or. (iden == 4)) then
            iret = nf90_inq_varid(sta%ncid, "salinity", &
                                  sta%v_station_data_id)
            call check_err(iret)
            if (MNPROC == 1) then
               iret = nf90_put_var(sta%ncid, sta%v_station_data_id, &
                                   descript2%array2D, start3D, kount3D)
            else
               iret = nf90_put_var(sta%ncid, sta%v_station_data_id, &
                                   descript2%array2D_g, start3D, kount3D)
            end if

            call check_err(iret)
         end if
         if ((iden == 3) .or. (iden == 4)) then
            iret = nf90_inq_varid(sta%ncid, "temperature", &
                                  sta%w_station_data_id)
            call check_err(iret)
            if (MNPROC == 1) then
               iret = nf90_put_var(sta%ncid, sta%w_station_data_id, &
                                   descript3%array2D, start3D, kount3D)
            else
               iret = nf90_put_var(sta%ncid, sta%w_station_data_id, &
                                   descript3%array2D_g, start3D, kount3D)
            end if
            call check_err(iret)
         end if
      case (42)
         iret = nf90_inq_varid(sta%ncid, "u-vel3D", sta%u_station_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(sta%ncid, "v-vel3D", sta%v_station_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(sta%ncid, "w-vel3D", sta%w_station_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id, &
                                descript1%array2D, start3D, kount3D)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%v_station_data_id, &
                                descript2%array2D, start3D, kount3D)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%w_station_data_id, &
                                descript3%array2D, start3D, kount3D)
            call check_err(iret)
         else
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id, &
                                descript1%array2D_g, start3D, kount3D)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%v_station_data_id, &
                                descript2%array2D_g, start3D, kount3D)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%w_station_data_id, &
                                descript3%array2D_g, start3D, kount3D)
            call check_err(iret)
         end if

      case (43)
         iret = nf90_inq_varid(sta%ncid, "q20", sta%u_station_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(sta%ncid, "l", sta%v_station_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(sta%ncid, "ev", sta%w_station_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id, &
                                descript1%array2D, start3D, kount3D)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%v_station_data_id, &
                                descript2%array2D, start3D, kount3D)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%w_station_data_id, &
                                descript3%array2D, start3D, kount3D)
            call check_err(iret)
         else
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id, &
                                descript1%array2D_g, start3D, kount3D)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%v_station_data_id, &
                                descript2%array2D_g, start3D, kount3D)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%w_station_data_id, &
                                descript3%array2D_g, start3D, kount3D)
            call check_err(iret)
         end if
         ! WJP 02.20.2018 Adding fort.51-52.nc capabilities
      case (51)
         iret = nf90_inq_varid(sta%ncid, "amp", sta%station_ha_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(sta%ncid, "phs", sta%station_hg_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_put_var(sta%ncid, sta%station_ha_data_id, &
                                descript1%array2D)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%station_hg_data_id, &
                                descript2%array2D)
            call check_err(iret)
         else
            iret = nf90_put_var(sta%ncid, sta%station_ha_data_id, &
                                descript1%array2D_g)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%station_hg_data_id, &
                                descript2%array2D_g)
            call check_err(iret)
         end if

      case (52)
         iret = nf90_inq_varid(sta%ncid, "u-vel-amp", &
                               sta%ha_u_station_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(sta%ncid, "u-vel-phs", &
                               sta%hg_u_station_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(sta%ncid, "v-vel-amp", &
                               sta%ha_v_station_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(sta%ncid, "v-vel-phs", &
                               sta%hg_v_station_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_put_var(sta%ncid, sta%ha_u_station_data_id, &
                                descript1%array2D)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%hg_u_station_data_id, &
                                descript2%array2D)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%ha_v_station_data_id, &
                                descript3%array2D)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%hg_v_station_data_id, &
                                descript4%array2D)
            call check_err(iret)
         else
            iret = nf90_put_var(sta%ncid, sta%ha_u_station_data_id, &
                                descript1%array2D_g)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%hg_u_station_data_id, &
                                descript2%array2D_g)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%ha_v_station_data_id, &
                                descript3%array2D_g)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%hg_v_station_data_id, &
                                descript4%array2D_g)
            call check_err(iret)
         end if

      case (61)
         iret = nf90_inq_varid(sta%ncid, "zeta", sta%station_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_put_var(sta%ncid, sta%station_data_id, &
                                descript1%array, start, kount)
         else
            iret = nf90_put_var(sta%ncid, sta%station_data_id, &
                                descript1%array_g, start, kount)
         end if
         call check_err(iret)
      case (62)
         iret = nf90_inq_varid(sta%ncid, "u-vel", sta%u_station_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(sta%ncid, "v-vel", sta%v_station_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id, &
                                descript1%array, start, kount)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%v_station_data_id, &
                                descript1%array2, start, kount)
         else
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id, &
                                descript1%array_g, start, kount)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%v_station_data_id, &
                                descript1%array2_g, start, kount)
         end if
      case (71)
         iret = nf90_inq_varid(sta%ncid, "pressure", sta%station_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_put_var(sta%ncid, sta%station_data_id, &
                                descript1%array, start, kount)
         else
            iret = nf90_put_var(sta%ncid, sta%station_data_id, &
                                descript1%array_g, start, kount)
         end if
         call check_err(iret)
      case (72)
         iret = nf90_inq_varid(sta%ncid, "windx", sta%u_station_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(sta%ncid, "windy", sta%v_station_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id, &
                                descript1%array, start, kount)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%v_station_data_id, &
                                descript1%array2, start, kount)
         else
            iret = nf90_put_var(sta%ncid, sta%u_station_data_id, &
                                descript1%array_g, start, kount)
            call check_err(iret)
            iret = nf90_put_var(sta%ncid, sta%v_station_data_id, &
                                descript1%array2_g, start, kount)
         end if

      case (91)
         iret = nf90_inq_varid(sta%ncid, "iceaf", sta%station_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_put_var(sta%ncid, sta%station_data_id, &
                                descript1%array, start, kount)
         else
            iret = nf90_put_var(sta%ncid, sta%station_data_id, &
                                descript1%array_g, start, kount)
         end if
         call check_err(iret)

      case (109)
         iret = nf90_inq_varid(sta%ncid, "dynamicWaterlevelCorrection", sta%station_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_put_var(sta%ncid, sta%station_data_id, &
                                descript1%array, start, kount)
         else
            iret = nf90_put_var(sta%ncid, sta%station_data_id, &
                                descript1%array_g, start, kount)
         end if
         call check_err(iret)

      case DEFAULT
         write (scratchMessage, &
                '("No netCDF for station files with unit number ",i0,".")') lun
         call allMessage(ERROR, scratchMessage)
      end select

!     Close netCDF file
      call check_err(nf90_close(sta%ncid))

!-----------------------------------------------------------------------
   end subroutine writeStationData
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E    W R I T E  N O D A L   D A T A
!-----------------------------------------------------------------------
!     jgf49.17.02 Writes data to a full domain file.
!-----------------------------------------------------------------------
   subroutine writeNodalData(dat, lun, descript1, timesec, &
                             descript2, descript3, descript4)
      use SIZES, only: MNPROC, MYPROC
      use GLOBAL, only: OutputDataDescript_t, NODECODE, IDEN

      implicit none

      type(nodalData), intent(inout) :: dat
      integer, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: descript1
      type(OutputDataDescript_t), intent(in), optional :: descript2
      type(OutputDataDescript_t), intent(in), optional :: descript3
      type(OutputDataDescript_t), intent(in), optional :: descript4
      real(8), intent(in) :: timesec

      integer :: kount(2), start(2)
      integer :: kount3D(3), start3D(3)
      integer :: iret ! success or failure of the netcdf call
      integer, allocatable  :: tempIntArray(:)
      real(8), allocatable :: tempArray(:)
      character(1024) :: scratchMessage

      LOG_SCOPE_TRACED("writeNodalData", NETCDFIO_TRACING)

      write (scratchMessage, '("Rank ",i0," writing to ",a,".")') myproc, &
         trim(dat%myFile%filename)
      call allMessage(DEBUG, scratchMessage)

      iret = nf90_open(dat%myFile%FILENAME, NF90_WRITE, dat%ncid)
      call check_err(iret)
      dat%myTime%timenc(dat%myTime%timenc_len) = timesec

      call setRecordCounterAndStoreTime(dat%ncid, dat%myFile, &
                                        dat%myTime)

!     Set up the 2D netcdf data extents
      kount(1) = dat%myMesh%num_nodes
      kount(2) = dat%myTime%timenc_len
      start(1) = 1
      start(2) = dat%myFile%record_counter

!     Set up the 3D netcdf data extents
      kount3D(1) = dat%myMesh%num_nodes
      kount3D(2) = dat%myMesh%num_v_nodes
      kount3D(3) = dat%myTime%timenc_len
      start3D(1) = 1
      start3D(2) = 1
      start3D(3) = dat%myFile%record_counter

!     Grab the data ids
      select case (lun)

      case (44)
         iret = nf90_inq_varid(dat%ncid, "sigmat", dat%u_nodal_data_id)
         call check_err(iret)
         if ((iden == 2) .or. (iden == 4)) then
            iret = nf90_inq_varid(dat%ncid, "salinity", &
                                  dat%v_nodal_data_id)
            call check_err(iret)
         end if
         if ((iden == 3) .or. (iden == 4)) then
            iret = nf90_inq_varid(dat%ncid, "temperature", &
                                  dat%w_nodal_data_id)
            call check_err(iret)
         end if
      case (45)
         iret = nf90_inq_varid(dat%ncid, "u-vel3D", dat%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "v-vel3D", dat%v_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "w-vel3D", dat%w_nodal_data_id)
         call check_err(iret)
      case (46)
         iret = nf90_inq_varid(dat%ncid, "q20", dat%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "l", dat%v_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "ev", dat%w_nodal_data_id)
         call check_err(iret)
      case (47)
         iret = nf90_inq_varid(dat%ncid, "qsurfkp1", dat%nodal_data_id)
         call check_err(iret)
      case (48)
         iret = nf90_inq_varid(dat%ncid, "bpgx", dat%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "bpgy", dat%v_nodal_data_id)
         call check_err(iret)
      case (53)
         iret = nf90_inq_varid(dat%ncid, "amp", dat%ha_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "phs", dat%hg_nodal_data_id)
         call check_err(iret)
      case (54)
         iret = nf90_inq_varid(dat%ncid, "u_amp", dat%ha_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "u_phs", dat%hg_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "v_amp", dat%v_ha_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "v_phs", dat%v_hg_nodal_data_id)
         call check_err(iret)
      case (63)
         iret = nf90_inq_varid(dat%ncid, "zeta", dat%nodal_data_id)
         call check_err(iret)
      case (64)
         iret = nf90_inq_varid(dat%ncid, "u-vel", dat%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "v-vel", dat%v_nodal_data_id)
         call check_err(iret)
      case (73)
         iret = nf90_inq_varid(dat%ncid, "pressure", dat%nodal_data_id)
         call check_err(iret)
      case (74)
         iret = nf90_inq_varid(dat%ncid, "windx", dat%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "windy", dat%v_nodal_data_id)
         call check_err(iret)
      case (93)
         iret = nf90_inq_varid(dat%ncid, "iceaf", dat%nodal_data_id)
         call check_err(iret)
      case (173)
         iret = nf90_inq_varid(dat%ncid, "winddrag", dat%nodal_data_id)
         call check_err(iret)
      case (77)
         iret = nf90_inq_varid(dat%ncid, "weir_dz", dat%nodal_data_id)
         call check_err(iret)
      case (90)
         iret = nf90_inq_varid(dat%ncid, "tau0", dat%nodal_data_id)
         call check_err(iret)
      case (164)
         iret = nf90_inq_varid(dat%ncid, "radstress_x", dat%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "radstress_y", dat%v_nodal_data_id)
         call check_err(iret)
#ifdef CSWAN
      case (301)
         iret = nf90_inq_varid(dat%ncid, "swan_HS", dat%nodal_data_id)
         call check_err(iret)
      case (302)
         iret = nf90_inq_varid(dat%ncid, "swan_DIR", dat%nodal_data_id)
         call check_err(iret)
      case (303)
         iret = nf90_inq_varid(dat%ncid, "swan_TM01", dat%nodal_data_id)
         call check_err(iret)
      case (304)
         iret = nf90_inq_varid(dat%ncid, "swan_TPS", dat%nodal_data_id)
         call check_err(iret)
      case (305)
         iret = nf90_inq_varid(dat%ncid, "swan_windx", dat%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "swan_windy", dat%v_nodal_data_id)
         call check_err(iret)
      case (306)
         iret = nf90_inq_varid(dat%ncid, "swan_TM02", dat%nodal_data_id)
         call check_err(iret)
      case (307)
         iret = nf90_inq_varid(dat%ncid, "swan_TMM10", dat%nodal_data_id)
         call check_err(iret)
#endif

      case (311)
         iret = nf90_inq_varid(dat%ncid, "zeta_max", dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "time_of_zeta_max", &
                               dat%time_max_nodal_data_id)
         call check_err(iret)
      case (312)
         iret = nf90_inq_varid(dat%ncid, "vel_max", &
                               dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "time_of_vel_max", &
                               dat%time_max_nodal_data_id)
         call check_err(iret)
      case (313)
         iret = nf90_inq_varid(dat%ncid, "pressure_min", &
                               dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "time_of_pressure_min", &
                               dat%time_max_nodal_data_id)
         call check_err(iret)
      case (314)
         iret = nf90_inq_varid(dat%ncid, "wind_max", &
                               dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "time_of_wind_max", &
                               dat%time_max_nodal_data_id)
         call check_err(iret)
      case (315)
         iret = nf90_inq_varid(dat%ncid, "radstress_max", &
                               dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "time_of_radstress_max", &
                               dat%time_max_nodal_data_id)
         call check_err(iret)
#ifdef CSWAN
      case (316)
         iret = nf90_inq_varid(dat%ncid, "swan_HS_max", &
                               dat%max_nodal_data_id)
         call check_err(iret)
      case (317)
         iret = nf90_inq_varid(dat%ncid, "swan_DIR_max", &
                               dat%max_nodal_data_id)
         call check_err(iret)
      case (318)
         iret = nf90_inq_varid(dat%ncid, "swan_TM01_max", &
                               dat%max_nodal_data_id)
         call check_err(iret)
      case (319)
         iret = nf90_inq_varid(dat%ncid, "swan_TPS_max", &
                               dat%max_nodal_data_id)
         call check_err(iret)
      case (320)
         iret = nf90_inq_varid(dat%ncid, "swan_wind_max", &
                               dat%max_nodal_data_id)
         call check_err(iret)
      case (321)
         iret = nf90_inq_varid(dat%ncid, "swan_TM02_max", &
                               dat%max_nodal_data_id)
         call check_err(iret)
      case (322)
         iret = nf90_inq_varid(dat%ncid, "swan_TMM10_max", &
                               dat%max_nodal_data_id)
         call check_err(iret)
#endif

         ! inundationtime.63
      case (400)
         iret = nf90_inq_varid(dat%ncid, "inun_time", &
                               dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "onset_inun_time", &
                               dat%time_max_nodal_data_id)
         call check_err(iret)
         ! maxinundepth.63
      case (401)
         iret = nf90_inq_varid(dat%ncid, "inun_max", &
                               dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "time_of_inun_max", &
                               dat%time_max_nodal_data_id)
         call check_err(iret)
         ! initiallydry.63
      case (402)
         iret = nf90_inq_varid(dat%ncid, "initiallydry", dat%nodal_data_id)
         call check_err(iret)
         ! endrisinginun.63
      case (403)
         iret = nf90_inq_varid(dat%ncid, "endrisinginun", dat%nodal_data_id)
         call check_err(iret)
      case (404)
         iret = nf90_inq_varid(dat%ncid, "everdried", &
                               dat%max_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(dat%ncid, "time_of_everdried", &
                               dat%time_max_nodal_data_id)
         call check_err(iret)

      case (108)
         iret = nf90_inq_varid(dat%ncid, "dynamicWaterlevelCorrection", dat%nodal_data_id)
         call check_err(iret)

      case DEFAULT
         write (scratchMessage, &
                '("No netCDF for files with unit number ",i0,".")') lun
         call allMessage(ERROR, scratchMessage)
      end select

!     Write the array values
      select case (lun)
      case (44)
         if (MNPROC == 1) then
            iret = nf90_put_var(dat%ncid, dat%u_nodal_data_id, &
                                descript1%array2D, start3D, kount3D)
         else
            iret = nf90_put_var(dat%ncid, dat%u_nodal_data_id, &
                                descript1%array2D_g, start3D, kount3D)
         end if
         call check_err(iret)

         if ((iden == 2) .or. (iden == 4)) then
            if (MNPROC == 1) then
               iret = nf90_put_var(dat%ncid, dat%v_nodal_data_id, &
                                   descript2%array2D, start3D, kount3D)
            else
               iret = nf90_put_var(dat%ncid, dat%v_nodal_data_id, &
                                   descript2%array2D_g, start3D, kount3D)
            end if
            call check_err(iret)
         end if
         if ((iden == 3) .or. (iden == 4)) then
            if (MNPROC == 1) then
               iret = nf90_put_var(dat%ncid, dat%w_nodal_data_id, &
                                   descript3%array2D, start3D, kount3D)
            else
               iret = nf90_put_var(dat%ncid, dat%w_nodal_data_id, &
                                   descript3%array2D_g, start3D, kount3D)
            end if
            call check_err(iret)
         end if
      case (45, 46)
         if (MNPROC == 1) then
            iret = nf90_put_var(dat%ncid, dat%u_nodal_data_id, &
                                descript1%array2D, start3D, kount3D)
            call check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%v_nodal_data_id, &
                                descript2%array2D, start3D, kount3D)
            call check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%w_nodal_data_id, &
                                descript3%array2D, start3D, kount3D)
            call check_err(iret)
         else
            iret = nf90_put_var(dat%ncid, dat%u_nodal_data_id, &
                                descript1%array2D_g, start3D, kount3D)
            call check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%v_nodal_data_id, &
                                descript2%array2D_g, start3D, kount3D)
            call check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%w_nodal_data_id, &
                                descript3%array2D_g, start3D, kount3D)
            call check_err(iret)
         end if
      case (48)
         if (MNPROC == 1) then
            iret = nf90_put_var(dat%ncid, dat%u_nodal_data_id, &
                                descript1%array2D, start3D, kount3D)
            call check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%v_nodal_data_id, &
                                descript2%array2D, start3D, kount3D)
            call check_err(iret)
         else
            iret = nf90_put_var(dat%ncid, dat%u_nodal_data_id, &
                                descript1%array2D_g, start3D, kount3D)
            call check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%v_nodal_data_id, &
                                descript2%array2D_g, start3D, kount3D)
            call check_err(iret)
         end if
      case (53, 54)
         if (MNPROC == 1) then
            iret = nf90_put_var(dat%ncid, dat%ha_nodal_data_id, &
                                descript1%array2D)
            call check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%hg_nodal_data_id, &
                                descript2%array2D)
            call check_err(iret)
            if (lun == 54) then
               iret = nf90_put_var(dat%ncid, dat%v_ha_nodal_data_id, &
                                   descript3%array2D)
               call check_err(iret)
               iret = nf90_put_var(dat%ncid, dat%v_hg_nodal_data_id, &
                                   descript4%array2D)
               call check_err(iret)
            end if
         else
            iret = nf90_put_var(dat%ncid, dat%ha_nodal_data_id, &
                                descript1%array2D_g)
            call check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%hg_nodal_data_id, &
                                descript2%array2D_g)
            call check_err(iret)
            if (lun == 54) then
               iret = nf90_put_var(dat%ncid, dat%v_ha_nodal_data_id, &
                                   descript3%array2D_g)
               call check_err(iret)
               iret = nf90_put_var(dat%ncid, dat%v_hg_nodal_data_id, &
                                   descript4%array2D_g)
               call check_err(iret)
            end if
         end if
      case (316, 317, 318, 319, 321, 322) !swan max/min
         if (MNPROC == 1) then ! SERIAL
            if (descript1%ConsiderWetDry .eqv. .true.) then
               allocate (tempArray(1:dat%myMesh%num_nodes))
               tempArray = merge(descript1%array, &
                                 spread(descript1%alternate_value, 1, kount(1)), nodecode > 0)
               iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id, tempArray)
            else
               iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id, &
                                   descript1%array)
            end if
            call check_err(iret)
         else
            iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id, &
                                descript1%array_g)
            call check_err(iret)
         end if
      case (311, 313) !adcirc max/min only
         if (MNPROC == 1) then ! SERIAL
            if (descript1%ConsiderWetDry .eqv. .true.) then
               allocate (tempArray(1:dat%myMesh%num_nodes))
               tempArray = merge(descript1%array, &
                                 spread(descript1%alternate_value, 1, kount(1)), nodecode > 0)
               iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id, tempArray)
               call check_err(iret)
               tempArray = merge(descript1%array2, &
                                 spread(descript1%alternate_value, 1, kount(1)), nodecode > 0)
               iret = nf90_put_var(dat%ncid, dat%time_max_nodal_data_id, tempArray)
               call check_err(iret)
            else
               iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id, &
                                   descript1%array)
               call check_err(iret)
               iret = nf90_put_var(dat%ncid, dat%time_max_nodal_data_id, &
                                   descript1%array2)
               call check_err(iret)
            end if
         else
            iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id, &
                                descript1%array_g)
            call check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%time_max_nodal_data_id, &
                                descript1%array2_g)
            call check_err(iret)
         end if
      case (312, 314, 315, 400, 401, 404) !adcirc max/min
         if (MNPROC == 1) then ! SERIAL
            if (descript1%ConsiderWetDry .eqv. .true.) then
               allocate (tempArray(1:dat%myMesh%num_nodes))
               tempArray = merge(descript1%array, &
                                 spread(descript1%alternate_value, 1, kount(1)), nodecode > 0)
               iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id, tempArray)
               call check_err(iret)
               tempArray = merge(descript1%array2, &
                                 spread(descript1%alternate_value, 1, kount(1)), nodecode > 0)
               iret = nf90_put_var(dat%ncid, dat%time_max_nodal_data_id, tempArray)
               call check_err(iret)
            else
               iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id, &
                                   descript1%array)
               call check_err(iret)
               iret = nf90_put_var(dat%ncid, dat%time_max_nodal_data_id, &
                                   descript1%array2)
               call check_err(iret)
            end if
         else ! PARALLEL
            iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id, &
                                descript1%array_g)
            call check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%time_max_nodal_data_id, &
                                descript1%array2_g)
            call check_err(iret)
         end if
      case (320) !swan max/min
         if (MNPROC == 1) then ! SERIAL
            if (descript1%ConsiderWetDry .eqv. .true.) then
               allocate (tempArray(1:dat%myMesh%num_nodes))
               tempArray = merge(descript1%array, &
                                 spread(descript1%alternate_value, 1, kount(1)), nodecode > 0)
               iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id, tempArray)
               call check_err(iret)
            else
               iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id, &
                                   descript1%array)
               call check_err(iret)
            end if
         else ! PARALLEL
            iret = nf90_put_var(dat%ncid, dat%max_nodal_data_id, &
                                descript1%array_g)
            call check_err(iret)
         end if
      case (47, 63, 73, 77, 90, 93, 173, 301, 302, 303, 304, 306, 307, 108) ! GML added 93 20210727
         if (MNPROC == 1) then ! SERIAL
            if (descript1%ConsiderWetDry .eqv. .true.) then
               allocate (tempArray(1:dat%myMesh%num_nodes))
               tempArray = merge(descript1%array, &
                                 spread(descript1%alternate_value, 1, kount(1)), nodecode > 0)

               iret = nf90_put_var(dat%ncid, dat%nodal_data_id, tempArray, start, kount)
            else
               iret = nf90_put_var(dat%ncid, dat%nodal_data_id, &
                                   descript1%array, start, kount)
            end if
            call check_err(iret)
         else
            iret = nf90_put_var(dat%ncid, dat%nodal_data_id, &
                                descript1%array_g, start, kount)
            call check_err(iret)
         end if
      case (64, 74, 164, 305)
         if (MNPROC == 1) then ! SERIAL
            if (descript1%ConsiderWetDry .eqv. .true.) then
               allocate (tempArray(1:dat%myMesh%num_nodes))
               tempArray = merge(descript1%array, &
                                 spread(descript1%alternate_value, 1, kount(1)), nodecode > 0)
               iret = nf90_put_var(dat%ncid, dat%u_nodal_data_id, tempArray, start, kount)
               call check_err(iret)
               tempArray = merge(descript1%array2, &
                                 spread(descript1%alternate_value, 1, kount(1)), nodecode > 0)
               iret = nf90_put_var(dat%ncid, dat%v_nodal_data_id, tempArray, start, kount)
               call check_err(iret)
            else
               iret = nf90_put_var(dat%ncid, dat%u_nodal_data_id, &
                                   descript1%array, start, kount)
               call check_err(iret)
               iret = nf90_put_var(dat%ncid, dat%v_nodal_data_id, &
                                   descript1%array2, start, kount)
               call check_err(iret)
            end if
         else ! PARALLEL
            iret = nf90_put_var(dat%ncid, dat%u_nodal_data_id, &
                                descript1%array_g, start, kount)
            call check_err(iret)
            iret = nf90_put_var(dat%ncid, dat%v_nodal_data_id, &
                                descript1%array2_g, start, kount)
            call check_err(iret)
         end if
      case (402, 403) ! integers at nodes (e.g., initiallydry.63, endrisinginun.63)
         if (MNPROC == 1) then ! SERIAL
            if (descript1%ConsiderWetDry .eqv. .true.) then
               allocate (tempIntArray(1:dat%myMesh%num_nodes))
               tempIntArray = merge(descript1%iarray, &
                                    spread(int(descript1%alternate_value), 1, kount(1)), nodecode > 0)
               iret = nf90_put_var(dat%ncid, dat%nodal_data_id, tempIntArray, start, kount)
            else
               iret = nf90_put_var(dat%ncid, dat%nodal_data_id, &
                                   descript1%iarray, start, kount)
            end if
            call check_err(iret)
         else
            iret = nf90_put_var(dat%ncid, dat%nodal_data_id, &
                                descript1%iarray_g, start, kount)
            call check_err(iret)
         end if
      case DEFAULT
         write (scratchMessage, &
                '("No netCDF for files with unit number ",i0,".")') lun
         call allMessage(ERROR, scratchMessage)
      end select

!     Close netCDF file
      call check_err(nf90_close(dat%ncid))

!-----------------------------------------------------------------------
   end subroutine writeNodalData
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!      S E T   R E C O R D   C O U N T E R   A N D   S T O R E   T I M E
!-----------------------------------------------------------------------
!     jgf49.48.08 Compares the current ADCIRC simulation time with the
!     array of output times in the file, and if the simulation time
!     is before the end of the file, it sets the record counter to the
!     right place within the existing data. Data that occur after the
!     inserted data will remain, due to the inability of netcdf to
!     delete data from files.
!-----------------------------------------------------------------------
   subroutine setRecordCounterAndStoreTime(ncid, f, t)
      implicit none

      integer, intent(in) :: ncid
      type(fileData), intent(inout) :: f
      type(timeData), intent(inout) :: t

      real(8), allocatable :: storedTimes(:) ! array of time values in file
      logical :: timeFound ! true if current time is in array of stored times
      integer :: ndim ! number of dimensions in the netcdf file
      integer :: nvar ! number of variables in the netcdf file
      integer :: natt ! number of attributes in the netcdf file
      integer :: counti(1), starti(1)
      integer :: iret ! success or failure of netcdf call
      integer :: i ! loop counter
      character(1024) :: scratchMessage, scratchFormat

      LOG_SCOPE_TRACED("setRecordCounterAndStoreTime", NETCDFIO_TRACING)

!     Inquire time variable
      iret = nf90_inquire(ncid, ndim, nvar, natt, t%timenc_dim_id)
      call check_err(iret)
      iret = nf90_inquire_dimension(ncid, t%timenc_dim_id, &
                                    len=f%record_counter)
      call check_err(iret)
      iret = nf90_inq_varid(ncid, "time", t%timenc_id)
      call check_err(iret)

!     Determine the relationship between the current simulation time
!     and the time array stored in the netcdf file. Set the record
!     counter based on this relationship.
      if (f%record_counter /= 0) then
         allocate (storedTimes(f%record_counter))
         iret = nf90_get_var(ncid, t%timenc_id, storedTimes)
         call check_err(iret)
         timeFound = .false.
         do i = 1, f%record_counter
            if ((t%timenc(1) < storedTimes(i)) .or. &
                (abs(t%timenc(1) - storedTimes(i)) < 1.0d-10)) then
               timeFound = .true.
               exit
            end if
         end do
         if (timeFound .eqv. .false.) then
            ! Increment the record counter so that we can store data at the
            ! next location in the netcdf file (i.e., all of the times
            ! in the netcdf file were found to be earlier than the current
            ! adcirc simulation time).
            f%record_counter = f%record_counter + 1
         else
            ! jgf49.48.08: set the counter at the index that reflects the
            ! current time within the netcdf file (or is between two times
            ! found in the netcdf file).
            ! WARNING: all subsequent data will remain in the file, we
            ! are just overwriting it ... if we don't overwrite all of it,
            ! the pre-existing data will still be there, which is probably
            ! not what the user intended ... but apparently there is no
            ! way to delete data from netcdf files:
            ! http://www.unidata.ucar.edu/support/help/MailArchives/netcdf/msg02367.html
            scratchFormat = &
               '("Overwriting pre-existing data in netcdf file ",a,'// &
               '" for time=",f17.8,". '// &
               'Subsequent data in netcdf file remain unchanged.")'
            write (scratchMessage, scratchFormat) &
               trim(f%FILENAME), t%timenc(1)
            call allMessage(INFO, scratchMessage)
            f%record_counter = i
         end if
         deallocate (storedTimes)
      else
         ! set the counter at 1 so we can record our first time value
         f%record_counter = 1
      end if

!     Store simulation time in netcdf file
      starti(1) = f%record_counter
      counti(1) = t%timenc_len
      iret = nf90_put_var(ncid, t%timenc_id, t%timenc, starti, counti)
      call check_err(iret)

!-----------------------------------------------------------------------
   end subroutine setRecordCounterAndStoreTime
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E    I N I T   N E T C D F   H O T S T A R T
!-----------------------------------------------------------------------
!     jgf49.35 Sets up netCDF variables and writes mesh data into netcdf
!     hotstart file.
!-----------------------------------------------------------------------
   subroutine initNetCDFHotstart(lun, Elev1Descript, ncerror)
      use GLOBAL, only: OutputDataDescript_t, IM, IMHS, NE_G, NP_G
      use MESH, only: ICS
      implicit none

      integer, intent(in) :: lun
      type(OutputDataDescript_t), intent(inout) :: Elev1Descript
      logical, intent(out) :: ncerror

      integer :: iret ! success or failure of the netcdf call
      integer :: tempid

      LOG_SCOPE_TRACED("initNetCDFHotstart", NETCDFIO_TRACING)
      ncerror = .false.

!     Point to the hotstart file we want to work on.
      if (lun == 67) then
         hs => hs67
      else
         hs => hs68
      end if
      hs%myMesh => adcircMesh
      if (hs%myMesh%initialized .eqv. .false.) then
         hs%myMesh%num_nodes = NP_G
         hs%myMesh%num_elems = NE_G
         hs%myMesh%nface_len = 3
      end if
      if (hs%myTime%initialized .eqv. .false.) then
         allocate (hs%myTime%timenc(hs%myTime%timenc_len))
         hs%myTime%initialized = .true.
      end if

!     Initialize netCDF hotstart file, creating a new one
      Elev1Descript%lun = lun
      call createNetCDFOutputFile(hs%ncid, hs%myFile, hs%myTime, &
                                  Elev1Descript, ncerror, 1)
      ! return an error flag to the calling routine if something went
      ! wrong when we tried to create the netcdf file
      if (ncerror .eqv. .true.) then
         return
      end if
      if (hs%myMesh%initialized .eqv. .false.) then
         hs%myMesh%num_nodes = NP_G
         hs%myMesh%num_elems = NE_G
         hs%myMesh%nface_len = 3
         call initNetCDFCoord(hs%myMesh)
      end if
      call defineMeshVariables(hs%ncid, hs%myMesh, hs%myFile)

!     Z E T A 1
      hs%zeta1%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%zeta1%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      call define_netcdf_hotstart_variable(hs, 'zeta1', hs%zeta1%nodal_data_dims, NF90_DOUBLE, &
                                           hs%zeta1%nodal_data_id)

      iret = nf90_put_att(hs%ncid, hs%zeta1%nodal_data_id, &
                          'long_name', 'water surface elevation at previous time step')
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%zeta1%nodal_data_id, &
                          'standard_name', &
                          'water_surface_elevation_at_previous_time step')
      call check_err(iret)
      call putUnitsAttribute(hs%ncid, hs%zeta1%nodal_data_id, &
                             'm', 'ft')
      iret = nf90_put_att(hs%ncid, hs%zeta1%nodal_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%zeta1%nodal_data_id, &
                          'positive', 'up')
      call check_err(iret)

!     Z E T A 2
      hs%zeta2%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%zeta2%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      call define_netcdf_hotstart_variable(hs, 'zeta2', hs%zeta2%nodal_data_dims, NF90_DOUBLE, &
                                           hs%zeta2%nodal_data_id)

      iret = nf90_put_att(hs%ncid, hs%zeta2%nodal_data_id, &
                          'long_name', 'water surface elevation at current time step')
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%zeta2%nodal_data_id, &
                          'standard_name', &
                          'water_surface_elevation_at_current_time_step')
      call check_err(iret)
      call putUnitsAttribute(hs%ncid, hs%zeta2%nodal_data_id, &
                             'm', 'ft')
      iret = nf90_put_att(hs%ncid, hs%zeta2%nodal_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%zeta2%nodal_data_id, &
                          'positive', 'up')
      call check_err(iret)

!     Z E T A D
      hs%zetad%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%zetad%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      call define_netcdf_hotstart_variable(hs, 'zetad', hs%zetad%nodal_data_dims, NF90_DOUBLE, &
                                           hs%zetad%nodal_data_id)

      iret = nf90_put_att(hs%ncid, hs%zetad%nodal_data_id, &
                          'long_name', &
                          'water elevation at flux specified boundary')
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%zetad%nodal_data_id, &
                          'standard_name', &
                          'water_elevation_at_flux_specified_boundary')
      call check_err(iret)
      call putUnitsAttribute(hs%ncid, hs%zetad%nodal_data_id, &
                             'm', 'ft')
      iret = nf90_put_att(hs%ncid, hs%zetad%nodal_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%zetad%nodal_data_id, 'positive', &
                          'up')
      call check_err(iret)

!     H 1 : total water depth at the previous time step
      hs%htot1%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%htot1%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      call define_netcdf_hotstart_variable(hs, 'h1', hs%htot1%nodal_data_dims, NF90_DOUBLE, &
                                           hs%htot1%nodal_data_id)

      iret = nf90_put_att(hs%ncid, hs%htot1%nodal_data_id, &
                          'long_name', 'total water depth at previous time step')
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%htot1%nodal_data_id, &
                          'standard_name', &
                          'total_water_depth_at_previous_time step')
      call check_err(iret)
      call putUnitsAttribute(hs%ncid, hs%htot1%nodal_data_id, &
                             'm', 'ft')
      iret = nf90_put_att(hs%ncid, hs%htot1%nodal_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)

!     H 2 : total water depth at the current time step
      hs%htot2%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%htot2%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      call define_netcdf_hotstart_variable(hs, 'h2', hs%htot2%nodal_data_dims, NF90_DOUBLE, &
                                           hs%htot2%nodal_data_id)

      iret = nf90_put_att(hs%ncid, hs%htot2%nodal_data_id, &
                          'long_name', 'total water depth at current time step')
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%htot2%nodal_data_id, &
                          'standard_name', &
                          'total_water_depth_at_current_time_step')
      call check_err(iret)
      call putUnitsAttribute(hs%ncid, hs%htot2%nodal_data_id, &
                             'm', 'ft')
      iret = nf90_put_att(hs%ncid, hs%htot2%nodal_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)

!     U V E L
      hs%vel%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%vel%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      call define_netcdf_hotstart_variable(hs, 'u-vel', hs%vel%nodal_data_dims, NF90_DOUBLE, &
                                           hs%vel%u_nodal_data_id)

      if (ics /= 1) then
         iret = nf90_put_att(hs%ncid, hs%vel%u_nodal_data_id, &
                             'long_name', 'vertically averaged water velocity u-component')
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%vel%u_nodal_data_id, &
                             'positive', 'east')
      else
         iret = nf90_put_att(hs%ncid, hs%vel%u_nodal_data_id, &
                             'long_name', 'vertically averaged velocity in x-direction')
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%vel%u_nodal_data_id, &
                             'positive', 'right')
      end if
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%vel%u_nodal_data_id, &
                          'standard_name', 'u_velocity')
      call check_err(iret)
      call putUnitsAttribute(hs%ncid, hs%vel%u_nodal_data_id, &
                             'm s-1', 'ft s-1')
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%vel%u_nodal_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%vel%u_nodal_data_id, &
                          'dry_Value', doubleval)
      call check_err(iret)

!     V V E L
      call define_netcdf_hotstart_variable(hs, 'v-vel', hs%vel%nodal_data_dims, NF90_DOUBLE, &
                                           hs%vel%v_nodal_data_id)

      if (ics /= 1) then
         iret = nf90_put_att(hs%ncid, hs%vel%v_nodal_data_id, &
                             'long_name', 'vertically averaged water velocity v-component')
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%vel%v_nodal_data_id, &
                             'positive', 'north')
      else
         iret = nf90_put_att(hs%ncid, hs%vel%v_nodal_data_id, &
                             'long_name', 'vertically averaged water velocity in y-direction')
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%vel%v_nodal_data_id, &
                             'positive', '90 degrees counterclockwise from x water velocity')
      end if
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%vel%v_nodal_data_id, &
                          'standard_name', 'v_velocity')
      call check_err(iret)
      call putUnitsAttribute(hs%ncid, hs%vel%v_nodal_data_id, &
                             'm s-1', 'ft s-1')
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%vel%v_nodal_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%vel%v_nodal_data_id, &
                          'dry_Value', doubleval)
      call check_err(iret)

!     C H 1
      if ((IM == 10) .or. (IMHS == 10)) then
         hs%ch1%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
         hs%ch1%nodal_data_dims(2) = hs%myTime%timenc_dim_id
         call define_netcdf_hotstart_variable(hs, 'ch1', hs%ch1%nodal_data_dims, NF90_DOUBLE, &
                                              hs%ch1%nodal_data_id)
         iret = nf90_put_att(hs%ncid, hs%ch1%nodal_data_id, &
                             'long_name', 'concentration')
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%ch1%nodal_data_id, &
                             'standard_name', 'concentration')
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%ch1%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
      end if

!     N O D E C O D E
      hs%nodecodenc%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%nodecodenc%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      call define_netcdf_hotstart_variable(hs, 'nodecode', hs%nodecodenc%nodal_data_dims, NF90_INT, &
                                           hs%nodecodenc%nodal_data_id)
      att_text = "wet or dry state of node where 1 indicates that the" &
                 //" node is wet and 0 indicates that the node is dry"
      iret = nf90_put_att(hs%ncid, hs%nodecodenc%nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      att_text = "node_wet_or_dry"
      iret = nf90_put_att(hs%ncid, hs%nodecodenc%nodal_data_id, &
                          'standard_name', trim(att_text))
      call check_err(iret)

!     N O F F
      hs%noffnc%nodal_data_dims(1) = hs%myMesh%num_elems_dim_id
      hs%noffnc%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      call define_netcdf_hotstart_variable(hs, 'noff', hs%noffnc%nodal_data_dims, NF90_INT, &
                                           hs%noffnc%nodal_data_id)
      att_text = "wet or dry state of element where 1 indicates that" &
                 //" the element is wet and 0 indicates that it is dry"
      iret = nf90_put_att(hs%ncid, hs%noffnc%nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      att_text = "element_wet_or_dry"
      iret = nf90_put_att(hs%ncid, hs%noffnc%nodal_data_id, &
                          'standard_name', trim(att_text))
      call check_err(iret)

!     Radiation Stress
      hs%rs1%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%rs1%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      call define_netcdf_hotstart_variable(hs, 'rsx1', hs%rs1%nodal_data_dims, NF90_DOUBLE, &
                                           hs%rs1%u_nodal_data_id)
      call define_netcdf_hotstart_variable(hs, 'rsy1', hs%rs1%nodal_data_dims, NF90_DOUBLE, &
                                           hs%rs1%v_nodal_data_id)

      hs%rs2%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%rs2%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      call define_netcdf_hotstart_variable(hs, 'rsx2', hs%rs2%nodal_data_dims, NF90_DOUBLE, &
                                           hs%rs2%u_nodal_data_id)
      call define_netcdf_hotstart_variable(hs, 'rsy2', hs%rs2%nodal_data_dims, NF90_DOUBLE, &
                                           hs%rs2%v_nodal_data_id)

      att_text = "radiation stress vector x component at previous time step"
      iret = nf90_put_att(hs%ncid, hs%rs1%u_nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      att_text = "radiation stress vector y component at previous time step"
      iret = nf90_put_att(hs%ncid, hs%rs1%v_nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      att_text = "radiation stress vector x component at current time step"
      iret = nf90_put_att(hs%ncid, hs%rs2%u_nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      att_text = "radiation stress vector y component at current time step"
      iret = nf90_put_att(hs%ncid, hs%rs2%v_nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)

      att_text = "radiation_stress_x_component_at_previous_time_step"
      iret = nf90_put_att(hs%ncid, hs%rs1%u_nodal_data_id, &
                          'standard_name', trim(att_text))
      call check_err(iret)
      att_text = "radiation_stress_y_component_at_previous_time_step"
      iret = nf90_put_att(hs%ncid, hs%rs1%v_nodal_data_id, &
                          'standard_name', trim(att_text))
      call check_err(iret)
      att_text = "radiation_stress_x_component_at_current_time_step"
      iret = nf90_put_att(hs%ncid, hs%rs2%u_nodal_data_id, &
                          'standard_name', trim(att_text))
      call check_err(iret)
      att_text = "radiation_stress_y_component_at_current_time_step"
      iret = nf90_put_att(hs%ncid, hs%rs2%v_nodal_data_id, &
                          'standard_name', trim(att_text))
      call check_err(iret)

!     SWAN Radiation Stress
      hs%swan_rs1%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%swan_rs1%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      call define_netcdf_hotstart_variable(hs, 'swan_rsx1', hs%swan_rs1%nodal_data_dims, NF90_DOUBLE, &
                                           hs%swan_rs1%u_nodal_data_id)
      call define_netcdf_hotstart_variable(hs, 'swan_rsy1', hs%swan_rs1%nodal_data_dims, NF90_DOUBLE, &
                                           hs%swan_rs1%v_nodal_data_id)

      hs%swan_rs2%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%swan_rs2%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      call define_netcdf_hotstart_variable(hs, 'swan_rsx2', hs%swan_rs2%nodal_data_dims, NF90_DOUBLE, &
                                           hs%swan_rs2%u_nodal_data_id)
      call define_netcdf_hotstart_variable(hs, 'swan_rsy2', hs%swan_rs2%nodal_data_dims, NF90_DOUBLE, &
                                           hs%swan_rs2%v_nodal_data_id)

      att_text = "radiation stress vector x component at previous time step"
      iret = nf90_put_att(hs%ncid, hs%swan_rs1%u_nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      att_text = "radiation stress vector y component at previous time step"
      iret = nf90_put_att(hs%ncid, hs%swan_rs1%v_nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      att_text = "radiation stress vector x component at current time step"
      iret = nf90_put_att(hs%ncid, hs%swan_rs2%u_nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      att_text = "radiation stress vector y component at current time step"
      iret = nf90_put_att(hs%ncid, hs%swan_rs2%v_nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)

      att_text = "radiation_stress_x_component_at_previous_time_step"
      iret = nf90_put_att(hs%ncid, hs%swan_rs1%u_nodal_data_id, &
                          'standard_name', trim(att_text))
      call check_err(iret)
      att_text = "radiation_stress_y_component_at_previous_time_step"
      iret = nf90_put_att(hs%ncid, hs%swan_rs1%v_nodal_data_id, &
                          'standard_name', trim(att_text))
      call check_err(iret)
      att_text = "radiation_stress_x_component_at_current_time_step"
      iret = nf90_put_att(hs%ncid, hs%swan_rs2%u_nodal_data_id, &
                          'standard_name', trim(att_text))
      call check_err(iret)
      att_text = "radiation_stress_y_component_at_current_time_step"
      iret = nf90_put_att(hs%ncid, hs%swan_rs2%v_nodal_data_id, &
                          'standard_name', trim(att_text))
      call check_err(iret)

      ! Define the fill value to be 0.0
      iret = nf90_put_att(hs%ncid, hs%rs1%u_nodal_data_id, &
                          '_FillValue', 0.0d0)
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%rs1%v_nodal_data_id, &
                          '_FillValue', 0.0d0)
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%rs2%u_nodal_data_id, &
                          '_FillValue', 0.0d0)
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%rs2%v_nodal_data_id, &
                          '_FillValue', 0.0d0)
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%swan_rs1%u_nodal_data_id, &
                          '_FillValue', 0.0d0)
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%swan_rs1%v_nodal_data_id, &
                          '_FillValue', 0.0d0)
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%swan_rs2%u_nodal_data_id, &
                          '_FillValue', 0.0d0)
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%swan_rs2%v_nodal_data_id, &
                          '_FillValue', 0.0d0)
      call check_err(iret)

!     Define hotstart parameters
      call check_err(iret)
      iret = nf90_def_var(hs%ncid, 'imhs', NF90_INT, varid=tempid)
      call check_err(iret)
      att_text = 'model_type'
      iret = nf90_put_att(hs%ncid, tempid, 'long_name', 'model_type')
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name', &
                          'model_type')
      call check_err(iret)
!
      iret = nf90_def_var(hs%ncid, 'iths', NF90_INT, varid=tempid)
      call check_err(iret)
      att_text = &
         'model time step number since the beginning of the model run'
      iret = nf90_put_att(hs%ncid, tempid, 'long_name', trim(att_text))
      call check_err(iret)
      att_text = 'model_time_step'
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name', &
                          trim(att_text))
      call check_err(iret)

      iret = nf90_def_var(hs%ncid, 'nscoue', NF90_INT, varid=tempid)
      call check_err(iret)
      att_text = 'time step counter to determine when the' &
                 //' next entry will be written to the elevation time series at' &
                 //' specified elevation recording Stations output file'
      iret = nf90_put_att(hs%ncid, tempid, 'long_name', trim(att_text))

      call check_err(iret)
      att_text = 'time_step_counter_for_next_entry_elev_rec_stations'
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name', &
                          trim(att_text))
      call check_err(iret)
!
      iret = nf90_def_var(hs%ncid, 'nscouv', NF90_INT, varid=tempid)
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'long_name', &
                          'time step counter to determine when the' &
                          //' next entry will be written to the depth-averaged velocity' &
                          //' time series at specified velocity recording stations output' &
                          //' file.')
      call check_err(iret)

      if ((IM == 10) .or. (IMHS == 10)) then
         iret = nf90_def_var(hs%ncid, 'nscouc', NF90_INT, varid=tempid)
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, tempid, 'long_name', &
                             'time step counter to determine when the next ' &
                             //'entry will be written to the scalar concentration time series' &
                             //' at specified concentration recording stations output file')
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, tempid, 'standard_name', &
                             'time_step_counter_for_next_entry_conc_rec_stations')
         call check_err(iret)
      end if

      iret = nf90_def_var(hs%ncid, 'nscoum', NF90_INT, varid=tempid)
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'long_name', &
                          'time step counter to determine when the' &
                          //' next entry will be written to the atmospheric pressure time' &
                          //' series at specified meteorological recording stations and' &
                          //' wind velocity time series at specified meteorological' &
                          //' recording stations output files')
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name', &
                          'time_step_counter_of_atm_press_and_wind_vel_at_rec_stations')
      call check_err(iret)

      iret = nf90_def_var(hs%ncid, 'nscouge', NF90_INT, varid=tempid)
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'long_name', &
                          'time step counter to determine when the' &
                          //' next entry will be written to the elevation time series at' &
                          //' all nodes in the model grid output file')
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name', &
                          'time_step_counter_of_elev_at_model_nodes')
      call check_err(iret)

      iret = nf90_def_var(hs%ncid, 'nscougv', NF90_INT, varid=tempid)
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'long_name', &
                          'time step counter to determine when the' &
                          //' next entry will be written to the depth-averaged velocity' &
                          //' time series at all nodes in the model grid output file')
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name', &
                          'time_step_counter_of_vel_at_model_nodes')
      call check_err(iret)

      if ((IM == 10) .or. (IMHS == 10)) then
         iret = nf90_def_var(hs%ncid, 'nscougc', NF90_INT, varid=tempid)
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, tempid, 'long_name', &
                             'time step counter to determine when the' &
                             //' next entry will be written to the scalar concentration time' &
                             //' series at All Nodes in the model grid output file')
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, tempid, 'standard_name', &
                             'time_step_counter_of_conc_at_model_nodes')
         call check_err(iret)
      end if

      iret = nf90_def_var(hs%ncid, 'nscougw', NF90_INT, varid=tempid)
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'long_name', &
                          'time step counter to determine when the' &
                          //' next entry will be written to the atmospheric pressure time' &
                          //' series at all nodes in the model grid and wind stress or' &
                          //' velocity time series at all nodes in the model grid output' &
                          //' files')
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, tempid, 'standard_name', &
                          'time_step_counter_of_atm_press_and_wind_vel_at_model_nodes')
      call check_err(iret)

!     define time attributes
      call defineTimeAttributes(hs%ncid, hs%myTime)

!     define metadata and selected fort.15 parameters in netcdf file
      call defineMetaData(hs%ncid)

!     Leave define mode
      iret = nf90_enddef(hs%ncid)
      call check_err(iret)

!     write mesh to netcdf file
      call putMeshVariables(hs%ncid, hs%myMesh)

!     now close the initialized netcdf file
      iret = nf90_close(hs%ncid)
      call check_err(iret)

!-----------------------------------------------------------------------
   end subroutine initNetCDFHotstart
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! @brief Define a netCDF variable for hotstart data with the appropriate storage format
! @param hs Pointer to hotstartData structure
! @param varname Name of the variable to define
! @param nodal_data_dims Dimensions of the variable
! @param nc_data_type NetCDF data type (e.g., NF90_DOUBLE, NF90_INT)
! @param nodal_data_id Output variable ID
!
! This is a helper function to select contiguous storage for netCDF variables
! when using netCDF4 output. Defining the storage as contiguous allows us to
! avoid copy-on-write issues with hdf5 and continuously overwriting data
! which causes the file to grow in size.
!-----------------------------------------------------------------------
   subroutine define_netcdf_hotstart_variable(hs, varname, nodal_data_dims, nc_data_type, nodal_data_id)
      use netcdf_error, only: check_err
      implicit none
      type(hotstartData), pointer, intent(in) :: hs
      character(len=*), intent(in) :: varname
      integer, intent(in) :: nodal_data_dims(:)
      integer, intent(in) :: nc_data_type
      integer, intent(out) :: nodal_data_id

      if (hs%myFile%ncformat == NF90_CLOBBER) then
         call check_err(nf90_def_var(hs%ncid, varname, nc_data_type, &
                                     nodal_data_dims, nodal_data_id))
      else
         call check_err(nf90_def_var(hs%ncid, varname, nc_data_type, &
                                     nodal_data_dims, nodal_data_id, contiguous=.true.))
      end if

   end subroutine define_netcdf_hotstart_variable

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!        I N I T   N E T C D F   H O T S T A R T   H A R M O N I C
!-----------------------------------------------------------------------
!     jgf49.35 Sets up netCDF variables for hotstarting harmonic analysis.
!-----------------------------------------------------------------------
   subroutine initNetCDFHotstartHarmonic(lun, GLOELVDescript, STAELVDescript, STAULVDescript, err)
      use SIZES, only: MNHARF
      use HARM, only: NHASE, NHASV, NHAGE, NHAGV
      use GLOBAL, only: OutputDataDescript_t
      implicit none

      integer, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: GLOELVDescript
      type(OutputDataDescript_t), intent(in) :: STAELVDescript
      type(OutputDataDescript_t), intent(in) :: STAULVDescript
      logical, intent(out) :: err

      integer :: iret ! success or failure of the netcdf call
      character(1024) :: att_text

      LOG_SCOPE_TRACED("initNetCDFHotstartHarmonic", NETCDFIO_TRACING)
      err = .false.

!     Point to the hotstart file we want to work on.
      if (lun == 67) then
         hs => hs67
      else
         hs => hs68
      end if

!     Open existing NetCDF file
      iret = nf90_open(hs%myFile%FILENAME, NF90_WRITE, hs%ncid)
      call check_err(iret)

!     Enter "redefine" mode
      iret = NF90_REDEF(hs%ncid)
      call check_err(iret)

!     Inquire variables
      call getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)

!     Define dimensions
      iret = nf90_def_dim(hs%ncid, 'mnharf', MNHARF, &
                          hs%mnharf_dim_id)
      call check_err(iret)
      iret = nf90_def_dim(hs%ncid, 'mnharfx2', (MNHARF*2), &
                          hs%load_vector_dim_id)
      call check_err(iret)

!     Create station dimension and station name dimension
      iret = nf90_def_dim(hs%ncid, 'namefrlen', hs%namefr_len, &
                          hs%namefr_len_dim_id)
      call check_err(iret)

!     Define harmonic analysis frequency names array
      hs%namefr_dims(1) = hs%namefr_len_dim_id
      hs%namefr_dims(2) = hs%mnharf_dim_id
      iret = nf90_def_var(hs%ncid, 'namefr', NF90_CHAR, &
                          hs%namefr_dims, hs%namefr_id)
      call check_err(iret)

!     harmonic analysis components
      hs%component_dims(1) = hs%mnharf_dim_id
      iret = nf90_def_var(hs%ncid, 'hafreq', NF90_DOUBLE, &
                          hs%component_dims, hs%hafreq_id)
      call check_err(iret)
      att_text = "frequencies (rad/s) of harmonic analysis constituents"
      iret = nf90_put_att(hs%ncid, hs%hafreq_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      att_text = "harmonic analysis frequencies (rad/s)"
      iret = nf90_put_att(hs%ncid, hs%hafreq_id, &
                          'standard_name', trim(att_text))

      iret = nf90_def_var(hs%ncid, 'haff', NF90_DOUBLE, &
                          hs%component_dims, hs%haff_id)
      call check_err(iret)
      att_text = "nodal factors of harmonic analysis constituents"
      iret = nf90_put_att(hs%ncid, hs%haff_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      att_text = "harmonic analysis nodal factors"
      iret = nf90_put_att(hs%ncid, hs%haff_id, &
                          'standard_name', trim(att_text))

      iret = nf90_def_var(hs%ncid, 'haface', NF90_DOUBLE, &
                          hs%component_dims, hs%haface_id)
      call check_err(iret)
      att_text = &
         "equilibrium arguments (degrees) of harmonic analysis " &
         //"constituents"
      iret = nf90_put_att(hs%ncid, hs%haface_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      att_text = "equilibrium arguments (degrees)"
      iret = nf90_put_att(hs%ncid, hs%haface_id, &
                          'standard_name', trim(att_text))

      hs%ha_dims(1) = hs%load_vector_dim_id
      hs%ha_dims(2) = hs%load_vector_dim_id
      iret = nf90_def_var(hs%ncid, 'ha', NF90_DOUBLE, hs%ha_dims, hs%ha_id)
      call check_err(iret)
      att_text = "left hand side matrix for harmonic analysis"
      iret = nf90_put_att(hs%ncid, hs%ha_id, 'long_name', trim(att_text))
      call check_err(iret)
      att_text = "LHS for harmonic analysis"
      iret = nf90_put_att(hs%ncid, hs%ha_id, &
                          'standard_name', trim(att_text))

!     global elevation load vector
      if (NHAGE /= 0) then
         hs%gloelv%nodal_data_dims_3D(1) = hs%load_vector_dim_id
         hs%gloelv%nodal_data_dims_3D(2) = hs%myMesh%num_nodes_dim_id
         hs%gloelv%nodal_data_dims_3D(3) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'gloelv', NF90_DOUBLE, &
                             hs%gloelv%nodal_data_dims_3D, hs%gloelv%nodal_data_id)
         call check_err(iret)
         att_text = "full domain elevation load vector at each node"
         iret = nf90_put_att(hs%ncid, hs%gloelv%nodal_data_id, &
                             'long_name', trim(att_text))
         call check_err(iret)
         att_text = "full domain elevation load vector"
         iret = nf90_put_att(hs%ncid, hs%gloelv%nodal_data_id, &
                             'standard_name', trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%gloelv%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%gloelv%nodal_data_id, &
                             'positive', 'up')
         call check_err(iret)
      end if

!     global velocity load vectors
      if (NHAGV /= 0) then
         hs%glovellv%nodal_data_dims_3D(1) = hs%load_vector_dim_id
         hs%glovellv%nodal_data_dims_3D(2) = hs%myMesh%num_nodes_dim_id
         hs%glovellv%nodal_data_dims_3D(3) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'gloulv', NF90_DOUBLE, &
                             hs%glovellv%nodal_data_dims_3D, hs%glovellv%u_nodal_data_id)
         call check_err(iret)
         att_text = "full domain u velocity load vector at each node"
         iret = nf90_put_att(hs%ncid, hs%glovellv%u_nodal_data_id, &
                             'long_name', trim(att_text))
         call check_err(iret)
         att_text = "full domain u velocity load vector"
         iret = nf90_put_att(hs%ncid, hs%glovellv%u_nodal_data_id, &
                             'standard_name', trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%glovellv%u_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%glovellv%u_nodal_data_id, &
                             'positive', 'east')
         call check_err(iret)

         iret = nf90_def_var(hs%ncid, 'glovlv', NF90_DOUBLE, &
                             hs%glovellv%nodal_data_dims_3D, hs%glovellv%v_nodal_data_id)
         call check_err(iret)
         att_text = "full domain v velocity load vector at each node"
         iret = nf90_put_att(hs%ncid, hs%glovellv%v_nodal_data_id, &
                             'long_name', trim(att_text))
         call check_err(iret)
         att_text = "full domain v velocity load vector"
         iret = nf90_put_att(hs%ncid, hs%glovellv%v_nodal_data_id, &
                             'standard_name', trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%glovellv%v_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%glovellv%v_nodal_data_id, &
                             'positive', 'north')
         call check_err(iret)
      end if

!     station elevation load vector
      if (NHASE /= 0) then
         hs%staelv%num_stations = STAELVDescript%num_fd_records
         iret = nf90_def_dim(hs%ncid, 'elevstation', &
                             hs%staelv%num_stations, hs%staelv%num_sta_dim_id)
         call check_err(iret)
         hs%staelv%station_data_dims_3D(1) = hs%load_vector_dim_id
         hs%staelv%station_data_dims_3D(2) = hs%staelv%num_sta_dim_id
         hs%staelv%station_data_dims_3D(3) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'staelv', NF90_DOUBLE, &
                             hs%staelv%station_data_dims_3D, hs%staelv%station_data_id)
         call check_err(iret)
         att_text = "elevation load vector at each elevation station"
         iret = nf90_put_att(hs%ncid, hs%staelv%station_data_id, &
                             'long_name', trim(att_text))
         call check_err(iret)
         att_text = "station elevation load vector"
         iret = nf90_put_att(hs%ncid, hs%staelv%station_data_id, &
                             'standard_name', trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%staelv%station_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%staelv%station_data_id, &
                             'positive', 'up')
         call check_err(iret)
      end if

!     station velocity load vectors
      if (NHASV /= 0) then
         hs%stavellv%num_stations = STAULVDescript%num_fd_records
         iret = nf90_def_dim(hs%ncid, 'velstation', &
                             hs%stavellv%num_stations, hs%stavellv%num_sta_dim_id)
         call check_err(iret)
!        define dimension
         hs%stavellv%station_data_dims_3D(1) = hs%load_vector_dim_id
         hs%stavellv%station_data_dims_3D(2) = hs%stavellv%num_sta_dim_id
         hs%stavellv%station_data_dims_3D(3) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'staulv', NF90_DOUBLE, &
                             hs%stavellv%station_data_dims_3D, &
                             hs%stavellv%u_station_data_id)
         call check_err(iret)
         att_text = "u velocity load vector at each velocity station"
         iret = nf90_put_att(hs%ncid, hs%stavellv%u_station_data_id, &
                             'long_name', trim(att_text))
         call check_err(iret)
         att_text = "station u velocity load vector"
         iret = nf90_put_att(hs%ncid, hs%stavellv%u_station_data_id, &
                             'standard_name', trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, &
                             hs%stavellv%u_station_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%stavellv%u_station_data_id, &
                             'positive', 'east')
         call check_err(iret)

         iret = nf90_def_var(hs%ncid, 'stavlv', NF90_DOUBLE, &
                             hs%stavellv%station_data_dims_3D, &
                             hs%stavellv%v_station_data_id)
         call check_err(iret)
         att_text = "v velocity load vector at each velocity station"
         iret = nf90_put_att(hs%ncid, hs%stavellv%v_station_data_id, &
                             'long_name', trim(att_text))
         call check_err(iret)
         att_text = "station v velocity load vector"
         iret = nf90_put_att(hs%ncid, hs%stavellv%v_station_data_id, &
                             'standard_name', trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%stavellv%v_station_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%stavellv%v_station_data_id, &
                             'positive', 'north')
         call check_err(iret)
      end if

!     jgf50.44: Automatically turn on compression if we are using the
!     netcdf4 file format.
#ifdef NETCDF_CAN_DEFLATE
      if ((GLOELVDescript%specifier == 5) .or. &
          (GLOELVDescript%specifier == 567) .or. &
          (GLOELVDescript%specifier == 568)) then
         iret = nf90_def_var_deflate(hs%ncid, hs%ha_id, &
                                     1, 1, 2)
         call check_err(iret)
         if (NHAGE /= 0) then
            iret = nf90_def_var_deflate(hs%ncid, &
                                        hs%gloelv%nodal_data_id, 1, 1, 2)
            call check_err(iret)
         end if
         if (NHAGV /= 0) then
            iret = nf90_def_var_deflate(hs%ncid, &
                                        hs%glovellv%u_nodal_data_id, 1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(hs%ncid, &
                                        hs%glovellv%v_nodal_data_id, 1, 1, 2)
            call check_err(iret)
         end if
         if (NHASE /= 0) then
            iret = nf90_def_var_deflate(hs%ncid, &
                                        hs%staelv%station_data_id, 1, 1, 2)
            call check_err(iret)
         end if
         if (NHASV /= 0) then
            iret = nf90_def_var_deflate(hs%ncid, &
                                        hs%stavellv%u_station_data_id, 1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(hs%ncid, &
                                        hs%stavellv%v_station_data_id, 1, 1, 2)
            call check_err(iret)
         end if
      end if
#endif

!     Define harmonic analysis parameters
      call defineParameterWithText(hs%ncid, 'icha', NF90_INT, &
                                   "harmonic analysis spool counter", &
                                   "harmonic analysis spool counter")

      call defineParameterWithText(hs%ncid, 'nz', NF90_INT, &
                                   "set to 0 if a steady harmonic component is included", &
                                   "indicator of steady harmonic component")

      call defineParameterWithText(hs%ncid, 'nf', NF90_INT, &
                                   "set to 1 if a steady harmonic component is included", &
                                   "steady harmonic component number")

      call defineParameterWithText(hs%ncid, 'mm', NF90_INT, &
                                   "2x the number of harmonic frequencies plus any steady component", &
                                   "2x harmonic frequencies [plus 1]")

      call defineParameterWithText(hs%ncid, 'nstae', NF90_INT, &
                                   "number of elevation recording stations for harmonic analysis", &
                                   "number of elevation recording stations")

      call defineParameterWithText(hs%ncid, 'nstav', NF90_INT, &
                                   "number of velocity recording stations for harmonic analysis", &
                                   "number of velocity recording stations")

      call defineParameterWithText(hs%ncid, 'nhase', NF90_INT, &
                                   "indicator for perfomance and formatting of harmonic analysis " &
                                   //"of elevation station data", &
                                   "elevation station harmonic analysis indicator")

      call defineParameterWithText(hs%ncid, 'nhasv', NF90_INT, &
                                   "indicator for perfomance and formatting of harmonic analysis " &
                                   //"of velocity station data", &
                                   "velocity station harmonic analysis indicator")

      call defineParameterWithText(hs%ncid, 'nhage', NF90_INT, &
                                   "indicator for perfomance and formatting of harmonic analysis " &
                                   //"of full domain elevation data (at every node)", &
                                   "full domain elevation harmonic analysis indicator")

      call defineParameterWithText(hs%ncid, 'nhagv', NF90_INT, &
                                   "indicator for perfomance and formatting of harmonic analysis " &
                                   //"of full domain velocity data (at every node)", &
                                   "full domain velocity harmonic analysis indicator")

      call defineParameterWithText(hs%ncid, 'icall', NF90_INT, &
                                   "number of subroutine calls to update load vectors and left " &
                                   //"matrix for harmonic analysis", &
                                   "number of calls to update harmonic analysis")

      call defineParameterWithText(hs%ncid, 'nfreq', NF90_INT, &
                                   "number of frequencies under consideration in harmonic analysis" &
                                   //" not including a steady component, if any", &
                                   "number of frequencies for harmonic analysis")

      call defineParameterWithText(hs%ncid, 'timeud', NF90_DOUBLE, &
                                   "ADCIRC time at the most recent update of the load vectors for " &
                                   //"harmonic analysis", &
                                   "update time for load vectors")

      call defineParameterWithText(hs%ncid, 'itud', NF90_INT, &
                                   "ADCIRC time step at the most recent update of the load vectors" &
                                   //" for harmonic analysis", &
                                   "update time step for load vectors")

!     Leave define mode
      iret = nf90_enddef(hs%ncid)
      call check_err(iret)

!     now close the initialized netcdf file
      iret = nf90_close(hs%ncid)
      call check_err(iret)

!-----------------------------------------------------------------------
   end subroutine initNetCDFHotstartHarmonic
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!        I N I T   N E T C D F   H O T S T A R T   H A R M O N I C
!                                  M E A N S   V A R I A N C E S
!-----------------------------------------------------------------------
!     jgf49.43.14 Sets up netCDF variables for hotstarting harmonic analysis
!     means and variances calculations.
!-----------------------------------------------------------------------
   subroutine initNetCDFHotstartHarmonicMeansVariances(lun, ELAVDescript, reterror)
      use HARM, only: NHAGE, NHAGV
      use GLOBAL, only: OutputDataDescript_t
      implicit none

      integer, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: ELAVDescript
      logical, intent(out) :: reterror

      integer :: iret ! success or failure of the netcdf call
      character(1024) :: att_text

      LOG_SCOPE_TRACED("initNetCDFHotstartHarmonicMeansVariances", NETCDFIO_TRACING)
      reterror = .false.

!     Point to the hotstart file we want to work on. Memory allocation
!     was already done for means and variances by initNetCDFHotstartHarmonic.
      if (lun == 67) then
         hs => hs67
      else
         hs => hs68
      end if

!     Open existing NetCDF file
      iret = nf90_open(hs%myFile%FILENAME, NF90_WRITE, hs%ncid)
      call check_err(iret)

!     Enter "redefine" mode
      iret = nf90_redef(hs%ncid)
      call check_err(iret)

!     Inquire variables
      call getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)

!     elevation
      if (NHAGE /= 0) then
         ! ELAV
         hs%elav%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
         hs%elav%nodal_data_dims(2) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'elav', NF90_DOUBLE, &
                             hs%elav%nodal_data_dims, hs%elav%nodal_data_id)
         call check_err(iret)
         att_text = "sum of elevations computed by ADCIRC, at every " &
                    //"node in the model grid, over all time steps since harmonic " &
                    //"analysis means and variance checking has begun"
         iret = nf90_put_att(hs%ncid, hs%elav%nodal_data_id, &
                             'long_name', trim(att_text))
         call check_err(iret)
         att_text = "sum of elevations"
         iret = nf90_put_att(hs%ncid, hs%elav%nodal_data_id, &
                             'standard_name', trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%elav%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%elav%nodal_data_id, &
                             'positive', 'up')
         call check_err(iret)
         ! ELVA
         hs%elva%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
         hs%elva%nodal_data_dims(2) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'elva', NF90_DOUBLE, &
                             hs%elva%nodal_data_dims, hs%elva%nodal_data_id)
         call check_err(iret)
         att_text = "sum of squares of elevations computed by ADCIRC, " &
                    //"at every node in the model grid, over all time steps since " &
                    //"harmonic analysis means and variance checking has begun"
         iret = nf90_put_att(hs%ncid, hs%elva%nodal_data_id, &
                             'long_name', trim(att_text))
         call check_err(iret)
         att_text = "sum of squares of elevations"
         iret = nf90_put_att(hs%ncid, hs%elva%nodal_data_id, &
                             'standard_name', trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%elva%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%elva%nodal_data_id, &
                             'positive', 'up')
         call check_err(iret)
      end if

!     global velocity load vectors
      if (NHAGV /= 0) then
         ! XVELAV
         hs%xvelav%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
         hs%xvelav%nodal_data_dims(2) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'xvelav', NF90_DOUBLE, &
                             hs%xvelav%nodal_data_dims, &
                             hs%xvelav%nodal_data_id)
         call check_err(iret)
         att_text = "sum of depth-averaged u velocities computed by " &
                    //"ADCIRC, at every node in the model grid, over all time steps " &
                    //"since harmonic analysis means and variance checking has begun"
         iret = nf90_put_att(hs%ncid, hs%xvelav%nodal_data_id, &
                             'long_name', trim(att_text))
         call check_err(iret)
         att_text = "sum of depth averaged u velocities"
         iret = nf90_put_att(hs%ncid, hs%xvelav%nodal_data_id, &
                             'standard_name', trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%xvelav%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%xvelav%nodal_data_id, &
                             'positive', 'east')
         call check_err(iret)
         ! YVELAV
         hs%yvelav%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
         hs%yvelav%nodal_data_dims(2) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'yvelav', NF90_DOUBLE, &
                             hs%yvelav%nodal_data_dims, hs%yvelav%nodal_data_id)
         call check_err(iret)
         att_text = "sum of depth-averaged v velocities computed by " &
                    //"ADCIRC, at every node in the model grid, over all time steps " &
                    //"since harmonic analysis means and variance checking has begun"
         iret = nf90_put_att(hs%ncid, hs%yvelav%nodal_data_id, &
                             'long_name', trim(att_text))
         call check_err(iret)
         att_text = "sum of depth averaged v velocities"
         iret = nf90_put_att(hs%ncid, hs%yvelav%nodal_data_id, &
                             'standard_name', trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%yvelav%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%yvelav%nodal_data_id, &
                             'positive', 'north')
         call check_err(iret)

         ! XVELVA
         hs%xvelva%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
         hs%xvelva%nodal_data_dims(2) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'xvelva', NF90_DOUBLE, &
                             hs%xvelva%nodal_data_dims, &
                             hs%xvelva%nodal_data_id)
         call check_err(iret)
         att_text = "sum of squares of depth averaged u velocities " &
                    //"computed by ADCIRC, at every node in the model grid, over all" &
                    //" time steps since harmonic analysis means and variance " &
                    //"checking has begun"
         iret = nf90_put_att(hs%ncid, hs%xvelva%nodal_data_id, &
                             'long_name', trim(att_text))
         call check_err(iret)
         att_text = "sum of squares of depth averaged u velocities"
         iret = nf90_put_att(hs%ncid, hs%xvelva%nodal_data_id, &
                             'standard_name', trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%xvelva%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%xvelva%nodal_data_id, &
                             'positive', 'east')
         call check_err(iret)
         ! YVELVA
         hs%yvelva%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
         hs%yvelva%nodal_data_dims(2) = hs%myTime%timenc_dim_id
         iret = nf90_def_var(hs%ncid, 'yvelva', NF90_DOUBLE, &
                             hs%yvelva%nodal_data_dims, hs%yvelva%nodal_data_id)
         call check_err(iret)
         att_text = "sum of squares of depth averaged v velocities " &
                    //"computed by ADCIRC, at every node in the model grid, over " &
                    //"all time steps since harmonic analysis means and variance " &
                    //"checking has begun"
         iret = nf90_put_att(hs%ncid, hs%yvelva%nodal_data_id, &
                             'long_name', trim(att_text))
         call check_err(iret)
         att_text = "sum of squares of depth averaged v velocities"
         iret = nf90_put_att(hs%ncid, hs%yvelva%nodal_data_id, &
                             'standard_name', trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%yvelva%nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%yvelva%nodal_data_id, &
                             'positive', 'up')
         call check_err(iret)
      end if

!     jgf50.44: Automatically turn on compression if we are using the
!     netcdf4 file format.
#ifdef NETCDF_CAN_DEFLATE
      if ((ELAVDescript%specifier == 5) .or. &
          (ELAVDescript%specifier == 567) .or. &
          (ELAVDescript%specifier == 568)) then
         if (NHAGE /= 0) then
            iret = nf90_def_var_deflate(hs%ncid, &
                                        hs%elav%nodal_data_id, 1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(hs%ncid, &
                                        hs%elva%nodal_data_id, 1, 1, 2)
            call check_err(iret)
         end if
         if (NHAGV /= 0) then
            iret = nf90_def_var_deflate(hs%ncid, &
                                        hs%xvelav%nodal_data_id, 1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(hs%ncid, &
                                        hs%xvelva%nodal_data_id, 1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(hs%ncid, &
                                        hs%yvelav%nodal_data_id, 1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(hs%ncid, &
                                        hs%yvelva%nodal_data_id, 1, 1, 2)
            call check_err(iret)
         end if
      end if
#endif

!     Define harmonic analysis parameters
      call defineParameterWithText(hs%ncid, 'ntsteps', NF90_INT, &
                                   "number of time steps since start of means and variance", &
                                   "number of time steps since start of means and variance")

!     Leave define mode
      iret = nf90_enddef(hs%ncid)
      call check_err(iret)

!     now close the initialized netcdf file
      iret = nf90_close(hs%ncid)
      call check_err(iret)

!-----------------------------------------------------------------------
   end subroutine initNetCDFHotstartHarmonicMeansVariances
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!        I N I T   N E T C D F   H O T S T A R T   3 D
!-----------------------------------------------------------------------
!     jgf49.49.02 Sets up netCDF variables for hotstarting 3D.
!-----------------------------------------------------------------------
   subroutine initNetCDFHotstart3D(lun, netcdf_format)
      use SIZES, only: MNPROC
      use GLOBAL, only: OutputDataDescript_t, IDEN
      use GLOBAL_3DVS, only: NFEN
      implicit none

      integer, intent(in) :: lun
      integer, intent(in) :: netcdf_format ! whether netcdf3 or netcdf4
      ! format (netcdf4=hdf5)
      ! classic data model in any case

      logical :: err

      integer :: iret ! success or failure of the netcdf call
      character(1024) :: att_text
      integer :: tempid

      LOG_SCOPE_TRACED("initNetCDFHotstart3D", NETCDFIO_TRACING)
      err = .false.

!     Point to the hotstart file we want to work on.
      if (lun == 67) then
         hs => hs67
      else
         hs => hs68
      end if

!     nodes and elements
      hs%myMesh => adcircMesh
!     jgf50.06: If this is called during a parallel run, and all other
!     fulldomain 3D output is turned off, the number of vertical nodes will not
!     have been set.
      if (MNPROC > 1) then
         hs%myMesh%num_v_nodes = NFEN
      end if

!     Open existing NetCDF file
      iret = nf90_open(hs%myFile%FILENAME, NF90_WRITE, hs%ncid)
      call check_err(iret)

!     Enter "redefine" mode
      iret = nf90_redef(hs%ncid)
      call check_err(iret)

!     Inquire variables
      call getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)

!     DUU
      hs%duu%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%duu%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'duu', NF90_DOUBLE, &
                          hs%duu%nodal_data_dims, hs%duu%nodal_data_id)
      call check_err(iret)
      att_text = "velocity dispersion term"
      iret = nf90_put_att(hs%ncid, hs%duu%nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%duu%nodal_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)

!     DUV
      hs%duv%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%duv%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'duv', NF90_DOUBLE, &
                          hs%duv%nodal_data_dims, hs%duv%nodal_data_id)
      call check_err(iret)
      att_text = "velocity dispersion term"
      iret = nf90_put_att(hs%ncid, hs%duv%nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%duv%nodal_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)

!     DVV
      hs%dvv%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%dvv%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'dvv', NF90_DOUBLE, &
                          hs%dvv%nodal_data_dims, hs%dvv%nodal_data_id)
      call check_err(iret)
      att_text = "velocity dispersion term"
      iret = nf90_put_att(hs%ncid, hs%dvv%nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%dvv%nodal_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)

!     UU
      hs%uu%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%uu%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'uu', NF90_DOUBLE, &
                          hs%uu%nodal_data_dims, hs%uu%nodal_data_id)
      call check_err(iret)
      att_text = "vertically averaged velocity in east direction"
      iret = nf90_put_att(hs%ncid, hs%uu%nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%uu%nodal_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)

!     VV
      hs%vv%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%vv%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'vv', NF90_DOUBLE, &
                          hs%vv%nodal_data_dims, hs%vv%nodal_data_id)
      call check_err(iret)
      att_text = "vertically averaged velocity in north direction"
      iret = nf90_put_att(hs%ncid, hs%vv%nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%vv%nodal_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)

!     BSX
      hs%bsx%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%bsx%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'bsx', NF90_DOUBLE, &
                          hs%bsx%nodal_data_dims, hs%bsx%nodal_data_id)
      call check_err(iret)
      att_text = "bottom stress in east direction"
      iret = nf90_put_att(hs%ncid, hs%bsx%nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%bsx%nodal_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)

!     BSY
      hs%bsy%nodal_data_dims(1) = hs%myMesh%num_nodes_dim_id
      hs%bsy%nodal_data_dims(2) = hs%myTime%timenc_dim_id
      iret = nf90_def_var(hs%ncid, 'bsy', NF90_DOUBLE, &
                          hs%bsy%nodal_data_dims, hs%bsy%nodal_data_id)
      call check_err(iret)
      att_text = "bottom stress in north direction"
      iret = nf90_put_att(hs%ncid, hs%bsy%nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%bsy%nodal_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)

!     3D DENSITY
      if (IDEN /= 0) then
         hs%density3D%nodal_data_dims_3D(1) = hs%myMesh%num_nodes_dim_id
         hs%density3D%nodal_data_dims_3D(2) = hs%myMesh%num_v_nodes_dim_id
         hs%density3D%nodal_data_dims_3D(3) = hs%myTime%timenc_dim_id
      end if
      if (IDEN == 1) then
         iret = nf90_def_var(hs%ncid, 'sigt', NF90_DOUBLE, &
                             hs%density3D%nodal_data_dims_3D, hs%density3D%u_nodal_data_id)
         call check_err(iret)
         att_text = "sigma t density"
         iret = nf90_put_att(hs%ncid, hs%density3D%u_nodal_data_id, &
                             'long_name', trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%density3D%u_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
      end if
      if ((IDEN == 2) .or. (IDEN == 4)) then
         iret = nf90_def_var(hs%ncid, 'salinity', NF90_DOUBLE, &
                             hs%density3D%nodal_data_dims_3D, &
                             hs%density3D%v_nodal_data_id)
         call check_err(iret)
         att_text = "salinity"
         iret = nf90_put_att(hs%ncid, hs%density3D%v_nodal_data_id, &
                             'long_name', trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%density3D%v_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
      end if
      if ((IDEN == 3) .or. (IDEN == 4)) then
         iret = nf90_def_var(hs%ncid, 'temperature', NF90_DOUBLE, &
                             hs%density3D%nodal_data_dims_3D, &
                             hs%density3D%w_nodal_data_id)
         call check_err(iret)
         att_text = "salinity"
         iret = nf90_put_att(hs%ncid, hs%density3D%w_nodal_data_id, &
                             'long_name', trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(hs%ncid, hs%density3D%w_nodal_data_id, &
                             '_FillValue', doubleval)
         call check_err(iret)
      end if

!     3D VELOCITY
      hs%velocity3D%nodal_data_dims_3D(1) = hs%myMesh%num_nodes_dim_id
      hs%velocity3D%nodal_data_dims_3D(2) = hs%myMesh%num_v_nodes_dim_id
      hs%velocity3D%nodal_data_dims_3D(3) = hs%myTime%timenc_dim_id
!     u-vel3D
      iret = nf90_def_var(hs%ncid, 'u-vel3D', NF90_DOUBLE, &
                          hs%velocity3D%nodal_data_dims_3D, &
                          hs%velocity3D%u_nodal_data_id)
      call check_err(iret)
      att_text = "3D fulldomain velocity in east direction"
      iret = nf90_put_att(hs%ncid, hs%velocity3D%u_nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%velocity3D%u_nodal_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)
!     v-vel3D
      iret = nf90_def_var(hs%ncid, 'v-vel3D', NF90_DOUBLE, &
                          hs%velocity3D%nodal_data_dims_3D, &
                          hs%velocity3D%v_nodal_data_id)
      call check_err(iret)
      att_text = "3D fulldomain velocity in north direction"
      iret = nf90_put_att(hs%ncid, hs%velocity3D%v_nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%velocity3D%v_nodal_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)
!     w-vel3D
      iret = nf90_def_var(hs%ncid, 'w-vel3D', NF90_DOUBLE, &
                          hs%velocity3D%nodal_data_dims_3D, &
                          hs%velocity3D%w_nodal_data_id)
      call check_err(iret)
      att_text = "3D full domain velocity in the vertical direction"
      iret = nf90_put_att(hs%ncid, hs%velocity3D%w_nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%velocity3D%w_nodal_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)

!     3D TURBULENCE
      hs%turbulence3D%nodal_data_dims_3D(1) = hs%myMesh%num_nodes_dim_id
      hs%turbulence3D%nodal_data_dims_3D(2) = hs%myMesh%num_v_nodes_dim_id
      hs%turbulence3D%nodal_data_dims_3D(3) = hs%myTime%timenc_dim_id
!     Q20
      iret = nf90_def_var(hs%ncid, 'q20', NF90_DOUBLE, &
                          hs%turbulence3D%nodal_data_dims_3D, &
                          hs%turbulence3D%u_nodal_data_id)
      call check_err(iret)
      att_text = "3D fulldomain turbulence kinetic energy"
      iret = nf90_put_att(hs%ncid, hs%turbulence3D%u_nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%turbulence3D%u_nodal_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)
!     L
      iret = nf90_def_var(hs%ncid, 'l', NF90_DOUBLE, &
                          hs%turbulence3D%nodal_data_dims_3D, &
                          hs%turbulence3D%v_nodal_data_id)
      call check_err(iret)
      att_text = "3D fulldomain turbulence length scale"
      iret = nf90_put_att(hs%ncid, hs%turbulence3D%v_nodal_data_id, &
                          'long_name', trim(att_text))
      call check_err(iret)
      iret = nf90_put_att(hs%ncid, hs%turbulence3D%v_nodal_data_id, &
                          '_FillValue', doubleval)
      call check_err(iret)

!     jgf50.44: Automatically turn on compression if we are using the
!     netcdf4 file format.
#ifdef NETCDF_CAN_DEFLATE
      if ((netcdf_format == 5) .or. &
          (netcdf_format == 567) .or. &
          (netcdf_format == 568)) then
         iret = nf90_def_var_deflate(hs%ncid, hs%duu%nodal_data_id, &
                                     1, 1, 2)
         call check_err(iret)
         iret = nf90_def_var_deflate(hs%ncid, hs%dvv%nodal_data_id, &
                                     1, 1, 2)
         call check_err(iret)
         iret = nf90_def_var_deflate(hs%ncid, hs%uu%nodal_data_id, &
                                     1, 1, 2)
         call check_err(iret)
         iret = nf90_def_var_deflate(hs%ncid, hs%vv%nodal_data_id, &
                                     1, 1, 2)
         call check_err(iret)
         iret = nf90_def_var_deflate(hs%ncid, hs%bsx%nodal_data_id, &
                                     1, 1, 2)
         call check_err(iret)
         iret = nf90_def_var_deflate(hs%ncid, hs%bsy%nodal_data_id, &
                                     1, 1, 2)
         call check_err(iret)
         if (IDEN == 1) then
            iret = nf90_def_var_deflate(hs%ncid, &
                                        hs%density3D%u_nodal_data_id, 1, 1, 2)
            call check_err(iret)
         end if
         if ((IDEN == 2) .or. (IDEN == 4)) then
            iret = nf90_def_var_deflate(hs%ncid, &
                                        hs%density3D%v_nodal_data_id, 1, 1, 2)
            call check_err(iret)
         end if
         if ((IDEN == 3) .or. (IDEN == 4)) then
            iret = nf90_def_var_deflate(hs%ncid, &
                                        hs%density3D%w_nodal_data_id, 1, 1, 2)
            call check_err(iret)
         end if
         iret = nf90_def_var_deflate(hs%ncid, &
                                     hs%velocity3D%u_nodal_data_id, 1, 1, 2)
         call check_err(iret)
         iret = nf90_def_var_deflate(hs%ncid, &
                                     hs%velocity3D%v_nodal_data_id, 1, 1, 2)
         call check_err(iret)
         iret = nf90_def_var_deflate(hs%ncid, &
                                     hs%velocity3D%w_nodal_data_id, 1, 1, 2)
         call check_err(iret)
         iret = nf90_def_var_deflate(hs%ncid, &
                                     hs%turbulence3D%u_nodal_data_id, 1, 1, 2)
         call check_err(iret)
         iret = nf90_def_var_deflate(hs%ncid, &
                                     hs%turbulence3D%v_nodal_data_id, 1, 1, 2)
         call check_err(iret)

      end if
#endif

      iret = nf90_def_var(hs%ncid, 'n3dsd', NF90_INT, varid=tempid)
      call check_err(iret)
      iret = nf90_def_var(hs%ncid, 'i3dsdrec', NF90_INT, varid=tempid)
      call check_err(iret)
      iret = nf90_def_var(hs%ncid, 'n3dsv', NF90_INT, varid=tempid)
      call check_err(iret)
      iret = nf90_def_var(hs%ncid, 'i3dsvrec', NF90_INT, varid=tempid)
      call check_err(iret)
      iret = nf90_def_var(hs%ncid, 'n3dst', NF90_INT, varid=tempid)
      call check_err(iret)
      iret = nf90_def_var(hs%ncid, 'i3dstrec', NF90_INT, varid=tempid)
      call check_err(iret)
      iret = nf90_def_var(hs%ncid, 'n3dgd', NF90_INT, varid=tempid)
      call check_err(iret)
      iret = nf90_def_var(hs%ncid, 'i3dgdrec', NF90_INT, varid=tempid)
      call check_err(iret)
      iret = nf90_def_var(hs%ncid, 'n3dgv', NF90_INT, varid=tempid)
      call check_err(iret)
      iret = nf90_def_var(hs%ncid, 'i3dgvrec', NF90_INT, varid=tempid)
      call check_err(iret)
      iret = nf90_def_var(hs%ncid, 'n3dgt', NF90_INT, varid=tempid)
      call check_err(iret)
      iret = nf90_def_var(hs%ncid, 'i3dgtrec', NF90_INT, varid=tempid)
      call check_err(iret)

!     Leave define mode
      iret = nf90_enddef(hs%ncid)
      call check_err(iret)

!     now close the initialized netcdf file
      iret = nf90_close(hs%ncid)
      call check_err(iret)

!-----------------------------------------------------------------------
   end subroutine initNetCDFHotstart3D
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!        D E F I N E   P A R A M E T E R   W I T H   T E X T
!-----------------------------------------------------------------------
!     jgf49.44 Defines a variable in the netcdf file and associates
!     attribute text with it.
!-----------------------------------------------------------------------
   subroutine defineParameterWithText(ncid, param, varType, &
                                      longName, standardName)
      implicit none
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: param
      integer, intent(in) :: varType ! netcdf-defined data type
      character(len=*), intent(in) :: longName
      character(len=*), intent(in) :: standardName

      integer :: tempid ! variable id for attaching to text
      integer :: iret ! netcdf err indicator

      LOG_SCOPE_TRACED("defineParameterWithText", NETCDFIO_TRACING)
      iret = nf90_def_var(ncid, param, varType, varid=tempid)
      call check_err(iret)
      iret = nf90_put_att(ncid, tempid, 'long_name', trim(longName))
      call check_err(iret)
      iret = nf90_put_att(ncid, tempid, 'standard_name', &
                          trim(standardName))
      call check_err(iret)

!-----------------------------------------------------------------------
   end subroutine defineParameterWithText
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E    W R I T E   N E T C D F   H O T S T A R T
!-----------------------------------------------------------------------
!     jgf49.17.02 Writes data to hotstart file.
!-----------------------------------------------------------------------
   subroutine writeNetCDFHotstart(lun, Elev1Descript, &
                                  Elev2Descript, VelDescript, CH1Descript, EtaDiscDescript, &
                                  NodeCodeDescript, NOFFDescript, H1_HS_Descript, &
                                  H2_HS_Descript, RS1Descript, RS2Descript, SWANRS1Descript, &
                                  SWANRS2Descript, timesec, it)
      use SIZES, only: globaldir, mnproc
      use GLOBAL, only: OutputDataDescript_t, im, nscoue, imhs, &
                        nscouv, nscouc, nscoum, &
                        nscouge, nscougv, nscougc, nscougw, nrs
      implicit none

      integer, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: Elev1Descript
      type(OutputDataDescript_t), intent(in) :: Elev2Descript
      type(OutputDataDescript_t), intent(in) :: VelDescript
      type(OutputDataDescript_t), intent(in) :: CH1Descript
      type(OutputDataDescript_t), intent(in) :: EtaDiscDescript
      type(OutputDataDescript_t), intent(in) :: NodeCodeDescript
      type(OutputDataDescript_t), intent(in) :: NOFFDescript
      type(OutputDataDescript_t), intent(in) :: H1_HS_Descript
      type(OutputDataDescript_t), intent(in) :: H2_HS_Descript
      type(OutputDataDescript_t), intent(in) :: RS1Descript
      type(OutputDataDescript_t), intent(in) :: RS2Descript
      type(OutputDataDescript_t), intent(in) :: SWANRS1Descript
      type(OutputDataDescript_t), intent(in) :: SWANRS2Descript

      real(8), intent(in) :: timesec
      integer, intent(in) :: it ! current ADCIRC time step

      integer :: counti(1), starti(1)
      integer :: kount(2), start(2) ! for nodally based data
      integer :: elekount(2) ! for elementally based data
      integer :: iret ! success or failure of the netcdf call
      integer :: tempid
      character(len=10) :: fext

      LOG_SCOPE_TRACED("writeNetCDFHotstart", NETCDFIO_TRACING)

!     Point to the hotstart file we want to work on
      if (lun == 67) then
         hs => hs67
      else
         hs => hs68
      end if
      hs%myMesh => adcircMesh

      ! create file name
      write (fext, '(i0)') lun
      hs%myFile%filename = trim(globaldir)//'/fort.'//trim(fext)//'.nc'
      call logMessage(INFO, 'Opening hotstart file "' &
                      //trim(hs%myFile%filename)//'" for writing.')

      iret = nf90_open(hs%myFile%FILENAME, NF90_WRITE, hs%ncid)
      call check_err(iret)

!     Inquire variables
      call getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)
      if (hs%myTime%initialized .eqv. .false.) then
         allocate (hs%myTime%timenc(hs%myTime%timenc_len))
         hs%myTime%initialized = .true.
      end if

!     Don't increment the record counter so that we can store data at the
!     next location in the netcdf file -- the hotstart file is only intended
!     to have a single snapshot of data in it.
!     hs%myFile%record_counter = hs%myFile%record_counter + 1
      hs%myFile%record_counter = 1

!     Store time
      iret = nf90_inq_varid(hs%ncid, "time", hs%myTime%timenc_id)
      starti(1) = hs%myFile%record_counter
      counti(1) = hs%myTime%timenc_len
      hs%myTime%timenc(hs%myTime%timenc_len) = timesec
      iret = nf90_put_var(hs%ncid, hs%myTime%timenc_id, hs%myTime%timenc, &
                          starti, counti)
      call check_err(iret)

      kount(1) = hs%myMesh%num_nodes
      kount(2) = hs%myTime%timenc_len
      elekount(1) = hs%myMesh%num_elems
      elekount(2) = kount(2)
      start(1) = 1
      start(2) = hs%myFile%record_counter

!     Get the NetCDF IDs of the relevant variables from the file
      iret = nf90_inq_varid(hs%ncid, "zeta1", hs%zeta1%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "zeta2", hs%zeta2%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "zetad", hs%zetad%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "u-vel", hs%vel%u_nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "v-vel", hs%vel%v_nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nodecode", hs%nodecodenc%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "noff", hs%noffnc%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "h1", hs%htot1%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "h2", hs%htot2%nodal_data_id)
      call check_err(iret)
      if (im == 10 .or. imhs == 10) then
         iret = nf90_inq_varid(hs%ncid, "ch1", hs%ch1%nodal_data_id)
         call check_err(iret)
      end if
      if (nrs /= 0) then
         iret = nf90_inq_varid(hs%ncid, "rsx1", hs%rs1%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(hs%ncid, "rsy1", hs%rs1%v_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(hs%ncid, "rsx2", hs%rs2%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(hs%ncid, "rsy2", hs%rs2%v_nodal_data_id)
         call check_err(iret)
#ifdef CSWAN
         if (nrs == 3) then
            iret = nf90_inq_varid(hs%ncid, "swan_rsx1", hs%swan_rs1%u_nodal_data_id)
            call check_err(iret)
            iret = nf90_inq_varid(hs%ncid, "swan_rsy1", hs%swan_rs1%v_nodal_data_id)
            call check_err(iret)
            iret = nf90_inq_varid(hs%ncid, "swan_rsx2", hs%swan_rs2%u_nodal_data_id)
            call check_err(iret)
            iret = nf90_inq_varid(hs%ncid, "swan_rsy2", hs%swan_rs2%v_nodal_data_id)
            call check_err(iret)
         end if
#endif
      end if
!     Write the nodal data to the netcdf file
      if (MNPROC == 1) then
         iret = nf90_put_var(hs%ncid, hs%zeta1%nodal_data_id, &
                             Elev1Descript%array, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%zeta2%nodal_data_id, &
                             Elev2Descript%array, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%zetad%nodal_data_id, &
                             EtaDiscDescript%array, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%vel%u_nodal_data_id, &
                             VelDescript%array, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%vel%v_nodal_data_id, &
                             VelDescript%array2, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%nodecodenc%nodal_data_id, &
                             NodeCodeDescript%iarray, start, kount)
         call check_err(iret)

         if (nrs /= 0) then
            iret = nf90_put_var(hs%ncid, hs%rs1%u_nodal_data_id, &
                                rs1Descript%array, start, kount)
            call check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%rs1%v_nodal_data_id, &
                                rs1Descript%array2, start, kount)
            call check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%rs2%u_nodal_data_id, &
                                rs2Descript%array, start, kount)
            call check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%rs2%v_nodal_data_id, &
                                rs2Descript%array2, start, kount)
            call check_err(iret)
            if (nrs == 3) then
               iret = nf90_put_var(hs%ncid, hs%swan_rs1%u_nodal_data_id, &
                                   SwanRs1Descript%array, start, kount)
               call check_err(iret)
               iret = nf90_put_var(hs%ncid, hs%swan_rs1%v_nodal_data_id, &
                                   SwanRs1Descript%array2, start, kount)
               call check_err(iret)
               iret = nf90_put_var(hs%ncid, hs%swan_rs2%u_nodal_data_id, &
                                   SwanRs2Descript%array, start, kount)
               call check_err(iret)
               iret = nf90_put_var(hs%ncid, hs%swan_rs2%v_nodal_data_id, &
                                   SwanRs2Descript%array2, start, kount)
               call check_err(iret)
            end if
         end if

         !------ H1, H2
         iret = nf90_put_var(hs%ncid, hs%htot2%nodal_data_id, &
                             H2_HS_Descript%array, start, kount)
         call check_err(iret)

         iret = nf90_put_var(hs%ncid, hs%htot1%nodal_data_id, &
                             H1_HS_Descript%array, start, kount)
         call check_err(iret)

         iret = nf90_put_var(hs%ncid, hs%noffnc%nodal_data_id, &
                             NOFFDescript%iarray, start, elekount)
         call check_err(iret)

         if (im == 10 .or. imhs == 10) then
            iret = nf90_put_var(hs%ncid, hs%ch1%nodal_data_id, &
                                CH1Descript%array, start, kount)
            call check_err(iret)
         end if

      else
         iret = nf90_put_var(hs%ncid, hs%zeta1%nodal_data_id, &
                             Elev1Descript%array_g, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%zeta2%nodal_data_id, &
                             Elev2Descript%array_g, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%zetad%nodal_data_id, &
                             EtaDiscDescript%array_g, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%vel%u_nodal_data_id, &
                             VelDescript%array_g, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%vel%v_nodal_data_id, &
                             VelDescript%array2_g, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%nodecodenc%nodal_data_id, &
                             NodeCodeDescript%iarray_g, start, kount)
         call check_err(iret)

         if (nrs /= 0) then
            iret = nf90_put_var(hs%ncid, hs%rs1%u_nodal_data_id, &
                                rs1Descript%array_g, start, kount)
            call check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%rs1%v_nodal_data_id, &
                                rs1Descript%array2_g, start, kount)
            call check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%rs2%u_nodal_data_id, &
                                rs2Descript%array_g, start, kount)
            call check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%rs2%v_nodal_data_id, &
                                rs2Descript%array2_g, start, kount)
            call check_err(iret)
#ifdef CSWAN
            if (nrs == 3) then
               iret = nf90_put_var(hs%ncid, hs%swan_rs1%u_nodal_data_id, &
                                   SwanRs1Descript%array_g, start, kount)
               call check_err(iret)
               iret = nf90_put_var(hs%ncid, hs%swan_rs1%v_nodal_data_id, &
                                   SwanRs1Descript%array2_g, start, kount)
               call check_err(iret)
               iret = nf90_put_var(hs%ncid, hs%swan_rs2%u_nodal_data_id, &
                                   SwanRs2Descript%array_g, start, kount)
               call check_err(iret)
               iret = nf90_put_var(hs%ncid, hs%swan_rs2%v_nodal_data_id, &
                                   SwanRs2Descript%array2_g, start, kount)
               call check_err(iret)
            end if
#endif
         end if

         ! H1, H2
         iret = nf90_put_var(hs%ncid, hs%htot2%nodal_data_id, &
                             H2_HS_Descript%array_g, start, kount)
         call check_err(iret)

         iret = nf90_put_var(hs%ncid, hs%htot1%nodal_data_id, &
                             H1_HS_Descript%array_g, start, kount)
         call check_err(iret)

         iret = nf90_put_var(hs%ncid, hs%noffnc%nodal_data_id, &
                             NOFFDescript%iarray_g, start, elekount)
         call check_err(iret)

         if (im == 10 .or. imhs == 10) then
            iret = nf90_put_var(hs%ncid, hs%ch1%nodal_data_id, &
                                CH1Descript%array_g, start, kount)
            call check_err(iret)
         end if

      end if

!     Get each variable ID for the model parameters in the netcdf file
!     and immediately write the parameter value to that variable ID before
!     going on to the next one.
      iret = nf90_inq_varid(hs%ncid, "imhs", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, im)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "iths", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, it)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nscoue", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nscoue)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nscouv", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nscouv)
      call check_err(iret)
      if ((IM == 10) .or. (IMHS == 10)) then
         iret = nf90_inq_varid(hs%ncid, "nscouc", tempid)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, tempid, nscouc)
         call check_err(iret)
      end if
      iret = nf90_inq_varid(hs%ncid, "nscoum", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nscoum)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nscouge", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nscouge)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nscougv", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nscougv)
      call check_err(iret)
      if ((IM == 10) .or. (IMHS == 10)) then
         iret = nf90_inq_varid(hs%ncid, "nscougc", tempid)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, tempid, nscougc)
         call check_err(iret)
      end if
      iret = nf90_inq_varid(hs%ncid, "nscougw", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nscougw)
      call check_err(iret)

!     now close the netcdf file
      iret = nf90_close(hs%ncid)
      call check_err(iret)

      ! jgf52.17: If we started with this hotstart file, the
      ! metadata still reflects the fort.15 from the previous run,
      ! instead of this one, so we need to update it.
      call updateMetaData(hs%ncid, hs%myFile)

!-----------------------------------------------------------------------
   end subroutine writeNetCDFHotstart
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!        W R I T E   N E T C D F   H O T S T A R T   H A R M O N I C
!-----------------------------------------------------------------------
!     jgf49.44.03 Writes harmonic analysis data to hotstart file.
!-----------------------------------------------------------------------
   subroutine writeNetCDFHotstartHarmonic(lun, &
                                          GLOELVDescript, STAELVDescript, &
                                          GLOULVDescript, GLOVLVDescript, &
                                          STAULVDescript, STAVLVDescript)
      use SIZES, only: MNHARF, MNPROC
      use GLOBAL, only: OutputDataDescript_t, NSTAE_G, NSTAV_G
      use HARM, only: nz, nf, mm, nhase, nhasv, nhage, nhagv, icall, &
                      nfreq, timeud, itud, namefr, hafreq, haff, &
                      haface, ha, icha
      implicit none

      integer, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: GLOELVDescript
      type(OutputDataDescript_t), intent(in) :: STAELVDescript
      type(OutputDataDescript_t), intent(in) :: GLOULVDescript
      type(OutputDataDescript_t), intent(in) :: GLOVLVDescript
      type(OutputDataDescript_t), intent(in) :: STAULVDescript
      type(OutputDataDescript_t), intent(in) :: STAVLVDescript

      integer :: i
      integer :: kount(3), start(3) ! for nodally based data
      integer :: hakount(2), hastart(2) ! for lhs
      integer :: iret ! success or failure of the netcdf call
      integer :: tempid

      LOG_SCOPE_TRACED("writeNetCDFHotstartHarmonic", NETCDFIO_TRACING)

!     Point to the hotstart file we want to work on
      if (lun == 67) then
         hs => hs67
      else
         hs => hs68
      end if

      iret = nf90_open(hs%myFile%FILENAME, NF90_WRITE, hs%ncid)
      call check_err(iret)

!     Inquire variables (time)
      call getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)

      if (NHASE /= 0) then
!        elevation station dimension
         iret = nf90_inq_dimid(hs%ncid, "elevstation", &
                               hs%staelv%num_sta_dim_id)
         call check_err(iret)
         iret = nf90_inquire_dimension(hs%ncid, hs%staelv%num_sta_dim_id, &
                                       len=hs%staelv%num_stations)
         call check_err(iret)
      end if
      if (NHASV /= 0) then
!        velocity station dimension
         iret = nf90_inq_dimid(hs%ncid, "velstation", &
                               hs%stavellv%num_sta_dim_id)
         call check_err(iret)
         iret = nf90_inquire_dimension(hs%ncid, hs%stavellv%num_sta_dim_id, &
                                       len=hs%stavellv%num_stations)
         call check_err(iret)
      end if

!     Don't increment the record counter so that we can store data at the
!     next location in the netcdf file -- the hotstart file is only intended
!     to have a single snapshot of data in it.
      hs%myFile%record_counter = 1

      kount(1) = MNHARF*2 ! for load vector data
      kount(2) = hs%myMesh%num_nodes ! for nodal data
      kount(3) = hs%myTime%timenc_len
      start(1) = 1
      start(2) = 1
      start(3) = hs%myFile%record_counter

!     Get the NetCDF IDs of the relevant variables from the file
      if (NHAGE /= 0) then
         iret = nf90_inq_varid(hs%ncid, "gloelv", hs%gloelv%nodal_data_id)
         call check_err(iret)
      end if
      if (NHASE /= 0) then
         iret = nf90_inq_varid(hs%ncid, "staelv", hs%staelv%station_data_id)
         call check_err(iret)
      end if
      if (NHAGV /= 0) then
         iret = nf90_inq_varid(hs%ncid, "gloulv", hs%glovellv%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(hs%ncid, "glovlv", hs%glovellv%v_nodal_data_id)
         call check_err(iret)
      end if
      if (NHASV /= 0) then
         iret = nf90_inq_varid(hs%ncid, "staulv", &
                               hs%stavellv%u_station_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(hs%ncid, "stavlv", &
                               hs%stavellv%v_station_data_id)
         call check_err(iret)
      end if

!     Write the nodal data to the netcdf file
      if (MNPROC == 1) then
         if (NHAGE /= 0) then
            kount(2) = hs%myMesh%num_nodes ! for nodal data
            iret = nf90_put_var(hs%ncid, hs%gloelv%nodal_data_id, &
                                GLOELVDescript%array2D, start, kount)
            call check_err(iret)
         end if
         if (NHASE /= 0) then
            kount(2) = hs%staelv%num_stations ! for elevation stations
            iret = nf90_put_var(hs%ncid, hs%staelv%station_data_id, &
                                STAELVDescript%array2D, start, kount)
            call check_err(iret)
         end if
         if (NHAGV /= 0) then
            kount(2) = hs%myMesh%num_nodes ! for nodal data
            iret = nf90_put_var(hs%ncid, hs%glovellv%u_nodal_data_id, &
                                GLOULVDescript%array2D, start, kount)
            call check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%glovellv%v_nodal_data_id, &
                                GLOVLVDescript%array2D, start, kount)
            call check_err(iret)
         end if
         if (NHASV /= 0) then
            kount(2) = hs%stavellv%num_stations ! for velocity stations
            iret = nf90_put_var(hs%ncid, hs%stavellv%u_station_data_id, &
                                STAULVDescript%array2D, start, kount)
            call check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%stavellv%v_station_data_id, &
                                STAVLVDescript%array2D, start, kount)
            call check_err(iret)
         end if
      else
         if (NHAGE /= 0) then
            kount(2) = hs%myMesh%num_nodes ! for nodal data
            iret = nf90_put_var(hs%ncid, hs%gloelv%nodal_data_id, &
                                GLOELVDescript%array2D_g, start, kount)
            call check_err(iret)
         end if
         if (NHASE /= 0) then
            kount(2) = hs%staelv%num_stations ! for elevation stations
            iret = nf90_put_var(hs%ncid, hs%staelv%station_data_id, &
                                STAELVDescript%array2D_g, start, kount)
            call check_err(iret)
         end if
         if (NHAGV /= 0) then
            kount(2) = hs%myMesh%num_nodes ! for nodal data
            iret = nf90_put_var(hs%ncid, hs%glovellv%u_nodal_data_id, &
                                GLOULVDescript%array2D_g, start, kount)
            call check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%glovellv%v_nodal_data_id, &
                                GLOVLVDescript%array2D_g, start, kount)
            call check_err(iret)
         end if
         if (NHASV /= 0) then
            kount(2) = hs%stavellv%num_stations ! for velocity stations
            iret = nf90_put_var(hs%ncid, hs%stavellv%u_station_data_id, &
                                STAULVDescript%array2D_g, start, kount)
            call check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%stavellv%v_station_data_id, &
                                STAVLVDescript%array2D_g, start, kount)
            call check_err(iret)
         end if
      end if

!     Get each variable ID for the model parameters in the netcdf file
!     and immediately write the parameter value to that variable ID before
!     going on to the next one.
      iret = nf90_inq_varid(hs%ncid, "icha", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, icha)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nz", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nz)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nf", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nf)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "mm", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, mm)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nstae", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nstae_g)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nstav", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nstav_g)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nhase", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nhase)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nhasv", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nhasv)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nhage", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nhage)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nhagv", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nhagv)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "icall", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, icall)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nfreq", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, nfreq)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "timeud", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, timeud)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "itud", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, itud)
      call check_err(iret)

!     left hand side
      hakount(1) = 2*MNHARF
      hakount(2) = 2*MNHARF
      hastart(1) = 1
      hastart(2) = 1
      iret = nf90_inq_varid(hs%ncid, "ha", hs%ha_id)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, hs%ha_id, ha, hastart, hakount)
      call check_err(iret)
!     frequency names
      iret = nf90_inq_varid(hs%ncid, "namefr", hs%namefr_id)
      call check_err(iret)
      do i = 1, mnharf
         start(1) = 1
         start(2) = i
         kount(1) = len(namefr(i))
         kount(2) = 1
         iret = nf90_put_var(hs%ncid, hs%namefr_id, namefr(i), &
                             start, kount)
         call check_err(iret)
      end do
!     harmonic constituents
      start(1) = 1
      start(2) = 1
      kount(1) = MNHARF ! for constituents
      kount(2) = 1
      iret = nf90_inq_varid(hs%ncid, "hafreq", hs%hafreq_id)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, hs%hafreq_id, hafreq, start, kount)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "haff", hs%haff_id)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, hs%haff_id, haff, start, kount)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "haface", hs%haface_id)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, hs%haface_id, haface, start, kount)
      call check_err(iret)
!
!     now close the netcdf file
      iret = nf90_close(hs%ncid)
      call check_err(iret)
!
!-----------------------------------------------------------------------
   end subroutine writeNetCDFHotstartHarmonic
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!        W R I T E   N E T C D F   H O T S T A R T   H A R M O N I C
!                                     M E A N S   V A R I A N C E S
!-----------------------------------------------------------------------
!     jgf49.43.14 Writes harmonic analysis means and variance data to
!     hotstart file.
!-----------------------------------------------------------------------
   subroutine writeNetCDFHotstartHarmonicMeansVariances(lun, &
                                                        ELAVDescript, ELVADescript, &
                                                        XVELAVDescript, YVELAVDescript, &
                                                        XVELVADescript, YVELVADescript)
      use SIZES, only: MNPROC
      use GLOBAL, only: OutputDataDescript_t
      use HARM, only: nhage, nhagv, ntsteps
      implicit none

      integer, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: ELAVDescript
      type(OutputDataDescript_t), intent(in) :: ELVADescript
      type(OutputDataDescript_t), intent(in) :: XVELAVDescript
      type(OutputDataDescript_t), intent(in) :: YVELAVDescript
      type(OutputDataDescript_t), intent(in) :: XVELVADescript
      type(OutputDataDescript_t), intent(in) :: YVELVADescript

      integer :: kount(2), start(2) ! for nodally based data
      integer :: iret ! success or failure of the netcdf call
      integer :: tempid

      LOG_SCOPE_TRACED("writeNetCDFHotstartHarmonicMeansVariances", NETCDFIO_TRACING)

!     Point to the hotstart file we want to work on
      if (lun == 67) then
         hs => hs67
      else
         hs => hs68
      end if

      iret = nf90_open(hs%myFile%FILENAME, NF90_WRITE, hs%ncid)
      call check_err(iret)

!     Inquire variables
!     time dimension
      call getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)

!     Don't increment the record counter so that we can store data at the
!     next location in the netcdf file -- the hotstart file is only intended
!     to have a single snapshot of data in it.
      hs%myFile%record_counter = 1
!
      kount(1) = hs%myMesh%num_nodes ! for nodal data
      kount(2) = hs%myTime%timenc_len
      start(1) = 1
      start(2) = hs%myFile%record_counter

!     Get the NetCDF IDs of the relevant variables from the file
      if (NHAGE /= 0) then
         iret = nf90_inq_varid(hs%ncid, "elav", hs%elav%nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(hs%ncid, "elva", hs%elva%nodal_data_id)
         call check_err(iret)
      end if
      if (NHAGV /= 0) then
         iret = nf90_inq_varid(hs%ncid, "xvelav", hs%xvelav%nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(hs%ncid, "yvelav", hs%yvelav%nodal_data_id)
         call check_err(iret)

         iret = nf90_inq_varid(hs%ncid, "xvelva", hs%xvelva%nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(hs%ncid, "yvelva", hs%yvelva%nodal_data_id)
         call check_err(iret)
      end if

!     Write the nodal data to the netcdf file
      if (MNPROC == 1) then
         if (NHAGE /= 0) then
            iret = nf90_put_var(hs%ncid, hs%elav%nodal_data_id, &
                                ELAVDescript%array, start, kount)
            call check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%elva%nodal_data_id, &
                                ELVADescript%array, start, kount)
            call check_err(iret)
         end if
         if (NHAGV /= 0) then
            iret = nf90_put_var(hs%ncid, hs%xvelav%nodal_data_id, &
                                XVELAVDescript%array, start, kount)
            call check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%yvelav%nodal_data_id, &
                                YVELAVDescript%array, start, kount)
            call check_err(iret)

            iret = nf90_put_var(hs%ncid, hs%xvelva%nodal_data_id, &
                                XVELVADescript%array, start, kount)
            call check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%yvelva%nodal_data_id, &
                                YVELVADescript%array, start, kount)
            call check_err(iret)
         end if
      else
         if (NHAGE /= 0) then
            iret = nf90_put_var(hs%ncid, hs%elav%nodal_data_id, &
                                ELAVDescript%array_g, start, kount)
            call check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%elva%nodal_data_id, &
                                ELVADescript%array_g, start, kount)
            call check_err(iret)
         end if
         if (NHAGV /= 0) then
            iret = nf90_put_var(hs%ncid, hs%xvelav%nodal_data_id, &
                                XVELAVDescript%array_g, start, kount)
            call check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%yvelav%nodal_data_id, &
                                YVELAVDescript%array_g, start, kount)
            call check_err(iret)

            iret = nf90_put_var(hs%ncid, hs%xvelva%nodal_data_id, &
                                XVELVADescript%array_g, start, kount)
            call check_err(iret)
            iret = nf90_put_var(hs%ncid, hs%yvelva%nodal_data_id, &
                                YVELVADescript%array_g, start, kount)
            call check_err(iret)
         end if
      end if

!     Get each variable ID for the model parameters in the netcdf file
!     and immediately write the parameter value to that variable ID before
!     going on to the next one.
      iret = nf90_inq_varid(hs%ncid, "ntsteps", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, ntsteps)
      call check_err(iret)

!     now close the netcdf file
      iret = nf90_close(hs%ncid)
      call check_err(iret)

!-----------------------------------------------------------------------
   end subroutine writeNetCDFHotstartHarmonicMeansVariances
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!        W R I T E   N E T C D F   H O T S T A R T   3 D
!-----------------------------------------------------------------------
!     jgf49.49.02 Writes 3D data to hotstart file. Only includes arrays
!     of rank 1 (dimensioned by number of nodes, NP). The 3D variables
!     (i.e., of rank 2) will be written in a subsequent subroutine call.
!-----------------------------------------------------------------------
   subroutine writeNetCDFHotstart3D(lun, DUUDescript, &
                                    DUVDescript, DVVDescript, UUDescript, VVDescript, &
                                    BSXDescript, BSYDescript)
      use SIZES, only: MNPROC
      use GLOBAL, only: OutputDataDescript_t
      use GLOBAL_3DVS, only: n3dsd, i3dsdrec, n3dsv, i3dsvrec, &
                             n3dst, i3dstrec, n3dgd, i3dgdrec, n3dgv, i3dgvrec, &
                             n3dgt, i3dgtrec
      implicit none

      integer, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: DUUDescript
      type(OutputDataDescript_t), intent(in) :: DUVDescript
      type(OutputDataDescript_t), intent(in) :: DVVDescript
      type(OutputDataDescript_t), intent(in) :: UUDescript
      type(OutputDataDescript_t), intent(in) :: VVDescript
      type(OutputDataDescript_t), intent(in) :: BSXDescript
      type(OutputDataDescript_t), intent(in) :: BSYDescript

      integer :: kount(2), start(2) ! for nodally based data
      integer :: iret ! success or failure of the netcdf call
      integer :: tempid

      LOG_SCOPE_TRACED("writeNetCDFHotstart3D", NETCDFIO_TRACING)

!     Point to the hotstart file we want to work on
      if (lun == 67) then
         hs => hs67
      else
         hs => hs68
      end if

      iret = nf90_open(hs%myFile%FILENAME, NF90_WRITE, hs%ncid)
      call check_err(iret)

!     Inquire variables
!     time dimension
      call getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)

!     Don't increment the record counter so that we can store data at the
!     next location in the netcdf file -- the hotstart file is only intended
!     to have a single snapshot of data in it.
      hs%myFile%record_counter = 1

      kount(1) = hs%myMesh%num_nodes ! for nodal data
      kount(2) = hs%myTime%timenc_len
      start(1) = 1
      start(2) = hs%myFile%record_counter

!     Get the NetCDF IDs of the relevant variables from the file
      iret = nf90_inq_varid(hs%ncid, "duu", hs%duu%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "duv", hs%duv%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "dvv", hs%dvv%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "uu", hs%uu%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "vv", hs%vv%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "bsx", hs%bsx%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "bsy", hs%bsy%nodal_data_id)
      call check_err(iret)

!     Write the nodal data to the netcdf file
      if (MNPROC == 1) then
         iret = nf90_put_var(hs%ncid, hs%duu%nodal_data_id, &
                             DUUDescript%array, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%duv%nodal_data_id, &
                             DUVDescript%array, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%dvv%nodal_data_id, &
                             DVVDescript%array, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%uu%nodal_data_id, &
                             UUDescript%array, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%vv%nodal_data_id, &
                             VVDescript%array, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%bsx%nodal_data_id, &
                             BSXDescript%array, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%bsy%nodal_data_id, &
                             BSYDescript%array, start, kount)
         call check_err(iret)
      else
         iret = nf90_put_var(hs%ncid, hs%duu%nodal_data_id, &
                             DUUDescript%array_g, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%duv%nodal_data_id, &
                             DUVDescript%array_g, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%dvv%nodal_data_id, &
                             DVVDescript%array_g, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%uu%nodal_data_id, &
                             UUDescript%array_g, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%vv%nodal_data_id, &
                             VVDescript%array_g, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%bsx%nodal_data_id, &
                             BSXDescript%array_g, start, kount)
         call check_err(iret)
         iret = nf90_put_var(hs%ncid, hs%bsy%nodal_data_id, &
                             BSYDescript%array_g, start, kount)
         call check_err(iret)
      end if

!     Get each variable ID for the model parameters in the netcdf file
!     and immediately write the parameter value to that variable ID before
!     going on to the next one.
      iret = nf90_inq_varid(hs%ncid, "n3dsd", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, n3dsd)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "i3dsdrec", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, i3dsdrec)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "n3dsv", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, n3dsv)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "i3dsvrec", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, i3dsvrec)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "n3dst", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, n3dst)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "i3dstrec", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, i3dstrec)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "n3dgd", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, n3dgd)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "i3dgdrec", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, i3dgdrec)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "n3dgv", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, n3dgv)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "i3dgvrec", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, i3dgvrec)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "n3dgt", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, n3dgt)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "i3dgtrec", tempid)
      call check_err(iret)
      iret = nf90_put_var(hs%ncid, tempid, i3dgtrec)
      call check_err(iret)

!     now close the netcdf file
      iret = nf90_close(hs%ncid)
      call check_err(iret)

!-----------------------------------------------------------------------
   end subroutine writeNetCDFHotstart3D
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!        W R I T E   N E T C D F   H O T S T A R T   3 D   V A R
!-----------------------------------------------------------------------
!     jgf49.49.02 Writes 3D data to hotstart file for a single variable
!     of rank 2.
!-----------------------------------------------------------------------
   subroutine writeNetCDFHotstart3DVar(lun, descript)
      use SIZES, only: MNPROC
      use GLOBAL, only: OutputDataDescript_t
      implicit none

      integer, intent(in) :: lun
      type(OutputDataDescript_t), intent(in) :: descript

      integer :: kount(3), start(3) ! for nodally based data
      integer :: iret ! success or failure of the netcdf call

      LOG_SCOPE_TRACED("writeNetCDFHotstart3DVar", NETCDFIO_TRACING)

!     Point to the hotstart file we want to work on
      if (lun == 67) then
         hs => hs67
      else
         hs => hs68
      end if

      iret = nf90_open(hs%myFile%FILENAME, NF90_WRITE, hs%ncid)
      call check_err(iret)
!
!     Inquire variables
!     time dimension
      call getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)
!     vertical dimension
      iret = nf90_inq_dimid(hs%ncid, "num_v_nodes", &
                            hs%myMesh%num_v_nodes_dim_id)
      call check_err(iret)
      iret = nf90_inquire_dimension(hs%ncid, hs%myMesh%num_v_nodes_dim_id, &
                                    len=hs%myMesh%num_v_nodes)
      call check_err(iret)
!
!
!     Don't increment the record counter so that we can store data at the
!     next location in the netcdf file -- the hotstart file is only intended
!     to have a single snapshot of data in it.
      hs%myFile%record_counter = 1
!
      kount(1) = hs%myMesh%num_nodes
      kount(2) = hs%myMesh%num_v_nodes
      kount(3) = hs%myTime%timenc_len
      start(1) = 1
      start(2) = 1
      start(3) = hs%myFile%record_counter
!
!     Get the NetCDF ID of the relevant variable from the file
      select case (trim(descript%field_name))
      case ("SigmaT")
         iret = nf90_inq_varid(hs%ncid, "sigt", hs%density3D%u_nodal_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_put_var(hs%ncid, hs%density3D%u_nodal_data_id, &
                                descript%array2D, start, kount)
         else
            iret = nf90_put_var(hs%ncid, hs%density3D%u_nodal_data_id, &
                                descript%array2D_g, start, kount)
         end if
         call check_err(iret)
      case ("Salinity")
         iret = nf90_inq_varid(hs%ncid, "salinity", &
                               hs%density3D%v_nodal_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_put_var(hs%ncid, hs%density3D%v_nodal_data_id, &
                                descript%array2D, start, kount)
         else
            iret = nf90_put_var(hs%ncid, hs%density3D%v_nodal_data_id, &
                                descript%array2D_g, start, kount)
         end if
         call check_err(iret)
      case ("Temperature")
         iret = nf90_inq_varid(hs%ncid, "temperature", &
                               hs%density3D%nodal_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_put_var(hs%ncid, &
                                hs%density3D%w_nodal_data_id, &
                                descript%array2D, start, kount)
         else
            iret = nf90_put_var(hs%ncid, &
                                hs%density3D%w_nodal_data_id, &
                                descript%array2D_g, start, kount)
         end if
         call check_err(iret)
      case ("u-vel3D")
         iret = nf90_inq_varid(hs%ncid, "u-vel3D", &
                               hs%velocity3D%u_nodal_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_put_var(hs%ncid, hs%velocity3D%u_nodal_data_id, &
                                descript%array2D, start, kount)
         else
            iret = nf90_put_var(hs%ncid, hs%velocity3D%u_nodal_data_id, &
                                descript%array2D_g, start, kount)
         end if
         call check_err(iret)
      case ("v-vel3D")
         iret = nf90_inq_varid(hs%ncid, "v-vel3D", &
                               hs%velocity3D%v_nodal_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_put_var(hs%ncid, &
                                hs%velocity3D%v_nodal_data_id, descript%array2D, &
                                start, kount)
         else
            iret = nf90_put_var(hs%ncid, &
                                hs%velocity3D%v_nodal_data_id, descript%array2D_g, &
                                start, kount)
         end if
         call check_err(iret)
      case ("w-vel3D")
         iret = nf90_inq_varid(hs%ncid, "w-vel3D", &
                               hs%velocity3D%w_nodal_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_put_var(hs%ncid, hs%velocity3D%w_nodal_data_id, &
                                descript%array2D, start, kount)
         else
            iret = nf90_put_var(hs%ncid, hs%velocity3D%w_nodal_data_id, &
                                descript%array2D_g, start, kount)
         end if
         call check_err(iret)
      case ("q20")
         iret = nf90_inq_varid(hs%ncid, "q20", &
                               hs%turbulence3D%u_nodal_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_put_var(hs%ncid, &
                                hs%turbulence3D%u_nodal_data_id, &
                                descript%array2D, start, kount)
         else
            iret = nf90_put_var(hs%ncid, &
                                hs%turbulence3D%u_nodal_data_id, &
                                descript%array2D_g, start, kount)
         end if
         call check_err(iret)
      case ("l")
         iret = nf90_inq_varid(hs%ncid, "l", hs%turbulence3D%v_nodal_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_put_var(hs%ncid, &
                                hs%turbulence3D%v_nodal_data_id, &
                                descript%array2D, start, kount)
         else
            iret = nf90_put_var(hs%ncid, &
                                hs%turbulence3D%v_nodal_data_id, &
                                descript%array2D_g, start, kount)
         end if
         call check_err(iret)
      case DEFAULT

      end select
!
!     now close the netcdf file
      iret = nf90_close(hs%ncid)
      call check_err(iret)
!

!-----------------------------------------------------------------------
   end subroutine writeNetCDFHotstart3DVar
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E    R E A D   N E T C D F   H O T S T A R T
!-----------------------------------------------------------------------
!     jgf49.17.02 Reads data from the hotstart file.
!-----------------------------------------------------------------------
   subroutine readNetCDFHotstart(lun, timeLoc)
      use SIZES, only: globaldir, mnproc
      use GLOBAL, only: OutputDataDescript_t, imhs, iths, nscoue, &
                        nscouv, nscouc, nscoum, nscouge, nscougv, nscougc, &
                        nscougw, ETA1, ETA2, EtaDisc, &
                        UU2, VV2, NNODECODE, NOFF, &
                        IM, NP_G, NE_G, &
                        htot1 => H1, htot2 => H2, &
                        nrs, rsnx1, rsnx2, rsny1, rsny2
#ifdef CSWAN
      use GLOBAL, only: swan_rsnx1, swan_rsnx2, swan_rsny1, swan_rsny2
#endif
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none

      integer, intent(in) :: lun
      real(8), intent(out) :: timeLoc

      integer :: counti(1), starti(1)
      integer :: kount(2), start(2)
      integer :: iret ! success or failure of the netcdf call
      integer :: tempid

      character(len=10) :: fext
      character(len=1024) :: scratchMessage

      LOG_SCOPE_TRACED("readNetCDFHotstart", NETCDFIO_TRACING)
!
!     Point to the hotstart file we want to work on
      if (lun == 67) then
         hs => hs67
      else
         hs => hs68
      end if
!
      hs%myMesh => adcircMesh
      hs%myFile%fileFound = .false.
!     ! create file name
      write (fext, '(i0)') lun
      hs%myFile%filename = trim(globaldir)//'/fort.'//trim(fext)//'.nc'
      call logMessage(INFO, 'Opening hot start file "' &
                      //trim(hs%myFile%filename)//'" for reading.')
!
!     Open fulldomain file
      inquire (FILE=hs%myFile%FILENAME, EXIST=hs%myFile%fileFound)
      if (hs%myFile%fileFound .eqv. .false.) then
         write (scratchMessage, '(a,a,a)') 'The file ', &
            trim(adjustl(hs%myFile%FILENAME)), &
            ' was not found; ADCIRC terminating.'
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message=scratchMessage)
      else
         iret = nf90_open(hs%myFile%FILENAME, NF90_NOWRITE, hs%ncid)
         call check_err(iret)
      end if
!
!     Inquire variables
      call getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)
!
!     nodes and elements
      hs%myMesh => adcircMesh
      if (hs%myMesh%initialized .eqv. .false.) then
         hs%myMesh%num_nodes = NP_G
         hs%myMesh%num_elems = NE_G
         hs%myMesh%nface_len = 3
      end if
      if (hs%myTime%initialized .eqv. .false.) then
         allocate (hs%myTime%timenc(hs%myTime%timenc_len))
         hs%myTime%initialized = .true.
      end if
!
      hs%myFile%record_counter = 1
      kount(1) = hs%myMesh%num_nodes
      kount(2) = hs%myTime%timenc_len
      start(1) = 1
      start(2) = hs%myFile%record_counter
!
!     Get time
      iret = nf90_inq_varid(hs%ncid, "time", hs%myTime%timenc_id)
      call check_err(iret)
      starti(1) = hs%myFile%record_counter
      counti(1) = hs%myTime%timenc_len
      iret = nf90_get_var(hs%ncid, hs%myTime%timenc_id, &
                          hs%myTime%timenc, starti, counti)
      call check_err(iret)
!     ! set timeLoc in hstart.F to the current time in the hotstart file
      timeLoc = hs%myTime%timenc(hs%myFile%record_counter)
!     !
!     ! get array variable ids
      iret = nf90_inq_varid(hs%ncid, "zeta1", hs%zeta1%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "zeta2", hs%zeta2%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "zetad", hs%zetad%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "u-vel", hs%vel%u_nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "v-vel", hs%vel%v_nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nodecode", hs%nodecodenc%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "noff", hs%noffnc%nodal_data_id)
      call check_err(iret)
      if (nrs /= 0) then
         iret = nf90_inq_varid(hs%ncid, "rsx1", hs%rs1%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(hs%ncid, "rsy1", hs%rs1%v_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(hs%ncid, "rsx2", hs%rs2%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(hs%ncid, "rsy2", hs%rs2%v_nodal_data_id)
         call check_err(iret)
#ifdef CSWAN
         if (nrs == 3) then
            iret = nf90_inq_varid(hs%ncid, "swan_rsx1", hs%swan_rs1%u_nodal_data_id)
            call check_err(iret)
            iret = nf90_inq_varid(hs%ncid, "swan_rsy1", hs%swan_rs1%v_nodal_data_id)
            call check_err(iret)
            iret = nf90_inq_varid(hs%ncid, "swan_rsx2", hs%swan_rs2%u_nodal_data_id)
            call check_err(iret)
            iret = nf90_inq_varid(hs%ncid, "swan_rsy2", hs%swan_rs2%v_nodal_data_id)
            call check_err(iret)
         end if
#endif
      end if
      iret = nf90_inq_varid(hs%ncid, "h1", hs%htot1%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "h2", hs%htot2%nodal_data_id)
      call check_err(iret)

      ! serial
      if (MNPROC == 1) then
!        Elev1
         iret = nf90_get_var(hs%ncid, hs%zeta1%nodal_data_id, &
                             eta1, start, kount)
         call check_err(iret)
!        Elev2
         iret = nf90_get_var(hs%ncid, hs%zeta2%nodal_data_id, &
                             eta2, start, kount)
         call check_err(iret)
!        EtaDisc
         iret = nf90_get_var(hs%ncid, hs%zetad%nodal_data_id, &
                             EtaDisc, start, kount)
         call check_err(iret)
!        Vel
         iret = nf90_get_var(hs%ncid, hs%vel%u_nodal_data_id, &
                             uu2, start, kount)
         call check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%vel%v_nodal_data_id, &
                             vv2, start, kount)
         call check_err(iret)
!        NodeCode
         iret = nf90_get_var(hs%ncid, hs%nodecodenc%nodal_data_id, &
                             nnodecode, start, kount)
         call check_err(iret)

         if (nrs /= 0) then
            iret = nf90_get_var(hs%ncid, hs%rs1%u_nodal_data_id, &
                                rsnx1, start, kount)
            call check_err(iret)
            iret = nf90_get_var(hs%ncid, hs%rs1%v_nodal_data_id, &
                                rsny1, start, kount)
            call check_err(iret)
            iret = nf90_get_var(hs%ncid, hs%rs2%u_nodal_data_id, &
                                rsnx2, start, kount)
            call check_err(iret)
            iret = nf90_get_var(hs%ncid, hs%rs2%v_nodal_data_id, &
                                rsny2, start, kount)
            call check_err(iret)
#ifdef CSWAN
            if (nrs == 3) then
               iret = nf90_get_var(hs%ncid, hs%rs1%u_nodal_data_id, &
                                   swan_rsnx1, start, kount)
               call check_err(iret)
               iret = nf90_get_var(hs%ncid, hs%rs1%v_nodal_data_id, &
                                   swan_rsny1, start, kount)
               call check_err(iret)
               iret = nf90_get_var(hs%ncid, hs%rs2%u_nodal_data_id, &
                                   swan_rsnx2, start, kount)
               call check_err(iret)
               iret = nf90_get_var(hs%ncid, hs%rs2%v_nodal_data_id, &
                                   swan_rsny2, start, kount)
               call check_err(iret)
            end if
#endif
         end if

         ! hot start file contain the varible h1 and h2
         ! H 2
         call check_err(nf90_get_var(hs%ncid, hs%htot2%nodal_data_id, &
                                     htot2, start, kount))

         ! H 1
         call check_err(nf90_get_var(hs%ncid, hs%htot1%nodal_data_id, &
                                     htot1, start, kount))

!        NOFF
         start(1) = 1
         kount(1) = hs%myMesh%num_elems
         iret = nf90_get_var(hs%ncid, hs%noffnc%nodal_data_id, &
                             noff, start, kount)
         call check_err(iret)

      else ! parallel

         ! Create the fulldomain node and element index lists for this subdomain
         call createFullDomainIndexLists()

         ! get fulldomain data and map the data to this subdomain
         fullDomainIndexList => fullDomainNodeList
         call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                       hs%zeta1%nodal_data_id, subdomain_reals=eta1)
         call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                       hs%zeta2%nodal_data_id, subdomain_reals=eta2)
         call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                       hs%zetad%nodal_data_id, subdomain_reals=EtaDisc)
         call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                       hs%vel%u_nodal_data_id, subdomain_reals=uu2)
         call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                       hs%vel%v_nodal_data_id, subdomain_reals=vv2)
         call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                       hs%nodecodenc%nodal_data_id, subdomain_ints=nnodecode)

         if (nrs /= 0) then
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                          hs%rs1%u_nodal_data_id, subdomain_reals=rsnx1)
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                          hs%rs1%v_nodal_data_id, subdomain_reals=rsny1)
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                          hs%rs2%u_nodal_data_id, subdomain_reals=rsnx2)
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                          hs%rs2%v_nodal_data_id, subdomain_reals=rsny2)
#ifdef CSWAN
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                          hs%swan_rs1%u_nodal_data_id, subdomain_reals=swan_rsnx1)
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                          hs%swan_rs1%v_nodal_data_id, subdomain_reals=swan_rsny1)
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                          hs%swan_rs2%u_nodal_data_id, subdomain_reals=swan_rsnx2)
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                          hs%swan_rs2%v_nodal_data_id, subdomain_reals=swan_rsny2)
#endif
         end if

         ! H 2
         call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                       hs%htot2%nodal_data_id, subdomain_reals=htot2)

         ! H 1
         call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                       hs%htot1%nodal_data_id, subdomain_reals=htot1)

         fullDomainIndexList => fullDomainElementList
         call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_elems, &
                                       hs%noffnc%nodal_data_id, subdomain_ints=noff)

      end if

!     Read in model parameters to ADCIRC variables
      iret = nf90_inq_varid(hs%ncid, "imhs", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, imhs)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "iths", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, iths)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nscoue", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, nscoue)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nscouv", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, nscouv)
      call check_err(iret)
      if ((IM == 10) .or. (IMHS == 10)) then
         iret = nf90_inq_varid(hs%ncid, "nscouc", tempid)
         call check_err(iret)
         iret = nf90_get_var(hs%ncid, tempid, nscouc)
         call check_err(iret)
      end if
      iret = nf90_inq_varid(hs%ncid, "nscoum", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, nscoum)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nscouge", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, nscouge)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nscougv", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, nscougv)
      call check_err(iret)
      if ((IM == 10) .or. (IMHS == 10)) then
         iret = nf90_inq_varid(hs%ncid, "nscougc", tempid)
         call check_err(iret)
         iret = nf90_get_var(hs%ncid, tempid, nscougc)
         call check_err(iret)
      end if
      iret = nf90_inq_varid(hs%ncid, "nscougw", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, nscougw)
      call check_err(iret)

      deallocate (hs%myTime%timenc)
      hs%myTime%initialized = .false. !..zc50.93 - reinitialize if we write a hot start later

!     now close the netcdf file
      iret = nf90_close(hs%ncid)
      call check_err(iret)

!-----------------------------------------------------------------------
   end subroutine readNetCDFHotstart
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!                    S U B R O U T I N E
!                 R E A D   A N D   M A P   T O
!         S U B D O M A I N   M I N   M A X   N E T C D F
!-----------------------------------------------------------------------
!     jgf52.08.03: Read in the min/max values from a previous run
!     at hotstart. In order to read min/max data from netCDF, we also
!     have to be writing that min/max data in netCDF. TODO: relax this
!     limitation.
!-----------------------------------------------------------------------
   subroutine readAndMapToSubdomainMaxMinNetCDF(descript, timeloc)
      use sizes, only: mnproc, myproc
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      use global, only: OutputDataDescript_t, nfover
      use mod_logging, only: allMessage, logMessage
#ifdef CMPI
      use mesh, only: np
      use global, only: np_g
      use messenger, only: mapToSubdomainRealMPI, mapToSubdomainIntMPI, msg_lbcast
#endif
      implicit none
      type(OutputDataDescript_t), intent(inout) :: descript
      real(8), intent(in) :: timeloc ! adcirc time in seconds since cold start

      type(nodalData), pointer :: dat
      integer :: counti(1), starti(1)
      integer :: kount(2), start(2)
      integer :: iret ! success or failure of the netcdf call
      logical :: timestampsFound !.true. if file contains time of occurrence dataset
      logical :: no_early_return
      character(len=1024) :: vn ! min/max variable name in netcdf file
      character(len=1025) :: tvn ! time of min/max occurrence variable name in netcdf
      character(len=1024) :: scratchMessage

#ifdef CMPI
      real(8), allocatable :: tmp0(:), tmp1(:)
      logical :: ldmy(1)
#endif

      LOG_SCOPE_TRACED("readAndMapToSubdomainMaxMinNetCDF", NETCDFIO_TRACING)

      timestampsFound = .false.

      select case (descript%lun)
      case (311) ! maxele
         dat => EtaMax
         vn = 'zeta_max'
         tvn = 'time_of_zeta_max'
      case (312) ! maxvel
         dat => UMax
         vn = 'vel_max'
         tvn = 'time_of_vel_max'
      case (313) ! minpr
         dat => PrMin
         vn = 'pressure_min'
         tvn = 'time_of_pressure_min'
      case (314) ! maxwvel
         dat => WVMax
         vn = 'wind_max'
         tvn = 'time_of_wind_max'
      case (315) ! maxrs
         dat => RSMax
         vn = 'radstress_max'
         tvn = 'time_of_radstress_max'
      case (400) ! inundationtime
         dat => inTime
         vn = 'inun_time'
         tvn = 'onset_inun_time'
      case (401) ! maxinundepth
         dat => maxInDep
         vn = 'inun_max'
         tvn = 'time_of_inun_max'
      case (404) ! everdried
         dat => evrDry
         vn = 'everdried'
         tvn = 'time_of_everdried'
      case default
         call allMessage(ERROR, 'Files of the type "' &
                         //trim(descript%file_name)//'.nc" cannot be read by ADCIRC.')
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG, "Return.")
#endif
         no_early_return = .false.
         call communicateMapToSubdomainMaxMinNetcdfStatus(no_early_return)
         return
      end select

      if (myproc == 0) then
         dat%myFile%filename = trim(descript%file_name)//'.nc'
         call logMessage(INFO, 'Opening file "' &
                         //trim(dat%myFile%filename)//'" for reading.')

!       Open fulldomain file
         inquire (file=dat%myfile%filename, exist=dat%myfile%filefound)
         if (dat%myFile%fileFound .eqv. .false.) then
            call allMessage(INFO, 'The file "'//trim(dat%myFile%FILENAME)// &
                            '" was not found.')
            no_early_return = .false.
            call communicateMapToSubdomainMaxMinNetcdfStatus(no_early_return)
            return ! RETURN RETURN
         else
            ! the file was found, let's open it
            iret = nf90_open(dat%myFile%FILENAME, NF90_NOWRITE, dat%ncid)
            call check_err(iret)
         end if

!        nodes and elements
         dat%myMesh => adcircMesh

!        Inquire variables
         call getDimensions(dat%ncid, dat%myTime, dat%myMesh, dat%myFile)

         ! jgf52.08.21: If the netcdf min/max file contains no time records,
         ! it was just created by adcprep, and there are no data to load.
         if (dat%myFile%record_counter == 0) then
            call allMessage(INFO, 'The file '//trim(dat%myFile%filename)// &
                            'contains no data, '// &
                            'so the min/max record will be started anew.')

            ! now close the netcdf file
            iret = nf90_close(dat%ncid)
            call check_err(iret)
            no_early_return = .false.
            call communicateMapToSubdomainMaxMinNetcdfStatus(no_early_return)
            return ! EARLY RETURN
         else
            write (scratchMessage, 3333) trim(descript%file_name)
3333        format('Values from ', (a), &
                   ' will be reflected from the solution prior to this hotstart.')
            call allMessage(INFO, scratchMessage)
         end if
         if (dat%myTime%initialized .eqv. .false.) then
            allocate (dat%myTime%timenc(dat%myTime%timenc_len))
            dat%myTime%initialized = .true.
         end if

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
         iret = nf90_inq_varid(dat%ncid, "time", dat%myTime%timenc_id)
         call check_err(iret)
         starti(1) = 1 ! min/max files have one element in the time array
         counti(1) = 1
         iret = nf90_get_var(dat%ncid, dat%myTime%timenc_id, &
                             dat%myTime%timenc, starti, counti)
         call check_err(iret)
         if (dat%myTime%timenc(1) > timeloc) then
            call allMessage(ERROR, 'Max/min file '//trim(dat%myFile%filename)// &
                            ' contains data written after the hotstart time. '// &
                            ' It may have been produced using different input parameters, '// &
                            'perhaps during a previously attempted hot start run. '// &
                            ' Its values will not be read.')
            if (nfover == 0) then
               call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                              message='Execution terminated due to invalid '// &
                              trim(dat%myFile%filename)//' file.')
            else
               call allMessage(INFO, 'The record for '// &
                               trim(dat%myFile%filename)// &
                               ' will be started anew and exection will continue.')

               ! now close the netcdf file
               iret = nf90_close(dat%ncid)
               call check_err(iret)

               no_early_return = .false.
               call communicateMapToSubdomainMaxMinNetcdfStatus(no_early_return)
               return ! EARLY RETURN
            end if
         end if

         dat%myFile%record_counter = 1
         kount(1) = dat%myMesh%num_nodes
         kount(2) = dat%myTime%timenc_len
         start(1) = 1
         start(2) = dat%myFile%record_counter

         ! get array variable ids
         !call logMessage(DEBUG,'netcdf var is '//trim(vn))
         iret = nf90_inq_varid(dat%ncid, trim(vn), dat%max_nodal_data_id)
         call check_err(iret)
         timestampsFound = .true.
         iret = nf90_inq_varid(dat%ncid, trim(tvn), dat%time_max_nodal_data_id)
         if (iret /= NF90_NOERR) then
            call logMessage(INFO, 'The file "'//trim(dat%myFile%filename)// &
                            '" does not contain time of occurrence for the min/max data. '// &
                            'As a result, time of occurrence will be initialized to -99999.')
            timestampsFound = .false.
         end if

#ifdef CMPI
         if (mnproc > 1) then
            no_early_return = .true.
            call communicateMapToSubdomainMaxMinNetcdfStatus(no_early_return)
            ldmy(1) = timestampsFound
            call msg_lbcast(ldmy, 1)
         end if
#endif
      else
#ifdef CMPI
         no_early_return = .false.
         call communicateMapToSubdomainMaxMinNetcdfStatus(no_early_return)
         if (.not. no_early_return) return

         call msg_lbcast(ldmy, 1)
         timestampsFound = ldmy(1)
#endif
      end if

      ! serial
      if (mnproc == 1) then

         if (descript%isInteger .eqv. .true.) then
            ! integer
            iret = nf90_get_var(dat%ncid, dat%max_nodal_data_id, &
                                descript%iarray, start, kount)
            call check_err(iret)
         else
            ! float
            iret = nf90_get_var(dat%ncid, dat%max_nodal_data_id, &
                                descript%array, start, kount)
            call check_err(iret)
         end if
         ! time of occurrence
         if (timestampsFound .eqv. .true.) then
            iret = nf90_get_var(dat%ncid, dat%time_max_nodal_data_id, &
                                descript%array2, start, kount)
            call check_err(iret)
         end if

      else ! parallel
#ifdef CMPI
         ! Create the fulldomain node and element index lists for this subdomain
         ! if it has not already happened for the reading of the netcdf
         ! hotstart file.
         if (fullDomainIndexListsInitialized .eqv. .false.) then
            call createFullDomainIndexLists()
         end if

         ! get fulldomain data and map the data to this subdomain
         fullDomainIndexList => fullDomainNodeList
         if (descript%isInteger .eqv. .true.) then
            ! integer
            if (myproc == 0) then
               iret = nf90_get_var(dat%ncid, dat%max_nodal_data_id, descript%iarray_g)
               call check_err(iret)
            end if
            call mapToSubdomainIntMPI(dat%myMesh%num_nodes, np, descript%iarray, descript%imap, descript%iarray_g)
         else
            ! float
            if (myproc == 0) then
               iret = nf90_get_var(dat%ncid, dat%max_nodal_data_id, descript%array_g)
               call check_err(iret)
            end if

            if (allocated(tmp0)) deallocate (tmp0)
            if (allocated(tmp1)) deallocate (tmp1)
            allocate (tmp0(1:np_g))
            allocate (tmp1(1:np))
            if (myproc == 0) then
               tmp0(:) = descript%array_g(:)
            end if
            call mapToSubdomainRealMPI(np_g, np, tmp1, descript%imap, tmp0)
            descript%array(:) = tmp1(:)
         end if
         ! time of occurrence
         if (timestampsFound) then
            if (myproc == 0) then
               iret = nf90_get_var(dat%ncid, dat%time_max_nodal_data_id, descript%array2_g)
               call check_err(iret)
            end if
            if (allocated(tmp0)) deallocate (tmp0)
            if (allocated(tmp1)) deallocate (tmp1)
            allocate (tmp0(1:np_g))
            allocate (tmp1(1:np))
            if (myproc == 0) then
               tmp0(:) = descript%array2_g(:)
            end if
            call mapToSubdomainRealMPI(np_g, np, tmp1, descript%imap, tmp0)
            descript%array2(:) = tmp1(:)
         end if
#endif
      end if

      if (myproc == 0) then
         ! now close the netcdf file
         iret = nf90_close(dat%ncid)
         call check_err(iret)
      end if

!-----------------------------------------------------------------------
   end subroutine readAndMapToSubdomainMaxMinNetCDF
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Alert the other processors if processor 0 tries to duck out of the
! read early
!-----------------------------------------------------------------------
   subroutine communicateMapToSubdomainMaxMinNetcdfStatus(communicateStatus)
!-----------------------------------------------------------------------
#ifdef CMPI
      use messenger, only: msg_lbcast
#endif
      implicit none
      logical, intent(inout) :: communicateStatus
      logical :: tmp(1)
      tmp(1) = communicateStatus
#ifdef CMPI
      call msg_lbcast(tmp, 1)
      communicateStatus = tmp(1)
#endif
!-----------------------------------------------------------------------
   end subroutine communicateMapToSubdomainMaxMinNetcdfStatus
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!        C R E A T E   F U L L D O M A I N   I N D E X   L I S T S
!-----------------------------------------------------------------------
!     jgf52.08.03 Create fulldomain node and element index lists for
!     this subdomain. This subroutine uses code originally found in
!     readNetCDFHotstart; these lists need to be created for reading
!     min/max files as well as hotstart files. If a netCDF hotstart
!     file was already read, these index lists will have already been
!     initialized, so we wouldn't have to do it again.
!-----------------------------------------------------------------------
   subroutine createFullDomainIndexLists()
      use mesh, only: ne, np
      use global, only: nodes_lg, imap_el_lg
      implicit none
      integer :: sd_element_number
      integer :: sd_node_number
      LOG_SCOPE_TRACED("createFullDomainIndexLists", NETCDFIO_TRACING)
      ! make a list of full domain nodes that correspond to this
      ! subdomain's nodes
      allocate (fullDomainNodeList(np))
      ! loop over subdomain indexes to form a list of corresponding
      ! fulldomain indexes
      forall (sd_node_number=1:np)
         ! get the corresponding fulldomain indexes
         fullDomainNodeList(sd_node_number) &
            = abs(nodes_lg(sd_node_number))
      end forall
      ! make a list of full domain elements that correspond to this
      ! subdomain's elements
      allocate (fullDomainElementList(ne))
      ! loop over subdomain indexes to form a list of corresponding
      ! fulldomain indexes
      forall (sd_element_number=1:ne)
         ! get the corresponding fulldomain indexes
         fullDomainElementList(sd_element_number) &
            = abs(imap_el_lg(sd_element_number))
      end forall
      fullDomainIndexListsInitialized = .true.

!-----------------------------------------------------------------------
   end subroutine createFullDomainIndexLists
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!        M A P   F U L L D O M A I N   T O   S U B D O M A I N
!-----------------------------------------------------------------------
!     jgf50.60.03 Maps full domain data to subdomain.
!-----------------------------------------------------------------------
   subroutine mapFullDomainToSubdomain(ncid, fd_array_size, &
                                       data_id, subdomain_reals, subdomain_ints)
      implicit none
      integer, intent(in) :: ncid ! file id to pull data from
      integer, intent(in) :: fd_array_size ! highest index in fulldomain array
      integer, intent(in) :: data_id ! netcdf variable id in file
      real(8), optional, intent(out) :: subdomain_reals(:) ! we need
      integer, optional, intent(out) :: subdomain_ints(:) ! we need

      real(8), allocatable :: work_reals(:) ! holds fulldomain data
      integer, allocatable :: work_ints(:) ! holds fulldomain data
      integer :: iret ! success or failure of the netcdf call

      LOG_SCOPE_TRACED("mapFullDomainToSubDomain", NETCDFIO_TRACING)

      ! grab array of full domain data from netcdf file and pull out the
      ! data that is needed for this subdomain
      if (present(subdomain_reals)) then
         allocate (work_reals(fd_array_size))
         iret = nf90_get_var(ncid, data_id, work_reals)
         call check_err(iret)
         subdomain_reals(:) = work_reals(fullDomainIndexList)
         deallocate (work_reals)
      else
         allocate (work_ints(fd_array_size))
         iret = nf90_get_var(ncid, data_id, work_ints)
         call check_err(iret)
         subdomain_ints(:) = work_ints(fullDomainIndexList)
         deallocate (work_ints)
      end if

!-----------------------------------------------------------------------
   end subroutine mapFulldomainToSubdomain
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!     M A P   F U L L D O M A I N   T O   S U B D O M A I N  M   B Y   NP
!
!-----------------------------------------------------------------------
!     jgf50.60.04 Maps mutlidimensional data from full domain data to
!     subdomain; specifically, this subroutine is useful for harmonic
!     analysis data.
!-----------------------------------------------------------------------
   subroutine mapFullDomainToSubdomainMByNP(ncid, m, n, &
                                            data_id, subdomain_reals)
      implicit none
      integer, intent(in) :: ncid ! file id to pull data from
      integer, intent(in) :: m ! non-nodal dimension
      integer, intent(in) :: n ! number of horizontal nodes
      integer, intent(in) :: data_id ! netcdf variable id in file
      real(8), optional, intent(out) :: subdomain_reals(:, :) ! we need

      real(8), allocatable :: work_reals(:, :) ! holds fulldomain data
      integer :: iret ! success or failure of the netcdf call
      integer :: i ! loop counter

      LOG_SCOPE_TRACED("mapFullDomainToSubDomainMbyNP", NETCDFIO_TRACING)
      ! grab array of full domain data from netcdf file and pull out the
      ! data that is needed for this subdomain
      allocate (work_reals(m, n))
      iret = nf90_get_var(ncid, data_id, work_reals)
      call check_err(iret)
      do i = 1, size(subdomain_reals(1, :))
         subdomain_reals(:, i) = work_reals(:, fullDomainIndexList(i))
      end do
      deallocate (work_reals)

!-----------------------------------------------------------------------
   end subroutine mapFulldomainToSubdomainMByNP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!     M A P   F U L L D O M A I N   T O   S U B D O M A I N  NP   B Y   M
!
!-----------------------------------------------------------------------
!     jgf50.60.04 Maps mutlidimensional data from full domain data to
!     subdomain; specifically, this subroutine is useful for 3D data.
!-----------------------------------------------------------------------
   subroutine mapFullDomainToSubdomainNPByM(ncid, m, n, &
                                            data_id, subdomain_reals)
      implicit none
      integer, intent(in) :: ncid ! file id to pull data from
      integer, intent(in) :: m ! non-nodal dimension
      integer, intent(in) :: n ! number of horizontal nodes
      integer, intent(in) :: data_id ! netcdf variable id in file
      real(8), optional, intent(out) :: subdomain_reals(:, :) ! we need

      real(8), allocatable :: work_reals(:, :) ! holds fulldomain data
      integer :: iret ! success or failure of the netcdf call

      LOG_SCOPE_TRACED("mapFullDomainToSubDomainNPbyM", NETCDFIO_TRACING)
      ! grab array of full domain data from netcdf file and pull out the
      ! data that is needed for this subdomain
      allocate (work_reals(n, m))
      iret = nf90_get_var(ncid, data_id, work_reals)
      call check_err(iret)
      subdomain_reals(:, :) = work_reals(fullDomainIndexList, :)
      deallocate (work_reals)

!-----------------------------------------------------------------------
   end subroutine mapFulldomainToSubdomainNPByM
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E   G E T   D I M E N S I O N S
!-----------------------------------------------------------------------
!     jgf50.60.04 Pulls the number of nodes, number of elements, and
!     number of time records from the netcdf file, as well as the
!     dimension ids of those dimensions.
!-----------------------------------------------------------------------
   subroutine getDimensions(ncid, time_struct, mesh_struct, &
                            file_struct)
      implicit none
      integer, intent(in) :: ncid ! file id
      type(timeData), intent(inout) :: time_struct
      type(fileData), intent(inout) :: file_struct
      type(meshStructure), intent(inout) :: mesh_struct

      integer :: iret ! success or failure of netcdf call

      LOG_SCOPE_TRACED("getDimensions", NETCDFIO_TRACING)

!     Inquire variables
      time_struct%timenc_dim_id = 0
      iret = nf90_inq_dimid(ncid, "time", time_struct%timenc_dim_id)
      call check_err(iret)
      iret = nf90_inquire_dimension(ncid, time_struct%timenc_dim_id, &
                                    len=file_struct%record_counter)
      call check_err(iret)

      iret = nf90_inq_dimid(ncid, "node", mesh_struct%num_nodes_dim_id)
      call check_err(iret)
      iret = nf90_inquire_dimension(ncid, mesh_struct%num_nodes_dim_id, &
                                    len=mesh_struct%num_nodes)
      call check_err(iret)

      iret = nf90_inq_dimid(ncid, "nfaces", mesh_struct%num_elems_dim_id)
      call check_err(iret)
      iret = nf90_inquire_dimension(ncid, mesh_struct%num_elems_dim_id, &
                                    len=mesh_struct%num_elems)
      call check_err(iret)

!-----------------------------------------------------------------------
   end subroutine getDimensions
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!        R E A D   N E T C D F   H O T S T A R T   H A R M O N I C
!-----------------------------------------------------------------------
!     jgf49.44.11 Reads harmonic analysis data from the hotstart file.
!-----------------------------------------------------------------------
   subroutine readNetCDFHotstartHarmonic(lun)
      use SIZES, only: mnproc
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      use GLOBAL, only: OutputDataDescript_t, &
                        IMAP_STAE_LG, NSTAE, IMAP_STAV_LG, NSTAV
      use HARM, only: GLOELV, STAELV, GLOULV, GLOVLV, STAULV, STAVLV, &
                      NHASE, NHASV, NHAGE, NHAGV, MNHARF, ICHA, INZ, &
                      INZ, INF, IMM, INSTAE, INSTAV, INHASE, INHASV, &
                      INHASE, INHAGE, INHAGV, ICALL, INFREQ, TIMEUD, &
                      ITUD, HA, INAMEFR, INP, IFF, IFACE, IFREQ

      implicit none

      integer, intent(in) :: lun

      integer :: kount(3), start(3)
      integer :: hakount(2), hastart(2)
      integer :: sd_station_number
      integer :: iret ! success or failure of the netcdf call
      integer :: tempid
      integer :: i
      character(1024) :: scratchMessage

      LOG_SCOPE_TRACED("readNetCDFHotstartHarmonic", NETCDFIO_TRACING)

!     Point to the hotstart file we want to work on
      if (lun == 67) then
         hs => hs67
      else
         hs => hs68
      end if

      hs%myMesh => adcircMesh
      hs%myFile%fileFound = .false.

!     Open fulldomain file
      inquire (FILE=hs%myFile%FILENAME, EXIST=hs%myFile%fileFound)
      if (hs%myFile%fileFound .eqv. .false.) then
         write (scratchMessage, '(a,a,a)') "The file ", &
            trim(hs%myFile%FILENAME), " was not found; ADCIRC terminating."
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message=scratchMessage)
      else
         iret = nf90_open(hs%myFile%FILENAME, NF90_NOWRITE, hs%ncid)
         call check_err(iret)
      end if

!     Inquire variables
      call getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)

!     elevation station dimension
      if (NHASE /= 0) then
         iret = nf90_inq_dimid(hs%ncid, "elevstation", &
                               hs%staelv%num_sta_dim_id)
         call check_err(iret)
         iret = nf90_inquire_dimension(hs%ncid, hs%staelv%num_sta_dim_id, &
                                       len=hs%staelv%num_stations)
         call check_err(iret)
         if (MNPROC > 1) then
            ! make a list of full domain elevation stations that
            ! correspond to this subdomain's stations
            allocate (fullDomainElevationStationList(hs%staelv%num_stations))
            ! loop over subdomain indexes to form a list of corresponding
            ! fulldomain indexes
            forall (sd_station_number=1:nstae)
               ! get the corresponding fulldomain indexes
               fullDomainElevationStationList(sd_station_number) &
                  = abs(imap_stae_lg(sd_station_number))
            end forall
         end if
      end if
!     velocity station dimension
      if (NHASV /= 0) then
         iret = nf90_inq_dimid(hs%ncid, "velstation", &
                               hs%stavellv%num_sta_dim_id)
         call check_err(iret)
         iret = nf90_inquire_dimension(hs%ncid, hs%stavellv%num_sta_dim_id, &
                                       len=hs%stavellv%num_stations)
         call check_err(iret)
         if (MNPROC > 1) then
            ! make a list of full domain velocity stations that
            ! correspond to this subdomain's velocity stations
            allocate (fullDomainVelocityStationList &
                      (hs%stavellv%num_stations))
            ! loop over subdomain indexes to form a list of corresponding
            ! fulldomain indexes
            forall (sd_station_number=1:nstav)
               ! get the corresponding fulldomain indexes
               fullDomainVelocityStationList(sd_station_number) &
                  = abs(imap_stav_lg(sd_station_number))
            end forall
         end if
      end if

!     Point to the hotstart file we want to work on.
      if (lun == 67) then
         hs => hs67
      else
         hs => hs68
      end if

      hs%myFile%record_counter = 1

      kount(1) = MNHARF*2 ! for load vector data
      kount(2) = hs%myMesh%num_nodes ! for nodal data
      kount(3) = hs%myTime%timenc_len
      start(1) = 1
      start(2) = 1
      start(3) = hs%myFile%record_counter

!     Read in fulldomain load vector data
!     GLOELV - full domain elevation
      fullDomainIndexList => fullDomainNodeList
      if (NHAGE /= 0) then
         iret = nf90_inq_varid(hs%ncid, "gloelv", hs%gloelv%nodal_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_get_var(hs%ncid, hs%gloelv%nodal_data_id, &
                                GLOELV, start, kount)
            call check_err(iret)
         else
            call mapFullDomainToSubdomainMByNP(hs%ncid, &
                                               2*MNHARF, hs%myMesh%num_nodes, &
                                               hs%gloelv%nodal_data_id, subdomain_reals=GLOELV)
         end if
      end if
!     GLOULV - fulldomain u velocity
!     GLOVLV - fulldomain v velocity
      if (NHAGV /= 0) then
         iret = nf90_inq_varid(hs%ncid, "gloulv", hs%glovellv%u_nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(hs%ncid, "glovlv", hs%glovellv%v_nodal_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_get_var(hs%ncid, hs%glovellv%u_nodal_data_id, &
                                GLOULV, start, kount)
            call check_err(iret)
            iret = nf90_get_var(hs%ncid, hs%glovellv%v_nodal_data_id, &
                                GLOVLV, start, kount)
            call check_err(iret)
         else
            call mapFullDomainToSubdomainMByNP(hs%ncid, &
                                               2*MNHARF, hs%myMesh%num_nodes, &
                                               hs%glovellv%u_nodal_data_id, subdomain_reals=GLOULV)
            call mapFullDomainToSubdomainMByNP(hs%ncid, &
                                               2*MNHARF, hs%myMesh%num_nodes, &
                                               hs%glovellv%v_nodal_data_id, subdomain_reals=GLOVLV)
         end if
      end if
!     STAELV - station elevation
      if (NHASE /= 0) then
         fullDomainIndexList => fullDomainElevationStationList
         iret = nf90_inq_varid(hs%ncid, "staelv", hs%staelv%station_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            kount(2) = hs%staelv%num_stations ! for elevation stations
            iret = nf90_get_var(hs%ncid, hs%staelv%station_data_id, &
                                STAELV, start, kount)
            call check_err(iret)
         else
            ! jgf52.30.12: Fix segfault by checking number of stations
            ! (provided by Damrongsak Wirasaet)
            if (NSTAE > 0) then
               call mapFullDomainToSubdomainMByNP(hs%ncid, &
                                                  2*MNHARF, hs%staelv%num_stations, &
                                                  hs%staelv%station_data_id, subdomain_reals=STAELV)
            end if
         end if
      end if
!     STAULV/STAVLV - station u and v velocity

      if (NHASV /= 0) then

         fullDomainIndexList => fullDomainVelocityStationList
         kount(2) = hs%stavellv%num_stations ! for velocity stations
         iret = nf90_inq_varid(hs%ncid, "staulv", &
                               hs%stavellv%u_station_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(hs%ncid, "stavlv", &
                               hs%stavellv%v_station_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_get_var(hs%ncid, hs%stavellv%u_station_data_id, &
                                STAULV, start, kount)
            call check_err(iret)
            iret = nf90_get_var(hs%ncid, hs%stavellv%v_station_data_id, &
                                STAVLV, start, kount)
            call check_err(iret)
         else
            ! jgf52.30.12: Fix segfault by checking number of stations
            ! (provided by Damrongsak Wirasaet)
            if (NSTAV > 0) then
               call mapFullDomainToSubdomainMByNP(hs%ncid, &
                                                  2*MNHARF, hs%stavellv%num_stations, &
                                                  hs%stavellv%u_station_data_id, subdomain_reals=STAULV)
               call mapFullDomainToSubdomainMByNP(hs%ncid, &
                                                  2*MNHARF, hs%stavellv%num_stations, &
                                                  hs%stavellv%v_station_data_id, subdomain_reals=STAVLV)
            end if
         end if
      end if

!     Read in model parameters to ADCIRC variables
      iret = nf90_inq_varid(hs%ncid, "icha", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, icha)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nz", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, inz)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nf", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, inf)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "mm", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, imm)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nstae", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, instae)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nstav", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, instav)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nhase", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, inhase)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nhasv", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, inhasv)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nhage", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, inhage)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nhagv", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, inhagv)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "icall", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, icall)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "nfreq", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, infreq)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "timeud", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, timeud)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "itud", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, itud)
      call check_err(iret)

!     Variables that are used to check that harmonic analysis data match
!     the simulation that is reading this hotstart file.
      inp = hs%myMesh%num_nodes

!     left hand side
      hakount(1) = 2*MNHARF
      hakount(2) = 2*MNHARF
      hastart(1) = 1
      hastart(2) = 1
      iret = nf90_inq_varid(hs%ncid, "ha", hs%ha_id)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, hs%ha_id, ha, hastart, hakount)
      call check_err(iret)
!     frequency names
      iret = nf90_inq_varid(hs%ncid, "namefr", hs%namefr_id)
      call check_err(iret)
      do i = 1, mnharf
         start(1) = 1
         start(2) = i
         kount(1) = len(inamefr(i))
         kount(2) = 1
         iret = nf90_get_var(hs%ncid, hs%namefr_id, inamefr(i), start, kount)
         call check_err(iret)
      end do
!     harmonic constituents
      start(1) = 1
      start(2) = 1
      kount(1) = MNHARF ! for constituents
      kount(2) = 1
      iret = nf90_inq_varid(hs%ncid, "hafreq", hs%hafreq_id)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, hs%hafreq_id, ifreq, start, kount)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "haff", hs%haff_id)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, hs%haff_id, iff, start, kount)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "haface", hs%haface_id)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, hs%haface_id, iface, start, kount)
      call check_err(iret)

!     now close the netcdf file
      iret = nf90_close(hs%ncid)
      call check_err(iret)

!-----------------------------------------------------------------------
   end subroutine readNetCDFHotstartHarmonic
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!        R E A D   N E T C D F   H O T S T A R T   H A R M O N I C
!                                   M E A N S   V A R I A N C E S
!-----------------------------------------------------------------------
!     jgf49.44.11 Reads harmonic analysis data from the hotstart file.
!-----------------------------------------------------------------------
   subroutine readNetCDFHotstartHarmonicMeansVariances(lun)
      use SIZES, only: MNPROC
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      use GLOBAL, only: OutputDataDescript_t
      use HARM, only: ELAV, ELVA, XVELAV, YVELAV, XVELVA, YVELVA, &
                      NTSTEPS, NHAGE, NHAGV

      implicit none

      integer, intent(in) :: lun

      integer :: kount(2), start(2)
      integer :: iret ! success or failure of the netcdf call
      integer :: tempid
      character(1024) :: scratchMessage

      LOG_SCOPE_TRACED("readNetCDFHotstartHarmonicMeansVariances", NETCDFIO_TRACING)

!     Point to the hotstart file we want to work on
      if (lun == 67) then
         hs => hs67
      else
         hs => hs68
      end if

      hs%myMesh => adcircMesh
      hs%myFile%fileFound = .false.

!     Open fulldomain file
      inquire (FILE=hs%myFile%FILENAME, EXIST=hs%myFile%fileFound)
      if (hs%myFile%fileFound .eqv. .false.) then
         write (scratchMessage, '(a,a,a)') "The file ", &
            trim(hs%myFile%FILENAME), " was not found; ADCIRC terminating."
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message=scratchMessage)
      else
         iret = nf90_open(hs%myFile%FILENAME, NF90_NOWRITE, hs%ncid)
         call check_err(iret)
      end if

!     Inquire variables
      call getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)

      hs%myFile%record_counter = 1

      kount(1) = hs%myMesh%num_nodes ! for nodal data
      kount(2) = hs%myTime%timenc_len
      start(1) = 1
      start(2) = hs%myFile%record_counter

      fullDomainIndexList => fullDomainNodeList
      if (NHAGE /= 0) then
         iret = nf90_inq_varid(hs%ncid, "elav", hs%elav%nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(hs%ncid, "elva", hs%elva%nodal_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_get_var(hs%ncid, hs%elav%nodal_data_id, &
                                ELAV, start, kount)
            call check_err(iret)
            iret = nf90_get_var(hs%ncid, hs%elva%nodal_data_id, &
                                ELVA, start, kount)
            call check_err(iret)
         else
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                          hs%elav%nodal_data_id, subdomain_reals=ELAV)
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                          hs%elva%nodal_data_id, subdomain_reals=ELVA)
         end if
      end if
      if (NHAGV /= 0) then
         iret = nf90_inq_varid(hs%ncid, "xvelav", hs%xvelav%nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(hs%ncid, "yvelav", hs%yvelav%nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(hs%ncid, "xvelva", hs%xvelva%nodal_data_id)
         call check_err(iret)
         iret = nf90_inq_varid(hs%ncid, "yvelva", hs%yvelva%nodal_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_get_var(hs%ncid, hs%xvelav%nodal_data_id, &
                                XVELAV, start, kount)
            call check_err(iret)
            iret = nf90_get_var(hs%ncid, hs%yvelav%nodal_data_id, &
                                YVELAV, start, kount)
            call check_err(iret)
            iret = nf90_get_var(hs%ncid, hs%xvelva%nodal_data_id, &
                                XVELVA, start, kount)
            call check_err(iret)
            iret = nf90_get_var(hs%ncid, hs%yvelva%nodal_data_id, &
                                YVELVA, start, kount)
            call check_err(iret)
         else
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                          hs%xvelav%nodal_data_id, subdomain_reals=XVELAV)
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                          hs%xvelva%nodal_data_id, subdomain_reals=XVELVA)
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                          hs%yvelav%nodal_data_id, subdomain_reals=YVELAV)
            call mapFullDomainToSubdomain(hs%ncid, hs%myMesh%num_nodes, &
                                          hs%yvelva%nodal_data_id, subdomain_reals=YVELVA)
         end if
      end if

!     Read in model parameters to ADCIRC variables
      iret = nf90_inq_varid(hs%ncid, "ntsteps", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, ntsteps)
      call check_err(iret)

!     now close the netcdf file
      iret = nf90_close(hs%ncid)
      call check_err(iret)

!-----------------------------------------------------------------------
   end subroutine readNetCDFHotstartHarmonicMeansVariances
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!        R E A D   N E T C D F   H O T S T A R T   3 D
!-----------------------------------------------------------------------
!     jgf49.48.03 Reads 3D data from the hotstart file.
!-----------------------------------------------------------------------
   subroutine readNetCDFHotstart3D(lun)
      use SIZES, only: MNPROC
      use GLOBAL, only: OutputDataDescript_t, IDEN, &
                        duu1, duv1, dvv1, uu2, vv2, bsx1, bsy1
      use GLOBAL_3DVS, only: q, nfen, sigt, sal, temp, wz, q20, &
                             l, n3dsd, i3dsdrec, n3dsv, i3dsvrec, n3dst, i3dstrec, n3dgd, &
                             i3dgdrec, n3dgv, i3dgvrec, n3dgt, i3dgtrec
      use MESH, only: NP
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none

      integer, intent(in) :: lun

      integer :: kount(2), start(2)
      integer :: kount3D(3), start3D(3)
      real(8), allocatable :: rp(:, :) ! real part of Q (subdomain), i.e. u-vel
      real(8), allocatable :: ip(:, :) ! imag part of Q (subdomain), i.e. v-vel
      integer :: iret ! success or failure of the netcdf call
      integer :: tempid
      character(1024) :: scratchMessage

      LOG_SCOPE_TRACED("readNetCDFHotstart3D", NETCDFIO_TRACING)

!     Point to the hotstart file we want to work on
      if (lun == 67) then
         hs => hs67
      else
         hs => hs68
      end if

      hs%myMesh => adcircMesh
      hs%myFile%fileFound = .false.

!     Open fulldomain file
      inquire (FILE=hs%myFile%FILENAME, EXIST=hs%myFile%fileFound)
      if (hs%myFile%fileFound .eqv. .false.) then
         write (scratchMessage, '(a,a,a)') "The file ", &
            trim(hs%myFile%FILENAME), " was not found; ADCIRC terminating."
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message=scratchMessage)
      else
         iret = nf90_open(hs%myFile%FILENAME, NF90_NOWRITE, hs%ncid)
         call check_err(iret)
      end if

!     Inquire variables
      call getDimensions(hs%ncid, hs%myTime, hs%myMesh, hs%myFile)
!     vertical node dimension
      iret = nf90_inq_dimid(hs%ncid, "num_v_nodes", &
                            hs%myMesh%num_v_nodes_dim_id)
      call check_err(iret)
      iret = nf90_inquire_dimension(hs%ncid, hs%myMesh%num_v_nodes_dim_id, &
                                    len=hs%myMesh%num_v_nodes)
      call check_err(iret)

      hs%myFile%record_counter = 1

      kount(1) = hs%myMesh%num_nodes ! for nodal data
      kount(2) = hs%myTime%timenc_len
      start(1) = 1
      start(2) = hs%myFile%record_counter

      kount3D(1) = hs%myMesh%num_nodes
      kount3D(2) = hs%myMesh%num_v_nodes ! for 3D data
      kount3D(3) = hs%myTime%timenc_len
      start3D(1) = 1
      start3D(2) = 1
      start3D(3) = hs%myFile%record_counter

      iret = nf90_inq_varid(hs%ncid, "duu", hs%duu%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "duv", hs%duv%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "dvv", hs%dvv%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "uu", hs%uu%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "vv", hs%vv%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "bsx", hs%bsx%nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "bsy", hs%bsy%nodal_data_id)
      call check_err(iret)
      if (MNPROC == 1) then
         iret = nf90_get_var(hs%ncid, hs%duu%nodal_data_id, &
                             DUU1, start, kount)
         call check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%duv%nodal_data_id, &
                             DUV1, start, kount)
         call check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%dvv%nodal_data_id, &
                             DVV1, start, kount)
         call check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%uu%nodal_data_id, &
                             UU2, start, kount)
         call check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%vv%nodal_data_id, &
                             VV2, start, kount)
         call check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%bsx%nodal_data_id, &
                             BSX1, start, kount)
         call check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%bsy%nodal_data_id, &
                             BSY1, start, kount)
         call check_err(iret)
      else
         fullDomainIndexList => fullDomainNodeList
         call mapFullDomainToSubdomain(hs%ncid, &
                                       hs%myMesh%num_nodes, hs%duu%nodal_data_id, &
                                       subdomain_reals=DUU1)
         call mapFullDomainToSubdomain(hs%ncid, &
                                       hs%myMesh%num_nodes, hs%duv%nodal_data_id, &
                                       subdomain_reals=DUV1)
         call mapFullDomainToSubdomain(hs%ncid, &
                                       hs%myMesh%num_nodes, hs%dvv%nodal_data_id, &
                                       subdomain_reals=DVV1)
         call mapFullDomainToSubdomain(hs%ncid, &
                                       hs%myMesh%num_nodes, hs%uu%nodal_data_id, &
                                       subdomain_reals=UU2)
         call mapFullDomainToSubdomain(hs%ncid, &
                                       hs%myMesh%num_nodes, hs%vv%nodal_data_id, &
                                       subdomain_reals=VV2)
         call mapFullDomainToSubdomain(hs%ncid, &
                                       hs%myMesh%num_nodes, hs%bsx%nodal_data_id, &
                                       subdomain_reals=BSX1)
         call mapFullDomainToSubdomain(hs%ncid, &
                                       hs%myMesh%num_nodes, hs%bsy%nodal_data_id, &
                                       subdomain_reals=BSY1)
      end if
      !
      ! 3D Density
      if (abs(IDEN) == 1) then
         iret = nf90_inq_varid(hs%ncid, "sigt", hs%density3D%u_nodal_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_get_var(hs%ncid, &
                                hs%density3D%u_nodal_data_id, SIGT, start3D, kount3D)
            call check_err(iret)
         else
            call mapFullDomainToSubdomainNPByM(hs%ncid, &
                                               nfen, hs%myMesh%num_nodes, hs%density3D%u_nodal_data_id, &
                                               subdomain_reals=SIGT)
         end if
      end if
      if ((abs(IDEN) == 2) .or. (abs(IDEN) == 4)) then
         iret = nf90_inq_varid(hs%ncid, "salinity", &
                               hs%density3D%v_nodal_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_get_var(hs%ncid, hs%density3D%v_nodal_data_id, &
                                SAL, start3D, kount3D)
            call check_err(iret)
         else
            call mapFullDomainToSubdomainNPByM(hs%ncid, &
                                               nfen, hs%myMesh%num_nodes, hs%density3D%v_nodal_data_id, &
                                               subdomain_reals=SAL)
         end if
      end if
      if ((abs(IDEN) == 3) .or. (abs(IDEN) == 4)) then
         iret = nf90_inq_varid(hs%ncid, "temperature", &
                               hs%density3D%w_nodal_data_id)
         call check_err(iret)
         if (MNPROC == 1) then
            iret = nf90_get_var(hs%ncid, hs%density3D%w_nodal_data_id, &
                                TEMP, start3D, kount3D)
            call check_err(iret)
         else
            call mapFullDomainToSubdomainNPByM(hs%ncid, &
                                               nfen, hs%myMesh%num_nodes, hs%density3D%w_nodal_data_id, &
                                               subdomain_reals=TEMP)
         end if
      end if

      ! 3D velocity
      iret = nf90_inq_varid(hs%ncid, "u-vel3D", hs%velocity3D%u_nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "v-vel3D", hs%velocity3D%v_nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "w-vel3D", hs%velocity3D%w_nodal_data_id)
      call check_err(iret)
      allocate (rp(np, nfen), ip(np, nfen))
      if (MNPROC == 1) then
         iret = nf90_get_var(hs%ncid, hs%velocity3D%u_nodal_data_id, &
                             rp, start3D, kount3D)
         call check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%velocity3D%v_nodal_data_id, &
                             ip, start3D, kount3D)
         call check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%velocity3D%w_nodal_data_id, &
                             wz, start3D, kount3D)
         call check_err(iret)
      else
         call mapFullDomainToSubdomainNPByM(hs%ncid, &
                                            nfen, hs%myMesh%num_nodes, hs%velocity3D%u_nodal_data_id, &
                                            subdomain_reals=rp)
         call mapFullDomainToSubdomainNPByM(hs%ncid, &
                                            nfen, hs%myMesh%num_nodes, hs%velocity3D%v_nodal_data_id, &
                                            subdomain_reals=ip)
         call mapFullDomainToSubdomainNPByM(hs%ncid, &
                                            nfen, hs%myMesh%num_nodes, hs%velocity3D%w_nodal_data_id, &
                                            subdomain_reals=wz)
      end if
      q(:, :) = cmplx(rp(:, :), ip(:, :), 8) ! construct q from real and imaginary
      deallocate (rp, ip)

      ! 3D turbulence
      iret = nf90_inq_varid(hs%ncid, "q20", hs%turbulence3D%u_nodal_data_id)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "l", hs%turbulence3D%v_nodal_data_id)
      call check_err(iret)
      if (MNPROC == 1) then
         iret = nf90_get_var(hs%ncid, hs%turbulence3D%u_nodal_data_id, &
                             q20, start3D, kount3D)
         call check_err(iret)
         iret = nf90_get_var(hs%ncid, hs%turbulence3D%v_nodal_data_id, &
                             l, start3D, kount3D)
         call check_err(iret)
      else
         call mapFullDomainToSubdomainNPByM(hs%ncid, &
                                            nfen, hs%myMesh%num_nodes, hs%turbulence3D%u_nodal_data_id, &
                                            subdomain_reals=q20)
         call mapFullDomainToSubdomainNPByM(hs%ncid, &
                                            nfen, hs%myMesh%num_nodes, hs%velocity3D%v_nodal_data_id, &
                                            subdomain_reals=l)
      end if

!     Read in model parameters to ADCIRC variables
      iret = nf90_inq_varid(hs%ncid, "n3dsd", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, n3dsd)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "i3dsdrec", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, i3dsdrec)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "n3dsv", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, n3dsv)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "i3dsvrec", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, i3dsvrec)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "n3dst", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, n3dst)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "i3dstrec", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, i3dstrec)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "n3dgd", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, n3dgd)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "i3dgdrec", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, i3dgdrec)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "n3dgv", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, n3dgv)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "i3dgvrec", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, i3dgvrec)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "n3dgt", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, n3dgt)
      call check_err(iret)
      iret = nf90_inq_varid(hs%ncid, "i3dgtrec", tempid)
      call check_err(iret)
      iret = nf90_get_var(hs%ncid, tempid, i3dgtrec)
      call check_err(iret)

!     now close the netcdf file
      iret = nf90_close(hs%ncid)
      call check_err(iret)

!-----------------------------------------------------------------------
   end subroutine readNetCDFHotstart3D
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E   D E F I N E   M E S H   V A R I A B L E S
!-----------------------------------------------------------------------
!     jgf49.17.02 Defines data that are common to all netCDF files.
!-----------------------------------------------------------------------
   subroutine defineMeshVariables(ncid, myMesh, myFile)
      use GLOBAL, only: C3D
      use GLOBAL_3DVS, only: NFEN
      implicit none
      integer :: iret ! Error status return
      type(meshStructure), intent(inout) :: myMesh
      integer, intent(in) :: ncid ! netCDF id
      type(fileData), intent(in) :: myFile ! file format info
      character(1024) :: att_text ! reused to hold attribute text
!     -----------------
!      date_string variables for time attribute
!     -----------------
      integer :: dmy

      LOG_SCOPE_TRACED("defineMeshVariables", NETCDFIO_TRACING)

!     Define dimensions
      iret = nf90_def_dim(ncid, 'node', myMesh%num_nodes, &
                          myMesh%num_nodes_dim_id)
      call check_err(iret)
      iret = nf90_def_dim(ncid, 'nele', myMesh%num_elems, dmy)
      call check_err(iret)
      iret = nf90_def_dim(ncid, 'nfaces', myMesh%num_elems, &
                          myMesh%num_elems_dim_id)
      call check_err(iret)
      if (C3D .eqv. .true.) then
         myMesh%num_v_nodes = NFEN
         iret = nf90_def_dim(ncid, 'num_v_nodes', myMesh%num_v_nodes, &
                             myMesh%num_v_nodes_dim_id)
         call check_err(iret)
         myMesh%sigma_dims(1) = myMesh%num_v_nodes_dim_id
         iret = nf90_def_var(ncid, 'sigma', NF90_DOUBLE, &
                             myMesh%sigma_dims, myMesh%sigma_id)
         call check_err(iret)
      end if
      iret = nf90_def_dim(ncid, 'nvertex', 3, myMesh%nface_dim_id)
      call check_err(iret)
      if (myMesh%nopenc /= 0) then
         iret = nf90_def_dim(ncid, 'nope', myMesh%nopenc, &
                             myMesh%nopenc_dim_id)
         call check_err(iret)
!        WJP to make arrays smaller only have dimension of neta
         iret = nf90_def_dim(ncid, 'neta', myMesh%netanc, &
                             myMesh%netanc_dim_id)
         call check_err(iret)
         iret = nf90_def_dim(ncid, 'max_nvdll', myMesh%max_nvdllnc, &
                             myMesh%max_nvdllnc_dim_id)
         call check_err(iret)
      end if
      if (myMesh%nbounc /= 0) then
         iret = nf90_def_dim(ncid, 'nbou', &
                             myMesh%nbounc, myMesh%nbounc_dim_id)
         call check_err(iret)
!        WJP to make arrays smaller only have dimension of nvel
         iret = nf90_def_dim(ncid, 'nvel', myMesh%nvelnc, &
                             myMesh%nvelnc_dim_id)
         call check_err(iret)
         iret = nf90_def_dim(ncid, 'max_nvell', &
                             myMesh%max_nvellnc, myMesh%max_nvellnc_dim_id)
      end if

!     Define variables
!     Define X
      myMesh%X_dims(1) = myMesh%num_nodes_dim_id
      iret = nf90_def_var(ncid, 'x', NF90_DOUBLE, &
                          myMesh%X_dims, myMesh%X_id)
      call check_err(iret)
!     Define Y coordinate
      myMesh%Y_dims(1) = myMesh%num_nodes_dim_id
      iret = nf90_def_var(ncid, 'y', NF90_DOUBLE, &
                          myMesh%Y_dims, myMesh%Y_id)
      call check_err(iret)

!     Define elements
      myMesh%ELE_dims(1) = myMesh%nface_dim_id
      myMesh%ELE_dims(2) = myMesh%num_elems_dim_id
      iret = nf90_def_var(ncid, 'element', NF90_INT, &
                          myMesh%ELE_dims, myMesh%ELE_id)
      call check_err(iret)

      ! Corbitt: Define ADCIRC-Mesh Variable
      call check_err(nf90_def_var(ncid, 'adcirc_mesh', NF90_INT, myMesh%MESH_id))

      ! jgf50.44: Turn on compression if this is a netcdf4-formatted file.
#ifdef NETCDF_CAN_DEFLATE
      if (myFile%ncformat == ior(NF90_CLASSIC_MODEL, NF90_NETCDF4)) then
         iret = nf90_def_var_deflate(ncid, myMesh%X_id, 1, 1, 2)
         call check_err(iret)
         iret = nf90_def_var_deflate(ncid, myMesh%Y_id, 1, 1, 2)
         call check_err(iret)
         iret = nf90_def_var_deflate(ncid, myMesh%ELE_id, 1, 1, 2)
         call check_err(iret)
      end if
#endif

!     Define elevation specified boundary forcing segments information
      if (myMesh%nopenc /= 0) then
         myMesh%nvdllnc_dims(1) = myMesh%nopenc_dim_id
         iret = nf90_def_var(ncid, 'nvdll', NF90_INT, &
                             myMesh%nvdllnc_dims, myMesh%nvdllnc_id)
         call check_err(iret)
         myMesh%ibtypeenc_dims(1) = myMesh%nopenc_dim_id
         iret = nf90_def_var(ncid, 'ibtypee', NF90_INT, &
                             myMesh%ibtypeenc_dims, &
                             myMesh%ibtypeenc_id)
         call check_err(iret)
!        WJP to make arrays smaller only have dimension of neta
         myMesh%nbdvnc_dims(1) = myMesh%netanc_dim_id
         iret = nf90_def_var(ncid, 'nbdv', NF90_INT, &
                             myMesh%nbdvnc_dims(1), myMesh%nbdvnc_id)
         call check_err(iret)
!        !jgf50.44: Turn on compression if this is a netcdf4 formatted file.
#ifdef NETCDF_CAN_DEFLATE
         if (myFile%ncformat == ior(NF90_CLASSIC_MODEL, NF90_NETCDF4)) then
            iret = nf90_def_var_deflate(ncid, myMesh%nvdllnc_id, &
                                        1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(ncid, myMesh%ibtypeenc_id, &
                                        1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(ncid, myMesh%nbdvnc_id, &
                                        1, 1, 2)
            call check_err(iret)
         end if
#endif
      end if
!     Define normal flow boundary information
      if (myMesh%nbounc /= 0) then
         myMesh%nvellnc_dims(1) = myMesh%nbounc_dim_id
         iret = nf90_def_var(ncid, 'nvell', NF90_INT, &
                             myMesh%nvellnc_dims, myMesh%nvellnc_id)
         call check_err(iret)
         myMesh%ibtypenc_dims(1) = myMesh%nbounc_dim_id
         iret = nf90_def_var(ncid, 'ibtype', NF90_INT, &
                             myMesh%ibtypenc_dims, myMesh%ibtypenc_id)
         call check_err(iret)
!        WJP to make arrays smaller only have dimension of nvel
         myMesh%nbvvnc_dims(1) = myMesh%nvelnc_dim_id
         iret = nf90_def_var(ncid, 'nbvv', NF90_INT, &
                             myMesh%nbvvnc_dims(1), myMesh%nbvvnc_id)
         call check_err(iret)
         if (mymesh%hasWeirs) then
            iret = nf90_def_var(ncid, 'barsp', NF90_DOUBLE, &
                                myMesh%nbvvnc_dims(1), myMesh%barsp_id)
            call check_err(iret)
            iret = nf90_def_var(ncid, 'barht', NF90_DOUBLE, &
                                myMesh%nbvvnc_dims(1), myMesh%barht_id)
            if (myMesh%hasInternalWeirs) then
               call check_err(iret)
               iret = nf90_def_var(ncid, 'ibconn', NF90_INT, &
                                   myMesh%nbvvnc_dims(1), myMesh%ibconn_id)
               call check_err(iret)
               iret = nf90_def_var(ncid, 'barsb', NF90_DOUBLE, &
                                   myMesh%nbvvnc_dims(1), myMesh%barsb_id)
               call check_err(iret)
               if (myMesh%hasPipes) then
                  iret = nf90_def_var(ncid, 'pipeht', NF90_DOUBLE, &
                                      myMesh%nbvvnc_dims(1), myMesh%pipeht_id)
                  call check_err(iret)
                  iret = nf90_def_var(ncid, 'pipediam', NF90_DOUBLE, &
                                      myMesh%nbvvnc_dims(1), myMesh%pipediam_id)
                  call check_err(iret)
                  iret = nf90_def_var(ncid, 'pipecoef', NF90_DOUBLE, &
                                      myMesh%nbvvnc_dims(1), myMesh%pipecoef_id)
                  call check_err(iret)
               end if
            end if
         end if

!        !jgf50.44: Turn on compression if this is a netcdf4 formatted file.
#ifdef NETCDF_CAN_DEFLATE
         if (myFile%ncformat == ior(NF90_CLASSIC_MODEL, NF90_NETCDF4)) then
            iret = nf90_def_var_deflate(ncid, myMesh%nvellnc_id, &
                                        1, 1, 2)
            call check_err(iret)
            iret = nf90_def_var_deflate(ncid, myMesh%ibtypenc_id, &
                                        1, 1, 2)
            call check_err(iret)
            if (myMesh%hasWeirs) then
               iret = nf90_def_var_deflate(ncid, myMesh%barht_id, &
                                           1, 1, 2)
               call check_err(iret)
               iret = nf90_def_var_deflate(ncid, myMesh%nbvvnc_id, &
                                           1, 1, 2)
               call check_err(iret)
               iret = nf90_def_var_deflate(ncid, myMesh%barsp_id, &
                                           1, 1, 2)
               call check_err(iret)
               if (myMesh%hasInternalWeirs) then
                  iret = nf90_def_var_deflate(ncid, myMesh%ibconn_id, &
                                              1, 1, 2)
                  call check_err(iret)
                  iret = nf90_def_var_deflate(ncid, myMesh%barsb_id, &
                                              1, 1, 2)
                  call check_err(iret)
                  if (myMesh%hasPipes) then
                     iret = nf90_def_var_deflate(ncid, myMesh%pipeht_id, &
                                                 1, 1, 2)
                     call check_err(iret)
                     iret = nf90_def_var_deflate(ncid, myMesh%pipediam_id, &
                                                 1, 1, 2)
                     call check_err(iret)
                     iret = nf90_def_var_deflate(ncid, myMesh%pipecoef_id, &
                                                 1, 1, 2)
                     call check_err(iret)
                  end if
               end if
            end if
         end if
#endif
      end if
!     -------------------
!     Define Z coordinate
!     --------------------
      myMesh%DEPTH_dims(1) = myMesh%num_nodes_dim_id
      iret = nf90_def_var(ncid, 'depth', NF90_DOUBLE, &
                          myMesh%DEPTH_dims, myMesh%DEPTH_id)
      call check_err(iret)
#ifdef NETCDF_CAN_DEFLATE
      if (myFile%ncformat == ior(NF90_CLASSIC_MODEL, NF90_NETCDF4)) then
         iret = nf90_def_var_deflate(ncid, myMesh%DEPTH_id, &
                                     1, 1, 2)
         call check_err(iret)
      end if
#endif

!     Set coordinates as representing latitude or longitude, depending on
!     the value of ICS
      call defineCoordinateAttributes(ncid, myMesh%X_id, myMesh%Y_id)

!     Define depth attributes (Corbitt)
      iret = nf90_put_att(ncid, myMesh%DEPTH_id, 'long_name', 'distance&
&  below geoid')
      call check_err(iret)
      iret = nf90_put_att(ncid, myMesh%DEPTH_id, 'standard_name', &
                          'depth below geoid')
      call check_err(iret)
      iret = nf90_put_att(ncid, myMesh%DEPTH_id, 'coordinates', &
                          'time y x')
      call check_err(iret)
      iret = nf90_put_att(ncid, myMesh%DEPTH_id, 'location', 'node')
      call check_err(iret)
      iret = nf90_put_att(ncid, myMesh%DEPTH_id, 'mesh', 'adcirc_mesh')
      call check_err(iret)
      iret = nf90_put_att(ncid, myMesh%DEPTH_id, 'units', 'm')
      call check_err(iret)
!      iret = nf90_put_att(ncid, myMesh%DEPTH_id, 'positive', 'down')
!      CALL check_err(iret) Corbitt: Not CF Compliant

!     Define Element Attributes (Corbitt)
      iret = nf90_put_att(ncid, myMesh%ELE_id, 'long_name', 'element')
      call check_err(iret)
!Corbitt 11/12/13 - Obsolete CF-UGRID Convention
!      iret = nf90_put_att(ncid, myMesh%ELE_id,'standard_name',
!     & 'face_node_connectivity')
!      CALL check_err(iret)
      iret = nf90_put_att(ncid, myMesh%ELE_id, 'cf_role', &
                          'face_node_connectivity')
      call check_err(iret)
      iret = nf90_put_att(ncid, myMesh%ELE_id, 'start_index', 1)
      call check_err(iret)
      iret = nf90_put_att(ncid, myMesh%ELE_id, 'units', 'nondimensional')
      call check_err(iret)

!     Define text attributes for boundary segments
      if (myMesh%nopenc /= 0) then
!        NOPE
!         att_text = "number of elevation specified boundary &
!     &forcing segments"
!         iret = nf90_put_att(ncid, myMesh%nopenc_id, 'long_name',
!     &          len_trim(att_text), trim(att_text))
!         CALL check_err(iret)
!         iret = nf90_put_att(ncid, myMesh%nopenc_id, 'units', 14,
!     &                     'nondimensional')
!         CALL check_err(iret)
!        NVDLL
         att_text = "number of nodes in each elevation specified " &
                    //"boundary segment"
         iret = nf90_put_att(ncid, myMesh%nvdllnc_id, 'long_name', &
                             trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(ncid, myMesh%nvdllnc_id, 'units', &
                             'nondimensional')
         call check_err(iret)
!        IBTYPEE
         att_text = "elevation boundary type"
         iret = nf90_put_att(ncid, myMesh%ibtypeenc_id, 'long_name', &
                             trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(ncid, myMesh%ibtypeenc_id, 'units', &
                             'nondimensional')
         call check_err(iret)
!        NBDV
         att_text = "node numbers on each elevation specified boundary " &
                    //"segment"
         iret = nf90_put_att(ncid, myMesh%nbdvnc_id, 'long_name', &
                             trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(ncid, myMesh%nbdvnc_id, 'units', &
                             'nondimensional')
         call check_err(iret)
      end if
      if (myMesh%nbounc /= 0) then
!        IBTYPE
         att_text = "type of normal flow (discharge) boundary"
         iret = nf90_put_att(ncid, myMesh%ibtypenc_id, 'long_name', &
                             trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(ncid, myMesh%ibtypenc_id, 'units', &
                             'nondimensional')
         call check_err(iret)
!        NVELL
         att_text = 'number of nodes in each normal flow ' &
                    //'specified boundary segment'
         iret = nf90_put_att(ncid, myMesh%nvellnc_id, 'long_name', &
                             trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(ncid, myMesh%nvellnc_id, 'units', &
                             'nondimensional')
         call check_err(iret)
!        NBVV
         att_text = "node numbers on normal flow boundary segment"
         iret = nf90_put_att(ncid, myMesh%nbvvnc_id, 'long_name', &
                             trim(att_text))
         call check_err(iret)
         iret = nf90_put_att(ncid, myMesh%nbvvnc_id, 'units', &
                             'nondimensional')
         call check_err(iret)
         if (myMesh%hasWeirs) then
!           BARHT
            att_text = &
               "external/internal barrier boundary height above geoid"
            iret = nf90_put_att(ncid, myMesh%barht_id, 'long_name', &
                                trim(att_text))
            iret = nf90_put_att(ncid, myMesh%barht_id, 'units', 'm')
            call check_err(iret)
!           BARSP
            att_text = "external/internal barrier boundary"// &
                       " weir coefficient of supercritical flow"
            iret = nf90_put_att(ncid, myMesh%barsp_id, "long_name", &
                                trim(att_text))
            call check_err(iret)
            iret = nf90_put_att(ncid, mymesh%barsp_id, 'units', &
                                'nondimensional')
            call check_err(iret)
            if (myMesh%hasInternalWeirs) then
!               IBCONN
               att_text = "internal barrier nodal connections"
               iret = nf90_put_att(ncid, myMesh%ibconn_id, 'long_name', &
                                   trim(att_text))
               call check_err(iret)
               iret = nf90_put_att(ncid, mymesh%ibconn_id, 'units', &
                                   'nondimensional')
               call check_err(iret)
!               BARSB
               att_text = "internal barrier boundary"// &
                          " weir coefficient of subcritical flow"
               iret = nf90_put_att(ncid, myMesh%barsb_id, "long_name", &
                                   trim(att_text))
               call check_err(iret)
               iret = nf90_put_att(ncid, mymesh%barsb_id, 'units', &
                                   'nondimensional')
               call check_err(iret)
               if (myMesh%hasPipes) then
!                   PIPEHT
                  att_text = &
                     "internal barrier boundary cross-barrier pipe height"
                  iret = nf90_put_att(ncid, myMesh%pipeht_id, "long_name", &
                                      trim(att_text))
                  iret = nf90_put_att(ncid, myMesh%pipeht_id, 'units', 'm')
                  call check_err(iret)
!                   PIPEDIAM
                  att_text = &
                     "internal barrier boundary cross-barrier pipe diameter"
                  iret = nf90_put_att(ncid, myMesh%pipediam_id, "long_name", &
                                      trim(att_text))
                  iret = nf90_put_att(ncid, myMesh%pipediam_id, 'units', 'm')
                  call check_err(iret)
!                   PIPECOEF
                  att_text = &
                     "internal barrier boundary cross-barrier pipe coefficient"
                  iret = nf90_put_att(ncid, myMesh%pipecoef_id, "long_name", &
                                      trim(att_text))
                  iret = nf90_put_att(ncid, mymesh%pipecoef_id, 'units', &
                                      'nondimensional')
                  call check_err(iret)
               end if
            end if
         end if
      end if
!Corbitt: Define ADCIRC-Mesh Attributes
      iret = nf90_put_att(ncid, myMesh%mesh_id, 'long_name', &
                          'mesh_topology')
      call check_err(iret)
      iret = nf90_put_att(ncid, myMesh%mesh_id, 'cf_role', &
                          'mesh_topology')
      call check_err(iret)
!Corbitt: 11/12/13 - Obsolete CF-UGRID Convention
!         iret = nf90_put_att(ncid, myMesh%mesh_id, 'standard_name',
!     &                      'mesh_topology')
!         CALL check_err(iret)
!         iret = nf90_put_att(ncid, myMesh%mesh_id, 'dimension',
!     &                      2)
!         CALL check_err(iret)
      iret = nf90_put_att(ncid, myMesh%mesh_id, 'topology_dimension', &
                          2)
      call check_err(iret)
      iret = nf90_put_att(ncid, myMesh%mesh_id, 'node_coordinates', &
                          'x y')
      call check_err(iret)
      iret = nf90_put_att(ncid, myMesh%mesh_id, &
                          'face_node_connectivity', 'element')
      call check_err(iret)

!-----------------------------------------------------------------------
   end subroutine defineMeshVariables
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!        D E F I N E  C O O R D I N A T E  A T T R I B U T E S
!-----------------------------------------------------------------------
!     jgf49.21 Defines coordinate attributes to identify coordinates as
!     either x and y in Cartesian length (feet or meters) or latitude
!     and longitude depending on the value of ICS.
!-----------------------------------------------------------------------
   subroutine defineCoordinateAttributes(ncid, xid, yid)
      use MESH, only: ICS
      implicit none
      integer, intent(in) :: ncid
      integer, intent(in) :: xid
      integer, intent(in) :: yid
      integer :: iret

      LOG_SCOPE_TRACED("defineCoordinateAttributes", NETCDFIO_TRACING)

      if (ics /= 1) then ! this indicates spherical coordinates

!        Define longitude attributes
         iret = nf90_put_att(ncid, xid, 'long_name', 'longitude')
         call check_err(iret)
         iret = nf90_put_att(ncid, xid, 'standard_name', 'longitude')
         call check_err(iret)
         iret = nf90_put_att(ncid, xid, 'units', 'degrees_east')
         call check_err(iret)
         iret = nf90_put_att(ncid, xid, 'positive', 'east')
         call check_err(iret)
!        Define latitude attributes
         iret = nf90_put_att(ncid, yid, 'long_name', 'latitude')
         call check_err(iret)
         iret = nf90_put_att(ncid, yid, 'standard_name', 'latitude')
         call check_err(iret)
         iret = nf90_put_att(ncid, yid, 'units', 'degrees_north')
         call check_err(iret)
         iret = nf90_put_att(ncid, yid, 'positive', 'north')
         call check_err(iret)
      else ! must be using Cartesian (x,y) coordinates
!        Define x-coordinate attributes
         iret = nf90_put_att(ncid, xid, 'long_name', &
                             'Cartesian coordinate x')
         call check_err(iret)
         iret = nf90_put_att(ncid, xid, 'standard_name', 'x_coordinate')
         call check_err(iret)
!        determine variable units
         call putUnitsAttribute(ncid, xid, 'm', 'ft')
         iret = nf90_put_att(ncid, xid, 'positive', 'right')
         call check_err(iret)
!        Define y-coordinate attributes
         iret = nf90_put_att(ncid, yid, 'long_name', &
                             'Cartesian coordinate y')
         call check_err(iret)
         iret = nf90_put_att(ncid, yid, 'standard_name', 'y_coordinate')
         call check_err(iret)
         call putUnitsAttribute(ncid, yid, 'm', 'ft')
         iret = nf90_put_att(ncid, yid, 'positive', &
                             '90 degrees counterclockwise from x')
         call check_err(iret)
      end if

!-----------------------------------------------------------------------
   end subroutine defineCoordinateAttributes
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E   P U T   M E S H   V A R I A B L E S
!-----------------------------------------------------------------------
!     jgf49.17.02 Writes data that are common to all netCDF files
!     to the specified file.
!-----------------------------------------------------------------------
   subroutine putMeshVariables(ncid, myMesh)
      use MESH, only: DP
      use GLOBAL, only: C3D
      use GLOBAL_3DVS, only: SIGMA
      implicit none
      integer, intent(in) :: ncid
      type(meshStructure), intent(inout) :: myMesh
      integer :: iret ! Error status return

      integer :: kount(2), start(2)

      LOG_SCOPE_TRACED("putMeshVariables", NETCDFIO_TRACING)

!     Store nodal coordinates
      iret = nf90_put_var(ncid, myMesh%X_id, myMesh%xnc)
      call check_err(iret)
      iret = nf90_put_var(ncid, myMesh%Y_id, myMesh%ync)
      call check_err(iret)
      if (C3D .eqv. .true.) then
         iret = nf90_put_var(ncid, myMesh%sigma_id, sigma)
         call check_err(iret)
      end if

!     Store depth
      iret = nf90_put_var(ncid, myMesh%DEPTH_id, DP)
      call check_err(iret)
!     Store elements
      kount(1) = myMesh%nface_len
      kount(2) = myMesh%num_elems
      start(1) = 1
      start(2) = 1
      iret = nf90_put_var(ncid, myMesh%ele_id, &
                          myMesh%element, start, kount)
      call check_err(iret)

!     Store elevation boundary information
      if (myMesh%nopenc /= 0) then
         iret = nf90_put_var(ncid, myMesh%nvdllnc_id, myMesh%nvdllnc)
         call check_err(iret)
         iret = nf90_put_var(ncid, myMesh%ibtypeenc_id, &
                             myMesh%ibtypeenc)
         call check_err(iret)
         iret = nf90_put_var(ncid, myMesh%nbdvnc_id, myMesh%nbdvnc)
         call check_err(iret)
      end if
!     Store normal flow boundary information
      if (myMesh%nbounc /= 0) then
         iret = nf90_put_var(ncid, myMesh%ibtypenc_id, myMesh%ibtypenc)
         call check_err(iret)
         iret = nf90_put_var(ncid, myMesh%nvellnc_id, myMesh%nvellnc)
         call check_err(iret)
         iret = nf90_put_var(ncid, myMesh%nbvvnc_id, myMesh%nbvvnc)
         call check_err(iret)
         if (myMesh%hasWeirs) then
            iret = nf90_put_var(ncid, myMesh%barht_id, myMesh%barht)
            call check_err(iret)
            iret = nf90_put_var(ncid, myMesh%barsp_id, myMesh%barsp)
            call check_err(iret)
            if (myMesh%hasInternalWeirs) then
               iret = nf90_put_var(ncid, myMesh%ibconn_id, myMesh%ibconnnc)
               call check_err(iret)
               iret = nf90_put_var(ncid, myMesh%barsb_id, myMesh%barsb)
               call check_err(iret)
               if (myMesh%hasPipes) then
                  iret = nf90_put_var(ncid, myMesh%pipeht_id, myMesh%pipeht)
                  call check_err(iret)
                  iret = nf90_put_var(ncid, myMesh%pipediam_id, myMesh%pipediam)
                  call check_err(iret)
                  iret = nf90_put_var(ncid, myMesh%pipecoef_id, myMesh%pipecoef)
                  call check_err(iret)
               end if
            end if
         end if
      end if

!-----------------------------------------------------------------------
   end subroutine putMeshVariables
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E   D E F I N E   T I M E   A T T R I B U T E S
!-----------------------------------------------------------------------
!     jgf50.38 Writes time data that are common to all netCDF files
!     to the specified file.
!-----------------------------------------------------------------------
   subroutine defineTimeAttributes(ncid, myTime)
      use GLOBAL, only: base_date
      implicit none
      integer, intent(in) :: ncid ! file id of netcdf file
      type(TimeData), intent(in) :: myTime ! time data for this netcdf file

      character :: date_string*60
      integer :: iret ! error status of netcdf call

      iret = nf90_put_att(ncid, myTime%timenc_id, &
                          'long_name', 'model time')
      call check_err(iret)
      iret = nf90_put_att(ncid, myTime%timenc_id, &
                          'standard_name', 'time')
      call check_err(iret)
      date_string = 'seconds since '//adjustl(trim(base_date))
      iret = nf90_put_att(ncid, myTime%timenc_id, 'units', &
                          trim(date_string))
      call check_err(iret)
      iret = nf90_put_att(ncid, myTime%timenc_id, &
                          'base_date', adjustl(trim(base_date)))
      call check_err(iret)
!-----------------------------------------------------------------------
   end subroutine defineTimeAttributes
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E   D E F I N E   M E T A D A T A
!-----------------------------------------------------------------------
!     jgf49.17.02 Writes data that are common to all netCDF files
!     to the specified file.
!-----------------------------------------------------------------------
   subroutine defineMetaData(ncid)
      use VERSION, only: ADC_VERSION, ADC_HASH
      use MESH, only: ICS, AGRID, SLAM0, SFEA0
      use ADC_CONSTANTS, only: RAD2DEG
      use GLOBAL, only: RUNDES, title, institution, source, &
                        history, references, comments, host, &
                        convention, contact, dtdp, ihot, &
                        nolifa, nolica, nolicat, &
                        ncor, ntip, nws, nramp, statim, &
                        reftim, rnday, dramp, h0, &
                        cori, ntif, nbfr, C3D, inunThresh, &
                        inundationOutput, runid
      use GLOBAL_3DVS, only: iden, islip, kp, z0s, z0b, theta1, theta2, &
                             ievc, evmin, evcon, alp1, alp2, alp3, igc, nlsd, nvsd, nltd, &
                             nvtd, alp4
      use GWCE, only: a00, b00, c00
      use NodalAttributes, only: nolibf, nwp, tau0, cf, eslm
      implicit none
      integer, intent(in) :: ncid
      integer :: iret ! success or failure of the netcdf call

      real(8) :: SLAM0DEG
      real(8) :: SFEA0DEG
!     date_string variables for time attribute
      character :: date_string*40
      character :: now_date*8
      character :: big_ben*10
      character :: zone*5
      integer :: values(8)

      LOG_SCOPE_TRACED("defineMetaData", NETCDFIO_TRACING)

!     Convert back to degrees ... the original input is in degrees,
!     but this gets converted to radians immediately and unfortunately
!     the values that were read in get overwritten ... need to go
!     back to degrees to write them back out
      SLAM0DEG = SLAM0*RAD2DEG
      SFEA0DEG = SFEA0*RAD2DEG
!     -----------------
!     Global attributes
!     -----------------
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'model', 'ADCIRC')
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'version', ADC_VERSION)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'git_hash', ADC_HASH)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'grid_type', 'Triangular')
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'description', &
                          trim(adjustl(rundes)))
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'agrid', &
                          trim(adjustl(agrid)))
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'rundes', &
                          trim(adjustl(rundes)))
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'runid', &
                          trim(adjustl(runid)))
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'title', &
                          trim(adjustl(title)))
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'institution', &
                          trim(adjustl(institution)))
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'source', &
                          trim(adjustl(source)))
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'history', &
                          trim(adjustl(history)))
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'references', &
                          trim(adjustl(references)))
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'comments', &
                          trim(adjustl(comments)))
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'host', &
                          trim(adjustl(host)))
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'convention', &
                          trim(adjustl(convention)))
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', &
                          'UGRID-0.9.0')
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'contact', &
                          trim(adjustl(contact)))
      call check_err(iret)
      call date_and_time(now_date, big_ben, zone, values)
      write (date_string, 71) values(1), values(2), values(3), &
         values(5), values(6), values(7), (values(4))/60
71    format(I4, '-', I2.2, '-', i2.2, ' ', i2, ':', i2.2, ':', i2.2, ' ' &
             , i3.2, ':00')
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'creation_date', &
                          trim(date_string))
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'modification_date', &
                          trim(date_string))
      call check_err(iret)

!     -------------------------------------------
!     writing global attributes from fort.15 file
!     -------------------------------------------
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'fort.15', &
                          '==== Input File Parameters (below) ====')
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'dt', dtdp)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'ihot', ihot)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'ics', ics)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'nolibf', nolibf)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'nolifa', nolifa)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'nolica', nolica)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'nolicat', nolicat)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'nwp', nwp)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'ncor', ncor)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'ntip', ntip)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'nws', nws)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'nramp', nramp)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'tau0', tau0)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'statim', statim)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'reftim', reftim)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'rnday', rnday)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'dramp', dramp)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'a00', a00)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'b00', b00)

      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'c00', c00)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'h0', h0)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'slam0', slam0deg)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'sfea0', sfea0deg)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'cf', cf)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'eslm', eslm)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'cori', cori)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'ntif', ntif)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'nbfr', nbfr)
      call check_err(iret)
      if (C3D .eqv. .true.) then
         iret = nf90_put_att(ncid, NF90_GLOBAL, 'iden', iden)
         call check_err(iret)
         iret = nf90_put_att(ncid, NF90_GLOBAL, 'islip', islip)
         call check_err(iret)
         iret = nf90_put_att(ncid, NF90_GLOBAL, 'kp', kp)
         call check_err(iret)
         iret = nf90_put_att(ncid, NF90_GLOBAL, 'z0s', z0s)
         call check_err(iret)
         iret = nf90_put_att(ncid, NF90_GLOBAL, 'z0b', z0b)
         call check_err(iret)
         iret = nf90_put_att(ncid, NF90_GLOBAL, 'alp1', alp1)
         call check_err(iret)
         iret = nf90_put_att(ncid, NF90_GLOBAL, 'alp2', alp2)
         call check_err(iret)
         iret = nf90_put_att(ncid, NF90_GLOBAL, 'alp3', alp3)
         call check_err(iret)
         iret = nf90_put_att(ncid, NF90_GLOBAL, 'igc', igc)
         call check_err(iret)
         iret = nf90_put_att(ncid, NF90_GLOBAL, 'ievc', ievc)
         call check_err(iret)
         iret = nf90_put_att(ncid, NF90_GLOBAL, 'evmin', evmin)
         call check_err(iret)
         iret = nf90_put_att(ncid, NF90_GLOBAL, 'evcon', evcon)
         call check_err(iret)
         if ((ievc == 50) .or. (ievc == 51)) then
            iret = nf90_put_att(ncid, NF90_GLOBAL, 'theta1', theta1)
            call check_err(iret)
            iret = nf90_put_att(ncid, NF90_GLOBAL, 'theta2', theta2)
            call check_err(iret)
         end if
         if (iden > 0) then
            iret = nf90_put_att(ncid, NF90_GLOBAL, 'nlsd', nlsd)
            call check_err(iret)
            iret = nf90_put_att(ncid, NF90_GLOBAL, 'nvsd', nvsd)
            call check_err(iret)
            iret = nf90_put_att(ncid, NF90_GLOBAL, 'nltd', nltd)
            call check_err(iret)
            iret = nf90_put_att(ncid, NF90_GLOBAL, 'nvtd', nvtd)
            call check_err(iret)
            iret = nf90_put_att(ncid, NF90_GLOBAL, 'alp4', alp4)
            call check_err(iret)
         end if
      end if
      if (inundationOutput .eqv. .true.) then
         iret = nf90_put_att(ncid, NF90_GLOBAL, 'inunThresh', inunThresh)
         call check_err(iret)
      end if

!-----------------------------------------------------------------------
   end subroutine defineMetaData
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     F U N C T I O N     M E T A   L E N G T H
!-----------------------------------------------------------------------
!     jgf49.29 Gets the length of the metadata line, depending on the
!     presence of a "!" in the line; if a "!" is present, it is used to
!     terminate the metadata ... if not, then the whole line is used (up
!     to 80 characters, or as declared in global.F).
!-----------------------------------------------------------------------
   function metalength(string)
      implicit none
      integer :: metalength
      character(*), intent(in) :: string
      metalength = index(string, "!") ! use the "!" as terminator if present
      if (metalength == 0) then
         ! there is no embedded "!" in the metadata line -- use the full line
         metalength = len_trim(string)
      else
         ! trim space between end of metadata and embedded "!"
         metalength = len_trim(string(1:metalength - 1))
      end if
!-----------------------------------------------------------------------
   end function metalength
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E   U P D A T E   M E T A   D A T A
!-----------------------------------------------------------------------
!     jgf49.17.02 Updates data that are common to all netCDF files
!     in the specified file.
!-----------------------------------------------------------------------
   subroutine updateMetaData(ncid, myFile)
      use GLOBAL, only: ihot, nramp, rnday, dramp, h0, rundes, runid
      use NodalAttributes, only: tau0, cf, eslm
      use GWCE, only: a00, b00, c00
      implicit none
      integer, intent(out) :: ncid
      type(fileData), intent(inout) :: myFile

      integer :: iret ! Error status return
!     date_string variables for time attribute
      character :: date_string*40
      character :: now_date*8
      character :: big_ben*10
      character :: zone*5
      integer :: values(8)

      LOG_SCOPE_TRACED("updateMetaData", NETCDFIO_TRACING)
!     Open existing NetCDF file
      iret = nf90_open(myFile%FILENAME, NF90_WRITE, ncid)
      call check_err(iret)

      iret = nf90_redef(ncid)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'rundes', trim(adjustl(rundes)))
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'runid', trim(adjustl(runid)))
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'ihot', ihot)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'nramp', nramp)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'tau0', tau0)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'rnday', rnday)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'dramp', dramp)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'a00', a00)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'b00', b00)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'c00', c00)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'h0', h0)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'cf', cf)
      call check_err(iret)
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'eslm', eslm)
      call check_err(iret)
      call date_and_time(now_date, big_ben, zone, values)
      write (date_string, 71) values(1), values(2), values(3), &
         values(5), values(6), values(7), (values(4))/60
71    format(I4, '-', I2.2, '-', i2.2, ' ', i2, ':', i2.2, ':', i2.2, ' ' &
             , i3.2, ':00')
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'modification_date', &
                          date_string)
      call check_err(iret)
      iret = nf90_enddef(ncid)
      call check_err(iret)
!     now close the updated netcdf file
      iret = nf90_close(ncid)
      call check_err(iret)

!-----------------------------------------------------------------------
   end subroutine updateMetaData
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E  R E A D  M E T A  D A T A
!-----------------------------------------------------------------------
!     zc50.74.05 - Give the writer processors a way to read the metadata
!     they need to update in the NetCDF file.
!-----------------------------------------------------------------------
   subroutine ReadMetaData(ncid, myFile)
!-----------------------------------------------------------------------
      use GLOBAL, only: ihot, nramp, rnday, dramp, h0, np_g, ne_g, rundes, runid
      use NodalAttributes, only: tau0, cf, eslm
      use GWCE, only: a00, b00, c00
      implicit none
      integer :: NCID, dimid_node, dimid_nele
      type(fileData), intent(INOUT) :: myFile

      LOG_SCOPE_TRACED("readMetaData", NETCDFIO_TRACING)
      call CHECK_ERR(nf90_open(myFile%FILENAME, NF90_NOWRITE, NCID))
      call CHECK_ERR(nf90_get_att(NCID, NF90_GLOBAL, 'rundes', rundes))
      call CHECK_ERR(nf90_get_att(NCID, NF90_GLOBAL, 'runid', runid))
      call CHECK_ERR(nf90_get_att(NCID, NF90_GLOBAL, 'ihot', ihot))
      call CHECK_ERR(nf90_get_att(NCID, NF90_GLOBAL, 'nramp', nramp))
      call CHECK_ERR(nf90_get_att(NCID, NF90_GLOBAL, 'tau0', tau0))
      call CHECK_ERR(nf90_get_att(NCID, NF90_GLOBAL, 'rnday', rnday))
      call CHECK_ERR(nf90_get_att(NCID, NF90_GLOBAL, 'dramp', dramp))
      call CHECK_ERR(nf90_get_att(NCID, NF90_GLOBAL, 'a00', a00))
      call CHECK_ERR(nf90_get_att(NCID, NF90_GLOBAL, 'b00', b00))
      call CHECK_ERR(nf90_get_att(NCID, NF90_GLOBAL, 'c00', c00))
      call CHECK_ERR(nf90_get_att(NCID, NF90_GLOBAL, 'h0', h0))
      call CHECK_ERR(nf90_get_att(NCID, NF90_GLOBAL, 'cf', cf))
      call CHECK_ERR(nf90_get_att(NCID, NF90_GLOBAL, 'eslm', eslm))
      call CHECK_ERR(nf90_inq_dimid(NCID, "node", dimid_node))
      call CHECK_ERR(nf90_inq_dimid(NCID, "nfaces", dimid_nele))
      call CHECK_ERR(nf90_inquire_dimension(ncid, dimid_node, len=NP_G))
      call CHECK_ERR(nf90_inquire_dimension(ncid, dimid_nele, len=NE_G))
      call CHECK_ERR(nf90_close(NCID))
!-----------------------------------------------------------------------
   end subroutine ReadMetaData
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E   P U T   U N I T S   A T T R I B U T E
!-----------------------------------------------------------------------
!     jgf49.17.02 Puts the right units label based on whether ADCIRC was
!     run with English units or SI units.
!-----------------------------------------------------------------------
   subroutine putUnitsAttribute(ncid, var_id, metric, english)
      use ADC_CONSTANTS, only: G
      implicit none
      integer, intent(in) :: ncid
      integer, intent(in) :: var_id
      character(*), intent(in) :: metric
      character(*), intent(in) :: english
      integer :: iret ! success or failure of netcdf call

      LOG_SCOPE_TRACED("putUnitsAttribute", NETCDFIO_TRACING)

      if (G < 11.d0) then
         iret = nf90_put_att(ncid, var_id, 'units', metric)
      else
         iret = nf90_put_att(ncid, var_id, 'units', english)
      end if
      call check_err(iret)

!-----------------------------------------------------------------------
   end subroutine putUnitsAttribute
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E   S E T   A D C I R C   P A R A M E T E R S
!-----------------------------------------------------------------------
!     jgf49.31 Called by adcprep to populate the global.F data normally
!     collected by read_input.F in adcirc. This is somewhat
!     convoluted, since the data are collected by adcprep in read14()
!     and read15() and used to populate pre_global.F, and then we
!     call this subroutine in the netcdf module to populate the global
!     and nodalattributes modules with the same data.
!
!     This twisted logic is the result of the unfortunate overlap between the
!     pre_global and global modules, among other things. Once adcprep
!     becomes integrated into adcirc, this silly subroutine will no longer
!     be needed. At the very least, adcprep should be made to populate
!     global, not pre_global with these data, in other words, both
!     adcirc and adcprep should use the global module.
!-----------------------------------------------------------------------
   subroutine setADCIRCParameters( &
      base_date_p, NE_G_p, &
      NBOU_p, NVEL_p, NOPE_p, NP_G_p, SLAM0_p, SFEA0_p, NBVV_p, &
      NVDLL_p, NBDV_p, NVELL_p, X_p, Y_p, &
      IBTYPE_p, IBTYPEE_p, SLAM_p, SFEA_p, &
      NODECODE_p, G_p, FileFmtRev_p, FileFmtMinor_p, &
      FileFmtMajor_p, im_p, nscoue_p, nscouv_p, &
      nscouc_p, nscoum_p, nscouge_p, nscougv_p, nscougc_p, &
      nscougw_p, NM_p, DP_p, RUNDES_p, AGRID_p, title_p, &
      institution_p, source_p, history_p, references_p, comments_p, &
      host_p, convention_p, contact_p, dtdp_p, ihot_p, ics_p, &
      nolifa_p, nolica_p, nolicat_p, ncor_p, ntip_p, nws_p, nramp_p, &
      statim_p, reftim_p, rnday_p, dramp_p, a00_p, b00_p, c00_p, &
      h0_p, cori_p, ntif_p, nbfr_p, myProc_p, nolibf_p, &
      nwp_p, tau0_p, cf_p, eslm_p, neta_p, &
      nabout_p, nscreen_p, &
      nfen_p, iden_p, islip_p, kp_p, z0s_p, z0b_p, theta1_p, theta2_p, &
      ievc_p, evmin_p, evcon_p, alp1_p, alp2_p, alp3_p, igc_p, nlsd_p, &
      nvsd_p, nltd_p, nvtd_p, alp4_p, C3D_p, RUNID_p)
      use SIZES, only: myproc
      use MESH, only: X, Y, SLAM, SFEA, NM, DP, ICS, &
                      SLAM0, SFEA0, AGRID
      use BOUNDARIES, only: NBOU, NVEL, NOPE, NBVV, NVDLL, NBDV, NVELL, &
                            IBTYPE, IBTYPEE, NETA
      use VERSION, only: FileFmtRev, FileFmtMinor, FileFmtMajor
      use ADC_CONSTANTS, only: G
      use GLOBAL, only: base_date, NODECODE, &
                        NP_G, NE_G, im, nscoue, nscouv, &
                        nscouc, nscoum, nscouge, nscougv, nscougc, nscougw, &
                        RUNDES, title, institution, source, history, &
                        references, comments, host, convention, contact, dtdp, ihot, &
                        nolifa, nolica, nolicat, ncor, ntip, nws, nramp, statim, &
                        reftim, rnday, dramp, h0, cori, ntif, nbfr, &
                        C3D, runid
      use GLOBAL_3DVS, only: &
         nfen, iden, islip, kp, z0s, z0b, theta1, theta2, &
         ievc, evmin, evcon, alp1, alp2, alp3, igc, nlsd, nvsd, nltd, &
         nvtd, alp4
      use GWCE, only: a00, b00, c00
      use NodalAttributes, only: nolibf, nwp, tau0, cf, eslm
      use mod_logging, only: nscreen, nabout, t_log_level
      implicit none

!     Declare the argument variables coming in from adcprep.
      character(80), intent(in) :: base_date_p
      integer, intent(in) :: NE_G_p
      integer, intent(in) :: NBOU_p
      integer, intent(in) :: NVEL_p
      integer, intent(in) :: NOPE_p
      integer, intent(in) :: NP_G_p
      real(8), intent(in) :: SLAM0_p
      real(8), intent(in) :: SFEA0_p
      integer, intent(in) :: NBVV_p(:, :)
      integer, intent(in) :: NVDLL_p(:)
      integer, intent(in) :: NBDV_p(:, :)
      integer, intent(in) :: NVELL_p(:)
      real(8), intent(in) :: X_p(:)
      real(8), intent(in) :: Y_p(:)
      integer, intent(in) :: IBTYPE_p(:)
      integer, intent(in) :: IBTYPEE_p(:)
      real(8), intent(in) :: SLAM_p(:)
      real(8), intent(in) :: SFEA_p(:)
      integer, intent(in) :: NODECODE_p(:)
      real(8), intent(in) :: G_p
      integer, intent(in) :: FileFmtRev_p
      integer, intent(in) :: FileFmtMinor_p
      integer, intent(in) :: FileFmtMajor_p
      integer, intent(in) :: im_p
      integer, intent(in) :: nscoue_p
      integer, intent(in) :: nscouv_p
      integer, intent(in) :: nscouc_p
      integer, intent(in) :: nscoum_p
      integer, intent(in) :: nscouge_p
      integer, intent(in) :: nscougv_p
      integer, intent(in) :: nscougc_p
      integer, intent(in) :: nscougw_p
      integer, intent(in) :: NM_p(:, :)
      real(8), intent(in) :: DP_p(:)
      character(80), intent(in) :: RUNDES_p
      character(80), intent(in) :: RUNID_p
      character(80), intent(in) :: AGRID_p
      character(80), intent(in) :: title_p
      character(80), intent(in) :: institution_p
      character(80), intent(in) :: source_p
      character(80), intent(in) :: history_p
      character(80), intent(in) :: references_p
      character(80), intent(in) :: comments_p
      character(80), intent(in) :: host_p
      character(80), intent(in) :: convention_p
      character(80), intent(in) :: contact_p
      real(8), intent(in) :: dtdp_p
      integer, intent(in) :: ihot_p
      integer, intent(in) :: ics_p
      integer, intent(in) :: nolifa_p
      integer, intent(in) :: nolica_p
      integer, intent(in) :: nolicat_p
      integer, intent(in) :: ncor_p
      integer, intent(in) :: ntip_p
      integer, intent(in) :: nws_p
      integer, intent(in) :: nramp_p
      real(8), intent(in) :: statim_p
      real(8), intent(in) :: reftim_p
      real(8), intent(in) :: rnday_p
      real(8), intent(in) :: dramp_p
      real(8), intent(in) :: a00_p
      real(8), intent(in) :: b00_p
      real(8), intent(in) :: c00_p
      real(8), intent(in) :: h0_p
      real(8), intent(in) :: cori_p
      integer, intent(in) :: ntif_p
      integer, intent(in) :: nbfr_p
      integer, intent(in) :: myProc_p
      integer, intent(in) :: nolibf_p
      integer, intent(in) :: nwp_p
      real(8), intent(in) :: tau0_p
      real(8), intent(in) :: cf_p
      real(8), intent(in) :: eslm_p
      integer, intent(in) :: neta_p
      integer, intent(in) :: nabout_p
      integer, intent(in) :: nscreen_p

      integer, intent(in) :: nfen_p
      integer, intent(in) :: iden_p
      integer, intent(in) :: islip_p
      real(8), intent(in) :: kp_p
      real(8), intent(in) :: z0s_p
      real(8), intent(in) :: z0b_p
      real(8), intent(in) :: theta1_p
      real(8), intent(in) :: theta2_p
      integer, intent(in) :: ievc_p
      real(8), intent(in) :: evmin_p
      real(8), intent(in) :: evcon_p
      real(8), intent(in) :: alp1_p
      real(8), intent(in) :: alp2_p
      real(8), intent(in) :: alp3_p
      integer, intent(in) :: igc_p
      real(8), intent(in) :: nlsd_p
      real(8), intent(in) :: nvsd_p
      real(8), intent(in) :: nltd_p
      real(8), intent(in) :: nvtd_p
      real(8), intent(in) :: alp4_p
      logical, intent(in) :: C3D_p

      LOG_SCOPE_TRACED("setADCIRCParameters", NETCDFIO_TRACING)

      base_date = base_date_p
      NE_G = NE_G_p
      NBOU = NBOU_p
      NVEL = NVEL_p
      NOPE = NOPE_p
      NP_G = NP_G_p
      SLAM0 = SLAM0_p
      SFEA0 = SFEA0_p
      allocate (NBVV(NBOU_p, 0:NVEL_p))
      NBVV = NBVV_p
      allocate (NVDLL(NOPE_p))
      NVDLL = NVDLL_p
      allocate (NBDV(NOPE_p, NETA_p))
      NBDV = NBDV_p
      allocate (NVELL(NBOU_p))
      NVELL = NVELL_p
      allocate (X(NP_G_p))
      X = X_p
      allocate (Y(NP_G_p))
      Y = Y_p
      allocate (IBTYPEE(NOPE_p))
      IBTYPEE = IBTYPEE_p
      allocate (IBTYPE(NBOU_p))
      IBTYPE = IBTYPE_p
      allocate (SLAM(NP_G_p))
      SLAM = SLAM_p
      allocate (SFEA(NP_G_p))
      SFEA = SFEA_p
      allocate (NODECODE(NP_G_p))
      NODECODE = NODECODE_p
      G = G_p
      FileFmtRev = FileFmtRev_p
      FileFmtMinor = FileFmtMinor_p
      FileFmtMajor = FileFmtMajor_p
      im = im_p
      nscoue = nscoue_p
      nscouv = nscouv_p
      nscouc = nscouc_p
      nscoum = nscoum_p
      nscouge = nscouge_p
      nscougv = nscougv_p
      nscougc = nscougc_p
      nscougw = nscougw_p
      allocate (NM(NE_G_p, 3))
      NM = NM_p
      allocate (DP(NP_G_p))
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
      c00 = c00_p
      h0 = h0_p
      cori = cori_p
      ntif = ntif_p
      nbfr = nbfr_p
      myProc = myProc_p
      nolibf = nolibf_p
      nwp = nwp_p
      tau0 = tau0_p
      cf = cf_p
      eslm = eslm_p
      neta = neta_p
      nabout = t_log_level(nabout_p)
      nscreen = nscreen_p
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

!-----------------------------------------------------------------------
   end subroutine setADCIRCParameters
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E   F R E E    N E T C D F   C O O R D
!-----------------------------------------------------------------------
!     jgf49.43.05 Frees memory allocated for mesh and boundary data
!     in NetCDF files.
!-----------------------------------------------------------------------
   subroutine freeNetCDFCoord()
      implicit none

      LOG_SCOPE_TRACED("freeNetCDFCoord.", NETCDFIO_TRACING)

      if (adcircMesh%initialized .eqv. .true.) then
         deallocate (adcircMesh%xnc)
         deallocate (adcircMesh%ync)
         deallocate (adcircMesh%nvdllnc)
         deallocate (adcircMesh%ibtypeenc)
         deallocate (adcircMesh%ibtypenc)
         deallocate (adcircMesh%nvellnc)
         deallocate (adcircMesh%nbdvnc)
         deallocate (adcircMesh%nbvvnc)
         deallocate (adcircMesh%ibconnnc)
         deallocate (adcircMesh%element)
         deallocate (adcircMesh%nmnc)
         adcircMesh%initialized = .false.
      end if

!-----------------------------------------------------------------------
   end subroutine freeNetCDFCoord
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module NETCDFIO
!-----------------------------------------------------------------------
